"""
UROPA Visualization Module — Genomic Location of Features

Generates summary plots for UROPA annotation results:
- Annotation status (Annotated vs Not annotated)
- Genomic location distribution (bar plots)
- Split by feature, query, and feature x query matrix
"""

import json
import logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LOCATION_ORDER = [
    "PeakInsideFeature",
    "FeatureInsidePeak",
    "Upstream",
    "Downstream",
    "OverlapStart",
    "OverlapEnd",
]

# Auto-generate colorblind-friendly palette from seaborn
_cb_palette = sns.color_palette("colorblind", len(LOCATION_ORDER))
LOCATION_COLORS = {loc: col for loc, col in zip(LOCATION_ORDER, _cb_palette)}
LOCATION_COLORS["Other"] = "#CCCCCC"

STATUS_COLORS = {
    "Annotated": "#55A868",
    "Not annotated": "#C44E52",
}

FIGURE_DPI = 150
MAX_SUBPLOT_COLS = 3
DEFAULT_OTHER_THRESHOLD = 0.1  # percent


# ---------------------------------------------------------------------------
# Data Loading
# ---------------------------------------------------------------------------

def _load_hits(filepath):
    """Read a UROPA hits TSV file, return DataFrame with key columns.

    Parameters
    ----------
    filepath : str
        Path to a ``_finalhits.txt`` or ``_allhits.txt`` file.

    Returns
    -------
    pd.DataFrame
        Columns: ``feature``, ``relative_location``, ``name``.
    """
    df = pd.read_csv(
        filepath, sep="\t",
        usecols=["feature", "relative_location", "name"],
        dtype=str, na_filter=False,
    )
    return df


def _load_config(config_file):
    """Read the UROPA JSON config file.

    Returns
    -------
    dict
        Parsed configuration dictionary.
    """
    with open(config_file) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Data Transformation
# ---------------------------------------------------------------------------

def _compute_annotation_status(df):
    """Compute annotated vs not-annotated counts.

    Parameters
    ----------
    df : pd.DataFrame
        Full DataFrame including NA rows.

    Returns
    -------
    pd.DataFrame
        Columns: ``status``, ``count``, ``percentage``.
    """
    total = len(df)
    na_count = (df["feature"] == "NA").sum()
    annotated_count = total - na_count

    records = [
        {"status": "Annotated", "count": annotated_count,
         "percentage": round(annotated_count / total * 100, 1) if total > 0 else 0.0},
        {"status": "Not annotated", "count": na_count,
         "percentage": round(na_count / total * 100, 1) if total > 0 else 0.0},
    ]
    return pd.DataFrame(records)


def _filter_annotated(df):
    """Remove rows where feature is NA (unannotated peaks).

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame without NA rows.
    """
    return df[df["feature"] != "NA"].copy()


def _compute_location_counts(df, group_cols=None):
    """Count occurrences of each relative_location, optionally grouped.

    Parameters
    ----------
    df : pd.DataFrame
        Annotated-only DataFrame (no NA rows).
    group_cols : list of str or None
        Grouping columns. ``None`` for overall, e.g. ``["feature"]``,
        ``["name"]``, or ``["feature", "name"]``.

    Returns
    -------
    pd.DataFrame
        Columns: ``[*group_cols, "relative_location", "count", "percentage"]``.
        Percentages sum to 100 within each group.
    """
    if group_cols:
        grouped = (df.groupby(group_cols + ["relative_location"])
                     .size()
                     .reset_index(name="count"))
        group_totals = grouped.groupby(group_cols)["count"].transform("sum")
        grouped["percentage"] = (grouped["count"] / group_totals * 100).round(1)
    else:
        grouped = (df.groupby("relative_location")
                     .size()
                     .reset_index(name="count"))
        total = grouped["count"].sum()
        grouped["percentage"] = (
            (grouped["count"] / total * 100).round(1) if total > 0 else 0.0
        )

    return grouped


def _apply_other_threshold(counts_df, threshold=DEFAULT_OTHER_THRESHOLD, group_cols=None):
    """Group small categories into 'Other'.

    Categories with percentage below *threshold* within their group are
    merged into a single 'Other' row.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Output of :func:`_compute_location_counts`.
    threshold : float
        Minimum percentage for a category to remain individual (default 0.1).
    group_cols : list of str or None
        Grouping columns (same as used in ``_compute_location_counts``).

    Returns
    -------
    main_df : pd.DataFrame
        Same structure, small categories replaced by one 'Other' row per group.
    other_details_df : pd.DataFrame
        The original small categories that were merged.  Empty DataFrame if
        none were below threshold.
    """
    if group_cols:
        main_parts = []
        other_parts = []

        for group_key, group_df in counts_df.groupby(group_cols):
            small = group_df[group_df["percentage"] < threshold]
            large = group_df[group_df["percentage"] >= threshold].copy()

            if not small.empty:
                other_parts.append(small.copy())
                # Build the "Other" row with correct group columns
                if isinstance(group_key, tuple):
                    other_row = {col: val for col, val in zip(group_cols, group_key)}
                else:
                    other_row = {group_cols[0]: group_key}
                other_row["relative_location"] = "Other"
                other_row["count"] = int(small["count"].sum())
                other_row["percentage"] = round(small["percentage"].sum(), 1)
                large = pd.concat(
                    [large, pd.DataFrame([other_row])], ignore_index=True,
                )

            main_parts.append(large)

        main_df = (pd.concat(main_parts, ignore_index=True)
                   if main_parts else counts_df.copy())
        other_details_df = (pd.concat(other_parts, ignore_index=True)
                            if other_parts else pd.DataFrame())
    else:
        small = counts_df[counts_df["percentage"] < threshold]
        large = counts_df[counts_df["percentage"] >= threshold].copy()

        if not small.empty:
            other_details_df = small.copy()
            other_row = {
                "relative_location": "Other",
                "count": int(small["count"].sum()),
                "percentage": round(small["percentage"].sum(), 1),
            }
            main_df = pd.concat(
                [large, pd.DataFrame([other_row])], ignore_index=True,
            )
        else:
            main_df = large
            other_details_df = pd.DataFrame()

    return main_df, other_details_df


# ---------------------------------------------------------------------------
# Plot Rendering — Single Plots
# ---------------------------------------------------------------------------

def _get_bar_colors(categories):
    """Return a list of colors matching *categories* using LOCATION_COLORS."""
    return [LOCATION_COLORS.get(cat, "#CCCCCC") for cat in categories]


def _annotate_bars(ax, orientation="horizontal"):
    """Add count + percentage labels to bars.

    Expects bar containers on *ax*.  Labels are placed to the right of
    horizontal bars or on top of vertical bars.
    """
    for container in ax.containers:
        if orientation == "horizontal":
            ax.bar_label(container, fmt="%g", padding=3, fontsize=8)
        else:
            ax.bar_label(container, fmt="%g", padding=3, fontsize=8)


def _order_counts(counts_df, order=None):
    """Sort *counts_df* rows so that ``relative_location`` follows *order*.

    Categories not in *order* (e.g. 'Other') are appended at the end.
    Missing categories from *order* are skipped.
    """
    if order is None:
        order = LOCATION_ORDER
    # Build ordered list: keep only categories present in data
    present = counts_df["relative_location"].tolist()
    ordered = [loc for loc in order if loc in present]
    # Append anything not in LOCATION_ORDER (e.g. "Other")
    ordered += [loc for loc in present if loc not in ordered]
    cat_type = pd.CategoricalDtype(categories=ordered, ordered=True)
    counts_df = counts_df.copy()
    counts_df["relative_location"] = counts_df["relative_location"].astype(cat_type)
    return counts_df.sort_values("relative_location").reset_index(drop=True)


def _plot_annotation_status(status_df, title="Annotation Status"):
    """Horizontal bar plot showing Annotated vs Not annotated.

    Parameters
    ----------
    status_df : pd.DataFrame
        Output of :func:`_compute_annotation_status`.
    title : str
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(8, 3))
    colors = [STATUS_COLORS[s] for s in status_df["status"]]

    bars = ax.barh(status_df["status"], status_df["count"], color=colors,
                   edgecolor="white", height=0.5)

    # Add labels: count (percentage%)
    for bar, row in zip(bars, status_df.itertuples()):
        label = f"{row.count:,}  ({row.percentage}%)"
        ax.text(bar.get_width() + ax.get_xlim()[1] * 0.01, bar.get_y() + bar.get_height() / 2,
                label, va="center", ha="left", fontsize=10)

    ax.set_xlabel("Number of peaks")
    ax.set_title(title, fontweight="bold", fontsize=12)
    ax.set_xlim(right=ax.get_xlim()[1] * 1.25)  # space for labels
    sns.despine(ax=ax, left=True)
    ax.tick_params(left=False)
    fig.tight_layout()
    return fig


def _plot_genomic_location(counts_df, other_details_df, title="Genomic Location Distribution"):
    """Horizontal bar plot of relative_location categories.

    If *other_details_df* is non-empty, a secondary subplot shows the
    breakdown of the 'Other' category.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Location counts with 'Other' row (output of :func:`_apply_other_threshold`).
    other_details_df : pd.DataFrame
        Small categories that were merged into 'Other'.
    title : str
        Plot title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    counts_df = _order_counts(counts_df)
    has_other = not other_details_df.empty

    if has_other:
        fig, (ax_main, ax_other) = plt.subplots(
            1, 2, figsize=(12, 5),
            gridspec_kw={"width_ratios": [3, 1]},
        )
    else:
        fig, ax_main = plt.subplots(figsize=(8, 5))

    # Main bar plot
    colors = _get_bar_colors(counts_df["relative_location"])
    bars = ax_main.barh(
        counts_df["relative_location"].astype(str),
        counts_df["count"],
        color=colors, edgecolor="white",
    )

    # Labels: count (percentage%)
    for bar, row in zip(bars, counts_df.itertuples()):
        label = f"{row.count:,}  ({row.percentage}%)"
        ax_main.text(
            bar.get_width() + ax_main.get_xlim()[1] * 0.01,
            bar.get_y() + bar.get_height() / 2,
            label, va="center", ha="left", fontsize=9,
        )

    ax_main.set_xlabel("Count")
    ax_main.set_title(title, fontweight="bold", fontsize=12)
    ax_main.set_xlim(right=ax_main.get_xlim()[1] * 1.3)
    ax_main.invert_yaxis()  # top-to-bottom order
    sns.despine(ax=ax_main, left=True)
    ax_main.tick_params(left=False)

    # Other breakdown subplot
    if has_other:
        other_details_df = _order_counts(other_details_df)
        other_colors = _get_bar_colors(other_details_df["relative_location"])
        other_bars = ax_other.barh(
            other_details_df["relative_location"].astype(str),
            other_details_df["count"],
            color=other_colors, edgecolor="white",
        )
        for bar, row in zip(other_bars, other_details_df.itertuples()):
            label = f"{row.count:,}  ({row.percentage}%)"
            ax_other.text(
                bar.get_width() + ax_other.get_xlim()[1] * 0.01,
                bar.get_y() + bar.get_height() / 2,
                label, va="center", ha="left", fontsize=8,
            )
        ax_other.set_title("Other Breakdown", fontsize=10)
        ax_other.set_xlabel("Count")
        ax_other.set_xlim(right=ax_other.get_xlim()[1] * 1.4)
        ax_other.invert_yaxis()
        sns.despine(ax=ax_other, left=True)
        ax_other.tick_params(left=False)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Plot Rendering — Group & Matrix Plots
# ---------------------------------------------------------------------------

def _plot_single_barh(ax, group_counts, fontsize=8):
    """Draw a horizontal bar plot on *ax* for a single group's location counts.

    Helper used by ``_plot_genomic_location_by_group`` and
    ``_plot_genomic_location_matrix``.
    """
    group_counts = _order_counts(group_counts)
    colors = _get_bar_colors(group_counts["relative_location"])
    bars = ax.barh(
        group_counts["relative_location"].astype(str),
        group_counts["count"],
        color=colors, edgecolor="white",
    )
    for bar, row in zip(bars, group_counts.itertuples()):
        label = f"{row.count:,} ({row.percentage}%)"
        ax.text(
            bar.get_width() + ax.get_xlim()[1] * 0.01,
            bar.get_y() + bar.get_height() / 2,
            label, va="center", ha="left", fontsize=fontsize,
        )
    ax.set_xlim(right=ax.get_xlim()[1] * 1.4)
    ax.invert_yaxis()
    sns.despine(ax=ax, left=True)
    ax.tick_params(left=False)


def _plot_genomic_location_by_group(counts_df, other_details_df, group_col,
                                     title="Genomic Location by Group"):
    """Grid of bar plots, one subplot per unique value of *group_col*.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Location counts grouped by *group_col* (+ 'Other' rows).
    other_details_df : pd.DataFrame
        Small categories that were merged into 'Other'.
    group_col : str
        Column name to split on (``"feature"`` or ``"name"``).
    title : str
        Overall figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    groups = sorted(counts_df[group_col].unique())
    has_other = not other_details_df.empty

    n_panels = len(groups) + (1 if has_other else 0)
    ncols = min(n_panels, MAX_SUBPLOT_COLS)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(5 * ncols, 3.5 * nrows),
        squeeze=False,
    )
    fig.suptitle(title, fontweight="bold", fontsize=13, y=1.02)

    # Plot one bar chart per group
    for idx, group_val in enumerate(groups):
        row_i, col_i = divmod(idx, ncols)
        ax = axes[row_i][col_i]
        subset = counts_df[counts_df[group_col] == group_val].copy()
        _plot_single_barh(ax, subset, fontsize=7)
        ax.set_title(group_val, fontsize=10, fontweight="bold")

    # Other breakdown panel
    if has_other:
        row_i, col_i = divmod(len(groups), ncols)
        ax_other = axes[row_i][col_i]
        other_details_df = _order_counts(other_details_df)
        # Show grouped: one bar per original small category, color-coded
        other_colors = _get_bar_colors(other_details_df["relative_location"])
        labels = [
            f"{r.relative_location} ({r[group_col]})"
            for _, r in other_details_df.iterrows()
        ]
        bars = ax_other.barh(labels, other_details_df["count"],
                             color=other_colors, edgecolor="white")
        for bar, row in zip(bars, other_details_df.itertuples()):
            ax_other.text(
                bar.get_width() + ax_other.get_xlim()[1] * 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{row.count:,} ({row.percentage}%)",
                va="center", ha="left", fontsize=6,
            )
        ax_other.set_xlim(right=ax_other.get_xlim()[1] * 1.5)
        ax_other.invert_yaxis()
        ax_other.set_title("Other Breakdown", fontsize=10)
        sns.despine(ax=ax_other, left=True)
        ax_other.tick_params(left=False)

    # Hide unused axes
    for idx in range(n_panels, nrows * ncols):
        row_i, col_i = divmod(idx, ncols)
        axes[row_i][col_i].set_visible(False)

    fig.tight_layout()
    return fig


def _plot_genomic_location_matrix(counts_df, other_details_df,
                                   row_col, col_col,
                                   title="Genomic Location Matrix"):
    """Matrix of bar plots: *row_col* as rows, *col_col* as columns.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Location counts grouped by ``[row_col, col_col]``.
    other_details_df : pd.DataFrame
        Small categories merged into 'Other'.
    row_col : str
        Column for matrix rows (typically ``"feature"``).
    col_col : str
        Column for matrix columns (typically ``"name"``).
    title : str
        Overall figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    row_vals = sorted(counts_df[row_col].unique())
    col_vals = sorted(counts_df[col_col].unique())
    nrows = len(row_vals)
    ncols = len(col_vals)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(4 * ncols, 3 * nrows),
        squeeze=False,
    )
    fig.suptitle(title, fontweight="bold", fontsize=13, y=1.02)

    for ri, rv in enumerate(row_vals):
        for ci, cv in enumerate(col_vals):
            ax = axes[ri][ci]
            subset = counts_df[
                (counts_df[row_col] == rv) & (counts_df[col_col] == cv)
            ].copy()

            if subset.empty:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        fontsize=10, color="gray",
                        transform=ax.transAxes)
                ax.set_facecolor("#f5f5f5")
                ax.set_xticks([])
                ax.set_yticks([])
                sns.despine(ax=ax, left=True, bottom=True)
            else:
                _plot_single_barh(ax, subset, fontsize=6)

            # Row label on leftmost column
            if ci == 0:
                ax.set_ylabel(rv, fontsize=9, fontweight="bold")
            # Column label on top row
            if ri == 0:
                ax.set_title(cv, fontsize=9, fontweight="bold")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Output Helpers
# ---------------------------------------------------------------------------

def _save_figure(fig, pdf_pages, png_path, dpi=FIGURE_DPI):
    """Save *fig* to the open *pdf_pages* and as a PNG, then close it."""
    pdf_pages.savefig(fig, bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _create_cover_page(config, total_peaks, annotated_count):
    """Create a text-based cover page summarizing the UROPA run.

    Parameters
    ----------
    config : dict
        Parsed UROPA JSON config.
    total_peaks : int
        Total number of peaks.
    annotated_count : int
        Number of annotated peaks.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")

    lines = ["UROPA Summary", ""]
    lines.append(f"Total peaks: {total_peaks:,}")
    lines.append(f"Annotated: {annotated_count:,}  "
                 f"({round(annotated_count / total_peaks * 100, 1) if total_peaks else 0}%)")
    lines.append("")

    queries = config.get("queries", [])
    lines.append(f"Queries: {len(queries)}")
    for i, q in enumerate(queries):
        name = q.get("name", f"query_{i}")
        feature = q.get("feature", "any")
        distance = q.get("distance", "default")
        lines.append(f"  [{i+1}] {name}: feature={feature}, distance={distance}")

    text = "\n".join(lines)
    ax.text(0.05, 0.95, text, transform=ax.transAxes,
            fontsize=11, verticalalignment="top", fontfamily="monospace")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main Orchestrator
# ---------------------------------------------------------------------------

def _generate_plots_for_hits(df, hit_type, output_prefix, pdf, threshold, logger):
    """Generate all genomic location plots for one hit type (finalhits or allhits).

    Parameters
    ----------
    df : pd.DataFrame
        Full DataFrame (including NA rows) for this hit type.
    hit_type : str
        ``"finalhits"`` or ``"allhits"``.
    output_prefix : str
        Base path for output files.
    pdf : PdfPages
        Open PDF file to append pages to.
    threshold : float
        Other-category threshold percentage.
    logger : logging.Logger or None

    Returns
    -------
    bool
        True if at least one plot was generated.
    """
    annotated = _filter_annotated(df)

    if annotated.empty:
        if logger:
            logger.info("No annotated peaks in %s — skipping location plots.", hit_type)
        # Text-only page
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, f"No annotated peaks in {hit_type}.",
                ha="center", va="center", fontsize=14)
        _save_figure(fig, pdf, f"{output_prefix}_{hit_type}_no_annotations.png")
        return False

    n_features = annotated["feature"].nunique()
    n_queries = annotated["name"].nunique()

    if logger:
        logger.info("Generating %s plots: %d features, %d queries.",
                     hit_type, n_features, n_queries)

    # Plot: Overall genomic location
    counts = _compute_location_counts(annotated)
    main, other = _apply_other_threshold(counts, threshold)
    fig = _plot_genomic_location(main, other,
                                  title=f"Genomic Location Distribution ({hit_type})")
    _save_figure(fig, pdf, f"{output_prefix}_{hit_type}_location_overall.png")

    # Plot: By feature (only if >1 feature)
    if n_features > 1:
        counts_f = _compute_location_counts(annotated, group_cols=["feature"])
        main_f, other_f = _apply_other_threshold(counts_f, threshold, group_cols=["feature"])
        fig = _plot_genomic_location_by_group(
            main_f, other_f, "feature",
            title=f"Genomic Location by Feature ({hit_type})",
        )
        _save_figure(fig, pdf, f"{output_prefix}_{hit_type}_location_by_feature.png")

    # Plot: By query (only if >1 query)
    if n_queries > 1:
        counts_n = _compute_location_counts(annotated, group_cols=["name"])
        main_n, other_n = _apply_other_threshold(counts_n, threshold, group_cols=["name"])
        fig = _plot_genomic_location_by_group(
            main_n, other_n, "name",
            title=f"Genomic Location by Query ({hit_type})",
        )
        _save_figure(fig, pdf, f"{output_prefix}_{hit_type}_location_by_query.png")

    # Plot: Matrix feature x query (only if >1 of each)
    if n_features > 1 and n_queries > 1:
        counts_b = _compute_location_counts(annotated, group_cols=["feature", "name"])
        main_b, other_b = _apply_other_threshold(
            counts_b, threshold, group_cols=["feature", "name"],
        )
        fig = _plot_genomic_location_matrix(
            main_b, other_b, "feature", "name",
            title=f"Genomic Location: Feature x Query ({hit_type})",
        )
        _save_figure(fig, pdf, f"{output_prefix}_{hit_type}_location_matrix.png")

    return True


def generate_summary(finalhits_file, allhits_file=None, config_file=None,
                     output_prefix="uropa_summary", other_threshold=DEFAULT_OTHER_THRESHOLD,
                     logger=None):
    """Generate UROPA summary visualization.

    Produces a multi-page PDF and individual PNG files.

    Parameters
    ----------
    finalhits_file : str
        Path to ``_finalhits.txt``.
    allhits_file : str or None
        Path to ``_allhits.txt`` (optional).
    config_file : str or None
        Path to UROPA JSON config (for cover page).
    output_prefix : str
        Base path for output files. PDF will be ``{output_prefix}_summary.pdf``.
    other_threshold : float
        Minimum percentage for a category to appear individually (default 0.1).
    logger : logging.Logger or None
    """
    if logger:
        logger.info("Creating summary visualization...")

    # Load data
    df_final = _load_hits(finalhits_file)
    df_all = _load_hits(allhits_file) if allhits_file else None

    config = _load_config(config_file) if config_file else {}

    # Compute annotation status from finalhits (one row per peak)
    status_df = _compute_annotation_status(df_final)
    total_peaks = len(df_final)
    annotated_count = int(status_df.loc[status_df["status"] == "Annotated", "count"].iloc[0])

    pdf_path = f"{output_prefix}_summary.pdf"

    with PdfPages(pdf_path) as pdf:
        # Cover page
        if config:
            fig_cover = _create_cover_page(config, total_peaks, annotated_count)
            _save_figure(fig_cover, pdf, f"{output_prefix}_summary_cover.png")

        # Annotation status (always from finalhits)
        fig_status = _plot_annotation_status(status_df, title="Annotation Status")
        _save_figure(fig_status, pdf, f"{output_prefix}_summary_annotation_status.png")

        # Finalhits plots
        _generate_plots_for_hits(
            df_final, "finalhits", output_prefix, pdf, other_threshold, logger,
        )

        # Allhits plots (if available)
        if df_all is not None:
            _generate_plots_for_hits(
                df_all, "allhits", output_prefix, pdf, other_threshold, logger,
            )

    if logger:
        logger.info("Summary saved to %s", pdf_path)
