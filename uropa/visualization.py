"""
visualization.py: Distance visualization for UROPA summary output.

Generates distance distribution plots (KDE or histogram) from UROPA
annotation output files. Called from uropa.py when --summary is set.
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

COLORS = plt.cm.tab10.colors


def load_uropa_output(filepath):
    """Read UROPA TSV output and drop rows where distance is NA.

    Parameters
    ----------
    filepath : str
        Path to a UROPA output TSV file (e.g. ``*_allhits.txt`` or
        ``*_finalhits.txt``).

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ``distance``, ``peak_start``, ``peak_end``,
        ``feat_start``, and ``feat_end`` cast to float. Rows with missing
        or ``'NA'`` distance values are removed.

    Raises
    ------
    FileNotFoundError
        If *filepath* does not exist (raised by ``pd.read_csv``).
    KeyError
        If the required column ``distance`` is missing from the file.
    """
    df = pd.read_csv(filepath, sep='\t')
    df = df[df['distance'].notna() & (df['distance'] != 'NA')]
    df['distance'] = df['distance'].astype(float)
    df['peak_start'] = df['peak_start'].astype(float)
    df['peak_end'] = df['peak_end'].astype(float)
    df['feat_start'] = df['feat_start'].astype(float)
    df['feat_end'] = df['feat_end'].astype(float)
    return df


def reconstruct_signed_distance(df):
    """Add a ``signed_distance`` column to the DataFrame.

    Negative values indicate upstream (5') and positive values indicate
    downstream (3') relative to the feature strand. The anchor position is
    resolved using the same strand-aware logic as ``annotation.py:67-69``.

    Parameters
    ----------
    df : pandas.DataFrame
        UROPA output DataFrame. Must contain the columns ``peak_start``,
        ``peak_end``, ``feat_start``, ``feat_end``, ``feat_strand``, and
        ``feat_anchor``.

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with an added ``signed_distance`` column (float).

    Notes
    -----
    For minus-strand features the ``start`` anchor maps to ``feat_end`` and
    the ``end`` anchor maps to ``feat_start``, matching the convention in
    ``annotation.py``.
    """
    peak_center = (df['peak_start'] + df['peak_end']) / 2.0
    is_minus = df['feat_strand'] == '-'

    anchor_pos = np.where(
        df['feat_anchor'] == 'center',
        (df['feat_start'] + df['feat_end']) / 2.0,
        np.where(
            df['feat_anchor'] == 'start',
            np.where(is_minus, df['feat_end'], df['feat_start']),
            np.where(is_minus, df['feat_start'], df['feat_end'])  # 'end'
        )
    )

    raw_distance = peak_center - anchor_pos
    df['signed_distance'] = np.where(is_minus, -raw_distance, raw_distance)
    return df


def plot_distance_kde(df, facet_col, distance_col, log_scale, ax, xlim=None):
    """Plot a KDE density curve with colored overlays per facet group.

    Each unique value in *facet_col* is drawn as a separate filled KDE curve.
    Groups with fewer than 2 data points or singular covariance matrices are
    silently skipped.

    Parameters
    ----------
    df : pandas.DataFrame
        UROPA output DataFrame containing at least *facet_col* and
        *distance_col*.
    facet_col : str
        Column name used to split the data into groups (e.g. ``'feature'``,
        ``'feat_anchor'``, or ``'name'``).
    distance_col : str
        Column name containing the distance values to plot. Typically
        ``'distance'`` (absolute) or ``'signed_distance'`` (directional).
    log_scale : bool
        If ``True``, apply symmetric-log scale (``symlog``) to the y-axis
        with ``linthresh=1e-6``.
    ax : matplotlib.axes.Axes
        Axes object to draw on.
    xlim : tuple of (float, float) or None, optional
        If provided, the KDE is evaluated over this ``(x_min, x_max)``
        range instead of the data range. Default is ``None``.

    Returns
    -------
    None
        Draws directly on the provided *ax*.
    """
    # Use xlim for KDE evaluation range if provided, otherwise use full data range
    all_data = df[distance_col].dropna().values
    if xlim is not None:
        x_min, x_max = xlim
    else:
        x_min, x_max = all_data.min(), all_data.max()
    x_range = np.linspace(x_min, x_max, 500)

    groups = sorted(df[facet_col].unique())
    for i, group in enumerate(groups):
        data = df.loc[df[facet_col] == group, distance_col].dropna().values
        if len(data) < 2:
            continue
        try:
            kde = gaussian_kde(data)
        except np.linalg.LinAlgError:
            continue
        ax.fill_between(x_range, kde(x_range), alpha=0.4,
                        color=COLORS[i % len(COLORS)], label=str(group))
        ax.plot(x_range, kde(x_range), color=COLORS[i % len(COLORS)], linewidth=1)
    ax.set_ylabel('Density')
    ax.legend()
    if log_scale:
        ax.set_yscale('symlog', linthresh=1e-6)


def plot_distance_hist(df, facet_col, distance_col, log_scale, ax, xlim=None):
    """Plot a histogram with colored overlays per facet group.

    Each unique value in *facet_col* is drawn as a separate step-filled
    histogram layer.

    Parameters
    ----------
    df : pandas.DataFrame
        UROPA output DataFrame containing at least *facet_col* and
        *distance_col*.
    facet_col : str
        Column name used to split the data into groups (e.g. ``'feature'``,
        ``'feat_anchor'``, or ``'name'``).
    distance_col : str
        Column name containing the distance values to plot. Typically
        ``'distance'`` (absolute) or ``'signed_distance'`` (directional).
    log_scale : bool
        If ``True``, apply symmetric-log scale (``symlog``) to the y-axis
        with ``linthresh=1``.
    ax : matplotlib.axes.Axes
        Axes object to draw on.
    xlim : tuple of (float, float) or None, optional
        If provided, used as the ``range`` parameter for the histogram
        binning. Default is ``None``.

    Returns
    -------
    None
        Draws directly on the provided *ax*.
    """
    groups = sorted(df[facet_col].unique())
    data_list = [df.loc[df[facet_col] == g, distance_col].dropna().values for g in groups]
    colors = [COLORS[i % len(COLORS)] for i in range(len(groups))]
    hist_range = tuple(xlim) if xlim is not None else None
    ax.hist(data_list, bins='auto', histtype='stepfilled', alpha=0.5,
            color=colors, label=[str(g) for g in groups], range=hist_range)
    ax.set_ylabel('Count')
    ax.legend()
    if log_scale:
        ax.set_yscale('symlog', linthresh=1)


def generate_summary(output_prefix, summary_input, plot_type, facet_by,
                     distance_mode, log_scale, logger):
    """Generate a summary PDF with distance distribution plots.

    This is the main entry point called from ``uropa.py`` when ``--summary``
    is set. All logic is wrapped in a try/except so that visualization
    failures are logged but never crash the UROPA pipeline.

    Parameters
    ----------
    output_prefix : str
        Path prefix for UROPA output files (e.g. ``'results/my_run'``).
        The input file is derived as ``output_prefix + '_finalhits.txt'`` or
        ``output_prefix + '_allhits.txt'``, and the PDF is saved as
        ``output_prefix + '_summary.pdf'``.
    summary_input : {'finalhits', 'allhits'}
        Which UROPA output to visualize. ``'finalhits'`` produces a single
        plot; ``'allhits'`` produces one subplot per query.
    plot_type : {'kde', 'hist'}
        Visualization type. ``'kde'`` draws kernel density estimates,
        ``'hist'`` draws histograms.
    facet_by : {'feature', 'anchor', 'query'}
        How to group (color) data within each plot. Mapped internally to
        DataFrame columns: ``'feature'`` -> ``'feature'``,
        ``'anchor'`` -> ``'feat_anchor'``, ``'query'`` -> ``'name'``.
    distance_mode : {'absolute', 'directional'}
        ``'absolute'`` plots unsigned distances. ``'directional'`` reconstructs
        strand-relative signed distances (negative = upstream / 5',
        positive = downstream / 3') and centers the x-axis symmetrically
        around zero.
    log_scale : bool
        If ``True``, apply symmetric-log scale to the y-axis.
    logger : logging.Logger
        Logger instance for info/warning/error messages.

    Returns
    -------
    None
        Saves a PDF to ``output_prefix + '_summary.pdf'``. Returns
        ``None`` in all cases, including on error.

    Notes
    -----
    The function returns early (with a warning) if the input file is missing
    or contains no annotated peaks. Any unexpected exception is caught,
    logged at ERROR level, and a full traceback is logged at DEBUG level.
    """
    suffix = '_finalhits.txt' if summary_input == 'finalhits' else '_allhits.txt'
    filepath = output_prefix + suffix

    try:
        if not os.path.exists(filepath):
            logger.warning("Summary input file not found: %s", filepath)
            return

        df = load_uropa_output(filepath)
        if df.empty:
            logger.warning("No annotated peaks found in %s — skipping summary.", filepath)
            return

        # Distance column
        logger.info("Summary params: distance_mode=%s, plot_type=%s, facet_by=%s, summary_input=%s",
                     distance_mode, plot_type, facet_by, summary_input)
        if distance_mode == 'directional':
            df = reconstruct_signed_distance(df)
            distance_col = 'signed_distance'
            x_label = '\u2190 upstream (5\')     Distance     downstream (3\') \u2192'
        else:
            distance_col = 'distance'
            x_label = 'Absolute Distance'
        logger.info("Using distance_col='%s', data range: [%s, %s]",
                     distance_col, df[distance_col].min(), df[distance_col].max())

        # Facet column mapping
        facet_map = {'feature': 'feature', 'anchor': 'feat_anchor', 'query': 'name'}
        facet_col = facet_map[facet_by]

        if facet_col not in df.columns:
            logger.error("Facet column '%s' not found in output. "
                         "Available columns: %s", facet_col, list(df.columns))
            return

        plot_fn = plot_distance_kde if plot_type == 'kde' else plot_distance_hist

        # Symmetric x-axis limits centered on 0 for directional mode
        if distance_mode == 'directional':
            max_abs = df[distance_col].abs().max()
            xlim = (-max_abs * 1.05, max_abs * 1.05)
            logger.info("Directional mode: max_abs=%s, xlim=%s", max_abs, xlim)
        else:
            xlim = None
            logger.info("Absolute mode: no xlim override")

        # Layout: allhits -> one subplot per query; finalhits -> single plot
        if summary_input == 'allhits':
            queries = sorted(df['name'].unique())
            n = len(queries)
            ncols = min(n, 2)
            nrows = (n + ncols - 1) // ncols
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(7 * ncols, 5 * nrows), squeeze=False)
            for idx, query in enumerate(queries):
                ax = axes[idx // ncols][idx % ncols]
                sub_df = df[df['name'] == query]
                plot_fn(sub_df, facet_col, distance_col, log_scale, ax, xlim=xlim)
                ax.set_title(query)
                ax.set_xlabel(x_label)
                if xlim is not None:
                    ax.set_xlim(xlim)
            for idx in range(n, nrows * ncols):
                axes[idx // ncols][idx % ncols].set_visible(False)
            fig.suptitle('Distance Distribution (all hits)', fontsize=14)
        else:
            fig, ax = plt.subplots(figsize=(8, 5))
            plot_fn(df, facet_col, distance_col, log_scale, ax, xlim=xlim)
            ax.set_xlabel(x_label)
            if xlim is not None:
                ax.set_xlim(xlim)
            ax.set_title('Distance Distribution (final hits)')

        # Log final axis limits for debugging
        if summary_input == 'allhits':
            for idx, query in enumerate(queries):
                ax = axes[idx // ncols][idx % ncols]
                logger.info("Subplot '%s' final xlim: %s", query, ax.get_xlim())
        else:
            logger.info("Final xlim: %s", ax.get_xlim())

        fig.tight_layout()
        output_path = output_prefix + '_summary.pdf'
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        logger.info("Summary plot saved to %s", output_path)

    except Exception as e:
        logger.error("Failed to generate summary plot: %s", str(e))
        logger.debug("Summary error details:", exc_info=True)
