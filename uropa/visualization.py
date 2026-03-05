"""
visualization.py: Enrichment/depletion analysis and plotting for UROPA.

Compares observed peak annotation proportions against genomic background
derived from the reference GTF. Uses binomial test with BH correction.
"""

import os
import csv
import logging
import numpy as np
from collections import defaultdict
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

logger = logging.getLogger("uropa")

# Biotype attribute names in order of detection priority
BIOTYPE_ATTRIBUTES = ["gene_type", "gene_biotype", "biotype"]


def _parse_gtf_attributes(attr_string):
    """Parse GTF column 9 attribute string into a dict."""
    attrs = {}
    for entry in attr_string.strip().split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(" ", 1)
        if len(parts) == 2:
            key = parts[0]
            value = parts[1].strip('"')
            attrs[key] = value
    return attrs


def _merge_intervals(intervals):
    """Merge overlapping intervals and return total non-redundant bp.

    Parameters
    ----------
    intervals : list of (int, int)
        List of (start, end) tuples, 1-based inclusive GTF coordinates.

    Returns
    -------
    int
        Total non-redundant base pairs covered.
    """
    if not intervals:
        return 0

    intervals.sort()
    merged_bp = 0
    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_bp += current_end - current_start + 1
            current_start, current_end = start, end

    merged_bp += current_end - current_start + 1
    return merged_bp


def compute_genomic_background(gtf_path, genome_size_mb=None):
    """Parse GTF to compute genomic coverage per feature type and biotype.

    Parameters
    ----------
    gtf_path : str
        Path to the full (unsubsetted) GTF file.
    genome_size_mb : float or None
        Total genome size in Mb. If None, inferred from max GTF coordinates.

    Returns
    -------
    dict
        {
            "feature_type": {name: {"coverage_bp": int, "proportion": float}},
            "biotype": {name: {"coverage_bp": int, "proportion": float}},
            "biotype_attribute": str or None,
            "genome_size_bp": int
        }
    """
    # Collect intervals per feature type and per biotype
    # Structure: {category: {chrom: [(start, end), ...]}}
    feat_intervals = defaultdict(lambda: defaultdict(list))
    biotype_intervals = defaultdict(lambda: defaultdict(list))
    chrom_max = defaultdict(int)

    # Detect biotype attribute from first gene-level feature
    biotype_attr = None

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            columns = line.rstrip().split("\t")
            if len(columns) < 9:
                continue

            chrom = columns[0]
            feature_type = columns[2]
            start = int(columns[3])
            end = int(columns[4])

            # Track max coordinate per chromosome
            if end > chrom_max[chrom]:
                chrom_max[chrom] = end

            # Collect feature type intervals
            feat_intervals[feature_type][chrom].append((start, end))

            # Collect biotype intervals (gene-level features only)
            if feature_type == "gene":
                attrs = _parse_gtf_attributes(columns[8])

                # Auto-detect biotype attribute on first gene
                if biotype_attr is None:
                    for candidate in BIOTYPE_ATTRIBUTES:
                        if candidate in attrs:
                            biotype_attr = candidate
                            break

                if biotype_attr and biotype_attr in attrs:
                    biotype_val = attrs[biotype_attr]
                    biotype_intervals[biotype_val][chrom].append((start, end))

    # Compute genome size
    if genome_size_mb is not None:
        genome_size_bp = int(genome_size_mb * 1e6)
    else:
        genome_size_bp = sum(chrom_max.values())

    if genome_size_bp == 0:
        logger.warning("Genome size is 0 — cannot compute enrichment background")
        return None

    # Merge intervals and compute coverage per feature type
    feat_background = {}
    for feat_type, chrom_dict in feat_intervals.items():
        total_bp = 0
        for chrom, intervals in chrom_dict.items():
            total_bp += _merge_intervals(intervals)
        feat_background[feat_type] = {
            "coverage_bp": total_bp,
            "proportion": total_bp / genome_size_bp
        }

    # Merge intervals and compute coverage per biotype
    bio_background = {}
    for biotype_val, chrom_dict in biotype_intervals.items():
        total_bp = 0
        for chrom, intervals in chrom_dict.items():
            total_bp += _merge_intervals(intervals)
        bio_background[biotype_val] = {
            "coverage_bp": total_bp,
            "proportion": total_bp / genome_size_bp
        }

    return {
        "feature_type": feat_background,
        "biotype": bio_background,
        "biotype_attribute": biotype_attr,
        "genome_size_bp": genome_size_bp
    }


def compute_observed_proportions(hits_path, biotype_attribute=None,
                                 deduplicate=False):
    """Parse a hits file to count peaks per feature type and biotype.

    Parameters
    ----------
    hits_path : str
        Path to a finalhits or allhits .txt output file.
    biotype_attribute : str or None
        Name of the biotype column to look for (e.g. "gene_type").
    deduplicate : bool
        If True, count each peak only once per category (for allhits where
        a peak can appear multiple times). Uses peak_id for deduplication.

    Returns
    -------
    dict
        {
            "feature_type": {name: count},
            "biotype": {name: count},
            "total_peaks": int,
            "query_feature_type": {query_name: {feature_type: count}},
            "query_biotype": {query_name: {biotype: count}},
        }
    """
    feat_counts = defaultdict(int)
    bio_counts = defaultdict(int)
    query_feat_counts = defaultdict(lambda: defaultdict(int))
    query_bio_counts = defaultdict(lambda: defaultdict(int))

    # For deduplication: track which (peak_id, category) pairs we've seen
    seen_feat = set()
    seen_bio = set()
    seen_query_feat = set()
    seen_query_bio = set()
    all_peak_ids = set()

    with open(hits_path) as f:
        header = f.readline().rstrip().split("\t")

        feat_idx = header.index("feature")
        name_idx = header.index("name")
        peak_id_idx = header.index("peak_id")
        bio_idx = None
        if biotype_attribute and biotype_attribute in header:
            bio_idx = header.index(biotype_attribute)

        for line in f:
            cols = line.rstrip().split("\t")
            peak_id = cols[peak_id_idx]
            all_peak_ids.add(peak_id)

            feature = cols[feat_idx]
            query_name = cols[name_idx]

            if feature == "NA":
                continue

            if deduplicate:
                key_feat = (peak_id, feature)
                if key_feat not in seen_feat:
                    seen_feat.add(key_feat)
                    feat_counts[feature] += 1

                key_qf = (peak_id, query_name, feature)
                if key_qf not in seen_query_feat:
                    seen_query_feat.add(key_qf)
                    query_feat_counts[query_name][feature] += 1

                if bio_idx is not None:
                    biotype_val = cols[bio_idx]
                    if biotype_val != "NA":
                        key_bio = (peak_id, biotype_val)
                        if key_bio not in seen_bio:
                            seen_bio.add(key_bio)
                            bio_counts[biotype_val] += 1

                        key_qb = (peak_id, query_name, biotype_val)
                        if key_qb not in seen_query_bio:
                            seen_query_bio.add(key_qb)
                            query_bio_counts[query_name][biotype_val] += 1
            else:
                feat_counts[feature] += 1
                query_feat_counts[query_name][feature] += 1

                if bio_idx is not None:
                    biotype_val = cols[bio_idx]
                    if biotype_val != "NA":
                        bio_counts[biotype_val] += 1
                        query_bio_counts[query_name][biotype_val] += 1

    total_peaks = len(all_peak_ids)

    return {
        "feature_type": dict(feat_counts),
        "biotype": dict(bio_counts),
        "total_peaks": total_peaks,
        "query_feature_type": {k: dict(v) for k, v in query_feat_counts.items()},
        "query_biotype": {k: dict(v) for k, v in query_bio_counts.items()},
    }


def enrichment_test(observed_counts, background, total_peaks):
    """Binomial test for enrichment/depletion with BH correction.

    Parameters
    ----------
    observed_counts : dict
        {category_name: peak_count}
    background : dict
        {category_name: {"proportion": float, "coverage_bp": int}}
    total_peaks : int
        Total number of peaks (including unannotated).

    Returns
    -------
    list of dict
        Each dict has: category, observed_count, observed_prop, expected_prop,
        fold_enrichment, log2_fc, p_value, p_adj, significance.
        Sorted by log2_fc descending.
    """
    results = []

    # Collect all categories present in either observed or background
    all_categories = set(observed_counts.keys()) | set(background.keys())

    for category in all_categories:
        obs_count = observed_counts.get(category, 0)
        obs_prop = obs_count / total_peaks if total_peaks > 0 else 0.0

        bg_info = background.get(category, {})
        exp_prop = bg_info.get("proportion", 0.0)

        # Skip categories with zero background proportion (can't compute enrichment)
        if exp_prop == 0.0:
            continue

        fold_enrichment = obs_prop / exp_prop
        log2_fc = np.log2(fold_enrichment) if fold_enrichment > 0 else -np.inf

        # Binomial test: is obs_count significantly different from expected?
        result = binomtest(obs_count, total_peaks, exp_prop, alternative='two-sided')
        p_value = result.pvalue

        results.append({
            "category": category,
            "observed_count": obs_count,
            "observed_prop": obs_prop,
            "expected_prop": exp_prop,
            "fold_enrichment": fold_enrichment,
            "log2_fc": log2_fc,
            "p_value": p_value,
        })

    # BH correction for multiple testing
    if results:
        p_values = [r["p_value"] for r in results]
        _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")
        for r, padj in zip(results, p_adj):
            r["p_adj"] = padj
            # Significance labels
            if padj < 0.001:
                r["significance"] = "***"
            elif padj < 0.01:
                r["significance"] = "**"
            elif padj < 0.05:
                r["significance"] = "*"
            else:
                r["significance"] = ""

    # Sort by log2_fc descending
    results.sort(key=lambda r: r["log2_fc"], reverse=True)
    return results


TSV_COLUMNS = [
    "category", "observed_count", "observed_prop", "expected_prop",
    "fold_enrichment", "log2_fc", "p_value", "p_adj", "significance"
]

TSV_METADATA_PREFIX = "#"


def write_enrichment_tsv(enrichment_results, tsv_path, total_peaks=None,
                         genome_size_bp=None):
    """Write enrichment results to a TSV file.

    The TSV contains all information needed for plotting. Metadata
    (total_peaks, genome_size_bp) is stored in comment lines at the top.

    Parameters
    ----------
    enrichment_results : list of dict
        Output from enrichment_test().
    tsv_path : str
        Output TSV file path.
    total_peaks : int or None
        Total number of peaks.
    genome_size_bp : int or None
        Genome size in bp.
    """
    with open(tsv_path, "w", newline="") as f:
        # Metadata as comment lines
        if total_peaks is not None:
            f.write(f"{TSV_METADATA_PREFIX}total_peaks={total_peaks}\n")
        if genome_size_bp is not None:
            f.write(f"{TSV_METADATA_PREFIX}genome_size_bp={genome_size_bp}\n")

        writer = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for r in enrichment_results:
            row = dict(r)
            # Format floats for readability; keep inf as string
            row["log2_fc"] = f"{r['log2_fc']:.6f}" if np.isfinite(r["log2_fc"]) else "-inf"
            row["observed_prop"] = f"{r['observed_prop']:.6f}"
            row["expected_prop"] = f"{r['expected_prop']:.6f}"
            row["fold_enrichment"] = f"{r['fold_enrichment']:.6f}"
            row["p_value"] = f"{r['p_value']:.6e}"
            row["p_adj"] = f"{r['p_adj']:.6e}"
            writer.writerow(row)

    logger.info(f"Enrichment TSV saved: {tsv_path}")


def read_enrichment_tsv(tsv_path):
    """Read enrichment results and metadata from a TSV file.

    Parameters
    ----------
    tsv_path : str
        Path to TSV file written by write_enrichment_tsv().

    Returns
    -------
    tuple (list of dict, dict)
        (enrichment_results, metadata) where metadata has keys like
        "total_peaks" and "genome_size_bp".
    """
    metadata = {}
    results = []

    with open(tsv_path) as f:
        # Read metadata comment lines
        lines = f.readlines()

    data_lines = []
    for line in lines:
        if line.startswith(TSV_METADATA_PREFIX) and "=" in line:
            kv = line[len(TSV_METADATA_PREFIX):].strip()
            key, val = kv.split("=", 1)
            metadata[key] = int(val) if val.isdigit() else val
        else:
            data_lines.append(line)

    reader = csv.DictReader(data_lines, delimiter="\t")
    for row in reader:
        results.append({
            "category": row["category"],
            "observed_count": int(row["observed_count"]),
            "observed_prop": float(row["observed_prop"]),
            "expected_prop": float(row["expected_prop"]),
            "fold_enrichment": float(row["fold_enrichment"]),
            "log2_fc": float(row["log2_fc"]) if row["log2_fc"] != "-inf" else -np.inf,
            "p_value": float(row["p_value"]),
            "p_adj": float(row["p_adj"]),
            "significance": row["significance"],
        })

    return results, metadata


def plot_enrichment_from_tsv(tsv_path, output_path, title, xlim_abs=None,
                            min_background_pct=0.0):
    """Create horizontal bar plot of log2(fold enrichment) from a TSV file.

    Parameters
    ----------
    tsv_path : str
        Path to enrichment TSV (written by write_enrichment_tsv).
    output_path : str
        Path to save PNG file.
    title : str
        Plot title.
    xlim_abs : float or None
        Symmetric x-axis limit. If None, computed from data.
    min_background_pct : float
        Minimum background proportion (in %) to include a category in the plot.
        E.g. 0.5 means categories covering < 0.5% of the genome are excluded.
        Default 0.0 (no filtering).
    """
    results, metadata = read_enrichment_tsv(tsv_path)
    if not results:
        return

    total_peaks = metadata.get("total_peaks")
    genome_size_bp = metadata.get("genome_size_bp")

    # Filter: remove categories with 0 observed peaks
    results = [r for r in results if r["observed_count"] > 0]

    # Filter: remove categories below minimum background proportion
    if min_background_pct > 0:
        min_prop = min_background_pct / 100.0
        results = [r for r in results if r["expected_prop"] >= min_prop]

    if not results:
        return

    # Extract data, capping -inf values
    finite_vals = [r["log2_fc"] for r in results if np.isfinite(r["log2_fc"])]
    if not finite_vals:
        return

    # Determine axis limit
    if xlim_abs is None:
        xlim_abs = max(abs(v) for v in finite_vals) * 1.15
    xlim_abs = max(xlim_abs, 0.5)  # minimum axis range

    # Color mapping: intensity reflects significance, gray for n.s.
    _COLOR_MAP = {
        # enriched (red gradient): darker = more significant
        ("enriched", "***"): "#922b21",
        ("enriched", "**"):  "#c0392b",
        ("enriched", "*"):   "#e74c3c",
        ("enriched", ""):    "#bdc3c7",
        # depleted (blue gradient): darker = more significant
        ("depleted", "***"): "#1a5276",
        ("depleted", "**"):  "#2980b9",
        ("depleted", "*"):   "#5dade2",
        ("depleted", ""):    "#bdc3c7",
    }

    # Prepare plot data
    categories = []
    log2_fcs = []
    colors = []
    significance = []
    for r in results:
        categories.append(r["category"])
        val = r["log2_fc"] if np.isfinite(r["log2_fc"]) else -xlim_abs
        log2_fcs.append(val)
        sig = r.get("significance", "")
        direction = "enriched" if val >= 0 else "depleted"
        colors.append(_COLOR_MAP[(direction, sig)])
        significance.append(sig)

    # Reverse for bottom-to-top ordering (highest enrichment at top)
    categories = categories[::-1]
    log2_fcs = log2_fcs[::-1]
    colors = colors[::-1]
    significance = significance[::-1]

    # Create figure
    n_bars = len(categories)
    fig_height = max(3, 0.4 * n_bars + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))

    y_pos = np.arange(n_bars)
    ax.barh(y_pos, log2_fcs, color=colors, edgecolor="none", height=0.7)

    # Significance markers
    for i, (val, sig) in enumerate(zip(log2_fcs, significance)):
        if sig:
            x_offset = 0.05 * xlim_abs if val >= 0 else -0.05 * xlim_abs
            ha = "left" if val >= 0 else "right"
            ax.text(val + x_offset, i, sig, ha=ha, va="center",
                    fontsize=9, fontweight="bold", color="#333333")

    # Reference line
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8)

    # Labels and formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories, fontsize=9)
    ax.set_xlim(-xlim_abs, xlim_abs)
    ax.set_xlabel("log$_2$(fold enrichment)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")

    # Annotation text
    annot_parts = []
    if total_peaks is not None:
        annot_parts.append(f"n = {total_peaks:,} peaks")
    if genome_size_bp is not None:
        annot_parts.append(f"genome = {genome_size_bp / 1e6:.1f} Mb")
    if annot_parts:
        annot_text = " | ".join(annot_parts)
        ax.text(0.99, 0.01, annot_text, transform=ax.transAxes, fontsize=8,
                ha="right", va="bottom", color="gray")

    # Significance legend
    ax.text(0.01, 0.01, "* p<0.05  ** p<0.01  *** p<0.001 (BH-adjusted)",
            transform=ax.transAxes, fontsize=7, ha="left", va="bottom", color="gray")

    sns.despine(ax=ax, left=True)
    ax.tick_params(axis="y", length=0)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Enrichment plot saved: {output_path}")


def plot_dotplot_from_tsv(tsv_path, output_path, title, xlim_abs=None,
                          min_background_pct=0.0):
    """Create dot plot of enrichment results from a TSV file.

    Dot size = observed count, color = adjusted p-value (log10 scale),
    x-axis = log2(fold enrichment). Shares the same xlim as barplots.

    Parameters
    ----------
    tsv_path : str
        Path to enrichment TSV (written by write_enrichment_tsv).
    output_path : str
        Path to save PNG file.
    title : str
        Plot title.
    xlim_abs : float or None
        Symmetric x-axis limit. If None, computed from data.
    min_background_pct : float
        Minimum background proportion (in %) to include a category in the plot.
    """
    results, metadata = read_enrichment_tsv(tsv_path)
    if not results:
        return

    total_peaks = metadata.get("total_peaks")
    genome_size_bp = metadata.get("genome_size_bp")

    # Filter: remove categories with 0 observed peaks
    results = [r for r in results if r["observed_count"] > 0]

    # Filter: remove categories below minimum background proportion
    if min_background_pct > 0:
        min_prop = min_background_pct / 100.0
        results = [r for r in results if r["expected_prop"] >= min_prop]

    if not results:
        return

    finite_vals = [r["log2_fc"] for r in results if np.isfinite(r["log2_fc"])]
    if not finite_vals:
        return

    # Determine axis limit
    if xlim_abs is None:
        xlim_abs = max(abs(v) for v in finite_vals) * 1.15
    xlim_abs = max(xlim_abs, 0.5)

    # Prepare data
    categories = []
    log2_fcs = []
    counts = []
    padj_vals = []
    for r in results:
        categories.append(r["category"])
        val = r["log2_fc"] if np.isfinite(r["log2_fc"]) else -xlim_abs
        log2_fcs.append(val)
        counts.append(r["observed_count"])
        padj_vals.append(r["p_adj"])

    # Reverse for bottom-to-top ordering (highest enrichment at top)
    categories = categories[::-1]
    log2_fcs = log2_fcs[::-1]
    counts = counts[::-1]
    padj_vals = padj_vals[::-1]

    log2_fcs = np.array(log2_fcs)
    counts = np.array(counts, dtype=float)
    padj_vals = np.array(padj_vals)

    # Size scaling: map observed counts to dot area
    min_size, max_size = 30, 400
    if counts.max() == counts.min():
        sizes = np.full_like(counts, (min_size + max_size) / 2)
    else:
        sizes = min_size + (counts - counts.min()) / (counts.max() - counts.min()) * (max_size - min_size)

    # Color: -log10(p_adj), capped for display
    neg_log10_p = -np.log10(np.clip(padj_vals, 1e-300, 1.0))
    cmap = plt.cm.Reds
    vmin = 0
    vmax = max(neg_log10_p.max(), 1.3)  # at least -log10(0.05)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Create figure
    n_dots = len(categories)
    fig_height = max(3, 0.4 * n_dots + 1.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))

    y_pos = np.arange(n_dots)
    scatter = ax.scatter(log2_fcs, y_pos, s=sizes, c=neg_log10_p,
                         cmap=cmap, norm=norm, edgecolors="#333333",
                         linewidths=0.5, zorder=3)

    # Reference line
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, zorder=1)

    # Labels and formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories, fontsize=9)
    ax.set_xlim(-xlim_abs, xlim_abs)
    ax.set_xlabel("log$_2$(fold enrichment)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02, aspect=30)
    cbar.set_label("$-\\log_{10}$(p$_{adj}$)", fontsize=10)

    # Size legend
    unique_counts = sorted(set(counts))
    if len(unique_counts) > 3:
        legend_counts = [int(counts.min()), int(np.median(counts)), int(counts.max())]
    else:
        legend_counts = [int(c) for c in unique_counts]
    legend_counts = sorted(set(legend_counts))
    legend_handles = []
    for c in legend_counts:
        if counts.max() == counts.min():
            s = (min_size + max_size) / 2
        else:
            s = min_size + (c - counts.min()) / (counts.max() - counts.min()) * (max_size - min_size)
        legend_handles.append(
            ax.scatter([], [], s=s, c="gray", edgecolors="#333333",
                       linewidths=0.5, label=str(c))
        )
    ax.legend(handles=legend_handles, title="Observed", loc="lower left",
              fontsize=8, title_fontsize=9, framealpha=0.9, labelspacing=1.2)

    # Annotation text
    annot_parts = []
    if total_peaks is not None:
        annot_parts.append(f"n = {total_peaks:,} peaks")
    if genome_size_bp is not None:
        annot_parts.append(f"genome = {genome_size_bp / 1e6:.1f} Mb")
    if annot_parts:
        annot_text = " | ".join(annot_parts)
        ax.text(0.99, 0.01, annot_text, transform=ax.transAxes, fontsize=8,
                ha="right", va="bottom", color="gray")

    sns.despine(ax=ax, left=True)
    ax.tick_params(axis="y", length=0)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Dot plot saved: {output_path}")


def _compute_shared_xlim(all_results):
    """Compute a shared symmetric x-axis limit across multiple enrichment results.

    Parameters
    ----------
    all_results : list of list of dict
        List of enrichment_test() outputs.

    Returns
    -------
    float
        Symmetric x-axis limit (positive value).
    """
    max_abs = 0.5  # minimum
    for results in all_results:
        for r in results:
            if np.isfinite(r["log2_fc"]):
                max_abs = max(max_abs, abs(r["log2_fc"]))
    return max_abs * 1.15


def _run_enrichment_for_hits(observed, background, output_dir, title_suffix,
                            min_background_pct=0.0):
    """Run enrichment tests, write TSVs, and generate plots for one hits file.

    Parameters
    ----------
    observed : dict
        Output from compute_observed_proportions().
    background : dict
        Output from compute_genomic_background().
    output_dir : str
        Directory for output files.
    title_suffix : str
        Suffix appended to plot titles (e.g. "(finalhits)" or "(allhits)").
    min_background_pct : float
        Minimum background proportion (%) to include a category in plots.

    Returns
    -------
    list of list of dict
        All enrichment results (for shared xlim computation).
    """
    total_peaks = observed["total_peaks"]
    genome_size_bp = background["genome_size_bp"]
    tsv_metadata = dict(total_peaks=total_peaks, genome_size_bp=genome_size_bp)

    os.makedirs(output_dir, exist_ok=True)

    # Collect (tsv_path, title, results) tuples
    tsv_entries = []

    # Global feature type
    results_ft = enrichment_test(
        observed["feature_type"], background["feature_type"], total_peaks
    )
    tsv_ft = os.path.join(output_dir, "enrichment_feature_type.tsv")
    write_enrichment_tsv(results_ft, tsv_ft, **tsv_metadata)
    tsv_entries.append((
        tsv_ft, f"Feature Type Enrichment {title_suffix}", results_ft
    ))

    # Global biotype
    if background["biotype"]:
        results_bio = enrichment_test(
            observed["biotype"], background["biotype"], total_peaks
        )
        tsv_bio = os.path.join(output_dir, "enrichment_biotype.tsv")
        write_enrichment_tsv(results_bio, tsv_bio, **tsv_metadata)
        tsv_entries.append((
            tsv_bio,
            f"Biotype Enrichment ({background['biotype_attribute']}) {title_suffix}",
            results_bio
        ))

    # Per-query results
    for query_name, q_feat_counts in observed["query_feature_type"].items():
        safe_name = "".join(
            c if c.isalnum() or c in "-_" else "_" for c in query_name
        )

        q_ft = enrichment_test(
            q_feat_counts, background["feature_type"], total_peaks
        )
        tsv_q_ft = os.path.join(output_dir, f"enrichment_{safe_name}_feature_type.tsv")
        write_enrichment_tsv(q_ft, tsv_q_ft, **tsv_metadata)
        tsv_entries.append((
            tsv_q_ft,
            f"Feature Type Enrichment — {query_name} {title_suffix}",
            q_ft
        ))

        if background["biotype"] and query_name in observed["query_biotype"]:
            q_bio = enrichment_test(
                observed["query_biotype"][query_name],
                background["biotype"],
                total_peaks
            )
            tsv_q_bio = os.path.join(output_dir, f"enrichment_{safe_name}_biotype.tsv")
            write_enrichment_tsv(q_bio, tsv_q_bio, **tsv_metadata)
            tsv_entries.append((
                tsv_q_bio,
                f"Biotype Enrichment — {query_name} {title_suffix}",
                q_bio
            ))

    return tsv_entries


def run_enrichment_analysis(gtf_path, finalhits_path, output_prefix,
                            genome_size_mb=None, min_background_pct=0.0):
    """Run full enrichment/depletion analysis and generate plots.

    Pipeline: compute background + observed -> statistical tests -> write TSV
    files -> read TSVs back and generate PNG plots. The TSV files are the
    canonical intermediate output and contain all data needed for plotting.

    Produces two output folders:
    - {output_prefix}_enrichment_finalhits/ — based on best-hit-per-peak
    - {output_prefix}_enrichment_allhits/  — based on all valid annotations
      (deduplicated: each peak counted once per category)

    Parameters
    ----------
    gtf_path : str
        Path to the full (unsubsetted) GTF file.
    finalhits_path : str
        Path to the finalhits.txt output file.
    output_prefix : str
        Prefix for output directories.
    genome_size_mb : float or None
        Total genome size in Mb. If None, inferred from GTF.
    min_background_pct : float
        Minimum background proportion (%) to include a category in plots.
        Default 0.0 (no filtering). TSV files always contain all categories.
    """
    logger.info("Computing genomic background from GTF for enrichment analysis")
    background = compute_genomic_background(gtf_path, genome_size_mb)
    if background is None:
        logger.warning("Skipping enrichment analysis — could not compute background")
        return

    # Derive allhits path from finalhits path
    allhits_path = finalhits_path.replace("_finalhits.txt", "_allhits.txt")

    # Compute observed proportions for both hit types
    observed_final = compute_observed_proportions(
        finalhits_path, biotype_attribute=background["biotype_attribute"],
        deduplicate=False
    )
    if observed_final["total_peaks"] == 0:
        logger.warning("Skipping enrichment analysis — no peaks found")
        return

    all_tsv_entries = []

    # Finalhits enrichment
    finalhits_dir = output_prefix + "_enrichment_finalhits"
    entries_final = _run_enrichment_for_hits(
        observed_final, background, finalhits_dir, "(finalhits)",
        min_background_pct
    )
    all_tsv_entries.extend(entries_final)

    # Allhits enrichment (if file exists)
    if os.path.isfile(allhits_path):
        observed_all = compute_observed_proportions(
            allhits_path, biotype_attribute=background["biotype_attribute"],
            deduplicate=True
        )
        allhits_dir = output_prefix + "_enrichment_allhits"
        entries_all = _run_enrichment_for_hits(
            observed_all, background, allhits_dir, "(allhits)",
            min_background_pct
        )
        all_tsv_entries.extend(entries_all)
    else:
        logger.warning("Allhits file not found: %s — skipping allhits enrichment",
                        allhits_path)

    # Shared xlim across ALL plots (both finalhits and allhits)
    all_results = [entry[2] for entry in all_tsv_entries]
    xlim_abs = _compute_shared_xlim(all_results)

    # Generate plots from TSV files (barplot + dotplot for each)
    for tsv_path, title, _ in all_tsv_entries:
        png_path = os.path.splitext(tsv_path)[0] + ".png"
        plot_enrichment_from_tsv(tsv_path, png_path, title, xlim_abs=xlim_abs,
                                min_background_pct=min_background_pct)
        dot_path = os.path.splitext(tsv_path)[0] + "_dotplot.png"
        plot_dotplot_from_tsv(tsv_path, dot_path, title, xlim_abs=xlim_abs,
                              min_background_pct=min_background_pct)

    logger.info("Enrichment analysis complete")
