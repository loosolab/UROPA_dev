import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


REQUIRED_COLUMNS = [
    "peak_start", "peak_end",
    "feat_start", "feat_end", "feat_strand", "feat_anchor",
    "distance"
]


def load_annotation_data(filepath):
    """Read a UROPA output file (finalhits.txt or allhits.txt) into a DataFrame.

    Filters out unannotated peaks (rows where feature == 'NA').

    Args:
        filepath: Path to a UROPA output TSV file.

    Returns:
        pandas.DataFrame with annotated peaks only.

    Raises:
        FileNotFoundError: If filepath does not exist.
        ValueError: If required columns are missing.
    """
    df = pd.read_csv(filepath, sep="\t", keep_default_na=False, na_values=[])

    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {filepath}: {', '.join(missing)}"
        )

    df = df[df["feature"] != "NA"].copy()

    return df


def compute_signed_distance(df):
    """Reconstruct feature-centric signed distance from output columns.

    Convention: negative = upstream of feature, positive = downstream.
    This is independent of genomic strand — the sign always reflects
    the feature's orientation.

    Args:
        df: DataFrame from load_annotation_data() with columns peak_start,
            peak_end, feat_start, feat_end, feat_strand, feat_anchor.

    Returns:
        pandas.Series of signed distances (float).
    """
    peak_start = df["peak_start"].astype(float)
    peak_end = df["peak_end"].astype(float)
    feat_start = df["feat_start"].astype(float)
    feat_end = df["feat_end"].astype(float)

    peak_center = (peak_start + peak_end) / 2

    # Compute anchor position (strand-aware, matching annotation.py logic)
    is_minus = df["feat_strand"] == "-"
    is_start = df["feat_anchor"] == "start"
    is_center = df["feat_anchor"] == "center"
    is_end = df["feat_anchor"] == "end"

    anchor_pos = np.where(
        is_start,
        np.where(is_minus, feat_end, feat_start),
        np.where(
            is_center,
            (feat_start + feat_end) / 2,
            np.where(is_minus, feat_start, feat_end)  # is_end
        )
    )

    raw_distance = peak_center - anchor_pos

    # Feature-centric: flip sign for minus strand so upstream is always negative
    signed_distance = np.where(is_minus, -raw_distance, raw_distance)

    return pd.Series(signed_distance, index=df.index, dtype=float)


def plot_distance_distribution(prefix, mode="relative", split_by=None, source="finalhits"):
    """Generate a KDE plot of peak-to-feature distances.

    Args:
        prefix: UROPA output prefix. Derives input from {prefix}_{source}.txt
                and output to {prefix}_distance_distribution.png.
        mode: "relative" (absolute distance) or "stranded" (feature-centric
              signed distance where upstream is negative).
        split_by: None (combined KDE), "feature" (split by feature type),
                  or "query" (split by query name).
        source: "finalhits" or "allhits".
    """
    filepath = f"{prefix}_{source}.txt"
    output = f"{prefix}_distance_distribution.png"

    df = load_annotation_data(filepath)

    # Prepare distance data based on mode
    if mode == "stranded":
        df["_distance"] = compute_signed_distance(df)
    else:
        df["_distance"] = df["distance"].astype(float)

    # Determine hue column for splitting
    hue_col = None
    if split_by == "feature":
        hue_col = "feature"
    elif split_by == "query":
        hue_col = "name"

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 5))

    sns.kdeplot(
        data=df,
        x="_distance",
        hue=hue_col,
        fill=True,
        alpha=0.3,
        ax=ax,
    )

    # Labels and styling
    if mode == "stranded":
        ax.set_xlabel("\u2190 Upstream | Downstream \u2192 (bp)")
        ax.axvline(x=0, color="grey", linestyle="--", linewidth=0.8)
    else:
        ax.set_xlabel("Distance (bp)")

    ax.set_ylabel("Density")
    ax.set_title("Distance to Feature")

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)
