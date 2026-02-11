# IMPLEMENTATION.md — Distance Distribution Plot

**Date:** 2026-02-11
**Status:** Phase 3 — Implementation Plan

---

## Summary of Decisions

| Topic | Decision |
|-------|----------|
| Signed distance | Reconstruct from output columns (no schema change) |
| Sign convention | Feature-centric (upstream = negative) |
| Data source | Parameter `source`, default `"finalhits"` |
| Plot type | KDE via Seaborn `kdeplot` |
| Grouping | Optional `split_by` parameter, default combined |
| CLI | No new flags — `--summary` triggers the plot with defaults |
| R summary | Comment out but keep the code |
| Dependencies | Add `pandas`, `matplotlib`, `seaborn` |

---

## Steps

### Step 1: Create `uropa/visualization.py` with `load_annotation_data()`

Create the new module with the data loading function.

```python
def load_annotation_data(filepath):
    """Read a UROPA output file (finalhits.txt or allhits.txt) into a DataFrame.
    Filters out unannotated peaks (feature == 'NA')."""
```

- Read TSV with `pandas.read_csv(filepath, sep="\t")`
- Filter out rows where `feature == "NA"`
- Validate that required columns exist: `peak_start`, `peak_end`, `feat_start`, `feat_end`, `feat_strand`, `feat_anchor`, `distance`
- Return DataFrame

**Checkpoint:** Function reads `test_data/` output files correctly.

---

### Step 2: Add `compute_signed_distance()`

Add the signed distance reconstruction function to `visualization.py`.

```python
def compute_signed_distance(df):
    """Reconstruct feature-centric signed distance from output columns.
    Returns a Series: negative = upstream, positive = downstream."""
```

Reconstruction logic (vectorized with pandas):

1. `peak_center = (peak_start + peak_end) / 2`
2. Compute `anchor_pos` based on `feat_anchor` and `feat_strand`:
   - `start` anchor: `feat_start` if strand != `-`, else `feat_end`
   - `center` anchor: `(feat_start + feat_end) / 2`
   - `end` anchor: `feat_end` if strand != `-`, else `feat_start`
3. `raw_distance = peak_center - anchor_pos`
4. Feature-centric flip: multiply by `-1` where `feat_strand == "-"`

**Checkpoint:** Unit-test with known peak/feature coordinates and strands.

---

### Step 3: Add `plot_distance_distribution()`

Add the main plotting function to `visualization.py`.

```python
def plot_distance_distribution(prefix, mode="relative", split_by=None, source="finalhits"):
    """Generate a KDE plot of peak-to-feature distances.

    Args:
        prefix: UROPA output prefix. Derives input from {prefix}_{source}.txt
                and output to {prefix}_distance_distribution.png
        mode: "relative" (absolute distance) or "stranded" (feature-centric signed distance)
        split_by: None (combined), "feature", or "query"
        source: "finalhits" or "allhits"
    """
```

Implementation:

1. Derive paths: `filepath = f"{prefix}_{source}.txt"`, `output = f"{prefix}_distance_distribution.png"`
2. Call `load_annotation_data(filepath)`
3. Choose distance data:
   - `mode == "relative"`: use the `distance` column directly (absolute values)
   - `mode == "stranded"`: call `compute_signed_distance(df)`
4. Determine `hue` parameter for seaborn:
   - `split_by is None`: no hue
   - `split_by == "feature"`: `hue="feature"`
   - `split_by == "query"`: `hue="name"`
5. Create figure with `matplotlib.pyplot`
6. Plot with `seaborn.kdeplot(data=df, x=distance_col, hue=hue_col, fill=True)`
7. Set labels:
   - X: `"Distance (bp)"` (relative) or `"← Upstream | Downstream → (bp)"` (stranded)
   - Y: `"Density"`
8. Add vertical line at x=0 for stranded mode
9. Save figure to output path
10. Close figure

**Checkpoint:** Generates a plot file from test data in both modes.

---

### Step 4: Update `setup.py` — add dependencies

Add `pandas`, `matplotlib`, and `seaborn` to `install_requires` in `setup.py:43-47`.

```python
install_requires=[
    'pysam',
    'psutil',
    'numpy',
    'pandas',
    'matplotlib',
    'seaborn'
],
```

**Checkpoint:** `pip install -e .` succeeds.

---

### Step 5: Integrate into `uropa.py` `--summary` block

Modify `uropa.py:574-590`:

1. Comment out the existing R summary call (lines 577-590), keeping the code intact
2. Add import and call to `plot_distance_distribution(output_prefix)`
3. Wrap in try/except with logger warning on failure (consistent with existing R fallback pattern)

```python
if args.summary:
    logger.info("Creating the Summary graphs of the results...")

    # R-based summary (commented out, kept for reference)
    # summary_script = "uropa_summary.R"
    # ...

    # Python distance distribution plot
    try:
        from uropa.visualization import plot_distance_distribution
        plot_distance_distribution(output_prefix)
        logger.info("Distance distribution plot created.")
    except Exception as e:
        logger.warning("Could not create distance distribution plot: %s", str(e))
```

**Checkpoint:** `uropa -b test.bed -g test.gtf --summary` produces the distance plot.

---

### Step 6: Test with test data

Run UROPA with `test_data/` to verify end-to-end:

1. Run annotation: `uropa -b test_data/genomic_regions.bed -g test_data/gencode.v29.annotation.chr19.gtf --summary`
2. Verify `*_distance_distribution.png` is created
3. Visually inspect the plot
4. Test standalone usage:
   ```python
   from uropa.visualization import plot_distance_distribution
   plot_distance_distribution("path/to/existing_output_prefix", mode="stranded")
   plot_distance_distribution("path/to/existing_output_prefix", split_by="feature")
   ```

---

## File Summary

| File | Action |
|------|--------|
| `uropa/visualization.py` | **Create** — Steps 1, 2, 3 |
| `setup.py` | **Edit** — Step 4 |
| `uropa/uropa.py` | **Edit** — Step 5 |

---

## Revision History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-11 | Initial implementation plan |
