# ARCHITECTURE.md — Distance Distribution Plot

**Date:** 2026-02-11
**Status:** Phase 2 — Architecture (Draft)

---

## Design Decisions (from Analysis)

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Signed distance source | Reconstruct from output columns | No schema change to output files |
| Sign convention | Feature-centric (upstream = negative) | Biologically meaningful |
| Data source | Parameter, default `finalhits.txt` | Flexible, sensible default |
| Overlapping peaks | Use distance as-is | UROPA's feat_anchor already accounts for this |
| Plot type | KDE (kernel density estimate) | Smooth distribution visualization |
| Grouping | Optional split by feature/query, default combined | Flexibility without clutter |
| Integration | Function in package + called from `--summary` | Reusable and integrated |

---

## Module Boundaries

### New file: `uropa/visualization.py`

A single module (not a subpackage) containing the distance distribution plot function(s). This keeps things simple and can be refactored into a subpackage later if more plots are added.

### Modified files

| File | Change |
|------|--------|
| `uropa/uropa.py` | Call visualization from `--summary` block (lines 574-590) |
| `setup.py` | Add `matplotlib` and `pandas` as dependencies |

---

## Data Flow

```
                    Output files (finalhits.txt / allhits.txt)
                                    │
                                    ▼
                    ┌──────────────────────────────┐
                    │  load_annotation_data()       │
                    │  Read TSV into pandas DataFrame│
                    └──────────────┬───────────────┘
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │  compute_signed_distance()    │
                    │  Reconstruct feature-centric  │
                    │  signed distance from:        │
                    │  - peak_start, peak_end       │
                    │  - feat_start, feat_end       │
                    │  - feat_strand, feat_anchor   │
                    └──────────────┬───────────────┘
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │  plot_distance_distribution() │
                    │  Parameters:                  │
                    │  - mode: "relative"|"stranded"│
                    │  - split_by: None|"feature"|  │
                    │             "query"            │
                    │  KDE plot via matplotlib       │
                    └──────────────┬───────────────┘
                                   │
                                   ▼
                              PNG / PDF file
```

---

## Key Functions / Contracts

### `load_annotation_data(filepath: str) -> pd.DataFrame`

Reads a UROPA output file (finalhits.txt or allhits.txt) into a DataFrame.

- **Input:** Path to TSV file with UROPA header
- **Output:** DataFrame with columns as per UROPA output schema
- **Filters out:** Rows where `feature == "NA"` (unannotated peaks)

### `compute_signed_distance(df: pd.DataFrame) -> pd.Series`

Reconstructs the feature-centric signed distance from existing columns.

**Reconstruction logic:**

```
1. peak_center = (peak_start + peak_end) / 2

2. anchor_pos = based on feat_anchor and feat_strand:
   - feat_anchor == "start" → feat_start if feat_strand != "-" else feat_end
   - feat_anchor == "center" → (feat_start + feat_end) / 2
   - feat_anchor == "end"   → feat_end if feat_strand != "-" else feat_start

3. raw_distance = peak_center - anchor_pos

4. Feature-centric sign:
   - feat_strand == "+" or ".": signed_distance = raw_distance
     (negative = upstream, positive = downstream)
   - feat_strand == "-":       signed_distance = -raw_distance
     (flip sign so upstream is still negative)
```

- **Input:** DataFrame with columns `peak_start`, `peak_end`, `feat_start`, `feat_end`, `feat_strand`, `feat_anchor`
- **Output:** Series of signed distances (int/float)

### `plot_distance_distribution(prefix, mode="relative", split_by=None, source="finalhits")`

Main entry point for generating the distance distribution plot.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `prefix` | str | — | UROPA output prefix (same as used by the run). Used to derive both input and output paths. |
| `mode` | str | `"relative"` | `"relative"` = absolute distance; `"stranded"` = feature-centric signed distance |
| `split_by` | str or None | `None` | `None` = combined; `"feature"` = split by feature type; `"query"` = split by query name |
| `source` | str | `"finalhits"` | Which output file to read: `"finalhits"` or `"allhits"` |

**Path derivation from prefix:**

- **Input file:** `{prefix}_{source}.txt` (e.g., `results_finalhits.txt`)
- **Output plot:** `{prefix}_distance_distribution.png`

**Behavior by mode:**

- **`"relative"`:** X-axis shows absolute distance (≥ 0). Single-sided KDE.
- **`"stranded"`:** X-axis shows signed distance (upstream ← 0 → downstream). Two-sided KDE centered on 0.

**Behavior by split_by:**

- **`None`:** Single KDE curve on one axes.
- **`"feature"`:** One KDE curve per unique value in the `feature` column, overlaid with legend.
- **`"query"`:** One KDE curve per unique value in the `name` column, overlaid with legend.

**Plot details:**

- Library: `matplotlib` (with optional `scipy.stats.gaussian_kde` for KDE)
- X-axis label: `"Distance (bp)"` (relative) or `"← Upstream | Downstream → (bp)"` (stranded)
- Y-axis label: `"Density"`
- Title: `"Distance to Feature"` (adjustable)
- Vertical line at x=0 for stranded mode
- Legend when split_by is used
- Output format inferred from file extension (.png, .pdf)

---

## Integration into `--summary`

Current flow in `uropa.py:574-590`:

```python
if args.summary:
    # ... calls uropa_summary.R ...
```

**Proposed change:** Add Python distance plot generation alongside or as replacement:

```python
if args.summary:
    # Existing R summary (keep for now, can deprecate later)
    # ...

    # New: Python distance distribution plot
    from uropa.visualization import plot_distance_distribution
    plot_distance_distribution(output_prefix)
    # All other parameters use defaults:
    # mode="relative", split_by=None, source="finalhits"
```

**Standalone usage** (outside CLI, e.g. in a notebook or script):

```python
from uropa.visualization import plot_distance_distribution

# Minimal — just prefix
plot_distance_distribution("/path/to/results/my_experiment")

# Customized
plot_distance_distribution(
    "/path/to/results/my_experiment",
    mode="stranded",
    split_by="feature",
    source="allhits"
)
```

Only `filepath` needs to be passed — `output`, `mode`, and `split_by` all use sensible defaults. This keeps the `--summary` integration minimal.

For now, the `mode` and `split_by` are hardcoded when called from `--summary`. Future work could expose these as CLI parameters (e.g., `--distance-mode stranded`).

---

## Dependencies

| Package | Purpose | New? |
|---------|---------|------|
| `pandas` | Data loading and manipulation | Yes |
| `matplotlib` | Base plotting backend | Yes |
| `seaborn` | KDE plotting via `kdeplot` | Yes |

**Note:** `pandas`, `matplotlib`, and `seaborn` are standard scientific Python packages with no heavy installation burden. Seaborn provides `kdeplot` with built-in support for `hue`-based grouping, which maps directly to the `split_by` parameter.

---

## Testing Strategy

1. **Unit tests for `compute_signed_distance()`:**
   - + strand feature, peak upstream → negative distance
   - + strand feature, peak downstream → positive distance
   - − strand feature, peak upstream → negative distance
   - − strand feature, peak downstream → positive distance
   - Overlapping peak/feature → correct small distance
   - Each feat_anchor type (start, center, end)

2. **Integration test for `plot_distance_distribution()`:**
   - Reads test_data output files
   - Generates plot without error in both modes
   - Output file is created and non-empty

3. **Visual inspection:**
   - Generate plots from `test_data/` and review manually

---

## Open Questions for Discussion

1. **Should `mode` and `split_by` be exposed as CLI arguments now?** Current proposal hardcodes defaults when called from `--summary`. We could add `--distance-mode` and `--distance-split` flags.

2. **Keep R summary alongside Python plot?** The proposal adds the Python plot in addition to the R summary. Should the R call be left as-is, or should we begin replacing it?

3. **scipy dependency:** Should we use `scipy.stats.gaussian_kde` for more control, or rely on matplotlib's built-in KDE capabilities (e.g., via `matplotlib`'s `Axes.violinplot` internals or manual KDE)? Matplotlib alone can do basic KDE via `Axes.hist(..., density=True)` + smoothing, but `scipy` gives cleaner KDE out of the box.

---

## Revision History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-11 | Initial architecture draft |
