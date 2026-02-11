# ANALYSIS.md — Distance Distribution Plot

**Date:** 2026-02-11
**Status:** Phase 1 — Analysis

---

## Problem Restatement

Implement a new Python-based plot that visualizes the **distance of annotated features to peaks**. The plot should:

1. Show a **distribution** of distances (e.g., histogram, KDE, or combined)
2. Support a **parameter** to toggle between two modes:
   - **Relative distance** — absolute distance without strand consideration
   - **Stranded distance** — signed distance accounting for feature strand (upstream = negative, downstream = positive)
3. Be implemented in **Python** (no R dependency)

---

## Relevant Files / Modules

| File | Relevance |
|------|-----------|
| `uropa/annotation.py:79-98` | `distance_to_peak_center()` — calculates both `raw_distance` (signed) and `distance` (absolute) |
| `uropa/annotation.py:118-158` | `get_relative_location()` — upstream/downstream classification, strand-aware |
| `uropa/annotation.py:62-69` | `anchor_pos` — strand-aware anchor positions (start/end flipped for `-` strand) |
| `uropa/annotation.py:249-258` | Distance validation — strand-aware upstream/downstream distance checks |
| `uropa/uropa.py:420-427` | Output column schema — `distance` column is absolute only; `raw_distance` is **not** written to output files |
| `uropa/uropa.py:575-590` | Current summary trigger — calls R script |
| `utils/uropa_summary.R:233-248` | Current R distance density plot |
| `utils/uropa_summary.R:307-339` | Current R distance histograms per query |
| `VISUALIZATION_ROADMAP.md:168-194` | Roadmap section 1.2 — Enhanced Distance Distribution Plot |

---

## Key Data Available

### In annotation dict (in-memory, during annotation)

| Field | Type | Description |
|-------|------|-------------|
| `raw_distance` | int | **Signed** distance: peak_center − anchor_pos. Negative = peak is upstream of feature on + strand |
| `distance` | int | **Absolute** distance (always ≥ 0) |
| `feat_strand` | str | `+` or `-` — needed to interpret raw_distance direction |
| `feat_anchor` | str | Which anchor was used: `start`, `center`, or `end` |
| `relative_location` | str | Upstream / Downstream / PeakInsideFeature / etc. |

### In output files (finalhits.txt / allhits.txt)

| Column | Available |
|--------|-----------|
| `distance` | Yes — absolute distance |
| `raw_distance` | **No** — not currently written to output |
| `feat_strand` | Yes |
| `feat_anchor` | Yes |
| `relative_location` | Yes |
| `feature` | Yes — for grouping by feature type |
| `name` (query) | Yes — for grouping by query |

---

## Risks / Unknowns / Assumptions

### Risk 1: `raw_distance` not in output files
The signed distance (`raw_distance`) is calculated during annotation but **not written** to the output TSV files. Only the absolute `distance` is stored.

**Options:**
- **(a)** Add `raw_distance` as a new output column — requires changing output schema (breaking change for downstream tools)
- **(b)** Reconstruct signed distance from existing columns (`distance`, `relative_location`, `feat_strand`) — no schema change, but reconstruction may be imperfect for overlapping peaks
- **(c)** Generate the plot in-memory during the annotation run (before data is written) — no schema change, but plot can't be generated post-hoc from output files

**Assumption:** This needs a design decision. Option (a) is cleanest but most disruptive.

### Risk 2: Stranded distance semantics
For the "stranded" mode, we need to define the sign convention clearly:
- **Convention A (feature-centric):** Negative = upstream of feature, Positive = downstream → requires strand-flipping for `-` strand features
- **Convention B (genomic):** Negative = peak is left of anchor, Positive = peak is right → no strand flipping

The existing `raw_distance` uses genomic coordinates (peak_center − anchor_pos), which is Convention B. For biological interpretation, Convention A (feature-centric) is more meaningful.

**Assumption:** We should use Convention A (feature-centric), where upstream/downstream is relative to the feature's orientation.

### Risk 3: Which data source?
- `finalhits.txt` — best annotation per peak (cleaner, less bias)
- `allhits.txt` — all valid annotations (more data points, but biased toward peaks with many annotations)

**Assumption:** The plot should support both, but default to `finalhits.txt`.

### Risk 4: Visualization package structure
The `VISUALIZATION_ROADMAP.md` proposes a `uropa/visualization/` subpackage. No such package exists yet. This is the first plot to be implemented.

**Assumption:** We need to decide whether to:
- Create the full package structure now (future-proof but more initial work)
- Start with a single module (simpler start, refactor later)

### Risk 5: Overlapping features (distance = 0)
When peaks fully overlap features (or vice versa), the distance may be 0 and the directional interpretation is ambiguous. These cases correspond to `PeakInsideFeature`, `FeatureInsidePeak`, `OverlapStart`, `OverlapEnd`.

**Assumption:** These should still be included in the distribution at their calculated distance values.

---

## Proposed Approach Options

### Option 1: Post-hoc plot from output files (add `raw_distance` column)

- Add `raw_distance` to the output TSV schema
- Create a Python plotting module that reads finalhits/allhits files
- Plot can be generated independently of the annotation run

**Pros:** Decoupled from annotation, reproducible, shareable
**Cons:** Schema change, breaks backwards compatibility

### Option 2: Post-hoc plot with reconstructed signed distance

- Keep output schema unchanged
- Reconstruct signed distance from `distance` + `relative_location` + `feat_strand`
- Create a Python plotting module that reads finalhits/allhits files

**Pros:** No schema change, post-hoc generation possible
**Cons:** Reconstruction is imprecise for overlapping peaks (distance=0 cases), adds logic complexity

### Option 3: In-memory plot during annotation + optional `raw_distance` column

- Generate the plot during the annotation run from in-memory data
- Optionally add `raw_distance` to output for post-hoc regeneration
- New `--distance-plot` flag or integrate into `--summary`

**Pros:** Clean access to all data, can also support post-hoc
**Cons:** More integration work with the annotation pipeline

---

## Open Design Questions

1. **Should `raw_distance` be added to the output files?** (enables post-hoc plotting)
2. **Stranded distance convention:** feature-centric (upstream=negative) or genomic (left=negative)?
3. **Default data source:** finalhits only, allhits only, or user-selectable?
4. **Package structure:** Full `visualization/` subpackage or single module for now?
5. **Plot type preference:** Histogram, KDE, violin, ridge plot, or combination?
6. **Grouping:** Should the plot support faceting by feature type / query name?
7. **Integration:** Standalone script, part of `--summary`, or new CLI flag?

---

## Revision History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-11 | Initial analysis |
