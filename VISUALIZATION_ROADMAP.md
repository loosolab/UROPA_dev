# UROPA Visualization Roadmap

**Document Version:** 1.0
**Date:** 2026-01-22
**Status:** Draft for Revision

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current Architecture Overview](#current-architecture-overview)
3. [Available Data for Visualization](#available-data-for-visualization)
4. [Current Visualization Capabilities](#current-visualization-capabilities)
5. [Proposed Future Visualizations](#proposed-future-visualizations)
6. [Integration Architecture](#integration-architecture)
7. [Implementation Priorities](#implementation-priorities)

---

## Executive Summary

UROPA (Universal RObust Peak Annotator) is a mature tool for annotating genomic peaks to features from GTF files. The current version (4.0.3) includes R-based summary visualizations, but there is significant opportunity to expand the visualization capabilities to provide deeper insights into annotation results.

This document proposes a comprehensive set of visualizations organized into three tiers:
- **Core Statistics** - Essential metrics for every annotation run
- **Comparative Analysis** - Multi-query and multi-sample comparisons
- **Advanced Analytics** - Feature enrichment, quality metrics, and genomic context

---

## Current Architecture Overview

### Core Data Flow

```
BED (peaks) + GTF (features) + Config (queries)
                    ↓
            [Annotation Engine]
            (Python/multiprocessing)
                    ↓
    ┌───────────────┴───────────────┐
    ↓                               ↓
allhits.txt                   finalhits.txt
(all candidates)              (best per peak)
    ↓                               ↓
    └───────────────┬───────────────┘
                    ↓
            [uropa_summary.R]
                    ↓
              summary.pdf
```

### Key Modules

| Module | Location | Responsibility |
|--------|----------|----------------|
| `uropa.py` | `uropa/uropa.py` | Main orchestrator, CLI, multiprocessing |
| `annotation.py` | `uropa/annotation.py` | Core annotation logic, distance/overlap calculations |
| `utils.py` | `uropa/utils.py` | File I/O, validation, logging |
| `uropa_summary.R` | `utils/uropa_summary.R` | Current visualization (R-based) |

### Existing Experimental Work

A `pyviz` branch exists with initial Python visualization code (`visualization.py`) that provides:
- Distribution plots (histogram, boxplot, violin)
- Count plots (pie, bar)
- Peak count plots (stacked bar)
- UpSet plot skeleton

This work can serve as a foundation for the proposed visualizations.

---

## Available Data for Visualization

### Per-Annotation Metrics (from annotation dict)

| Field | Type | Description | Visualization Potential |
|-------|------|-------------|------------------------|
| `distance` | int | Absolute distance to feature anchor | Histograms, density plots, boxplots |
| `raw_distance` | int | Signed distance (upstream negative) | Directional distance plots |
| `relative_location` | str | Spatial relationship (6 categories) | Pie charts, stacked bars |
| `feat_ovl_peak` | float | % of feature overlapped by peak | Overlap distribution plots |
| `peak_ovl_feat` | float | % of peak overlapped by feature | Overlap distribution plots |
| `feature` | str | Feature type (gene, exon, etc.) | Category distributions |
| `query_name` | str | Query identifier | Grouping variable |
| `best_hit` | int | Whether this is the best annotation | Filtering |
| `peak_score` | float | Peak score/p-value | Score-based filtering, correlation |
| `feat_strand` / `peak_strand` | str | Strand information | Strand concordance analysis |

### Aggregatable Statistics

| Statistic | Derivation | Use Case |
|-----------|------------|----------|
| Annotation rate | annotated peaks / total peaks | QC metric |
| Distance distribution | histogram of distances | Feature proximity analysis |
| Feature type distribution | counts per feature | Genomic context |
| Query hit distribution | peaks per query | Query effectiveness |
| Overlap statistics | mean/median overlap % | Binding site characterization |
| Strand concordance | same vs opposite strand hits | Regulatory direction |

### Multi-Run Comparisons

When multiple annotation runs exist (different samples, conditions, or parameters):
- Peak overlap between samples
- Annotation consistency
- Feature distribution shifts

---

## Current Visualization Capabilities

### Existing Plots (uropa_summary.R)

1. **Cover Page** - Configuration summary, run parameters
2. **Distance Density Plot** - Overall distance distribution
3. **Relative Location Pie Charts** - Per-feature breakdown of spatial relationships
4. **Feature Distribution Bar Chart** - Counts per feature type
5. **Query-wise Histograms** - Distance distributions per query
6. **Pairwise Venn Diagrams** - Query overlap (2 queries)
7. **Chow-Ruskey/UpSet Diagrams** - Multi-query overlap (3-5 queries)

### Current Limitations

- R dependency required
- Fixed plot layouts
- No interactive exploration
- Limited customization
- No batch/comparison mode across samples
- No enrichment statistics

---

## Proposed Future Visualizations

### Tier 1: Core Statistics (Essential)

#### 1.1 Annotation Summary Dashboard

**Purpose:** Single-page overview of annotation results

**Components:**
- Annotation rate gauge (% peaks annotated)
- Total peaks vs annotated peaks bar
- Top 5 features by count
- Distance summary statistics (median, IQR)

**Data Source:** `finalhits.txt`

**Mockup:**
```
┌─────────────────────────────────────────────────────────────┐
│                    ANNOTATION SUMMARY                        │
├───────────────┬─────────────────────────────────────────────┤
│  ANNOTATED    │  Feature Distribution        Distance Stats  │
│    ████████   │  gene      ████████ 45%     Median: 1,234 bp│
│     87.3%     │  promoter  █████    28%     Mean:   2,456 bp│
│               │  exon      ███      15%     IQR:    890 bp  │
│  12,456 / 14,271│ enhancer ██       12%                      │
└───────────────┴─────────────────────────────────────────────┘
```

**Integration:** New function in visualization module, called from main workflow

---

#### 1.2 Distance Distribution Plot (Enhanced)

**Purpose:** Understand proximity of peaks to features

**Enhancements over current:**
- Separate upstream/downstream distributions
- Log-scale option for wide distance ranges
- Kernel density estimate with rug plot
- Annotate median/percentiles
- Facet by feature type

**Data Source:** `allhits.txt` (distance, raw_distance columns)

**Mockup:**
```
         Upstream                    Downstream
    ←──────────────────|──────────────────→
         ▲                          ▲
        ▓▓▓                        ▓▓
       ▓▓▓▓▓                      ▓▓▓▓
      ▓▓▓▓▓▓▓▓                  ▓▓▓▓▓▓▓
    ▓▓▓▓▓▓▓▓▓▓▓▓▓▓          ▓▓▓▓▓▓▓▓▓▓▓▓▓
    ─────────────────────────────────────
    -10kb      -5kb       TSS       +5kb      +10kb
                         median ↓
```

---

#### 1.3 Relative Location Breakdown

**Purpose:** Visualize spatial relationships between peaks and features

**Plot Types:**
- Stacked bar chart (by feature type or query)
- Alluvial/Sankey diagram (feature → location flow)
- Treemap for hierarchical view

**Categories:**
- `PeakInsideFeature` - Peak fully contained in feature
- `FeatureInsidePeak` - Feature fully contained in peak
- `Upstream` - Peak upstream of feature
- `Downstream` - Peak downstream of feature
- `OverlapStart` - Peak overlaps feature start
- `OverlapEnd` - Peak overlaps feature end

**Data Source:** `relative_location` column

---

#### 1.4 Feature Type Distribution

**Purpose:** Show which genomic features peaks annotate to

**Plot Types:**
- Horizontal bar chart (sorted by count)
- Waffle chart (proportional squares)
- Donut chart with legend

**Grouping Options:**
- Overall
- Per query
- Per chromosome

---

### Tier 2: Comparative Analysis

#### 2.1 Query Comparison Matrix

**Purpose:** Compare annotation patterns across multiple queries

**Components:**
- Heatmap of feature counts (queries × features)
- Shared/unique peak counts
- Distance comparison boxplots

**Mockup:**
```
                   Query 1   Query 2   Query 3
           gene      45        32        28
        promoter     128       156       89
           exon      23        18        31
       enhancer      67        45        102

       Color intensity = count, normalized per row or column
```

**Data Source:** `allhits.txt` grouped by query_name

---

#### 2.2 Multi-Sample Comparison

**Purpose:** Compare annotation results across different samples/conditions

**Use Case:** ChIP-seq time course, treatment vs control, different cell types

**Components:**
- Annotation rate comparison
- Feature distribution shift (stacked bars, side-by-side)
- Distance distribution comparison (overlaid densities)
- Peak overlap UpSet plot

**Data Source:** Multiple `finalhits.txt` files from different runs

**Integration:** New CLI mode or separate utility script

---

#### 2.3 Query Overlap Visualization (Enhanced)

**Purpose:** Visualize which peaks are captured by which queries

**Enhancements:**
- UpSet plot (Python native, no R dependency)
- Interactive hover showing peak IDs
- Exportable peak lists per intersection
- Proportional Venn for 2-3 queries

**Data Source:** Peak IDs grouped by query

---

#### 2.4 Query Effectiveness Metrics

**Purpose:** Evaluate how well each query performs

**Metrics:**
- Unique hit rate (peaks only found by this query)
- Redundancy rate (overlap with other queries)
- Distance quality (median distance)
- Coverage (% of peaks with hits)

**Plot:** Radar chart or parallel coordinates

---

### Tier 3: Advanced Analytics

#### 3.1 Feature Enrichment Analysis

**Purpose:** Determine if peaks are enriched for specific features vs background

**Method:**
- Compare observed feature distribution to expected (based on GTF composition)
- Calculate fold enrichment
- Statistical testing (chi-square, Fisher's exact)

**Plot:**
- Volcano plot (fold enrichment vs p-value)
- Bar chart with significance markers
- Forest plot for multiple comparisons

**Mockup:**
```
Feature Enrichment Analysis

promoter  ████████████████████  3.2x ***
enhancer  ██████████████        2.1x **
gene      ████████              1.2x ns
exon      ██████                0.8x ns
intron    ████                  0.5x *

                   Fold Enrichment →
*** p < 0.001, ** p < 0.01, * p < 0.05
```

**Data Requirements:** GTF feature counts as background

---

#### 3.2 Chromosome Distribution

**Purpose:** Visualize annotation patterns across chromosomes

**Plot Types:**
- Bar chart of peaks per chromosome
- Ideogram with peak density overlay
- Circus plot for feature-chromosome relationships

**Insights:**
- Identify chromosome-specific binding patterns
- Detect annotation biases
- QC for expected distribution

---

#### 3.3 Peak Score Correlation

**Purpose:** Analyze relationship between peak quality and annotation

**Analyses:**
- Peak score vs distance scatter plot
- Score distribution for annotated vs unannotated peaks
- Score threshold impact on annotation rate

**Plot:**
```
    Peak Score
        ▲
    10  │    ·  ·
        │   · · · · ·
     5  │  · · · · · · · ·
        │ · · · · · · · · · ·
     0  └─────────────────────→
        0     5kb    10kb   15kb
              Distance to Feature
```

---

#### 3.4 Strand Concordance Analysis

**Purpose:** Analyze strand relationships between peaks and features

**Metrics:**
- Same strand vs opposite strand ratio
- Strand preference by feature type
- Antisense binding patterns

**Plot:**
- Diverging bar chart (same vs opposite)
- Heatmap (feature × strand relationship)

---

#### 3.5 Overlap Distribution Analysis

**Purpose:** Characterize the extent of peak-feature overlaps

**Plots:**
- Dual histogram (feat_ovl_peak and peak_ovl_feat)
- 2D density plot (both overlap metrics)
- Categorized: full overlap, partial, minimal

**Insights:**
- Binding site precision
- Feature boundary relationships
- Peak width optimization

---

#### 3.6 Genomic Context Heatmap

**Purpose:** Show peak-feature relationships across genomic windows

**Method:**
- Aggregate peaks relative to feature anchors
- Create metaprofile of binding around features
- Similar to ChIP-seq heatmaps but for annotations

**Plot:**
```
            -5kb    TSS    +5kb
    Peak 1   ░░░░░░░████░░░░░░░
    Peak 2   ░░░░░████████░░░░░
    Peak 3   ░░░░░░░░████░░░░░░
    Peak 4   ░░░░░░████░░░░░░░░
    ...

    Average  ▁▂▃▅███████▅▃▂▁
```

---

#### 3.7 Annotation Quality Report

**Purpose:** Comprehensive QC report for annotation runs

**Components:**
- Annotation rate by chromosome
- Distance distribution per feature type
- Multi-mapping rate (peaks with multiple annotations)
- Query coverage heatmap
- Suggested parameter adjustments

**Output:** Multi-page PDF or HTML report

---

### Tier 4: Interactive Visualizations (Future)

#### 4.1 Interactive Dashboard

**Technology:** Plotly Dash or Streamlit

**Features:**
- Filter by query, feature, distance
- Zoom and pan on distributions
- Click to see peak details
- Export filtered data

#### 4.2 Genome Browser Track Export

**Purpose:** Visualize annotations in IGV/UCSC

**Outputs:**
- BED track colored by feature type
- BigWig for distance heatmap
- BED with annotation metadata

#### 4.3 Report Builder

**Purpose:** Custom report generation

**Features:**
- Select plots to include
- Add custom annotations
- Export to PDF/HTML/PNG

---

## Integration Architecture

### Proposed Module Structure

```
uropa/
├── __init__.py
├── uropa.py              # Main entry point
├── annotation.py         # Core annotation logic
├── utils.py              # Utilities
└── visualization/        # NEW: Visualization subpackage
    ├── __init__.py
    ├── core.py           # Base classes, data loading
    ├── distributions.py  # Distance, overlap distributions
    ├── categories.py     # Feature, location breakdowns
    ├── comparisons.py    # Query/sample comparisons
    ├── enrichment.py     # Statistical enrichment
    ├── quality.py        # QC metrics and reports
    ├── export.py         # PDF/HTML/PNG generation
    └── interactive.py    # Dash/Streamlit apps (optional)
```

### Data Flow with Visualization

```
                    uropa.py
                       │
                       ▼
              ┌────────────────┐
              │  annotation.py │
              └────────┬───────┘
                       │
                       ▼
              ┌────────────────┐
              │   Output Files │
              │  (txt/bed/json)│
              └────────┬───────┘
                       │
          ┌────────────┼────────────┐
          │            │            │
          ▼            ▼            ▼
    ┌──────────┐ ┌──────────┐ ┌──────────┐
    │ summary  │ │ compare  │ │ enrich   │
    │ (--summary)│ (--compare)│(--enrich) │
    └────┬─────┘ └────┬─────┘ └────┬─────┘
         │            │            │
         └────────────┼────────────┘
                      ▼
              ┌────────────────┐
              │ visualization/ │
              │    package     │
              └────────┬───────┘
                       │
                       ▼
              ┌────────────────┐
              │  Output (PDF/  │
              │   HTML/PNG)    │
              └────────────────┘
```

### CLI Integration Options

#### Option A: Extend Current `--summary` Flag

```bash
# Current behavior (R-based)
uropa -i config.json --summary

# Enhanced with visualization level
uropa -i config.json --summary basic     # Tier 1 only
uropa -i config.json --summary full      # Tier 1 + 2
uropa -i config.json --summary advanced  # Tier 1 + 2 + 3
```

#### Option B: Separate Visualization Command

```bash
# Run annotation
uropa -i config.json -o results/

# Generate visualizations separately
uropa-viz results/finalhits.txt --output report.pdf
uropa-viz results/ --compare sample2/ --output comparison.pdf
uropa-viz results/allhits.txt --enrich background.gtf
```

#### Option C: Python API

```python
from uropa.visualization import AnnotationReport

# Load results
report = AnnotationReport("results/finalhits.txt")

# Generate specific plots
report.distance_distribution(save="distance.png")
report.feature_breakdown(save="features.png")

# Generate full report
report.generate_pdf("full_report.pdf", level="advanced")

# Compare multiple runs
from uropa.visualization import compare_runs
compare_runs(["sample1/", "sample2/"], output="comparison.pdf")
```

### Recommended Approach

**Hybrid: Option A + C**

1. Keep `--summary` for quick visualization (replace R with Python)
2. Add Python API for programmatic access
3. Add `--viz-level` flag for controlling detail level
4. Future: Add standalone `uropa-viz` command for post-hoc analysis

---

## Implementation Priorities

### Phase 1: Foundation (Replace R)

**Goal:** Remove R dependency, establish Python visualization base

**Tasks:**
1. Port `uropa_summary.R` functionality to Python
2. Implement Tier 1 core statistics
3. Create `visualization/` package structure
4. Update `--summary` flag to use Python

**Dependencies:** matplotlib, seaborn, pandas

**Timeline Estimate:** Core refactoring work

---

### Phase 2: Enhanced Analytics

**Goal:** Add statistical depth to visualizations

**Tasks:**
1. Implement Tier 2 comparison tools
2. Add enrichment analysis (Tier 3.1)
3. Create multi-sample comparison mode
4. Add chromosome distribution plots

**Dependencies:** scipy (statistics), upsetplot

---

### Phase 3: Quality & Reporting

**Goal:** Comprehensive QC and reporting

**Tasks:**
1. Implement quality metrics (Tier 3.7)
2. Add genome browser track export
3. Create customizable report builder
4. Add batch processing mode

---

### Phase 4: Interactive (Optional)

**Goal:** Interactive exploration capabilities

**Tasks:**
1. Create Streamlit/Dash dashboard
2. Add IGV.js integration
3. Implement live filtering

**Dependencies:** streamlit or dash, plotly

---

## Technical Considerations

### Performance

- Use pandas for efficient data manipulation
- Implement lazy loading for large datasets
- Cache computed statistics
- Support chunked processing for memory efficiency

### Compatibility

- Python 3.7+ (match current UROPA requirements)
- Optional dependencies for advanced features
- Graceful fallback if visualization deps missing

### Testing

- Unit tests for statistical calculations
- Visual regression tests for plots
- Integration tests with real annotation outputs

### Documentation

- Docstrings for all public functions
- Gallery of example plots
- Tutorial notebooks

---

## Appendix: Data Schema Reference

### finalhits.txt Columns

```
peak_chr        str     Chromosome of peak
peak_start      int     Peak start coordinate (0-based)
peak_end        int     Peak end coordinate
peak_id         str     Peak identifier
peak_score      str     Peak score (or '.')
peak_strand     str     Peak strand (+/-/.)
feature         str     GTF feature type
feat_start      int     Feature start coordinate
feat_end        int     Feature end coordinate
feat_strand     str     Feature strand
feat_anchor     str     Anchor point used (start/center/end)
distance        int     Distance to feature anchor
relative_location str   Spatial relationship category
feat_ovl_peak   float   Fraction of feature overlapped by peak
peak_ovl_feat   float   Fraction of peak overlapped by feature
[attributes]    str     GTF attributes (gene_id, gene_name, etc.)
query_name      str     Query identifier
```

### allhits.txt Additional Columns

Same as finalhits.txt, includes all valid annotations (not just best hit)

### Configuration JSON Schema

```json
{
  "queries": [
    {
      "name": "string",
      "feature": ["string"],
      "distance": [int, int],
      "feature_anchor": ["string"],
      "strand": "string",
      "relative_location": ["string"],
      "filter_attribute": "string",
      "attribute_values": ["string"],
      "internals": float
    }
  ],
  "show_attributes": ["string"] | "all",
  "priority": boolean,
  "gtf": "path",
  "bed": "path",
  "prefix": "string",
  "outdir": "path",
  "threads": int
}
```

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-22 | Claude | Initial draft |

---

*This document is a living roadmap and should be updated as requirements evolve and implementations progress.*
