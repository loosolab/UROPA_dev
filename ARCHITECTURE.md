# ARCHITECTURE — Feature Enrichment/Depletion Visualization

## Overview

Two enrichment/depletion plots triggered by `--summary`:
1. **Feature-type level** — enrichment by GTF feature type (gene, exon, transcript, CDS, ...)
2. **Biotype level** — enrichment by GTF attribute (gene_type: protein_coding, lncRNA, ...)

Both compare observed peak annotation proportions against genomic background proportions derived from the GTF.

## Module Boundaries

### New: `uropa/visualization.py`

All visualization code lives here. No annotation logic, no file I/O beyond reading inputs.

```
visualization.py
├── compute_genomic_background(gtf_path, genome_size_mb=None) -> dict
│   Parses full GTF, merges overlapping intervals per feature type,
│   returns {feature_type: coverage_bp} and total_genome_bp.
│   Also extracts biotype coverage from GTF attributes.
│
├── compute_observed_proportions(finalhits_path, total_peaks) -> dict
│   Parses finalhits file, counts peaks per feature type (and biotype),
│   returns {feature_type: count} with total_peaks as denominator.
│
├── enrichment_test(observed_counts, background_proportions, total_peaks) -> DataFrame
│   Binomial test per feature/biotype. Returns fold_enrichment,
│   log2_fc, p_value, p_adjusted (BH).
│
├── plot_enrichment(enrichment_df, output_path, title) -> None
│   Bar plot: log2(fold enrichment) on y-axis, feature types on x-axis.
│   Colored by enrichment (red) / depletion (blue).
│   Significance asterisks from adjusted p-values.
│
└── run_enrichment_analysis(gtf_path, finalhits_path, output_prefix,
│                           genome_size_mb=None, query_names=None) -> None
    Orchestrator: calls the above functions, produces PNG files.
    Creates: {prefix}_enrichment_feature_type.png
             {prefix}_enrichment_biotype.png
             {prefix}_enrichment_{query_name}.png (per query)
```

### Modified: `uropa/uropa.py`

- Add `--genome-size` CLI argument (optional, in Mb)
- In the `--summary` block, call `run_enrichment_analysis()` after annotation

## Data Flow

```
┌──────────────┐     ┌─────────────────────┐
│  Full GTF    │────>│ compute_genomic_     │
│  (original)  │     │ background()         │
└──────────────┘     │                      │
                     │  Per feature type:   │
                     │  - merge intervals   │
                     │  - sum coverage bp   │
                     │                      │
                     │  Per biotype:        │
                     │  - merge intervals   │
                     │  - sum coverage bp   │
                     │                      │
                     │  Genome size:        │
                     │  - max(chr coords)   │
                     │    OR user-provided  │
                     └──────────┬──────────┘
                                │
                     ┌──────────▼──────────┐
                     │ enrichment_test()    │
                     │                      │
                     │  For each category:  │
                     │  - obs_prop = k/n    │
                     │  - exp_prop = bg/G   │◄──── compute_observed_
                     │  - binom_test(k,n,p) │      proportions()
                     │  - BH correction     │      (from finalhits)
                     └──────────┬──────────┘
                                │
                     ┌──────────▼──────────┐
                     │ plot_enrichment()    │
                     │                      │
                     │  log2(FC) bar plot   │
                     │  + significance      │
                     │  -> PNG output       │
                     └─────────────────────┘
```

## Key Types / Contracts

### Background dict
```python
{
    "feature_type": {
        "gene": {"coverage_bp": 45000000, "proportion": 0.35},
        "exon": {"coverage_bp": 12000000, "proportion": 0.09},
        ...
    },
    "biotype": {
        "protein_coding": {"coverage_bp": 30000000, "proportion": 0.23},
        "lncRNA": {"coverage_bp": 8000000, "proportion": 0.06},
        ...
    },
    "genome_size_bp": 128000000  # inferred or user-provided
}
```

### Enrichment DataFrame
```
feature_type | observed_count | observed_prop | expected_prop | fold_enrichment | log2_fc | p_value | p_adj
gene         | 350            | 0.70          | 0.35          | 2.0             | 1.0     | 1e-15   | 5e-15
exon         | 30             | 0.06          | 0.09          | 0.67            | -0.59   | 0.04    | 0.06
```

### Interval merging
For computing non-redundant coverage per feature type:
```python
# Input: list of (chrom, start, end) tuples for one feature type
# Output: total non-redundant bp
# Algorithm: sort by (chrom, start), merge overlapping intervals
```

## Genome Size Inference

Default: sum of max coordinates per chromosome from GTF.
```python
chrom_max = {}  # {chr1: 248956422, chr2: 242193529, ...}
genome_size_bp = sum(chrom_max.values())
```

Override: `--genome-size 3100` (in Mb) → `genome_size_bp = 3100 * 1e6`

## Biotype Attribute Detection

GTF attributes vary across sources. Detection priority:
1. `gene_type` (GENCODE)
2. `gene_biotype` (Ensembl)
3. `biotype` (older Ensembl)
4. Skip biotype plot if none found

## Per-Query Enrichment

For each named query:
1. Filter finalhits rows matching that query name
2. Use the same full-genome background
3. Denominator = total peaks (not just peaks matching the query)
4. Produce a separate PNG per query

## Plot Design

- Horizontal bar plot (feature names can be long)
- x-axis: log2(fold enrichment), centered at 0
- Bars colored: enrichment gradient (red) / depletion gradient (blue)
- Significance markers: `*` p<0.05, `**` p<0.01, `***` p<0.001 (BH-adjusted)
- Gray dashed line at x=0
- Title indicates "Feature Type Enrichment" or "Biotype Enrichment"
- Annotation: total peaks count, genome size used

## Dependencies

- `matplotlib` — plotting (already available in typical conda envs)
- `scipy.stats.binomtest` — binomial test
- `numpy` — numerical operations (already imported)
- `pandas` — DataFrame for enrichment results (optional, could use dicts)

## Design Decisions

1. **Biotype background**: Use only "gene"-level features for biotype background computation. Genes are the natural unit for biotype classification and avoid double-counting from nested features (exons inside transcripts inside genes).

2. **Minimum count threshold**: Show all feature types / biotypes — no filtering. The binomial test p-values will indicate which are significant.

3. **Plot library**: seaborn (with matplotlib backend).

4. **Axis symmetry**: x-axis limits are symmetric around 0 (same magnitude for + and -) and consistent across all plots in a run.
