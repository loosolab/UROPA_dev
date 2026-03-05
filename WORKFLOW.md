# Enrichment/Depletion Visualization — Workflow Overview

## Pipeline

```
GTF ──> compute_genomic_background() ──> background proportions
                                              │
                                              ▼
finalhits.txt ──> compute_observed_proportions() ──> enrichment_test() ──> write TSV ──> plot PNG
  allhits.txt ──> compute_observed_proportions(deduplicate=True) ──> enrichment_test() ──> write TSV ──> plot PNG
```

## What it does

Compares **how often peaks land on each feature type** vs **how much of the genome that feature type covers**. A feature covering 10% of the genome but annotated to 30% of peaks is 3x enriched.

## Key steps

1. **Background**: Parse GTF, merge overlapping intervals per feature type/biotype, compute genome coverage proportions
2. **Observed**: Count peaks per category from finalhits (best hit) and allhits (deduplicated)
3. **Test**: Binomial test per category + Benjamini-Hochberg FDR correction
4. **Export**: Write results to TSV (intermediate data file)
5. **Plot**: Read TSV, generate horizontal bar plot (log2 fold enrichment, symmetric axes)

## Output structure

```
{prefix}_enrichment_finalhits/
    enrichment_feature_type.tsv + .png
    enrichment_biotype.tsv + .png
    enrichment_{query}_feature_type.tsv + .png
    enrichment_{query}_biotype.tsv + .png

{prefix}_enrichment_allhits/
    (same structure)
```

## CLI parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--summary` | Triggers enrichment plots (existing flag) | off |
| `--genome-size` | Total genome size in Mb | inferred from GTF |
| `--min-background` | Hide features covering < N% of genome from plots | 0 (show all) |

## Pseudocode

```
function run_enrichment_analysis(gtf, finalhits, output_prefix):

    # 1. Background: parse GTF, collect intervals per feature_type and biotype (genes only)
    #    Merge overlapping intervals per category per chrom, compute proportion = merged_bp / genome_size

    # 2. Observed: for each hits_file in [finalhits, allhits]:
    #    Count peaks per feature_type, biotype, and per (query, category)
    #    allhits: deduplicate by (peak_id, category); total_peaks = unique peak_ids

    # 3. Test: for each category: binomtest(k, n, background_prop, two-sided)
    #    BH-correct p-values, compute log2(fold_enrichment), write sorted results to TSV

    # 4. Plot: shared xlim = max(|log2fc|) * 1.15 across all TSVs
    #    For each TSV: filter 0-count and min-background, plot horizontal bars (red/blue), add significance stars
```

## Files changed

- **New**: `uropa/visualization.py` — all enrichment logic and plotting
- **Modified**: `uropa/uropa.py` — added CLI args, calls visualization in `--summary` block
