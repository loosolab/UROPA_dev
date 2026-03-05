# Prompt: Feature Enrichment/Depletion Visualization for UROPA

## Task

Add an enrichment/depletion visualization to UROPA that compares the proportion of annotated peaks against the genomic background derived from the reference GTF. Normalize by the length of features in the genome. Show both enrichment and depletion with statistical significance.

## Requirements

### Two analysis levels (both produced by default)
1. **Feature type level** — GTF column 3 (gene, transcript, exon, CDS, ...)
2. **Biotype level** — GTF attribute (gene_type / gene_biotype / biotype, auto-detected). Use only gene-level features for biotype background to avoid double-counting nested features.

### Background computation
- Parse the full (unsubsetted) GTF to compute genomic coverage per category.
- Merge overlapping intervals of the same type per chromosome before summing coverage (no double-counting).
- Genome size: infer from max GTF coordinates per chromosome by default. Add `--genome-size` CLI parameter (in Mb) as alternative override.

### Observed computation
- Parse both the finalhits file (best hit per peak) and allhits file (all valid annotations).
- For allhits: deduplicate — count each peak only once per feature type/biotype (use peak_id).
- Denominator = total peaks including unannotated.

### Statistical test
- Binomial test (two-sided) per category: is the observed count significantly different from expected given background proportion?
- Benjamini-Hochberg FDR correction for multiple testing.
- Significance markers on plots: `*` p<0.05, `**` p<0.01, `***` p<0.001 (BH-adjusted).

### Per-query enrichment
- UROPA supports multiple named queries. Produce enrichment plots per query in addition to global plots.
- Same full-genome background for all; denominator = total peaks for all.

### Output structure
Produce two output folders per run:
- `{prefix}_enrichment_finalhits/` — from best-hit-per-peak table
- `{prefix}_enrichment_allhits/` — from all-hits table (deduplicated)

Each folder contains TSV + PNG pairs:
- `enrichment_feature_type.tsv/.png` — global feature type
- `enrichment_biotype.tsv/.png` — global biotype
- `enrichment_{query}_feature_type.tsv/.png` — per query
- `enrichment_{query}_biotype.tsv/.png` — per query

### TSV intermediate format
The TSV is the canonical data output. Metadata stored as comment lines:
```
#total_peaks=1000
#genome_size_bp=58599801
category	observed_count	observed_prop	expected_prop	fold_enrichment	log2_fc	p_value	p_adj	significance
exon	298	0.298000	0.119447	2.494826	1.318939	2.92e-51	2.92e-51	***
```
The plot function reads from the TSV, not from in-memory data.

### Plot design
- Horizontal bar plot using seaborn (matplotlib backend).
- x-axis: log2(fold enrichment), centered at 0.
- Bars colored: red (#c0392b) for enrichment, blue (#2980b9) for depletion.
- x-axis limits symmetric around 0 and consistent across ALL plots in a run (finalhits + allhits + all queries share the same limits).
- Categories with 0 observed peaks: always hidden from plots (kept in TSV).
- Categories with -inf log2_fc (0 observed): capped at the negative axis limit for display.
- `--min-background` CLI parameter (in %, default 0): optionally hide categories covering less than N% of the genome from plots (TSV always contains all).
- Annotation text: total peaks count and genome size.
- Significance legend in bottom-left.

### Integration
- Triggered by existing `--summary` flag.
- New code goes into `uropa/visualization.py`.
- New CLI args in `uropa/uropa.py`: `--genome-size` (float, Mb) and `--min-background` (float, %).
- Called in the `--summary` block of `uropa.py` before the legacy R summary.
- Dependencies: scipy (binomtest), statsmodels (multipletests), seaborn, matplotlib, numpy.

## Architecture

### `uropa/visualization.py` functions

```
compute_genomic_background(gtf_path, genome_size_mb=None) -> dict
    Parse GTF, merge overlapping intervals per feature type and per biotype
    (gene-level only). Return coverage proportions and genome size.

compute_observed_proportions(hits_path, biotype_attribute=None, deduplicate=False) -> dict
    Parse finalhits or allhits file. Count peaks per feature type, biotype,
    and per query. If deduplicate=True, count each peak only once per category
    using peak_id. Total peaks = unique peak IDs (includes unannotated).

enrichment_test(observed_counts, background, total_peaks) -> list of dict
    Binomial test per category + BH correction. Returns sorted list with
    category, observed_count/prop, expected_prop, fold_enrichment, log2_fc,
    p_value, p_adj, significance.

write_enrichment_tsv(results, tsv_path, total_peaks, genome_size_bp)
    Write results to TSV with metadata comment lines.

read_enrichment_tsv(tsv_path) -> (results, metadata)
    Read TSV back into list of dicts + metadata dict.

plot_enrichment_from_tsv(tsv_path, output_path, title, xlim_abs, min_background_pct)
    Read TSV, filter (0-count, min background), generate horizontal bar plot PNG.

_run_enrichment_for_hits(observed, background, output_dir, title_suffix, min_background_pct)
    Run enrichment for one hits file: test all categories, write TSVs.
    Returns list of (tsv_path, title, results) tuples.

run_enrichment_analysis(gtf_path, finalhits_path, output_prefix, genome_size_mb, min_background_pct)
    Orchestrator: compute background once, run for finalhits + allhits in
    separate folders, compute shared xlim across all results, generate all plots.
```

### `uropa/uropa.py` changes
- Add `--genome-size` and `--min-background` to argument parser.
- In `--summary` block: `from .visualization import run_enrichment_analysis` and call it with gtf_path, finalhits_path, output_prefix, and the new CLI args.
