# ANALYSIS — Feature Enrichment/Depletion Plot

## 1. Problem Restatement

Add a visualization showing whether annotated peaks are **enriched** or **depleted** for specific feature types relative to what would be expected from the genomic background. The background is derived from the reference GTF: if feature X covers 50% of the genome and 70% of peaks are annotated to it, feature X is enriched (1.4x).

### Constraints
- Must parse the reference GTF to compute background proportions
- Normalize by feature length in genome
- Show both enrichment and depletion
- Integrates with existing `--summary` flag
- New code goes into `uropa/visualization.py`

## 2. Relevant Files/Modules

| File | Relevance |
|------|-----------|
| `uropa/annotation.py` | Annotation output format; `feature` column (GTF col 3) carries the feature type |
| `uropa/utils.py` | `subset_gtf()` filters GTF by feature type; `parse_bedlines()` reads BED; GTF chr-prefix handling |
| `uropa/uropa.py` | Orchestrator; `--summary` triggers visualization; writes finalhits/allhits output files |
| `uropa/visualization.py` | Target file for new visualization code (to be created) |
| `test_data/gencode.v29.annotation.chr19.gtf` | Reference GTF for testing |

## 3. Risks / Unknowns / Assumptions

### Overlapping features in the GTF
GTF features overlap hierarchically: genes contain transcripts which contain exons. Naively summing lengths per feature type **double-counts** shared genomic space. For example, an exon is always inside a transcript which is inside a gene.

### Genome size
The GTF doesn't contain total chromosome lengths — only annotated features. We need either:
- (a) Chromosome sizes from an external source (e.g., chromInfo, genome FASTA index)
- (b) Use only the annotated portion of the genome as the denominator
- (c) Derive approximate chromosome boundaries from max feature coordinates in the GTF

### Feature granularity
GTF column 3 has broad types (`gene`, `transcript`, `exon`, `CDS`, `UTR`, `start_codon`, `stop_codon`). Users may also care about **biotype** enrichment (e.g., `protein_coding` vs `lncRNA` genes), which lives in GTF attributes (`gene_type`, `transcript_type`).

### Statistical significance
A simple ratio (observed/expected) shows fold enrichment but doesn't indicate whether the difference is statistically significant. Options: binomial test, Fisher's exact test, chi-squared test.

### Unannotated peaks
Some peaks have no annotation (NA). These reduce the denominator for "proportion of peaks annotated to feature X" or should they be counted separately?

### Multiple queries
UROPA supports multiple queries with different feature types. The enrichment could be computed:
- Per-query (each query has its own set of target features)
- Globally across all queries (using finalhits which picks the best)

## 4. Proposed Approach

### Background computation
1. Parse the full GTF file
2. For each feature type (GTF column 3), merge overlapping intervals per chromosome to get non-redundant coverage in bp
3. Sum across chromosomes to get total bp covered by each feature type
4. Compute proportions relative to total genome size

### Observed computation
1. Parse the finalhits output file
2. Count how many peaks are annotated to each feature type
3. Compute proportions relative to total annotated peaks (or total peaks)

### Enrichment metric
- **Fold enrichment** = (observed proportion) / (expected proportion)
- Values > 1 = enriched, < 1 = depleted
- Log2 scale for symmetric visualization

### Visualization
- Bar plot with log2(fold enrichment) on y-axis, feature types on x-axis
- Bars colored by enrichment (red) vs depletion (blue)
- Optional: statistical significance markers (asterisks or p-value annotations)
- Optional: confidence intervals

## 5. Open Design Questions (for discussion)

1. **Feature granularity**: Should enrichment be computed at GTF feature-type level (gene/transcript/exon/CDS) or also at biotype level (protein_coding/lncRNA via gene_type attribute)? Or both?

2. **Genome denominator**: How to determine total genome size? Options:
   - (a) Require chromosome sizes file as input
   - (b) Infer from GTF max coordinates (approximate)
   - (c) Use only annotated genome portion

3. **Overlap handling for background**: Should overlapping features of the same type be merged (interval union) before computing coverage? (Yes seems correct, but adds complexity)

4. **Unannotated peaks**: Include in denominator or exclude?

5. **Statistical test**: Include significance testing, or just show fold enrichment?

6. **Multiple queries**: Per-query enrichment, global enrichment, or both?

7. **Output format**: Add to existing summary PDF, or separate output file?
