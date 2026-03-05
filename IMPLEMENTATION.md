# IMPLEMENTATION PLAN — Feature Enrichment/Depletion Visualization

## Step 1: Add `--genome-size` CLI argument to `uropa.py`

**File**: `uropa/uropa.py`

- Add `--genome-size` argument to the `additional` argument group
- Type: float, optional, in megabases (Mb)
- Help text: "Total genome size in Mb for enrichment analysis (default: inferred from GTF)"
- Pass through to visualization call

## Step 2: Create `uropa/visualization.py` — background computation

**File**: `uropa/visualization.py` (new)

Implement `compute_genomic_background(gtf_path, genome_size_mb=None)`:

1. Parse full GTF line by line (skip comment lines)
2. Collect intervals per feature type: `{feature_type: {chrom: [(start, end), ...]}`
3. Collect intervals per biotype (from "gene"-level features only):
   - Auto-detect attribute: try `gene_type`, then `gene_biotype`, then `biotype`
   - `{biotype_value: {chrom: [(start, end), ...]}}`
4. Track max coordinate per chromosome: `{chrom: max_coord}`
5. Merge overlapping intervals per category per chromosome (sort + merge algorithm)
6. Sum non-redundant bp per category
7. Compute genome size: `sum(chrom_max.values())` or `genome_size_mb * 1e6`
8. Return background dict with proportions

**Checkpoint**: Unit-testable with test GTF. Run on test data and print background proportions.

## Step 3: Create `uropa/visualization.py` — observed proportions

Implement `compute_observed_proportions(finalhits_path)`:

1. Read finalhits TSV, parse header to find column indices for `feature`, `name`, and the biotype attribute column (if present in show_attributes)
2. Count peaks per feature type (exclude NA rows from feature counts but include in total)
3. Count peaks per biotype (from the biotype attribute column, exclude NA)
4. Count total peaks (including unannotated)
5. If biotype column not in finalhits, also accept the detected biotype attribute name to look for
6. Return counts dict and total_peaks

**Checkpoint**: Run on test finalhits output.

## Step 4: Create `uropa/visualization.py` — statistical testing

Implement `enrichment_test(observed_counts, background_proportions, total_peaks)`:

1. For each category (feature type or biotype):
   - `k` = observed count
   - `n` = total_peaks
   - `p` = background proportion (expected)
   - Run `scipy.stats.binomtest(k, n, p, alternative='two-sided')`
   - Compute fold enrichment = `(k/n) / p`
   - Compute log2(fold enrichment)
2. Collect p-values, apply Benjamini-Hochberg correction (`scipy.stats.false_discovery_control` or manual)
3. Return list of dicts (or DataFrame) with: category, observed_count, observed_prop, expected_prop, fold_enrichment, log2_fc, p_value, p_adj

**Checkpoint**: Print enrichment table for test data.

## Step 5: Create `uropa/visualization.py` — plotting

Implement `plot_enrichment(enrichment_results, output_path, title)`:

1. Create horizontal bar plot using seaborn
2. x-axis: log2(fold enrichment), centered at 0
3. y-axis: feature type / biotype names, sorted by log2_fc
4. Bar colors: red gradient for enrichment, blue gradient for depletion
5. Significance markers at bar ends: `*` p<0.05, `**` p<0.01, `***` p<0.001
6. Gray dashed vertical line at x=0
7. Symmetric x-axis limits: `xlim = (-max_abs, +max_abs)`
8. Annotate with total peaks count and genome size
9. Save as PNG

**Checkpoint**: Visual inspection of test plot.

## Step 6: Create `uropa/visualization.py` — orchestrator

Implement `run_enrichment_analysis(gtf_path, finalhits_path, output_prefix, genome_size_mb=None, biotype_attribute=None)`:

1. Call `compute_genomic_background()`
2. Call `compute_observed_proportions()`
3. Determine shared x-axis limits across all plots
4. **Global feature-type enrichment**:
   - `enrichment_test()` on feature type counts vs background
   - `plot_enrichment()` → `{prefix}_enrichment_feature_type.png`
5. **Global biotype enrichment** (if biotype detected):
   - `enrichment_test()` on biotype counts vs background
   - `plot_enrichment()` → `{prefix}_enrichment_biotype.png`
6. **Per-query enrichment** (if multiple queries):
   - For each query name, filter observed counts to that query
   - `enrichment_test()` + `plot_enrichment()` for feature type
   - `enrichment_test()` + `plot_enrichment()` for biotype
   - → `{prefix}_enrichment_{query_name}_feature_type.png`
   - → `{prefix}_enrichment_{query_name}_biotype.png`

## Step 7: Integrate into `uropa.py`

**File**: `uropa/uropa.py`

1. Import `run_enrichment_analysis` from visualization module
2. In the `--summary` block (after line 575), add call:
   ```python
   from .visualization import run_enrichment_analysis
   run_enrichment_analysis(
       gtf_path=cfg_dict["gtf"],
       finalhits_path=output_prefix + "_finalhits.txt",
       output_prefix=output_prefix,
       genome_size_mb=args.genome_size,
   )
   ```
3. Keep existing R summary call intact (if available)

## Step 8: Test on test data

- Run UROPA with `--summary` on test data
- Verify PNG outputs are created
- Inspect plots for correctness
- Test with `--genome-size` override
- Test with single query (no per-query plots) and multiple queries
