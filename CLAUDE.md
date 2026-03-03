# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

UROPA (Universal RObust Peak Annotator) is a Python CLI tool for genomic region annotation. It annotates BED-format peaks against GTF reference files using configurable queries with parameters for feature type, distance, strand, anchor position, and attribute filtering. Supports multiprocessing via chunked peak annotation with pysam/tabix.

## Enviroment

For testing use the conda env uropa-dev

## Running UROPA

```bash
# With config file (multi-query)
uropa -i sample_config.json -t 4

# Direct CLI (single query)
uropa --bed test_data/genomic_regions.bed --gtf test_data/gencode.v29.annotation.chr19.gtf

# With visualization summary
uropa --bed test_data/genomic_regions.bed --gtf test_data/gencode.v29.annotation.chr19.gtf --summary
```

## Architecture

### Core Modules (under `uropa/`)

- **`uropa.py`** — CLI entry point and orchestrator. Parses args/config, validates input, prepares GTF (subset + tabix index), runs multiprocessed annotation via `mp.Pool`, writes output through a sorted queue writer, and optionally generates summary plots.

- **`annotation.py`** — Core annotation logic. `annotate_single_peak()` fetches GTF hits from tabix within distance range, validates each hit against query criteria (feature type, strand, distance, overlap, attributes, relative location), and returns valid annotations. `annotate_peaks()` processes a chunk of peaks and pushes results to the output queue. Key functions: `create_anno_dict()`, `distance_to_peak_center()`, `calculate_overlap()`, `get_relative_location()`.

- **`utils.py`** — Configuration parsing (`format_config`, `parse_json`), file validation, BED parsing (`parse_bedlines`), GTF utilities (`check_chr`, `subset_gtf`), sorted multiprocess file writer (`sorted_file_writer`), and `UROPALogger` (custom logger with queue-based multiprocess logging and custom DEBUG2 level).

### Data Flow

1. **Input**: BED file (peaks) + GTF file (annotations) + config (queries)
2. **Preparation**: GTF is optionally subsetted by feature type, then sorted/compressed/indexed with pysam tabix
3. **Annotation**: BED file is read in chunks; each chunk is annotated in parallel via `mp.Pool.apply_async`. Peak annotation uses tabix fetch within max distance, then validates hits against all query criteria
4. **Output**: Results are written via a sorted queue writer to maintain chunk order. Produces `{prefix}_allhits.txt/.bed` (all valid annotations) and `{prefix}_finalhits.txt/.bed` (best hit per peak). Optionally per-query output files with `--output-by-query`
5. **Visualization**: Optional `--summary` 

### Key Conventions

- **Strand-aware anchors**: For minus-strand features, "start" anchor maps to `feat_end` and "end" anchor maps to `feat_start` (see `annotation.py:67-69`).
- **Config format**: JSON with `queries` array, each query has keys like `feature`, `feature_anchor`, `distance` (array of [upstream, downstream]), `strand`, `relative_location`, `internals`, `filter_attribute`/`attribute_values`, `name`.
- **Output columns**: `peak_chr, peak_start, peak_end, peak_id, peak_score, peak_strand, feature, feat_start, feat_end, feat_strand, feat_anchor, distance, relative_location, feat_ovl_peak, peak_ovl_feat` + show_attributes + query name.
- **Multiprocessing**: Uses `mp.Pool` with memory-aware job throttling. Output ordering is maintained via indexed queue and `sorted_file_writer`.

### Test Data

`test_data/` contains `genomic_regions.bed` and `gencode.v29.annotation.chr19.gtf` for testing. No automated test suite exists; testing is done by running UROPA on test data and inspecting output.

Ignore existing R utilities.

# Default workflow

- Do not implement code immediately.
- Always start with:
(1) Planning phase (Codebase understanding if available, Architecture / design proposal), in exchange with the user (discussing design questions).
(2) Implementation plan.

- All new findings and decisions are saved into markdown files.

## Phase 1 — Analysis (required if codebase existing)
- Produce an ANALYSIS.md (or a Markdown section in chat) that includes:
1. Problem restatement + constraints
2. Relevant files/modules to inspect
3. Risks/unknowns + assumptions
4. Proposed approach options (if applicable) with tradeoffs

## Phase 2 — Architecture (required)
- Propose the target architecture and interfaces before touching code:
1. Module boundaries
2. Data flow
3. Key types/contracts
4. Testing strategy
- Write an ARCHITECTURE.md with the proposed architecture.
- Discuss important design questions with the user and update the ARCHITECTURE.md.

## Phase 3 — Plan then implement
- Write a step-by-step implementation plan with checkpoints, store the implementation plan in IMPLEMENTATION.md.

## PHASE 4 - Implementation
- Implement the IMPLEMENTATION.md step by step. Between each step the user will review the code and trigger the next step.

## New visualization related code:
- goes into uropa/visualization.py
- are executed with the --summary param
