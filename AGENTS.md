# AGENTS.md — commaKit AI Context

> This file makes any AI coding tool smarter about commaKit. It works
> with Claude Code, Cursor, Copilot, Codex, Gemini CLI, and any tool
> that reads repo-level context files.

## What is commaKit?

**commaKit** (Comparative Microbial Methylomics Analysis Kit) is an
R/Bioconductor-style package for analyzing bacterial DNA methylation
from Oxford Nanopore sequencing data. It supports genome-wide
characterization of methylation patterns across multiple modification
types (6mA, 5mC, 4mC), annotation of methylation sites relative to
genomic features, and identification of differentially methylated sites
between conditions.

- **R package namespace**: `commaKit`
- **GitHub repo**: `carl-stone/commaKit`
- **Version**: 0.2.0
- **R**: \>= 4.3.0 (CI pinned to R 4.5)
- **License**: MIT

## Core Architecture

### Data structure: `commaData` (S4 class)

`commaData` extends `RangedSummarizedExperiment`. Genomic positions are
stored as `GRanges` in `rowRanges()` (1-bp ranges per site). Methylation
fractions and coverage counts are in assay matrices (rows = sites,
columns = samples).

**Key accessors:** - `methylation(object)` — methylation fractions
matrix (0–1) - `coverage(object)` — coverage count matrix (non-negative
integers) - `modTypes(object)` — unique modification types (e.g., “6mA”,
“5mC”) - `modContexts(object)` — unique modification+context units
(e.g., “6mA:GATC”, “5mC:CCWGG”) - `motifs(object)` — sequence motifs per
site - `siteInfo(object)` — flat DataFrame with chrom, position, strand,
mod_type, motif, and computed `site_key` - `sampleInfo(object)` — sample
metadata DataFrame - `caller(object)` — which caller produced the data
(“modkit”, “dorado”, “megalodon”)

**Important:** `mod_context` is **modification + context**, not just
“modification context.” `mod_type` is the chemical modification (6mA,
5mC, 4mC); motif/context is the sequence context (GATC, CCWGG, CpG);
`mod_context` is the two-part biological unit (e.g., `6mA:GATC`). The
package is modification-type agnostic — never assume a single mod_type.

**Alignment:** All alignment is by genomic position (GRanges) via
`findOverlaps()`, not by row names or string keys. `site_key` is a
computed display column (e.g., `chr1:512:+:6mA:GATC`), never used for
alignment.

### Analysis pipeline

1.  [`commaData()`](https://carl-stone.github.io/commaKit/reference/commaData.md)
    — construct from modkit/Dorado/Megalodon output
2.  [`annotateSites()`](https://carl-stone.github.io/commaKit/reference/annotateSites.md)
    — annotate sites to genomic features (list-columns:
    CharacterList/IntegerList/NumericList)
3.  [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
    — differential methylation (3 backends: beta-binomial/quasi_f,
    limma, methylKit)
4.  [`results()`](https://carl-stone.github.io/commaKit/reference/results.md)
    /
    [`filterResults()`](https://carl-stone.github.io/commaKit/reference/filterResults.md)
    — extract and filter differential methylation results
5.  [`enrichMethylation()`](https://carl-stone.github.io/commaKit/reference/enrichMethylation.md)
    — ORA/GSEA on differentially methylated genes
6.  [`methylomeSummary()`](https://carl-stone.github.io/commaKit/reference/methylomeSummary.md)
    — per-sample QC summary
7.  [`coverageDepth()`](https://carl-stone.github.io/commaKit/reference/coverageDepth.md)
    /
    [`varianceByDepth()`](https://carl-stone.github.io/commaKit/reference/varianceByDepth.md)
    — quality diagnostics
8.  [`slidingWindow()`](https://carl-stone.github.io/commaKit/reference/slidingWindow.md)
    — genome-wide smoothed profiles
9.  [`mValues()`](https://carl-stone.github.io/commaKit/reference/mValues.md)
    — M-value transformation for PCA
10. [`writeBED()`](https://carl-stone.github.io/commaKit/reference/writeBED.md)
    — export to BED

### Input parsers

- **modkit** (primary) — `parse_modkit.R` — modkit pileup BED output
- **Dorado** — `parse_dorado.R` — Dorado BAM with modification tags
- **Megalodon** (legacy) — `parse_megalodon.R` — Megalodon output

### Plot functions (8 total)

All return ggplot/patchwork except
[`plot_heatmap()`](https://carl-stone.github.io/commaKit/reference/plot_heatmap.md)
(ComplexHeatmap):
[`plot_coverage()`](https://carl-stone.github.io/commaKit/reference/plot_coverage.md),
[`plot_methylation_distribution()`](https://carl-stone.github.io/commaKit/reference/plot_methylation_distribution.md),
[`plot_pca()`](https://carl-stone.github.io/commaKit/reference/plot_pca.md),
[`plot_genome_track()`](https://carl-stone.github.io/commaKit/reference/plot_genome_track.md),
[`plot_metagene()`](https://carl-stone.github.io/commaKit/reference/plot_metagene.md),
[`plot_tss_profile()`](https://carl-stone.github.io/commaKit/reference/plot_tss_profile.md),
[`plot_volcano()`](https://carl-stone.github.io/commaKit/reference/plot_volcano.md),
[`plot_heatmap()`](https://carl-stone.github.io/commaKit/reference/plot_heatmap.md)

### diffMethyl backends

- **quasi_f** (default) — beta-binomial quasi-likelihood F-test. Most
  robust for small samples.
- **limma** — limma-voom with empirical Bayes moderation
- **methylKit** — wraps methylKit’s logistic regression

[`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
loops by `mod_context`, not `mod_type` — this prevents spurious pooling
across different sequence contexts. Effect sizes are always reported on
the beta scale (0–1), not M-value scale. Multiple testing correction is
genome-wide across all mod_contexts.

## R/Bioconductor Conventions

### Naming

| Category | Convention | Examples |
|----|----|----|
| S4 class | camelCase, lowercase first | `commaData` |
| Analysis functions | verbNoun() camelCase | [`annotateSites()`](https://carl-stone.github.io/commaKit/reference/annotateSites.md), [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md) |
| Plot functions | plot_noun() snake_case | [`plot_volcano()`](https://carl-stone.github.io/commaKit/reference/plot_volcano.md), [`plot_metagene()`](https://carl-stone.github.io/commaKit/reference/plot_metagene.md) |
| Internal functions | `.` prefix | `.parseBetaValues()`, [`.circularIndex()`](https://carl-stone.github.io/commaKit/reference/dot-circularIndex.md) |
| Arguments | snake_case | `mod_type`, `min_coverage` |
| Test files | test-functionName.R | `test-annotateSites.R` |

### Hard rules

- Every exported function accepts `commaData` as primary input
- Use
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  for interval overlap — never nested for-loops
- Return tidy data frames or updated commaData
- All `plot_*()` return ggplot/patchwork objects (except `plot_heatmap`
  → ComplexHeatmap)
- Genome size from `seqlengths(object)` (Seqinfo), never hardcoded
- Document every exported function with roxygen2: `@param`, `@return`,
  `@examples`
- Write testthat tests for every exported function
- Import `dplyr`, `tidyr` individually — never import `tidyverse`
  (Bioconductor requirement)
- Do not use
  [`purrr::map_dfr()`](https://purrr.tidyverse.org/reference/map_dfr.html)
  (superseded) — use `map(...) |> list_rbind()`
- Do not use
  [`purrr::map_dbl()`](https://purrr.tidyverse.org/reference/map.html)
  (superseded) — use [`vapply()`](https://rdrr.io/r/base/lapply.html)

### Known R gotchas in this codebase

1.  [`S4Vectors::rename()`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html)
    masks
    [`dplyr::rename()`](https://dplyr.tidyverse.org/reference/rename.html)
    — always use
    [`dplyr::rename()`](https://dplyr.tidyverse.org/reference/rename.html)
    explicitly
2.  [`matrixStats::count()`](https://rdrr.io/pkg/matrixStats/man/rowCounts.html)
    masks
    [`dplyr::count()`](https://dplyr.tidyverse.org/reference/count.html)
    — always use
    [`dplyr::count()`](https://dplyr.tidyverse.org/reference/count.html)
    explicitly
3.  [`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) /
    [`mclust::map()`](https://mclust-org.github.io/mclust/reference/map.html)
    collision — use [`lapply()`](https://rdrr.io/r/base/lapply.html) +
    [`purrr::list_rbind()`](https://purrr.tidyverse.org/reference/list_c.html)
    when mclust might be loaded
4.  [`diag()`](https://rdrr.io/r/base/diag.html) scalar trap: `diag(x)`
    creates x\*x identity matrix when x is scalar. Use `diag(x, nrow=1)`
    or `matrix(x, 1, 1)`
5.  methylKit crashes on zero-variance or all-zero-coverage sites — must
    filter before testing
6.  Non-ASCII characters (e.g., x) cause R CMD check notes
7.  CI runs with `--run-donttest`, so `\donttest{}` examples needing
    user-provided files will fail. Use `\dontrun{}` for examples needing
    modkit BED, GFF3, FASTA, BSgenome.
8.  `org.EcK12.eg.db` requires `::` syntax in examples

## Test data

- 588 sites: 393 x 6mA (GATC), 195 x 5mC (CCWGG)
- 6 samples: ctrl_1-3, treat_1-3
- Genome: chr_sim, 100 kb
- Ground truth: 30 of 393 6mA sites are differentially methylated
- Created by `data-raw/create_example_data.R` (`set.seed(1312)`)

## Development commands

``` bash
# Run tests
Rscript -e "devtools::test()"

# Regenerate documentation
Rscript -e "devtools::document()"

# Full R CMD check
Rscript -e "devtools::check()"

# R CMD check without vignettes (if pandoc unavailable)
Rscript -e "devtools::check(build_args = c('--no-build-vignettes'))"

# Install locally
Rscript -e "devtools::install()"
```

## Key design decisions

1.  **diffMethyl loops by mod_context, not mod_type** — prevents
    spurious pooling across different sequence contexts
2.  **commaData extends RangedSummarizedExperiment** — genomic positions
    in `rowRanges()` (GRanges), not rowData columns. Use
    [`siteInfo()`](https://carl-stone.github.io/commaKit/reference/siteInfo.md)
    for backward-compatible flat DataFrame access
3.  **Effect sizes always on beta scale (0–1)**, not M-value scale
4.  **Multiple testing correction is genome-wide** across all
    mod_contexts
5.  **annotateSites uses list-columns**
    (CharacterList/IntegerList/NumericList) — do NOT revert to
    single-match
6.  **Enrichment supports gene_role = target/regulator/both**
7.  **KEGG offline path** uses 2 bulk API calls then caches to RDS
8.  **Constructor uses `findOverlaps()` + mod_type/motif matching** for
    merge alignment — no string-key rownames on assay matrices
9.  **`.validateModType()`** internal helper validates mod_type values
    before subsetting — all exported functions with a `mod_type`
    argument use this
10. **`site_key` is a computed display column**, never used for
    alignment. Separator is `\001` (ASCII unit separator), not `\r`
11. **`plot_heatmap` row indexing**: use `as.integer(rownames(top_res))`
    to get original row indices, NOT `which(!is.na(...))`
12. **diffMethyl formula parameter**: accepts one-sided formulas with 2+
    levels. Multi-level formulas error with a clear message

## Current project status

- **Schema v2 milestone**: COMPLETE (merged to main 2026-05-21 as
  v0.2.0). All 14 issues closed.
- **Test Quality**: In progress. Audit done, PRs \#132/#133 strengthened
  tests. Plot tests still mostly smoke tests.
- **Code Quality Audits**: In progress. Thermonuclear review filed ~30
  findings as issues \#135-#163, index \#164. PRs \#165 and \#166
  merged.
- **Circle Ops** (issue \#122): Complete — circular genome boundary
  behavior audited and documented
- **Layered Assays** (issue \#118): Proposed for v0.3.0 — assay key
  system for multiple analysis runs
- **Technical rename to commaKit** (#168-#173): Complete — package/repo
  identity renamed from `comma`/`CoMMA` to `commaKit`

## Durable knowledge in the repo

- `dev/knowledge/test-quality.md` — full test audit with categorization
- `dev/knowledge/design-decisions.md` — architectural decision records
- `dev/knowledge/known-issues.md` — known issues and workarounds
- `dev/knowledge/git-discipline.md` — branching and versioning
  conventions
- `dev/knowledge/branching-releases.md` — release strategy
- `dev/ROADMAP.md` — strategic roadmap and milestone sequence
