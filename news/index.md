# Changelog

## commaKit 0.2.0

### Current package baseline

- The canonical package version is `0.2.0`.
- The technical package identity has changed from `comma` to `commaKit`;
  repository, documentation, and installation URLs now use
  `carl-stone/commaKit`. The exported `commaData` class and
  `comma_example_data` dataset keep their names for API continuity.

### New features

- **`mod_context` rowData column** — every `commaData` object now stores
  a `mod_context` column in `rowData` that combines `mod_type` and
  `motif` (e.g. `"6mA_GATC"`, `"5mC_CCWGG"`). When `motif` is `NA`
  (Dorado/Megalodon callers), the fallback is the `mod_type` string
  alone (never `"6mA_NA"`). `mod_context` is required and enforced by
  `setValidity()`.

- **[`modContexts()`](https://carl-stone.github.io/commaKit/reference/modContexts.md)**
  — new exported S4 accessor that returns sorted unique modification
  context strings from a `commaData` object.

- **`expected_mod_contexts` constructor parameter** —
  [`commaData()`](https://carl-stone.github.io/commaKit/reference/commaData.md)
  now accepts `expected_mod_contexts`, a named list mapping mod types to
  allowed motifs
  (e.g. `list("6mA" = c("GATC", "ACCACC"), "5mC" = "CCWGG")`). Sites
  with unexpected mod_type × motif combinations are dropped at
  construction time with an informative message per mod type.

- **`subset(object, mod_context = ...)` filter** —
  [`subset()`](https://rdrr.io/r/base/subset.html) gains a `mod_context`
  parameter for filtering by modification context.

- **[`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md)
  loops by `mod_context`** — differential methylation is now computed
  independently for each `mod_context` group (e.g. <6mA@GATC> and
  <6mA@ACCACC> are tested separately). A `mod_context` parameter is
  added to
  [`diffMethyl()`](https://carl-stone.github.io/commaKit/reference/diffMethyl.md),
  [`results()`](https://carl-stone.github.io/commaKit/reference/results.md),
  and
  [`filterResults()`](https://carl-stone.github.io/commaKit/reference/filterResults.md)
  for context-specific extraction.

- **`mod_context` parameter added throughout** —
  [`methylomeSummary()`](https://carl-stone.github.io/commaKit/reference/methylomeSummary.md),
  [`slidingWindow()`](https://carl-stone.github.io/commaKit/reference/slidingWindow.md),
  [`mValues()`](https://carl-stone.github.io/commaKit/reference/mValues.md),
  [`writeBED()`](https://carl-stone.github.io/commaKit/reference/writeBED.md),
  and all eight `plot_*()` functions gain a `mod_context = NULL`
  parameter for context-level filtering.
  [`plot_tss_profile()`](https://carl-stone.github.io/commaKit/reference/plot_tss_profile.md)
  additionally supports `color_by = "mod_context"` and
  `facet_by = "mod_context"`.

### Breaking changes

- **Old `commaData` objects are invalid** — `rowData` must include a
  `mod_context` character column. Objects created with earlier versions
  will fail `validObject()`. Re-create from source files using the
  updated constructor.

### Bug fixes

- `diffMethyl(method = "methylkit")` no longer crashes with “object of
  type ‘closure’ is not subsettable” when a modification context
  contains sites where all samples have zero coverage after filtering.
  [`methylKit::unite()`](https://rdrr.io/pkg/methylKit/man/unite-methods.html)
  retains such sites in the united object; `calculateDiffMeth()` then
  calls `glm.fit` with all-zero weights and all-NaN response, which
  crashes. The wrapper now filters those sites out before calling
  `calculateDiffMeth` and assigns them `p = 1` (consistent with the null
  hypothesis). Regression test added.

### Package improvements

- Added explicit `Author` and `Maintainer` fields to DESCRIPTION for R
  4.6.0 compatibility.
- Replaced non-ASCII character in
  [`writeBED()`](https://carl-stone.github.io/commaKit/reference/writeBED.md)
  documentation with ASCII equivalent.
- Updated `.Rbuildignore` to exclude development files (`.claude/`,
  `.codex/`, `.letta/`, `.lteams/`, `AGENTS.md`, `PRD.md`, `VISION.md`,
  `ROADMAP.md`, `SPECS.md`).
- Renamed test files to match `test-functionName.R` convention:
  `test-coverageAnalysis.R` → `test-coverageDepth.R` +
  `test-varianceByDepth.R`, `test-enrichment.R` →
  `test-enrichMethylation.R`, `test-find_motif_sites.R` →
  `test-findMotifSites.R`, `test-load_annotation.R` →
  `test-loadAnnotation.R`, `test-m_values.R` → `test-mValues.R`,
  `test-plot_distribution.R` → `test-plot_methylation_distribution.R`.
- Updated README: corrected example data size (588 sites), added
  enrichment analysis and TSS profile workflow sections, updated feature
  list (eight plot functions, four DM backends, enrichment analysis),
  and expanded roadmap table.

### Historical pre-reset changelog entries

The entries below predate the `0.2.0` Schema v2 reset and are retained
for provenance. They are not the current package version sequence.
