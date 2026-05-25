# Resolve the two-level differential methylation design contract

Internal helper shared by diffMethyl() and all statistical backends.
comma currently supports one two-level contrast per diffMethyl() call.
Multi-level primary variables must be modeled in a future
explicit-contrast API.

## Usage

``` r
.resolveDiffMethylDesign(coldata, formula, ref_level = NULL)
```

## Arguments

- coldata:

  Sample metadata as a data.frame-like object.

- formula:

  One-sided design formula.

- ref_level:

  Optional reference level.

## Value

A list containing primary_var, ref_level, treat_level, cond_levels,
cond, and group_idx.
