# Compute per-group beta means and treatment-reference delta beta

Compute per-group beta means and treatment-reference delta beta

## Usage

``` r
.computeDiffMethylGroupStats(methyl_mat, design)
```

## Arguments

- methyl_mat:

  Sites x samples methylation matrix.

- design:

  Output from .resolveDiffMethylDesign().

## Value

A list with group_means matrix and delta_beta vector.
