# Deprecated coverage accessor for commaData objects

`coverage(commaData)` is deprecated because `coverage()` is an
established IRanges/GenomicRanges generic for computing genomic
coverage, not for retrieving an assay matrix. Use
[`siteCoverage`](https://carl-stone.github.io/comma/reference/siteCoverage.md)
instead.

## Usage

``` r
# S4 method for class 'commaData'
coverage(x, shift = 0L, width = NULL, weight = 1L, ...)
```

## Arguments

- x:

  A `commaData` object.

- shift, width, weight, ...:

  Inherited from
  [`IRanges::coverage`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).
  These arguments are not meaningful for the commaData assay accessor
  and must be left at their defaults.

## Value

An integer matrix with rows corresponding to methylation sites and
columns corresponding to samples.
