# GenomeInfoDb genome accessor compatibility method for commaData

Historically, `genome(commaData)` returned chromosome sizes. New code
should use
[`genomeSizes`](https://carl-stone.github.io/commaKit/reference/genomeSizes.md)
for chromosome lengths and `GenomeInfoDb::genome(seqinfo(object))` for
genome build/name metadata. This method preserves the historical
size-vector behavior for compatibility.

## Usage

``` r
# S4 method for class 'commaData'
genome(x)
```

## Arguments

- x:

  A `commaData` object.

## Value

A named integer vector of chromosome sizes, or `NULL` if no size
information was provided at construction.

## Examples

``` r
data(comma_example_data)
genome(comma_example_data)
#> chr_sim 
#>  100000 
```
