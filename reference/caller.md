# Accessor for the methylation caller

Returns the name of the methylation caller that produced the data (e.g.,
`"modkit"`, `"megalodon"`, or `"dorado"`). The caller is stored in
`metadata(object)` at construction time.

## Usage

``` r
caller(object)

# S4 method for class 'commaData'
caller(object)
```

## Arguments

- object:

  A `commaData` object.

## Value

A character string naming the caller, or `NA` if not stored (e.g.,
objects created before caller storage was implemented).

## Examples

``` r
data(comma_example_data)
caller(comma_example_data)
#> [1] "modkit"
```
