# Get j-th column with names

When indexing a matrix or data.frame, rownames are lost. This function
will return a column of a matrix, as a vector (dropping rank so to
speak), but the vector will retain the rownames of the matrix

## Usage

``` r
column(m, j = 1)
```

## Arguments

- m:

  a matrix

- j:

  a column index

## Value

a named vector

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
u <- column(m$Parameter,"unit")
```
