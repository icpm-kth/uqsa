# Modifies a value

This function does the same as `x <- x + sign*value`, but without
repeating `x`. The expression `modify(x) <- rnorm(length(x),0,1)` will
add Gaussian noise to it. This is meant as a replacement for the
`x += 1` syntax of C, it exists only for aesthetic reasons.

## Usage

``` r
modify(x, i = seq(NROW(x)), j = seq(NCOL(x)), sgn = +1) <- value
```

## Arguments

- x:

  a numeric value to be modified

- i:

  row-indices of `x` to be modified

- j:

  column-indices of `x` to be modified

- sgn:

  modification

- value:

  a numeric value of appropriate size, depending on `i` and `j`

## Value

The value of x is modified additively ihn place: `x <- x + sgn*value`

## Details

Specifically, this function should work for matrices, and it is possible
to supply row and column index vectors: `x[i,j]` will be modified.

This function is quite useful if `x` has a very long name, e.g.
`experiments[[1]]$func`.

## Examples

``` r
 x <- matrix(seq(12),3,4)
 modify(x,seq(2),seq(2)) <- 10
 print(x)
#>      [,1] [,2] [,3] [,4]
#> [1,]   11   14    7   10
#> [2,]   12   15    8   11
#> [3,]    3    6    9   12
```
