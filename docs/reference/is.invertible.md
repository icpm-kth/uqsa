# checks whether a given matrix is a valid, invertible fisherInformation

This matrix has to be symmetric and invertible. But, because the matrix
has a perhaps sketchy origin, it could be defective in all possible
ways.

## Usage

``` r
is.invertible(G = NULL, abs_tol = 1e-11)
```

## Arguments

- G:

  a matrix

- abs_tol:

  absolute tolerance for the reciprocal condition number of G

## Value

TRUE or FALSE
