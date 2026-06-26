# highCor returns ordered index-pairs of high to low correlation

This function uses the correlation matrix C of a sample X, orders all
values from the upper triangle of C (excluding the diagonal) from
highest to lowest correlation value and returns the indices as a
data.frame.

## Usage

``` r
highCor(C)
```

## Arguments

- C:

  the correlation matrix of a sample, no attributes need be present
  other than dim.

## Value

data.frame with columns i and j, representing the rows and columns of
high to low correlation pairs.

## Details

When truncated, the result can be used to plot only pairs with high
correlation.

## Examples

``` r
A <- matrix(
  c(
     1,  -1, 0.1,
    -1,   1, 0.4,
   0.1, 0.4,   1
  ),3,3
)
print(highCor(A))
#>   i j
#> 1 1 2
#> 2 2 3
#> 3 1 3
#> 4 1 1
#> 5 2 1
#> 6 3 1
#> 7 2 2
#> 8 3 2
#> 9 3 3
```
