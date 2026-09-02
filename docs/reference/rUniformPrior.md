# rUniformPrior returns a random vector generator

The return value is a function that generates random vectors of the same
size as ll and ul from a uniform distribution within the limits defined
by ul and ll. The random vectors are returned as n rows of a matrix,
where n is the only argument of the returned function.

## Usage

``` r
rUniformPrior(ll, ul)
```

## Arguments

- ll:

  lower limit of the random variables (a vector)

- ul:

  upper limit of the random variables (same size vector as ll)

## Value

a uniform random vector generating function: runiform(n), where n is the
requested number of vectors (rows)

## Examples

``` r
rup<-rUniformPrior(ll=c(0,1,2),ul=c(1,2,3))
rup(12)
#>             [,1]     [,2]     [,3]
#>  [1,] 0.83861256 1.315263 2.851363
#>  [2,] 0.23035291 1.542582 2.571535
#>  [3,] 0.80322417 1.418249 2.220343
#>  [4,] 0.79697265 1.304960 2.420876
#>  [5,] 0.48053315 1.636662 2.345794
#>  [6,] 0.75043391 1.106935 2.662378
#>  [7,] 0.45698796 1.614150 2.499770
#>  [8,] 0.01997855 1.473466 2.670255
#>  [9,] 0.12349001 1.260695 2.374025
#> [10,] 0.24126926 1.766407 2.733486
#> [11,] 0.53543044 1.185453 2.163538
#> [12,] 0.01318925 1.638889 2.928599
```
