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
#>  [1,] 0.99320374 1.645509 2.176973
#>  [2,] 0.09472838 1.863167 2.426114
#>  [3,] 0.84327161 1.252896 2.887468
#>  [4,] 0.25783174 1.468557 2.412131
#>  [5,] 0.57244787 1.324868 2.288161
#>  [6,] 0.37764145 1.458114 2.399713
#>  [7,] 0.22708202 1.372173 2.102753
#>  [8,] 0.39634991 1.118303 2.427391
#>  [9,] 0.44124259 1.693120 2.181922
#> [10,] 0.89493297 1.745078 2.011168
#> [11,] 0.30456573 1.787181 2.735519
#> [12,] 0.38429417 1.266346 2.083381
```
