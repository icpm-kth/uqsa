# rNormalPrior returns a random vector generator

The return value is a function that generates random vectors of the same
size as mean and sd from a multivariate normal distribution with
independent components with mean "mean" and standard deviation "sd". The
random vectors are returned as n rows of a matrix, where n is the only
argument of the returned function.

## Usage

``` r
rNormalPrior(mean, sd)
```

## Arguments

- mean:

  mean of the random variables (a vector)

- sd:

  standard deviation of the random variables (same size vector as mean)

## Value

an independentent multivariate normal random vector generating function:
rprior(n), where n is the requested number of vectors (rows)

## Examples

``` r
rnp<-rNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
rnp(12)
#>             [,1]       [,2]      [,3]
#>  [1,] -1.6674932 -0.0426894  1.444180
#>  [2,]  0.2003264  3.5639663  4.206825
#>  [3,] -0.4741987 -1.6850613  4.071467
#>  [4,] -0.0361994 -0.3343590 -0.753566
#>  [5,]  0.2784966  0.5178827  8.137944
#>  [6,] -0.6123089  0.2345288 -3.326241
#>  [7,]  0.2653067  1.9973628  4.111856
#>  [8,] -1.1648482  0.5261603  3.222167
#>  [9,]  0.6811271  1.2913329  1.848706
#> [10,]  0.7157739  1.9016031 -1.432438
#> [11,]  1.7167973  0.6066698  4.219390
#> [12,]  0.7866607  2.5848295  3.947625
```
