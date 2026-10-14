# Returns a random vector generator

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

an independent multivariate normal random vector generating function:
rprior(n), where n is the requested number of vectors (rows)

## Examples

``` r
rnp<-rNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
rnp(12)
#>              [,1]       [,2]       [,3]
#>  [1,] -1.61202347  1.5378328  1.6958883
#>  [2,] -1.07835464 -0.1603457 -0.4584329
#>  [3,] -0.12222376  0.2729641  5.2928250
#>  [4,]  1.07169584  0.5730665 -1.1282050
#>  [5,]  0.31435000  4.6946895  2.4037720
#>  [6,]  1.39894874  0.1543286  4.8904265
#>  [7,] -0.52140003 -1.7570553 -1.6423737
#>  [8,] -0.79180154 -1.4568923  1.3120742
#>  [9,] -0.02373642 -4.0431584 -0.5145879
#> [10,]  0.18570900  1.9083450 -2.5021438
#> [11,]  0.55674870  1.8623681  4.9973316
#> [12,]  0.49183051  2.1719005  5.6278598
```
