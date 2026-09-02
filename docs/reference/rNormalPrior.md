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

an independent multivariate normal random vector generating function:
rprior(n), where n is the requested number of vectors (rows)

## Examples

``` r
rnp<-rNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
rnp(12)
#>              [,1]       [,2]       [,3]
#>  [1,]  1.23651533  3.0715204  0.4955232
#>  [2,]  0.01481459  0.7907957  0.9250191
#>  [3,]  2.23591826  1.3295411  8.1114178
#>  [4,]  0.73027659  1.2443280  5.6798020
#>  [5,] -0.63376135  4.6383517 -0.2788240
#>  [6,]  0.93135425  3.5743076  2.9020095
#>  [7,]  1.24999877 -1.0599455  6.7854998
#>  [8,]  2.59822929  5.2735788 -1.0380100
#>  [9,] -1.58289504  2.9588575  5.7162921
#> [10,] -1.01788411  0.7161362 -1.8499366
#> [11,] -0.64225651  1.6527338  4.0023513
#> [12,]  2.58065149 -1.3193510 -3.8654205
```
