# gNormalPrior creates the gradient function of a multivariate normal distribution with independent components, in log-space

The returned density function takes vectors of the same size as mean and
sd. It returns the gradient of the logarithm of the multivariate normal
distribution, with mean `mean` and standard deviation `sd`.

## Usage

``` r
gNormalPrior(mean, sd)
```

## Arguments

- mean:

  mean of the random variables (a vector)

- sd:

  standard deviation of the random variables (same size vector as mean)

## Value

a probability density function on vectors with the same length as mean
and sd.

## Examples

``` r
gnp <- gNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
gnp(c(0.5,1.5,2.5))
#> [1] -0.50000000 -0.12500000 -0.05555556
```
