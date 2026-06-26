# dNormalPrior creates the density function of a multivariate normal distribution with independent components

The returned density function takes vectors of the same size as mean and
sd. It returns the product of the components' one-dimensional normal
distribution, with mean "mean" and standard deviation "sd".

## Usage

``` r
dNormalPrior(mean, sd)
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
dnp<-dNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
dnp(c(0.5,1.5,2.5))
#> [1] 0.008926651
```
