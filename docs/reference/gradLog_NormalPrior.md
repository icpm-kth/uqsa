# Gradient of the logarithm of a normal prior

This makes a function that returns grad(log(dprior(x))) The returned
function implictly remembers the parameters of the normal distribution.

## Usage

``` r
gradLog_NormalPrior(mean = 0, sd = 1)
```

## Arguments

- mean:

  vector of mu values

- sd:

  vector of standard deviation values

## Value

g(x) a closure that remembers mean and sd from its creation
