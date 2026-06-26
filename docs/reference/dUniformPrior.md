# dUniformPrior creates a uniform density function

The returned denisty function takes vectors of the same size as ll and
ul. It returns the product of the component's one-dimensional uniform
distribtions.

## Usage

``` r
dUniformPrior(ll, ul)
```

## Arguments

- ll:

  lower limit of the random variables (a vector)

- ul:

  upper limit of the random variables (same size vector as ll)

## Value

a probability density function on vectors withthe same length as ll and
ul.

## Examples

``` r
dup<-dUniformPrior(ll=c(0,1,2),ul=c(1,2,3))
dup(c(0.5,1.5,2.5))
#> [1] 1
```
