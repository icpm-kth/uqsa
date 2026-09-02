# Creates a prior probability density function

This function accepts the return list of
[fitCopula](https://icpm-kth.github.io/uqsa/reference/fitCopula.md) and
creates a density function from it.

## Usage

``` r
dCopulaPrior(Copula)
```

## Arguments

- Copula:

  a list, as returned by
  [fitCopula](https://icpm-kth.github.io/uqsa/reference/fitCopula.md)

## Value

a function that maps parameters (a vector) to probability density values
(scalar)

## Examples

``` r
x<-rnorm(300,mean=1,sd=2)
X<-matrix(x,100,3)
C<-fitCopula(X)
d<-dCopulaPrior(C)
print(d(c(1,2,3)))
#> [1] 0.003934713
print(prod(sapply(c(1,2,3),FUN=dnorm,mean=1,sd=2)))
#> [1] 0.004248212
```
