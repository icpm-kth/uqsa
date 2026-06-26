# dCopulaPrior creates a prior probability density function

This function accepts the return list of fitCopula() or
makeIndepCopula() and creates a density function from it.

## Usage

``` r
dCopulaPrior(Copula)
```

## Arguments

- Copula:

  a list, as returned by fitCopula() or makeIndepCopula

## Value

a function that maps parameters (a vector) to probability density values
(scalar)

## Examples

``` r
x<-rnorm(300,mean=1,sd=2)
X<-matrix(x,100,3)
C<-fitCopula(X)
#> Loading required namespace: ks
d<-dCopulaPrior(C)
print(d(c(1,2,3)))
#> [1] 0.00277532
print(prod(sapply(c(1,2,3),FUN=dnorm,mean=1,sd=2)))
#> [1] 0.004248212
```
