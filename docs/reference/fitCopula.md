# Makes a Probability Density Estimate (from a sample)

Given a sample (from some probability distribution) this function makes
a Copula fit to the source distribution using the VineCopula package.

## Usage

``` r
fitCopula(X)
```

## Arguments

- X:

  sample that characterizes the traget distribution (rows)

## Value

as list: vineCop, U, Z, and Y where U are marginal probability samples,
Z are cummulative density values for U, and Y are the probability
density values of U.

## Examples

``` r
rprior <- rNormalPrior(c(1,2,3),c(4,5,6))
X <- rprior(1000)
C <- fitCopula(X)
rCopula <- rCopulaPrior(C)
Z <- rCopula(1000)
print(norm(cov(X) - cov(Z),"2")/norm(X,"2"))
#> [1] 0.009099276
print(abs(sum(colMeans(X) - colMeans(Z)))/sum(colMeans(X)))
#> [1] 0.01314942
```
