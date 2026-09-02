# Kullback Leibler Divergence

This function calculates the Kullback Leibler Divergence \$D(P\\Q)\$
value between two distributions \$P\$ and \$Q\$, represented by their
two samples, `X` and `Y`. Both samples will have their density inferred.
The intended use-case is to compare 2D and 3D densities, e.g.: to find
interesting pairs of parameters within a bigger distribution.

## Usage

``` r
KLD(X, Y, de = c("copula", "ks", "mvtnorm"))
```

## Arguments

- X:

  sample from distribution P

- Y:

  sample from distribution Q

- de:

  density estimation mechanism (character scalar)

## Value

D, a scalar value, the Kullback Leibler Divergence

## Details

This estimate requires a method of density estimation, by default we use
the copula based methods
[fitCopula](https://icpm-kth.github.io/uqsa/reference/fitCopula.md) and
[dCopulaPrior](https://icpm-kth.github.io/uqsa/reference/dCopulaPrior.md)
(which has a fairly high accuracy, but can be quite slow).

Effects of setting `de` to

- `"ks"`: density is estimated via
  [`ks::kde()`](https://mvstat.net/ks/reference/kde.html)

- `"mvtnorm"`: density is estimated via manual kernel density
  estimation, using the multivariate Gaussian density
  [`mvtnorm::dmvnorm`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)
