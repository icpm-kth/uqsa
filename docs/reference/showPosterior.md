# showPosterior makes a pairs plot for a sample

This function will display the difference between the posterior and
prior by plotting the posterior as shaded density plots and the prior as
contour lines of level sets. If the two are identical, the lines will be
invisible as they blend into the density plot. Otherwise the contour
lines will show up as a distinct feature.

## Usage

``` r
showPosterior(posterior, prior, ...)
```

## Arguments

- posterior:

  a matrix, each row is a sample member

- prior:

  a matrix of the same size as the posterior

- ...:

  passed to [`graphics::pairs()`](https://rdrr.io/r/graphics/pairs.html)

## Value

pairs plot object

## Examples

``` r
# \donttest{
rprior <- rNormalPrior(c(-1,0,1),c(1,2,3))
A <- matrix(rnorm(9),3,3)
A <- (A + t(A))^2/norm(A)^2
X <- rprior(1000)
Z <- X %*% A
colnames(Z) <- letters[seq(3)]
colnames(X) <- letters[seq(3)]
showPosterior(Z,X) # this can take a while

# }
```
