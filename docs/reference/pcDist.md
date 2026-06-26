# plots a sample in parallel coordinates

This function makes a plot that is quite similar to parallel
coordinates. It includes information about the prior as error-bars,
centered around th eprior's median.

## Usage

``` r
pcDist(posterior, prior, color = rgb(0.5, 0.5, 0.5, 0.05), ...)
```

## Arguments

- posterior:

  a matrix, with N rows (sample-members), and M columns (different model
  parameters). The columns must be named.

- prior:

  a data.frame with at least \$median, and \$stdv columns. This
  data.frame may also include the fields: color, and colorOutline to
  change the prior error-bars.

- color:

  the color of the sample lines, should have some transparency.

- ...:

  parameters are passed to matplot.

## Value

produces a plot

## Examples

``` r
rprior <- rNormalPrior(c(-1,0,1),c(1,2,3))
A <- matrix(rnorm(9),3,3)
A <- (A + t(A))^2/norm(A)^2
X <- rprior(1000)
Z <- X %*% A
colnames(Z) <- letters[seq(3)]
pr <- data.frame(median=apply(X,2,median),stdv=apply(X,2,sd))
pcDist(Z,pr)
```
