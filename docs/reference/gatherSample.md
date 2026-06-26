# gatherSample collects all sample points, from all files, with the given temperature

This function assumes that each supplied RDS file contains a matrix of
model MCMC parameters, with an attribute called "beta" that lists the
temperature of each row.

## Usage

``` r
gatherSample(files, beta = 1, size = NA)
```

## Arguments

- files:

  a list of file names

- beta:

  the inverse temperture to extract sample for

- size:

  a size the is smaller than the actual sample size, if left unchanged,
  all sampled points are returned

## Value

a matrix of sampled points, all with the same temperature

## Details

This function selects and collects all rows, from all files with the
same (given) temperature.

This function should be used if you need to inspect only one of the
temperatures, not all of them. This function is similar to
[loadSample_mpi](https://icpm-kth.github.io/uqsa/reference/loadSample_mpi.md),
which returns all temperatures. But, whearas
[loadSample_mpi](https://icpm-kth.github.io/uqsa/reference/loadSample_mpi.md)
returns a list, this function returns the sample-matrix itself (because
the result of this function is conceptually similar to sampling on one
node, with one temperature).

## Examples

``` r
rprior <- rNormalPrior(seq(3),seq(4,5)) # some nonsense
N <- 100
f <- c(tempfile(),tempfile())

## first fake sample
X <- rprior(N)
attr(X,"beta") <- sample(1/seq(2)^2,N,replace=TRUE)
attr(X,"acceptanceRate") <- 0.23
attr(X,"swapRate") <- 0.1
attr(X,"logLikelihood") <- rnorm(N,-100,30)
saveRDS(X,file=f[1])

## second fake sample
X <- rprior(N)
attr(X,"beta") <- sample(1/seq(2)^2,N,replace=TRUE)
attr(X,"acceptanceRate") <- 0.23
attr(X,"swapRate") <- 0.1
attr(X,"logLikelihood") <- rnorm(N,-100,30)
saveRDS(X,file=f[2])

Z <- gatherSample(f,beta=1)
print(N)
#> [1] 100
print(dim(Z)) ## should be c(2*N,3)
#> [1] 101   3
print(names(attributes(Z)))
#> [1] "dim"           "beta"          "logLikelihood" "stepSize"     
```
