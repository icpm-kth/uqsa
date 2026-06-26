# Collect statistical Replicas

`gatherReplicas` collects all sample-points, from all files, which are
assumed to be exact replicas. Replicas hav different random number seeds
(and possibly sample sizes).

## Usage

``` r
gatherReplicas(files)
```

## Arguments

- files:

  a list of file names

## Value

a Sample matrix, with effective sample size (auto-correlation thinned)

## Details

This function uses mclapply to process the files, which may be quicker
than `gatherSample`. The temperature `beta` is disregarded, assuming
that no parallel tempering was used. To facilitate the loading of a very
big sample, this function will analyse the auto-correltation within each
file and returned a thinned sub-sample of size `N/(2*tau_int)`
(returning the effective sample size). The value of tau_int is
calculated on the likelihood values, either with the hadron package, or
the bultin `acf` function. There is no need to further reduce the
result.

For small samples, it is better to load the entire sample and analyse it
in full. This function is intended for samples that are so big that they
challenge the memory of the machine.

This function is quicker if you have used trivial parallelism, *without*
MPI communication between the ranks (or another method of obtaining
several replicas, like forking or sequential repetition).

This function assumes that each supplied RDS file contains a matrix of
model MCMC parameters. The returned value `X` will be similar to effect
of `Reduce(...,rbind)` of all the smaller samples contained in the
individual files. The value `X` will have several attributes attached to
it:

- logLikelihood: log(likelihood(X\[i,\])), one value per row of X

- stepSize: the MCMC step size used in each given file#'

## Examples

``` r
rprior <- rNormalPrior(seq(3),seq(4,5)) # some nonsense
N <- 100
f <- c(tempfile(),tempfile())

## first fake sample
X <- rprior(N)
attr(X,"acceptanceRate") <- 0.23
## fake auto-correlation
attr(X,"logLikelihood") <- sqrt(seq(N)) + rnorm(N,-100,3)
saveRDS(X,file=f[1])

## second fake sample
X <- rprior(N)
attr(X,"acceptanceRate") <- 0.23
attr(X,"logLikelihood") <- sqrt(seq(N)) + rnorm(N,-100,3)
saveRDS(X,file=f[2])

Z <- gatherReplicas(f)
print(N)
#> [1] 100
print(dim(Z))
#> [1] 49  3
print(names(attributes(Z)))
#> [1] "dim"           "logLikelihood" "stepSize"      "tau"          
```
