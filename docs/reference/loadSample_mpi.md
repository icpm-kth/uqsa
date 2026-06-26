# This function merges mpi-samples into one

When using MPI, we save the sample immediately into a file, each rank
saves to its own file. This function is basically a wrapper with several
calls to `Reduce`, it collects all of these smaller samples into one.

## Usage

``` r
loadSample_mpi(files)
```

## Arguments

- files:

  the rds files where the individual samples are stored

## Value

a list of named items, with `$Sample` representing one matrix where all
file-samples are concatenated (with rbind).

## Details

The samples should have been saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). This function
extracts the attributes that MPI sampling typically attaches to a
sample. The sample itself and all of these attributes are returned as a
list.

If the samples contain different temperatures, then no attempt is made
to untangle or sort them.

NOTE: If the big result-sample doesn't fit into memory, this function
will crash. Samples can be quite large, depending on the problem size.

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

Z <- loadSample_mpi(f)
#> loading sample files with acceptances:
#> [1] 0.23 0.23
print(dim(Z$Sample))
#> [1] 200   3
print(names(Z))
#> [1] "Sample"         "beta"           "acceptanceRate" "swapRate"      
#> [5] "logLikelihood"  "betaSelection"  "uB"            
```
