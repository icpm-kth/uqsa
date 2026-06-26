# This function merges mpi-samples into one

When using MPI, we save the sample immediately into a file, each rank
saves to its own file. This function collects all of these smaller
samples into one. The samples should be saved with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

## Usage

``` r
loadSubSample_mpi(
  files,
  size = NA,
  selection = NA,
  mc.cores = parallel::detectCores()
)
```

## Arguments

- files:

  the files where the individual samples are stored

- size:

  sub sample size, if not set, the whole sample is returned

- selection:

  integer index vector or logical vector indicating which temperatures
  to return: beta\[selection\] is returned, in decreasing order of beta.

- mc.cores:

  defaults to the total number of cores, but can be reduced with this
  option.

## Value

a list of matrices, by temperature, concatenated.
