# The MPI version of the mcmc function

this version of the MCMC function returns a Markov chain closure that
assumes that it is bein run in an MPI context: R was launched in an MPI
context, e.g. using

    mpirun -H localhost:8 -N 8 Rscript ...

and the pbdMPI package is installed. The chains shall communicate using
the provided `comm` object.

## Usage

``` r
mcmc_mpi(
  update,
  comm,
  swapDelay = 0,
  swapFunc = pbdMPI_bcast_reduce_temperatures
)
```

## Arguments

- update:

  an update function

- comm:

  an mpi comm which this function will use for send/receive operations

- swapDelay:

  swaps will be attempted every 2\*swapDelay+1 iterations
  [deprecated](https://rdrr.io/r/base/Deprecated.html)

- swapFunc:

  can be a custom function that does the MPI communication and decides
  whether or nopt to swap temperatures

## Value

an mcmc closure m(parMCMC,N,eps) that implicitly uses the supplied
update function

## Details

This function is intended for use within a parallel tempering approach
and MPI. For trivial parallelization (many chains), this is not at all
required, only a random number seed for each worker.

It is possible to supply a custom swap function, with the interface:

    swapFunc <- function(i, B, LL, H, r, comm, cs)

where `i` is the current iteration (for round robin rank choices), `B`
is the current beta value, `LL` the current log-likelihood (scalar) and
`H` the current step-size (scalar); `r`, `comm`, and `cs` are the MPI
rank, comm, and comm-size. The swap function returnsa list:
`list(B=,LL=,H=)` with the updated values (after swapping) or the old
values if the swap was rejected.

## Examples

``` r
## works in an MPI context
## similar to mcmc without _mpi prefix
if (FALSE) { # \dontrun{
  ## prepare the update function
  pt_mcmc <- mcmc_mpi(update, comm, swapDelay=0, swapFunc=pbdMPI_bcast_reduce_temperatures)
} # }
```
