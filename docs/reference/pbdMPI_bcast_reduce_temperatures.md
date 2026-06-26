# Broadcast to other ranks and swap temperatures with any of them

Using this function, at most two ranks will swap.

## Usage

``` r
pbdMPI_bcast_reduce_temperatures(i, B, LL, H, r, comm, cs)
```

## Arguments

- i:

  MCMC iteration

- B:

  inverse temperature (parallel tempering)

- LL:

  log-likelihood value of current point

- H:

  algorithm's step size (often called epsilon in literature)

- r:

  MPI rank

- comm:

  MPI communicator

- cs:

  MPI comm size

## Details

Given a current log-likelihood, temperature and step-size, this function
will broadcast a log-likelihood value to all other ranks and they can
each decide to swap temperatures with the root process. Root is cycled
around all ranks (round-robin).

Each other rank is allowed to make the offer to swap. The root process
decides which rank to swap with.
