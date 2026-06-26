# Communicate with other ranks and swap beta

Given a current log-likelihood, temperature and step-size, this decides
whether to send or receive the same variables from a neighboring process
and swap temperatures with them.

## Usage

``` r
Rmpi_swap_temperatures(i, B, LL, H, r, comm, cs)
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
