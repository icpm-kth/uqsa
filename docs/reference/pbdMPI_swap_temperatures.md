# Communicate with other ranks and swap beta

Given a current log-likelihood, temperature and step-size, this decides
whether to send or receive the same variables from a neighboring process
and swap temperatures with them.

## Usage

``` r
pbdMPI_swap_temperatures(i, B, LL, H, r, comm, cs)
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

Only one of the two communicating ranks is allowed to make the decision
to swap, because a random variable is used to make this decision. Which
rank will make this swap decision needs to be determined somehow.
Currently we alternate this reposibility based on the current iteration
`i`.
