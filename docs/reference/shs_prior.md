# Outpts the random sample on which to perform the Sobol-Homma-Saltelli global sensitivity analysis

Each parameter vector has length nPars, The sample consists of two
random (nSamples x nPars) matrices M1, M2 and a third (nSamples x nPars
x nPars) array N. N consists of nPars copies of M2, except that in each
M2-matrix one column has been replaced by the corresponding column of
M1. M1 and M2 consists of random numbers from a normal distribution.

## Usage

``` r
shs_prior(nSamples, rprior)
```

## Arguments

- nSamples:

  number of rows to return

- rprior:

  a function that samples from the prior distribution

## Value

a list with the components M1, M2 (both matrices) and N (a 3D-array).

## Details

These matrices provide prior distribution samples to be further
processed by the simulator, similar to this:

    sim <- simulator.c(experiments,modelName)
    fM1 <- t(sim(t(M1))[[1]]$state[,ti,])       # or similar

For details see: Halnes, Geir, et al. J. comp. neuroscience 27.3 (2009):
471.
