# Markov Chain Monte Carlo

This function creates an MCMC function for a given set of experiments.
The Markov chains have no communication between them if more than one is
created using this mechanism.

## Usage

``` r
mcmc(update)
```

## Arguments

- update:

  and update function

## Value

M(initiPar,N), a function of initial starting values and number of
Markov chain steps

## Details

The algorithm is entirely determined by the update function. Any
intermediate values that updates requires aside from simulation results
have to be attributes of the MCMC variable: parMCMC.

The update function: update(parGiven) -\> parUpdate depends only on the
given parameters, all other dependencies have to be either implicit (as
a closure) or attributes of parGiven.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- model_from_tsv(uqsa_example("AKAP79"))
rwm <- high_level_metropolis(m) # "random walk", metropolis algorithm
p <- rwm %@% "init"             # a valid starting point
smallSample <- rwm(rwm %@% "init",600,1e-4)
pairs(smallSample[,seq(6)])
} # }
```
