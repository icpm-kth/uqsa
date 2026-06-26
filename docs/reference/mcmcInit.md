# Initialize the Markov chain

This function must append all required attributes to the MCMC varible,
for the Markov chain to update correctly.

## Usage

``` r
mcmcInit(
  beta,
  parMCMC,
  simulate,
  logLikelihood,
  dprior,
  gradLogLikelihood = NULL,
  gprior = NULL,
  fisherInformation = NULL
)
```

## Arguments

- beta:

  inverse temperature for the Markov chain (parallel tempering)

- parMCMC:

  a plain starting value for the Markov chain

- logLikelihood:

  a function that maps simulations to logLikelihood values

- gradLogLikelihood:

  the gradient function of the logLikelihood (optional) -- only if the
  algorithm requires it

- fisherInformation:

  a function that calculates the Fisher Information matrix

## Value

the same starting parameter vector, but with attributes.
