# SMMALA Extract the gradient of the log-likelihood from the simfi solver's return value

This function extracts the approximate gradient of the log-likelihood
and transforms the gradient using the Jacobian of the parameter map
between Markov chain variables and model parameters.

## Usage

``` r
simfiGaussianGradLogLikelihood(
  ParMapJac = function(x) {
     diag(1, length(x))
 }
)
```

## Arguments

- ParMapJac:

  Jacobian of the parMap function

## Value

the gradient of the Gaussian log-likelihood with respect to the MCMC
variable

## Details

The simfi() gradient is with respect to the raw model parameters
