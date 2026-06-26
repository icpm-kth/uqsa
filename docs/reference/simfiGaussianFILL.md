# SMMALA Extract the approximate Fisher infomration from the simfi results

This function extracts the approximate Fisher information matrix G of
the log-likelihood and transforms it using the Jacobian of the parameter
map between Markov chain variables and model parameters.

## Usage

``` r
simfiGaussianFILL(
  ParMapJac = function(x) {
     diag(1, length(x))
 }
)
```

## Arguments

- ParMapJac:

  Jacobian of the parMap function

## Value

the approximate Fisher information of the Gaussian log-likelihood with
respect to the MCMC variable, useful for SMMALA

## Details

The simfi() values are with respect to the raw model parameters, while
this function rephrases them in terms of the Markov chain's position.
