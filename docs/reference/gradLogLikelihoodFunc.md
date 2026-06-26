# Default log-likelihood function, gradient

This returns a function g(x,simulations), which maps simulation results
and the MCMC variables x to the gradient of log(likelihood) values withj
respect to x. The experiments are used implicitly; simulations is a list
as returned by rgsl::r_gsl_odeiv2_outer().

## Usage

``` r
gradLogLikelihoodFunc(
  experiments,
  parMap = identity,
  parMapJac = function(x) {
     diag(1, length(x))
 }
)
```

## Arguments

- experiment:

  will be compared tp the simulation results
