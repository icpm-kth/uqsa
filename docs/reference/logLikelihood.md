# Default log-likelihood function

This returns a function f(simulations), which maps simulation results to
log(likelihood) values. The experiments are used implicitly; simulations
is a list as returned by rgsl::r_gsl_odeiv2_outer().

## Usage

``` r
logLikelihood(experiments)
```

## Arguments

- experiment:

  will be compared tp the simulation results
