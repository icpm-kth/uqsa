# Default log-likelihood function

This returns a function f(simulations), which maps simulation results to
log(likelihood) values. The experiments are used implicitly; simulations
is a list as returned by rgsl::r_gsl_odeiv2_outer().

## Usage

``` r
logLikelihoodFunc(experiments, perExpLLF = NULL, simpleUserLLF = NULL)
```

## Arguments

- experiments:

  will be compared tp the simulation results

- perExpLLF:

  (optional) a user supplied function with the interface
  `perExpLLF(p,s,e)`, where `p` are the parameters, `s` are the
  simulations and `e` are the experiments (with data). Supply this
  function if some of your experiments need to be normalized by the
  other experiments (and other complex cases).

- simpleUserLLF:

  (optional) a user supplied function that is used instead of the
  default sum of ((y-h)/stdv)^2 terms. The interface is:
  `simpleUserLLF(y,h,stdv,name=NULL)`, where each of them is an
  N-M-matrix where N is the dimensionality of the model output and M the
  number of data time-points. Here, `y` is `t(experiments[[i]]$data)`
  and may contain NA values. This function should also accept an
  optional *name* argument (this is the name of the experiment this
  function is currently called for).

## Value

`llf(parMCMC)`, a closure (function) of the mcmc-variable: parMCMC;
returns a scalar log-likelihood value. Alternatively, the user can
define such a function: `parMCMC -> log(Likelihood(parMCMC))`, and use
that during sampling. A test simulation of `p`: `y <- simulate(p)` will
reveal which values the simulator produces. These values will be
attached to p during sampling, as an attribute.
[mcmc_init](https://icpm-kth.github.io/uqsa/reference/mcmc_init.md) will
attach the same values for the initial Markov chain state. The
log-likelihood function can use these attributes.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simulator.c(ex,o,omit=0)
p <- values(m$Parameter)
attr(p,"simulations") <- s(p)
## this function is fairly flexible and accepts some user settings
llf <- logLikelihoodFunc(ex)
#> experiments contain 675 non-missing values
print(llf(p))
#> Warning: In 'Ops' : non-'errors' operand automatically coerced to an 'errors' object with no uncertainty
#> [1] -2491.245
## this function uses the values from the solver:
print(ll(p))
#> [1] -2491.245
```
