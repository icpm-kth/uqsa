# This creates a closure that simulates the model

This function will use the
[parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html) to do the
simulations simultaneously. Set `options(mc.cores=detectCores())` or a
similar sensible value: `options(mc.cores=length(experiments))`

## Usage

``` r
simulator.c(
  experiments,
  modelName,
  parMap = identity,
  noise = FALSE,
  omit = 3,
  method = 0,
  time.out = 1,
  num.steps = 0
)
```

## Arguments

- experiments:

  a list of experiments to simulate: inital values, inputs, time
  vectors, initial times

- modelName:

  a string (with optional comment indicating an .so file) which points
  out the model to simulate if modelName is a cme object, the simulation
  will be done stochasitcally

- parMap:

  the model will be called with parMap(parABC); so any parameter
  transformation can happen there.

- noise:

  boolean variable. If `noise=TRUE`, Gaussian noise is added to the
  output of the simulations. The standard deviation of the Gaussian
  noise is equal to the measurement error. If `noise=FALSE` the output
  is the deterministic solution of the ODE system. noise and sensitivity
  calculations are mutually exclusive.

- omit:

  `omit=0` returns all optional return values form the simulator,
  `omit=1` will not calculate the fisher information (and thus not
  return it), `omit=2` will omit the gradient of the log-likelihood, and
  `omit=3` will omit the likelihood calculations alltogether. Omission
  is cumulative: `omit=3` omits all the previously mentioned optional
  quantities.

- method:

  an integer offset, integration method (for ODE models), see
  [method](https://icpm-kth.github.io/uqsa/reference/method.md) and
  [name_method](https://icpm-kth.github.io/uqsa/reference/name_method.md)

- time.out:

  in seconds, for early stops.

- num.steps:

  maximum number of steps taken by the integrator (in the case of ODEs),
  or maximum number of total reaction-steps performed by the Gillespie
  algorithm (over time) for stochastic models.

## Value

a closure that returns the model's output for a given parameter vector

## Details

It returns a closure around: - experiments, - the model, and - parameter
mapping

The returned function depends only on the parameter vector (or matrix if
more than one simulation per experiment is desired). The parameter
vector this simnulator accepts is probably derived from the sampling
space of a Bayesian method \\\theta\\, so in the list of arguments, it
is called `parABC` or (parMCMC would also have been a valid choice).
These sampling parameters can be mapped to values the simulator can use
via `parMap`. `parModel <- parMap(parABC)`, where the ODE model is
expected to work with `parModel`. The model can be specified by name
(with a comment indicating a file location)

Some return values are optional and omiting them saves time.

## Examples

``` r
# \donttest{
  requireNamespace("errors")
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  o <- as_ode(m)
  ex <- experiments(m,o)
  C <- generate_code(o)
  c_path(o) <- write_c_code(C)
  so_path(o) <- shlib(o)
  s <- simulator.c(ex,o)
  y <- s(values(m$Parameter))
# }
```
