# This creates a closure that simulates the model, similar to simulator.c

This is a shorter alternative to simulator.c (C backend). It also
returns the log-likelihood, Fisher Information, and the gradient of the
log-likelihood, under the assumption that the measurement error is
Gaussian. No attempt is made to parallelize this call, all simulations
will be done in sequence.

## Usage

``` r
simfi(
  experiments,
  odeModel,
  parMap = identity,
  method = 0,
  omit = 0,
  time.out = 1,
  num.steps = 0
)
```

## Arguments

- experiments:

  a list of experiments to simulate: inital values, inputs, time
  vectors, initial times

- odeModel:

  Either the ode object created by
  [as_ode](https://icpm-kth.github.io/uqsa/reference/as_ode.md) (with a
  shared library field inserted), or a string (with a comment indicating
  an .so file) which points out the model to simulate

- parMap:

  the model will be called with parMap(parABC); so any parameter
  transformation can happen there.

- method:

  the integration method as an integer (higher numbers are simpler
  methods, lower numbers are more advanced methods, 0 maps to 'msbdf')

- omit:

  integer, omit optional return values, in this order: Fisher
  Information, gradient of the log-likelihood, the log-likelihood,
  output functions. Omission includes all previous entries. `omit = 1`
  omits only the Fisher Information, `omit=3`, omits FI, grad-ll, and
  log-likelihood calculations.

- time.out:

  (in seconds); simulations are aborted at a time greater than this.

- num.steps:

  unlimited by default, setting this to a finite value can help to stop
  very stiff simulations early.

## Value

a closure that returns the model's output for a given parameter vector,
and approximate sensitivity matrices, for each state variable, function,
time-point, and parameter vector.

## Details

It returns a closure around: - experiments, - the model, and - parameter
mapping

The returned function depends only on parABC (the sampling parameters).

This version of the function does not use the parallel package at all
and cannot add noise to the simulations (unlike simulator.c).

A hopeless simulation can be stopped early using the settings
`num.steps` and `time.out`. The value of `num.steps` applies to every
continuous simulation stretch (e.g. betwee two events), teh count of
steps is reset whenever an event occurs or one simulation ends (between
different parameters and different experiments).

The `time.out` is given in seconds and can trigger at measurement time
points (when `t_wallclock > time.out`), not between points. How much
time has passed is checked when the integrator is stopped to record the
state.

The limit on the number of steps, on the other hand, is a feature of the
GSL ODE solvers and can trigger precisely.

## Examples

``` r
# \donttest{
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  o <- as_ode(m)
  ex <- experiments(m,o)
  C <- generate_code(o)
  c_path(o) <- write_c_code(C)
  so_path(o) <- shlib(o)
  s <- simfi(ex,o)
  y <- s(values(m$Parameter)) # simulates
  print(y)
#> number of simulation experiments: 3
#>                                      400nM 
#> ------------------------------------------ 
#>               cpuSeconds: 0.00033
#>                 numSteps: 269
#>                   status: 0
#>                    state: 2, 225, 1 (dim)
#>                     func: 1, 225, 1 (dim)
#>            logLikelihood: -864.384
#>        gradLogLikelihood: 5, 1 (dim)
#>        FisherInformation: 5, 5, 1 (dim)
#> 
#>                                      100nM 
#> ------------------------------------------ 
#>               cpuSeconds: 0.000335
#>                 numSteps: 248
#>                   status: 0
#>                    state: 2, 225, 1 (dim)
#>                     func: 1, 225, 1 (dim)
#>            logLikelihood: -840.397
#>        gradLogLikelihood: 5, 1 (dim)
#>        FisherInformation: 5, 5, 1 (dim)
#> 
#>                                       25nM 
#> ------------------------------------------ 
#>               cpuSeconds: 0.000314
#>                 numSteps: 237
#>                   status: 0
#>                    state: 2, 225, 1 (dim)
#>                     func: 1, 225, 1 (dim)
#>            logLikelihood: -786.464
#>        gradLogLikelihood: 5, 1 (dim)
#>        FisherInformation: 5, 5, 1 (dim)
#> 
#> experiments:  400nM, 100nM, 25nM 
# }
```
