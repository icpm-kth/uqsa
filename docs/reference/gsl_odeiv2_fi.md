# simulates an ode model with extra work

This function calls a C function which solves an initial value problem,
calculates the sensitivity of the solution, log-likelihood value `ll`,
gradient of `ll` amd Fisher-Information.

## Usage

``` r
gsl_odeiv2_fi(
  odeModel,
  experiments,
  p,
  abs.tol = 1e-06,
  rel.tol = 1e-05,
  initial.step.size = 0.001,
  method = 0,
  omit = 0,
  time.out = 1,
  num.steps = 0
)
```

## Arguments

- odeModel:

  the name of the ODE model to simulate (a shared library of the same
  name will be dynamically loaded and needs to be created first).
  Alternatively this can be the ode object created by
  [as_ode](https://icpm-kth.github.io/uqsa/reference/as_ode.md), with a
  shared library path attached to it.

- experiments:

  a list of `N` simulation experiments (time, parameters, initial value,
  events)

- p:

  a matrix of parameters with M columns

- abs.tol:

  absolute tolerance, real scalar

- rel.tol:

  relative tolerance, real scalar

- initial.step.size:

  initial value for the step size; the step size will adapt to a value
  that observes the tolerances, real scalar

- method:

  integration method (see
  [method](https://icpm-kth.github.io/uqsa/reference/method.md) and
  [name_method](https://icpm-kth.github.io/uqsa/reference/name_method.md)).

- omit:

  an integer that indicates how many of these to omit in this order:
  fisher information, gradient of the log-likelihood, log-likelihood

- time.out:

  in seconds (early rejection due to long simulation time). This can
  trigger at measurement times (outputTime).

- num.steps:

  maximum number of steps the integration method is permitted to do;
  early rejection. This condition can trigger at any point during the
  integration.

## Value

a list of the solution trajectories `y(t;p)` for all experiments (named
like the experiments), as well as the output functions

## Details

The model is always simulated using a shared library. The path to the
shared library can be passed in three different ways:

1.  Character vector: `odeModel <- c("AKAKR4","/tmp/path/AKAR4.so")`

2.  A comment: `comment(odeModel) <- "/tmp/path/AKAR4.so"`

3.  As part of the ode object:

    odeModel <- as_ode(m)`
    so_path(odeModel) <- "/tmp/path/AKAR4.so"

The shared library needs to be created first. Either with `R CMD SHLIB`,
[shlib](https://icpm-kth.github.io/uqsa/reference/shlib.md), or manually
on the system's command line (bash, zsh, etc.).

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
  print(o)
#>                 Model name : AKAR4
#>                     C file : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.c [2026-06-26 13:34:58.400742]
#>             shared library : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.so [2026-06-26 13:34:58.400742]
#>  Number of state variables : 2
#>       Number of parameters : 5
#>          Number of outputs : 1
#>          Conservation laws : 2
#>            Transformations : no
  y <- gsl_odeiv2_fi(o,ex,values(m$Parameter))
  print(length(y))
#> [1] 3
  print(names(y[[1]]))
#> [1] "cpuSeconds"        "numSteps"          "status"           
#> [4] "state"             "func"              "logLikelihood"    
#> [7] "gradLogLikelihood" "FisherInformation"
  par(mfrow=c(length(ex),1))
  for (i in seq_along(y)){
      plot(
          errors::as.errors(ex[[i]]$outputTimes),
          ex[[i]]$data,
          xlab="time",
          ylab=rownames(y[[i]]$data)[1],
          main=names(ex)[i],
          ylim=c(100,200)
      )
      lines(ex[[i]]$outputTimes,drop(y[[i]]$func),col='red')
  }

# }
```
