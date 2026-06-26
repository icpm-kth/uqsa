# High Level Metropolis function

This function uses default assumption everywhere and returns a function
that will sample from the given model. This funciton will generate code,
compile the code, create an ODE solver for it, infer the sampling space
from the scale of the parameters, create all necessary functions to move
in parameter space (gradients of likelihood and prior), as well as
Fisher Information functions.

## Usage

``` r
high_level_metropolis(
  m,
  o = as_ode(m, cla = FALSE),
  ex = experiments(m, o),
  x = values(m$Parameter),
  beta = 1
)
```

## Arguments

- m:

  the model's TSV representation read via `model_from_tsv`

- o:

  (optional) ode representation of `m`

- ex:

  experiments of `m`, with simulation instructions for `o`.

- x:

  initial point of the markov chain, pre in itialized to have the right
  attributes.

- beta:

  for parallel tempering, the log-likelihood will have a factor of
  `beta` applied to it

## Value

`smmala` a function of three arguments: p0, N, eps; where p0 is the
starting point, N is the desired sample-size, and eps is the step size.
This function has an attribute called "init", with a pre-initialized
starting point.

## Examples

``` r
# \donttest{
m <- model_from_tsv(uqsa_example("AKAP79"))
rwm <- high_level_metropolis(m) # "random walk", metropolis algorithm
#> The parameters are given in log10-scale, so the simulator will do the reverse transformation: 10^p.
p <- rwm %@% "init"             # a valid starting point
N <- 200
smallSample <- rwm(rwm %@% "init",N,1e-6)
plot(
  smallSample %@% "logLikelihood",
  type="l",
  main=sprintf("%i iterations",N),
  xlab="iterations",
  ylab="log-likelihood"
)

# }
```
