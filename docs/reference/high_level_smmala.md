# High Level SMMALA function

This function uses default assumption everywhere and returns a function
that will sample from the given model. This function will generate code,
compile the code, create an ODE solver for it, infer the sampling space
from the scale of the parameters, create all necessary functions to move
in parameter space (gradients of likelihood and prior), as well as
Fisher Information functions.

## Usage

``` r
high_level_smmala(
  m,
  o = as_ode(m, cla = TRUE),
  ex = experiments(m, o),
  x = values(m$Parameter),
  verbose = interactive()
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

  initial point of the Markov chain, pre in initialized to have the
  right attributes.

- verbose:

  prints extra messages when TRUE

## Value

`smmala` a function of three arguments: p0, N, eps; where p0 is the
starting point, N is the desired sample-size, and eps is the step size.
This function has an attribute called "init", with a pre-initialized
starting point.

## Examples

``` r
  m <- model_from_tsv(uqsa_example("AKAP79"))
  rwm <- high_level_smmala(m) # "random walk", metropolis algorithm
#> Warning: CONFLICT: This model seems to have event based transformations and conservation laws. These two concepts clash with one another if a compound is conserved, but also changed by scheduled events.
  p <- rwm %@% "init"          # a valid starting point
  N <- 100
  if (interactive()){
    smallSample <- rwm(rwm %@% "init",N,1e-4)
    plot(
      smallSample %@% "logLikelihood",
      type='l',
      main=sprintf("%i iterations",N),
      xlab="iterations",
      ylab="log-likelihood"
    )
  } else {
    smallSample <- rwm(rwm %@% "init",1,1e-4)
  }
```
