# LOG10 parameter mapping used by the MCMC module

This map is used by the simulator to transform sampling variables into
ODE-model parameters. This function is an example for the `parMap` slot
in sampling functions. A `parMap` function, like this one, must
transform an MCMC variable (vector) to a parameter vector that the
scientific model we simulate can work with.

## Usage

``` r
log10ParMap(parMCMC)
```

## Arguments

- parMCMC:

  the sampling variables (numeric vector)

## Value

a numeric vector intended for the simulator.
