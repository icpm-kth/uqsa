# Updates Parameter Values

under valid ABC update conditions (successful simulation) the parameters
are updated to new values.

## Usage

``` r
parUpdate(
  objectiveFunction,
  curPar,
  canPar,
  curDelta,
  curPrior,
  delta,
  dprior,
  parAcceptable
)
```

## Arguments

- objectiveFunction:

  function that, given a (vectorial) parameter as input, simulates the
  model, and outputs the distance between experimental data and data
  simulated from the model with the parameter provided in input

- curPar:

  current parameter values (as ABC samples them)

- canPar:

  candidate parameter values (for MCMC)

- curDelta:

  current distance between data and simulation, if the MCMC chain has
  not yet reached any point where this is below the threshold (delta),
  this can be accepted as the new current state for the chain.

- curPrior:

  current Prior values given curPar

- delta:

  distance threshold for ABC

- dprior:

  prior probability density function

- parAcceptable:

  user-specified constraint function, must return a scalar Boolean

## Value

updated values for curPar, curDelta, and curPrior
