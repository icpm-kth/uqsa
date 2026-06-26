# ABC acceptance of currently sampled values given old data (Prior)

The prior probability density model using copulas and vines is not
perfect, so values sampled from an imperfect prior estimate can be
checked against old data.

## Usage

``` r
checkFitWithPreviousExperiments(draws, objectiveFunction, delta)
```

## Arguments

- draws:

  matrix of sampled values (to be filtered).

- objectiveFunction:

  function that, given a (vectorial) parameter as input, simulates the
  model, and outputs the distance between experimental data and data
  simulated from the model with the parameter provided in input

- delta:

  the acceptance threshold.

## Value

a filtered subset of acceptable parameter draws

## Examples

``` r
if (FALSE) { # \dontrun{
  posterior <- ABCMCMC(...)
  passed <- checkFitWithPreviousExperiments(posterior$draws, objFunc, delta=0.5)
  posterior$draws <- passed
} # }
```
