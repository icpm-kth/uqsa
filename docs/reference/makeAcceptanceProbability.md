# Ths function Creates an acceptanceProbability function

The returned closure needs only the sampling variables (parABC) as
inputand calculates a probability of accepting Markov chainmoves.

## Usage

``` r
makeAcceptanceProbability(
  experiments,
  modelName,
  getAcceptanceProbability,
  parMap = identity
)
```

## Arguments

- experiments:

  a list of experiments

- modelName:

  an annotated string, with the model name and model file as comment

- getAcceptanceProbability:

  an R function that mape the results of a simulation and experimental
  data to an acceptance probability

- parMap:

  an optional mapping between sampling parameters (parABC) and model
  parameters (e.g. rescaling,re-ordering).

## Value

a function that calculates probabilities given only parABC as input; it
implicitly uses all the argiments to this function.
