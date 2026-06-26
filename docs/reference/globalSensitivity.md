# Global Sensitivity Analysis

This function performs a binning based estimation of the global
sensitivity of a model's output with respect to the model's parameters.
The output can be a prediction of the model's behaviour in a scenario of
interest (parameters, input, intial values, boundary conditions,
scheduled events etc.). The output models a potentially measurable value
(the "observable"). The sample-rows and the output rows must correspond
(they must be from the same model simulation).

## Usage

``` r
globalSensitivity(parSample, outputSample, nBins = "Sturges")
```

## Arguments

- parSample:

  a matrix of parameter vectors (rows)

- outputSample:

  a matrix, with rows of outputs (row-index is the sample index)

- nBins:

  number of bins, if unset defaults to the default of the hist function

## Value

sensitivity Si,j of outputi with respect to parameterj
