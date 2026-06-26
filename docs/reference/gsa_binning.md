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
gsa_binning(parSample, outputSample, nBins = "Sturges")
```

## Arguments

- parSample:

  a matrix of parameter vectors (rows)

- outputSample:

  a matrix, with rows of outputs (row-index is the sample index)

- nBins:

  number of bins, if unset defaults to the default of the hist function

## Value

sensitivity `S[i,j]` of `output[i]` with respect to `parameter[j]`

## Examples

``` r
  rprior <- rNormalPrior(c(-1,0,1),c(1,2,3))
  X <- rprior(10000)
  colnames(X) <- LETTERS[seq(3)]
  Z <- exp(X[,1,drop=FALSE]+X[,2,drop=FALSE])
  colnames(Z) <- "alpha"
  GSA <- gsa_binning(X,Z)
  print(GSA)
#>               A         B            C
#> alpha 0.2412941 0.1565027 0.0003389213
  cat("global sensitivity of alpha with respect to B: ",GSA['alpha','B'],"\n")
#> global sensitivity of alpha with respect to B:  0.1565027 
```
