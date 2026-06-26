# SMMALA -- The default Extractor of the log-likelihood computed by the simfi solver

This function will only extract the log-likelihood value from the
solution via simfi. Simfi returns the likelihhod value based on the
assumption that the data provided in the list of experiments has a
Gaussian standard error. The value includes the normalising factor
1/sqrt(2*pi*sigma^2) for each measured value (the logarithm).

## Usage

``` r
simfiGaussianLogLikelihood(init = 0)
```

## Arguments

- init:

  a base value that will be added to the log-likelihood

## Value

a log-likelihood value for the set of experiments the solver was set up
with.
