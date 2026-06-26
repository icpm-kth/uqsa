# (for testing) A non-Copula sampling function as fallback

If the sample is not suited to infer a Copula (fitCopula fails), this
fuction uses much simpler rules to re-draw a new sample from an older
sample with some added noise.

## Usage

``` r
sampleWithNoise(X, sdf = 0.01, ...)
```

## Arguments

- X:

  an N×M matrix (N is the sample size), M is the number of variables
  (MCMC or ABC vars)

- sdf:

  factor to increase or decrease the standard deviation of the added
  noise

- ...:

  passed to [`base::sample.int()`](https://rdrr.io/r/base/sample.html)

- size:

  size of returned sample (passed to
  [`sample.int()`](https://rdrr.io/r/base/sample.html))

## Value

a matrix with `size` rows and `nrow(X)` columns.

## Details

This can be used during testing, in cases where the acceptance was very
low and we have to deal with a very low quality sample This function
should work like [`base::sample`](https://rdrr.io/r/base/sample.html),
but adds small noise. Missing values are always removed.
