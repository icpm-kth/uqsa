# Weighted Sum of Bin-specific variances

This function calculates the variance sum of a vector valued observable.

## Usage

``` r
sum.of.bin.variance(hst, binMeans, totalMean)
```

## Arguments

- hst:

  the histogram of the parameter sample

- binMeans:

  the means of the observable within each bin (rows of means)

- totalMean:

  the mean of the observable over the entire sample (vector)

## Value

The weighted sum of square differences between the binMean and the
totalMean
