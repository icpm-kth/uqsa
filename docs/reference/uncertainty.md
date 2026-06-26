# Find the uncertainty of values in a data.frame that is derived from a tsv file or similar

given a data.frame, this function will look for a column that contains
some kind of standard error and retrieve it. The returned numeric vector
will be named. This function is not intended for data, for data, the
[values](https://icpm-kth.github.io/uqsa/reference/values.md) function
will retrieve both the value and the standard error if it was specified.

## Usage

``` r
uncertainty(df)
```

## Arguments

- df:

  a data frame with a "value" column

## Value

a named numeric vector

## Details

This function is for the case that the table specifies a distribution
with a mean and an range (of some sort). The type of uncertainty found
will be attached as a comment to the returned value: "sd" standard
deviation for normal distribution, "se" standard error (for a normal
prior), and "range" for a uniform prior. Other priors are not recognised
yet.

The distinction between standard-error and standard-deviation doesn't
matter much here: either the value is some kind of mean and the
*uncertainty* is the standard-error or standard-deviation of the mean,
or it is a raw data-point (not averaged) and we know the standard
deviation (noise) of the device that measured it, then *uncertainty* is
the standard deviation of the noise distribution. In either case, the
value will be taken at face value and the uncertainty is used as sigma
in the default log-likelihood function.

Any entry of prior.distribution other than "uniform", will start a
search for some kind of standard deviation or standard error (or sigma).
As more priors are added, this function will look for the parameters of
those distributions.

This function makes many assumptions specifically that all variables in
the table have the same type of prior distribution (but not identically
distributed).
