# The mean value of an observable value given a parameter bin

The model parameters were binned with a histogram function. In each bin
one of the parameters is almost fixed (it varies much less than the
other parameters \[full range\]). This function returns the mean of the
observable for each bin, as a vector.

## Usage

``` r
observable.mean.in.bin(id, outputSample)
```

## Arguments

- id:

  is an integer vector that identifies the bin the parameter vector
  falls into to create the same row of the outputSample (the output
  stems from a model simulation with parameters). id has the same length
  as outputSample rows.

- outputSample:

  a matrix of output values, one output vector per row (different rows
  are results at different parameter values)

## Value

`M\[i,j\]` the mean of each `observable\[j\]` in `bin\[i\]`
