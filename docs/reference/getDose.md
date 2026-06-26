# Get the values of the input for a series of dose response experiments

This function finds the vector of inputs that varies among a series of
experiments that are part of a dose response experiment

## Usage

``` r
getDose(experiments)
```

## Arguments

- experiments:

  list of experimental data from the same dose response experiment

## Value

vector of inputs (i.e. dose) that varies among the experiments provided
to the function. The name of the input is saved in "comment(dose)".
