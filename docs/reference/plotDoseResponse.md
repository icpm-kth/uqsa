# Plot dose response simulations with experimental data

This function plots simulations of one dose response experiment and
plots them against experimental data.

## Usage

``` r
plotDoseResponse(simulations, experiments, dose, show.plot = TRUE)
```

## Arguments

- simulations:

  list of simualtions as output from the simulator

- experiments:

  list of experimental data from the same dose response experiment

- dose:

  vector of dose values to plot on the x axis

- show.plot:

  boolean variable. Set `show.plot=TRUE` to display plots when running
  the funcion, `FALSE` otherwise

## Value

plot with simulations and experimental data
