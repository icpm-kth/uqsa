# Plot time series simulations with experimental data

This function plots simulations of time series experiments and plots
them against experimental data. The input in the provided experiments
must differ only in one vector component.

## Usage

``` r
plotTimeSeries(simulations, experiments, show.plot = TRUE)
```

## Arguments

- simulations:

  list of simualtions as output from the simulator

- experiments:

  list of experiments

- show.plot:

  boolean variable. Set `show.plot=TRUE` to display plots when running
  the funcion, FALSE otherwise

## Value

list of plots with simulations and experimental data
