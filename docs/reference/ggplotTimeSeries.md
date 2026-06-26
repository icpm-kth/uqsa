# Plot time series simulations with experimental data

This function plots simulations of time series experiments and plots
them against experimental data. The input in the provided experiments
must differ only in one vector component.

## Usage

``` r
ggplotTimeSeries(
  simulations,
  experiments,
  nrow = NULL,
  ncol = NULL,
  plot.state = FALSE
)
```

## Arguments

- simulations:

  list of simualtions as output from the simulator

- experiments:

  list of experiments

- show.plot:

  boolean variable. Set show.plot=TRUE to display plots when running the
  funcion, FALSE otherwise

## Value

list of plots with simulations and experimental data
