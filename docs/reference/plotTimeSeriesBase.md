# plot function for experiments

This function does not use ggplot2, only base plot functions like
lines() and arrows().

## Usage

``` r
plotTimeSeriesBase(
  simulations,
  experiments,
  nmax = NULL,
  by = 1,
  ylimit = NULL
)
```

## Arguments

- simulations:

  simulation results (a list)

- experiments:

  experiment setup (a list)

- nmax:

  maximum index of lines to plot from sample, between 1 and size of
  sample

- by:

  skip this many sampled lines in between plotted lines

- ylimit:

  y-axis limits

## Value

plot object
