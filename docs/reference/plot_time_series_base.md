# plot function for experiments

This function does not use ggplot2, only base plot functions like
lines() and arrows().

## Usage

``` r
plot_time_series_base(simulations, experiments, by = 1, ylimit = NULL)
```

## Arguments

- simulations:

  simulation results (a list)

- experiments:

  experiment setup (a list)

- by:

  skip this many sampled lines in between plotted lines

- ylimit:

  y-axis limits

## Value

plot object

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- write_and_compile(as_ode(m))
ex <- experiments(m,o)
s <- simulator.c(ex,o)
p0 <- values(m$Parameter)
y <- s(p0)
plot_time_series_base(y,ex)

#> NULL
```
