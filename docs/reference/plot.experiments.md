# plot function for experiments

This function uses plot.errors and the base plot functions like
[matplot](https://rdrr.io/r/graphics/matplot.html).

## Usage

``` r
# S3 method for class 'experiments'
plot(x, y, ...)
```

## Arguments

- x:

  experiment setup (a list)

- y:

  simulation results (a list)

- ...:

  forwarded to the more specific plot function
  [plot.errors](https://r-quantities.github.io/errors/reference/plot.errors.html)

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
plot(ex,y)



#> NULL
```
