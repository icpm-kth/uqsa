# creates Objective functions from ingredients

the returned objective function has only one argument: the ABC variables
that shall be mapped to ODE-model parameters.

## Usage

``` r
makeObjective(experiments, simulate, distance = defaultDistance)
```

## Arguments

- experiments:

  a list of simulation experiments

- simulate:

  closure that simulates the model

- distance:

  a function that calculates ABC scores (distance between data and
  simulations)

## Value

an objective function

## Details

The user supplied distance function should accept three arguments:
distance(SIM, DATA, STDV), all three matrices. SIM is the model output
(simulation), DATA is the measured data, while STDV represents the
standard error of that measurement. All three have the same size: N×M,
where N is the number of observables (outputs), and M is th enumber of
measurement time-points (length of the time-series).

## Examples

``` r
# \donttest{
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  o <- as_ode(m)
  ex <- experiments(m,o)
  C <- generate_code(o)
  c_path(o) <- write_c_code(C)
  so_path(o) <- shlib(o)
  s <- simulator.c(ex,o)
  objFunc <- makeObjective(ex,s)
  print(objFunc(values(m$Parameter)))
#>             [,1]
#> 400nM 0.16439580
#> 100nM 0.13721115
#> 25nM  0.09171478
# }
```
