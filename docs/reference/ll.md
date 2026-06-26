# Default Log-likelihood Function

Extracts the `logLikelihood` value from the simulations attribute of the
parMCMC argument, requires:

- parMCMC has simulations attribute

- simulations list includes logLikelihood values (omit\<3)

## Usage

``` r
ll(parMCMC)
```

## Arguments

- parMCMC:

  a numeric vector, with attributes for MCMC, specifically smmala

## Value

a scalar value: log(likelihood(data\|parMCMC))

## Details

This function will take the log-likelihood-values claculated by the ode
solver in this package, and return the sum of those values over all
experiments. The value the simulator returns is calculated with the
assumption of a normal distribution on measurement errors.

This function does *almost no work*, it merely sums up the values
calculated during simulation.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simulator.c(ex,o,omit=2) # not 3
p <- values(m$Parameter)
attr(p,"simulations") <- s(p)
print(ll(p))
#> [1] -2491.245
```
