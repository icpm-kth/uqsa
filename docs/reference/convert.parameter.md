# Convert ODE parameter to Gillespie parameter

ODE parameters usually have a different unit of measurement than the
parameters we need for stochastic simulators. ODEs have fluxes, which
are in multiples of `M/s` (M is mol/liter), same unit as the first
derivative of the state variables.

## Usage

``` r
convert.parameter(k, n = 0, LV = 602214076)
```

## Arguments

- k:

  the ODE reaction rate coefficient (mandatory)

- n:

  multiplicity of each reactant, if any (order \> 0); omit for
  zero-order

- LV:

  `L*V` -- product of *Avogadro's number* and *volume* defaults to
  6.02214076e+8

## Value

rescaled parameter for stochastic simulation with a comment of how to
re-scale it

## Details

The reaction rate coefficients of mass action kinetics, kf and kb have
units that are compatible with the flux units, depending on the order of
the reaction (the order is related to the reaction's stoichiometry).

## Examples

``` r
# reaction: "2 A + B -> C"
k <- 1.0
attr(k,'unit') <- "µM/s"
n <- c(2,1)
reactants <- c('A','B')
uqsa:::convert.parameter(k,n)
#> [1] 2.75739e-24
#> attr(,"unit")
#> [1] "µM/s"
#> attr(,"conversion")
#> [1] 2.75739e-24
```
