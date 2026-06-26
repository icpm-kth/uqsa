# Generate C Code to solve a model stochastically

This function tries to generate code for a stochastic solver, but
assumming that the SBtab model is written in terms of concentrations and
rate coefficients.

## Usage

``` r
generateGillespieCode(sm, LV = 602214076)
```

## Arguments

- LV:

  Avogadro's Constant \* volume (in litres)

- sb:

  a stochatsic Gillespie model obtained via `makeGillespieModel`

## Value

character vector with code

## Details

The model is interpreted by `stochasticModel()`, parameters found in
reaction kinetics are converted to stochastic parameters. Input
parameters are not converted.

The default system size is 1 femtolitre.
