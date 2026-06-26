# Attempt to find multiplicative reaction rate coefficients

This function assumes that the kinetic Law is mass action kinetics. This
function helps in converting the units of an ODE into units that work in
stochastic simulations. Converting the units of a general formula
(Michaelis Menten, Hill kinetics, etc.) is difficult and dubious in
stochastic simulations.

## Usage

``` r
parameter.from.kinetic.law(kineticLaw, tab)
```

## Arguments

- kineticLaw:

  a string with a mathematical formula

- tab:

  an SBtab document, as returned by
  [`SBtabVFGEN::sbtab_from_tsv()`](https://rdrr.io/pkg/SBtabVFGEN/man/sbtab_from_tsv.html)

## Value

a parameter with value, and unit as attribute

## Details

For these reasons, this functions assumes: `kf*A*B*[...]` with the first
word `kf` representing the reaction rate coefficient.

Given an SBtab document, this function finds the value (if any) and unit
of this coefficient. The coefficient itself can be defined as a fixed
constant, a parameter, or an algebraic expression in the document. The
most important attribute here is the unit.
