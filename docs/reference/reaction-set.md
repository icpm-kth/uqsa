# Add a Reaction to an ODE

This function operates on an ODE vector field. It makes it possible to
construct the right hand side one reaction at a time:
`reaction(vf,r,p) <- flux` will add a flux term to the ODE's right hand
side (vector field), adding the given `flux` (reaction kinetic) in the
right places, informed by the reactants and products.

## Usage

``` r
reaction(vf, r, p) <- value
```

## Arguments

- vf:

  a named character vector of the right length (number of state
  variables)

- r:

  a named vector of stoichiometric coefficients for the reactants

- p:

  a named vector of stoichiometric coefficients for the products

- value:

  a reaction rate (string)

## Value

updated vf

## Details

The reactants, e.g.: `c(A=2,B=1)` name which state variables are
affected by the flux (negatively, for reactants) The products, e.g.:
`c(C=1)` name which state variables are affected positively by the flux.

## Examples

``` r
if (FALSE) { # \dontrun{
Reaction <- "A + B <=> C"
r <- c(A=1,B=1)
p <- c(C=1)
vf <- c(A="",B="",C="") # empty
reaction(vf,r,p) <- "A*B-C"
print(vf)
} # }
```
