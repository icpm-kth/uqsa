# Find the stoichiometry for a given parameter name

Given a list representation of stoichiometry (list of named integer
vectors), and a character vector of kinetic laws, this function returns
the correct stoichiometry for the given parameter.

## Usage

``` r
find_parameter_reaction(p, reactants, products, kinetic.law)
```

## Arguments

- p:

  parameter name

- reactants:

  stoichiometry of the reaction's left side

- products:

  stoichiometry of the reaction's right side

- kinetic.law:

  a character matrix of fluxes, column 1 for forward reaction, column2
  for backward reactions

## Value

the stoichiometry entry that belongs to the given parameter
