# Find the stoichiometry for a given parameter name

Given a list representation of stoichiometry (list of named integer
vectors), and a character vector of kinetic laws, this function returns
the correct stoichiometry for the given parameter.

## Usage

``` r
parameter_stoichiometry(p, reactants, products, kinetic.law)
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

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
r <- stoichiometry(strsplit(m$Reaction$reactants,"+",fixed=TRUE))
p <- stoichiometry(strsplit(m$Reaction$products ,"+",fixed=TRUE))
k <- kinetic_law_matrix(m$Reaction)
j <- parameter_stoichiometry("k5_1",r,p,k)
#> Error in parameter_stoichiometry("k5_1", r, p, k): could not find function "parameter_stoichiometry"
print(j)
#> Error in eval(expr, envir, enclos): object 'j' not found
```
