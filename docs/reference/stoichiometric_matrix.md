# The stoichiometric matrix of a reaction network

Given a model, described in tabular form (`m` is a list of data-frames).
The stoichiometric matrix is the linear map between the model's flux
vector and the ODE's right-gand-side vector field. If the flux vector is
`rr <- flux(t,x,p)`, which maps the state variables `x` and parameters
`p` to the reaction rate `rr` of each reaction. The stoichiometric
matrix `nu` (\\\nu\\), will map the reaction rates to the rate of change
of the state variables: dx/dt := nu %\*% flux(t,x,p).

## Usage

``` r
stoichiometric_matrix(m, compound.names = rownames(m$Compound))
```

## Arguments

- m:

  list of data frames with at least the 'Reaction' table, and the
  'Compound' table

- compound.names:

  all names of the reacting compounds

## Value

the stoichimetric matrix, with some additional attributes.

## Details

The matrix is usually sparse, but not extremely big. This function
attaches a sparse version of the same information as attributes to the
return-value, as two lists, for convenience.

## Examples

``` r
the_reaction <- "A + B <=> C"
m <- list(
    Reaction=data.frame(reactants=c("A+B"),products=c("C"))
)
nu <- stoichiometric_matrix(m,c("A","B","C"))
```
