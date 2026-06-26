# Function that creates a function that reads the SBtab expression table and computes the expressions (which are functions of parameters) given the parameters passed in input

Function that creates a function that reads the SBtab expression table
and computes the expressions (which are functions of parameters) given
the parameters passed in input

## Usage

``` r
parameters_from_expressions_func(model.tab)
```

## Arguments

- model.tab:

  an SBtab file imported in R, as it is returned by the function
  SBtabVFGEN::sbtab_from_tsv

## Value

a function that given a vector of (sampling) parameters returns a vector
of parameters corresponding to the evaluations of the expressions in the
SBtab expression table (which are functions of the parameters)
