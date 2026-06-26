# Splits a formula into a left and right side

This function splits a reaction formulka apart into its parts, removing
whitespace on each side: `"A + 2 B <=> AB2"` will be split into a list
with two entries

## Usage

``` r
parse.formula(reactionFormula)
```

## Arguments

- string:

  - reactionFormula

## Value

a named list with a forward component and a backward component, each
entry contains a character vector

## Details

    list$reactants == c("A","2*B")
    list$products == c("AB2")
