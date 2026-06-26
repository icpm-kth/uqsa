# Returns only the names in a reaction formula

This is the companion function to onlyCoefficients. It returns the names
of reactants, without the stoichiometry.

## Usage

``` r
onlyNames(formulaList)
```

## Arguments

- formulaList:

  a list of strings like: "2 B" or "45 X"

## Value

a list of name vectors

## Examples

``` r
print(onlyNames("12 A"))
#> [[1]]
#> [1] " A"
#> 
```
