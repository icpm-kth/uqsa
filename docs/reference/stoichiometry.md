# This function returns a list of named stoichiometric vectors

Given an already split list of entries such as c("3 A","B"), this
function returns a numeric vector c(3,1) with names c("A","B").

## Usage

``` r
stoichiometry(formulaList)
```

## Arguments

- formulaList:

  reaction formulae, either as a pre-split list or character vector

## Value

named numeric vector of stoichiometric coefficients

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
r <- stoichiometry(m$Reaction$reactants)
print(head(r))
#> [[1]]
#> Rii_C 
#>     1 
#> 
#> [[2]]
#> RiiP    C 
#>    1    1 
#> 
#> [[3]]
#> RiiP_C   cAMP 
#>      1      1 
#> 
#> [[4]]
#> cAMP RiiP 
#>    1    1 
#> 
#> [[5]]
#> RiiP_cAMP         C 
#>         1         1 
#> 
#> [[6]]
#> cAMP  Rii 
#>    1    1 
#> 
print(tail(r))
#> [[1]]
#> RiiP  CaN 
#>    1    1 
#> 
#> [[2]]
#>       CaN RiiP_cAMP 
#>         1         1 
#> 
#> [[3]]
#> RiiP_CaN 
#>        1 
#> 
#> [[4]]
#> RiiP_cAMP_CaN 
#>             1 
#> 
#> [[5]]
#>     C AKAR4 
#>     1     1 
#> 
#> [[6]]
#> AKAR4_C 
#>       1 
#> 
```
