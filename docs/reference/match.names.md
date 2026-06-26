# Find the variable names in a formula

A reaction formula has reactants and products, separated by `<=>`, with
reactants on the left and products on the right (by convention). Each of
those is a plus separated list of reacting compounds and modifiers, with
optional coefficients, e.g.: `A + 2 B <=> AB2`

## Usage

``` r
match.names(chrv)
```

## Arguments

- chrv:

  a character vector as returned by parse.formula

## Value

coefficients, as a vector

## Details

Once the formula is split into left and right side, this function
determines the names. For the above example, this function returns
`c("A","B")` for the left side and `"AB2"` for the right side.

## Examples

``` r
lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.names)
#> $reactants
#> [1] "A"   "2*B"
#> 
#> $products
#> [1] "AB2"
#> 
lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.names)
#> $reactants
#> [1] "A"   "2*B"
#> 
#> $products
#> [1] "AB2"
#> 
```
