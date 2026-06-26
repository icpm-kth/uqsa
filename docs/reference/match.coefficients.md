# find the coefficients in a formula

A reaction formula has reactants and products, separated by \<=\>, with
reactants on the left and products on the right (by convention). Each of
those is a plus separated list of reacting compounds and modifiers, with
optional coefficients, e.g.: `A + 2 B <=> AB2`

## Usage

``` r
match.coefficients(chrv)
```

## Arguments

- chrv:

  a character vector as returned by parse.formula

## Value

coefficients, as a vector

## Details

Once the formula is split into left and right side, this function
determines the coefficients. For the above example, this function
returns `c(1,2)` for the left side and 1 for the right side.

## Examples

``` r
lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.coefficients)
#> $reactants
#> [1] 1 2
#> 
#> $products
#> [1] 1
#> 
```
