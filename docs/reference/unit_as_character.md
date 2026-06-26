# converts a unit data.frame into a printable string

This is a crude function to make a printable representation of a unit
data.frame, with very explicit parentheses and exponents.

## Usage

``` r
unit_as_character(unit)
```

## Arguments

- unit:

  a data.frame created by unit.from.string()

## Value

a string representation of that data.frame purely for printing

## Details

This function is similar to unit.info.

## Examples

``` r
u <- unit.from.string("s^-1")
str <- unit_as_character(u)
print(str)
#> [1] "(second*10^(0))^(-1)"
unit.info("s^-1")
#> «s^-1» has been interpreted as the product of: 
#> (1 × second × 10^(0))^(-1)
#> [1] "(1 × second × 10^(0))^(-1)"
```
