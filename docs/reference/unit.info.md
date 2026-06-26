# Prints an interpretation string of a unit

given a string describing a unit of measurement, this function prints
the interpretation on screen, rather than returning it as a data.frame

## Usage

``` r
unit.info(unit.str, unit = unit.from.string(unit.str))
```

## Arguments

- unit.str:

  unit string

- unit:

  optionally, the data.frame that describes the unit

## Examples

``` r
print(unit.info("km/h",unit.from.string("km/h")))
#> «km/h» has been interpreted as the product of: 
#> (1 × metre × 10^(3))^(1)
#> (60 × second × 10^(0))^(-1)
#> [1] "(1 × metre × 10^(3))^(1)"    "(60 × second × 10^(0))^(-1)"
```
