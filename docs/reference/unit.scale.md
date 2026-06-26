# Unit scale from SI prefix

This function reads a prefix from a string and returns the exponent
(base-10) that this prefix represents.

## Usage

``` r
unit.scale(prefix)
```

## Arguments

- prefix:

  a string, e.g.: "M", "mega", "m", "milli", "µ", "micro", etc.

## Value

an integer that corresponds to the prefix, defaults to 0.

## Examples

``` r
print(unit.scale("M"))
#> Error in unit.scale("M"): could not find function "unit.scale"
```
