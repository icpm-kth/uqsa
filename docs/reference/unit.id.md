# Converts a unit to a string that works as an identifier

Some formats require a name for a unit definition. This functions
creates a name from a unit, converting math/symbols to text. The
returned value should work as an SBML unit id.

## Usage

``` r
unit.id(unit.str, prnt = FALSE)
```

## Arguments

- unit.str:

  the original string representastion of that unit

- prnt:

  logical switch: if TRUE, the name will be printed.

## Value

unit.id string

## Examples

``` r
print(unit.id("s^9"))
#> [1] "s_to_the_power_of_9"
print(unit.id("cm^2"))
#> [1] "cm_square"
print(unit.id("1/s"))
#> [1] "one_over_s"
```
