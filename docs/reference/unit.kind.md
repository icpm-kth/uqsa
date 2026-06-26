# find unit category

This function reads a unit, without SI prefix, and returns a string from
a smaller subset of unit kinds, similar to what is defined in SBML. This
normalizes the various ways to write the same unit: "meter", "m" and
"metre".

## Usage

``` r
unit.kind(kind)
```

## Arguments

- kind:

  the unnormalized string that humans use to write a unit (but without
  prefix)

## Value

normalized category name of the unit kind: litre, mole, metre, kilogram,
gram, ampere, candela, second, kelvin, hour, molarity, dimensionless.
defaults to dimensionless.

## Details

The unit kind of "m" is "metre", the kind of "g" is "gram".

## Examples

``` r
print(unit.kind("meter"))
#> Error in unit.kind("meter"): could not find function "unit.kind"
```
