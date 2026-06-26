# This function calculates the conversion from moles to particle counts

Given some molar quantity with a unit (character vector), this function
calculates a conversion factor based on the SI prefixes used, as well as
what the exponent of LV (Avogadro's constant \* Volume) applies in this
case. This is an internal function.

## Usage

``` r
moleCountConversion(unit)
```

## Arguments

- unit:

  character vector, will be parsed for SI prefixes and to determine the
  mole component

## Value

a data.frame with an lvpower and factor component, with names like the
names of the unit vector

## Examples

``` r
f <- moleCountConversion("nM") # nanomoles/litre
#> Error in moleCountConversion("nM"): could not find function "moleCountConversion"
```
