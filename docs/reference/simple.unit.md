# Simple unit from string

This function takes a simple, human readable unit (without '\*' or '/'),
from a string and returns a data.frame with the unit's meaning.

## Usage

``` r
simple.unit(u = NULL)
```

## Arguments

- u:

  a unit with no fractions or products

## Value

a data.frame with the unit's properties

## Details

In this context, a simple unit is just a prefix, a unit kind, and an
exponent, e.g. cm^2 A not-simple unit is: m/s, kg*m/s^2, kg*h
