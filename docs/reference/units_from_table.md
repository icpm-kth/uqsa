# Get units from a data.frame column

Given a data.frame this funciton retrieves the strings in the unit
column named: unit, Unit, units (partial matching disregarding
capitalization).

## Usage

``` r
units_from_table(df, default = "1")
```

## Arguments

- df:

  a data.frame

- default:

  default value if no unit column exists

## Value

a character vector of units with names

## Details

The returned value uses the row names of the data.frame as names of the
character vector of units.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
u <- units_from_table(m$Compound)
```
