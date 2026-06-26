# %as% is a binary operator on strings with units in them

The function calls the units utility and converts the string on the left
into the unit on the right, e.g.: "cm" %as% "inches", both units can
contain numbers. Any input that is accepted by the units utility is
acceptable, as long as it makes sense with the command line arguments:
`units --strict --compact -1 "$originalUnit" "$targetUnit"`

## Usage

``` r
txtUnit %as% target
```

## Arguments

- txtUnit:

  a string with numeric values, including units, e.g. "3 cm", can be a
  character vector

- target:

  string, target unit, e.g. "m", must be scalar

## Value

a numeric value y: val*originalUnit = y*targetUnit, the target unit is
attached to the returned value, as a comment.

## Examples

``` r
if (FALSE) { # \dontrun{
  ## needs `unit` utility (system utility)
  y <- "21 cm" %as% "inches"
  y <- "12 nmol/L" %as% "mol/L"
  print(comment(y))
  y <- "12 mol/m^3" %as% "mmol/L"
} # }
```
