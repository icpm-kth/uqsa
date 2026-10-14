# A binary operator on strings with units in them

The function calls the units utility and converts the string on the left
into the unit on the right, e.g.: "cm" %as% "inches", both units can
contain numbers. Any input that is accepted by the units utility is
acceptable, as long as it makes sense with the command line arguments:
`units -t "$originalUnit" "$targetUnit"`

## Usage

``` r
txtUnit %as% target
```

## Arguments

- txtUnit:

  a string with numeric values, including units, e.g. "3 cm"

- target:

  string, target unit, e.g. "m", must be scalar

## Value

a numeric value y: val*originalUnit = y*targetUnit, the target unit is
attached to the returned value, as a comment.

## Details

FreeBSD and MacOS have a buggy `units` command: it fails to signal
failure, so this function does some additional checks to make sure that
the conversion worked.

## Examples

``` r
  ## needs `units` utility (a system utility)
  if (nzchar(Sys.which("units"))){
      y <- "21 cm" %as% "inches"
      print(y)
      y <- "12 nmol/L" %as% "mol/L"
      print(y)
      y <- "12 mol/m^3" %as% "mmol/L"
      print(y)
  }
#>       21 cm 
#> "8.2677165" 
#> attr(,"unit")
#> [1] "inches"
#> 12 nmol/L 
#> "1.2e-08" 
#> attr(,"unit")
#> [1] "mol/L"
#> 12 mol/m^3 
#>       "12" 
#> attr(,"unit")
#> [1] "mmol/L"
```
