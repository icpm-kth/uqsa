# Returns a list of reaction coefficients

This function maps c("A","2 B") to c(1,2)

## Usage

``` r
onlyCoefficients(formulaList)
```

## Arguments

- formulaList:

  a list of character vectors, derived from the left or right side of a
  reaction formula:

## Value

a list of numeric coefficient vectors

## Details

This is a C function because it is much easier to write in C. C has the
strtod() function which expects a leading number and stops when the
numbers end. as.character() returns NA if the input contains any dirt.

The reaction formula is as tring like this: "A + 2 B \<=\> C", when
split at `<=>` and then later at `+`, we get the strings that must be
parsed: "A" and "2 B" for the left side and "C" for the right side. The
numbers are the stoichiometric constants, or coefficients.

## Examples

``` r
print(onlyCoefficients("12 A"))
#> [[1]]
#> [1] 12
#> 
```
