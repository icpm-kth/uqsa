# Read Concise Error Notation

Convert a vector of strings of the form: c("1.2(3)E-4","1.2(3)E-2") to a
matrix with two rows:

1.  values,

2.  uncertainties.

## Usage

``` r
parse_concise(v, use.errors = requireNamespace("errors"), na = c(NA, NA))
```

## Arguments

- v:

  a character vector of numbers in concise error notation

- use.errors:

  if TRUE, the errors package will be used to retrun an object of type
  "errors" (from that package). Otherwise, the errors will be attached
  as an attribute (also called "errors" to be consistent with the errors
  package)

- na:

  a two element vector which will replace NA values, e.g. c(NA,NA);
  na=c(0,Inf) means *infinite uncertainty* for missing values

## Value

either a numeric object with class errors (with the same dimensions as
`v`), or a numeric matrix of values and uncertainties (2 rows),
dimensions of original object are lost

## Details

If the errors package is available, then an *errors* object is returned
instead (uncertainties are an attribute). In that case the dimensions of
`v` are preserved on output. You can override this choice using the
second argument `use.errors`.

Concise error notation means that a floating point number is followed by
an integer in parentheses which indicates the uncertainty of the last
digits of the value: \$\$1.2345(12) = 1.2345 \pm 0.0012\$\$.

If the errors package is installed, then it will be used to represent
the return value.

## Examples

``` r
x <- parse_concise(c("1.23(4)","0.51099895069(16)","1.25663706127(20)e-6","1.3±1.6","5;1"))
print(as.data.frame(x))
#>                    x
#> 1            1.23(4)
#> 2    0.5109989507(2)
#> 3 1.2566370613(2)e-6
#> 4               1(2)
#> 5               5(1)
```
