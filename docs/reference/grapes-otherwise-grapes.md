# This function can be used to specify default values

When attributes are missing, the
[`base::attr()`](https://rdrr.io/r/base/attr.html) function returns
NULL. In those cases this function can be used to find an alternative
value in one expression: `attr(x,"dim") %otherwise% length(x)`

## Usage

``` r
a %otherwise% b
```

## Arguments

- a:

  value to check for NULL

- b:

  value to substitute

## Value

a, or b if a is NULL

## Examples

``` r
x <- numeric(10)
l <- dim(x) %otherwise% c(length(x),1)
## example with attributes:
attr(x,"logLikelihood") <- -980
## elsewhere:
logLF <- attr(x,"logLikelihood") %otherwise% -Inf
```
