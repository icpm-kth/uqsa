# Find a column that contains some kind of mathematic expression in a data.frame

Given a data.frame that should contain a column that assigns a math
expression to a name (in row names), this function returns a named
character vector with the expressions. The formula column shoul dbe
named "formula" (if it exists, only this column will be used). But some
other spellings will also work as fallback. As a fallback "value" is
acceptable as well, because it makes sense to say "the value of x is
'y/2+1'", even though it is not an atomic value (but an expression).

## Usage

``` r
formulae(df)
```

## Arguments

- df:

  a data.frame with a "formula" column

## Value

character vector with names taken from the row names of df

## Examples

``` r
df <- data.frame(formula=c("exp(x)","10^x","2*x + 3"),row.names=c("f1","f2",'f3'))
formulae(df)
#>        f1        f2        f3 
#>  "exp(x)"    "10^x" "2*x + 3" 
#> attr(,"unit")
#>  f1  f2  f3 
#> "1" "1" "1" 
df <- data.frame(value=c("exp(x)","10^x","2*x + 3"),row.names=c("f1","f2",'f3'))
formulae(df)
#>        f1        f2        f3 
#>  "exp(x)"    "10^x" "2*x + 3" 
#> attr(,"unit")
#>  f1  f2  f3 
#> "1" "1" "1" 
```
