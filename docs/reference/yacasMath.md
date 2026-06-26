# yacasMath converts math to Ryacas compatible math

Given a string like `"exp(2*x)"` this function returns a string that
yacas can process: `"Exp(2*x)"`

## Usage

``` r
yacasMath(v, reverse = FALSE)
```

## Arguments

- v:

  a chcarcter vector with math expressions

- reverse:

  do the reverse operation

## Value

a character vector with yacas math
