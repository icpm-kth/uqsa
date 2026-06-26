# replace_powers

replace_powers

## Usage

``` r
replace_powers(v)
```

## Arguments

- v:

  a character vector

## Examples

``` r
print(replace_powers(c("2^3.1","10^-6","x^2")))
#> [1] "pow(2, 3.1)"  "pow(10, -6)"  "gsl_pow_2(x)"
```
