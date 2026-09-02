# replace_powers does string manipulation

This function takes a string argument with human readable math (e.g. R
code), and replaces the power operator `z^n` with C-compatible function
calls: `pow(x,n)`, it counts parentheses to determine the base and
exponent automatically.

## Usage

``` r
replace_powers(v)
```

## Arguments

- v:

  a character vector

## Value

a string where all occurrences of `^` have been replaced by function
calls like `pow()`

## Details

This functions assumes that gsl functions can be used, the GNU
Scientific Library includes powers of small integers. These functions
may be faster than always calling `pow` from `math.h`.

This is necessary because in C the `^` operator means something else
(exclusive bitwise xor for integers). No attempt will be made to cast
the numbers to `float` or `double`.

## Examples

``` r
print(replace_powers(c("2^3.1","10^-6","x^2","(1+(1+x))^(n-0.5)")))
#> [1] "pow(2, 3.1)"         "pow(10, -6)"         "gsl_pow_2(x)"       
#> [4] "pow(1+(1+x), n-0.5)"
```
