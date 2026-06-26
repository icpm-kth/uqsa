# Convert to linear space

A number given in some logarithmic space can be transformed back to
linear space A call like `base(x) <- 10` means that x was provided in
decadic logarithm form. This will adjust `x` so that it is now in linear
space.

## Usage

``` r
base(x, i = seq_along(x)) <- value
```

## Arguments

- x:

  a numeric vector

- i:

  a subset of values in x, defaults to all values of x

- value:

  the base of the logarithm x was provided in

## Value

x will be changed to be in linear space

## Details

If `x` was provided in logarithmic space, then it is an exponent to the
given base (value).

## Examples

``` r
 x <- 2
 base(x) <- 10
 print(x)
#> [1] 100
```
