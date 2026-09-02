# Prints an interpretation string of a unit

The unit object is a tagged data frame, with these columns:

- multiplier

- kind

- scale

- exponent

## Usage

``` r
# S3 method for class 'unit_of_measurement'
print(x, ...)
```

## Arguments

- x:

  an object of type 'unit_of_measurement'

- ...:

  required by the generic print function.

## Value

called for the side-effect; no value.

## Details

The interpretation is the same as in SBML units. This function also
prints an inferred unit id: a string that has no special characters in
it and can be used in places where such characters are not allowed (e.g.
SBML unit `id` attribute).

The original string that a unit was derived from is attached to the unit
object as a [comment](https://rdrr.io/r/base/comment.html).

Units are produced by the function
[unit.from.string](https://icpm-kth.github.io/uqsa/reference/unit.from.string.md).

## Examples

``` r
lapply(lapply(c("km/h","s^-2","1/s"),unit.from.string),print)
#> «km/h» has been interpreted as
#> the product of: 
#>                       km_per_h
#> ==============================
#> (1 × metre × 10^(3))^(1)
#> (60 × second × 10^(0))^(-1)
#> «s^-2» has been interpreted as
#> the product of: 
#>   s_to_the_power_of_2_inverted
#> ==============================
#> (1 × second × 10^(0))^(-2)
#> «1/s» has been interpreted as
#> the product of: 
#>                     one_over_s
#> ==============================
#> (1 × dimensionless × 10^(0))^(1)
#> (1 × second × 10^(0))^(-1)
#> [[1]]
#> «km/h» has been interpreted as
#> the product of: 
#>                       km_per_h
#> ==============================
#> (1 × metre × 10^(3))^(1)
#> (60 × second × 10^(0))^(-1)
#> 
#> [[2]]
#> «s^-2» has been interpreted as
#> the product of: 
#>   s_to_the_power_of_2_inverted
#> ==============================
#> (1 × second × 10^(0))^(-2)
#> 
#> [[3]]
#> «1/s» has been interpreted as
#> the product of: 
#>                     one_over_s
#> ==============================
#> (1 × dimensionless × 10^(0))^(1)
#> (1 × second × 10^(0))^(-1)
#> 
```
