# Fetch an Attribute

This function differs from `rlang::%@%` in that it stops if the
attribute doesn't exist.

## Usage

``` r
x %@% a
```

## Arguments

- x:

  an R object (variable with attributes)

- a:

  the name of an attribute

## Value

the value of the attribute: `attr(x,a)`

## Details

This function tries to find a similarly named attribute disregarding
capitalization and using partila matching.

The only way from this function to return NULL is when `x` is null (the
object that supposedly has the attribute). For the purposes of this
function , NULL objects are treated as optional things, and thus their
attributes do not matter. Non-NULL objects that should have an
attribute, but don't are considered erroneous.

## Examples

``` r
x <- 1
attr(x,"unit") <- "m"
print(x %@% "unit")
#> [1] "m"
```
