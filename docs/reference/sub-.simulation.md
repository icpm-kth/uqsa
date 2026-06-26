# Subset simulations with preserved class

The normal list subset operation would drop the class from the
simulation object (which is fine in theory). With this override, the
class is preserved.

## Usage

``` r
# S3 method for class 'simulation'
x[i, ...]
```

## Arguments

- x:

  an object with class "simulation"

- i:

  an index-set

- ...:

  passed on the list-\[ function

## Value

subset of experiments
