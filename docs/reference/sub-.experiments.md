# Subset experiments with preserved class

The normal list subset operation would drop the class from the
experiment object (which is fine in theory). With this override, the
class is preserved.

## Usage

``` r
# S3 method for class 'experiments'
x[i, ...]
```

## Arguments

- x:

  an object with class "experiment"

- i:

  an index-set

- ...:

  passed on the list-\[ function

## Value

subset of experiments

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
x <- experiments(m)
class(x)
#> [1] "experiments"
class(x[seq(2)])
#> [1] "experiments"
```
