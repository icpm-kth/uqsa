# Initial count of reacting Compounds

This function returns the molecule count apart from LV. This returned
number must be multiplied by Avogadro's constant and volume.

## Usage

``` r
initialCount(m)
```

## Arguments

- m:

  the model data.frames obtained from
  [model_from_tsv](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md)

## Value

numeric vector or character vector, depending on how the initial
concentration was provided

## Details

The multiplication with LV isn't done here, because this allows the user
to change the volume of the system in the C-file, without re-generating
it from a model. The same c-file can be re-used that way.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
i <- initialCount(m) # as multiples of LV
#> Error in initialCount(m): could not find function "initialCount"
l <- i>0
#> Error in eval(expr, envir, enclos): object 'i' not found
cat(
  sprintf("initial count of %s is %i (%g * LV)",rownames(m$Compound)[l],round(i[l]*6e8),i[l]),
  sep="\n"
)
#> Error in eval(expr, envir, enclos): object 'l' not found
```
