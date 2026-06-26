# Print a Summary about the CME model

This information printed on screen omits the deatls about the
interactions, only the lengths of the vectors included in the data
structure `CME`.

## Usage

``` r
# S3 method for class 'cme'
print(x, ...)
```

## Arguments

- x:

  a model created by
  [as_cme](https://icpm-kth.github.io/uqsa/reference/as_cme.md)

- ...:

  requirement of print generic, not used.

## Value

Nil

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
cmeModel <- as_cme(m)
print(cmeModel)
#>                       Name : AKAR4
#> 
#> 
#>  Number of state variables : 4
#>       Number of parameters : 3
#>          Number of outputs : 1
#>        Number of constants : 0
```
