# Print a summary about the ode

An ODE model was crteated by `as_ode` can be summarized here, including
information about the compiled version of the model.

## Usage

``` r
# S3 method for class 'ode'
print(x, ...)
```

## Arguments

- x:

  the ode

- ...:

  requirement of print generic, not used.

## Details

The ode model is for the most part a list of named vectors and matrices
which together encode the mathematical structure of the ode.

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
o <- as_ode(m)
print(o)
#>                 Model name : AKAR4
#> 
#> 
#>  Number of state variables : 2
#>       Number of parameters : 5
#>          Number of outputs : 1
#>          Conservation laws : 2
#>            Transformations : no
```
