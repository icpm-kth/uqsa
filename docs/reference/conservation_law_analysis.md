# Reduce the size of the system

Given a stoichiometric matrix, this function performs model reduction
via linear algebra operations, with
[pracma::null](https://rdrr.io/pkg/pracma/man/nullspace.html).

## Usage

``` r
conservation_law_analysis(nu, iv, verbose = FALSE)
```

## Arguments

- nu:

  stoichiometric matrix

- iv:

  initial values

- verbose:

  if `TRUE`, this function will print the conservation laws on screen

## Value

a list of conservation laws

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
nu <- stoichiometric_matrix(m)
CL <- conservation_law_analysis(nu,values(m$Compound))
print(names(CL))
#> [1] "Eliminates"   "value"        "Constant"     "ConstantName" "Formula"     
print(CL[,c('value','Formula')])
#>         value                              Formula
#> AKAR4_C   0.0        AKAR4_C_ConservedConst - (+C)
#> AKAR4     0.2 AKAR4_ConservedConst - (+AKAR4p-1*C)
```
