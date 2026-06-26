# checks whether a variable has the named attributes

checks whether a variable has the named attributes

## Usage

``` r
var %has% attrNames
```

## Arguments

- var:

  a variable to check for attributes

- attrNames:

  named attributes

## Value

TRUE if all attributes are present

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79"))
x <- values(m$Compound)
x %has% "unit"
#> [1] TRUE
print(x %@% "unit")
#>           Rii          cAMP          RiiP         Rii_C     RiiP_cAMP 
#>          "µM"          "µM"          "µM"          "µM"          "µM" 
#>        RiiP_C   RiiP_C_cAMP             C      Rii_cAMP    Rii_C_cAMP 
#>          "µM"          "µM"          "µM"          "µM"          "µM" 
#>           CaN      RiiP_CaN RiiP_cAMP_CaN         AKAR4       AKAR4_C 
#>          "µM"          "µM"          "µM"          "µM"          "µM" 
#>        AKAR4p 
#>          "µM" 
```
