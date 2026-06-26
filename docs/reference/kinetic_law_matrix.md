# Split Kinetic Law

This function performs a very simplified split of a kinetic law into a
forward part and a backward part, if it isn't pre-split in the file.

## Usage

``` r
kinetic_law_matrix(r)
```

## Arguments

- r:

  the reaction table (data.frame)

## Value

a character matrix with a forward and backward column

## Details

If the data.frame contains separate forward and backward rates these
will be returned instead. Instead of using this function,
`m$Reaction[,c("fwd","bwd")]` would accomplish a very similar thing.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAP79")) # not pre-split
print(colnames(m$Reaction))
#> [1] "kinetic.law"   "reactants"     "arw"           "products"     
#> [5] "is.reversible"
k <- kinetic_law_matrix(m$Reaction)
print(k)
#>     fwd                   bwd                  
#> r51 "k5_1*Rii_C"          "0"                  
#> r14 "k4_1*RiiP*C"         "k1_4*RiiP_C"        
#> r12 "k1_2*RiiP_C*cAMP"    "k2_1*RiiP_C_cAMP"   
#> r43 "k4_3*cAMP*RiiP"      "k3_4*RiiP_cAMP"     
#> r23 "k3_2*RiiP_cAMP*C"    "k2_3*RiiP_C_cAMP"   
#> r78 "k8_7*cAMP*Rii"       "k7_8*Rii_cAMP"      
#> r56 "k5_6*Rii_C*cAMP"     "k6_5*Rii_C_cAMP"    
#> r76 "k7_6*Rii_cAMP*C"     "k6_7*Rii_C_cAMP"    
#> r62 "k6_2*Rii_C_cAMP"     "0"                  
#> r58 "k8_5*Rii*C"          "k5_8*Rii_C"         
#> r44 "k4_4p*RiiP*CaN"      "k4p_4*RiiP_CaN"     
#> r33 "k3_3p*CaN*RiiP_cAMP" "k3p_3*RiiP_cAMP_CaN"
#> r48 "k4p_8*RiiP_CaN"      "0"                  
#> r37 "k3p_7*RiiP_cAMP_CaN" "0"                  
#> r1  "kf_C_AKAR4*C*AKAR4"  "kb_C_AKAR4*AKAR4_C" 
#> r2  "kcat_AKARp*AKAR4_C"  "0"                  
```
