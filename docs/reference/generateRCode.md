# Write R code

This function expects a list of character vectors, as returned by
`as_ode`. This list describes an ODE model (initial values, default
parameters, transformation events, output functions). This function uses
this information, calculates Jacobians via Ryacas and returns a
character vector with R source code for the deSolve package.

## Usage

``` r
generateRCode(odeModel)
```

## Arguments

- odeModel:

  a list of named character vectors with math expressions (which work as
  R code)

## Value

a character vector with the generated code, one vector-element is one
line of code.

## Details

The value can be written to a file:
`cat(generateRCode(odeModel),sep="\n",file=...)`. This file can be
sourced later.

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
generateRCode(as_ode(m))
#>   [1] "# ODE vector field: y' = f(t,y;p)"                                                                                               
#>   [2] "AKAR4_vf <- function(t, state, parameters) {"                                                                                    
#>   [3] "##\tconstants"                                                                                                                   
#>   [4] "##\tparameter values"                                                                                                            
#>   [5] "\tkf_C_AKAR4 = parameters[1];"                                                                                                   
#>   [6] "\tkb_C_AKAR4 = parameters[2];"                                                                                                   
#>   [7] "\tkcat_AKARp = parameters[3];"                                                                                                   
#>   [8] "\tAKAR4_C_ConservedConst = parameters[4];"                                                                                       
#>   [9] "\tAKAR4_ConservedConst = parameters[5];"                                                                                         
#>  [10] "##\tstate variables"                                                                                                             
#>  [11] "\tAKAR4p = state[1];"                                                                                                            
#>  [12] "\tC = state[2];"                                                                                                                 
#>  [13] "##\texpressions"                                                                                                                 
#>  [14] "\tAKAR4_C = AKAR4_C_ConservedConst - (+C);"                                                                                      
#>  [15] "\tAKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);"                                                                                 
#>  [16] "\treaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;"                                                                         
#>  [17] "\treaction_2 = kcat_AKARp*AKAR4_C;"                                                                                              
#>  [18] "\tf_ <- numeric(2)"                                                                                                              
#>  [19] "\tf_[1] <- +reaction_2"                                                                                                          
#>  [20] "\tf_[2] <- -reaction_1+reaction_2"                                                                                               
#>  [21] "\tnames(f_) <- c(\"AKAR4p\", \"C\")"                                                                                             
#>  [22] "\treturn(list(f_))"                                                                                                              
#>  [23] "}"                                                                                                                               
#>  [24] ""                                                                                                                                
#>  [25] "# ODE Jacobian: df(t,y;p)/dy"                                                                                                    
#>  [26] "AKAR4_jac <- function(t, state, parameters) {"                                                                                   
#>  [27] "##\tconstants"                                                                                                                   
#>  [28] "##\tparameter values"                                                                                                            
#>  [29] "\tkf_C_AKAR4 = parameters[1];"                                                                                                   
#>  [30] "\tkb_C_AKAR4 = parameters[2];"                                                                                                   
#>  [31] "\tkcat_AKARp = parameters[3];"                                                                                                   
#>  [32] "\tAKAR4_C_ConservedConst = parameters[4];"                                                                                       
#>  [33] "\tAKAR4_ConservedConst = parameters[5];"                                                                                         
#>  [34] "##\tstate variables"                                                                                                             
#>  [35] "\tAKAR4p = state[1];"                                                                                                            
#>  [36] "\tC = state[2];"                                                                                                                 
#>  [37] "##\texpressions"                                                                                                                 
#>  [38] "\tAKAR4_C = AKAR4_C_ConservedConst - (+C);"                                                                                      
#>  [39] "\tAKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);"                                                                                 
#>  [40] "\treaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;"                                                                         
#>  [41] "\treaction_2 = kcat_AKARp*AKAR4_C;"                                                                                              
#>  [42] "\tjac_ <- matrix(0.0,2,2)"                                                                                                       
#>  [43] "\tjac_[2,1] <- kf_C_AKAR4*C"                                                                                                     
#>  [44] "\tjac_[1,2] <- -kcat_AKARp"                                                                                                      
#>  [45] "\tjac_[2,2] <- -(kcat_AKARp+kb_C_AKAR4+kf_C_AKAR4*C+kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4p-C)))"                               
#>  [46] "## no names for return value"                                                                                                    
#>  [47] "\treturn(jac_)"                                                                                                                  
#>  [48] "}"                                                                                                                               
#>  [49] ""                                                                                                                                
#>  [50] "# ODE parameter Jacobian: df(t,y;p)/dp"                                                                                          
#>  [51] "AKAR4_jacp <- function(t, state, parameters) {"                                                                                  
#>  [52] "##\tconstants"                                                                                                                   
#>  [53] "##\tparameter values"                                                                                                            
#>  [54] "\tkf_C_AKAR4 = parameters[1];"                                                                                                   
#>  [55] "\tkb_C_AKAR4 = parameters[2];"                                                                                                   
#>  [56] "\tkcat_AKARp = parameters[3];"                                                                                                   
#>  [57] "\tAKAR4_C_ConservedConst = parameters[4];"                                                                                       
#>  [58] "\tAKAR4_ConservedConst = parameters[5];"                                                                                         
#>  [59] "##\tstate variables"                                                                                                             
#>  [60] "\tAKAR4p = state[1];"                                                                                                            
#>  [61] "\tC = state[2];"                                                                                                                 
#>  [62] "##\texpressions"                                                                                                                 
#>  [63] "\tAKAR4_C = AKAR4_C_ConservedConst - (+C);"                                                                                      
#>  [64] "\tAKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);"                                                                                 
#>  [65] "\treaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;"                                                                         
#>  [66] "\treaction_2 = kcat_AKARp*AKAR4_C;"                                                                                              
#>  [67] "\tjacp_ <- matrix(0.0,2,5)"                                                                                                      
#>  [68] "\tjacp_[2,1] <- -C*(AKAR4_ConservedConst-(AKAR4p-C))"                                                                            
#>  [69] "\tjacp_[2,2] <- AKAR4_C_ConservedConst-C"                                                                                        
#>  [70] "\tjacp_[1,3] <- AKAR4_C_ConservedConst-C"                                                                                        
#>  [71] "\tjacp_[2,3] <- AKAR4_C_ConservedConst-C"                                                                                        
#>  [72] "\tjacp_[1,4] <- kcat_AKARp"                                                                                                      
#>  [73] "\tjacp_[2,4] <- kb_C_AKAR4+kcat_AKARp"                                                                                           
#>  [74] "\tjacp_[2,5] <- -kf_C_AKAR4*C"                                                                                                   
#>  [75] "## no names for return value"                                                                                                    
#>  [76] "\treturn(jacp_)"                                                                                                                 
#>  [77] "}"                                                                                                                               
#>  [78] ""                                                                                                                                
#>  [79] "# Output Function (Observables)"                                                                                                 
#>  [80] "AKAR4_func <- function(t, state, parameters) {"                                                                                  
#>  [81] "##\tconstants"                                                                                                                   
#>  [82] "##\tparameter values"                                                                                                            
#>  [83] "\tkf_C_AKAR4 = parameters[1];"                                                                                                   
#>  [84] "\tkb_C_AKAR4 = parameters[2];"                                                                                                   
#>  [85] "\tkcat_AKARp = parameters[3];"                                                                                                   
#>  [86] "\tAKAR4_C_ConservedConst = parameters[4];"                                                                                       
#>  [87] "\tAKAR4_ConservedConst = parameters[5];"                                                                                         
#>  [88] "##\tstate variables"                                                                                                             
#>  [89] "\tAKAR4p = state[1];"                                                                                                            
#>  [90] "\tC = state[2];"                                                                                                                 
#>  [91] "##\texpressions"                                                                                                                 
#>  [92] "\tAKAR4_C = AKAR4_C_ConservedConst - (+C);"                                                                                      
#>  [93] "\tAKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);"                                                                                 
#>  [94] "\treaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;"                                                                         
#>  [95] "\treaction_2 = kcat_AKARp*AKAR4_C;"                                                                                              
#>  [96] "\tfunc_ <- numeric(1)"                                                                                                           
#>  [97] "\tfunc_[1] <- 108 + 380*AKAR4p"                                                                                                  
#>  [98] "\tnames(func_) <- c(\"AKAR4pOUT\")"                                                                                              
#>  [99] "\treturn(func_)"                                                                                                                 
#> [100] "}"                                                                                                                               
#> [101] ""                                                                                                                                
#> [102] "AKAR4_default <- function(t) {"                                                                                                  
#> [103] "##\tconstants"                                                                                                                   
#> [104] "\tparameters <- numeric(5)"                                                                                                      
#> [105] "\tparameters[1] <- 0.018"                                                                                                        
#> [106] "\tparameters[2] <- 0.106"                                                                                                        
#> [107] "\tparameters[3] <- 10.2"                                                                                                         
#> [108] "\tparameters[5] <- 0.2"                                                                                                          
#> [109] "\tnames(parameters) <- c(\"kf_C_AKAR4\", \"kb_C_AKAR4\", \"kcat_AKARp\", \"AKAR4_C_ConservedConst\", \"AKAR4_ConservedConst\")"  
#> [110] "\treturn(parameters)"                                                                                                            
#> [111] "}"                                                                                                                               
#> [112] ""                                                                                                                                
#> [113] "AKAR4_init <- function(t, parameters) {"                                                                                         
#> [114] "##\tconstants"                                                                                                                   
#> [115] "##\tparameter values"                                                                                                            
#> [116] "\tkf_C_AKAR4 = parameters[1];"                                                                                                   
#> [117] "\tkb_C_AKAR4 = parameters[2];"                                                                                                   
#> [118] "\tkcat_AKARp = parameters[3];"                                                                                                   
#> [119] "\tAKAR4_C_ConservedConst = parameters[4];"                                                                                       
#> [120] "\tAKAR4_ConservedConst = parameters[5];"                                                                                         
#> [121] "\tstate <- numeric(2)"                                                                                                           
#> [122] "\tnames(state) <- c(\"AKAR4p\", \"C\")"                                                                                          
#> [123] "\treturn(state)"                                                                                                                 
#> [124] "}"                                                                                                                               
#> [125] ""                                                                                                                                
#> [126] "# a variable that collects all functions into one list:"                                                                         
#> [127] "model <- list(vf=AKAR4_vf, jac=AKAR4_jac, jacp=AKAR4_jacp, default=AKAR4_default, init=AKAR4_init, func=AKAR4_func,name='AKAR4')"
```
