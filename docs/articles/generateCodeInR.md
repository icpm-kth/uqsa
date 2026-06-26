# Generate Model Code

``` r
library(uqsa)
```

It is possible to generate C code for the model in R. The currently
available computer algebra system for this task is the
[Ryacas](https://cran.r-project.org/package=Ryacas) package that
interfaces between R and [yacas](http://www.yacas.org/) – both need to
be installed.

The function `generateCode` accepts the return value of
[`SBtabVFGEN::sbtab_to_vfgen()`](https://rdrr.io/pkg/SBtabVFGEN/man/sbtab_to_vfgen.html)
(a list of ode properties) and returns a character vector of C code:

``` r
 f <- uqsa::uqsa_example("AKAP79")    # paths to the TSV files
sb <- SBtabVFGEN::sbtab_from_tsv(f)   # loads the TSV files
#> [tsv] file[1] «AKAP79_Compound.tsv» belongs to Document «AKAP79»
#>  I'll take this as the Model Name.
#> AKAP79_Compound.tsv  AKAP79_Experiments.tsv  AKAP79_Expression.tsv  AKAP79_Input.tsv  AKAP79_Output.tsv  AKAP79_Parameter.tsv  AKAP79_Reaction.tsv  X0uM_cAMPCaN_AKAP79_0_nM_cAMP.tsv  X0uM_cAMPCaN_only_0_nM_cAMP.tsv  X0uM_cAMPno_CaN_0_nM_cAMP.tsv  X1000nM_cAMPCaN_AKAP79_1_uM_cAMP.tsv  X1000nM_cAMPCaN_only_1_uM_cAMP.tsv  X1000nM_cAMPno_CaN_1_uM_cAMP.tsv  X100nM_cAMPCaN_AKAP79_100_nM_cAMP.tsv  X100nM_cAMPCaN_only_100_nM_cAMP.tsv  X100nM_cAMPno_CaN_100_nM_cAMP.tsv  X2000nM_cAMPCaN_AKAP79_2_uM_cAMP.tsv  X2000nM_cAMPCaN_only_2_uM_cAMP.tsv  X2000nM_cAMPno_CaN_2_uM_cAMP.tsv  X200nM_cAMPCaN_AKAP79_200_nM_cAMP.tsv  X200nM_cAMPCaN_only_200_nM_cAMP.tsv  X200nM_cAMPno_CaN_200_nM_cAMP.tsv  X500nM_cAMPCaN_AKAP79_500_nM_cAMP.tsv  X500nM_cAMPCaN_only_500_nM_cAMP.tsv  X500nM_cAMPno_CaN_500_nM_cAMP.tsv
```

Create the ODE files:

``` r
odeModel <- SBtabVFGEN::sbtab_to_vfgen(sb)
```

The first two steps above are unrelated to code generation. They just
find, and load the biological model. Then we convert it to an ordinary
differential equation (and create intermediate files, various formats).

``` r
C <- uqsa::generateCode(odeModel)
```

This returns the generated code as a character vector, it can be viewed
(and written to a file using the `cat` function). Here are a few lines
of the generated code:

``` r
cat(head(C,25),sep="\n")
#> #include <stdlib.h>
#> #include <math.h>
#> #include <string.h>
#> #include <gsl/gsl_errno.h>
#> #include <gsl/gsl_odeiv2.h>
#> #include <gsl/gsl_math.h>
#> 
#> /* Enums will be used for indexing purposes.   */
#> enum stateVariable { _RiiP, _RiiP_cAMP, _RiiP_C, _RiiP_C_cAMP, _C, _Rii_cAMP, _Rii_C_cAMP, _RiiP_CaN, _RiiP_cAMP_CaN, _AKAR4_C, _AKAR4p, numStateVar };
#> enum param { _kf_Rii_C__RiiP_C, _kf_RiiP_CxcAMP__RiiP_C_cAMP, _kf_RiiP_cAMPxC__RiiP_C_cAMP, _kb_RiiP_cAMPxC__RiiP_C_cAMP, _kb_RiiPXcAMP__RiiP_cAMP, _kf_RiiPXcAMP__RiiP_cAMP, _kf_RiiPxC__RiiP_C, _kb_RiiPxC__RiiP_C, _kf_cAMPxRii__Rii_cAMP, _kb_cAMPxRii__Rii_cAMP, _kf_Rii_CxcAMP__Rii_C_cAMP, _kb_Rii_CxcAMP__Rii_C_cAMP, _kf_RiixC__Rii_C, _kf_Rii_cAMPxC__Rii_C_cAMP, _kb_Rii_cAMPxC__Rii_C_cAMP, _kf_Rii_C_cAMP__RiiP_C_cAMP, _kb_RiixC__Rii_C, _AKAPoff_1, _AKAPoff_3, _AKAPon_1, _AKAPon_3, _kf_C_AKAR4, _kb_C_AKAR4, _kcat_AKARp, _kmOFF, _kmON, _KD_T, _b_AKAP, _AKAR4_ConservedConst, _CaN_ConservedConst, _Rii_C_ConservedConst, _cAMP_ConservedConst, _Rii_ConservedConst, numParam };
#> enum func { _AKAR4pOUT, numFunc };
#> 
#> /* The error codes indicate how many values a function returns.                             */
#> /* Each function expects the output buffer to be allocated with at least that many values   */
#> 
#> /* ODE vector field: y' = f(t,y;p)   */
#> int AKAP79_vf(double t, const double y_[], double *f_, void *par){
#>  double *p_=par;
#>  if (!y_ || !f_) return 11;
#> /*   constants   */
#> /*   parameter values   */
#>  double kf_Rii_C__RiiP_C = p_[_kf_Rii_C__RiiP_C];                    /* [  0] */
#>  double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP];       /* [  1] */
#>  double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  2] */
#>  double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP];       /* [  3] */
```

This command will print the code, use the `file=` argument of `cat` to
print it into a file.

## Code for deSolve

It is also possible to create R code, intended for use with
[deSolve](https://cran.r-project.org/package=deSolve). Combining all
previous function calls into one (writing to file immediately)

``` r
rFile <- paste0(comment(odeModel),".R")
cat(uqsa::generateRCode(odeModel),sep="\n",file=rFile)
source(rFile)

# example call:
p <- AKAP79_default(0.0)
y0 <- AKAP79_init(0.0,p)

print(y0)
#>          RiiP     RiiP_cAMP        RiiP_C   RiiP_C_cAMP             C 
#>             0             0             0             0             0 
#>      Rii_cAMP    Rii_C_cAMP      RiiP_CaN RiiP_cAMP_CaN       AKAR4_C 
#>             0             0             0             0             0 
#>        AKAR4p 
#>             0
```

A test-simulation:

``` r
if (require(deSolve)){
    t_ <- seq(0,50,length.out=128)
    y_ <- ode(y0,t_,AKAP79_vf,p)
    plot(t_,y_[,1],type='l',main="AKAP79",xlab="time")
}
#> Loading required package: deSolve
```

![](generateCodeInR_files/figure-html/simulation-1.png)
