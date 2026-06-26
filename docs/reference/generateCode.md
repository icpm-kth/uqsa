# Write C code

This function expects a list of character vectors, as returned by
[`as_ode()`](https://icpm-kth.github.io/uqsa/reference/as_ode.md). This
list describes an ODE model (initial values, default parameters,
transformation events, output functions). This function uses this
information, calculates Jacobians via Ryacas and returns a character
vector with C source code for the solvers in the GNU Scientific Library
(GSL).

## Usage

``` r
generateCode(odeModel)
```

## Arguments

- odeModel:

  a list that represents an ODE

## Value

a character vector with the generated code, one vector-element is one
line of code.

## Details

The value can be written to a file:
`cat(generateCode(odeModel),sep="\n",file=...)`. This file can be
compiled into a shared library.

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
C <- generateCode(as_ode(m))
#> Warning: This function will start a background yacas process via Ryacas.
#> There is currently no working way to reset/restart that process.
#>  It is therefore not advisable to generate the code for two different models in the same R session.
#> The definitions for the two models will be mixed up. 
cat(head(C,12),sep='\n')
#> #include <stdlib.h>
#> #include <math.h>
#> #include <string.h>
#> #include <gsl/gsl_errno.h>
#> #include <gsl/gsl_odeiv2.h>
#> #include <gsl/gsl_math.h>
#> 
#> /* Enums will be used for indexing purposes.   */
#> enum stateVariable { _AKAR4p, _C, numStateVar };
#> enum param { _kf_C_AKAR4, _kb_C_AKAR4, _kcat_AKARp, _AKAR4_C_ConservedConst, _AKAR4_ConservedConst, numParam };
#> enum func { _AKAR4pOUT, numFunc };
#> 
```
