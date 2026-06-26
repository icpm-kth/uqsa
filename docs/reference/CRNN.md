# CRNN creates C code for a chemical reaction neural network

This function creates a very general ODE (c source code), that can be
compiled and simulated using the UQSA package.

## Usage

``` r
CRNN(numReactions, initialValues, funcValues, model.name = "CRNN")
```

## Arguments

- numReactions:

  the number of reversible mass action law reactions

- initialValues:

  named vector of initial values, names will be used as the neames of
  the reacting compounds.

- funcValues:

  named character vector, can be any valid C expression (one line) of
  the available state variables (the names can be used literally).

- model.name:

  the prefix of all created model functions: CRNN_vf, CRNN_jac, ...

## Value

a character vector suitable for writing to a file (.c)

## Details

Example: A + B \<=\> C numReactions: n \<- 1 initialValues: x \<-
c(A=2,B=3,C=0) funcValues: f \<- c("A+B","log(A)")

The above definition would create a CRNN inspired ODE, where A+B and
log(A) are treated as observable (measureable) values (functions of the
state variables).

Note: In addition to the state variables, the function values can also
reference the log-parameters of the model as `l[i]`, where `i` is a
0-based offset to the reaction `i`; the backward rate, is stored at
position `l[i+numRct]`.

## Examples

``` r
C <- CRNN(4,c(A=1,B=2,C=3),c(out="A+B+C"),model.name="testmodel")
cat(head(C),sep='\n')
#> #include <stdlib.h>
#> #include <math.h>
#> #include <string.h>
#> #include <gsl/gsl_errno.h>
#> #include <gsl/gsl_odeiv2.h>
#> #include <gsl/gsl_math.h>
```
