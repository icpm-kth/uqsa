# Construct Code

Interpret the first argument and generate code in the specified language
for the model type.

## Usage

``` r
generate_code(Model, language = "C", LV = 602214076)
```

## Arguments

- Model:

  either CME or ODE model

- language:

  either C or R

- LV:

  Avogadro's constant multiplied by the system's volume in litres, only
  used for CME models

## Value

a character vector with the code

## Details

Whenever the model `Model` is of type `"cme"`, the LV parameter is used
to determine the actual number of molecules in the system. Otherwise it
is ignored.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
C <- generate_code(o)
cat(head(C),sep="\n")
#> #include <stdlib.h>
#> #include <math.h>
#> #include <string.h>
#> #include <gsl/gsl_errno.h>
#> #include <gsl/gsl_odeiv2.h>
#> #include <gsl/gsl_math.h>
```
