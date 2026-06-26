# Add information about compiled code

Adds the path of the shared library (.so file) to the ODE model.

## Usage

``` r
so_path(o) <- value
```

## Arguments

- o:

  the ODE (list of named arrays and matrices), or CME model

- value:

  the path to the compiled model

## Value

modified o, with information about compiled code

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
print(o)
#>                 Model name : AKAR4
#>                     C file : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.c [2026-06-26 13:35:35.036577]
#>             shared library : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.so [2026-06-26 13:35:35.036577]
#>  Number of state variables : 2
#>       Number of parameters : 5
#>          Number of outputs : 1
#>          Conservation laws : 2
#>            Transformations : no
```
