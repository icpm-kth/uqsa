# Compile C code to shared library

Calls `R CMD SHLIB` to create the model's shared library.

## Usage

``` r
shlib(file)
```

## Arguments

- file:

  the c file that is to be compiled, OR an ODE/CME object with a c.file
  defined and recorded in it.

## Value

the path of the created shared library

## Details

The first argument can be a raw character scalar with just the path of
the c code to be compiled, or alternatively an object that has this
information stored within it. The models returned by
[as_cme](https://icpm-kth.github.io/uqsa/reference/as_cme.md) and
[as_ode](https://icpm-kth.github.io/uqsa/reference/as_ode.md) can both
carry this information, attach it via c_path\<-.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
C <- generate_code(o)
c_path(o) <- write_c_code(C)
so_path(o) <- shlib(o)
print(o)
#>                 Model name : AKAR4
#>                     C file : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.c [2026-06-26 13:35:31.176118]
#>             shared library : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.so [2026-06-26 13:35:31.176118]
#>  Number of state variables : 2
#>       Number of parameters : 5
#>          Number of outputs : 1
#>          Conservation laws : 2
#>            Transformations : no
if (file.exists(so_path(o))) cat("shared library exists.\n")
#> shared library exists.
```
