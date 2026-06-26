# Writes code to file and compiles

This function accepts an ode model, or cme model, generates code,
compiles it to a shared library, and returns a changed object. possibly
changed by the user. It writes the contents to a c file named
'modelName_gvf.c'. This file is compiled to './modelName.so' using
normal command line tools, not `R CMD SHLIB`

## Usage

``` r
write_and_compile(M)
```

## Arguments

- M:

  ode or cme Model for which code is generated and written to a file

## Value

a copy of `o` with file paths added to it

## Details

This entire function can be replaced with a call to
[`cat()`](https://rdrr.io/r/base/cat.html) and then compiling the
written file in the system's shell.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- write_and_compile(as_ode(m))
print(o)
#>                 Model name : AKAR4
#>                     C file : /tmp/RtmpiY0iWI/RtmpQMQrCy/AKAR4_cb4825f036ad8.c [2026-06-26 13:35:48.657317]
#>             shared library : /tmp/RtmpiY0iWI/RtmpQMQrCy/AKAR4_cb4825f036ad8.so [2026-06-26 13:35:48.657317]
#>  Number of state variables : 2
#>       Number of parameters : 5
#>          Number of outputs : 1
#>          Conservation laws : 2
#>            Transformations : no
```
