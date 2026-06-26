# Write C code

This function generates C code, expecting the information about the ODE
to come from the files that SBtabVFGEN::sbtab_to_vfgen() writes. It
exists for compatibility with the scripts in the RPN-derivative package.

## Usage

``` r
generateCodeFromFile(fileList)
```

## Arguments

- ode:

  a list of data.frames with the text representation of the ordinary
  differential equation

## Value

a character vector with the generated code, one element is one line.

## Details

sbtab_to_vfgen() also returns a variable. This variable (a list) has the
same information as the files and can be used directly (by a different
function).
