# Write the C code to a file

This function does not compile the code, it only writes it to a file in
a temporary location (tempdir). By default, the name of the file will
contain the hash of the entire code.

## Usage

``` r
write_c_code(C, model.name = comment(C), file = NULL)
```

## Arguments

- C:

  the code to write, as a character array.

- model.name:

  a string with no special characters, will be used in the file name

- file:

  override the default file name (based on hashing)

## Value

the path of the written file

## Details

If instead of a character vector, an ode or cme object is passed, this
function will generate code from it with default options.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
C <- generate_code(o)
c_path(o) <- write_c_code(C)
print(o)
#>                 Model name : AKAR4
#>                     C file : /tmp/RtmpiY0iWI/RtmpQMQrCy/adf9204aaf2748b8/AKAR4.c [2026-06-26 13:35:49.123137]
#> 
#>  Number of state variables : 2
#>       Number of parameters : 5
#>          Number of outputs : 1
#>          Conservation laws : 2
#>            Transformations : no
if (file.exists(c_path(o))) cat("c file exists.\n")
#> c file exists.
```
