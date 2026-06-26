# Load an ODE model from a file

The input can be either a list of files, or one compressed archive file
like "model.tar.gz". This function determines the name of the model from
the name of the archive, or from the parent directory of uncompressed
tsv files. This auto-determined name is attached as a comment of the
return value. To override this choice, replace the comment. This
function exists for compatibility with the scripts in the RPN-derivative
repository and loads the zip and tar.gz files that
SBtabVFGEN::sbtab_to_vfgen writes to disk.

## Usage

``` r
loadODE(fileList)
```

## Arguments

- fileList:

  list of file paths

## Value

a list of data.frames describing the model
