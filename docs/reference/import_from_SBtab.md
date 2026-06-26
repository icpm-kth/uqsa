# Import Systems Biology Models given in Tabular Form

This function uses SBtabVFGEN to read a series of tsv files, each
containing one systems biology table that together create a model
(Reactions, Parameters, etc.). The model is converted to a vfgen
compatible file. This file is processed by vfgen through a system call
to create source code for R (deSolve) and C (GSL solvers).

## Usage

``` r
import_from_SBtab(SBtabDir)
```

## Arguments

- SBtabDir:

  the directory that contains \`.tsv\` files (with SBtab content)

## Value

a model as a list of data.frames, one per tsv file, the Document title
is attached as a comment attribute: comment(model) = Document Title

## Details

This requires vfgen to be installed
(https://github.com/WarrenWeckesser/vfgen) - not an R package.

SBtab is a particular convention on how to structure the tables
(sbtab.net)

## Examples

``` r
 model <- import_from_SBtab("./model")
#> Warning: cannot open file './model/': No such file or directory
#> Error in file(con, "r"): cannot open the connection
 comment(model)
#> Error in comment(model): object 'model' not found
 source("model.R")
#> Warning: cannot open file 'model.R': No such file or directory
#> Error in file(filename, "r", encoding = encoding): cannot open the connection
```
