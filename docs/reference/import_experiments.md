# Reads the Data and Model Contained in an SBtab Document (tsv)

An SBtab Document is a set of tables that represent reactions,
compounds, parameters, and measured data that correspond to simulations
of the model under certain input conditions and initial values.

## Usage

``` r
import_experiments(modelName = NULL, SBtabDir)
```

## Arguments

- modelName:

  (string) the functions of the model have this prefix

- SBtabDir:

  (string) a local directory that contains tsv files (with SBtab
  content)

## Value

list of simulation experiments (and the data corresponding to that
simulation)

## Details

This function assumes that this information is stored in a series of tsv
files. The content is imported using the SBtabVFGEN package.

The data contents are reorganized into a list of simulation experiments
(initial values, measurement time points, etc.)
