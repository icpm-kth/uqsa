# Load an example model for this package

This function finds the path to an example model, given by name. In the
SBtab format, model and data travel together (in different tables, but
the same documents).

## Usage

``` r
uqsa_example(
  modelName = NULL,
  full.names = TRUE,
  pattern = "[.]tsv$",
  f = NULL
)
```

## Arguments

- modelName:

  name of model, e.g.: "AKAR4", "AKAP79", "CaMKII"; if empty, this
  function lists all available examples.

- full.names:

  return full paths to files - defaults to TRUE

- pattern:

  pattern to find specific files; if `NULL`, this function returns the
  directory of the example

- f:

  file ending, search for file endings in `f`, alternative to `pattern`

## Value

The location of the examples in the current environment if called with
no arguments, the paths to the model files if a modelName was provided
or the full path to the example if the file pattern *pattern* is unset

## Details

By default this function returns the names of the tsv files belonging to
the named model. If no modelName is provided it returns possible names
(contents of the top-level example directory).

## Examples

``` r
uqsa_example()
#> [1] "AKAP79"    "AKAR4"     "CaMKII"    "CaMKIIs"   "README.md" "Spike"    
uqsa_example("AKAR4",full.names=FALSE)
#> [1] "100nM.tsv"       "25nM.tsv"        "400nM.tsv"       "Compound.tsv"   
#> [5] "Experiments.tsv" "Output.tsv"      "Parameter.tsv"   "Reaction.tsv"   
uqsa_example("AKAP79",f='R',full.names=FALSE)
#> character(0)
uqsa_example("AKAP79",pat="^run.*R$")
#> character(0)
```
