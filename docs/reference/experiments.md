# Extract Measured Data and Simulation Experiment Instructions

This function accepts the model obtained via `model_fromt_tsv` or a
similar function. It finds the data tables for this model (if any are
present), and finds the simulation instructions to reproduce these data
sets using the model.

## Usage

``` r
experiments(m, o = NULL)
```

## Arguments

- m:

  the model (with data), as obtained via
  [`model_from_tsv()`](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md),
  or similar.

- o:

  the ode derived from `m`, only necessary if the experiments need to
  take conservation laws into account

## Value

a list of simulation instructions

## Details

This function requires that the files the model is stored as contains
measurements (data) that can be interpreted fairly easily. Each data
file needs columns that are named like the observable quantities listed
in the Output table.

If the data is very indirectly related to the model, then we don't
interpret the data files themselves and the user needs to write a
specialized likelihood function to relate the raw data in the files with
something that the model does. In such cases, don't use this function.

The simulation experiments returned here, include model input
parameters. Whenever conservation law analysis is perfomed, the
conserved constants are set as input parameters, because the conserved
amount can differ between experiments. For this reason the Experiment
table is interpreted differently in the presence of conservation laws.
Otherwise (no conservation laws), the `o` parameter can be omitted.

The instructions must be organised in a table called Experiment(s).

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
o <- as_ode(m)
ex <- experiments(m,o)
print(names(ex))
#> [1] "400nM" "100nM" "25nM" 
print(ex[[1]]$input)
#> AKAR4_C_ConservedConst   AKAR4_ConservedConst 
#>                    0.4                   -0.2 
```
