# Interpret a model as an ODE

This function accepts a list generated from a collection of TSV files
(or a similar format) and interprets the contents as an ordinary
differential equation (ODE).

## Usage

``` r
as_ode(m, cla = requireNamespace("pracma"))
```

## Arguments

- m:

  a list of data.frames, each corresponding to a TSV file or sheet in a
  spreadsheet.

- cla:

  a Boolean value indicating whether conservation law analysis should be
  performed.

## Value

a list that contains a summary of this model interpreted as an ODE,
crucially, the list contains the element `vf`, the right-hand-side
(vector field) of the ODE, this is the main result of this function.

## Details

The argument `m` can be obtained via
[`model_from_tsv()`](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md).
It has the components:

- `m$Constant`

- `m$Parameter`

- `m$Input`

- `m$Expression`

- `m$Compound`

- `m$Reaction`

- `m$Experiment`

There can be additional components describing measured data for this
model.

## Examples

``` r
f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
o <- as_ode(m)
#> Loading required namespace: pracma
print(names(o))
#>  [1] "vf"                    "const"                 "par"                  
#>  [4] "var"                   "exp"                   "func"                 
#>  [7] "stoichiometric_matrix" "conservationLaws"      "tf"                   
#> [10] "name"                  "c_path"                "c.date"               
#> [13] "so_path"               "so.date"              
print(o$vf)
#>                   AKAR4p                        C 
#>            "+reaction_2" "-reaction_1+reaction_2" 
```
