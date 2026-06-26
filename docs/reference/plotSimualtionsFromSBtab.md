# Simulate and plot Data and Simulation

This function imports model and experimental data saved as SBtab files,
it simualtes the model with the same initial conditions and input as in
the corresponding experiments, and it plots the simulations together
with the corresponding experimental data. It currently works with one
dimensional output.

## Usage

``` r
plotSimualtionsFromSBtab(
  SBtabDir,
  paramVal,
  plotDir = NULL,
  width = 15,
  heigth = 10
)
```

## Arguments

- SBtabDir:

  the directory that contains `.tsv` files (with SBtab content)

- paramVal:

  parameter vector with which the model has to be simulated

- plotDir:

  the directory where the plots will be saved as .pdf and .RData
  variables

- width:

  width of the plot window

- heigth:

  heigth of the plot window

## Value

a vector of R plot objects
