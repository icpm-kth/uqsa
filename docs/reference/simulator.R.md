# This creates a closure that simulates the model

This is a shorter alternative to the runModel function (R, deSolve
backend).

## Usage

``` r
simulator.R(experiments, model, parMap = identity)
```

## Arguments

- experiments:

  a list of experiments to simulate: inital values, inputs, time
  vectors, initial times

- parMap:

  the model will be called with parMap(parABC); so any parameter
  transformation can happen there.

- modelName:

  a string (with optional comment indicating an .so file) which points
  out the model to simulate

- parABC:

  the parameters for the model, subject to change by parMap.

## Value

a closure that returns the model's output for a given parameter vector

## Details

It returns a closure around: - experiments, - the model, and - parameter
mapping

The returned function depends only on parABC (the sampling parameters).

## Examples

``` r
   model.sbtab <- SBtabVFGEN::sbtab_from_tsv(dir(pattern="[.]tsv$"))
#> Error in file(con, "r"): invalid 'description' argument
   experiments <- SBtabVFGEN::sbtab.data(model.sbtab)
#> Error in is.factor(x): object 'model.sbtab' not found
   parABC <- SBtabVFGEN::sbtab.quantity(model.sbtab$Parameter)
#> Error: 'sbtab.quantity' is not an exported object from 'namespace:SBtabVFGEN'

   source("<model name>.R") # this defines the `model` variable
#> Warning: cannot open file '<model name>.R': No such file or directory
#> Error in file(filename, "r", encoding = encoding): cannot open the connection
   simulate <- simulator.R(experiments, model,  parABC)
#> Error in simulator.R(experiments, model, parABC): object 'experiments' not found
   yf <- sim(parABC)
#> Error in sim(parABC): could not find function "sim"
```
