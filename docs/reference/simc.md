# This creates a closure that simulates the model, similar to simulator.c

This is a shorter alternative to simulator.c (C backend).

## Usage

``` r
simc(experiments, modelName, parMap = identity, method = 0)
```

## Arguments

- experiments:

  a list of experiments to simulate: inital values, inputs, time
  vectors, initial times

- modelName:

  a string (with optional comment indicating an .so file) which points
  out the model to simulate

- parMap:

  the model will be called with parMap(parABC); so any parameter
  transformation can happen there.

- parABC:

  the parameters for the model, subject to change by parMap.

## Value

a closure that returns the model's output for a given parameter vector,
and approximate sensitivity matrices, for each state variable, function,
time-point, and parameter vector.

## Details

It returns a closure around: - experiments, - the model, and - parameter
mapping

The returned function depends only on parABC (the sampling parameters).
The simulation will be done suing the rgsl backend.

This version of the function does not use the parallel package at all
and cannot add noise to the simulations.

## Examples

``` r
 #  model.sbtab <- SBtabVFGEN::sbtab_from_tsv(dir(pattern="[.]tsv$"))
 #  experiments <- SBtabVFGEN::sbtab.data(model.sbtab)
 #  parABC <- SBtabVFGEN::sbtab.quantity(model.sbtab$Parameter)

 #  modelName <- checkModel("<insert_model_name>_gvf.c")
 #  simulate <- simc(experiments, modelName,  parABC)
 #  yf <- simulate(parABC)
```
