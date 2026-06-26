# NATURAL LOG parameter mapping used by the MCMC module

This map is used by the simulator to transform sampling variables into
ODE-model porameters.

## Usage

``` r
logParMap(parMCMC)
```

## Arguments

- parMCMC:

  the sampling variables (numeric vector)

## Examples

``` r
p <- c(-1,0,1)
parMap <- logParMap
print(parMap(p))
#> [1] 0.3678794 1.0000000 2.7182818
```
