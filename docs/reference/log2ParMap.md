# LOG2 parameter mapping used by the MCMC module

This map is used by the simulator to transform sampling variables into
ODE-model porameters.

## Usage

``` r
log2ParMap(parMCMC)
```

## Arguments

- parMCMC:

  the sampling variables (numeric vector)

## Examples

``` r
p <- c(-1,0,1)
parMap <- log2ParMap
print(parMap(p))
#> [1] 0.5 1.0 2.0
```
