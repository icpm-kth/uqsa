# NATURAL LOG parameter mapping, jacobian

This map is used by the simulator to transform sampling variables into
ODE-model porameters. As we often calculate sensitivites, we alos need
the jacobian of the map, due to the chain rule of differentiation.

## Usage

``` r
logParMapJac(parMCMC)
```

## Arguments

- parMCMC:

  the sampling variables (numeric vector)

## Examples

``` r
p <- c(-1,0,1)
parMap <- logParMap
parMpJ <- logParMapJac
print(parMap(p))
#> [1] 0.3678794 1.0000000 2.7182818
print(parMpJ(p))
#>           [,1] [,2]     [,3]
#> [1,] 0.3678794    0 0.000000
#> [2,] 0.0000000    1 0.000000
#> [3,] 0.0000000    0 2.718282
```
