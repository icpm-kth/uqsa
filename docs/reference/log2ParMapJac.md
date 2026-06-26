# LOG2 parameter mapping, jacobian

This map is used by the simulator to transform sampling variables into
ODE-model porameters. As we often calculate sensitivites, we alos need
the jacobian of the map, due to the chain rule of differentiation.

## Usage

``` r
log2ParMapJac(parMCMC)
```

## Arguments

- parMCMC:

  the sampling variables (numeric vector)

## Examples

``` r
p <- c(-1,0,1)
parMap <- log2ParMap
parMpJ <- log2ParMapJac
print(parMap(p))
#> [1] 0.5 1.0 2.0
print(parMpJ(p))
#>           [,1]      [,2]     [,3]
#> [1,] 0.3465736 0.0000000 0.000000
#> [2,] 0.0000000 0.6931472 0.000000
#> [3,] 0.0000000 0.0000000 1.386294
```
