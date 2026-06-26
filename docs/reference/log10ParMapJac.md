# LOG10 parameter mapping, jacobian

This map is used by the simulator to transform sampling variables into
ODE-model porameters. As we often calculate sensitivites, we alos need
the jacobian of the map, due to the chain rule of differentiation.

## Usage

``` r
log10ParMapJac(parMCMC)
```

## Arguments

- parMCMC:

  the sampling variables (numeric vector)

## Examples

``` r
p <- c(-1,0,1)
parMap <- log10ParMap
parMpJ <- log10ParMapJac
print(parMap(p))
#> [1]  0.1  1.0 10.0
print(parMpJ(p))
#>           [,1]     [,2]     [,3]
#> [1,] 0.2302585 0.000000  0.00000
#> [2,] 0.0000000 2.302585  0.00000
#> [3,] 0.0000000 0.000000 23.02585
```
