# Clear Yacas variables

The `reset` operation doesn't work in yacas, so this function wipes
every variable one by one.

## Usage

``` r
clear_yacas_environment()
```

## Value

list of variables cleared as a character array

## Examples

``` r
Ryacas::yac_str("y := 2*x")
#> [1] "2*x"
Ryacas::yac_str("restart")
#> [1] "restart"
print(Ryacas::yac_str("D(x) y^2"))
#> [1] "8*x"
clear_yacas_environment()
#> [1] "y"
print(Ryacas::yac_str("D(x) y^2"))
#> [1] "0"
```
