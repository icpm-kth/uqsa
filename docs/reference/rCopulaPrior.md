# rCopulaPrior returns a function that generates random values from the copula model

The returned function generates n random vectors, as rows of a matrix.

## Usage

``` r
rCopulaPrior(Copula)
```

## Arguments

- Copula:

  the return value of fitCopula()

## Value

a matrix of random values

## Examples

``` r
rprior <- rNormalPrior(c(-1,0,1),c(1,2,3))
C <- fitCopula(rprior(1000))
D <- rCopulaPrior(C)
print(cov(D(100)))
#>            [,1]       [,2]       [,3]
#> [1,]  0.8335288  0.3255351 -0.2465015
#> [2,]  0.3255351  4.0484247 -0.4511446
#> [3,] -0.2465015 -0.4511446  9.5039137
print(D(10))
#>             [,1]       [,2]       [,3]
#>  [1,] -1.1764805 -3.4237210  2.3613460
#>  [2,] -2.0822643 -1.2496824  0.9324375
#>  [3,] -1.8590347 -1.0796811 -1.0975961
#>  [4,] -1.9594211  1.6908897  0.9703570
#>  [5,]  0.2076257  4.5649173 -3.1332637
#>  [6,] -1.5208706  0.2375049  2.5788466
#>  [7,] -0.1620578  0.1140025 -2.7497282
#>  [8,] -0.1321217  2.3784682  2.8464257
#>  [9,] -0.9231011 -2.9415876  0.8929962
#> [10,] -2.2243977  1.5347527  0.8935569
```
