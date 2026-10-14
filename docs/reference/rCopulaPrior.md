# Returns a function that generates random values from the copula model

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
#> [1,]  1.1488690 -0.2509360 -0.2321524
#> [2,] -0.2509360  3.6417069 -0.6788716
#> [3,] -0.2321524 -0.6788716  8.9359975
print(D(10))
#>             [,1]       [,2]       [,3]
#>  [1,] -0.2462738 -0.4576989  2.6709661
#>  [2,] -1.8516084 -3.1100752 -0.8035401
#>  [3,] -1.5372687 -1.5350163 -1.5315393
#>  [4,] -2.1734932  3.1669473  4.2324968
#>  [5,]  0.7266450 -2.2483499  1.0641963
#>  [6,] -1.0876107  1.0389677 -1.5518637
#>  [7,] -0.9799121 -1.1294310 -2.0050341
#>  [8,]  0.1881597 -0.1782131  2.4098788
#>  [9,] -2.8288944  0.1447864  3.5936022
#> [10,]  1.0513657 -2.7416260 -2.9970855
```
