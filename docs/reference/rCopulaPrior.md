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
#>             [,1]       [,2]        [,3]
#> [1,]  1.01670956 -0.1359325 -0.05631105
#> [2,] -0.13593247  3.7574837 -0.80323033
#> [3,] -0.05631105 -0.8032303  9.79196158
print(D(10))
#>             [,1]       [,2]       [,3]
#>  [1,] -0.4615189  0.5237721  4.4344848
#>  [2,]  0.3114071 -1.0177968 -0.9907598
#>  [3,] -1.1048411  2.5744604 -1.0946116
#>  [4,] -1.8492980  3.5376979 -0.4918001
#>  [5,]  0.8443637  3.6610423  3.2282745
#>  [6,]  0.2632791 -1.2016008  1.6495253
#>  [7,] -2.0098311 -2.7904786 -1.8639339
#>  [8,] -1.1475635  0.6969277 -2.9795569
#>  [9,] -2.5214668  2.7463249 -1.5826019
#> [10,] -0.6532845  0.3971487 -0.7200967
```
