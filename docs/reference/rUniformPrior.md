# Returns a random vector generator

The return value is a function that generates random vectors of the same
size as ll and ul from a uniform distribution within the limits defined
by ul and ll. The random vectors are returned as n rows of a matrix,
where n is the only argument of the returned function.

## Usage

``` r
rUniformPrior(ll, ul)
```

## Arguments

- ll:

  lower limit of the random variables (a vector)

- ul:

  upper limit of the random variables (same size vector as ll)

## Value

a uniform random vector generating function: runiform(n), where n is the
requested number of vectors (rows)

## Examples

``` r
rup<-rUniformPrior(ll=c(0,1,2),ul=c(1,2,3))
rup(12)
#>               [,1]     [,2]     [,3]
#>  [1,] 2.119115e-02 1.112956 2.087808
#>  [2,] 8.372698e-02 1.517037 2.193110
#>  [3,] 2.006098e-01 1.304639 2.643416
#>  [4,] 3.109665e-01 1.198347 2.945671
#>  [5,] 9.962794e-01 1.989742 2.097866
#>  [6,] 5.427049e-05 1.234473 2.671795
#>  [7,] 2.311328e-01 1.268707 2.998162
#>  [8,] 6.490393e-01 1.024905 2.459072
#>  [9,] 9.478592e-01 1.959817 2.000005
#> [10,] 2.069484e-01 1.900305 2.335639
#> [11,] 7.698646e-01 1.313786 2.661171
#> [12,] 8.494464e-01 1.138105 2.438771
```
