# prints the simulation results

The results, if accidentally printed, are difficult to read. This
function prevents these accidental prints. It summarizes the results
instead.

## Usage

``` r
# S3 method for class 'simulation'
print(x, ...)
```

## Arguments

- x:

  simulation results

- ...:

  requirement of print generic, not used.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
ex <- experiments(m,o)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
s <- simfi(ex,o)
y <- s(values(m$Parameter))
print(y)
#> number of simulation experiments: 3
#>                                      400nM 
#> ------------------------------------------ 
#>               cpuSeconds: 0.000311
#>                 numSteps: 269
#>                   status: 0
#>                    state: 2, 225, 1 (dim)
#>                     func: 1, 225, 1 (dim)
#>            logLikelihood: -864.384
#>        gradLogLikelihood: 5, 1 (dim)
#>        FisherInformation: 5, 5, 1 (dim)
#> 
#>                                      100nM 
#> ------------------------------------------ 
#>               cpuSeconds: 0.000321
#>                 numSteps: 248
#>                   status: 0
#>                    state: 2, 225, 1 (dim)
#>                     func: 1, 225, 1 (dim)
#>            logLikelihood: -840.397
#>        gradLogLikelihood: 5, 1 (dim)
#>        FisherInformation: 5, 5, 1 (dim)
#> 
#>                                       25nM 
#> ------------------------------------------ 
#>               cpuSeconds: 0.0003
#>                 numSteps: 237
#>                   status: 0
#>                    state: 2, 225, 1 (dim)
#>                     func: 1, 225, 1 (dim)
#>            logLikelihood: -786.464
#>        gradLogLikelihood: 5, 1 (dim)
#>        FisherInformation: 5, 5, 1 (dim)
#> 
#> experiments:  400nM, 100nM, 25nM 
```
