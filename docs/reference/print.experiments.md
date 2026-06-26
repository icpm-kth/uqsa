# prints the simulation experiments

The experiments, if accidentally printed, are difficult to read. This
function prevents these accidental prints. It summarizes the data and
simulation experiments instead.

## Usage

``` r
# S3 method for class 'experiments'
print(x, ...)
```

## Arguments

- x:

  simulation experiments with data

- ...:

  ignored.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
ex <- experiments(m,o)
print(ex)
#> number of simulation experiments: 3
#>                                      400nM 
#> ------------------------------------------ 
#>             measurements: 2 columns (data.frame)
#>                     data: 1, 225 (dim)
#>                    input: 2 (length)
#>              initialTime: -15
#>             initialState: 2 (length)
#>              outputTimes: 225 (length)
#>                   events: NULL (class), NULL (type)
#> 
#>                                      100nM 
#> ------------------------------------------ 
#>             measurements: 2 columns (data.frame)
#>                     data: 1, 225 (dim)
#>                    input: 2 (length)
#>              initialTime: -15
#>             initialState: 2 (length)
#>              outputTimes: 225 (length)
#>                   events: NULL (class), NULL (type)
#> 
#>                                       25nM 
#> ------------------------------------------ 
#>             measurements: 2 columns (data.frame)
#>                     data: 1, 225 (dim)
#>                    input: 2 (length)
#>              initialTime: -15
#>             initialState: 2 (length)
#>              outputTimes: 225 (length)
#>                   events: NULL (class), NULL (type)
#> 
#> experiments:  400nM, 100nM, 25nM 
```
