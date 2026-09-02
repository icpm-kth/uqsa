# Interprets a character vector as names of logarithms

The values in `x` are possibly given in a logarithmic space. The
parameter `str_scale` gives this logarithmic scale (provided in a
language agnostic form), by a human. An empty string causes no
transformations. Similarly, providing no scale at all causes no
transformations.

## Usage

``` r
linear_scale(x, str_scale = attr(x, "scale"))
```

## Arguments

- x:

  values

- str_scale:

  character vector

## Value

a copy of `x`, transformed into linear space

## Details

The words in `str_scale` name a logarithm, e.g. "log10". Currently
understood scales:

- log10

- log2, ld

- ln, log

## Examples

``` r
x <- c(1,2,3,1,1,1)
attr(x,"scale") <- c("log10","log2","log","ln","ld","log5")
print(linear_scale(x))
#> [1] 10.000000  4.000000 20.085537  2.718282  2.000000  5.000000
#> attr(,"scale")
#> [1] "log10" "log2"  "log"   "ln"    "ld"    "log5" 
```
