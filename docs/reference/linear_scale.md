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

a copy of `x`, transformed in to linear space

## Details

The words in `str_scale` name a logarithm, e.g. "log10". Currently
understood scales:

- log10

- log2, ld

- ln, log

## Examples

``` r
if (FALSE) { # \dontrun{
x <- c(1,2,3)
attr(x,"scale") <- c("log10","log2","log")
print(linear_scale(x))
} # }
```
