# Standard Error Matrix from an errors object

If a matrix has an `errors` attribute, it is usually a vector. THis
function returns the values of this attribute as a matrix (it preserves
the dimensions of the host matrix).

## Usage

``` r
standard_error_matrix(M)
```

## Arguments

- M:

  a matrix with errors (uncertainties)

## Value

A matrix similar to E, with standard error values

## Examples

``` r
M <- matrix(seq(12),3,4,dimnames=list(letters[seq(3)],LETTERS[seq(4)]))
errors::errors(M) <- abs(M*0.1 + 0.1)
E <- standard_error_matrix(M)
print(E)
#>     A   B   C   D
#> a 0.2 0.5 0.8 1.1
#> b 0.3 0.6 0.9 1.2
#> c 0.4 0.7 1.0 1.3
```
