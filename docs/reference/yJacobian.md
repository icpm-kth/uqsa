# Jacobian of string-math

Given a named character array of math expressions and a vector of
independent variables, this function calculates the Jacobian matrix of
the math expressions with respect to the variables.

## Usage

``` r
yJacobian(f, x)
```

## Arguments

- f:

  a character vector of length n

- x:

  a character vector of length m

## Value

a character matrix (n×m) with derivatives `df[i]/dx[j]`

## Examples

``` r
f <- c("2*x*y","exp(-k*x)")
x <- c("x","y")
J <- yJacobian(f,x)
print(J)
#>      [,1]           [,2] 
#> [1,] "2*y"          "2*x"
#> [2,] "-k*exp(-k*x)" "0"  
```
