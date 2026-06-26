# Copula Formulation for Uniform Prior Distributions

Covers the (simpler) special case where the `prior(x)` is iid uniform.
The return value has the same structure as the value of
[`fitCopula()`](https://icpm-kth.github.io/uqsa/reference/fitCopula.md).

## Usage

``` r
makeIndepCopula(ll, ul)
```

## Arguments

- ll:

  `ll[i]` is the lower limit of random variable `x[i]`

- ul:

  upper limit, analogous to ll.

## Value

list with: copula, U, Z, and Y entries.
