# Find Integer

Given a ODE solver name (from the GSL solver module odeiv2), return an
integer offset `{0..10}`. This integer can be passed as the "method"
argument for all ODE simulator functions
([simulator.c](https://icpm-kth.github.io/uqsa/reference/simulator.c.md),
[simfi](https://icpm-kth.github.io/uqsa/reference/simfi.md))

## Usage

``` r
method(name)
```

## Arguments

- name:

  character scalar, name of the method

## Value

an integer that is acceptable to
[simfi](https://icpm-kth.github.io/uqsa/reference/simfi.md) and
[simulator.c](https://icpm-kth.github.io/uqsa/reference/simulator.c.md)
