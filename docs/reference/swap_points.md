# Swap the end-points of two Markov chains

This is a conditional swap, according to the rules of parallel
tempering. This function is only useful if the Markov chains have
returned to the global scope and one process will make the decision and
perform the swap.

## Usage

``` r
swap_points(parMCMC)
```

## Arguments

- parMCMC:

  a list of Markov chain end points

## Value

a list with some members swapped
