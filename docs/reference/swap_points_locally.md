# Swap the end-points of two Markov chains

This is a conditional swap, according to the rules of parallel
tempering. This function is only useful if the Markov chains have
returned to the global scope and one process will make the decision and
perform the swap, i.e.: the current state of each chain is locally
available.

## Usage

``` r
swap_points_locally(parMCMC)
```

## Arguments

- parMCMC:

  a list of Markov chain end points, each entry annotated with a
  temperature attribute `attr(parMCMC[[i]],"beta")`

## Value

a list with some members swapped
