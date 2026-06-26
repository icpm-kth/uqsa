# Functions to construct and run the stochastic simulation using GillespieSSA2 package

This translates the Reaction network into the specific form required by
GillespieSSA2

## Usage

``` r
importReactionsSSA(model.tab, compile = TRUE)
```

## Arguments

- compile:

  boolean value: if TRUE (default value), the GillespieSSA2 reactions
  are compiled (recommended for faster simulations)

- model:

  the model, represented by a list of data.frames with SBtab content

## Value

a list of GillespieSSA2::reaction items (if compile==TRUE, the reactions
are compiled)
