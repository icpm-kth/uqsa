# Find the Reaction a parameter appears in

Usually Parameters are a list of names and values, Reactions have
reaction kinetics, which can contain any number of parameters. To make
matters worse, a parameter may influence two or more reactions. The goal
is to decide on the parameter's conversion factor. We shall assume that
if the parameter appears in more than one reaction, it is the same kind
of context, and thus it must be converted exactly the sam eway. This
function finds the first reaction a given parameter appears in.

## Usage

``` r
parameterReaction(reactionKinetic, parNames)
```

## Arguments

- reactionKinetic:

  a character vector with all of the kinetic law expressions

- parNames:

  parameter names

## Value

list with two named components: fwd, bwd; each is a named vector of
integers, th enames are taken from parNames, the integer indicates a
reaction.

## Details

Note that the entire conversion may not work at all for complicated
kinetics, we use very simple rules here.
