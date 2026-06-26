# propensity creates a propensity formula

given the custom math expressions needed to calculate a propensity, the
propensity coefficient and the kinetic law of the reaction, this
function makes a string that can be used with GillespieSSA2.

## Usage

``` r
propensity(conv.coeff, kinetic.law, rExpressions)
```

## Arguments

- conv.coeff:

  propensity conversion coefficient:
  `conv.coeff*kinetic.law = propensity function`

- kinetic.law:

  the kinetic law of this reaction (as used with ODEs)

- rExpressions:

  named math expressions that appear in the kinetic.law of this reaction

## Value

a string representation of the propensity function

## Details

The propensity coefficient translates between
