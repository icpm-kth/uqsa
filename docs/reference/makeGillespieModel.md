# makeGillespieModel interprets the provided SBtab file as a stochastic model

The SBtab file is assumed to describe a reaction network. The systems
biology information is assumed to be concentrations and rate
coefficients.

## Usage

``` r
makeGillespieModel(m)
```

## Arguments

- m:

  list of data.frames, obtained via
  [`model_from_tsv()`](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md)

- LV:

  Avogadro's constant L multiplied by the system's volume V.

## Value

a list containing the interpreted model.

## Details

With the information provided with the rate coefficient units and a
volume, this function tries to convert everything to Gillespie rate
constants.
