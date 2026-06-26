# Interprets the provided model as a stochastic model

The chemical master equation can be simulated as a Markov jump process
(or continuous time Markov chain). One of the stochastic solver
algorithms is the Gillespie algorithm. This function return sa data
structure that can be used to generate code for the Gillespie solver in
this package.

## Usage

``` r
as_cme(m)
```

## Arguments

- m:

  list of data.frames, obtained via
  [`model_from_tsv()`](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md)

## Value

a list containing the interpreted model.

## Details

This function interprets the contionuous model `m` as a discrete state
model with molecule counts and propensities. For this reason, we need to
specify a volume for the simulations to take place in.

The model `m` is assumed to describe a reaction network, as a list of
data.frames (as retuned by
[model_from_tsv](https://icpm-kth.github.io/uqsa/reference/model_from_tsv.md)).
The systems biology information in the file is assumed to be
concentrations and rate coefficients, regardless of the interpretation
this function will derive from it. This is to make the model format of
the TSV file fairly uniform and independent of how we want to solve the
derived equations, be it ODE or CME.

Like the ode object, the returned object can also store the paths of
files we create for this model, with: c_path\<-, and so_path\<-

With the information provided with the rate coefficient units and a
volume, this function tries to convert everything to Gillespie rate
constants.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
stochasticModel <- as_cme(m)
print(stochasticModel)
#>                       Name : AKAR4
#> 
#> 
#>  Number of state variables : 4
#>       Number of parameters : 3
#>          Number of outputs : 1
#>        Number of constants : 0
```
