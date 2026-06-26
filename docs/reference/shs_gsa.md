# Outputs the global sensitivity scores SI and SIT, calculated by the Sobol-Homma-Saltelli method

M1, M2, and N are matrices prepared by
[`uqsa::shs_prior()`](https://icpm-kth.github.io/uqsa/reference/shs_prior.md).
The parameters (rows) from these matrices need to be simulated (using
any method), to obtain fM1, fM2 and fN.

## Usage

``` r
shs_gsa(fM1, fM2, fN, subtractMean = TRUE)
```

## Arguments

- fM1:

  output (f)unction values for M1, nSamples × nOuts

- fM2:

  output (f)unction values for M2, nSamples × nOuts

- fN:

  output (f)unction values for N, nSamples × nOuts × nPars

## Value

a list with sensitivity indices \$SI and total sensitivities \$SIT

## Details

These matrices are shaped similarly to M1, M2 and N respectively, but
now the parameters are replaced by the effects they have on a observable
of interest (the output). It can be the vector valued output at a
specific (single) time-point or a scalar output at different
time-points.

See Geir Halnes et al. (Halnes, Geir, et al. J. comp. neuroscience 27.3
(2009): 471.
