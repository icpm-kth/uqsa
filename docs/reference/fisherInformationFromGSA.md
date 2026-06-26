# Calculate Global Fisher Information

Given a sample, this performs global sensitivity analysis, and then
squares the sensitivity.

## Usage

``` r
fisherInformationFromGSA(Sample, yf = NULL, E)
```

## Arguments

- Sample:

  an MCMC sample, or ABC sample

- yf:

  the simulations of the sample

- E:

  experiments list
