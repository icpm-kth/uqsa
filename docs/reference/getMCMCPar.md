# Selects MCMC scheme specific setup parameters

The MCMC scheme uses a transition kernel. This function returns the
parameters of that transition kernel. Better parameters make the Markov
chain perform better (i.e. lower auto-correlation).

## Usage

``` r
getMCMCPar(prePar, preDelta, p = 0.05, sfactor = 0.1, delta = 0.01, num = 1)
```

## Arguments

- prePar:

  a sample of parameters from pre-Calibration

- preDelta:

  distance values (scores) for those parameters

- p:

  fraction (top scoring) of sampled points to base Sigma on

- sfactor:

  scales Sigma up or down

- delta:

  ABC threshold

- num:

  number of different starting parameter vectors.

## Value

Sigma and startPar (matrix with `num` rows) as a list
