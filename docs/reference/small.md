# Reduces the sample to its effective size

When plotting, we want to show only a few representative lines or points
derived from a sample. This function will determine the auto-correlation
length very roughly and use that number to thin out the sample to a
minimal size that still represents the original sample well.

## Usage

``` r
small(
  S,
  L = attr(S, "logLikelihood"),
  verbose = getOption("uqsa.verbose", default = interactive())
)
```

## Arguments

- S:

  an MCMC sample

- L:

  the log-likelihood values of S

- verbose:

  when TRUE the acf plot option is set to TRUE, and the found
  auto-correlation length is printed.

## Value

a smaller version of S

## Examples

``` r
S <- matrix(rnorm(300),100,3)
## the next line fakes auto-correlation:
attr(S,"logLikelihood") <- cos(seq(0,1,length.out=100)) + rnorm(100,sd=0.05)
print(dim(S))
#> [1] 100   3
print(dim(small(S)))
#> [1] 3 3
```
