# Multivariate Normal Distribution KLD

Like the [KLD](https://icpm-kth.github.io/uqsa/reference/KLD.md)
function, this function calculates KLD values, but for the specific case
of multivariate normal distributions \\\\mathcal{N}\_{A}\\ and
\\\\mathcal{N}\_{B}\\. The two distributions are specified using
\\\\mu\\ and \\\\Sigma\\ values (mean and covariance).

## Usage

``` r
exact_normal_kld(muA, SigmaA, muB, SigmaB)
```

## Arguments

- muA:

  mean of distribution A

- SigmaA:

  Covariance of distribution A

- muB:

  mean of distribution B

- SigmaB:

  Covariance of distribution B

## Value

the KLD value D(A\|B)
