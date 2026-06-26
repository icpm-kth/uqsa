# SMMALA move

The Simiplified Manifold Metropolis Adjusted Langevin Algorithm uses a
move instriction that uses a Gaussian kernel that is shifted away from
the current point

## Usage

``` r
smmala_move(beta, parGiven, fisherInformationPrior, eps = 0.01)
```

## Arguments

- beta:

  inverse temperature (parallel tempering)

- parGiven:

  given point

## Value

SMMALA proposal point
