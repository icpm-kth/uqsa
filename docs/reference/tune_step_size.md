# Find a good Step-Size for a given MCMC Algorithm

Given a closure `MCMC(p,N,eps)`, where `p` is the initial Markov-chain
position, `N` a sample-size, and `eps` a step-size, this function finds
a good value for eps.

## Usage

``` r
tune_step_size(
  MCMC,
  parMCMC = attr(MCMC, "init"),
  target_acceptance = 0.25,
  iter.max = 6,
  h = 1e-04
)
```

## Arguments

- MCMC:

  a Markov chain Monte Carlo closure (function)

- parMCMC:

  initial position of the Markov chain, has to be initialized with
  [mcmc_init](https://icpm-kth.github.io/uqsa/reference/mcmc_init.md).

- target_acceptance:

  a scalar value for the desired acceptance rate, some algorithms are
  most efficient with 20% to 30% acceptance, some work well with a very
  high acceptance.

- iter.max:

  maximum number of iterations until the function has to return.

- h:

  initial guess for the MCMC step size

## Value

optimal step size

## Details

It will take 100 sample points repeatedly, until an acceptance of
`target_acceptance` is reached (defaults to 25%). The step-size is
decreased if acceptance is very low and increased when it is too high.

This function will do at most

## Examples

``` r
# \donttest{
m <- model_from_tsv(uqsa_example("AKAP79"))
rwm <- high_level_metropolis(m) # "random walk", metropolis algorithm
#> The parameters are given in log10-scale, so the simulator will do the reverse transformation: 10^p.
p <- rwm %@% "init"             # a valid starting point
h <- tune_step_size(rwm,p)
#> acceptance rate: 0.14, step-size: 0.0001;
#> acceptance rate: 0.2, step-size: 4.77467e-05;
#> acceptance rate: 0.25, step-size: 3.72657e-05;
N <- 200
smallSample <- rwm(rwm %@% "init",N,h)
print(h)
#> [1] 3.726568e-05
plot(
  smallSample %@% "logLikelihood",
  type="l",
  main=sprintf("step size: %g",h),
  xlab="iterations",
  ylab="log-likelihood"
)

# }
```
