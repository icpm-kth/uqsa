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
  h = 1e-04,
  N = 100,
  verbose = interactive()
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

- N:

  size of test-samples for acceptance rate estimate

- verbose:

  when TRUE, this function prints a progress bar, 'a:' reports current
  acceptance rate, and 'h:' reports the current step-size.

## Value

optimal step size

## Details

It will take 100 sample points repeatedly, until an acceptance of
`target_acceptance` is reached (defaults to 25%). The step-size is
decreased if acceptance is very low and increased when it is too high.

When verbose This function will do at most

## Examples

``` r
  opt <- options(mc.cores=2)
  m <- model_from_tsv(uqsa_example("AKAP79"))
  rwm <- high_level_metropolis(m) # "random walk", metropolis algorithm
  p <- rwm %@% "init"             # a valid starting point
  N <- 100
  if (interactive()){
    h <- tune_step_size(rwm,p)
    smallSample <- rwm(rwm %@% "init",N,h)
    print(h)
    plot(
      smallSample %@% "logLikelihood",
      type="l",
      main=sprintf("step size: %g",h),
      xlab="iterations",
      ylab="log-likelihood"
    )
  } else {
    h <- tune_step_size(rwm,p,N=20,iter.max=1)
  }
  options(opt)
```
