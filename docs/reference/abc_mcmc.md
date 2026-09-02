# Performs and Approximate Bayesian Computation Sampling of Model Parameters

ABC replaces the need for an exact likelihood function and uses a
distance function instead: the distance between data and simulation.
This distance function is very similar to a likelihood function but
lacks a statistical justification. Nevertheless, this distance function,
like the likelihood function of a deterministic model, performs a
simulation of the scientific model, be it fully stochastic or a
stochastically embedded, but deterministic in its core.

## Usage

``` r
abc_mcmc(
  objectiveFunction,
  startPar,
  N,
  burnIn = ceiling(sqrt(N)),
  Sigma0 = cov(t(startPar)),
  dprior = NULL,
  deltaSpan = NULL,
  batchSize = 100 * NROW(Sigma0),
  parAcceptable = function(p) {
     all(is.finite(p))
 },
  verbose = TRUE
)
```

## Arguments

- objectiveFunction:

  function that, given a parameter matrix as input, simulates the model,
  and outputs the distance between experimental data and data simulated
  from the model with the parameter provided in input. This has to be a
  closure that contains the experimental data within itself. The closure
  must be vectorized over the columns of its matrix argument.

- startPar:

  starting values for the parameter vector, can (and should) be a matrix
  with n columns, where each column is a valid parameter vector; the
  number of columns determines the batch size.

- N:

  requested number of batches to return, the sample will be of size
  `batchSize*N` (batch size is the number of parallel Markov chains).

- burnIn:

  number of batches where the transition kernel will be adjusted to
  achieve an acceptance rate of below 10%.

- Sigma0:

  multivariate normal covariance of Markov chain transition kernel,
  defaults to the covariance of the initial parameters. If startPar is
  one vector, this matrix must be provided explicitly.

- dprior:

  a function that returns prior probability density values.

- deltaSpan:

  either an initial and final value for the ABC threshold delta, or a
  fixed value for delta that will never change.

- batchSize:

  the size of each batch, this should be a number that could be
  sufficient to calculate the covariance of in the given parameter
  space.

- parAcceptable:

  a function that can reject a parameter vector early based on
  user-requirements. Has to return a scalar Boolean. Use this to test
  for inequalities that you find difficult to encode in the prior.

- verbose:

  when TRUE, a progress bar is printed during burn-in and actual
  sampling.

## Value

a list containing a sample matrix and a vector of scores (values of
delta for each sample)

## Details

The distance of the ABC setting is compared to a threshold value
\\\\delta\\. The threshold doesn't need to be explicitly provided. You
can however provide a span of acceptable values in any order, the
smaller value will be used as a lower bound, the larger value will be
used initially.

The ABC procedure will attempt to converge first, using the initial
delta value, and decreasing it slowly using observed distance values.

This function always operates on a bundle of Markov chains. The size of
this bundle can be determined through `startPar` (a matrix of column
vectors). Each column will be used as the initial point of a Markov
chain. The chains will be resampled at each step during the convergence
phase, and decouple from one another once the burn-in is complete.

ABC methods (distance function, threshold delta) can be combined with
several other methods (like particle filters). Here we use several
parallel Markov chains to sample from the approximate posterior.

## Examples

``` r
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  o <- as_ode(m)
  ex <- experiments(m,o)
  C <- generate_code(o)
  c_path(o) <- write_c_code(C)
  so_path(o) <- shlib(o)
  s <- simulator.c(ex,o,parMap=log10ParMap)
  objFunc <- makeObjective(ex,s)
  p0 <- log10(values(m$Parameter))
  lowerBound <- p0 - 3
  upperBound <- p0 + 3
  dprior <- dUniformPrior(lowerBound,upperBound)
  rprior <- rUniformPrior(lowerBound,upperBound)
  X <- rprior(96)    # increase this number ...
  ## this always takes more than 5s to run:
  if (interactive()){
      abcSample <- abc_mcmc(
      objFunc,
      startPar=t(X),
      N=128,           # ... and this number
      Sigma0=cov(X),
      dprior=dprior
    )
  }
```
