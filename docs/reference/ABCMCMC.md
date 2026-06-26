# Performs and Approximate Bayesian Computation Sampling of Model Parameters

Given a set of simulation experiments (list), a model, parameter
boundaries, this function will draw a sample of parameters from the
posterior probability density of the given problem.

## Usage

``` r
ABCMCMC(
  objectiveFunction,
  startPar,
  nSims,
  Sigma0,
  delta,
  dprior,
  batchSize = 100,
  parAcceptable = function(p) {
     all(is.finite(p))
 },
  allow.reg = FALSE
)
```

## Arguments

- objectiveFunction:

  function that, given a (vectorial) parameter as input, simulates the
  model, and outputs the distance between experimental data and data
  simulated from the model with the parameter provided in input

- startPar:

  starting value for the parameter vector

- nSims:

  requested sample size

- Sigma0:

  multivariate normal covariance of Markov chain transition kernel

- delta:

  ABC acceptance threshold

- dprior:

  a function that returns prior probability density values

- batchSize:

  number of points to produce to record one point

- parAcceptable:

  a function that can reject a parameter vector early based on
  user-requirements. Has to return a scalar boolean.

- allow.reg:

  allow regularization (logical), if TRUE, then Sigma will be made
  smaller once a very low acceptance rate is detected: one accepted
  update per batch

## Value

a list containing a sample matrix and a vector of scores (values of
delta for each sample)

## Details

Normally, ABC would produce a highly auto-correlated sample, wasting
lots of disk-space; `batchSize` can be used to thin out the sample,
recording every batchSize-th point: with `batchSize=100`, we perform 100
updates to the Markov chain variable and then save the state to the
sample. Higher batchSize numbers improve the apparent quality of the
sample, but create more work per sampled point.

## Examples

``` r
if (FALSE) { # \dontrun{
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  o <- as_ode(m)
  ex <- experiments(m,o)
  C <- generate_code(o)
  c_path(o) <- write_c_code(C)
  so_path(o) <- shlib(o)
  s <- simulator.c(ex,o)
  objFunc <- makeObjective(ex,s)
  startPar <- values(m$Parameter)
  lowerBound <- startPar - m$Paramster$stdv # m$Parameter$min
  upperBound <- startPar + m$Paramster$stdv # m$Parameter$max
  abcSample <- ABCMCMC(
    objFunc,
    startPar=values(m$Parameter),
    100,
    cov(rprior(1000)),
    delta=1,
    dprior=dUniformPrior(lowerBound,upperBound),
    batchSize = 100
  )
} # }
```
