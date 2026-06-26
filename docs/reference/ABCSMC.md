# Performs and Approximate Bayesian Computation as a Particle Filter

Given a set of simulation experiments (list), a model, parameter
boundaries, this function will draw a sample of parameters from the
posterior probability density of the given problem.

## Usage

``` r
ABCSMC(
  objectiveFunction,
  startPar,
  Sigma = 2 * cov(startPar),
  dprior,
  delta = c(2, 0.5),
  parAcceptable = function(p) {
     all(is.finite(p))
 },
  messages = FALSE
)
```

## Arguments

- objectiveFunction:

  a function that can simulate the model for a batch of parameter
  vectors provided as a matrix of columns (batches)

- startPar:

  a matrix that has the same shape as the desired sample, but
  transposed, this can be a sample from the prior or a pre-conditioned
  sample that approximates the posterior, e.g.: t(rprior(1000))

- Sigma:

  multivariate normal covariance of Markov chain transition kernel

- dprior:

  a function that returns prior probability density values

- delta:

  ABC acceptance threshold, either a scalar, then it is the initial
  value of delta, or a pair of values, then it is the starting value and
  the final value of delta: `c(initialDelta,finalDelta)`

- parAcceptable:

  is a rejection-shortcut function; if `parAcceptable(p)` returns
  `FALSE` for a specific value of `p`, it means that simulations
  shouldn't even be attempted.

- messages:

  a logical value indicating whether log messages should be printed

## Value

a list containing a sample matrix and a vector of scores (values of
delta for each sample)

## Details

This is a variant of ABC where the entire batch is simulated with one
call to the simulator. startPar is the initial batch to be simulated: it
is a matrix where columns are different parameter vectors (e.g. prior
sample members). In other words: `startPar[,i]` must be a valid argument
for the objectiveFunction.

The objective function is a closure

The Objective-Function `objectiveFuntion(P)` should return a matrix with
`n` rows, where `n` is the number of simulation experiments (and thus
data-sets), and `m` columns, where `m` is the number of
parameterizations `NCOL(P)`.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(parallel)
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  ex <- experiments(m,as_ode(m,cla=FALSE))
  G <- as_cme(m)         # for Gillespie solver
  C <- generate_code(G)
  c_path(G) <- write_c_code(C)
  so_path(G) <- shlib(G)
  muX <- m$Parameter$value
  sdX <- m$Parameter$stdv
  rprior <- rNormalPrior(log(muX^2/(muX^2+sdX^2)),sqrt(log(1+sdX^2/muX^2)))
  dprior <- dNormalPrior(log(muX^2/(muX^2+sdX^2)),sqrt(log(1+sdX^2/muX^2)))
  s <- simstoch(ex,G,logParMap)
  O <- makeObjective(ex,s)
  X <- rprior(1000)
  colnames(X) <- rownames(m$Parameter)
  posterior <- ABCSMC(O,t(X),Sigma=cov(X),dprior=dprior,delta=c(0.5,1.5))
} # }
```
