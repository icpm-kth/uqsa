# This function proposes an MCMC candidate variable, and either accepts or rejects the candidate

This function receives a current MCMC variable, then calculates a
possible successor and returns it in the case of acceptance. It returns
the (old) current state upon rejection of the candidate.

## Usage

``` r
automatic_update(
  simulate,
  logLikelihood = ll,
  dprior = function(x) prod(dnorm(x)),
  gradLogLikelihood = NULL,
  gprior = NULL,
  fisherInformation = NULL,
  fisherInformationPrior = NULL,
  Sigma = NULL,
  parAcceptable = function(p) {
     TRUE
 }
)
```

## Arguments

- simulate:

  a function that simulates the model for a given parMCMC

- logLikelihood:

  a function that calculates log-likelihood values for given parMCMC

- dprior:

  prior density function

- gradLogLikelihood:

  a function that calculates the gradient of the log-likelihood for
  given parMCMC

- gprior:

  gradient of the prior density

- fisherInformation:

  a function that calculates approximates Fisher information matrices

- fisherInformationPrior:

  a constant matrix, the prior distributions fisher information

- Sigma:

  alternatively, `Sigma=solve(fisherInformationPrior)`, \[the inverse to
  fisherInformationPrior\] can be specified for the metropolis algorithm

- parAcceptable:

  a user shaped function that returns a Boolean scalar indicating
  whether a proposed parameter satisfies any user chosen constraints. A
  value of FALSE will trigger an early return and the model will not be
  simulated.

## Value

a function that returns possibly updated states of the Markov chain

## Details

The Markov chain has a current state (the MCMC variable, often x in
literature), but in the context of sampling the MCMC variables are used
as the parameters to a scientific model of some sort (and these often
have state variables, also x, or y). This is why we call the variables
parMCMC (parABC), or parCurrent\|Given\|Proposal, and similar.
