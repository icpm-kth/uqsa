# SMMALA Update is an MCMC update function

During Markov chain Monte Carlo a given parameter needs to be updated,
the model needs to be simulated at the updated point.

## Usage

``` r
smmala_update(
  simulate,
  logLikelihood = ll,
  dprior = function(x) prod(dnorm(x)),
  gradLogLikelihood = gllf(log10ParMapJac),
  gprior = function(x) (-x),
  fisherInformation = fi(log10ParMapJac),
  fisherInformationPrior = 0,
  parAcceptable = function(p) all(is.finite(p))
)
```

## Arguments

- simulate:

  a function that simulates the model

- logLikelihood:

  a function that returns the log-likelihood value given the paramegter
  value, with simulations attached to the parameter as an attreibute
  (probably a closure)

- dprior:

  a function that returns the prior density of the given parameter
  vector

- gradLogLikelihood:

  any function that calculates or estimates the gradient of the
  log-likelihood function, for the chosen parameter mapping. Function
  must take one argument (the MCMC variable)

- gprior:

  a function that returns the gradient of the log-prior distribution.

- fisherInformation:

  a function that estimates the Fisher Information for a given MCMC
  variable (parMCMC).

- fisherInformationPrior:

  a constant fisherInformation of the prior distribution (or rather, the
  precision of the prior)

- parAcceptable:

  a function that can be used to reject a proposal based on the values
  of the parameters alone (shortcut to rejection, sans simulation)

## Details

Using the simulations, and an acceptance rule, the proposed update is
either accepted or rejected.

This function returns a closure `smmala`, with only `parMCMC` as it's
sole argument: `parProposal <- smmala(parGiven)`

An optional argument to this function is `parAcceptable`, during
sampling, when `metropolis` is called as the update function, and
`parAcceptable(parProposal)` returns `FALSE`, then metropolis shortcuts
to `retrun(parGiven)` without performing simulations.

This function can be used to weed out parameter combinations that would
result in obviously nonsensical simulations without wasting CPU-time.

The argument `fisherInformationPrior` is really the precision of the
prior (a constant matrix). It's role is additive to the
fisherInformation and is used to regularize the *final* Fisher
Information Matrix (makes it invertible).

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simulator.c(ex,o,omit=0)
### without parameter transformations
gll <- gllf()
FI <- fi()
dprior <- dNormalPrior(values(m$Parameter),m$Parameter$stdv)
gprior <- dNormalPrior(values(m$Parameter),m$Parameter$stdv)
p <- mcmc_init(1.0,values(m$Parameter),s,ll,dprior,gll,gprior,FI)
UP <- smmala_update(s,ll,dprior=dprior,gll,gprior=gprior,FI,solve(diag(m$Parameter$stdv)))
p2 <- UP(p)
## updated value:
print(p2)
#> [1]  0.01801046  0.10893704 10.21378666
#>              simulations: 3 (length)
#>            logLikelihood: -2491.24
#>                    prior: 10.5776
#>        gradLogLikelihood: 3 (length)
#>             gradLogPrior: 10.5776
#>        fisherInformation: 9 (length)
print(sum(abs(p2-p)))
#> [1] 0.01673416
```
