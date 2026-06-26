# Metropolis Update is an MCMC update function

During Markov chain Monte Carlo a given parameter needs to be updated,
the model needs to be simulated at the updated point.

## Usage

``` r
metropolis_update(
  simulate,
  logLikelihood = ll,
  dprior = function(x) prod(dnorm(x)),
  Sigma = NULL,
  parAcceptable = function(p) {
     all(is.finite(p))
 }
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

- Sigma:

  the transition kernel's covariance matrix.

- parAcceptable:

  a function that can be used to reject a proposal based on the values
  of the parameters alone (shortcut to rejection, sans simulation)

## Details

Using the simulations, and an acceptance rule, the proposed update is
either accepted or rejected.

This function returns a closure `metropolis`, with only `parMCMC` as
it's sole argument: `parProposal <- metropolis(parGiven)`

An optional argument to this function is `parAcceptable`, during
sampling, when `metropolis` is called as the update function, and
`parAcceptable(parProposal)` returns `FALSE`, then metropolis shortcuts
to `retrun(parGiven)` without performing simulations.

This function can be used to weed out parameter combinations that would
result in obviously nonsensical simulations without wasting CPU-time.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simulator.c(ex,o,omit=0)
dprior <- dNormalPrior(values(m$Parameter),m$Parameter$stdv)
p <- mcmc_init(1.0,values(m$Parameter),s,ll,dprior)
UP <- metropolis_update(s,ll,dprior=dprior,Sigma=diag(m$Parameter$stdv)^2)
p2 <- UP(p)
## updated value:
print(p2)
#> kf_C_AKAR4 kb_C_AKAR4 kcat_AKARp 
#>      0.018      0.106     10.200 
#>              simulations: 3 (length)
#>            logLikelihood: -2491.25
#>                    prior: 10.5823
print(sum(abs(p2-p)))
#> [1] 0
```
