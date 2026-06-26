# Initialize the Markov chain

This function must append all required attributes to the MCMC varible,
for the Markov chain to update correctly.

## Usage

``` r
mcmc_init(
  beta,
  parMCMC,
  simulate,
  logLikelihood = ll,
  dprior = function(x) prod(rnorm(x)),
  gradLogLikelihood = NULL,
  gprior = NULL,
  fisherInformation = NULL
)
```

## Arguments

- beta:

  inverse temperature for the Markov chain (parallel tempering)

- parMCMC:

  a plain starting value for the Markov chain

- simulate:

  a closure that maps the MCMC variable to simulation results (the
  simulation experiments are enclosed in this function).

- logLikelihood:

  a function that maps simulations to logLikelihood values

- dprior:

  density of the prior distribution

- gradLogLikelihood:

  the gradient function of the logLikelihood (optional) – only if the
  algorithm requires it

- gprior:

  the gradient pf the log-prior (for SMMALA and similar algorithms).

- fisherInformation:

  a function that calculates the Fisher Information matrix

## Value

the same starting parameter vector, but with attributes.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
p0 <- values(m$Parameter)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simfi(ex,o)
dprior <- dNormalPrior(p0,m$Parameter$stdv)
p <- mcmc_init(1.0,p0,s,dprior=dprior)
print(names(attributes(p))) ## now has attributes necessary for MCMC
#> [1] "names"         "unit"          "beta"          "simulations"  
#> [5] "prior"         "logLikelihood" "class"        
```
