# simulates a CRNN ode model with extra work

This function calls a C function which solves an initial value problem,
derived from a CRNN.

## Usage

``` r
gsl_odeiv2_CRNN(
  name,
  experiments,
  l,
  nu,
  m,
  abs.tol = 1e-06,
  rel.tol = 1e-05,
  initial.step.size = 0.001,
  method = 0,
  time.out = 1,
  nstep = 0
)
```

## Arguments

- name:

  either the name of a file (shared library file) or the name of an ODE
  model to simulate (a shared library of the same name will be
  dynamically loaded and needs to be created first). If the name of the
  model is given, then the so file must have the same name in the
  current directory or a comment indicates its location.

- experiments:

  a list of `N` simulation experiments (time, parameters, initial value,
  events).

- l:

  a matrix of parameters with M columns, in log-space.

- nu:

  a stoichiometry matrix (N × R) where N is the number of state
  variables and R the number of reactions, all reactions are assumed to
  be reversible.

- m:

  modifiers – similar to stoichiometry, but indicates whether the
  species takes part in the reaction without being consumed.

- abs.tol:

  absolute tolerance, real scalar.

- rel.tol:

  relative tolerance, real scalar.

- initial.step.size:

  initial value for the step size; the step size will adapt to a value
  that observes the tolerances, real scalar.

- method:

  one of the integration methods bundled with GSL (see
  [method](https://icpm-kth.github.io/uqsa/reference/method.md) and
  [name_method](https://icpm-kth.github.io/uqsa/reference/name_method.md)).

- time.out:

  time limit in seconds, checked at every measurement time-point (in the
  data).

- nstep:

  maximum number of ODE integrator steps, checked at every step,
  defaults to unlimited (0).

## Value

a list of the solution trajectories y(t;p) for all experiments (named
like the experiments), as well as the output functions.

## Examples

``` r
if (FALSE) { # \dontrun{
  f <- uqsa_example("AKAR4")
  m <- model_from_tsv(f)
  ex <- experiments(m,as_ode(m,cla=FALSE))
  nu <- stoichiometric_matrix(m)
  l <- matrix(c(log(values(m$Parameter)),0),2,2,dimnames=list(rownames(m$Reaction),c("fwd","bwd")))
  C <- CRNN(NCOL(nu),initialValues=values(m$Compound),funcValues=formulae(m$Output))
  c.file <- tempfile("AKAR4_",fileext=".c")
  cat(C,file=c.file,sep='\n')
  so.file <- shlib(c.file,model.name="AKAR4")
  y <- gsl_odeiv2_CRNN(so.file,ex,l,nu,nu*0)
} # }
```
