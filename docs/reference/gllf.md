# Default gradient-log-likelihood Function

Extracts the `FisherInformation` values from the simulations attribute
of the parMCMC argument, requires:

- parMCMC has simulations attribute

- simulations list includes gradient values (omit \<2)

## Usage

``` r
gllf(parMapJac = function(x) diag(1, length(x), length(x)))
```

## Arguments

- parMapJac:

  a function; maps parameter vectors to the Jacobian of the parameter
  transformation.

## Value

a numeric vector: grad(log(likelihood(data\|parMCMC)))

## Details

This function will take the log-likelihood gradient values claculated by
the ode solver in this package, and return the sum of those vectors over
all experiments. The gll-value the simulator returns is calculated with
the assumption of a normal distribution on measurement errors, and uses
the "identity" map between MCMC parameters and model-parameters by
default (i.e. no transformation).

Like [ll](https://icpm-kth.github.io/uqsa/reference/ll.md) this function
does almost no work, it merely sums up the gradient values calculated
during simulation, but it also performs a transformation of the gradient
vector, taking the parameter-mapping between the sampling-space and
model-parameter-space into account.

The returned function takes one argument, the MCMC variable `parMCMC` (a
numeric vector). This variable requires all smmala specific attributes.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simulator.c(ex,o,omit=1) # not 3
p <- values(m$Parameter)
attr(p,"simulations") <- s(p)
print(ll(p))
#> [1] -2491.245
trivialJac <- \(x) diag(1,length(x),length(x)) # the default
gll <- gllf(parMapJac=trivialJac)
print(gll(p))
#> [1] 1643.53085001   -2.87257643    0.03029852
```
