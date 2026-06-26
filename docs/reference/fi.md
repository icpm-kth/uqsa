# Default gradient-Log-likelihood Function

Extracts the `gradLogLikelihood` values from the simulations attribute
of the parMCMC argument, requires:

- parMCMC has simulations attribute

- simulations list includes Fisher Information values (omit=0)

## Usage

``` r
fi(parMapJac = function(x) diag(1, length(x), length(x)))
```

## Arguments

- parMapJac:

  a function; maps parameter vectors to the Jacobian of the parameter
  transformation.

## Value

a scalar value: log(likelihood(data\|parMCMC))

## Details

This function will take the Fisher-Information-matrices claculated by
the ode solver in this package, and return the sum of those values over
all experiments. The gll-value the simulator returns is calculated with
the assumption of a normal distribution on measurement errors, and uses
the `identity` map between the MCMC variable and the model's parameters
by default (i.e. no transformation).

Like [ll](https://icpm-kth.github.io/uqsa/reference/ll.md) and
[gllf](https://icpm-kth.github.io/uqsa/reference/gllf.md) this function
does almost no work, it merely sums up the FI values calculated during
simulation, but it also performs a transformation of the Fisher
Information Matrix, taking the parameter-mapping between the
sampling-space and model-parameter-space into account.

The only argument is a function that takes the current MCMC variable,
`parMCMC` (a numeric vector), with all necessary attributes for smmala
to work (e.g. through initialization).

## Examples

``` r
# \donttest{
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
s <- simulator.c(ex,o,omit=0)
p <- values(m$Parameter)
attr(p,"simulations") <- s(p)
### without parameter transformations
gll <- gllf()
FI <- fi()
print(ll(p))
#> [1] -2491.245
print(gll(p))
#> [1] 1643.53085001   -2.87257643    0.03029852
print(FI(p))
#>               [,1]          [,2]          [,3]
#> [1,] 2628400.11438 -4.591588e+03 50.1228420741
#> [2,]   -4591.58828  8.021109e+00 -0.0875605903
#> [3,]      50.12284 -8.756059e-02  0.0009571775
# }
```
