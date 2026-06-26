# print information about the mcmc variable

Some mcmc variables have many attributes, which clutter the screen when
accidentally printed. This function prevents these long printouts.

## Usage

``` r
# S3 method for class 'mcmcVariable'
print(x, ...)
```

## Arguments

- x:

  the variable

- ...:

  requirement of print generic, not used.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
dprior <- dNormalPrior(values(m$Parameter),m$Parameter$stdv)
s <- simfi(ex,o)
p <- mcmc_init(1.0,values(m$Parameter),s,dprior=dprior)
print(p)
#> kf_C_AKAR4 kb_C_AKAR4 kcat_AKARp 
#>      0.018      0.106     10.200 
#>              simulations: 3 (length)
#>            logLikelihood: -2491.25
#>                    prior: 10.5823
```
