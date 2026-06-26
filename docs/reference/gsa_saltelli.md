# Outputs the global sensitivity scores SI and SIT, calculated by the Sobol-Homma-Saltelli method

M1, M2, and N are matrices prepared by
[`uqsa::saltelli_prior()`](https://icpm-kth.github.io/uqsa/reference/saltelli_prior.md).
The parameters (rows) from these matrices need to be simulated (using
any method), to obtain `fM1`, `fM2` and `fN`.

## Usage

``` r
gsa_saltelli(fM1, fM2, fN, subtract.mean = TRUE)
```

## Arguments

- fM1:

  output (f)unction values for `M1`, \\n_S \times n_O\\

- fM2:

  output (f)unction values for `M2`, \\n_S \times n_O\\

- fN:

  output (f)unction values for `N`, \\n_S \times n_O \times n_P\\

- subtract.mean:

  whether or not to subtract the column-means from all matrices/arrays

## Value

a list with sensitivity indices `$SI` and total sensitivities `$SIT`

## Details

These matrices are shaped similarly to `M1`, `M2` and `N` respectively,
but now the parameters are replaced by the effects they have on a
observable of interest (the output). It can be the vector valued output
at a specific (single) time-point or a scalar output at different
time-points.

See Geir Halnes et al. (Halnes, Geir, et al. J. comp. neuroscience 27.3
(2009): 471.

## Examples

``` r
m <- model_from_tsv(uqsa_example("AKAR4"))
o <- write_and_compile(as_ode(m))
ex <- experiments(m,o)
s <- simulator.c(ex[1],o)
p0 <- values(m$Parameter)
rprior <- rUniformPrior(p0/2,p0*2)
SP <- saltelli_prior(1000,rprior)
fM1 <- t(s(t(SP$M1))[[1]]$func[1,,])
fM2 <- t(s(t(SP$M2))[[1]]$func[1,,])
fN <- lapply(asplit(SP$N,3),\(N) t(s(t(N))[[1]]$func[1,,]))
fN <- simplify2array(fN)
GSA <- gsa_saltelli(fM1,fM2,fN)
print(names(GSA))
#> [1] "SI"  "SIT"
cat(
  "average relative senitivity S(p1) / S(p2): ",
  mean(abs(GSA$SI[,1]/GSA$SI[,2]),na.rm=TRUE)
)
#> average relative senitivity S(p1) / S(p2):  14864.47
```
