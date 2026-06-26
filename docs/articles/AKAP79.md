# Sample AKAP79 with ABC

We use `icpm-kth/SBtabVFGEN` to read the model files and `icpm-kth/rgsl`
as an interface to `gsl_odeiv2` solvers (for initial value problems) to
simulate the model:

``` r
library(uqsa)
library(parallel)
```

Next, we import the model’s data tables from the SBtab files:

``` r
m <- model_from_tsv(uqsa_example("AKAP79")
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)
ex <- experiments(m,o)
```

The function `shlib` compiles the model to a *shared library* using
`R CMD SHLIB`. On a fairly normal system you can also make a shared
library using any c compiler with the options `cc -shared -fPIC ...`.

``` r
stopifnot(all(m$Parameter$scale == "log10"))
parABC <- values(m$Parameter) # initial value
```

The parameters are given in logarithmic space (log10). So, the natural
approach is to sample in logarithmic space: the ABC algorithm will
sample the logarithms of the model parameters. This way, the actual
parameters `10^parABC` are guaranteed to be positive.

We limit the sampling to ranges, as appropriate for each parameter:

``` r
ll <- m$Parameter$median - 2*m$Parameter$stdv
ul <- m$Parameter$median + 2*m$Parameter$stdv
```

The sampling procedure takes data in the list of experiments `ex` and
compares the data points to the solution that `gsl_odeiv2` returns. We
define a list of experiments, subdividing them into smaller groups and
processing the groups in sequence. Between the rounds of sampling the
posterior of every result is used as the prior distribution of the next
round. This mimics the arrival of data sets in sequence (from the lab).
The intermediate distributions will be modelled using `VineCopula`.

``` r
chunks <- list(c(3, 12,18, 9), c(2, 11, 17, 8), c(1, 10, 16, 7))
```

Problem Size and core distribution:

``` r
N <- 3000                    # sample size
```

Build a random number generator, and density function for the intial
prior:

``` r
rprior <- rUniformPrior(ll, ul)
dprior <- dUniformPrior(ll, ul)
```

The sampling loop:

``` r
set.seed(7619201)
options(mc.cores=detectCores())

start_time = Sys.time()
X <- rprior(N) # sample from the prior

for (i in seq_along(chunks)){
    I <- chunks[[i]]
    s <- simulator.c(experiments,o,log10ParMap,omit=3)
    Obj <- makeObjective(experiments[I],s)
    Z <- ABCSMC(
        Obj,
        t(X),
        Sigma=cov(X),
        dprior=dprior,
        delta=c(1,3) # either a range or c(initial, final)
    )
    C <- fitCopula(ABCMCMCoutput$draws)
    rprior <- rCopulaPrior(C)
    dprior <- dCopulaPrior(C)
    X <- Z$draws
}

end_time = Sys.time()
time_span = end_time - start_time

cat("Total time:\n",time_span)
```

The result is a collection of intermediate samples and a final posterior
sample. This article was build on this CPU:

``` sh
grep "model name" /proc/cpuinfo | head -n 1
head -n 1 /proc/meminfo
```

Let’s display the difference between the posterior and prior
distribution for the last chunk iteration:

``` r
posterior <- Z$draws
prior <- rprior(NROW(posterior))
uqsa::showPosterior(posterior[,seq(6)],prior[,seq(6)])
```

Alternatively:

``` r
if (require(hexbin)){
    hexbin::hexplom(posterior[,seq(6)])
} else {
    pairs(posterior[,seq(6)])
}
```
