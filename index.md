# Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

This is an R package that performs *parameter estimation*,
*uncertainty quantification*, and *global sensitivity analysis* using
Bayesian methods and ABC-MCMC sampling.

* **Source code:** https://github.com/icpm-kth/uqsa/

## Documentation

The Articles on this page are a user guide to this package. The
function reference is also accessible through R (`?uqsa::ABCMCMC`).

The package itself includes example models and a *run* script. We
advise users to copy the run scripts either from github or the package
and adapt them. The location of a *run* file can be found like this:

```R
uqsa::uqsa_example("AKAP79",pat="^run.*R$")
```
Or by browsing the folders on GitHub (`inst/extdata/`)

These examples are also covered by the articles: AKAR4, AKAP79, and CaMKII.

## Installation

The first two packages are optional in your own work, but you will
need them for our examples:

```R
remotes::install_github("icpm-kth/rgsl")       # requires gsl in your OS
remotes::install_github("icpm-kth/SBtabVFGEN") # if you plan to use SBtab
remotes::install_github("icpm-kth/uqsa")       # this package
```

- `rgsl` is an interface between R and the ODE solvers in the GNU Scientific Library
- `SBtabVFGEN` is an R package that works with models written in the SBtab format

### GNU Scientific Library (optional)

For the simulation backend using the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (gsl) needs to be installed in your OS.

|   Operating System | command                  |
|-------------------:|:-------------------------|
| Debian and Ubuntu: | `apt install libgsl-dev` |
|      Alpine Linux: | `apk add gsl`            |
|               Guix | `guix install gsl`       |
|         Arch Linux | `pacman -S gsl`          |
|             Gentoo | `emerge sci-libs/gsl`    |
|              MACOX | `brew install gsl`       |



