# Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

This is an R package that performs *parameter estimation*,
*uncertainty quantification*, and *global sensitivity analysis* using
Bayesian methods and ABC-MCMC sampling.

* **Source code:** https://github.com/icpm-kth/uqsa/

## Documentation

Some topics are covered in Articles on this page. There is also a function reference.

The package includes example models. These examples are also covered by the articles: AKAR4, AKAP79, and CaMKII.

We also have a GitHub Wiki.

* **Wiki:** https://github.com/icpm-kth/uqsa/wiki/Home

## Installation

```R
remotes::install_github("icpm-kth/SBtabVFGEN")
remotes::install_github("icpm-kth/uqsa")
```

### GNU Scientific Library (optional)

For the simulation backend using the GNU Scientific Library, the gsl needs to be installed in your OS.

Debian and Ubuntu:
```sh
apt install libgsl-dev
```

Alpine Linux:
```sh
apk add gsl
```

Guix:
```sh
guix install gsl
```

Arch Linux:
```sh
pacman -S gsl
```

Gentoo:
```sh
emerge sci-libs/gsl
```

And finally, in R:

```R
remotes::install_github("icpm-kth/rgsl")
```

