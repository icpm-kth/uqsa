# Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

This is an R package that performs *parameter estimation*,
*uncertainty quantification*, and *global sensitivity analysis* using
Bayesian methods and ABC-MCMC sampling.

* **Source code:** https://github.com/icpm-kth/uqsa/

The Articles on this page are a user guide to this package. As always,
the function reference is also accessible within R (`?uqsa::ABCMCMC`),
after installation. See *Get Started* for detailed installation
instructions.

```R
if (!require(remotes)) install.packages("remotes")
remotes::install_github("icpm-kth/uqsa",dependencies=TRUE)
```

We have built a toolset for _model building_, _automated code
generation_ and _simulation_ around this package. The functinality of
the entire toolset is covered on this page, not only this R package.

