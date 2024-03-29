% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{gradLogLikelihood}
\alias{gradLogLikelihood}
\title{Default log-likelihood function, gradient}
\usage{
gradLogLikelihood(model, experiments, parMap = identity, parMapJac = 1)
}
\arguments{
\item{experiment}{will be compared tp the simulation results}
}
\description{
This returns a function g(x,simulations), which maps simulation
results and the MCMC variables x to the gradient of log(likelihood)
values withj respect to x. The experiments are used implicitly;
simulations is a list as returned by rgsl::r_gsl_odeiv2_outer().
}
