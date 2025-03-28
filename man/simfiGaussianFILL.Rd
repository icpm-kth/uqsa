% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{simfiGaussianFILL}
\alias{simfiGaussianFILL}
\title{\link{SMMALA} Extract the approximate Fisher infomration from the simfi results}
\usage{
simfiGaussianFILL(
  ParMapJac = function(x) {
     diag(1, length(x))
 }
)
}
\arguments{
\item{ParMapJac}{Jacobian of the parMap function}
}
\value{
the approximate Fisher information of the Gaussian log-likelihood with respect to
the MCMC variable, useful for SMMALA
}
\description{
This function extracts the approximate Fisher information matrix G
of the log-likelihood and transforms it using the
Jacobian of the parameter map between Markov chain variables and
model parameters.
}
\details{
The simfi() values are with respect to the raw model parameters,
while this function rephrases them in terms of the Markov chain's
position.
}
