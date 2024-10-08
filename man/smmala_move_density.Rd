% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{smmala_move_density}
\alias{smmala_move_density}
\title{SMMALA transition kernel density}
\usage{
smmala_move_density(
  beta,
  parProposal,
  parGiven,
  fisherInformationPrior,
  eps = 0.01
)
}
\arguments{
\item{beta}{inverse temperature (parallel tempering)}

\item{parGiven}{given point}
}
\value{
SMMALA proposal point
}
\description{
The Simiplified Manifold Metropolis Adjusted Langevin Algorithm uses a
move instriction that uses a Gaussian kernel that is shifted away
from the current point.
}
