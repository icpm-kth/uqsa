% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gillespie.R
\name{simstoch}
\alias{simstoch}
\title{Simulate stochastic model}
\usage{
simstoch(experiments, model.so, parMap = identity)
}
\arguments{
\item{experiments}{list of expperiments, same as for the deterministic solvers.}

\item{model.so}{compiled C code for the model, this has to be a path with at least one slash in it, e.g.: ./model.so}

\item{parameters}{a numeric vector of appropriate size}
}
\value{
simulation result list
}
\description{
Simulate a stochastic model generated with
\code{uqsa::generateGillespieModel()}, using the solver in this package.
}
\details{
This will simulate all experimental conditions included in the list of experiments, including applying the inputs:
\code{u <- experiments[[i]]$input} - the input will be copied to the end of the model's internal parameter vector.

Like for deterministic models, we assume that there is a vector of
unknown parameter (a Markov chain variable, a vector of
optimization variables) and also known parameters (aka the input
parameters). The model itself does not distinguish between the two,
but one is the same between the experiments and one is different
between different experiments: \code{modelParam <- c(mcmcParam, inputParam)}
}
