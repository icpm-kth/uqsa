% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior.R
\name{rCopulaPrior}
\alias{rCopulaPrior}
\title{rCopulaPrior returns a function that generates random values from the copula model}
\usage{
rCopulaPrior(Copula)
}
\arguments{
\item{Copula}{the return value of fitCopula()}
}
\value{
a matrix of random values
}
\description{
The returned function generates n random vectors, as rows of a matrix.
}
