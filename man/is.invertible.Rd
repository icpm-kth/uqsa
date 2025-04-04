% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{is.invertible}
\alias{is.invertible}
\title{checks whether a given matrix is a valid, invertible fisherInformation}
\usage{
is.invertible(G = NULL, abs_tol = 1e-11)
}
\arguments{
\item{G}{a matrix}

\item{abs_tol}{absolute tolerance for the reciprocal condition number of G}
}
\value{
TRUE or FALSE
}
\description{
This matrix has to be symmetric and invertible. But, because the
matrix has a perhaps sketchy origin, it could be defective in all
possible ways.
}
