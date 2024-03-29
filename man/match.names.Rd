% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runModel_StochasticSimulation.R
\name{match.names}
\alias{match.names}
\title{Find the variable names in a formula}
\usage{
match.names(chrv)
}
\arguments{
\item{chrv}{a character vector as returned by parse.formula}
}
\value{
coefficients, as a vector
}
\description{
A reaction formula has reactants and products, separated by \verb{<=>},
with reactants on the left and products on the right (by convention).
Each of those is a plus separated list of reacting compounds and
modifiers, with optional coefficients, e.g.: \verb{A + 2 B <=> AB2}
}
\details{
Once the formula is split into left and right side, this function
determines the names. For the above example, this function
returns \code{c("A","B")} for the left side and \code{"AB2"} for the right side.
}
\examples{
lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.names)
lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.names)
}
