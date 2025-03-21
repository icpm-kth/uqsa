% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runModel_StochasticSimulation.R
\name{importReactionsSSA}
\alias{importReactionsSSA}
\title{Functions to construct and run the stochastic simulation using GillespieSSA2 package}
\usage{
importReactionsSSA(model.tab, compile = TRUE)
}
\arguments{
\item{model}{the model, represented by a list of data.frames with SBtab content}
}
\value{
a list of compiled GillespieSSA2::reaction items
}
\description{
This translates the Reaction network into the specific form required by GillespieSSA2
}
