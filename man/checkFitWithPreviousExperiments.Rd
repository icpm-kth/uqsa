% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABCMCMCFunctions.R
\name{checkFitWithPreviousExperiments}
\alias{checkFitWithPreviousExperiments}
\title{ABC acceptance of currently sampled values given old data (Prior)}
\usage{
checkFitWithPreviousExperiments(draws, objectiveFunction, delta)
}
\arguments{
\item{draws}{matrix of sampled values (to be filtered).}

\item{delta}{the acceptance threshold.}

\item{modelName}{name, is used to find the file and model functions
therein (comment(modelName) can contain a file-name if it
differs from `modelName.R`).}

\item{experiments}{a list of experiments (all of them, or up to
currentExpSet).}

\item{parMap}{optional remapping function:
parModel<-parMap(parABC); the ABC variables will be transformed
to make the parameter vector acceptable to the model. This is
necessary whenever ABC sampling happens in a different space
than the model parameter domain for whatever reason.}

\item{getScore}{scoring function.}

\item{nCores}{number of cores to use in parallel::mclapply() calls.}
}
\value{
a filtered subset of acceptable parameter draws
}
\description{
The prior probability density model using copulas and vines is not
perfect, so values sampled from an imperfect prior estimate can be
checked against old data.
}