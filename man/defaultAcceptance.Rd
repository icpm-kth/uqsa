% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runModel.R
\name{defaultAcceptance}
\alias{defaultAcceptance}
\title{default ABC acceptance probability function for one experiment}
\usage{
defaultAcceptance(funcSim, dataVAL, dataERR = max(dataVAL))
}
\arguments{
\item{funcSim}{a matrix, contains model solution (output values),
columns of output vectors}

\item{dataVAL}{a data.frame of experimental data}

\item{dataERR}{a data.frame of measurement errors, if available,
defaults to the maximum data value.}
}
\description{
if each experiment corresponds to one simulation and is fully
quanitified by itself, then calculating the overall distance
between data and experiment can be done one by one. This function
describes the default way a simulation is compared to data.
}
\details{
If the data is more complex, and two or more simulations are needed
to calculate one distance value then the objective-Function needs
to be entirely user-supplied. This is the case with experiments
that have a "control" -- this is needed when the measurement is in
arbitrary units and only makes sense comparatively to a secondary
(control) scenario.

This function will be used if none is provided by the user.

The funcSim values need to be supplied as a matrix of size N×T with
N the length of the model's output vectors and T the amount of
measurement times (this is how the rgsl package returns the
simulation results).
}
