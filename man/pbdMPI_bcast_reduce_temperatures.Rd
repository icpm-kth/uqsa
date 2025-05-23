% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{pbdMPI_bcast_reduce_temperatures}
\alias{pbdMPI_bcast_reduce_temperatures}
\title{Broadcast to other ranks and swap temperatures with any of them}
\usage{
pbdMPI_bcast_reduce_temperatures(i, B, LL, H, r, comm, cs)
}
\arguments{
\item{i}{MCMC iteration}

\item{B}{inverse temperature (parallel tempering)}

\item{LL}{log-likelihood value of current point}

\item{H}{algorithm's step size (often called epsilon in literature)}

\item{r}{MPI rank}

\item{comm}{MPI communicator}

\item{cs}{MPI comm size}
}
\description{
Using this function, at most two ranks will swap.
}
\details{
Given a current log-likelihood, temperature and step-size, this
function will broadcast a log-likelihood value to all other ranks
and they can each decide to swap temperatures with the root
process. Root is cycled around all ranks (round-robin).

Each other rank is allowed to make the
offer to swap. The root process decides which rank to swap with.
}
