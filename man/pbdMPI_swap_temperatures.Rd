% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{pbdMPI_swap_temperatures}
\alias{pbdMPI_swap_temperatures}
\title{Communicate with other ranks and swap beta}
\usage{
pbdMPI_swap_temperatures(i, B, LL, H, r, comm, cs)
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
Given a current log-likelihood, temperature and step-size, this
decides whether to send or receive the same variables from a
neighboring process and swap temperatures with them.
}
\details{
Only one of the two communicating ranks is allowed to make the
decision to swap, because a random variable is used to make this
decision. Which rank will make this swap decision needs to be
determined somehow. Currently we alternate this reposibility based
on the current iteration \code{i}.
}
