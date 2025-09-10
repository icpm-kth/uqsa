# Uncertainty Quantification: Copula functions for ABC-MCMC
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#' Makes a Probability Density Estimate (from a sample)
#'
#' Given a sample (from some probability distribution) this function
#' makes a Copula fit to the source distribution using the VineCopula
#' package.
#'
#' @importFrom VineCopula RVineStructureSelect
#' @export
#' @param X sample that characterizes the traget distribution (rows)
#' @return as list: vineCop, U, Z, and Y where U are marginal
#'     probability samples, Z are cummulative density values for U,
#'     and Y are the probability density values of U.
fitCopula <- function(X){
  stopifnot(is.matrix(X))
  ncx <- ncol(X)
  ns <- nrow(X)
  eps <- 0.1
  npoints <- 5000
  # randomly pick sample points
  if(ns > npoints){
    I <- sample(1:ns, npoints, replace=FALSE)
  }else{
    I <- 1:ns
  }
  # add max and min
  I <- c(I, apply(X, 2, which.max))
  I <- c(I, apply(X, 2, which.min))
  I <- unique(I)
  Z <- U <- Y <-  matrix(NA, length(I), ncx)
  # must evaluate in real datapoints to
  # keep connection between params
  # this is a normal kernel, looks similar
  # to using the ecdf function

  for(i in 1:ncx){
    minx <- min(X[,i])
    maxx <- max(X[,i])
    ls <- minx-eps
    us <- maxx+eps
    U[,i] <- X[I,i]
		if (requireNamespace("ks")){
			Z[,i] = ks::kcde(X[,i], xmin=ls, xmax=us, eval.points = X[I,i])$estimate
			Y[,i] = ks::kde(X[,i], xmin=ls, xmax=us, eval.points = X[I,i])$estimate
		} else {
			stop("ks:: is required for this function.")
		}
  }

  # fit copula
  vineCop <- RVineStructureSelect(Z,indeptest = T)
  return(list(copula=vineCop, U=U, Z=Z, Y=Y))
}

#' Copula Formulation for Uniform Prior Distributions
#'
#' Covers the (simpler) special case where the `prior(x)` is iid uniform.
#' The return value has the same structure as the value of `fitCopula()`.
#'
#' @importFrom VineCopula RVineStructureSelect
#' @export
#' @param ll `ll[i]` is the lower limit of random variable `x[i]`
#' @param ul upper limit, analogous to ll.
#' @return list with: copula, U, Z, and Y entries.
makeIndepCopula <- function(ll, ul){
  npoints <- 5000
  np <- length(ll)
  Z <- U <- Y <- matrix(NA, npoints, np)
  for(i in 1:np){
     minx <- ll[i]
     maxx <- ul[i]
     U[,i] <- seq(minx, maxx, length.out = npoints)
     Z[,i] <- seq(0,1, length.out=npoints)
     Y[,i] <- rep(1/(maxx-minx), npoints)
  }
  vineCop <- RVineStructureSelect(Z, family=0)
  return(list(copula=vineCop, U=U, Z=Z, Y=Y))
}

#' (for testing) A non-Copula sampling function as fallback
#'
#' If the sample is not suited to infer a Copula (fitCopula fails),
#' this fuction uses much simpler rules to re-draw a new sample from
#' an older sample with some added noise.
#'
#' This can be used during testing, in cases where the acceptance was
#' very low and we have to deal with a very low quality sample This
#' function should work like `base::sample`, but adds small
#' noise. Missing values are always removed.
#'
#' @export
#' @param X an N×M matrix (N is the sample size), M is the number of variables (MCMC or ABC vars)
#' @param sdf factor to increase or decrease the standard deviation of the added noise
#' @param size size of returned sample (passed to `sample.int()`)
#' @param ... passed to `base::sample.int()`
#' @return a matrix with `size` rows and `nrow(X)` columns.
sampleWithNoise <- function(X,sdf=1e-2,...){
	i <- sample.int(NROW(X),...)
	m <- colMeans(X,na.rm=TRUE)
	C <- var(X,na.rm=TRUE)
	Y <- matrix(rnorm(length(i)*NCOL(X),mean=as.numeric(X[i,]),sd=(abs(m)+sqrt(diag(C))+1.0)*sdf),length(i),NCOL(X))
	return(Y)
}
