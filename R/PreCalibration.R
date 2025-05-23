# Uncertainty Quantification: Precalibration for ABC-MCMC
# Federica Milinanni (fedmil@kth.se)
# (based on: Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com))

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#' Determine a starting value for ABC's delta
#'
#' In ABC settings a model solution is compared to data with an
#' acceptance threshold: delta. This-pre calibration function attempts
#' to adjust this delta value.
#'
#' @export
#' @param objectiveFunction function that, given a (vectorial)
#'     parameter as input, (1) simulates the model with the given parameter, and (2) outputs the
#'     distance between experimental data and simulated data simulated.
#' @param npc sample size of pre-calibration.
#' @param rprior a function that generates random ABC variables,
#'     distributed according to the prior.
#' @param rep number of repetitions of the preCalibration process.
#' @param p fraction (top scoring) of sampled points to base Sigma on (Sigma is the covariance matrix for the moves proposed in the ABCMCMC algorithm).
#' @param sfactor scales Sigma up or down (Sigma is the covariance matrix for the moves proposed in the ABCMCMC algorithm).
#' @param delta ABC threshold.
#' @param num number of different starting parameter vectors (initial states of the chains) to generate. Usually, num is equal to the number of chain that will be run in the sampling procedure.
#' @return list with entries prePar (sampled parameters),  preDelta (distances between experimental data and trajectories produced with each of the parameters in prePar), Sigma (covariance matrix for the moves proposed in the ABCMCMC algortihm) and startPar (starting parameters for the ABCMCMC chains)
preCalibration <- function(objectiveFunction, npc=1000, rprior, rep = 1, p=0.05, sfactor=0.1, delta=0.01, num=1){
	nCores <- unlist(options("mc.cores"))
	if (is.null(nCores)){
		nCores <- parallel::detectCores()
		options(mc.cores=nCores)
	}
	# make npc a multiple of cores
	npc <- (max(2,ceiling(npc/nCores)))*nCores
	prePar <- t(rprior(npc))
	n <- dim(prePar)
	# split work
	dim(prePar) <- c(n[1],n[2]/nCores,nCores)
	preDelta<-unlist(mclapply(1:nCores, function(j) {apply(objectiveFunction(prePar[,,j]),2,max)}))
	dim(prePar) <- n
	# With the following loop we repeat the process and keep the npc best parameters and deltas
	for( i in 1:rep){
		p<-t(rprior(npc))
		newPrePar <- cbind(prePar, p)
		dim(p)<-c(n[1],n[2]/nCores,nCores)
		d <- unlist(mclapply(1:nCores,function(j) {apply(objectiveFunction(p[,,j]),2,max)}))
		newPreDelta <- c(preDelta,d)
		dim(p) <- n
		ix <- order(newPreDelta)[1:npc]
		preDelta <- newPreDelta[ix]
		prePar <- newPrePar[,ix]
	}
	
	M <- getMCMCPar(prePar, preDelta, p=p, sfactor=sfactor, delta = delta, num=num)
	return(list(prePar=prePar, preDelta=preDelta, Sigma=M$Sigma, startPar=M$startPar))
}

#' Selects MCMC scheme specific setup parameters
#'
#' The MCMC scheme uses a transition kernel. This function returns the
#' parameters of that transition kernel. Better parameters make the
#' Markov chain perform better (i.e. lower auto-correlation).
#'
#' @param prePar a sample of parameters from pre-Calibration
#' @param preDelta distance values (scores) for those parameters
#' @param p fraction (top scoring) of sampled points to base Sigma on
#' @param sfactor scales Sigma up or down
#' @param delta ABC threshold
#' @param num number of different starting parameter vectors.
#' @return Sigma and startPar (matrix with `num` rows) as a list
getMCMCPar <- function(prePar, preDelta, p=0.05, sfactor=0.1, delta=0.01, num=1){
	if (all(is.na(preDelta)) || is.null(preDelta))
		stop("no usable pre-calibration parameters.")
	if (sum(is.finite(preDelta)) < num){
		cat(sprintf("There are %i valid (finite) distance scores in the pre-calibration sample (%i starting positions requested).\n",sum(is.finite(preDelta)),num))
		stop("The number of valid points is too small to make MCMC starting parameters.")
	}

	prePar <- prePar[,!is.na(preDelta)]
	preDelta <- preDelta[!is.na(preDelta)]
	nk <- max(ceiling(ncol(prePar)*p), num)
	pick1 <- which(preDelta <= delta)   # pick all pars that meet threshold
	pick2 <- head(order(preDelta, decreasing = FALSE),nk) # pick top p percent
	if(length(pick1)>length(pick2)){
		pick <- pick1
	}else{
		pick <- pick2
		warning(sprintf("distances between experiment and simulation are too big; selecting the best (%i) parameter vectors.\n",length(pick)))
	}
	Scorr <- cor(t(prePar[,pick]))
	diag(Scorr) <- 1
	sdv <- apply(t(prePar[,pick]), 2, sd)
	Sigma <- sfactor * Scorr * tcrossprod(sdv)
	startPar <- prePar[,sample(pick, num, replace = FALSE),drop=FALSE]
	list(Sigma=Sigma, startPar=startPar)
}

