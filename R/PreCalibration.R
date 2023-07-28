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
#' @param experiments a list of simulation experiments. Same as for
#'     ABCMCMC.
#' @param modelName this name will be used to find the file and
#'     functions within the file according to naming conventions (the
#'     model name and file can differ, the file-name can be attached
#'     to the name as a comment).
#' @param parMap a remapping function that takes ABC sampling
#'     variables and returns valid model parameters
#' @param npc sample size of pre-calibration.
#' @param rprior a function that generates random ABC variables,
#'     distributed according to the prior.
#' @param getScore a function that maps the model's output values to
#'     ABC score values (in comparison to data).
#' @param rep number of repetitions of the preCalibration process
#' @return list with entries preDelta and prePar, final values of
#'     calibration run
preCalibration <- function(objectiveFunction, npc=1000, rprior, rep = 1){
	nCores <- options()$mc.cores
	# make npc a multiple of cores
	npc <- ceiling(npc/nCores)*nCores
	prePar <- t(rprior(npc))
	n <- dim(prePar)
	# split work
	dim(prePar) <- c(n[1],n[2]/nCores,nCores)
	preDelta<-apply(Reduce(rbind,mclapply(1:nCores,function(j) {objectiveFunction(prePar[,,j])})),2,mean)
	dim(prePar) <- n
	# With the following loop we repeat the process and keep the npc best parameters and deltas
	for( i in 1:rep){
		newPrePar <- cbind(prePar, p<-t(rprior(npc)))
		dim(p)<-c(n[1],n[2]/nCores,nCores)
		d <- apply(Reduce(rbind,mclapply(1:nCores,function(j) {objectiveFunction(p[,,j])})),2,mean)
		newPreDelta <- c(preDelta,d)
		dim(p) <- n
		ix <- order(newPreDelta)[1:npc]
		preDelta <- newPreDelta[ix]
		prePar <- newPrePar[ix,]
	}
	return(list(preDelta=preDelta, prePar=prePar))
}

#' Selects MCMC scheme specific setup parameters
#'
#' The MCMC scheme uses a transition kernel. This function returns the
#' parameters of that transition kernel. Better parameters make the
#' Markov chain perform better (i.e. lower auto-correlation).
#'
#' @export
#' @param prePar a sample of parameters from pre-Calibration
#' @param preDelta distance values (scores) for those parameters
#' @param p fraction (top scoring) of sampled points to base Sigma on [default 0.05]
#' @param sfactor scales Sigma up or down [default 0.1]
#' @param delta ABC threshold [default 1e-2]
#' @param num number of different starting parameter vectors [default 1].
#' @return Sigma and startPar (matrix with `num` rows) as a list
getMCMCPar <- function(prePar, preDelta, p=0.05, sfactor=0.1, delta=0.01, num=1){
	if (all(is.na(preDelta)) || is.null(preDelta))
		stop("no usable pre-calibration parameters.")
	if (sum(is.finite(preDelta)) < num){
		cat(sprintf("There are %i valid (finite) distance scores in the pre-calibration sample (%i starting positions requested).\n",sum(is.finite(preDelta)),num))
		stop("The number of valid points is too small to make MCMC starting parameters.")
	}

	prePar <- prePar[!is.na(preDelta),]
	preDelta <- preDelta[!is.na(preDelta)]
	nk <- ceiling(nrow(prePar)*p)
	pick1  <- which(preDelta <= delta)   # pick all pars that meet threshold
	pick2 <- order(preDelta, decreasing = FALSE)[1:nk] # pick top p percent
	if(length(pick1)>length(pick2)){
		pick <- pick1
	}else{
		pick <- pick2
		warning(sprintf("distances between experiment and simulation are too big; selecting the best (%i) parameter vectors.\n",length(pick)))
	}
	Scorr <- cor(prePar[pick,])
	diag(Scorr) <- 1
	sdv <- apply(prePar[pick,], 2, sd)
	Sigma <- sfactor * Scorr * tcrossprod(sdv)
	startPar <- prePar[sample(pick, num, replace = FALSE),]
	list(Sigma=Sigma, startPar=startPar)
}

