# Uncertainty Quantification: ABC-MCMC with copulas
# Federica Milinanni (fedmil@kth.se)
# (based on: Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# and based on modifications: 2021 by Joao Antunes (joaodgantunes@gmail.com) and Olivia Eriksson (olivia@kth.se))

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#' Performs and Approximate Bayesian Computation Sampling of Model Parameters
#'
#' Given a set of simulation experiments (list), a model, parameter
#' boundaries, this function will draw a sample of parameters from the
#' posterior probability density of the given problem.
#'
#' Initially this function performs a similar job as an optimizer, and
#' then transitions to MCMC sampling.
#'
#' @export
#' @param experiments a list of experiments
#' @param modelName model functions will be assumed to have this
#'     prefix (comment(modelName) can contain a file-name)
#' @param startPar starting value for the parameter vector
#' @param parMap re-mapping function between the ABC MCMC variables and model
#'     parameters
#' @param nSims requested sample size
#' @param Sigma0 multivariate normal covariance of Markov chain
#'     transition kernel
#' @param delta ABC acceptance threshold
#' @param dprior a function that returns prior probability density
#'     values
#' @param getScore a function(model_output,experiment) that returns
#'     scores
#' @param nCores setting for multicore package
#' @return a sample matrix
ABCMCMC <- function(experiments, modelName, startPar, parMap, nSims, Sigma0, delta, dprior, getScore, nCores=detectCores()){
  cat("Started chain.\n")
  Sigma1 <- 0.25*diag(diag(Sigma0))
  curDelta <- Inf
  np <- length(startPar)
  scount <- 1
  curPar  <- startPar
  numExperiments <- length(experiments)
  out <- runModel(experiments, modelName, curPar, parMap, nCores)
  curDelta <- mclapply(1:length(out),
                       function(i) getScore(out[[i]], experiments[[i]][["outputValues"]]),
                       mc.preschedule = FALSE,
                       mc.cores = nCores)
  curDelta <- unlist(curDelta)

  #Similarly to what we did in the preCalibration, we average the score obtained with (the same) startPar applied to all the simulations (corresponding to different experiments setup)
  #As in preCalibration, we can use - for instance - the sum of squares
  curDelta <- mean(curDelta)

  if(is.na(curDelta)){
    cat("\n*** [parUpdate] curDelta is NA. Replacing it with Inf ***\n")
    curDelta <- Inf
  }
  curPrior <- dprior(curPar)
  draws <- matrix(0, nSims,np)
  n <- 0
  while (n < nSims){
    if(runif(1)<=0.95){
      canPar <- mvrnorm(n=1, curPar, Sigma0)
    }else{
      canPar <- mvrnorm(n=1, curPar, Sigma1)
    }
    out <- parUpdate(experiments, modelName, parMap, curPar, canPar, curDelta, curPrior, delta, dprior, getScore, nCores)
    curPar <- out$curPar
    curDelta <- out$curDelta
    curPrior <- out$curPrior
    scount <- ifelse(!all(curPar == canPar), scount+1, 1)
    if(!is.na(curDelta <= delta & all(curPar == canPar))){
      if (curDelta <= delta & all(curPar == canPar)){
        n=n+1
        draws[n,] = curPar
      }
    } else {
      cat('NA when evaluating: "curDelta <= delta & all(curPar == canPar)"')
    }

    if(scount>500){ #terminate chain if stuck
      cat('Aborted chain.\n')
      return(draws)
    }
  }
  cat("Finished chain.\n")
  return(draws)
}

#' Updates Parameter Values
#'
#' under valid ABC update conditions (successful simulation) the
#' parameters are updated to new values.
#'
#' @export
#' @param experiments list of simulation experiments
#' @param modelName (characters), this will be used to find the model
#'     file and the functions in that file
#' @param parMap optional remapping function set for passing parameters to the
#'     model: parModel<-parMap(parABC)
#' @param curPar current parameter values (as ABC samples them)
#' @param canPar candidate parameter values (for MCMC)
#' @param curDelta current distance between data and simulation, if
#'     the MCMC chain has not yet reached any point where this is
#'     below the threshold (delta), this can be accepted as the new
#'     current state for the chain.
#' @param curPrior current Prior values given curPar
#' @param delta distance threshold for ABC
#' @param dprior prior probability density function
#' @param getScore a scoring function
#' @param nCores number of cores to use in mclapply().
#' @return updated values for curPar, curDelta, and curPrior
parUpdate <- function(experiments, modelName, parMap, curPar, canPar, curDelta, curPrior, delta, dprior, getScore, environment, nCores=detectCores()){
  numExperiments <- length(experiments)
  invisible(capture.output(out <- runModel(experiments, modelName, parABC=canPar, parMap, mc.cores=nCores)))

  if(is.null(out)){
    canDelta <- Inf
    canPrior <- 0
  }else{
    canDelta <- mclapply(1:length(out), function(i) getScore(out[[i]], experiments[[i]][["outputValues"]]), mc.preschedule = FALSE, mc.cores = nCores)
    canDelta <- unlist(canDelta)
    #Similarly to what we did in the preCalibration and in ABCMCMC, we average the score obtained with (the same) startPar applied to all the simulations (corresponding to different experiments setup)
    #As in preCalibration and ABCMCMC, we can use - for instance - the sum of squares
    canDelta <- mean(canDelta)
    if(is.na(canDelta))
    {
      cat("\n*** [ABCMCMC] canDelta is NA. Replacing it with Inf ***")
      canDelta <- Inf
    }
    canPrior <- dprior(canPar)
  }

  if (canDelta <= max(delta,curDelta)){
    if (runif(1) <= canPrior/curPrior){
      curDelta <- canDelta
      curPrior <- canPrior
      curPar <- canPar
    }
  }
  return(list(curPar=curPar, curDelta=curDelta, curPrior=curPrior))
}

#' ABC acceptance of currently sampled values given old data (Prior)
#'
#' The prior probability density model using copulas and vines is not
#' perfect, so values sampled from an imperfect prior estimate can be
#' checked against old data.
#'
#' @export
#' @param modelName name, is used to find the file and model functions
#'     therein (comment(modelName) can contain a file-name if it
#'     differs from `modelName.R`).
#' @param draws matrix of sampled values (to be filtered).
#' @param experiments a list of experiments (all of them, or up to
#'     currentExpSet).
#' @param parMap optional remapping function:
#'     parModel<-parMap(parABC); the ABC variables will be transformed
#'     to make the parameter vector acceptable to the model. This is
#'     necessary whenever ABC sampling happens in a different space
#'     than the model parameter domain for whatever reason.
#' @param getScore scoring function.
#' @param delta the acceptance threshold.
#' @param nCores number of cores to use in parallel::mclapply() calls.
#' @return a filtered subset of acceptable parameter draws
checkFitWithPreviousExperiments <- function(modelName, draws, experiments, parMap=identity(), getScore, delta, nCores=detectCores()){
	numExperiments <- length(experiments)
	cat("-Checking fit with previous data\n")
	nDraws = dim(draws)[1]
	outputTimes_list <- list()
	outputFunctions_list <- list()
	for(k in 1:numExperiments){
		outputTimes_list <- c(outputTimes_list, replicate(nDraws, list(experiments[[k]][["outputTimes"]])))
		outputFunctions_list <- c(outputFunctions_list,replicate(nDraws, list(experiments[[k]][["outputFunction"]])))
	}

	output_yy <- runModel(experiments, modelName, t(draws), parMap, nCores)
	scores <- mclapply(seq(length(output_yy)),function(k) getScore(output_yy[[k]], experiments[[((k-1) %/% nDraws)+1]][["outputValues"]]), mc.preschedule = FALSE, mc.cores = nCores)
	scores <- unlist(scores)
	dim(scores) <- c(nDraws,numExperiments)
	acceptable <- apply(scores <= delta,1,all)
	stopifnot(length(acceptable)==nDraws)
	if (any(acceptable)){
		draws <- draws[acceptable,]
		nPickedDraws <- nrow(draws)
		nonFits <- nDraws - nPickedDraws;
		cat("-- ", nonFits, " samples  did not fit previous datasets")
	} else {
		print(scores)
		warning("none of the draws have been accepted.")
	}
	return(draws)
}
