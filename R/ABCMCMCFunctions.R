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
#' @return a list containing a sample matrix and a vector of scores (values of delta for each sample)
ABCMCMC <- function(objectiveFunction=NULL, startPar, nSims, Sigma0, delta, dprior, acceptanceProbability=NULL){
  if(is.null(objectiveFunction) && is.null(acceptanceProbability)){
    error("Provide objectiveFunction or acceptanceProbability.")
  }
  
  cat("Started chain.\n")
  Sigma1 <- 0.25*diag(diag(Sigma0))
  curDelta <- Inf
  np <- length(startPar)
  curPar  <- startPar
  curDelta <- NA
  if(is.null(acceptanceProbability)){
    curDelta <- max(objectiveFunction(curPar))
    if(is.na(curDelta)){
      cat("*** [ABCMCMC] curDelta is NA. Replacing it with Inf ***\n")
      curDelta <- Inf
    }
  }
  curPrior <- dprior(curPar)
  draws <- matrix(NA, nSims,np)
  scores <- rep(NA, nSims)

  n <- 0
  acceptedSamples <- 0
  nRegularizations <- 0
  batchSize <- 100
  while (n/batchSize < nSims){
    if(n %% batchSize == 0 && acceptedSamples<0.0005*n){
      nRegularizations <- nRegularizations + 1
      if(nRegularizations >= 5){
        timeStr <- Sys.time()
        timeStr <- gsub(":","_", timeStr)
        timeStr <- gsub(" ","_", timeStr)
        save(draws, file = paste0("AbortedChainAfterRegularization_",timeStr,".RData"))
        warning(paste0("Stuck chain (nRegularizations = ", nRegularizations,")"))
        return(list(draws = c(), scores = c(), acceptanceRate = c(), nRegularizations = nRegularizations))
      }
      cat(paste0("Regularization of proposal covariance matrix (nRegularizations = ", nRegularizations,")"))

      Sigma0 <- solve(solve(Sigma0)+solve(0.1*norm(Sigma0)*diag(1,np,np)))
      Sigma1 <- 0.25*diag(diag(Sigma0))
      draws <- matrix(NA, nSims,np)
      scores <- rep(NA, nSims)
      n <- 0
      acceptedSamples <- 0
    }
    #else if (n %% 100 == 0 && acceptedSamples > 0.1*n) {
    #  delta <- delta * .9
    #}

    if(runif(1)<=0.95){
      canPar <- MASS::mvrnorm(n=1, curPar, Sigma0)
    }else{
      canPar <- MASS::mvrnorm(n=1, curPar, Sigma1)
    }

    if(is.null(acceptanceProbability)){
      out <- parUpdate(objectiveFunction, curPar, canPar, curDelta, curPrior, delta, dprior)
      curDelta <- out$curDelta
    }else{
      out <- parUpdate_ProbabilisticAcceptance(acceptanceProbability, curPar, canPar, curPrior, dprior)
    }
    
    curPar <- out$curPar
    curDelta <- out$curDelta
    curPrior <- out$curPrior
    acceptedSamples <- acceptedSamples + out$acceptance

    n <- n+1
    if(n %% batchSize == 0){
      draws[n/batchSize,]  <- curPar
      scores[n/batchSize] <- curDelta
    }
    
    if(n %% 10000 == 0){
      cat("n =", n)
      print(gc())
    }
  }
  cat("Finished chain.\n")
  return(list(draws = draws, scores = scores, acceptanceRate = acceptedSamples/n, nRegularizations = nRegularizations))
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
parUpdate <- function(objectiveFunction, curPar, canPar, curDelta, curPrior, delta, dprior){
  canDelta <- max(objectiveFunction(canPar))

  if(is.na(canDelta)){
    cat("\n*** [parUpdate] canDelta is NA. Replacing it with Inf ***")
    canDelta <- Inf
  }
  canPrior <- dprior(canPar)

  if (canDelta <= max(delta, curDelta)){
    acceptance <- (runif(1) <= canPrior/curPrior)
    if (acceptance){
      curDelta <- canDelta
      curPrior <- canPrior
      curPar <- canPar
    }
  } else {
	  # curPar, curDelta, and curPrior remain unchanged
	  acceptance <- FALSE
  }
  return(list(curPar=curPar, curDelta=curDelta, curPrior=curPrior, acceptance=acceptance))
}



parUpdate_ProbabilisticAcceptance <- function(acceptanceProbability, curPar, canPar, curPrior, dprior){
  canPrior <- dprior(canPar)
  acceptance <- (runif(1) <= acceptanceProbability(canPar)*min(1,canPrior/curPrior))
  if (acceptance){
    curPrior <- canPrior
    curPar <- canPar
  }
  return(list(curPar=curPar, curPrior=curPrior, acceptance=acceptance))
}

#' ABC acceptance of currently sampled values given old data (Prior)
#'
#' The prior probability density model using copulas and vines is not
#' perfect, so values sampled from an imperfect prior estimate can be
#' checked against old data.
#'
#' @export
#' @param draws matrix of sampled values (to be filtered).
#' @param objectiveFunction function that, given a (vectorial) parameter as input, simulated the model and outputs the distance between experimental data and data simulated from the model with the parameter provided in input
#' @param delta the acceptance threshold.
#' @return a filtered subset of acceptable parameter draws
checkFitWithPreviousExperiments <- function(draws, objectiveFunction, delta){
  cat("\n-Checking fit with previous data\n")
  nDraws = dim(draws)[1]

  scores <- objectiveFunction(t(draws))

  dim(scores) <- c(nDraws,length(scores)/nDraws)
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
