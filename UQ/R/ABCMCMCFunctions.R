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
#' boundaries, this function will sample the parameters to
#' characterize the Bayesian posterior probability for the given problem
#'
#' Initially this function performs a similar job as an optimizer, and
#' then transitions to MCMC sampling.
#'
#' @export
#' @param experiments a list of experiments
#' @param modelName model functions will be assumed to have this prefix
#' @param startPar starting value for the parameter vector
#' @param parIdx re-mapping between the MCMC variables and model parameters
#' @param parDefVal default values for parameters, useful when some are not subject to estimation
#' @param nSims requested sample size
#' @param Sigma0 multivariate normal covariance of Markov chain transition kernel
#' @param delta ABC acceptance threshold
#' @param U sampled values for marginal distribution
#' @param Z (CDE) marginal cumulative probability density estimate
#' @param Y (PDE) marginal probability density estimate
#' @param copula An RVineMatrix() as returned by RVineStructureSelect
#' @param ll lower limit of model parameters
#' @param ul upper limit of model parameters
#' @param getScore a function(model_output,experiment) that returns scores
#' @param nCores setting for multicore package
#' @param environment (string) "R" or "C", selects solvers
#' @return a sample matrix
ABCMCMC <- function(experiments, modelName, startPar, parIdx, parDefVal, nSims, Sigma0, delta, U, Z, Y, copula, ll, ul, getScore, nCores, environment){

  cat("Started chain.\n")
  Sigma1 <- 0.25*diag(diag(Sigma0))
  curDelta <- Inf
  np <- length(startPar)
  scount <- 1

  curPar  <- startPar

  numExperiments <- length(experiments)

  currPar <- parDefVal
  currPar[parIdx] <- 10^(startPar)

  tmp_list <- mclapply(experiments,
                       function(exper) c(currPar, exper[["input"]]),
                       mc.preschedule = FALSE,
                       mc.cores = nCores)
  params_inputs <- do.call(cbind, tmp_list)

  tmp_list <- mclapply(experiments,
                       function(exper) exper[["initialState"]],
                       mc.preschedule = FALSE,
                       mc.cores = nCores)
  y0 <- do.call(cbind, tmp_list)

  outputTimes_list <- mclapply(experiments,
                               function(exper) exper[["outputTimes"]],
                               mc.preschedule = FALSE,
                               mc.cores = nCores)
  outputFunctions_list <- mclapply(experiments,
                                   function(exper) exper[["outputFunction"]],
                                   mc.preschedule = FALSE,
                                   mc.cores = nCores)

  out <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)

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

  curPrior <- dprior(curPar, U, Z, Y, copula, ll, ul)
  draws <- matrix(0, nSims,np)

  n <- 0

  while (n < nSims){
    if(runif(1)<=0.95){
      canPar <- mvrnorm(n=1, curPar, Sigma0)
    }else{
      canPar <- mvrnorm(n=1, curPar, Sigma1)
    }
    out <- parUpdate(experiments, modelName, parIdx, parDefVal, curPar, canPar, curDelta, curPrior, delta, U, Z, Y, copula, ll, ul, getScore, environment, nCores)
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

#' Prior Probability Density Value from Copula Fit
#'
#' This function calculates a probability density value of an argument value x (vector) given a copula
#' obtained from a sample from the prior (from preceding data sets).
#'
#' Returns PDF(x)
#'
#' @export
#' @param inx PDF argument (value of random variable)
#' @param U samples of marginal distributions
#' @param Z (CDE) cumulative density estimates of U
#' @param Y (PDE) probability density estimates of U
#' @param copula as returned by copula estimator
#' @param ll lower limit of inx
#' @param ul upper limit of inx
#' @return probability density function value
dprior <- function(inx, U, Z, Y, copula, ll, ul){
  if(all(!is.na(inx))){
    np <- length(inx)
    ed <- sapply(1:np, function(i) approx(U[,i], Z[,i], xout=inx[i])$y)
    mpdf <- sapply(1:np, function(i) approx(U[,i], Y[,i], xout=inx[i])$y)

    if(any(is.na(ed)) || any(is.na(mpdf))){ # outside of copula defined limits
      jpdf <- 0
    }else if(!(all(inx >=ll) && all(inx<=ul))){ # outside of prior
      jpdf <- 0
    }else{
      jpdf <- RVinePDF(ed, copula, verbose = TRUE)*prod(mpdf)
    }
  } else {
    jpdf <- 0
  }
  return(jpdf)
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
#' @param parIdx remapping index set for passing parameters to the
#'     model
#' @param parDefVal default values for the parameters (for cases in
#'     which some never change)
#' @param curPar current parameter values (as ABC samples them)
#' @param canPar candidate parameter values (for MCMC)
#' @param curDelta current distance between data and simulation, if
#'     the MCMC chain has not yet reached any point where this is
#'     below the threshold (delta), this can be accepted as the new
#'     current state for the chain.
#' @param curPrior current Prior values given curPar
#' @param delta distance threshold for ABC
#' @param U sampled marginal parameter values
#' @param Z (CDE) cumulative probability estimate of U
#' @param Y (PDE) probability density estimate values of U
#' @param copula as returned by copula estimators
#' @param ll lower limit of parameters
#' @param ul upper limit of parameters
#' @param getScore a scoring function
#' @param environment "C" selects GSL solvers, "R" selects the deSolve
#'     as backend
#' @param nCores number of cores to use in mclapply().
#' @return updated values for curPar, curDelta, and curPrior
parUpdate <- function(experiments, modelName, parIdx, parDefVal, curPar, canPar, curDelta, curPrior, delta, U, Z, Y, copula, ll, ul, getScore, environment, nCores){

  numExperiments <- length(experiments)

  par <- parDefVal
  par[parIdx] <- 10^(canPar)
  tmp_list <- mclapply(experiments,
                       function(exper) c(par, exper[["input"]]),
                       mc.preschedule = FALSE,
                       mc.cores = nCores)
  params_inputs <- do.call(cbind, tmp_list)

  tmp_list <- mclapply(experiments,
                       function(exper) exper[["initialState"]],
                       mc.preschedule = FALSE,
                       mc.cores = nCores)
  y0 <- do.call(cbind, tmp_list)

  outputTimes_list <- mclapply(experiments,
                               function(exper) exper[["outputTimes"]],
                               mc.preschedule = FALSE,
                               mc.cores = nCores)
  outputFunctions_list <- mclapply(experiments,
                                   function(exper) exper[["outputFunction"]],
                                   mc.preschedule = FALSE,
                                   mc.cores = nCores)

  invisible(capture.output(out <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)))

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
    canPrior <- dprior(canPar, U, Z, Y, copula, ll, ul)
  }

  if (canDelta <= max(delta,curDelta)){
    if (canPrior==0){
      h <- 0
    }else{
      h <- min(1,canPrior/curPrior)
    }

    if (runif(1) <= h){
      curDelta <- canDelta
      curPrior <- canPrior
      curPar <- canPar
    }
  }
  list(curPar=curPar, curDelta=curDelta, curPrior=curPrior)
}

#' ABC acceptance of currently sampled values given old data (Prior)
#'
#' The prior probability density model using copulas and vines is not
#' perfect, so values sampled from an imperfect prior estimate can be
#' checked against old data.
#'
#' @export
#' @param currentExpSet an index, all experiments before that index
#'     will be evaluated for scores and acceptance.
#' @param experimentsIndices can be used to filter simulation
#'     experiments (to exclude some), e.g. 2:4, to exclude 1.
#' @param modelName name (prefix), is used to find the file and model functions therein.
#' @param draws matrix of sampled values (to be filtered).
#' @param experiments a list of experiments (all of them, or up to currentExpSet).
#' @param parVal not used other than for size?
#' @param parIdx remapping index passed to runModel().
#' @param getScore scoring function.
#' @param delta the acceptance threshold.
#' @param environment passed to runModel(), selects solver.
#' @param nCores number of cores to use in mclapply() calls.
#' @param nChains number of parallel Markov chains (unused?).
checkFitWithPreviousExperiments <- function(currentExpSet, experimentsIndices, modelName, draws, experiments, parVal, parIdx, getScore, delta, environment, nCores, nChains){
  if(currentExpSet>1){
    for(j in 1:(currentExpSet-1)){
      filtInd <- experimentsIndices[[j]]
      cat("-Checking fit with dataset", filtInd, "\n")
      nDraws = dim(draws)[1]


      tmp_list <- mclapply(experiments[filtInd],
                           function(x) replicate(nDraws, c(parVal,x[["input"]])),
                           mc.preschedule = FALSE,
                           mc.cores = nCores*nChains)
      params_inputs <- do.call(cbind, tmp_list)
      params_inputs[parIdx,] <- 10^t(draws)

      tmp_list <- mclapply(experiments[filtInd],
                           function(x) replicate(nDraws, x[["initialState"]]),
                           mc.preschedule = FALSE,
                           mc.cores = nCores*nChains)
      y0 <- do.call(cbind, tmp_list)

      outputTimes_list <- list()
      outputFunctions_list <- list()
      for(k in 1:length(filtInd)){
        outputTimes_list <- c(outputTimes_list,
                              replicate(nDraws, list(experiments[[k]][["outputTimes"]])))
        outputFunctions_list <- c(outputFunctions_list,
                                  replicate(nDraws, list(experiments[[k]][["outputFunction"]])))
      }

      output_yy <- runModel(y0, modelName, params_inputs, outputTimes_list,
                            outputFunctions_list, environment, nCores*nChains)
      scores <- mclapply(1:length(output_yy),
                         function(k) getScore(output_yy[[k]],
                                              experiments[[filtInd[(k-1)%/%nDraws+1]]][["outputValues"]]),
                         mc.preschedule = FALSE,
                         mc.cores = nCores*nChains)
      scores <- unlist(scores)

      scores <-  scores[!is.na(scores)]

      pick <- scores <= delta
      draws <- draws[pick,];
      nPickedDraws <- nrow(draws)
      nonFits <-  nDraws - nPickedDraws;
      cat("-- ", nonFits, " samples of posterior after datasets ",
          expInd, " did not fit dataset ", filtInd)
    }
  }
  return(draws)
}
