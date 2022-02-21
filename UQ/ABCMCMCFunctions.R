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
  
  tmp_list <- mclapply(experiments, function(exper) c(currPar, exper[["input"]]), mc.preschedule = FALSE, mc.cores = nCores)
  params_inputs <- do.call(cbind, tmp_list)
  
  tmp_list <- mclapply(experiments, function(exper) exper[["initialState"]],  mc.preschedule = FALSE, mc.cores = nCores)
  y0 <- do.call(cbind, tmp_list)
  
  outputTimes_list <- mclapply(experiments, function(exper) exper[["outputTimes"]], mc.preschedule = FALSE, mc.cores = nCores)
  outputFunctions_list <- mclapply(experiments, function(exper) exper[["outputFunction"]], mc.preschedule = FALSE, mc.cores = nCores)
  
  invisible(capture.output(out <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)))
  
  curDelta <- mclapply(1:length(out), function(i) getScore(out[[i]], experiments[[i]][["outputValues"]]), mc.preschedule = FALSE, mc.cores = nCores)
  curDelta <- unlist(curDelta)
  
  #Similarly to what we did in the preCalibration, we average the score obtained with (the same) startPar applied to all the simulations (corresponding to different experiments setup)
  #As in preCalibration, we can use - for instance - the sum of squares
  curDelta <- mean(curDelta^2)
  
  curPrior <- dprior(curPar, U, Z, Y, copula, ll, ul)
  draws <- matrix(0, nSims,np)

  n <- 0
  
  while (n < nSims){
    if(runif(1)<=0.95){
      canPar <- mvrnorm(n=1, curPar, Sigma0)
    }else{
      canPar <- mvrnorm(n=1, curPar, Sigma1)
    }
    out <- parUpdate(experiments, modelName, parIdx, parDefVal, curPar, canPar, curDelta, curPrior, delta, U, Z, Y, copula, ll, ul, environment, nCores)
    curPar <- out$curPar
    curDelta <- out$curDelta
    curPrior <- out$curPrior
    
    scount <- ifelse(!all(curPar == canPar), scount+1, 1)
    
    if (curDelta <= delta & all(curPar == canPar)){
      n=n+1
      draws[n,] = curPar
    }
    
    if(scount>500){ #terminate chain if stuck
      cat('Aborted chain.\n')
      return(draws)
    }		
  }
  cat("Finished chain.\n")
  return(draws)
}



dprior <- function(inx, U, Z, Y, copula, ll, ul){
  np <- length(inx)
  
  ed <- sapply(1:np, function(i) approx(U[,i], Z[,i], xout=inx[i])$y)
  mpdf <- sapply(1:np, function(i) approx(U[,i], Y[,i], xout=inx[i])$y)
  
  if(any(is.na(ed))|any(is.na(mpdf))){ # outside of copula defined limits
    jpdf <- 0
  }else if(!(all(inx >=ll) & all(inx<=ul))){ # outside of prior
    jpdf <- 0
  }else{
    jpdf <- RVinePDF(ed, copula, verbose = TRUE)*prod(mpdf)
  }
  
  return(jpdf)
}


parUpdate <- function(experiments, modelName, parIdx, parDefVal, curPar, canPar, curDelta, curPrior, delta, U, Z, Y, copula, ll, ul, getScore, environment, nCores){
  #browser()
  
  numExperiments <- length(experiments)
  
  par <- parDefVal
  par[parIdx] <- 10^(canPar)
  tmp_list <- mclapply(experiments, function(exper) c(par, exper[["input"]]), mc.preschedule = FALSE, mc.cores = nCores)
  params_inputs <- do.call(cbind, tmp_list)
  
  tmp_list <- mclapply(experiments, function(exper) exper[["initialState"]],  mc.preschedule = FALSE, mc.cores = nCores)
  y0 <- do.call(cbind, tmp_list)
  
  outputTimes_list <- mclapply(experiments, function(exper) exper[["outputTimes"]], mc.preschedule = FALSE, mc.cores = nCores)
  outputFunctions_list <- mclapply(experiments, function(exper) exper[["outputFunction"]], mc.preschedule = FALSE, mc.cores = nCores)
  
  invisible(capture.output(out <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)))
  
  if(is.null(out)){
    canDelta <- Inf
    canPrior <- 0
  }else{  	
    canDelta <- mclapply(1:length(out), function(i) getScore(out[[i]], experiments[[i]][["outputValues"]]), mc.preschedule = FALSE, mc.cores = nCores)
    canDelta <- unlist(canDelta)
    
    #Similarly to what we did in the preCalibration and in ABCMCMC, we average the score obtained with (the same) startPar applied to all the simulations (corresponding to different experiments setup)
    #As in preCalibration and ABCMCMC, we can use - for instance - the sum of squares
    canDelta <- mean(canDelta^2)
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
