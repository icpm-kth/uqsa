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

preCalibration <- function(experiments, modelName, parDefVal, parIdx, npc, copula, U, Z, getScore, nCores, environment){
  
  numExperiments <- length(experiments)
  np <- length(parIdx)
  
  R <- RVineSim(npc, copula)
  prePar <- matrix(0, npc, np)
  for(i in 1:np){
    prePar[,i] = spline(Z[,i],U[,i],xout=R[,i])$y
  }
  
  tmp_list <- mclapply(experiments, function(x) replicate(npc, c(parDefVal,x[["input"]])),  mc.preschedule = FALSE, mc.cores = nCores)
  params_inputs <- do.call(cbind, tmp_list)
  params_inputs[parIdx,] <- 10^t(prePar)
  
  tmp_list <- mclapply(experiments, function(x) replicate(npc, x[["initialState"]]),  mc.preschedule = FALSE, mc.cores = nCores)
  y0 <- do.call(cbind, tmp_list)
  
  outputTimes_list <- list()
  outputFunctions_list <- list()
  for(i in 1:numExperiments){
    outputTimes_list <- c(outputTimes_list, replicate(npc, list(experiments[[i]][["outputTimes"]])))
    outputFunctions_list <- c(outputFunctions_list, replicate(npc, list(experiments[[i]][["outputFunction"]])))
  }
  
  invisible(capture.output(output_yy <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)))
  preDelta <- mclapply(1:length(output_yy), function(i) getScore(output_yy[[i]], experiments[[(i-1)%/%npc+1]][["outputValues"]]), mc.preschedule = FALSE, mc.cores = nCores)
  preDelta <- unlist(preDelta)
  
  #preDelta is a vector of length npc*numExperiments.
  #It is obtained using npc different parameter vectors, each of them tested on all the experiment.
  #In particular, parameter i (in 1:npc) was used in the generation of preDelta[i+j*npc] (j in 0:numExperiments-1)
  #Hence, to get an estimate of the delta that sums up the goodness of a certain parameter on the chosen experiments, we can use - for instance - the mean of squares
  preDelta <- sapply(1:npc, function(i) sum((preDelta[i+seq(0,npc*(numExperiments-1), npc)])^2)/numExperiments)
  
  return(list(preDelta=preDelta, prePar=prePar))
}


getMCMCPar <- function(prePar, preDelta, p, sfactor, delta, nChains){
  prePar <- prePar[!is.na(preDelta),]
  preDelta <- preDelta[!is.na(preDelta)]
  nk <- nrow(prePar)*p
  pick1  <- grep(TRUE, (preDelta <= delta)) # pick all pars that meet threshold
  pick2 <- order(preDelta, decreasing = F)[1:nk] # pick top p percent
  
  if(length(pick1)>length(pick2)){
    pick <- pick1
  }else{
    pick <- pick2
  }
  
  Scorr <- cor(prePar[pick,])*0.8 #tone down corrs
  diag(Scorr) <- 1
  sdv <- apply(prePar[pick,], 2, sd)
  Sigma <- sfactor * Scorr * tcrossprod(sdv)
  
  startPar <- prePar[sample(pick, nChains, replace = FALSE),]
  list(Sigma=Sigma, startPar=startPar)
}

