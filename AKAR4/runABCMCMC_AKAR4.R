setwd('/Users/federicamilinanni/Documents/PhD/UQSA/uqsa/AKAR4')
source('../import_from_SBtab.R')
source('../UQ/copulaFunctions.R')
source('../UQ/runModel.R')
source('../UQ/PreCalibration.R')
source('../UQ/ABCMCMCFunctions.R')
source('../UQ/ScoringFunction.R')
library(parallel)
library(VineCopula)
library(MASS)
library(ks)
library(R.utils)
library(R.matlab)

SBtabDir <- getwd()
model = import_from_SBtab(SBtabDir)

modelName <- "AKAR4"
source(paste(SBtabDir,"/",modelName,".R",sep=""))

parVal <- model[["Parameter"]][["!DefaultValue"]]
parNames <- model[["Parameter"]][["!Name"]]

# load experiments
experiments <- import_experiments(modelName, SBtabDir)

# scale to determine prior values
defRange <- 1000

# Parameter set on which we perform uncertainty quantification. The other parameters are fixed at their default value.
parIdx <- 1:length(parVal)
parNames <- parNames[parIdx]

# Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- c(parVal/defRange)
ul <- c(parVal*defRange)
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale


# Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
experimentsIndices <- list(1,2,3)

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 10000 # no of samples required from each ABC-MCMC chain #WAS 1000
npc <- 50000 # pre-calibration  WAS 50.000

# Define ABC-MCMC Settings
p <- 0.01     # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values
nChains <- 12 # For the ABC-MCMC Process: Nr of Parallel Chains; 
delta <- 0.01 

# Define the number of Cores for the parallelization
nCores <- 10

set.seed(2022)

# Define the score function to compare simulated data with experimental data
maxVal<-numeric(0)
minVal<-numeric(0)
for(i in 1:3){
  maxVal[i] <- (max(experiments[[i]]$outputValues))
  minVal[i] <- (min(experiments[[i]]$outputValues))
}
maxVal<-max(maxVal)
minVal<-min(minVal)

getScore  <- function(yy_sim, yy_exp){
  
  yy_sim <- (yy_sim-0)/(0.2-0.0)
  yy_exp <- (yy_exp-minVal)/(maxVal-minVal)
  distance <- mean((yy_sim-yy_exp)^2)
  
  return(distance)
}

environment <- "C"


# Loop through the Different Experimental Settings
start_time = Sys.time()
for (i in 1:length(experimentsIndices)){
  
  expInd <- experimentsIndices[[i]]
  
  cat("#####Starting run for Experiments ", expInd, "######\n")
  
  ## If First Experimental Setting, Create an Independente Colupla
  if(i==1){
    cat(sprintf("-Fitting independent Copula \n"))
    out <- makeIndepCopula(ll, ul)
    ## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
  } else {
    cat(sprintf("-Fitting Copula based on previous MCMC runs\n"))
    out <- fitCopula(draws, ll, ul, nChains)
  }
  
  copula <- out$copula
  U <- out$U
  Z <- out$Z
  Y <- out$Y
  
  ## Run Pre-Calibration Sampling
  cat(sprintf("-Precalibration \n"))
  out1 <- preCalibration(experiments[expInd], modelName, parVal, parIdx, npc, copula, U, Z, getScore, nCores = nCores, environment)
  
  sfactor <- 0.1 # scaling factor 
  
  ## Get Starting Parameters from Pre-Calibration
  out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
  Sigma <- out2$Sigma
  startPar <- out2$startPar
  
  ## Run ABC-MCMC Sampling
  cat(sprintf("-Running MCMC chains \n"))
  cl <- makeCluster(nChains, outfile="out-log.txt")
  clusterEvalQ(cl, c(library(parallel), library(VineCopula), library(MASS), source('../UQ/ABCMCMCFunctions.R')))
  clusterExport(cl, list("runModel", "modelName", "getScore", "delta", "experiments", "parVal", "parIdx", "ns", "ll", "ul", "nCores", "environment","Sigma", "startPar", "expInd", "copula","Z", "U", "Y", "minVal", "maxVal"))
  
  # run outer loop
  draws <- parLapply(cl, 1:nChains, function(k) ABCMCMC(experiments[expInd], modelName, startPar[k,], parIdx, parVal, ns, Sigma, delta, U, Z, Y, copula, ll, ul, getScore, nCores, environment))
  stopCluster(cl)
  rm(cl)
  
  # put together
  draws <- do.call("rbind", draws)
  pick <- !apply(draws, 1, function(rw) all(rw==0))
  draws <- draws[pick,]
  if(i>1){
    for(j in 1:(i-1)){
      filtInd <- experimentsIndices[[j]]
      cat("-Checking fit with dataset", filtInd, "\n")
      nDraws = dim(draws)[1]
      
      
      tmp_list <- mclapply(experiments[filtInd], function(x) replicate(nDraws, c(parVal,x[["input"]])),  mc.preschedule = FALSE, mc.cores = nCores)
      params_inputs <- do.call(cbind, tmp_list)
      params_inputs[parIdx,] <- 10^t(draws)
      
      tmp_list <- mclapply(experiments[filtInd], function(x) replicate(nDraws, x[["initialState"]]),  mc.preschedule = FALSE, mc.cores = nCores)
      y0 <- do.call(cbind, tmp_list)
      
      outputTimes_list <- list()
      outputFunctions_list <- list()
      for(k in 1:length(filtInd)){
        outputTimes_list <- c(outputTimes_list, replicate(nDraws, list(experiments[[k]][["outputTimes"]])))
        outputFunctions_list <- c(outputFunctions_list, replicate(nDraws, list(experiments[[k]][["outputFunction"]])))
      }
      
      output_yy <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)
      scores <- mclapply(1:length(output_yy), function(k) getScore(output_yy[[k]], experiments[[filtInd[(k-1)%/%nDraws+1]]][["outputValues"]]), mc.preschedule = FALSE, mc.cores = nCores)
      scores <- unlist(scores)
      
      pick <- scores <= delta
      draws <- draws[pick,];
      nPickedDraws <- nrow(draws)
      nonFits <-  nDraws - nPickedDraws;
      cat("-- ", nonFits, " samples of posterior after datasets ", expInd, " did not fit dataset ", filtInd)
    }
  }
  
  # Save Resulting Samples to MATLAB and R files.
  cat("-Saving sample \n")
  outFile <- paste(unlist(experimentsIndices[1:i]), collapse="_")
  timeStr <- Sys.time()
  timeStr <- gsub(":","_", timeStr)
  timeStr <- gsub(" ","_", timeStr)
  outFileR <- paste0("DrawsExperiments_",modelName,outFile,timeStr,".RData",collapse="_")
  outFileM <- paste0("DrawsExperiments_",modelName,outFile,timeStr,".mat",collapse="_")
  save(draws, parNames, file=outFileR)
  writeMat(outFileM, samples=10^draws)
  
}
end_time = Sys.time()
time_ = end_time - start_time


