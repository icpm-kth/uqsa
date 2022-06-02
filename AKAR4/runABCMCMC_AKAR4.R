## source('../UQ/copulaFunctions.R')
## source('../UQ/runModel.R')
## source('../UQ/PreCalibration.R')
## source('../UQ/ABCMCMCFunctions.R')
## source('../UQ/ScoringFunction.R')

#library(parallel)
#library(VineCopula)
#library(MASS)
#library(ks)
#library(R.utils)
library(UQ)
library(rgsl)
library(SBtabVFGEN)

modelName <- "AKAR4"
SBtabDir <- getwd()
model = import_from_SBtab(SBtabDir)

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
ns <- 1000 # no of samples required from each ABC-MCMC chain 
npc <- 500 # pre-calibration 

# Define ABC-MCMC Settings
p <- 0.01     # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values
nChains <- 1 # For the ABC-MCMC Process: Nr of Parallel Chains; 
delta <- 0.01 

# Define the number of Cores for the parallelization
nCores <- 1

set.seed(2022)

# Define the score function to compare simulated data with experimental data
maxVal <- max(unlist(lapply(experiments, function(x) max(x[["outputValues"]]))))
minVal <- min(unlist(lapply(experiments, function(x) min(x[["outputValues"]]))))

getScore  <- function(yy_sim, yy_exp){
  
  yy_sim <- (yy_sim-0)/(0.2-0.0)
  ifelse(!is.na(yy_exp), yy_exp <- (yy_exp-minVal)/(maxVal-minVal), Inf)
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
  out1 <- preCalibration(experiments[expInd], modelName, parVal, parIdx, npc, copula, U, Z, getScore, nCores = nCores*nChains, environment)
  
  sfactor <- 0.1 # scaling factor 
  
  ## Get Starting Parameters from Pre-Calibration
  out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
  Sigma <- out2$Sigma
  startPar <- out2$startPar
  
  ## Run ABC-MCMC Sampling
  #cat(sprintf("-Running MCMC chains \n"))
  #cl <- makeCluster(nChains, type="FORK", outfile=paste("out-log_",modelName,"_",environment,"_ns",ns,".txt",sep=""))
 
  # run outer loop
  #draws <- parLapply(cl, 1:nChains, function(k) ABCMCMC(experiments[expInd], modelName, startPar[k,], parIdx, parVal, ns, Sigma, delta, U, Z, Y, copula, ll, ul, getScore, nCores, environment))
  #stopCluster(cl)
  #rm(cl)
  draws <- ABCMCMC(experiments[expInd], modelName, startPar, parIdx, parVal, ns, Sigma, delta, U, Z, Y, copula, ll, ul, getScore, nCores, environment)
  # put together
  #draws <- do.call("rbind", draws)
  pick <- !apply(draws, 1, function(rw) all(rw==0))
  draws <- draws[pick,]
  draws <- checkFitWithPreviousExperiments(i, experimentsIndices, modelName, draws, experiments, parVal, parIdx, getScore, delta, environment, nCores, nChains)
  
  # Save Resulting Samples to MATLAB and R files.
  cat("-Saving sample \n")
  outFile <- paste(unlist(experimentsIndices[1:i]), collapse="_")
  timeStr <- Sys.time()
  timeStr <- gsub(":","_", timeStr)
  timeStr <- gsub(" ","_", timeStr)
  outFileR <- paste0("../PosteriorSamples/Draws",modelName,"_",environment,"_ns",ns,"_npc",npc,"_",outFile,timeStr,".RData",collapse="_")
  outFileM <- paste0("../PosteriorSamples/Draws",modelName,"_",environment,"_ns",ns,"_npc",npc,"_",outFile,timeStr,".mat",collapse="_")
  save(draws, parNames, file=outFileR)
  #writeMat(outFileM, samples=10^draws)
  
}
end_time = Sys.time()
time_ = end_time - start_time

#### PLOT RESTULTS FOR AKAR4 
# par(mfrow=c(2,3))
# for(i in 1:3){
#   hist(draws[,i], main=parNames[i], xlab = "Value in log scale")
# }
# combinePar <- list(c(1,2), c(1,3), c(2,3))
# for(i in combinePar){
#   plot(draws[,i[1]], draws[,i[2]], xlab = parNames[i[1]], ylab = parNames[i[2]])
# }
# 
# library(plotly)
# df = as.data.frame(draws)
# colnames(df) <- parNames
# plot_ly(dat = df, x = ~kf_C_AKAR4, y = ~kb_C_AKAR4, z = ~kcat_AKARp, type="scatter3d", mode="markers", marker=list(size = 1, color = "red"))
