setwd('/Users/federicamilinanni/Documents/PhD/UQSA/uqsa')  
source('import_from_SBtab.R')
source('UQ/copulaFunctions.R')
source('../newUQSA/runModel_multiExper.R')
source('../newUQSA/newPreCalibration_multiExper.R')
source('../newUQSA/ABCMCMCFunctions2_multiExper.R')
source('../newUQSA/ScoringFunction.R')
library(parallel)
library(VineCopula)
library(MASS)

SBtabDir <- "/Users/federicamilinanni/Documents/PhD/UQSA/uqsa/AKAP79"
model = import_from_SBtab(SBtabDir)

modelName <- "AKAP79"
source(paste(SBtabDir,"/",modelName,".R",sep=""))

#compoundNames <- model[["Compound"]][["!Name"]]

parVal <- model[["Parameter"]][["!DefaultValue"]]
parNames <- model[["Parameter"]][["!Name"]]

# load experiments
experiments <- import_experiments(modelName, SBtabDir)

# scale to determine prior values
scale <- 1000


# Parameter set on which we perform uncertainty quantification. The other parameters are fixed at their default value.
parIdx <- 1:length(parVal)
parNames <- parNames[parIdx]

# prior
ll = parVal[parIdx]/scale
ul = parVal[parIdx]*scale
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale
#NOTE THAT JOÃO USED A SCALE OF 1.25, 1.5 AND 2 FOR SOME OF THE PARAMETERS

normalizeAKAR4p <- function(x) ((x/1.08 - 100) / ( 171.67 - 100)) * 0.2

#FIX EXPERIMENTS (THIS INFORMATION SHOULD BE WRITTEN DIRECTLY IN THE SBTAB FILES, 
#I.E. NOT ADMIT NA VALUES IN THE SBTABS
# HAVING THE CORRECT INITIAL CONDITIONS IN THE SBTABS
# ...)
for(i in 1:length(experiments)){
  experiments[[i]][["outputValues"]] <- normalizeAKAR4p(experiments[[i]][["outputValues"]])
}

experimentsIndices <- list(c(2,4,7),c(1,5,9,12),c(3,6), c(8,10,11,13,14,15))

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 10 # no of samples required from each ABC-MCMC chain #WAS 1000
npc <- 50 # pre-calibration  WAS 50.000

# Define ABC-MCMC Settings
p <- 0.01 # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values
nChains <- 19 # For the ABC-MCMC Process: Nr of Parallel Chains; 
delta <- 0.1 #JOÃO USES DIFFERENT VALUES FOR DIFFERENT EXPERIMENTS (IN SOME CASES, EVEN -1)

nCores <- 8

set.seed(7619201)

# Loop through the Different Experimental Settings
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
  
  #start_time = Sys.time()
  out1 <- preCalibration(experiments[expInd], modelName, parVal, parIdx, npc, copula, U, Z, nCores = nCores, environment = "C")
  #end_time = Sys.time()
  #timeC_preCal <- end_time - start_time
   
  sfactor <- 0.1 # scaling factor 
  
  ## Get Starting Parameters from Pre-Calibration
  out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
  Sigma <- out2$Sigma
  startPar <- out2$startPar
  
  ## Run ABC-MCMC Sampling
  cat(sprintf("-Running MCMC chains \n"))
  
  #start_time = Sys.time()
  draws <- mclapply(1:nChains, function(k) ABCMCMC(experiments[expInd], modelName, startPar[k,], parIdx, parVal, ns, Sigma, delta, U, Z, Y, copula, ll, ul, nCores, "C"), mc.preschedule = FALSE, mc.cores = nChains);
  #end_time = Sys.time()
  #timeC_ABCMCMC2 <- end_time - start_time
  
  }

