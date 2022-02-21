source('import_from_SBtab.R')
source('UQ/copulaFunctions.R')
source('UQ/runModel.R')
source('UQ/PreCalibration.R')
source('UQ/ABCMCMCFunctions.R')
source('UQ/ScoringFunction.R')
library(parallel)
library(VineCopula)
library(MASS)

SBtabDir <- paste(getwd(),"/AKAP79",sep="")
model = import_from_SBtab(SBtabDir)

modelName <- "AKAP79"
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
ll <- c(parVal[1:19]/defRange, parVal[20]/1.9, parVal[21]/defRange, parVal[22:24]/1.25, par[25:26]/1.5, par[27]/2)
ul <- c(parVal[1:19]*defRange, parVal[20]*1.9, parVal[21]*defRange, parVal[22:24]*1.25, par[25:26]*1.5, par[27]*2)
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale


# Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
experimentsIndices <- list(3, 12, -1, 9, 2, 11, -1, 8, 1, 10, -1, 7)
#1uM cAMP is missing

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 10 # no of samples required from each ABC-MCMC chain #WAS 1000
npc <- 50 # pre-calibration  WAS 50.000

# Define ABC-MCMC Settings
p <- 0.01     # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values
nChains <- 19 # For the ABC-MCMC Process: Nr of Parallel Chains; 
delta <- 0.01 

# Define the number of Cores for the parallelization
nCores <- 20

set.seed(7619201)

# Define the score function to compare simulated data with experimental data
getScore  <- function(yy_sim, yy_exp){
  distance <- 100
  yy_sim <- (yy_sim-0)/(0.2-0.0)
  yy_exp <- (yy_exp-100)/(171.67-100)
  
  distance <- mean((yy_sim-yy_exp)^2)
  
  return(distance)
}

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
  out1 <- preCalibration(experiments[expInd], modelName, parVal, parIdx, npc, copula, U, Z, getScore, nCores = nCores, environment = "C")
  
  sfactor <- 0.1 # scaling factor 
  
  ## Get Starting Parameters from Pre-Calibration
  out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
  Sigma <- out2$Sigma
  startPar <- out2$startPar
  
  ## Run ABC-MCMC Sampling
  cat(sprintf("-Running MCMC chains \n"))
  
  draws <- mclapply(1:nChains, function(k) ABCMCMC(experiments[expInd], modelName, startPar[k,], parIdx, parVal, ns, Sigma, delta, U, Z, Y, copula, ll, ul, getScore, nCores, "C"), mc.preschedule = FALSE, mc.cores = nChains);
  
}

