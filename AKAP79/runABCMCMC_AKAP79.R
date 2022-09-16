remotes::install_github("a-kramer/rgsl")
remotes::install_github("a-kramer/SBtabVFGEN")
library(rgsl)
library(SBtabVFGEN)
library(UQ)

SBtabDir <- getwd()
model = import_from_SBtab(SBtabDir)
#modelName <- checkModel(comment(model),paste0(comment(model),'.R'))
modelName <- checkModel(comment(model),paste0(comment(model),'_gvf.c'))

#source(paste(SBtabDir,"/",modelName,".R",sep=""))

parVal <- model[["Parameter"]][["!DefaultValue"]]
parNames <- model[["Parameter"]][["!Name"]]

# load experiments
experiments <- import_experiments(modelName, SBtabDir)

# scale to determine prior values
defRange <- 1000

# Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- c(parVal[1:19]/defRange, parVal[20]/1.9, parVal[21]/defRange, parVal[22:24]/1.25, parVal[25:26]/1.5, parVal[27]/2)
ul <- c(parVal[1:19]*defRange, parVal[20]*1.9, parVal[21]*defRange, parVal[22:24]*1.25, parVal[25:26]*1.5, parVal[27]*2)
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale


# Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
experimentsIndices <- c(3, 12, 18, 9, 2, 11, 17, 8, 1, 10, 16, 7)

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 50 # no of samples required from each ABC-MCMC chain
npc <- 5000 # pre-calibration

# Define ABC-MCMC Settings
delta <- 8 #0.01

# Define the number of Cores for the parallelization
nCores <- parallel::detectCores() %/% 2

set.seed(7619201)

# Define the score function to compare simulated data with experimental data
getScore	<- function(yy_sim, yy_exp, yy_expErr){
	yy_sim <- (yy_sim-0)/(0.2-0.0)
	ifelse(!is.na(yy_exp), yy_exp <- (yy_exp-100)/(171.67-100), Inf)
	distance <- mean(((yy_sim-yy_exp)/(yy_expErr/(171.67-100)))^2)
	
	#When output function is fixed:
	#distance <- mean((yy_sim-yy_exp)/(yy_expErr)^2)
	return(distance)
}

parMap <- function(parABC){
	return(10^parABC)
}
  
  
start_time = Sys.time()
for (i in 1:length(experimentsIndices)){
  
	expInd <- experimentsIndices[i]
	objectiveFunction <- makeObjective(experiments[expInd], modelName, getScore, parMap)
	
	cat("#####Starting run for Experiments ", expInd, "######\n")
	## If First Experimental Setting, Create an Independente Colupla
	if(i==1){
		message("- Initial Prior: uniform product distribution")
		rprior <- rUniformPrior(ll, ul)
		dprior <- dUniformPrior(ll, ul)
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		message("- New Prior: fitting Copula based on previous MCMC runs")
		C <- fitCopula(draws, nCores)
		rprior <- rCopulaPrior(C)
		dprior <- dCopulaPrior(C)
	}
	## Run Pre-Calibration Sampling
	message("- Precalibration")
	pC <- preCalibration(objectiveFunction, npc, rprior)
	## Get Starting Parameters from Pre-Calibration
	M <- getMCMCPar(pC$prePar, pC$preDelta, delta=delta)
	stopifnot(dprior(M$startPar)>0)
	## Run ABC-MCMC Sampling
	cat(sprintf("-Running MCMC chains \n"))
	# run outer loop
	out_ABCMCMC <- ABCMCMC(objectiveFunction, M$startPar, ns, M$Sigma, delta, dprior)
	
	draws <- out_ABCMCMC$draws
	scores <- out_ABCMCMC$scores
	acceptanceRate <- out_ABCMCMC$acceptanceRate
	nRegularizations <- out_ABCMCMC$nRegularizations
	
	if (i>1){
		precursors <- experimentsIndices[1:(i-1)]
		objectiveFunction <- makeObjective(experiments[precursors], modelName, getScore, parMap, nCores)
		draws <- checkFitWithPreviousExperiments(draws, objectiveFunction, delta)
	}
	# Save Resulting Samples to MATLAB and R files.
	cat("-Saving sample \n")
	outFile <- paste(experimentsIndices[1:i], collapse="_")
	timeStr <- Sys.time()
	timeStr <- gsub(":","_", timeStr)
	timeStr <- gsub(" ","_", timeStr)
	outFileR <- paste("../PosteriorSamples/Draws",modelName,"ns",ns,"npc",npc,outFile,timeStr,".RData",collapse="_",sep="_")
	#save(draws, parNames, file=outFileR)
	#if (require(R.matlab)){
	#	outFileM <- paste("../PosteriorSamples/Draws",modelName,"ns",ns,"npc",npc,outFile,timeStr,".mat",collapse="_",sep="_")
	#	writeMat(outFileM, samples=10^draws)
	#}
}
end_time = Sys.time()
time_ = end_time - start_time


