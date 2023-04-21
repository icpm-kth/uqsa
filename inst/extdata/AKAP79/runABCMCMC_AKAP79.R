require(rgsl)
require(SBtabVFGEN)
library(uqsa)
library(parallel)

SBtabDir <- getwd()
model = import_from_SBtab(SBtabDir)
print(comment(model))
modelName <- checkModel(comment(model),paste0(SBtabDir,'/',comment(model),'_gvf.c'))

parVal <- model[["Parameter"]][["!DefaultValue"]]
names(parVal)<-model[["Parameter"]][["!Name"]]
parNames <- model[["Parameter"]][["!Name"]]

# load experiments
experiments <- import_experiments(modelName, SBtabDir)

parMap <- function(parABC){
	return(10^parABC)
}

# scale to determine prior values
defRange <- 1000

# Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- c(parVal[1:19]/defRange, parVal[20]/1.9, parVal[21]/defRange, parVal[22:24]/1.25, parVal[25:26]/1.5, parVal[27]/2)
ul <- c(parVal[1:19]*defRange, parVal[20]*1.9, parVal[21]*defRange, parVal[22:24]*1.25, parVal[25:26]*1.5, parVal[27]*2)
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale

# Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
experimentsIndices <- list(c(3, 12,18, 9, 2, 11, 17, 8, 1, 10, 16, 7))

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 25000 # Size of the sub-sample from each chain
npc <- 5000 # pre-calibration sample size
nChains <- 4
n <- ns*nChains

# Define ABC-MCMC Settings
delta <- 7 #0.01

# Define the number of Cores for the parallelization


nCores <- parallel::detectCores() %/% nChains

set.seed(7619201)

getScore	<- function(yy_sim, yy_exp=Inf, yy_expErr=Inf){
	distance <- mean(((yy_sim-yy_exp)/yy_expErr)^2, na.rm=TRUE)
	return(distance)
}

getAcceptanceProbability <- function(yy_sim, yy_exp, yy_expErr){
	yy_sim <- (yy_sim-0)/(0.2-0.0)
	ifelse(!is.na(yy_exp), yy_exp <- (yy_exp-100)/(171.67-100), Inf)
	yy_expErr <- yy_expErr/(171.67-100)

	return(exp(-sum((yy_sim-yy_exp)^2/(2*yy_expErr),na.rm = TRUE)))
}

start_time = Sys.time()
for (i in 1:length(experimentsIndices)){

	expInd <- experimentsIndices[[i]]
	objectiveFunction <- makeObjective(experiments[expInd], modelName, getScore, parMap)
	acceptanceProbability <- makeAcceptanceProbability(experiments[expInd], modelName, getAcceptanceProbability, parMap)

	cat("#####Starting run for Experiments ", expInd, "######\n")
	## If First Experimental Setting, Create an Independente Colupla
	if(i==1){
		message("- Initial Prior: uniform product distribution")
		rprior <- rUniformPrior(ll, ul)
		dprior <- dUniformPrior(ll, ul)
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		start_time_fitCopula <- Sys.time()
		message("- New Prior: fitting Copula based on previous MCMC runs")
		C <- fitCopula(draws, nCores)
		rprior <- rCopulaPrior(C)
		dprior <- dCopulaPrior(C)
		cat("\nFitting copula:")
		print(Sys.time()-start_time_fitCopula)
	}

	## Run Pre-Calibration Sampling
	message("- Precalibration")
	start_time_preCalibration <- Sys.time()
	pC <- preCalibration(objectiveFunction, npc, rprior)
	cat("\nPreCalibration:")
	print(Sys.time()-start_time_preCalibration)

	## Get Starting Parameters from Pre-Calibration
	M <- getMCMCPar(pC$prePar, pC$preDelta, delta=delta, num = nChains)
	M$startPar <- matrix(M$startPar, nChains)
	for(j in 1 : nChains){
		stopifnot(dprior(M$startPar[j,])>0)
	}

	## Run ABC-MCMC Sampling
	cat(sprintf("-Running MCMC chains \n"))
	start_time_ABC = Sys.time()
	cl <- makeForkCluster(nChains)
	clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "dprior", "acceptanceProbability"))
	out_ABCMCMC <- parLapply(cl, 1:nChains, function(j) ABCMCMC(objectiveFunction, M$startPar[j,], ns, M$Sigma, delta, dprior, acceptanceProbability))
	stopCluster(cl)
	draws <- c()
	scores <- c()
	acceptanceRate <- c()
	nRegularizations <- c()
	for(j in 1:nChains){
		draws <- rbind(draws, out_ABCMCMC[[j]]$draws)
		scores <- c(scores, out_ABCMCMC[[j]]$scores)
		acceptanceRate <- c(acceptanceRate, out_ABCMCMC[[j]]$acceptanceRate)
		nRegularizations <- c(nRegularizations, out_ABCMCMC[[j]]$nRegularizations)
	}
	end_time = Sys.time()
	time_ = end_time - start_time_ABC
	cat("\nABCMCMC for experimental set",i,":")
	print(time_)
	cat("\nRegularizations:", nRegularizations)
	cat("\nAcceptance rate:", acceptanceRate)
	# if (i>1){
	#   precursors <- experimentsIndices[1:(i-1)]
	#   objectiveFunction <- makeObjective(experiments[precursors], modelName, getScore, parMap, nCores)
	#   draws <- checkFitWithPreviousExperiments(draws, objectiveFunction, delta)
	# }

	cat("\nNumber of draws after fitting with previous experiments:",dim(draws)[1])

	# Save Resulting Samples to MATLAB and R files.
	cat("\n-Saving sample \n")
	outFile <- paste(experimentsIndices[[1]], collapse="_")
	timeStr <- Sys.time()
	timeStr <- gsub(":","_", timeStr)
	timeStr <- gsub(" ","_", timeStr)
	if (!dir.exists("./PosteriorSamples")) {
		dir.create("./PosteriorSamples")
	}

	outFileR <- paste("./PosteriorSamples/Draws",modelName,"nChains",nChains,"ns",ns,"npc",npc,outFile,timeStr,".RData",collapse="_",sep="_")
	save(draws, parNames, file=outFileR)
}
end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)


