require(rgsl)
require(SBtabVFGEN)
library(uqsa)
library(parallel)

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 1000 # Size of the sub-sample from each chain
npc <- 5000 # pre-calibration sample size
nChains <- 4
# Define ABC-MCMC Settings
delta <- 0.01

model.tsv <- uqsa_example("AKAP79")
model.tab <- sbtab_from_tsv(model.tsv)

# source all R functions for this model, this also loads a variable called "model"
source(uqsa_example("AKAP79",pat="^AKAP79[.]R$"))

experiments <- sbtab.data(model.tab)

print(comment(model.tab))
modelName <- checkModel(comment(model.tab),uqsa_example("AKAP79",pat="_gvf.c")))

numPar <- nrow(model.tab$Parameter)
parVal <- model$par()[1:numPar]
parNames <- row.names(model.tab$Parameter)

# load experiments
experiments <- import_experiments(model.tab)

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


# Define the number of Cores for the parallelization
nCores <- parallel::detectCores() %/% nChains

set.seed(7619201)

distanceMeasure <- function(funcSim, dataExpr=Inf, dataErr=Inf){
  if (all(is.finite(funcSim))){
    distance <- mean((1/71.67*(funcSim-as.matrix(dataExpr))/as.matrix(dataErr))^2, na.rm=TRUE)
  } else {
    distance <- Inf
  }
  return(distance)
}

save_sample <- function(ind,ABCMCMCoutput){
	I <- paste(ind, collapse="_")
	timeStr <- gsub("[ :]","_", Sys.time())
	if (!dir.exists("./PosteriorSamples")) {
		dir.create("./PosteriorSamples")
	}
	outFileR <- paste("./PosteriorSamples/Draws",modelName,"nChains",nChains,"ns",ns,"npc",npc,I,timeStr,".RData",collapse="_",sep="_")
	save(ABCMCMCoutput, file=outFileR)
}


getAcceptanceProbability <- function(yy_sim, yy_exp, yy_expErr){
  return(exp(-getScore(yy_sim,yy_exp,yy_expErr)/(delta)))

}

start_time = Sys.time()
for (i in 1:length(experimentsIndices)){

	expInd <- experimentsIndices[[i]]
	simulate <- simulator.c(experiments[expInd], modelName, parMap)
	objectiveFunction <- makeObjective(experiments[expInd], modelName, distanceMeasure, parMap, simulate)
	
	#acceptanceProbability <- makeAcceptanceProbability(experiments[expInd], modelName, getAcceptanceProbability, parMap,simulate)
  acceptanceProbability <- NULL
  
	cat("#####Starting run for Experiments ", expInd, "######\n")
	## If First Experimental Setting, Create an Independente Colupla
	if(i==1){
		message("- Initial Prior: uniform product distribution")
		#rprior <- rUniformPrior(ll, ul)
		rprior <- rNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5) 
		# with this choice of sd, each component has 98.8% probability
		# of being in its interval [ll,ul], 
		# and the vector has (98.8%)^27 = 71.2% probability of being in the hyperrectangle
		# defined by ll and ul
		#dprior <- dUniformPrior(ll, ul)
		dprior <- dNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		start_time_fitCopula <- Sys.time()
		message("- New Prior: fitting Copula based on previous MCMC runs")
		C <- fitCopula(ABCMCMCoutput$draws, nCores)
		rprior <- rCopulaPrior(C)
		dprior <- dCopulaPrior(C)
		cat("\nFitting copula:")
		print(Sys.time()-start_time_fitCopula)
	}

	## Run Pre-Calibration Sampling
	message("- Precalibration")
	start_time_preCalibration <- Sys.time()
	pC <- preCalibration(objectiveFunction, npc, rprior, rep = 5)
	cat("\nPreCalibration:")
	print(Sys.time()-start_time_preCalibration)

	## Get Starting Parameters from Pre-Calibration
	M <- getMCMCPar(pC$prePar, pC$preDelta, delta=delta, num = nChains)
	options(mc.cores=parallal::detectCores %/% nChains)
	for(j in 1 : nChains){
		stopifnot(all(dprior(M$startPar[,j])>0))
	  cat("Chain", j, "\n")
	  cat("\tMin distance of starting parameter for chain",j," = ", min(objectiveFunction(M$startPar[,j])),"\n")
	  cat("\tMean distance of starting parameter for chain",j," = ", mean(objectiveFunction(M$startPar[,j])),"\n")
	  cat("\tMax distance of starting parameter for chain",j," = ", max(objectiveFunction(M$startPar[,j])),"\n")
	}

	## Run ABC-MCMC Sampling
	cat(sprintf("-Running MCMC chains \n"))
	start_time_ABC = Sys.time()
	cl <- makeForkCluster(nChains)
	clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "dprior", "acceptanceProbability"))
	out_ABCMCMC <- parLapply(cl, 1:nChains, function(j) ABCMCMC(objectiveFunction, M$startPar[,j], ns, M$Sigma, delta, dprior, acceptanceProbability))
	stopCluster(cl)
	
	ABCMCMCoutput <- do.call(Map, c(rbind,out_ABCMCMC))
	ABCMCMCoutput$scores <- as.numeric(t(ABCMCMCoutput$scores))
	end_time = Sys.time()
	time_ = end_time - start_time_ABC
	cat("\nABCMCMC for experimental set",i,":")
	print(time_)
	cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
	cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)
	
	# Save Resulting Samples to MATLAB and R files.
	save_sample(expInd,ABCMCMCoutput)
}
end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)


