require(rgsl)
require(SBtabVFGEN)
library(uqsa)
library(parallel)

model.tsv <- uqsa_example("CaMKII",full.names=TRUE)
model.tab <- sbtab_from_tsv(model.tsv)

# source all R functions for this model
source(uqsa_example("CaMKII",pat="^CaMKII.*R$",full.names=TRUE))

experiments <- sbtab.data(model.tab)

modelName <- checkModel(comment(model.tab), uqsa_example("CaMKII",pat="_gvf[.]c$"))

## this is a function from the CaMKIIs.R file,
## it takes the "!Scale" mentioned in SBtab into account:
parVal <- model$par()


parMap <- function(parABC){
	return(10^parABC)
}

## scale to determine prior values
defRange <- 1000

## Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- log10(parVal/defRange)
ul <- log10(parVal*defRange)

## Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
experimentsIndices <- list(
 which(startsWith(names(experiments),"E0")),
 which(startsWith(names(experiments),"E1")),
 which(startsWith(names(experiments),"E2")),
 which(startsWith(names(experiments),"E3")),
 which(startsWith(names(experiments),"E4")),
 which(startsWith(names(experiments),"E5"))
)


## Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 10000 # Size of the sub-sample from each chain
npc <- 50000 # pre-calibration sample size

# Define ABC-MCMC Settings
delta <- 0.01

# Define the number of Cores for the parallelization
nChains <- 4
nCores <- parallel::detectCores() %/% nChains

set.seed(7619201)

Score <- function(yy_sim, yy_exp=Inf, yy_expErr=Inf){
	distance <- mean(((yy_sim-yy_exp)/yy_expErr)^2, na.rm=TRUE)
	return(distance)
}

save_sample <- function(ind,ABCMCMCoutput){
	I <- paste(ind, collapse="_")
	timeStr <- Sys.time()
	timeStr <- gsub(":","_", timeStr)
	timeStr <- gsub(" ","_", timeStr)
	if (!dir.exists("./PosteriorSamples")) {
		dir.create("./PosteriorSamples")
	}
	outFileR <- paste("./PosteriorSamples/Draws",modelName,"nChains",nChains,"ns",ns,"npc",npc,I,timeStr,".RData",collapse="_",sep="_")
	save(ABCMCMCoutput, file=outFileR)
}


start_time = Sys.time()
for (i in 1:length(experimentsIndices)){
	expInd <- experimentsIndices[[i]]
	objectiveFunction <- makeObjective(experiments[expInd], modelName, Score, parMap)
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
	M$startPar <- matrix(M$startPar, nChains)
	for(j in 1 : nChains){
		stopifnot(dprior(M$startPar[j,])>0)
		cat("Chain", j, "\n")
		cat("\tMin distance of starting parameter for chain",j," = ", min(objectiveFunction(M$startPar[j,])),"\n")
		cat("\tMean distance of starting parameter for chain",j," = ", mean(objectiveFunction(M$startPar[j,])),"\n")
		cat("\tMax distance of starting parameter for chain",j," = ", max(objectiveFunction(M$startPar[j,])),"\n")
	}
	## Run ABC-MCMC Sampling
	cat(sprintf("-Running MCMC chains \n"))
	start_time_ABC = Sys.time()
	cl <- makeForkCluster(nChains)
	clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "dprior", "acceptanceProbability"))
	out_ABCMCMC <- parLapply(cl, 1:nChains, function(j) ABCMCMC(objectiveFunction, M$startPar[j,], ns, M$Sigma, delta, dprior, acceptanceProbability))
	stopCluster(cl)

	ABCMCMCoutput <- do.call(Map, c(rbind,out_ABCMCMC))
	ABCMCMCoutput$scores <- as.numeric(t(ABCMCMCoutput$scores))
	end_time = Sys.time()
	time_ = end_time - start_time_ABC
	cat("\nABCMCMC for experimental set",i,":")
	print(time_)
	cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
	cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)

	cat("\n-Saving sample \n")
	save_sample(experimentsIndices[[1]],ABCMCMCoutput)
}
end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)


