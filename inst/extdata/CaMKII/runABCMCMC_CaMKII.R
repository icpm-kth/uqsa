require(rgsl)
require(SBtabVFGEN)
require(parallel)
library(uqsa)

nChains <- 4
options(mc.cores=parallel::detectCores() %/% nChains)

model.tsv <- uqsa_example("CaMKII",full.names=TRUE)
model.tab <- sbtab_from_tsv(model.tsv)

# source all R functions for this model, this also loads a variable called "model"
source(uqsa_example("CaMKII",pat="^CaMKIIs[.]R$",full.names=TRUE))

experiments <- sbtab.data(model.tab)

# we take the model name as inferred from the sbtab document:
modelName <- checkModel(comment(model.tab), uqsa_example("CaMKII",pat="_gvf[.]c$"))
## Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
<<<<<<< HEAD
ns <- 500 # Size of the sub-sample from each chain
npc <- 1000 # pre-calibration sample size
=======
ns <- 5000 # Size of the sub-sample from each chain
npc <- 5000 # pre-calibration sample size
>>>>>>> e16f55d (re-created AKAP79.R)

## model is a list variable defined in CaMKIIs.R, model$par() is the
## CaMKII_default() function from the same file.  But, model$par() is
## a more generic name that is defined for any model made by the same
## tool (RPN-derivative/sh/ode.sh).  These functions takes the
## "!Scale" mentioned in SBtab into account and always return the
## parameters in linear scale. So, if the SBtab gives us the
## logarithmic values, this function has the converted values:
numPar <- nrow(model.tab$Parameter)
parVal <- model$par()[1:numPar]
## There is one caveat: in our SBtab files we make a distinction
## between *model parameters* and *input parameters*. The input
## parameters represent something we did to the system during the
## experiment (known values), while the model parameters represent
## internal properties of the system (unknown). The ODE model file
## makes no such distinction (it just has ODE parameters). So, the
## parameters the model accepts are more numerous than what we need to
## estimate using ABC. We assume the internal model parameters to come
## first and input parameters last. During sampling (or optimization)
## we only want to consider the not-so-well-known internal parameters:
## that's why we do [1:numPar].

# we sample in logarithmic space, so the model gets a transformation
# function to compensate:
parMap <- function(parABC){
	return(10^parABC)
}

## scale to determine prior values
defRange <- 100

## Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- log10(parVal/defRange)
ul <- log10(parVal*defRange)
## parameter 23 should not vary quite as much:
ll[23] <- log10(parVal[23]/10)
ul[23] <- log10(parVal[23]*10)

## Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
experimentsIndices <- list(
 which(startsWith(names(experiments),"E0")),
 which(startsWith(names(experiments),"E1")),
 which(startsWith(names(experiments),"E2")),
 which(startsWith(names(experiments),"E3")),
 which(startsWith(names(experiments),"E4")),
 which(startsWith(names(experiments),"E5"))
)

experimentsIndices <- list(1:99)

# Define ABC-MCMC Settings
delta <- 1.0

set.seed(7619201)

distanceMeasure <- function(funcSim, dataExpr=Inf, dataErr=Inf){
	if (all(is.finite(funcSim))){
		distance <- mean(abs((funcSim-as.matrix(dataExpr))/as.matrix(dataErr)), na.rm=TRUE)
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

start_time = Sys.time()
for (i in 1:length(experimentsIndices)){
	expInd <- experimentsIndices[[i]]
	simulate <- simulator.c(experiments[expInd], modelName, parMap)
	objectiveFunction <- makeObjective(experiments[expInd], modelName, distanceMeasure, parMap, simulate)

	cat("#####Starting run for Experiments ", expInd, "######\n")
	if(i==1){
		message("- Initial Prior: uniform product distribution")
		rprior <- rUniformPrior(ll, ul)
		dprior <- dUniformPrior(ll, ul)
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		start_time_fitCopula <- Sys.time()
		message("- New Prior: fitting Copula based on previous MCMC runs")
		C <- fitCopula(ABCMCMCoutput$draws)
		rprior <- rCopulaPrior(C)
		dprior <- dCopulaPrior(C)
		cat("\nFitting copula:")
		print(Sys.time()-start_time_fitCopula)
	}

	## Run Pre-Calibration Sampling
	message("- Precalibration")
	start_time_preCalibration <- Sys.time()

	## Pre-Calibration uses only mclapply calls, we set up no cluster for this.
	options(mc.cores=parallel::detectCores())
	pC <- preCalibration(objectiveFunction, npc, rprior, rep=5)
	cat("\nPreCalibration:")
	print(Sys.time()-start_time_preCalibration)

	## Get Starting Parameters from Pre-Calibration
	M <- getMCMCPar(pC$prePar, pC$preDelta, delta=delta, num = nChains)
	for(j in 1 : nChains){
		stopifnot(dprior(M$startPar[,j])>0)
		cat("Chain", j, "\n")
		cat("\tMin distance of starting parameter for chain",j," = ", min(objectiveFunction(M$startPar[,j])),"\n")
		cat("\tMean distance of starting parameter for chain",j," = ", mean(objectiveFunction(M$startPar[,j])),"\n")
		cat("\tMax distance of starting parameter for chain",j," = ", max(objectiveFunction(M$startPar[,j])),"\n")
	}
	## Run ABC-MCMC Sampling
	cat(sprintf("-Running MCMC chains \n"))
	start_time_ABC = Sys.time()

	## here we set up a cluster and start n parallel MCMC chains
	## so, we reduce the number of cores per chain (the chains are on the same computing node)
	options(mc.cores=parallel::detectCores() %/% nChains)
	cl <- makeForkCluster(nChains, outfile="outputMessagesABCMCMC.txt")
	clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "dprior"))
	out_ABCMCMC <- parLapply(cl,
		1:nChains,
		function(j) {
			tryCatch(
				ABCMCMC(objectiveFunction, M$startPar[,j], ns, M$Sigma, delta, dprior, objectiveFunction),
				error=function(cond) {message("ABCMCMC crashed"); print(M); print(j); return(NULL)}
			)
		}
	)
	stopCluster(cl)

	ABCMCMCoutput <- do.call(Map, c(rbind,out_ABCMCMC))
	if (!is.matrix(ABCMCMCoutput$draws)){
		stop("All chains got stuck.")
	}
	if (!is.null(ABCMCMCoutput$scores)){
		ABCMCMCoutput$scores <- as.numeric(t(ABCMCMCoutput$scores))
	}
	end_time = Sys.time()
	time_ = end_time - start_time_ABC
	cat("\nABCMCMC for experimental set",i,":")
	print(time_)
	cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
	cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)

	cat("\n-Saving sample \n")
	save_sample(experimentsIndices[[1]],ABCMCMCoutput)
	cat("memory usage: \n")
	print(gc())
}
end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)


