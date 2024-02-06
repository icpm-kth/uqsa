require(rgsl)
require(SBtabVFGEN)
require(parallel)
library(uqsa)

nCores <- parallel::detectCores()
nChains <- max(1,nCores %/% 2)
coresPerChain <- max(1,nCores %/% nChains)

model.tsv <- uqsa_example("CaMKII",full.names=TRUE)
model.tab <- sbtab_from_tsv(model.tsv)

# source all R functions for this model, this also loads a variable called "model"
source(uqsa_example("CaMKII",pat="^CaMKIIs[.]R$",full.names=TRUE))

experiments <- sbtab.data(model.tab)

# we take the model name as inferred from the sbtab document:
modelName <- checkModel(comment(model.tab), uqsa_example("CaMKII",pat="_gvf[.]c$"))
## Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 500 # Size of the sub-sample from each chain
npc <- 100000 # pre-calibration sample size

numPar <- nrow(model.tab$Parameter)
parVal <- model$par()[1:numPar]
parMap <- function(parABC){
	return(10^parABC)
}

## scale to determine prior values
defRange <- 10

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

# Define ABC-MCMC Settings
delta <- 3.0

set.seed(7619201)

distanceMeasure <- function(funcSim, dataExpr=Inf, dataErr=Inf){
	if (all(is.finite(funcSim))){
		distance <- mean(abs((as.numeric(funcSim)-as.numeric(dataExpr))/as.numeric(dataErr)), na.rm=TRUE)
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
rprior <- rUniformPrior(ll, ul)
dprior <- dUniformPrior(ll, ul)


for (i in 1:length(experimentsIndices)){
	expInd <- experimentsIndices[[i]]
	simulate <- simulator.c(experiments[expInd], modelName, parMap)
	objectiveFunction <- makeObjective(experiments[expInd], modelName, distanceMeasure, parMap, simulate)

	cat("#####Starting run for Experiments ", expInd, "######\n")

	## Run Pre-Calibration Sampling
	cat("- Precalibration\n")
	start_time_preCalibration <- Sys.time()
	options(mc.cores=nCores)
	pC <- preCalibration(objectiveFunction, npc, rprior, rep=3)
	cat("\tPreCalibration time:")
	print(Sys.time()-start_time_preCalibration)
	cat("\n")
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
	cat("-Running MCMC chains \n")
	start_time_ABC <-  Sys.time()

	## here we set up a cluster and start n parallel MCMC chains
	## so, we reduce the number of cores per chain (the chains are on the same computing node)
	options(mc.cores=coresPerChain)
	cl <- makeForkCluster(nChains, outfile="outputMessagesABCMCMC.txt")
	clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "dprior"))
	out_ABCMCMC <- parLapply( # parallel block
		cl,
		1:nChains,
		function(j) {
			ABCMCMC(objectiveFunction, M$startPar[,j,drop=FALSE], ns, M$Sigma, delta, dprior)
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
	cat("ABCMCMC for experimental set",i,":")
	print(time_)
	cat("\n")

	cat("Regularizations:", ABCMCMCoutput$nRegularizations,"\n")
	cat("Acceptance rate:", ABCMCMCoutput$acceptanceRate,"\n")

	start_time_fitCopula <- Sys.time()
	message("- New Prior: fitting Copula based on previous MCMC runs")
	C <- fitCopula(ABCMCMCoutput$draws)
	rprior <- rCopulaPrior(C)
	dprior <- dCopulaPrior(C)
	cat("\nFitting copula:")
	print(Sys.time()-start_time_fitCopula)

	cat("\n-Saving sample \n")
	save_sample(experimentsIndices[[1]],ABCMCMCoutput)
	cat("memory usage: \n")
	print(gc())
}
end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)
