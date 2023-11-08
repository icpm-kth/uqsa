#!/usr/bin/env Rscript
require(rgsl)
require(SBtabVFGEN)
library(uqsa)
library(parallel)

## the config file contains default values and standard definitions
## it loads the model and creates simulation functions for it
source(uqsa_example("AKAP79",pat="^config.*R$"))
cat("these are the loaded variables:",ls())

start_time = Sys.time()

CA <- commandArgs(trailingOnly=TRUE)
if (length(CA)>0){
	KEY <- subset(CA,c(TRUE,FALSE))
	VAL <- subset(CA,c(FALSE,TRUE))
	for (i in seq(length(KEY))){
		k <- tolower(KEY[i])
		v <- as.numeric(VAL[i])
		switch(k,
			ns={ns <- v},
			npc={npc <- v},
			delta={delta <- v},
			{cat(sprintf("Option %s with value %g unknown.\n",k,v))}
		)
	}
}

#acceptanceProbability <- NULL

message("- Initial Prior: uniform product distribution")

#rprior <- rUniformPrior(ll, ul)
rprior <- rNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)
# with this choice of sd, each component has 98.8% probability
# of being in its interval [ll,ul],
# and the vector has (98.8%)^27 = 71.2% probability of being in the hyperrectangle
# defined by ll and ul

#dprior <- dUniformPrior(ll, ul)
dprior <- dNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)

## Run Pre-Calibration Sampling
message("- Precalibration")
start_time_preCalibration <- Sys.time()
options(mc.cores=parallel::detectCores())
pC <- preCalibration(objectiveFunction, npc, rprior, rep = 5)
cat("\nPreCalibration:")
print(Sys.time()-start_time_preCalibration)

## Get Starting Parameters from Pre-Calibration
M <- getMCMCPar(pC$prePar, pC$preDelta, delta=delta, num = nChains)
options(mc.cores=nCoresPerChain)
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
clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "dprior"))
out_ABCMCMC <- parLapply(cl, 1:nChains, function(j) ABCMCMC(objectiveFunction=objectiveFunction, M$startPar[,j], ns, M$Sigma, delta, dprior, acceptanceProbability=NULL))
stopCluster(cl)

ABCMCMCoutput <- do.call(Map, c(rbind,out_ABCMCMC))
ABCMCMCoutput$scores <- as.numeric(t(ABCMCMCoutput$scores))
end_time = Sys.time()
time_ = end_time - start_time_ABC
print(time_)
cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)

# Save Resulting Samples to MATLAB and R files.
save_sample(ABCMCMCoutput)

end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)


