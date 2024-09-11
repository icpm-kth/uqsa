## this project:
library(uqsa)
## our other packages:
require(rgsl)
require(SBtabVFGEN)

model.tsv <- uqsa_example("AKAR4",full.names=TRUE)
model.tab <- sbtab_from_tsv(model.tsv) # SBtabVFGEN
source(uqsa_example("AKAR4",pat="^AKAR4[.]R$"))

modelName <- checkModel(comment(model.tab),uqsa_example("AKAR4",pat="_gvf[.]c$")) # SBtabVFGEN

numPar <- nrow(model.tab$Parameter)
parNames <- row.names(model.tab$Parameter)
parVal <- model$par()[1:numPar]

# load experiments
experiments <- sbtab.data(model.tab)

# scale to determine prior values
defRange <- 1000

# a function that tansforms the ABC variables to acceptable model
# parameters, re-indexing could also happen here
parMap <- function (parABC=0) {
	return(10^parABC)
}

# Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- log10(c(parVal/defRange))
ul <- log10(c(parVal*defRange))

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 500 # no of samples required from each ABC-MCMC chain
npc <- 5000 # pre-calibration

delta <- 0.002

nChains <- 4

set.seed(2024)

## library(GauPro)

## gaussianProcessPrediction <- function(e) {
##   gp <- GauPro(X=e$outputTimes, Z=e$outputValues, parallel=FALSE)
##   e$outputValues<-gp$predict(e$outputTimes)
##   return(e)
## }
## experiments <- lapply(experiments, gaussianProcessPrediction)

distanceMeasure <- function(funcSim, dataExpr, dataErr = 1.0){
  distance <- mean(((funcSim-dataExpr$AKAR4pOUT)/max(dataExpr$AKAR4pOUT))^2,na.rm=TRUE)
  return(distance)
}

save_sample <- function(ABCMCMCoutput){
  timeStr <- gsub("[ :]","_", Sys.time())
  if (!dir.exists("./PosteriorSamples")) {
    dir.create("./PosteriorSamples")
  }
  outFileR <- paste("./PosteriorSamples/Draws",modelName,"nChains",nChains,"ns",ns,"npc",npc,timeStr,".RData",collapse="_",sep="_")
  save(ABCMCMCoutput, file=outFileR)
}

nCores <- parallel::detectCores()
# Loop through the Different Experimental Settings
start_time = Sys.time()

# work packages
#chunks <- list(c(1,2),3)
chunks <- list(c(1,2,3))

for (i in seq(length(chunks))){
	expInd <- chunks[[i]]
	simulate <- simulator.c(experiments[expInd],modelName,parMap)
	objectiveFunction <- makeObjective(experiments[expInd], modelName, distanceMeasure, parMap, simulate)
	
	cat("#####Starting run for Experiments ", expInd, "######\n")
	## If First Experimental Setting, Create an Independente Colupla
	if(i==1){
		cat(sprintf("- Starting with uniform prior \n"))
		
	  # uniform prior
	  #dprior <- dUniformPrior(ll, ul)
		#rprior <- rUniformPrior(ll, ul)
		
	  # normal prior
		dprior <- dNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5) 
		rprior <- rNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5) 
		
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		cat(sprintf("- Fitting Copula based on previous MCMC runs\n"))
		C <- fitCopula(mcmc$draws)
		dprior <- dCopulaPrior(C)
		rprior <- rCopulaPrior(C)
	}
	
	
	## Run Pre-Calibration Sampling
	message("- Precalibration")
	start_time_preCalibration <- Sys.time()
	options(mc.cores=parallel::detectCores())
	pC <- preCalibration(objectiveFunction, npc, rprior, rep = 5)
	cat("\nPreCalibration:")
	print(Sys.time()-start_time_preCalibration)
	

	
	## Get Starting Parameters from Pre-Calibration
	M <- getMCMCPar(pC$prePar, pC$preDelta, delta=delta, num = nChains)
	options(mc.cores=parallel::detectCores() %/% nChains)
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
	out_ABCMCMC <- parLapply(cl, 1:nChains, function(j) ABCMCMC(objectiveFunction, M$startPar[,j], ns, M$Sigma, delta, dprior))
	stopCluster(cl)
	
	ABCMCMCoutput <- do.call(Map, c(rbind,out_ABCMCMC))
	ABCMCMCoutput$scores <- as.numeric(t(ABCMCMCoutput$scores))
	end_time = Sys.time()
	time_ = end_time - start_time_ABC
	print(time_)
	cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
	cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)
	

	# if (i>1){
	#   simulate <- simulator.c(experiments[expInd],modelName,parMap)
	#   objectiveFunction <- makeObjective(experiments[expInd], modelName, distanceMeasure, parMap, simulate)
	#   
	# 	simulate <- simulator.c(experiments,modelName,parMap)
	# 	Obj <- makeObjective(experiments,modelName,distanceMeasure,parMap,simulate)
	# 	mcmc$draws <- checkFitWithPreviousExperiments(mcmc$draws, Obj, delta)
	# }
	
	
	save_sample(ABCMCMCoutput)
	
	## this section makes a little sensitivity plot:
	y<-simulate(t(ABCMCMCoutput$draws))
	f<-aperm(y[[1]]$func[1,,]) # aperm makes the sample-index (3rd) the first index of f, default permutation
	S<-globalSensitivity(ABCMCMCoutput$draws,f)
	S[1,]<-0 # the first index of S is time, and initially sensitivity is 0
	cuS<-t(apply(S,1,cumsum))
	plot.new()
	tm<-experiments[[3]]$outputTimes
	plot(tm,cuS[,3],type="l")
	colors <- c(rgb(98,152,210,max = 255),
	            rgb(255,190,0, max=255),
	            rgb(232,106,88,max = 255))
	for (si in dim(S)[2]:1){
		polygon(c(tm,rev(tm)),c(cuS[,si],numeric(length(tm))),col=colors[si])
	}
}
end_time = Sys.time()
time_ = end_time - start_time


