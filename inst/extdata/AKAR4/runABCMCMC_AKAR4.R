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
ll <- c(parVal/defRange)
ul <- c(parVal*defRange)
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 50000 # no of samples required from each ABC-MCMC chain
npc <- 50000 # pre-calibration

# Define ABC-MCMC Settings
p <- 0.01 # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values

delta <- 0.002
set.seed(2022)

## library(GauPro)

## gaussianProcessPrediction <- function(e) {
##   gp <- GauPro(X=e$outputTimes, Z=e$outputValues, parallel=FALSE)
##   e$outputValues<-gp$predict(e$outputTimes)
##   return(e)
## }
## experiments <- lapply(experiments, gaussianProcessPrediction)

maxVal <- 184
minVal <- 108
distanceMeasure <- function(funcSim, dataExpr, dataErr = 1.0){
  distance <- mean(((funcSim-dataExpr$AKAR4pOUT)/max(dataExpr$AKAR4pOUT))^2,na.rm=TRUE)
  return(distance)
}


nCores <- parallel::detectCores()
nChains <- nCores %/% 2

# Loop through the Different Experimental Settings
start_time = Sys.time()

# work packages
chunks <- list(c(1,2),3)

for (i in seq(length(chunks))){
	expInd <- chunks[[i]]
	simulate <- simulator.c(experiments[expInd],modelName,parMap)
	Obj <- makeObjective(experiments[expInd],modelName,distanceMeasure,parMap,simulate)
	cat("#####Starting run for Experiments ", expInd, "######\n")
	## If First Experimental Setting, Create an Independente Colupla
	if(i==1){
		cat(sprintf("- Starting with uniform prior \n"))
		priorPDF <- dUniformPrior(ll, ul)
		rprior <- rUniformPrior(ll, ul)
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		cat(sprintf("- Fitting Copula based on previous MCMC runs\n"))
		C <- fitCopula(mcmc$draws)
		priorPDF <- dCopulaPrior(C)
		rprior <- rCopulaPrior(C)
	}
	## Run Pre-Calibration Sampling
	cat(sprintf("- Precalibration \n"))

	time_pC <- Sys.time()
	options(mc.cores = nCores)
	pC <- preCalibration(Obj, npc, rprior, rep=3)
	time_pC <- Sys.time() - time_pC
	cat(sprintf("- time spent on precalibration: \n"))
	print(time_pC)

	## Get Starting Parameters from Pre-Calibration
	M <- getMCMCPar(pC$prePar, pC$preDelta, delta, num=nChains)
	## Run ABC-MCMC Sampling
	cat(sprintf("- Running MCMC\n"))
	time_ABC <- Sys.time()
	options(mc.cores=nCores %/% nChains)
	mcmc <- ABCMCMC(Obj, M$startPar, ns, M$Sigma, delta, priorPDF)
	time_ABC <- Sys.time() - time_ABC
	cat(sprintf("- time spent on ABC-MCMC: \n"))
	print(time_ABC)
	if (i>1){
		simulate <- simulator.c(experiments,modelName,parMap)
		Obj <- makeObjective(experiments,modelName,distanceMeasure,parMap,simulate)
		mcmc$draws <- checkFitWithPreviousExperiments(mcmc$draws, Obj, delta)
	}
	# Save Resulting Samples to MATLAB and R files.

	## this section makes a little sensitivity plot:
	y<-simulate(t(mcmc$draws))
	f<-aperm(y[[3]]$func[1,,]) # aperm makes the sample-index (3rd) the first index of f, default permutation
	S<-sensitivity(mcmc$draws,f)
	S[1,]<-0 # the first index of S is time, and initially sensitivity is 0
	cuS<-t(apply(S,1,cumsum))
	plot.new()
	tm<-experiments[[3]]$outputTimes
	plot(tm,cuS[,3],type="l")
	for (si in dim(S)[2]:1){
		polygon(c(tm,rev(tm)),c(cuS[,si],numeric(length(tm))),col=si+1)
	}
}
end_time = Sys.time()
time_ = end_time - start_time

#### PLOT RESTULTS FOR AKAR4
par(mfrow=c(2,3))
for(i in 1:3){
	 hist(mcmc$draws[,i], main=parNames[i], xlab = "Value in log scale")
}
combinePar <- list(c(1,2), c(1,3), c(2,3))
for(i in combinePar){
	 plot(mcmc$draws[,i[1]], mcmc$draws[,i[2]], xlab = parNames[i[1]], ylab = parNames[i[2]])
}


library(reshape2)
library(ggplot2)

for(i in 1:3){
  experiment <- experiments[i]
  stopifnot(length(experiment)==1)
  output_yy <- simulate(t(mcmc$draws))
  df_ <- mclapply(1:dim(output_yy[[1]][["func"]])[3], function(i) as.data.frame(x = list(output_yy[[1]][["func"]][1,,i]/0.2,experiment[[1]][['outputTimes']]), col.names = c("y","t")))
  df__ <- melt(df_,id=c("t","y"))
  yy_exp <- (experiments[[i]][["outputValues"]]-minVal)/(maxVal-minVal)
  dfExpA<- data.frame(t=experiment[[1]][["outputTimes"]], y=yy_exp)
  ggp <- ggplot(df__,aes(x=t, y=y, group=L1))+
    geom_line(color="blue")+
    geom_point(data=dfExpA, aes(x=t, y=y), inherit.aes=FALSE)
  show(ggp)
}

