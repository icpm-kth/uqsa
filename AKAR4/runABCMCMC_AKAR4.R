## source('../UQ/copulaFunctions.R')
## source('../UQ/runModel.R')
## source('../UQ/PreCalibration.R')
## source('../UQ/ABCMCMCFunctions.R')
## source('../UQ/ScoringFunction.R')

#library(parallel)
#library(VineCopula)
#library(MASS)
#library(ks)
#library(R.utils)
library(UQ)
library(rgsl)
library(SBtabVFGEN)
SBtabDir <- getwd()
model = import_from_SBtab(SBtabDir)
print(comment(model)) # this should be the name of the model, if everything works
modelName <- checkModel(comment(model),"./AKAR4_gvf.c")
#modelName <- checkModel(comment(model),"./AKAR4.R")

source(paste(SBtabDir,"/",modelName,".R",sep=""))

parVal <- model[["Parameter"]][["!DefaultValue"]]
parNames <- model[["Parameter"]][["!Name"]]

# load experiments
experiments <- import_experiments(modelName, SBtabDir)

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
ns <- 1000 # no of samples required from each ABC-MCMC chain
npc <- 500 # pre-calibration

# Define ABC-MCMC Settings
p <- 0.01		 # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values

delta <- 0.01

set.seed(2022)

# Define the score function to compare simulated data with experimental data
maxVal <- max(unlist(lapply(experiments, function(x) max(x[["outputValues"]]))))
minVal <- min(unlist(lapply(experiments, function(x) min(x[["outputValues"]]))))

getScore	<- function(yy_sim, yy_exp){
	yy_sim <- (yy_sim-0)/(0.2-0.0)
	ifelse(!is.na(yy_exp), yy_exp <- (yy_exp-minVal)/(maxVal-minVal), Inf)
	distance <- mean((yy_sim-yy_exp)^2)
	return(distance)
}

Obj <- makeObjective(experiments,modelName,getScore,parMap)

# Loop through the Different Experimental Settings
start_time = Sys.time()

# work packages
chunks <- list(c(1,2),3)

for (i in seq(length(chunks))){
	expInd <- chunks[[i]]
	cat("#####Starting run for Experiments ", expInd, "######\n")
	## If First Experimental Setting, Create an Independente Colupla
	if(i==1){
		cat(sprintf("- Starting with uniform prior \n"))
		priorPDF <- dUniformPrior(ll, ul)
		rprior <- rUniformPrior(ll, ul)
		## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
	} else {
		cat(sprintf("- Fitting Copula based on previous MCMC runs\n"))
		C<-fitCopula(draws)
		priorPDF <- dCopulaPrior(C)
		rprior <- rCopulaPrior(C)
	}
	## Run Pre-Calibration Sampling
	cat(sprintf("- Precalibration \n"))
	par <- rprior(1)
	opt <- optim(par,Obj,lower=ll,upper=ul)
	startPar <- opt$par
	
	out1 <- preCalibration(experiments[expInd], modelName, parMap, npc, rprior, getScore)
	sfactor <- 0.1 # scaling factor
	## Get Starting Parameters from Pre-Calibration
	out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta)
	Sigma <- out2$Sigma
	#startPar <- out2$startPar
	## Run ABC-MCMC Sampling
	cat(sprintf("- Running MCMC\n"))
	draws <- ABCMCMC(experiments[expInd], modelName, startPar, parMap, ns, Sigma, delta, dprior=priorPDF, getScore)

	if (i>1){
	 draws <- checkFitWithPreviousExperiments(modelName, draws, experiments[unlist(chunks[1:i-1])], parMap, getScore, delta)
	}
	# Save Resulting Samples to MATLAB and R files.
	cat("-Saving sample \n")
	outFile <- paste(seq(1,i), collapse="_")
	timeStr <- Sys.time()
	timeStr <- gsub(":","_", timeStr)
	timeStr <- gsub(" ","_", timeStr)
	outFileR <- paste0("../PosteriorSamples/Draws",modelName,"_",basename(comment(modelName)),"_ns",ns,"_npc",npc,"_",outFile,timeStr,".RData",collapse="_")
	if (require("R.matlab")){
		outFileM <- paste0("../PosteriorSamples/Draws",modelName,"_",basename(comment(modelName)),"_ns",ns,"_npc",npc,"_",outFile,timeStr,".mat",collapse="_")
	}
	save(draws, parNames, file=outFileR)
	#writeMat(outFileM, samples=10^draws)
}
end_time = Sys.time()
time_ = end_time - start_time

#### PLOT RESTULTS FOR AKAR4
par(mfrow=c(2,3))
for(i in 1:3){
	 hist(draws[,i], main=parNames[i], xlab = "Value in log scale")
}
combinePar <- list(c(1,2), c(1,3), c(2,3))
for(i in combinePar){
	 plot(draws[,i[1]], draws[,i[2]], xlab = parNames[i[1]], ylab = parNames[i[2]])
}

# library(plotly)
# df = as.data.frame(draws)
# colnames(df) <- parNames
# plot_ly(dat = df, x = ~kf_C_AKAR4, y = ~kb_C_AKAR4, z = ~kcat_AKARp, type="scatter3d", mode="markers", marker=list(size = 1, color = "red"))
