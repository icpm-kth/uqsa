# ABCMCMC with stochastic simulation
library(GillespieSSA2)
library(uqsa)
library(rgsl)
library(SBtabVFGEN)
library(parallel)
library(pracma)

nChains <- 4

#options(mc.cores=parallel::detectCores() %/% nChains)
options(mc.cores=4)

SBtabDir <- getwd()
#model = import_from_SBtab(SBtabDir)
model.tsv <- uqsa_example("AKAR4",full.names=TRUE) 
source(uqsa_example("AKAR4",pat = "^AKAR4[.]R",full.names=TRUE)) # model is loaded
model.tab <- SBtabVFGEN::sbtab_from_tsv(model.tsv)
modelName <- checkModel(comment(model.tab),"./AKAR4_gvf.c")

parVal <- model.tab[["Parameter"]][["!DefaultValue"]]
parNames <- model.tab[["Parameter"]][["!Name"]]
names(parVal) <- parNames

# load experiments
#experiments <- import_experiments(modelName, SBtabDir)
experiments <- SBtabVFGEN::sbtab.data(model.tab)

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
ns <- 100 # no of samples required from each ABC-MCMC chain
npc <- 100 # pre-calibration

# Define ABC-MCMC Settings
p <- 0.01		 # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values

delta <- 0.01

set.seed(2022)

# Function to compute the score (distance) between experimental and simulated data
distanceMeasure <- function(funcSim, dataExpr, dataErr = 1.0){
  distance <- mean(((funcSim-dataExpr)/max(dataExpr))^2,na.rm=TRUE)
  return(distance)
}

# Create reactions and parameters for stochastic simulations with the GillespieSSA2 package
reactions <- importReactionsSSA(model.tab)

AvoNum <- 6.022e23
unit <- 1e-6
vol <- 4e-15
Phi <- AvoNum * vol * unit

compiled_reactions <- GillespieSSA2::compile_reactions(
  reactions = reactions,
  state_ids = model.tab$Compound[["!Name"]],
  params = c(parVal, Phi=Phi)
)


# Loop through the Different Experimental Settings
start_time = Sys.time()

# work packages
#chunks <- list(c(1,2),3)
chunks <- list(c(1,2,3))

i<-1
#for (i in seq(length(chunks))){
  expInd <- chunks[[i]]
  cat("#####Starting run for Experiments ", expInd, "######\n")
  objectiveFunction <- makeObjectiveSSA(experiments[expInd], model, parNames,distanceMeasure,parMap, Phi, compiled_reactions, nStochSim = 3)
  ## If First Experimental Setting, Create an Independente Colupla
  if(i==1){
    cat(sprintf("- Starting with uniform prior \n"))
    priorPDF <- dUniformPrior(ll, ul)
    rprior <- rUniformPrior(ll, ul)
    ## Otherwise, Take Copula from the Previous Exp Setting and Use as a Prior
  } else {
    cat(sprintf("- Fitting Copula based on previous MCMC runs\n"))
    C<-fitCopula(draws$draws)
    priorPDF <- dCopulaPrior(C)
    rprior <- rCopulaPrior(C)
  }
  ## Run Pre-Calibration Sampling
  cat(sprintf("- Precalibration \n"))
  
  time_pC <- Sys.time()
  pC <- preCalibration(objectiveFunction, npc, rprior,rep = 1)
  time_pC <- Sys.time() - time_pC
  cat(sprintf("- time for precalibration: \n"))
  print(time_pC)
  
  
  sfactor <- 0.1 # scaling factor
  ## Get Starting Parameters from Pre-Calibration
  M <- getMCMCPar(pC$prePar, pC$preDelta, p, sfactor, delta, num = nChains)
  Sigma <- M$Sigma
  startPar <- M$startPar
  for(j in 1 : nChains){
    cat("Chain", j, "\n")
    cat("\tMin distance of starting parameter for chain",j," = ", min(objectiveFunction(M$startPar[,j])),"\n")
    cat("\tMean distance of starting parameter for chain",j," = ", mean(objectiveFunction(M$startPar[,j])),"\n")
    cat("\tMax distance of starting parameter for chain",j," = ", max(objectiveFunction(M$startPar[,j])),"\n")
  }
  
  ## Run ABC-MCMC Sampling
  cat(sprintf("- Running MCMC\n"))
  time_ABC <- Sys.time()
  cl <- makeForkCluster(nChains, outfile="outputMessagesABCMCMC.txt")
  clusterExport(cl, c("objectiveFunction", "M", "ns", "delta", "priorPDF"))
  out_ABCMCMC <- parLapply(cl,
                           1:nChains,
                           function(j) {
                             tryCatch(
                               ABCMCMC(objectiveFunction, M$startPar[,j], ns, M$Sigma, delta, priorPDF),
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
  
  #draws <- ABCMCMC(objectiveFunction, startPar, ns, Sigma, delta, priorPDF)
  
  time_ABC <- Sys.time() - time_ABC
  cat(sprintf("- time for ABCMCMC: \n"))
  print(time_ABC)
  
  cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
  cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)
  

  if (i>1){
    draws$draws <- checkFitWithPreviousExperiments(draws$draws, objectiveFunction, delta)
  }
  cat("-Saving sample \n")
  outFile <- paste(seq(1,i), collapse="_")
  timeStr <- Sys.time()
  timeStr <- gsub(":","_", timeStr)
  timeStr <- gsub(" ","_", timeStr)
  outFileR <- paste0("./PosteriorSamples/DrawsStochastic",modelName,"_",basename(comment(modelName)),"_ns",ns,"_npc",npc,"_",outFile,timeStr,".RData",collapse="_")
  save(draws, parNames, file=outFileR)
#}
end_time = Sys.time()
time_ = end_time - start_time
cat("Total time:")
print(time_)

#### PLOT RESTULTS FOR AKAR4
# par(mfrow=c(2,3))
# for(i in 1:3){
#   hist(ABCMCMCoutput$draws[,i], main=parNames[i], xlab = "Value in log scale")
# }
# combinePar <- list(c(1,2), c(1,3), c(2,3))
# for(i in combinePar){
#   plot(ABCMCMCoutput$draws[,i[1]], ABCMCMCoutput$draws[,i[2]], xlab = parNames[i[1]], ylab = parNames[i[2]])
# }


#### PLOT SIMULATIONS FROM DRAWS
simulateSSA <- function(e, param, nStochSim){
  avgOutput <- rep(0, length(e[["outputTimes"]]))
  for(i in 1:nStochSim){
    out_ssa <- GillespieSSA2::ssa(
      initial_state = ceil(e[["initialState"]]*Phi),
      reactions = compiled_reactions,
      params = c(parMap(param), Phi=Phi),
      final_time = max(e[["outputTimes"]]),
      method = ssa_exact(),
      verbose = FALSE,
      log_propensity = TRUE,
      log_firings = TRUE,
      census_interval = 0.001,
      sim_name = modelName)

    # out$state is a matrix of dimension (time points)x(num compounds)
    output <- apply(out_ssa$state/Phi, 1, function(state) model$func(t=0,state=state,parameters=param))
    interpOutput <- approx(out_ssa$time, output, e[["outputTimes"]])
    #interpOutput$y[is.na(interpOutput$y)] <- tail(output,1)
    avgOutput <- avgOutput + interpOutput$y
  }
  avgOutput <- avgOutput/nStochSim
  return(avgOutput)
}


exp.ind <- 2
par(mfrow=c(1,1))
plot(experiments[[exp.ind]][["outputTimes"]],experiments[[exp.ind]][["outputValues"]][[1]],ylim=c(90,250))
for(i in 1:100){
  param <- ABCMCMCoutput$draws[i,]
  #param <- pC$prePar[,i]
  names(param) <- parNames
  sim <- simulateSSA(experiments[[exp.ind]], param, nStochSim = 3)
  lines(experiments[[exp.ind]][["outputTimes"]],sim, col="blue")
}
points(experiments[[exp.ind]][["outputTimes"]],experiments[[exp.ind]][["outputValues"]][[1]],ylim=c(90,250))



library(plotly)
df = as.data.frame(ABCMCMCoutput$draws)
colnames(df) <- parNames
plot_ly(dat = df, x = ~kf_C_AKAR4, y = ~kb_C_AKAR4, z = ~kcat_AKARp, type="scatter3d", mode="markers", marker=list(size = 1, color = "red"))

