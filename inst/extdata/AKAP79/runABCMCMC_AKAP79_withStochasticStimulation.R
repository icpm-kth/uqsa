# ABCMCMC with stochastic simulation
library(GillespieSSA2)
library(uqsa)
library(rgsl)
library(SBtabVFGEN)
library(parallel)
library(pracma)

nChains <- 4

options(mc.cores=parallel::detectCores())
model.tsv <- uqsa_example("AKAP79")
model.tab <- sbtab_from_tsv(model.tsv)

# source all R functions for this model, this also loads a variable called "model"
source(uqsa_example("AKAP79",pat="^AKAP79[.]R$"))

# without conservation laws
experiments <- sbtab.data(model.tab)

# with conservation laws
#load("./inst/extdata/AKAP79/ConservationLaws.RData")
#experiments <- sbtab.data(model.tab,ConLaw)

print(comment(model.tab))
modelName <- checkModel(comment(model.tab),uqsa_example("AKAP79",pat="_gvf.c"))

numPar <- nrow(model.tab$Parameter)
parVal <- model$par()[1:numPar]
parNames <- row.names(model.tab$Parameter)


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

set.seed(7619201)


distanceMeasure <- function(funcSim, dataExpr=Inf, dataErr=Inf){
  if (all(is.finite(funcSim))){
    distance <- mean((1/71.67*(funcSim-as.matrix(dataExpr))/as.matrix(dataErr))^2, na.rm=TRUE)
  } else {
    distance <- Inf
  }
  return(distance)
}


save_sample <- function(ABCMCMCoutput){
  timeStr <- gsub("[ :]","_", Sys.time())
  if (!dir.exists("./PosteriorSamples")) {
    dir.create("./PosteriorSamples")
  }
  outFileR <- paste("./PosteriorSamples/Draws_StochSim_",modelName,"nChains",nChains,"ns",ns,"npc",npc,timeStr,".RData",collapse="_",sep="_")
  save(ABCMCMCoutput, file=outFileR)
}


# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
#ns <- 1000 or more on a cluster
ns <- 100 # no of samples required from each ABC-MCMC chain
#npc <- 1000 or more on a cluster
npc <- 100 # pre-calibration

# Threshold for distance between experimental data and simulations
delta <- 0.01

# Define ABC-MCMC Settings
p <- 0.01		 # For the Pre-Calibration: Choose Top 1% Samples with Shortest Distance to the Experimental Values

# Define volume where the reactions take place (e.g., volume of neuron synapse) and unit of measure (1e-6 = micrometer)
vol <- 4e-16
unit <- 1e-06

# Create reactions and parameters for stochastic simulations with the GillespieSSA2 package
reactions <- importReactionsSSA(model.tab)

# Loop through the Different Experimental Settings
start_time = Sys.time()

expInd <- c(16,17,18)
cat("#####Starting run for Experiments ", expInd, "######\n")
objectiveFunction <- makeObjectiveSSA(experiments[expInd],
                                      model.tab = model.tab, 
                                      parNames = parNames,
                                      distance = distanceMeasure,
                                      parMap = parMap, 
                                      outputFunction = model$func,
                                      vol = vol,
                                      unit = unit,
                                      reactions = reactions, 
                                      nStochSim = 3)

## If First Experimental Setting, Create an Independent Colupla
rprior <- rNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5) 
dprior <- dNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)

## Run Pre-Calibration Sampling
cat(sprintf("- Precalibration \n"))
  
time_pC <- Sys.time()
sfactor <- 0.1 # scaling factor
pC <- preCalibration(objectiveFunction, npc, rprior, rep = 1, sfactor = sfactor, num = nChains)
time_pC <- Sys.time() - time_pC
cat(sprintf("- time for precalibration: \n"))
print(time_pC)
  
  
## Get Starting Parameters from Pre-Calibration
Sigma <- pC$Sigma
startPar <- pC$startPar
for(j in 1 : nChains){
    cat("Chain", j, "\n")
    obFun <- objectiveFunction(startPar[,j])
    cat("\tMin distance of starting parameter for chain",j," = ", min(obFun),"\n")
    cat("\tMean distance of starting parameter for chain",j," = ", mean(obFun),"\n")
    cat("\tMax distance of starting parameter for chain",j," = ", max(obFun),"\n")
}
  
## Run ABC-MCMC Sampling
cat(sprintf("- Running MCMC\n"))
time_ABC <- Sys.time()
cl <- makeForkCluster(nChains, outfile="outputMessagesABCMCMC.txt")
out_ABCMCMC <- parLapply(cl,
                           1:nChains,
                           function(j) {
                             tryCatch(
                               ABCMCMC(objectiveFunction, startPar[,j], ns, Sigma, delta, dprior, batchSize = 1),
                               error=function(cond) {message("ABCMCMC crashed"); return(NULL)})
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
  
time_ABC <- Sys.time() - time_ABC
cat(sprintf("- time for ABCMCMC: \n"))
print(time_ABC)
  
cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)
  
#save_sample(ABCMCMCoutput)

end_time = Sys.time()
time_ = end_time - start_time
cat("Total time:")
print(time_)


##### PLOT RESULTS #####

# HISTOGRAM (posterior marginal densities)
for(i in 1:length(parVal)){
  hist(ABCMCMCoutput$draws[,i], main=par_names[i], xlab = "Value in log scale")
}


# PLOT DRAWS of first state
plot(ABCMCMCoutput$draws[,1])

# PLOT EXPERIMENT 
expInd <- 16
e <- experiments[[expInd]]

plot(e$outputTimes, e$outputValues$AKAR4pOUT)

n_draws <- dim(ABCMCMCoutput$draws)[1]
n_sample_param <- 5
idx_sample_param <- sample(1:n_draws, n_sample_param)
idx_min_delta <- head(sort(ABCMCMCoutput$scores, index.return = TRUE)$ix)
for(i in idx_min_delta){
  param <- ABCMCMCoutput$draws[i, ]
  names(param) <- parNames
  sim <- simulator.stoch(experiments,
                         model.tab = model.tab,
                         reactions = reactions, 
                         parMap = parMap, 
                         outputFunction = model$func, 
                         vol = vol, 
                         unit = unit, 
                         nStochSim = 3)
  # compute trajectory in all experimental conditions
  output <- sim(param)
  
  # trajectory with experimental conditions as in the (expInd)-th experiment
  output_expInd <- output[[expInd]]
  
  lines(output_expInd$time, output_expInd$output, col="blue")
}

