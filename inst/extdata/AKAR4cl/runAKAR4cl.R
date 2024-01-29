## this project:
devtools::load_all("~/uqsa")
## our other packages:
require(rgsl)
require(SBtabVFGEN)

model.tsv <- dir("AKAR4cl",pattern="[.]tsv$",full.names=TRUE)
model.tab <- sbtab_from_tsv(model.tsv) # SBtabVFGEN
source("./AKAR4cl.R")

modelName <- checkModel(comment(model.tab),"./AKAR4cl.so")

# load experiments
load("ConservationLaws.RData",verbose=TRUE) #ConLaw
experiments <- sbtab.data(model.tab,ConLaw)

n <- length(experiments[[1]]$input)
if (n>0) {
	parVal <- head(AKAR4cl_default(),-n)
} else {
	parVal <- AKAR4cl_default()
}

# scale to determine prior values
defRange <- 100

# Define Lower and Upper Limits for logUniform prior distribution for the parameters

start_time = Sys.time()
                                                # experiments, model, parMap=identity, parMapJac=1.0
sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox)
# test run
y <- rgsl::r_gsl_odeiv2_outer(modelName,experiments,matrix(parVal))

# more test run
y <- simulate(log10(parVal))

dprior <- dNormalPrior(mean=log10(parVal),sd=rep(2,length(parVal)))
rprior <- rNormalPrior(mean=log10(parVal),sd=rep(2,length(parVal)))
llf <- logLikelihood(experiments)
gradLL <- gradLogLikelihood(model,experiments,parMap=log10ParMap, parMapJac=log10ParMapJac)
fiIn <- fisherInformation(model, experiments,parMap=log10ParMap)
fiPrior <- solve(diag(2,length(parVal)))
update  <- mcmcUpdate(simulate=simulate,
	experiments=experiments,
	model=model,
	logLikelihood=llf,
	gradLogLikelihood=gradLL,
	fisherInformation=fiIn,
	fisherInformationPrior=fiPrior,
	dprior)

#parMCMC <- mcmcInit(log10(parVal),simulate,llf,gradLL,fiIn)
#yinint <- attr(parMCMC,"simulations")
m <- mcmc(update)
h <- 1e-1

# adjusting  step size h:
L <- function(a) { (1.0 / (1.0+exp(-(a-0.25)/0.1))) + 0.5}
parMCMC <- log10(parVal)

N <- 100

for (j in seq(10)){
 cat("adjusting step size: ",h," \n");
 parMCMC <- mcmcInit(parMCMC,simulate,llf,gradLL,fiIn)
 sample <- m(parMCMC,N,eps=h)
 a <- attr(sample,"acceptanceRate")
 cat("acceptance: ",a*100," %\n")
 h <- h * L(a)
 parMCMC <- sample[N,]
}
cat("final step size: ",h,"\n")
parMCMC <- mcmcInit(parMCMC,simulate,llf,gradLL,fiIn)
cat("finished adjusting after",Sys.time() - start_time," seconds\n")

sample <- m(parMCMC,N=10000,eps=h)
colnames(sample) <- names(parVal)

time_ = Sys.time() - start_time
cat("finished sampling after",time_," seconds\n")
hexbin::hexplom(sample)
