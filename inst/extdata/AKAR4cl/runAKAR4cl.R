## this project:
library(uqsa)
## our other packages:
require(rgsl)
require(SBtabVFGEN)
require(parallel)

modelFiles <- uqsa_example("AKAR4cl",pattern="[.]tsv$",full.names=TRUE)
SBtab <- SBtabVFGEN::sbtab_from_tsv(modelFiles)

# this will load a series of functions: AKAR4cl_vf, AKAR4cl_jac, etc.
# also loads "model", a list of the same functions with generic names:
source(uqsa_example("AKAR4cl",pat="^AKAR4cl[.]R$"))
names(model)

modelName <- checkModel("AKAR4cl",uqsa_example("AKAR4cl",pat="_gvf[.]c$"))

# load experiments
load(uqsa_example("AKAR4cl",pat="^ConservationLaws[.]RData$")) # loads "ConLaw" variable
experiments <- sbtab.data(SBtab,ConLaw)

#default values for the parameters:
n <- length(experiments[[1]]$input)
stopifnot(n>0)
parVal <- head(AKAR4cl_default(),-n)
# scale to determine prior values
defRange <- 100

# Define Lower and Upper Limits for logUniform prior distribution for the parameters

start_time = Sys.time()
sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox)
# test run
y <- rgsl::r_gsl_odeiv2_outer(modelName,experiments,matrix(parVal))

# more test run
y <- simulate(log10(parVal))

dprior <- dNormalPrior(mean=log10(parVal),sd=rep(2,length(parVal)))
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

m <- mcmc(update)
h <- 1e-1

# adjusting  step size h:
L <- function(a) { (1.0 / (1.0+exp(-(a-0.25)/0.1))) + 0.5}
parMCMC <- log10(parVal)

N <- 100

for (j in seq(5)){
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

options(mc.cores = parallel::detectCores() %/% 4)
cl <- parallel::makeForkCluster(4)

pL <- parLapply(cl, rep(list(parMCMC),4), m, N=4000, eps=h)
sample <- Reduce(rbind,pL)
colnames(sample) <- names(parVal)

time_ = Sys.time() - start_time
stopCluster(cl)
cat("finished sampling after",time_," seconds\n")
hexbin::hexplom(sample)
