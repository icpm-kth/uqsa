#!/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(hexbin)
library(SBtabVFGEN)
library(Rmpi)

start_time <- Sys.time()                         # measure sampling time
r <- mpi.comm.rank(comm=0)
cs <- mpi.comm.size(comm=0)
beta <- (1.0 - (r/cs))^4
nChains <- cs

a <- commandArgs(trailing=TRUE)
N <- 300 # default
if (!is.null(a)) {
	N <- as.numeric(a[1])
}

## ## ----label="SBtab content"----------------------------------------------------
modelFiles <- uqsa_example("AKAP79",pattern="[.]tsv$",full.names=TRUE)
SBtab <- SBtabVFGEN::sbtab_from_tsv(modelFiles)

## ----model--------------------------------------------------------------------
source(uqsa_example("AKAP79",pat="^AKAP79[.]R$"))
modelName <- checkModel("AKAP79","./AKAP79.so")

## ----ConservationLaws---------------------------------------------------------
load(uqsa_example("AKAP79",pat="^ConservationLaws[.]RData$"))
##print(ConLaw$Text)
experiments <- sbtab.data(SBtab,ConLaw)


## ----default------------------------------------------------------------------
n <- length(experiments[[1]]$input)
stopifnot(n>0)
parVal <- log10(head(AKAP79_default(),-n))

## ----range--------------------------------------------------------------------
defRange <- 2 # log-10 space
dprior <- dNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
rprior <- rNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))

## ----simulate-----------------------------------------------------------------
sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
simulate <- simc(experiments,modelName,log10ParMap,sensApprox)
y <- simulate(parVal)

## ----likelihood---------------------------------------------------------------
llf <- logLikelihood(experiments)
gradLL <- gradLogLikelihood(model,experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
fiIn <- fisherInformation(model, experiments, parMap=log10ParMap)
fiPrior <- solve(diag(defRange, length(parVal)))


## ----update-------------------------------------------------------------------
update <- mcmcUpdate(simulate=simulate,
		          experiments=experiments,
		          model=model,
		          logLikelihood=llf,
		          gradLogLikelihood=gradLL,
		          fisherInformation=fiIn,
		          fisherInformationPrior=fiPrior,
		          dprior=dprior)
## ----init---------------------------------------------------------------------
m <- mcmc(update)   # a Markov chain function
h <- 1e-3           # initial step size guess
## ----adjust-------------------------------------------------------------------
accTarget <- 0.25
L <- function(a) { (1.0 / (1.0+exp(-(a-accTarget)/0.1))) + 0.5 }

start_time <- Sys.time()
x <- parVal
nj <- 2
initFile <- sprintf("rmpi-init-rank-%i-of-%i.RData",r,cs)
if (file.exists(initFile)){
	load(initFile)
} else {
	for (j in seq(nj)){
		x <- mcmcInit(beta,x,simulate,dprior,llf,gradLL,fiIn)
		Sample <- m(x,100,eps=h)
		a <- attr(Sample,"acceptanceRate")
		h <- h * L(a)
		x <- attr(Sample,"lastPoint")
		cat(sprintf("iteration %02i/%02i for rank %02i/%02i,\th = %g\n",j,nj,r,cs,h))
	}
	save(x,h,beta,file=initFile)
	cat("final step size: ",h,"\n")
	cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n")
}
## ----sample-------------------------------------------------------------------
m <- mcmc_mpi(update,comm=0)    # m is now an mpi aware function, and passes messages on "comm"
s <- m(x,N,h) # the main amount of work is done here
colnames(s) <- names(parVal)

saveRDS(s,file=sprintf("Rmpi-testSample-rank%i-of%i.RData",r,cs))
time_ <- difftime(Sys.time(),start_time,units="sec")
mpi.finalize()
print(time_)
