#!/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(Rmpi)

start_time <- Sys.time()                         # measure sampling time
r <- mpi.comm.rank(comm=0)
cs <- mpi.comm.size(comm=0)
nChains <- cs # this number can be used to pretend like we have more than cs chains, for the calculation of beta

a <- commandArgs(trailing=TRUE)

if (!is.null(a) && length(a)>0) {
	if (length(a)==1) {
		Args <- c(N=as.numeric(a[1]),h=NA,nChains=cs)
	} else {
		a <- strsplit(a,"=")
		Args <- Reduce(\(a,b) {c(a,as.numeric(b[2]))},a,init=NULL)
		names(Args) <- Reduce(\(a,b) c(a,make.names(b[1])),a,init=NULL)
	}
} else {
	Args <- c(N=300,h=NA,nChains=cs)
}
print(Args)
N <- Args['N']
beta <- (1.0 - (r/nChains))^2

cat(sprintf("rank %i of %i workers will sample %i points.\n",r,cs,N))
## ----label="SBtab content"----------------------------------------------------
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
gprior <- gradLog_NormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
## ----simulate-----------------------------------------------------------------
sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
simulate <- simc(experiments,modelName,log10ParMap,sensApprox)

#options(mc.cores = 2)
#simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox=sensApprox)
y <- simulate(parVal)
#print(length(y))

## ----likelihood---------------------------------------------------------------
llf <- logLikelihoodFunc(experiments)
gradLL <- gradLogLikelihoodFunc(model,experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
fiIn <- fisherInformationFunc(model, experiments, parMap=log10ParMap)
fiPrior <- solve(diag(defRange, length(parVal)))

## ----update-------------------------------------------------------------------
smmala <- mcmcUpdate(simulate=simulate,
		          experiments=experiments,
		          model=model,
		          logLikelihood=llf,
		          dprior=dprior,
		          gradLogLikelihood=gradLL,
		          gprior=gprior,
		          fisherInformation=fiIn,
		          fisherInformationPrior=fiPrior)
## ----init---------------------------------------------------------------------
#m <- mcmc(smmala)                     # a serial Markov chain Monte Carlo function
ptsmmala <- mcmc_mpi(smmala,comm=0)    # MPI aware function, passes messages on "comm"
#smmala <- mcmc(smmala)
## ----adjust-------------------------------------------------------------------
accTarget <- 0.25
L <- function(a) { (1.0 / (1.0+exp(-(a-accTarget)/0.1))) + 0.5 }

start_time <- Sys.time()
x <- parVal
nj <- 20
h <- Args['h']
initFile <- sprintf("rmpi-init-rank-%i-of-%i.RData",r,cs)
if (file.exists(initFile)){
	load(initFile)
} else if (is.na(h)) {
	h <- 1e-3
	x <- mcmcInit(beta,x,simulate,llf,dprior,gradLL,gprior,fiIn)
	for (j in seq(nj)) {
		Sample <- ptsmmala(x,100,eps=h)
		a <- attr(Sample,"acceptanceRate")
		h <- h * L(a)
		x <- attr(Sample,"lastPoint")
		beta <- attr(x,"beta")
		cat(sprintf("iteration %02i/%02i for rank %02i/%02i,\th = %g,\tacceptance = %i %%\n",j,nj,r,cs,h,round(100*a)))
		flush.console()
	}
	save(x,h,beta,file=initFile)
	cat("final step size: ",h,"\n")
	cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n")
} else {
	x <- mcmcInit(beta,x,simulate,llf,dprior,gradLL,gprior,fiIn)
}

## ----sample-------------------------------------------------------------------
s <- ptsmmala(x,Args['N'],h) # the main amount of work is done here
colnames(s) <- names(parVal)
saveRDS(s,file=sprintf("Rmpi-testSample-rank%i-of%i.RData",r,cs))
x <- attr(s,"lastPoint")
beta <- attr(x,"beta")

cat(sprintf("rank %02i/%02i finished with acceptance rate of %02i %% and swap rate of %02i %%.\n",r,cs,round(100*attr(s,"acceptanceRate")),round(100*attr(s,"swapRate"))))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
mpi.finalize()

