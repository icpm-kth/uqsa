#!/usr/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(pbdMPI)

pbdMPI::init()

comm=0
r <- pbdMPI::comm.rank(comm=comm)
cs <- pbdMPI::comm.size(comm=comm)
attr(comm,"rank") <- r
attr(comm,"size") <- cs

start_time <- Sys.time()                         # measure sampling time
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
nChains <- Args['nChains']
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
print(ConLaw$Text)
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
simulate <- simcf(experiments,modelName,log10ParMap)

#options(mc.cores = 2)
#simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox=sensApprox)
y <- simulate(parVal)
stopifnot(is.list(y) && length(y)==length(experiments) && !all(c("state","func") %in% names(y[[1]])))

## ----likelihood---------------------------------------------------------------
llf <- logLikelihoodFunc(experiments)

## ----update-------------------------------------------------------------------
## metropolis <- metropolisUpdate(simulate, experiments, model, logLikelihood=llf, dprior=dprior)

## (simulate, experiments, model, logLikelihood, dprior, gradLogLikelihood=NULL, gprior=NULL, fisherInformation=NULL, fisherInformationPrior=NULL)
metropolis <- mcmcUpdate(simulate=simulate,
		          experiments=experiments,
		          model=model,
		          logLikelihood=llf,
		          dprior=dprior)
## ----init---------------------------------------------------------------------
#m <- mcmc(smmala)                     # a serial Markov chain Monte Carlo function
ptMetropolis <- mcmc_mpi(metropolis,comm=comm,swapDelay=0,swapFunc=pbdMPI_bcast_reduce_temperature)# MPI aware function, passes messages on "comm"
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
	x <- mcmcInit(beta,x,simulate,llf,dprior)
	for (j in seq(nj)) {
		Sample <- ptMetropolis(x,100,eps=h)
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
	x <- mcmcInit(beta,x,simulate,llf,dprior)
}
## ----sample-------------------------------------------------------------------

s <- ptMetropolis(x,Args['N'],h) # the main amount of work is done here
colnames(s) <- names(parVal)
x <- attr(s,"lastPoint")
beta <- attr(x,"beta")
saveRDS(s,file=sprintf("Rmpi-testSample-rank%i-of%i.RData",r,cs))
save(x,h,beta,file=initFile)
cat(sprintf("rank %02i/%02i finished with acceptance rate of %02i %% and swap rate of %02i %%.\n",r,cs,round(100*attr(s,"acceptanceRate")),round(100*attr(s,"swapRate"))))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
finalize()

