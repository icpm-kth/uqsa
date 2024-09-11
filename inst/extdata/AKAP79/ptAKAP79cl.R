#!/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(pbdMPI)

MPI <- 'pbdMPI'

start_time <- Sys.time()                         # measure sampling time
comm  <- 0
if (MPI == 'Rmpi'){
	r <- Rmpi::mpi.comm.rank(comm=comm)
	cs <- Rmpi::mpi.comm.size(comm=comm)
} else if (MPI == "pbdMPI"){
	pbdMPI::init()
	r <- pbdMPI::comm.rank(comm=comm)
	cs <- pbdMPI::comm.size(comm=comm)
} else {
	stop('unknown type of MPI')
}
attr(comm,"rank") <- r
attr(comm,"size") <- cs

Args=c(N=1000,h=NA,nChains=cs) # nChains can be used to pretend like we have more than cs chains, for the calculation of beta

a <- commandArgs(trailing=TRUE)

if (!is.null(a) && length(a)>0) {
	if (length(a)==1) {
		Args['N']=as.numeric(a[1])
	} else {
		a <- strsplit(a,"=")
		Args <- Reduce(\(a,b) {c(a,as.numeric(b[2]))},a,init=NULL)
		names(Args) <- Reduce(\(a,b) c(a,make.names(b[1])),a,init=NULL)
	}
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
if (0 == r) {
	cat("Conservation laws for this model:")
	print(ConLaw$Text)
}
experiments <- sbtab.data(SBtab,ConLaw)


## ----default------------------------------------------------------------------
n <- length(experiments[[1]]$input)
stopifnot(n>0)
parVal <- c( 1.84512,  -1.11594,  -3.28677,  -0.225118,  0.630559,  -1.44445,  -3.19814,  -4.59471,  -2.1549,  -4.4207,  -2.17395,  0.937206,  -1.53715,  0.163359,  -2.40479,  2.12073,  -3.34629,  -0.913875,  0.492607,  -0.532405,  1.75744,  -1.6767,  -1.02841,  0.982068,  1.93345,  0.118566,  -0.197953)

# alternatively:
if (n>0){
	defaultVal <- log10(head(AKAP79_default(),-n))
} else {
	defaultVal <- log10(AKAP79_default())
}
## ----range--------------------------------------------------------------------
defRange <- 2 # log-10 space
dprior <- dNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
rprior <- rNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
gprior <- gradLog_NormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
## ----simulate-----------------------------------------------------------------
simulate <- simc(experiments,modelName,log10ParMap)

#options(mc.cores = 2)
#simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox=sensApprox)
y <- simulate(parVal)
#print(names(y[[1]]))

## ----likelihood---------------------------------------------------------------
llf <- logLikelihoodFunc(experiments)
gradLL <- gradLogLikelihoodFunc(model,experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
fiIn <- fisherInformationFunc(model, experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
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
ptsmmala <- mcmc_mpi(smmala,comm=comm,swapDelay=0,swapFunc=pbdMPI_swap_temperatures)
## ptsmmala is a now an MPI aware function, passes messages on "comm"
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
	txtLog <- paste0("adjusting-step-size-for-rank-",r,"-of-",cs,"-",gsub(" ","T",Sys.time()),".txt")
	for (j in seq(nj)) {
		Sample <- ptsmmala(x,100,eps=h)
		a <- attr(Sample,"acceptanceRate")
		h <- h * L(a)
		x <- attr(Sample,"lastPoint")
		beta <- attr(x,"beta")
		cat(sprintf("iteration %02i/%02i for rank\t%02i/%02i,\th = %10g,\tacceptance = %i %%\n",j,nj,r,cs,h,round(100*a)),file=txtLog,append=TRUE)
		flush.console()
	}
	save(x,h,beta,file=initFile)
	cat("final step size: ",h,"\n",file=txtLog,append=TRUE)
	cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n",file=txtLog,append=TRUE)
} else {
	x <- mcmcInit(beta,x,simulate,llf,dprior,gradLL,gprior,fiIn)
}

## ----sample-------------------------------------------------------------------
sfile <- sprintf("Rmpi-AKAP79-rank%i-of%i.RData",r,cs)
if (file.exists(sfile)){ # continue from last sample
	s <- readRDS(sfile)
	x <- attr(s,"lastPoint")
}
s <- ptsmmala(x,Args['N'],h) # the main amount of work is done here
colnames(s) <- names(parVal)
saveRDS(s,file=sfile)
x <- attr(s,"lastPoint")
beta <- attr(x,"beta")

cat(sprintf("rank %02i/%02i finished with acceptance rate of %02i %% and swap rate of %02i %%.\n",r,cs,round(100*attr(s,"acceptanceRate")),round(100*attr(s,"swapRate"))))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

if (MPI == "Rmpi"){
	Rmpi::mpi.finalize()
} else {
	finalize()
}
