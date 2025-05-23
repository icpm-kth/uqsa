#!/usr/bin/env Rscript

# pt = parallel tempering
# mh = Metropolis Hastings algorithm
# vF = vesicularFusogenicity

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(pbdMPI)

start_time <- Sys.time()                         # measure sampling time
comm  <- 0
r <- pbdMPI::comm.rank(comm=comm)
cs <- pbdMPI::comm.size(comm=comm)
pbdMPI::init()
attr(comm,"rank") <- r
attr(comm,"size") <- cs

# random number seed:
M <- 7 # integration method
set.seed(137*r + 1337*cs)

N <- 10000
h <- 2e-2
cycl <- 7
PREFIX <- "main"
## possibly adjust sampling-settings
a <- commandArgs(trailingOnly=TRUE)
if (!is.null(a) && length(a)>0 && length(a) %% 2 == 0){
	for (i in seq(1,length(a),by=2)){
		key = a[i]
		val = a[i+1]
		switch(key,
			"-N"={N <- as.integer(val)},
			"-h"={h <- as.double(val)},
			"-c"={cycl <- as.integer(val)},
			"--prefix"={PREFIX=as.character(val)},
			{warning(sprintf("unknown option %s %s",key,val))}
		)
	}
}
R <- round(cs/2)
beta <- (1.0 - (r<R)*(r/R))^2
if (N %% 4 != 0) N <- round(N/4) * 4

message(sprintf("rank %i of %i will sample %i points.\n",r,cs,N))
## ----label="SBtab content"----------------------------------------------------
#modelFiles <- uqsa_example("CaMKII")
modelFiles <- dir(pattern="tsv$",full.names=TRUE)
sb <- SBtabVFGEN::sbtab_from_tsv(modelFiles)
modelName <- checkModel(comment(sb),sprintf("./%s.so",comment(sb)))
source("./CaMKIIs.R")
experiments <- sbtab.data(sb)[seq(5)]
options(mc.cores = length(experiments))
## ----default------------------------------------------------------------------
n <- length(experiments[[1]]$input)

stopifnot(all(trimws(sb$Parameter[["!Scale"]])=="natural logarithm"))
parMCMC <- sb$Parameter[["!DefaultValue"]]/log(10)  # this was in natural-log-space


if ("!Std" %in% names(sb$Parameter)) {
	stdv <- sb$Parameter[["!Std"]]/log(10)
} else {
	warning("no standard error («!Std» field) in 'SBtab$Parameter'")
	stdv <- parMCMC*0.5 + 0.5 + 0.5*max(parMCMC)
}
dprior <- dNormalPrior(mean=parMCMC,sd=stdv)
rprior <- rNormalPrior(mean=parMCMC,sd=stdv)
gprior <- gradLog_NormalPrior(mean=parMCMC,sd=stdv)

## ----simulate-----------------------------------------------------------------
simulate <- simc(experiments,modelName,log10ParMap,approximateSensitivity=TRUE,method=M)

y <- simulate(parMCMC) ## little test
stopifnot(all(c("state","func") %in% names(y[[1]])))

llf <- logLikelihoodFunc(experiments)
gradLL <- gradLogLikelihoodFunc(model, experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
fiIn <- fisherInformationFunc(model, experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
fiPrior <- solve(diag(stdv, length(parMCMC)))

X <- NULL

sampleSize <- round(exp(seq(log(100),log(N),length.out=cycl)))

## ----sample-------------------------------------------------------------------
	smmala <- mcmcUpdate(simulate=simulate,
		experiments=experiments,
		model=model,
		logLikelihood=llf,
		dprior=dprior,
		gradLogLikelihood=gradLL,
		gprior=gprior,
		fisherInformation=fiIn,
		fisherInformationPrior=fiPrior
	)
	mhmcmc <- mcmc(smmala) # no-swap version for step-size tuning
	ptSmmala <- mcmc_mpi(
		smmala,
		comm=comm,
		swapDelay=0,
		swapFunc=pbdMPI_bcast_reduce_temperatures
	)
	for (i in seq(cycl)){
		rm(X)
		sampleFile=sprintf("%s-Sample-cycle-%i-rank-%02i-of-%02i.RDS",PREFIX,i,r,cs)
		x <- mcmcInit(beta,parMCMC,simulate,llf,dprior,gradLL,gprior,fiIn)
		if (i < round(cycl/2)){
			for (k in seq(5)){
				s <- mhmcmc(x,300,h)
				ar <- attr(s,"acceptanceRate")
				message(sprintf("[r%ic%ik%i] h=%g, a=%g",as.integer(r),as.integer(i),k,h,ar))
				h <- h * (0.5 + ar^4/(0.25^4 + ar^4)) # 0.25 is the target acceptance value
			}
		}
		s <- ptSmmala(x,sampleSize[i],h) # the main amount of work is done here
		cat(sprintf("rank %02i/%02i finished %i iterations for cycle %i with acceptance rate of %02i %% and swap rate of %02i %% (step-size is %g).\n",r,cs,sampleSize[i],i,round(100*attr(s,"acceptanceRate")),round(100*attr(s,"swapRate")),h))
		saveRDS(s,file=sampleFile)
		rm(s)
		gc()
		pbdMPI::barrier()
		f <- dir(pattern=sprintf('^%s-Sample-.*RDS$',PREFIX))
		X <- uqsa::gatherSample(f,beta,size=min(sampleSize[i],40000))
		attr(X,"beta") <- beta
		parMCMC <- as.numeric(tail(X,1))
	}

saveRDS(X,file="%s-final-sample-%i-%i-lb100-%i.RDS",PREFIX,r,cs,round(log(beta)*100))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
finalize()
print(warnings())
