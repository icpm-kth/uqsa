#!/usr/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(SBtabVFGEN)
library(pbdMPI)

start_time <- Sys.time()                         # measure sampling time
comm  <- 0
if (require(pbdMPI)){
	MPI <- "pbdMPI"
	r <- pbdMPI::comm.rank(comm=comm)
	cs <- pbdMPI::comm.size(comm=comm)
	pbdMPI::init()
} else if (require(Rmpi)){
	MPI <- "Rmpi"
	r <- Rmpi::mpi.comm.rank(comm=comm)
	cs <- Rmpi::mpi.comm.size(comm=comm)
} else {
	stop('unknown type of MPI')
}
attr(comm,"rank") <- r
attr(comm,"size") <- cs

Args <- c(N=100,h=NA,nChains=cs)
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
beta <- (1.0 - (r/nChains))

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
## ----simulate-----------------------------------------------------------------
simulate <- simcf(experiments,modelName,log10ParMap) # no sensitivities

#options(mc.cores = 2)
#simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox=sensApprox)
y <- simulate(parVal)
stopifnot(is.list(y) && length(y)==length(experiments))
stopifnot(all(c("state","func") %in% names(y[[1]])))
#print(length(y))

## ----likelihood---------------------------------------------------------------
llf <- logLikelihoodFunc(experiments)

## ----update-------------------------------------------------------------------
metropolis <- mcmcUpdate(simulate=simulate,
		          experiments=experiments,
		          model=model,
		          logLikelihood=llf,
		          dprior=dprior)
## ----init---------------------------------------------------------------------
#m <- mcmc(smmala)                     # a serial Markov chain Monte Carlo function
ptMetropolis <- mcmc_mpi(metropolis,comm=comm,swapDelay=0,swapFunc=pbdMPI_bcast_reduce_temperatures)# MPI aware function, passes messages on "comm"
serialMetropolis <- mcmc(metropolis)
#smmala <- mcmc(smmala)
## ----adjust-------------------------------------------------------------------
accTarget <- 0.25
L <- function(a) { (1.0 / (1.0+exp(-(a-accTarget)/0.1))) + 0.5 }

start_time <- Sys.time()
x <- parVal
h <- Args['h']
initFile <- sprintf("%s-init-rank-%i-of-%i.RData",MPI,r,cs)
if (file.exists(initFile)){
	message(sprintf("loading initial value and step size from file: %s",initFile))
	load(initFile, verbose = TRUE)
} else if (is.na(h)) {
	h <- 1.0
	x <- mcmcInit(beta,x,simulate,llf,dprior)
	txtLog <- paste0("adjusting-step-size-for-rank-",r,"-of-",cs,"-",gsub(" ","T",Sys.time()),".txt")
	nj <- 60
	for (j in seq(nj)) {
		stopifnot(is.numeric(x) && all(is.finite(x)))
		if (j%%2==1) {
			Sample <- ptMetropolis(x,100,eps=h) # make some progress, don't change h
			Label <- MPI
			swapRate <- attr(Sample,"swapRate")
		} else {
			Sample <- serialMetropolis(x,100,eps=h) # evaluate acceptance rate, make changes to h
			Label <- "no MPI"
			swapRate <- NA
		}
		a <- attr(Sample,"acceptanceRate")
		x <- attr(Sample,"lastPoint")
		beta <- attr(x,"beta")
		h <- attr(Sample,"stepSize")

		if (j%%2==0) h <- h * L(a)              # here is where the change happens
		## logging and info print
		Msg <- sprintf("[%s] iteration %02i/%02i for rank\t%02i/%02i (β=%0.3f),\tlog10(h) = %5g,\tacceptance = %i %%, swap rate = %i %%\n",Label,j,nj,r,cs,beta,log10(h),round(100*a),round(swapRate*100))
		cat(Msg)
		cat(Msg,file=txtLog,append=TRUE)
	}
	save(x,h,beta,file=initFile)
	cat("final step size: ",h,"\n",file=txtLog,append=TRUE)
	cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n",file=txtLog,append=TRUE)
} else {
	x <- mcmcInit(beta,x,simulate,llf,dprior)
}
## ----sample-------------------------------------------------------------------

s <- ptMetropolis(x,Args['N'],h) # the main amount of work is done here
colnames(s) <- names(parVal)
saveRDS(s,file=sprintf("%s-testSample-rank%i-of%i.RData",MPI,r,cs))
x <- attr(s,"lastPoint")
beta <- attr(x,"beta")
h <- attr(x,"stepSize")

save(x,h,beta,file=initFile)
cat(sprintf("rank %02i/%02i finished with acceptance rate of %02i %% and swap rate of %02i %%.\n",r,cs,round(100*attr(s,"acceptanceRate")),round(100*attr(s,"swapRate"))))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)
##print(warnings())

if (MPI == "Rmpi"){
	Rmpi::mpi.finalize()
} else {
	finalize()
}
