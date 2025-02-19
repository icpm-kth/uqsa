#!/usr/bin/env Rscript
library(rgsl)
library(uqsa)
library(pbdMPI)

start_time <- Sys.time()
comm  <- 0
pbdMPI::init()
r <- pbdMPI::comm.rank(comm=comm)
cs <- pbdMPI::comm.size(comm=comm)
attr(comm,"rank") <- r
attr(comm,"size") <- cs

ca <- commandArgs(trailingOnly=TRUE)
if (length(ca)>0){
	N <- as.integer(ca[1])
} else {
	N <- 1000              # default sample size
}
h <- 1e-2                  # step size

beta <- (1.0 - (r/cs))^2   # PT: inverse temperature
modelName <- "AKAR4"
comment(modelName) <- "./AKAR4.so"
sb <- readRDS(file="AKAR4-sb.RDS")
ex <- readRDS(file="AKAR4-ex.RDS")

parMCMC <- log10(sb$Parameter[["!DefaultValue"]])
stdv <- rep(2,length(parMCMC))

dprior <- dNormalPrior(mean=parMCMC,sd=stdv)
rprior <- rNormalPrior(mean=parMCMC,sd=stdv)

## ----simulate-----------------------------------------------------------------
sim <- simcf(ex,modelName,log10ParMap) # or simulator.c

## ----likelihood---------------------------------------------------------------
logLH <- function(y,h,stdv,name){
	n <- sum(!is.na(stdv))
	llf_const <- sum(log(stdv),na.rm=TRUE) + 0.5*log(2*pi)*n
	llf_sq <- 0.5*sum(((y - h)/stdv)^2,na.rm=TRUE)
	return(-llf_const-llf_sq)
}

suppressMessages(
	llf <- logLikelihoodFunc(ex,simpleUserLLF=logLH)
)
## ----update-------------------------------------------------------------------
metropolis <- mcmcUpdate(
	simulate=sim,
	ex=ex,
	logLikelihood=llf,
	dprior=dprior)

## ----mcmc-method--------------------------------------------------------------
mh <- mcmc(metropolis)
ptMetropolis <- mcmc_mpi(metropolis,comm=comm,swapDelay=0,swapFunc=pbdMPI_bcast_reduce_temperatures)
## ----init---------------------------------------------------------------------
x <- mcmcInit(
	beta,
	parMCMC,
	simulate=sim,
	logLikelihood=llf,
	dprior)

A <- function(a) { # step-size adjuster
    return(0.5 + a^4/(0.25^4 + a^4))
}

## ----converge-and-adapt-------------------------------------------------------
for (i in seq(30)){
	s <- mh(x,200,h)           # evaluate rate of acceptance
	a <- attr(s,"acceptanceRate")
	h <- h*A(a)                # adjust h up or down
	x <- attr(s,"lastPoint")   # start next iteration from last point
	cat(sprintf("rank: %02i; iteration: %02i; a: %f; h: %g\n",r,i,a,h))
}

pbdMPI::barrier()
x <- mcmcInit(beta,x,simulate=sim,logLikelihood=llf,dprior)
s <- ptMetropolis(x,N,h) # the main amount of work is done here

saveRDS(s,file=sprintf("sample-rank-%i.RDS",r))
# ---- when all are done, we load the sampled points from the files but only for the right temperature:
pbdMPI::barrier()
f <- dir(pattern=sprintf('sample-rank-.*RDS$'))
X <- uqsa::gatherSample(f,beta)
saveRDS(X,file=sprintf("AKAR4-parameter-sample-for-rank-%i.RDS",r))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f %% and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
