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
N <- 50000                 # default sample size
h <- 1e-2                  # step size

beta <- (1.0 - (r/cs))^2   # PT: inverse temperature
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"
sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

p <- sb$Parameter[["!DefaultValue"]]
parMCMC <- log10(p)
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
	model=NULL,
	logLikelihood=llf,
	dprior=dprior)

## ----init---------------------------------------------------------------------
mh <- mcmc(metropolis)
ptMetropolis <- mcmc_mpi(metropolis,comm=comm,swapDelay=0,swapFunc=pbdMPI_bcast_reduce_temperatures)
## ----sample-------------------------------------------------------------------
x <- mcmcInit(
	beta,
	parMCMC,
	simulate=sim,
	logLikelihood=llf,
	dprior)

A <- function(a) { # step-size adjuster
    return(0.5 + a^4/(0.25^4 + a^4))
}

for (i in seq(10)){
	s <- mh(x,300,h)
	a <- attr(s,"acceptanceRate")
	h <- h*A(h)
	x <- attr(s,"lastPoint")
	if (r == 0) cat(sprintf("a: %f\n",a))
}

s <- ptMetropolis(x,N,h) # the main amount of work is done here
saveRDS(s,file=sprintf("sample-rank-%i.RDS",r))
# ---- when all are done, we load the sampled points from the files but only for the right temperature:
pbdMPI::barrier()
f <- dir(pattern=sprintf('sample-rank-.*RDS$'))
X <- uqsa::gatherSample(f,beta)
saveRDS(X,file=sprintf("AKAP79-parameter-sample-for-rank-%i.RDS",r))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %02i %% and swap rate of %f.\n",round(r),round(cs),round(attr(s,"acceptanceRate")),attr(s,"swapRate")))
finalize()
