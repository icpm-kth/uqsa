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
beta <- (1.0 - (r/cs))^2   # PT: inverse temperature

stepSize <- function(beta){
	return(7e-3*exp(-9*beta))
}

a <- commandArgs(trailingOnly=TRUE)
if (length(a)>0){
      N <- as.integer(a[1])
} else {
      N <- 50000           # default sample size
}

h <- 0.01                  # step size
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"
sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

for (i in seq_along(ex)){
	D <- t(ex[[i]]$outputValues)
	SD <- t(ex[[i]]$errorValues)
	D[is.na(D)] <- 0.0
	ex[[i]]$data <- D
	ex[[i]]$stdv <- D*0.05+apply(D,1,FUN=max,na.rm=TRUE)*0.05
	ex[[i]]$stdv[is.na(SD)] <- Inf
}

parMCMC <- log10(sb$Parameter[["!DefaultValue"]])
stopifnot("!Max" %in% names(sb$Parameter))
stopifnot("!Min" %in% names(sb$Parameter))
stdv <- 0.5*(log10(sb$Parameter[["!Max"]])-log10(sb$Parameter[["!Min"]]))
stopifnot(length(parMCMC)==length(stdv))

dprior <- dNormalPrior(mean=parMCMC,sd=stdv)
rprior <- rNormalPrior(mean=parMCMC,sd=stdv)
gprior <- \(p) {return(-1.0*(p-parMCMC)/stdv^2)}
## ----simulate-----------------------------------------------------------------
sim <- simfi(ex,modelName,log10ParMap) # or simulator.c

## log-likelihood
llf <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(Reduce(\(a,b) a + b$logLikelihood[1], attr(parMCMC,"simulations"), init = 0.0))
}

## gradient of the log-likelihood
gllf <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(as.numeric(Reduce(\(a,b) a + b$gradLogLikelihood[i,1], attr(parMCMC,"simulations"), init = 0.0) %*% J))
}

## Fisher Information
fi <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(t(J) %*% Reduce(\(a,b) a + b$FisherInformation[i,i,1], attr(parMCMC,"simulations"), init = 0.0) %*% J)
}

## ----update-------------------------------------------------------------------
smmala <- mcmcUpdate(
	simulate=sim,
	experiments=ex,
	logLikelihood=llf,
	dprior=dprior,
	gradLogLikelihood=gllf,
	gprior,
	fisherInformation=fi,
	fisherInformationPrior=diag(1.0/stdv^2)
)

## ----mcmc-method--------------------------------------------------------------
MC <- mcmc(smmala)
ptSMMALA <- mcmc_mpi(smmala,comm=comm,swapDelay=0,swapFunc=pbdMPI_bcast_reduce_temperatures)
## ----init---------------------------------------------------------------------
	x <- mcmcInit(
		beta,
		parMCMC,
		simulate=sim,
		logLikelihood=llf,
		dprior,
		gllf,
		gprior,
		fi
	)

## beta,parMCMC,simulate,logLikelihood,dprior,gradLogLikelihood=NULL,gprior=NULL,fisherInformation=NULL
A <- function(a) { # step-size adjuster
    return(0.5 + a^4/(0.25^4 + a^4))
}

## ----converge-and-adapt-------------------------------------------------------
if (r==0){
	cat(sprintf("%10s  %12s %16s %16s\n","rank","iteration","acceptance","step size"))
	cat(sprintf("%10s  %12s %16s %16s\n","----","---------","----------","---------"))
}

for (j in seq(2)){
	for (i in seq(6)){
		s <- MC(x,100,h)           # evaluate rate of acceptance
		a <- attr(s,"acceptanceRate")
		h <- h*A(a)                # adjust h up or down
		x <- attr(s,"lastPoint")   # start next iteration from last point
		cat(sprintf("%10i  %12i %16.4f %16.4g\n",r,i,a,h)) # cat(sprintf("rank: %02i; iteration: %02i; a: %f; h: %g\n",r,i,a,h))
	}
	saveRDS(s,file=sprintf("smmala-last-adaptation-sample-for-rank-%i.RDS",r))
	pbdMPI::barrier()
	## ---- here, the sampling happens
	s <- ptSMMALA(x,N,h) # the main amount of work is done here
	saveRDS(s,file=sprintf("AKAP79-smmala-sample-%i-rank-%i.RDS",j,r))
	## ---- when all are done, we load the sampled points from the files but only for the right temperature:
	pbdMPI::barrier()
	f <- dir(pattern=sprintf("^AKAP79-smmala-sample-%i-rank-.*RDS$",j))
	X <- uqsa::gatherSample(f,beta)
	saveRDS(X,file=sprintf("AKAP79-temperature-ordered-smmala-sample-%i-for-rank-%i.RDS",j,r))
	x <- mcmcInit(
		beta,
		as.numeric(tail(X,1)), # last row
		simulate=sim,
		logLikelihood=llf,
		dprior,
		gllf,
		gprior,
		fi
	)
}
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
