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

N <- 6000                 # default sample size
h <- stepSize(beta)       # step size
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"
sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

## here we create a secondary data structure that represents the data,
## but all NA values replaced with sensible values that lead to the sam eresults
## This way, the underlying C code won't need to filter out NAs
for (i in seq_along(ex)){
	D <- t(ex[[i]]$outputValues)
	SD <- t(ex[[i]]$errorValues)
	D[is.na(D)] <- 0.0
	SD[is.na(SD)] <- Inf
	ex[[i]]$data <- D
	ex[[i]]$stdv <- SD
}

parMCMC <- c(0.298079363763137, 1.92414841193426, -3.46777437505845, -0.526142972399993, -0.335094425404467, 0.0154698674191168, 0.0589600659244421, -3.61722729859103, -2.81780027714603, -0.6294158469834, -0.337804704712963, -0.321516588604802, 1.32711476805281, -2.46115660220286, -2.43983426718858, 2.23812031125173, -3.72472670364874, -0.436575562361984, 2.33085374648653, -0.185827208095606, -1.06956748735169, -0.371748694335462, -0.810380219155343, 2.76171905370505, 2.09976590481578, 1.72733346358505, -2.19495503976565)
stdv <- rep(2,length(parMCMC))

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
## beta,parMCMC,simulate,logLikelihood,dprior,gradLogLikelihood=NULL,gprior=NULL,fisherInformation=NULL
x <- mcmcInit(
	beta,
	parMCMC,
	simulate=sim,
	logLikelihood=llf,
	dprior,
	gllf,
	gprior,
	fi)

A <- function(a) { # step-size adjuster
    return(0.5 + a^4/(0.25^4 + a^4))
}

## ----converge-and-adapt-------------------------------------------------------
if (r==0){
	cat(sprintf("%10s  %12s %16s %16s\n","rank","iteration","acceptance","step size"))
	cat(sprintf("%10s  %12s %16s %16s\n","----","---------","----------","---------"))
}
for (i in seq(30)){
	s <- MC(x,200,h)           # evaluate rate of acceptance
	a <- attr(s,"acceptanceRate")
	h <- h*A(a)                # adjust h up or down
	x <- attr(s,"lastPoint")   # start next iteration from last point
	cat(sprintf("%10i  %12i %16.4f %16.4f\n",r,i,a,h)) # cat(sprintf("rank: %02i; iteration: %02i; a: %f; h: %g\n",r,i,a,h))
}
saveRDS(s,file=sprintf("smmala-last-adaptation-sample-for-rank-%i.RDS",r))
pbdMPI::barrier()
## ---- here, the sampling happens
s <- ptSMMALA(x,N,h) # the main amount of work is done here

saveRDS(s,file=sprintf("AKAP79-smmala-sample-rank-%i.RDS",r))
## ---- when all are done, we load the sampled points from the files but only for the right temperature:
pbdMPI::barrier()
f <- dir(pattern='^AKAP79-smmala-sample-rank-.*RDS$')
X <- uqsa::gatherSample(f,beta)
saveRDS(X,file=sprintf("AKAP79-tempertaure-ordered-smmala-sample-for-rank-%i.RDS",r))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f %% and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
