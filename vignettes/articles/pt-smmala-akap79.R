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

parMCMC <- c(0.298079363763137, 1.92414841193426, -3.46777437505845, -0.526142972399993, -0.335094425404467, 0.0154698674191168, 0.0589600659244421, -3.61722729859103, -2.81780027714603, -0.6294158469834, -0.337804704712963, -0.321516588604802, 1.32711476805281, -2.46115660220286, -2.43983426718858, 2.23812031125173, -3.72472670364874, -0.436575562361984, 2.33085374648653, -0.185827208095606, -1.06956748735169, -0.371748694335462, -0.810380219155343, 2.76171905370505, 2.09976590481578, 1.72733346358505, -2.19495503976565)
stdv <- rep(2,length(parMCMC))

dprior <- dNormalPrior(mean=parMCMC,sd=stdv)
rprior <- rNormalPrior(mean=parMCMC,sd=stdv)
gprior <- \(p) {return(-1.0*(p-parMCMC)/stdv^2)}
## ----simulate-----------------------------------------------------------------
sim <- simc(ex,modelName,log10ParMap) # or simulator.c

## ----likelihood---------------------------------------------------------------
logLH <- function(y,h,stdv,name){
	n <- sum(!is.na(stdv))
	llf_const <- sum(log(stdv),na.rm=TRUE) + 0.5*log(2*pi)*n
	llf_sq <- 0.5*sum(((y - h)/stdv)^2,na.rm=TRUE)
	return(-llf_const-llf_sq)
}

gllf <- gradLogLikelihoodFunc(ex,parMap=log10ParMap,parMapJac=log10ParMapJac)
fi <- fisherInformationFunc(ex, parMap=log10ParMap, parMapJac=log10ParMapJac)


suppressMessages(
	llf <- logLikelihoodFunc(ex,simpleUserLLF=logLH)
)
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
for (i in seq(1)){
	s <- MC(x,200,h)           
	a <- attr(s,"acceptanceRate") # evaluate rate of acceptance
	h <- h*A(a)                # adjust h up or down
	x <- attr(s,"lastPoint")   # start next iteration from last point
	cat(sprintf("rank: %02i; iteration: %02i; a: %f; h: %g\n",r,i,a,h))
}
saveRDS(s,file=sprintf("smmala-last-adaptation-sample-for-rank-%i.RDS",r))
pbdMPI::barrier()
## ---- here, the sampling happens
s <- ptSMMALA(x,N,h) # the main amount of work is done here

saveRDS(s,file=sprintf("AKAP79-smmala-sample-rank-%i.RDS",r))
## ---- when all are done, we load the sampled points from the files but only for the right temperature:
pbdMPI::barrier()
f <- dir(pattern=sprintf('AKAP79-smmala-sample-rank-.*RDS$'))
X <- uqsa::gatherSample(f,beta)
saveRDS(X,file=sprintf("AKAP79-ordered-smmala-sample-for-rank-%i.RDS",r))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f %% and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
