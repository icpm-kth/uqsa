#!/usr/bin/env Rscript

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
ex <- readRDS("AKAP79-ex.RDS")

f <- uqsa_example(modelName)
m <- model_from_tsv(f)

parMCMC <- values(m$Parameter)
mu <- m$Parameter$median
stdv <- m$Parameter$stdv

stopifnot(length(mu)==length(stdv))
stopifnot(length(parMCMC)==length(mu))

dprior <- dNormalPrior(mean=mu,sd=stdv)
rprior <- rNormalPrior(mean=mu,sd=stdv)
gprior <- gNormalPrior(mean=mu,sd=stdv)

## ----simulate-----------------------------------------------------------------
sim <- simulator.c(ex,modelName,log10ParMap)

## ----update-------------------------------------------------------------------
smmala <- smmala_update(
	simulate=sim,
	dprior=dprior,
	gradLogLikelihood=gprior,
	fisherInformation=fi,
	fisherInformationPrior=diag(1.0/stdv^2)
)

## ----mcmc-method--------------------------------------------------------------
MC <- mcmc(smmala)
ptSMMALA <- mcmc_mpi(
	smmala,
	comm=comm,
	swapDelay=0,
	swapFunc=pbdMPI_bcast_reduce_temperatures
)

h <- tune_step_size(MC)

for (j in seq(2)){
	pbdMPI::barrier()
	x <- mcmc_init(
		beta,
		as.numeric(tail(X,1)), # last row
		simulate=sim,
		logLikelihood=llf,
		dprior,
		gllf,
		gprior,
		fi
	)
	## ---- here, the sampling happens
	s <- ptSMMALA(x,N,h) # the main amount of work is done here
	saveRDS(s,file=sprintf("AKAP79-pt-smmala-sample-%i-rank-%i.RDS",j,r))
	pbdMPI::barrier()
	f <- dir(pattern=sprintf("^AKAP79-pt-smmala-sample-%i-rank-.*RDS$",j))
	X <- uqsa::gatherSample(f,beta)
	saveRDS(X,file=sprintf("AKAP79-temperature-ordered-pt-smmala-sample-%i-for-rank-%i.RDS",j,r))
}
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
