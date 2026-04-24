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

beta <- (1.0 - (r/cs))^2 # PT: inverse temperature

stepSize <- function(beta){
	return(7e-3*exp(-9*beta))
}

a <- commandArgs(trailingOnly=TRUE)
if (length(a)>0){
	N <- as.integer(a[1])
} else {
	N <- 50000 # default sample size
}

## /dev/shm doesn't exist on all OS types. On GNU/Linux it is a
## location that exists physically in RAM. That means that it doesn't
## survive reboots, which is great for out purposes.  If you want to
## adapt this file, then this has to be any predictable path, that
## works on your system, even the current directory `./`.
o <- readRDS("/dev/shm/AKAP79-ode.RDS")
ex <- readRDS("/dev/shm/AKAP79-ex.RDS")

f <- uqsa_example("AKAP79")
m <- model_from_tsv(f)

stopifnot(all(m$Parameter$scale == "log10"))
parMCMC <- values(m$Parameter)
mu <- m$Parameter$median
stdv <- m$Parameter$stdv

stopifnot(length(mu)==length(stdv))
stopifnot(length(parMCMC)==length(mu))

dprior <- dNormalPrior(mean=mu,sd=stdv)
rprior <- rNormalPrior(mean=mu,sd=stdv)
gprior <- gNormalPrior(mean=mu,sd=stdv)

## ----simulate-----------------------------------------------------------------
sim <- simfi(ex,o,log10ParMap)

## ----update-------------------------------------------------------------------
smmala <- smmala_update(
	simulate=sim,
	logLikelihood=ll,
	gradLogLikelihood=gllf(log10ParMapJac),
	dprior=dprior,
	gprior=gprior,
	fisherInformationPrior=diag(1.0/stdv^2),
	fisherInformation=fi(log10ParMapJac)
)

## ----mcmc-method--------------------------------------------------------------
MC <- mcmc(smmala)
ptSMMALA <- mcmc_mpi(
	smmala,
	comm=comm,
	swapDelay=0,
	swapFunc=pbdMPI_bcast_reduce_temperatures
)

x <- mcmc_init(
	beta,
	values(m$Parameter),
	simulate=sim,
	logLikelihood=ll,
	gradLogLikelihood=gllf(log10ParMapJac),
	dprior=dprior,
	gprior=gprior,
	fisherInformation=fi(log10ParMapJac)
)

#h <- tune_step_size(MC,x)
h <- 1e-6

for (j in seq(2)){
	pbdMPI::barrier()
	## ---- here, the sampling happens
	X <- ptSMMALA(x,N,h) # the main amount of work is done here
	saveRDS(X,file=sprintf("AKAP79-pt-smmala-sample-%i-rank-%i.RDS",j,r))
	pbdMPI::barrier()
	f <- dir(pattern=sprintf("^AKAP79-pt-smmala-sample-%i-rank-.*RDS$",j))
	Z <- uqsa::gatherSample(f,beta)
	saveRDS(Z,file=sprintf("AKAP79-temperature-ordered-pt-smmala-sample-%i-for-rank-%i.RDS",j,r))
	x <- mcmc_init(
		beta,
		tail(Z,1),
		simulate=sim,
		logLikelihood=ll,
		gradLogLikelihood=gllf(log10ParMapJac),
		dprior=dprior,
		gprior=gprior,
		fisherInformation=fi(log10ParMapJac)
	)
}
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
