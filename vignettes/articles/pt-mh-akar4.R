#!/usr/bin/env Rscript
library(parallel)
library(uqsa)
library(pbdMPI)

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

beta <- (1.0 - (r/cs))^2   # PT: inverse temperature

m <- readRDS(file="/dev/shm/AKAR4-m.RDS")
o <- readRDS(file="/dev/shm/AKAR4-o.RDS")
ex <- readRDS(file="/dev/shm/AKAR4-e.RDS")

parMCMC <- log10(values(m$Parameter))
stdv <- rep(2,length(parMCMC))
dprior <- dNormalPrior(mean=parMCMC,sd=stdv)
rprior <- rNormalPrior(mean=parMCMC,sd=stdv)

sim <- simfi(ex,o,log10ParMap,omit=2) # or simulator.c

metropolis <- metropolis_update(
	simulate=sim,
	dprior=dprior
)

rwm <- mcmc(metropolis) # non-parallel

ptMetropolis <- mcmc_mpi( # parallel
	metropolis,
	comm=comm
)

x <- mcmc_init(
	beta,
	parMCMC,
	simulate=sim,
	dprior=dprior
)

h <- tune_step_size(rwm,x,1e-3)
gc()
Z <- ptMetropolis(x,N,h) # the main amount of work is done here
saveRDS(Z,file=sprintf("/dev/shm/sample-rank-%i.RDS",r))
pbdMPI::barrier()
rm(Z)

f <- dir("/dev/shm",pattern=sprintf('sample-rank-.*RDS$'),full.names=TRUE)
X <- gatherSample(f,beta)
saveRDS(X,file=sprintf("/dev/shm/AKAR4-parameter-sample-for-rank-%i.RDS",r))
gc()

time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(
	sprintf(
		"rank %02i of %02i finished with an acceptance rate of %f %% and swap rate of %f.\n",
		round(r),
		round(cs),
		attr(Z,"acceptanceRate"),
		attr(Z,"swapRate")
	)
)

finalize()
