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
	comm=comm
)

p0 <- c(
	1.535276, 0.2129433, -2.337059, -0.2758802, -1.097471, -2.148527, -2.310619, -3.046435,
	-2.11772, -0.9887686, -0.32303, 0.1178787, -0.9677111, -0.9030976, -1.682916, 1.469942,
	-3.538294, -0.53798, 1.320372, -0.3500506, 0.2507297, -1.737913, -0.9742049, 1.008909,
	2.003937, 0.002032132, -0.1675174
) # values(m$Parameter)

h <- stepSize(beta) # initially

for (j in seq(2)){
	## 1. initialize Markov chain
	x <- mcmc_init(
		beta,
		p0,
		simulate=sim,
		logLikelihood=ll,
		gradLogLikelihood=gllf(log10ParMapJac),
		dprior=dprior,
		gprior=gprior,
		fisherInformation=fi(log10ParMapJac)
	)
	## 2. adjust step soze
	h <- tune_step_size(MC,x,h=h) # adjust using test chain
	pbdMPI::barrier()
	## 3. Prepare a location for sample files, use '/dev/shm' if available
	TMPDIR <- ifelse(dir.exists('/dev/shm'),'/dev/shm',Sys.getenv(TMPDIR,unset='/tmp'))
	rank_file.rds <- tempfile(
		pattern=sprintf("AKAP79-pt-smmala-sample-%i-rank-%i",j,r),
		tmpdir=ifelse(dir.exists('/dev/shm'),'/dev/shm',TMPDIR),
		fileext='.rds'
	)
	## 4. sampling
	X <- ptSMMALA(x,N,h) # the main amount of work is done here
	## 5. save the sample to previously determined location TMPDIR
	saveRDS(X,file=rank_file.rds)
	pbdMPI::barrier()
	## 6. from all written files, load only the rows with consistent beta values
	f <- dir(TMPDIR,pattern=sprintf("^AKAP79-pt-smmala-sample-%i-rank-.*rds$",j),full.names=TRUE)
	Z <- uqsa::gatherSample(f,beta) # <- we only want this value of beta
	print(dim(Z))                   # little sanity check
	print(head(Z))                  #
	## 6. save the temperature consistent sample (constant beta) to a different file
	saveRDS(Z,file=file.path(TMPDIR,sprintf("AKAP79-temperature-ordered-pt-smmala-sample-%i-for-rank-%i.rds",j,r)))
	## 7. remove the non-ordered file
	file.remove(rank_file.rds)
	## 8. set a new starting location for the next iteration
	p0 <- as.numeric(tail(Z,1))
}
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(
	sprintf(
		"rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",
		round(r),
		round(cs),
		attr(Z,"acceptanceRate"),
		attr(Z,"swapRate")
	)
)
finalize()
