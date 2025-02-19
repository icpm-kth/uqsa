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
N <- 600000                 # default sample size
h <- 1e-2                  # step size

beta <- (1.0 - (r/cs))^2   # PT: inverse temperature
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"
sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

parMCMC <- c(0.298079363763137, 1.92414841193426, -3.46777437505845, -0.526142972399993, -0.335094425404467, 0.0154698674191168, 0.0589600659244421, -3.61722729859103, -2.81780027714603, -0.6294158469834, -0.337804704712963, -0.321516588604802, 1.32711476805281, -2.46115660220286, -2.43983426718858, 2.23812031125173, -3.72472670364874, -0.436575562361984, 2.33085374648653, -0.185827208095606, -1.06956748735169, -0.371748694335462, -0.810380219155343, 2.76171905370505, 2.09976590481578, 1.72733346358505, -2.19495503976565)
##c(1.42128435379381, 1.07608222486463, -3.47763108283383, -0.561580401738087, -0.346454801930942, -0.498138250397358, -0.28269649369305, -3.87563598331343, -2.43857253090425, -0.362298924611361, -0.262577607410472, 0.0228356788871791, 1.01918914320521, -1.80055534720865, -2.08958368679738, 0.982054225368052, -4.33148775774964, -0.374742076645538, 1.39602114864247, -0.405627329177627, -1.09200949057625, -0.645750168550101, -0.986412740227694, 2.46204002900275, 2.20222596616132, 1.50140437474347, -1.91950096195717)
## c(1.55904088262934, 0.8557303195304, -3.20117169640499, -0.478110748350246, -0.299286851364919, -1.1437050857187, -0.689170405609668, -4.0284739833275, -2.00978901806362, -0.821962702575839, -0.276461476232757, 0.375784872580963, 0.535304235795207, -1.8431373767986, -1.98970077563667, 1.62296477655548, -3.65076632672971, -0.421460188026877, 1.49604603864883, -0.802882209110794, 0.353833845351874, -1.00930900172998, -1.07519136390445, 1.88239613138114, 2.01849600590187, 0.8710450448686, -0.927898760548087)
## log10(sb$Parameter[["!DefaultValue"]])
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
saveRDS(X,file=sprintf("AKAP79-parameter-sample-for-rank-%i.RDS",r))
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f %% and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
