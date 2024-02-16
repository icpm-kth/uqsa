#!/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(hexbin)
library(SBtabVFGEN)
library(doMPI)

devtools::load_all("~/uqsa")
## ----label="SBtab content"----------------------------------------------------
modelFiles <- uqsa_example("AKAP79",pattern="[.]tsv$",full.names=TRUE)
SBtab <- SBtabVFGEN::sbtab_from_tsv(modelFiles)


## ----model--------------------------------------------------------------------
source(uqsa_example("AKAP79",pat="^AKAP79[.]R$"))
names(model)
# compile
modelName <- checkModel("AKAP79",uqsa_example("AKAP79",pat="_gvf[.]c$"))


## ----ConservationLaws---------------------------------------------------------
load(uqsa_example("AKAP79",pat="^ConservationLaws[.]RData$"))
print(ConLaw$Text)
experiments <- sbtab.data(SBtab,ConLaw)


## ----default------------------------------------------------------------------
n <- length(experiments[[1]]$input)
stopifnot(n>0)
#load("goodMCMCvalue.RData",verbose=TRUE) # parStart
#parVal <- parStart
parVal <- log10(head(AKAP79_default(),-n))
print(parVal)


## ----range--------------------------------------------------------------------
defRange <- 2 # log-10 space
dprior <- dNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
rprior <- rNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))

## ----simulate-----------------------------------------------------------------
sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
simulate <- simc(experiments,modelName,log10ParMap,sensApprox)
y <- simulate(parVal)

plot(experiments[[1]]$outputTimes,as.numeric(y[[1]]$state[1,,1]),xlab='time',ylab='AKAR4p', main='state',type='l')


## ----likelihood---------------------------------------------------------------
llf <- logLikelihood(experiments)
gradLL <- gradLogLikelihood(model,experiments, parMap=log10ParMap, parMapJac=log10ParMapJac)
fiIn <- fisherInformation(model, experiments, parMap=log10ParMap)
fiPrior <- solve(diag(defRange, length(parVal)))


## ----update-------------------------------------------------------------------
update  <- mcmcUpdate(simulate=simulate,
	experiments=experiments,
	model=model,
	logLikelihood=llf,
	gradLogLikelihood=gradLL,
	fisherInformation=fiIn,
	fisherInformationPrior=fiPrior,
	dprior=dprior)


## ----init---------------------------------------------------------------------
m <- mcmc(update)   # a Markov chain function
h <- 1e-3           # initial step size guess



## ----initCluster--------------------------------------------------------------
n <- 15                                          # cluster size, apart from the main process
nChains <- n+1
options(mc.cores = 1)

#cl <- parallel::makeForkCluster(n)
cl <- startMPIcluster(n)
registerDoMPI(cl)
#parallel::clusterSetRNGStream(cl, 1337)          # seeding random numbers sequences

betaSchedule <- seq(1,0,length.out=nChains)^4

## ----adjust-------------------------------------------------------------------
accTarget <- 0.25
L <- function(a) { (1.0 / (1.0+exp(-(a-accTarget)/0.1))) + 0.5 }
N <- 100

start_time <- Sys.time()
x <- parVal
                              # do the adjustment of h a few times
options(mc.cores = 1)


parMCMC <- foreach(b=betaSchedule) %dopar% {
	set.seed(1337*Rmpi::mpi.comm.rank(comm=0))
	for (j in seq(3)){
		x <- mcmcInit(beta=b,x,simulate,dprior,llf,gradLL,fiIn)
		Sample <- m(x,N,eps=h)
		a <- attr(Sample,"acceptanceRate")
		h <- h * L(a)
		x <- attr(Sample,"lastPoint")
	}
	attr(x,"stepSize") <- h
	return(x)
}
cat("final step size: ",h,"\n")
cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n")
cat("final parameter vectors, for each chain: \n")
print(Reduce(rbind,parMCMC))
save(parMCMC,file="mcmcStart.RData")
## ----sample-------------------------------------------------------------------
stopifnot(is.list(parMCMC) && length(parMCMC)==nChains)

start_time <- Sys.time()                         # measure sampling time

m <- mcmc_mpi(update)                            # m is now an mpi aware function

Sample <- foreach (p=parMCMC, .combine="rbind") %dopar% {
	h <- attr(p,"stepSize")
	print(names(attributes(p)))
	cat(sprintf("starting a Markov chain with step size h = %g.\n",h))
	s <- m(p,N=100,h)
	b <- attr(s,"beta")
	b1 <- (abs(b-1.0) < 1e-6)
	return(s[b1,])
}

colnames(Sample) <- names(parVal)
save(Sample,"doMPI-testSample.RData")

time_ <- difftime(Sys.time(),start_time,units="sec")
#parallel::stopCluster(cl)
closeCluster(cl)
mpi.finalize()
print(time_)


## ----hexplom, fig.width = 12, fig.height = 12---------------------------------

print(tail(Sample,10))
#hexbin::hexplom(Sample)

