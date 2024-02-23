#!/bin/env Rscript

## ----setup--------------------------------------------------------------------
library(uqsa)
library(parallel)
library(rgsl)
library(hexbin)
library(SBtabVFGEN)
library(Rmpi)
start_time <- Sys.time()                         # measure sampling time
a <- commandArgs(trailing=TRUE)
n <- 8
if (length(a)>0) n <- a[1]
mpi.spawn.Rslaves(nslaves=n)

retSample <- mpi.remote.exec({
	library(uqsa)
	library(parallel)
	library(rgsl)
	library(SBtabVFGEN)
	r <- mpi.comm.rank()
	cs <- mpi.comm.size()
	beta <- (1.0 - (r-1)/cs)^4
	nChains <- cs-1
	devtools::load_all("~/uqsa")

	## ## ----label="SBtab content"----------------------------------------------------
	modelFiles <- uqsa_example("AKAP79",pattern="[.]tsv$",full.names=TRUE)
	SBtab <- SBtabVFGEN::sbtab_from_tsv(modelFiles)

	## ----model--------------------------------------------------------------------
	source(uqsa_example("AKAP79",pat="^AKAP79[.]R$"))
	names(model)
	modelName <- checkModel("AKAP79","./AKAP79.so")

	## ----ConservationLaws---------------------------------------------------------
	load(uqsa_example("AKAP79",pat="^ConservationLaws[.]RData$"))
	print(ConLaw$Text)
	experiments <- sbtab.data(SBtab,ConLaw)


	## ----default------------------------------------------------------------------
	n <- length(experiments[[1]]$input)
	stopifnot(n>0)
	parVal <- log10(head(AKAP79_default(),-n))

	## ----range--------------------------------------------------------------------
	defRange <- 2 # log-10 space
	dprior <- dNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
	rprior <- rNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))

	## ----simulate-----------------------------------------------------------------
	sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
	simulate <- simc(experiments,modelName,log10ParMap,sensApprox)
	y <- simulate(parVal)

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
	## ----adjust-------------------------------------------------------------------
	accTarget <- 0.25
	L <- function(a) { (1.0 / (1.0+exp(-(a-accTarget)/0.1))) + 0.5 }
	N <- 100

	start_time <- Sys.time()
	x <- parVal
	parMCMC <- list()
	for (j in seq(nChains)){      # it doesn't have to be nChains, but we also collect possible starting posiitons
			x <- mcmcInit(beta,x,simulate,dprior,llf,gradLL,fiIn)
			Sample <- m(x,N,eps=h)
			a <- attr(Sample,"acceptanceRate")
			h <- h * L(a)
			x <- attr(Sample,"lastPoint")
			parMCMC[[j]] <- x
	}

	cat("final step size: ",h,"\n")
	cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n")

	## ----sample-------------------------------------------------------------------


	m <- mcmc_mpi(update)                            # m is now an mpi aware function

	r <- mpi.comm.rank()
	j <- (r%%nChains)+1
	cs <- mpi.comm.size()
	p <- parMCMC[[j]]
	attr(p,"beta") <- betaSchedule[j]
	h <- attr(p,"stepSize")
	s <- m(p,N=100,h)
	colnames(s) <- names(parVal)
	return(s)
},ret = TRUE)

betaTrace <- Reduce(function(a,b) return(c(a,attr(b,"beta"))),retSample,init=NULL)
Sample <- Reduce(rbind,retSample)
mpi.close.Rslaves()
save(retSample,file="retSample.RData")
dim(Sample)

save(Sample,file="Rmpi-testSample.RData")

time_ <- difftime(Sys.time(),start_time,units="sec")
mpi.finalize()
print(time_)


## ----hexplom, fig.width = 12, fig.height = 12---------------------------------

print(tail(Sample,10))
#hexbin::hexplom(Sample)

