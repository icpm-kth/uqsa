## ----setup--------------------------------------------------------------------
#library(uqsa)
library(parallel)
library(rgsl)
library(hexbin)
library(SBtabVFGEN)

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
load("goodMCMCvalue.RData",verbose=TRUE) # parStart
parVal <- parStart
#parVal <- log10(head(AKAP79_default(),-n))
print(parVal)


## ----range--------------------------------------------------------------------
defRange <- 2 # log-10 space
dprior <- dNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))
rprior <- rNormalPrior(mean=parVal,sd=rep(defRange,length(parVal)))

## ----simulate-----------------------------------------------------------------
sensApprox <- sensitivityEquilibriumApproximation(experiments, model, log10ParMap, log10ParMapJac)
simulate <- simulator.c(experiments,modelName,log10ParMap,noise=FALSE,sensApprox)
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
m <- mcmc(update)   # a Markov chain
h <- 1e-3           # step size guess

## ----adjust-------------------------------------------------------------------
accTarget <- 0.25
L <- function(a) { (1.0 / (1.0+exp(-(a-accTarget)/0.1))) + 0.5 }
N <- 100

start_time <- Sys.time()
x <- parVal
                              # do the adjustment of h a few times
options(mc.cores = parallel::detectCores())
for (j in seq(8)){
 cat("adjusting step size: ",h," \n");
 x <- mcmcInit(1.0,x,simulate,dprior,llf,gradLL,fiIn)
 Sample <- m(x,N,eps=h)
 a <- attr(Sample,"acceptanceRate")
 cat("acceptance: ",a*100," %\n")
 h <- h * L(a)
 x <- attr(Sample,"lastPoint")
}
plot(attr(Sample,"logLikelihood"),xlab="iteration",ylab="log-likelihood",main="small Sample to find a good step size",type='l')
cat("final step size: ",h,"\n")
cat("finished adjusting after",difftime(Sys.time(),start_time,units="sec")," seconds\n")


## ----initCluster--------------------------------------------------------------
n <- 8                                          # cluster size
nChains <- 32
options(mc.cores = parallel::detectCores() %/% n)
cl <- parallel::makeForkCluster(n)
parallel::clusterSetRNGStream(cl, 1337)          # seeding random numbers sequences

betas <- seq(1,0,length.out=nChains)^4
parMCMC <- lapply(betas,mcmcInit,parMCMC=parVal,simulate=simulate,dprior=dprior,logLikelihood=llf,gradLogLikelihood=gradLL,fisherInformation=fiIn)


## ----sample-------------------------------------------------------------------

start_time <- Sys.time()                         # measure sampling time
Sample <- NULL
for (i in seq(100)){
 s <- parallel::parLapply(cl, parMCMC, m, N=100, eps=h)
 parMCMC <- lapply(s,attr,which="lastPoint")
 parMCMC <- swap_points(parMCMC)
 if (i>2) {
  Sample <- rbind(Sample,s[[1]])
 }
}

colnames(Sample) <- names(parVal)


time_ <- difftime(Sys.time(),start_time,units="sec")
parallel::stopCluster(cl)
print(time_)


## ----hexplom, fig.width = 12, fig.height = 12---------------------------------

print(tail(Sample,10))
hexbin::hexplom(Sample)

