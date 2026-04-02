#' Performs and Approximate Bayesian Computation as a Particle Filter
#'
#' Given a set of simulation experiments (list), a model, parameter
#' boundaries, this function will draw a sample of parameters from the
#' posterior probability density of the given problem.
#'
#' This is a variant of ABC where the entire batch is simulated with
#' one call to the simulator. startPar is the initial batch to be
#' simulated: it is a matrix where columns are different parameter
#' vectors (e.g. prior sample members). In other words: `startPar[,i]`
#' must be a valid argument for the objectiveFunction.
#'
#' The objective function is a closure
#'
#' The Objective-Function `objectiveFuntion(P)` should return a matrix
#' with `n` rows, where `n` is the number of simulation experiments
#' (and thus data-sets), and `m` columns, where `m` is the number of
#' parameterizations `NCOL(P)`.
#'
#' @export
#' @param objectiveFunction a function that can simulate the model for
#'     a batch of parameter vectors provided as a matrix of columns
#'     (batches)
#' @param startPar a matrix that has the same shape as the desired
#'     sample, but transposed, this can be a sample from the prior or a
#'     pre-conditioned sample that approximates the posterior, e.g.: t(rprior(1000))
#' @param Sigma multivariate normal covariance of Markov chain
#'     transition kernel
#' @param delta ABC acceptance threshold, either a scalar, then it is
#'     the initial value of delta, or a pair of values, then it is the
#'     starting value and the final value of delta:
#'     `c(initialDelta,finalDelta)`
#' @param dprior a function that returns prior probability density
#'     values
#' @param parAcceptable is a rejection-shortcut function; if
#'     `parAcceptable(p)` returns `FALSE` for a specific value of `p`,
#'     it means that simulations shouldn't even be attempted.
#' @param messages a logical value indicating whether log messages should be printed
#' @return a list containing a sample matrix and a vector of scores
#'     (values of delta for each sample)
#' @examples
#' \dontrun{
#'   library(parallel)
#'   f <- uqsa_example("AKAR4")
#'   m <- model_from_tsv(f)
#'   ex <- experiments(m,as_ode(m,cla=FALSE))
#'   G <- as_cme(m)         # for Gillespie solver
#'   C <- generate_code(G)
#'   c_path(G) <- write_c_code(C)
#'   so_path(G) <- shlib(G)
#'   options(mc.cores=detectCores())
#'   muX <- m$Parameter$value
#'   sdX <- m$Parameter$stdv
#'   rprior <- rNormalPrior(log(muX^2/(muX^2+sdX^2)),sqrt(log(1+sdX^2/muX^2)))
#'   dprior <- dNormalPrior(log(muX^2/(muX^2+sdX^2)),sqrt(log(1+sdX^2/muX^2)))
#'   s <- simstoch(ex,G,logParMap)
#'   O <- makeObjective(ex,s)
#'   X <- rprior(1000)
#'   colnames(X) <- rownames(m$Parameter)
#'   posterior <- ABCSMC(O,t(X),Sigma=cov(X),dprior=dprior,delta=c(0.5,1.5))
#' }
ABCSMC <- function(objectiveFunction, startPar, Sigma=2*cov(startPar), dprior, delta=c(2,0.5),  parAcceptable=\(p){all(is.finite(p))}, messages=FALSE){
	delta <- sort(delta,decreasing=TRUE) # in case someone enters a range for delta, e.g. c(0.1,0.9) rather than c(initial,final)
	initialDelta <- delta[1]
	if (length(delta)>1){
		finalDelta <- delta[2]
	} else {
		finalDelta <- 0 # automatic lower bound
	}
	deltaLowerBound <- finalDelta
	if (messages) message(sprintf("allowed range for delta: [%g,%g]",deltaLowerBound,initialDelta))
	## initial delta, subject to change:
	delta <- initialDelta
	## prepare starting values
	stopifnot(is.matrix(startPar))
	batchSize <- NCOL(startPar)
	## calculate initial distances
	curPar  <- startPar
	curDistance <- colMeans(objectiveFunction(curPar))
	stopifnot(all(is.finite(curDistance)))
	curPrior <- dprior(curPar)
	curWeight <- numeric(batchSize) + 1.0 # init to 1.0
	acceptanceRate <- 1.0                 # startPar has 100% acceptance, we don't reject any of them
	while (delta > deltaLowerBound && acceptanceRate > 0.03) {
		if (messages) message(sprintf("delta: %g",delta))
		newPar <- matrix(NA,NROW(startPar),0)
		newPrior <- numeric(0)
		newDistance <- numeric(0)
		accepted <- 0
		proposed <- 0
		# This while loop aggregates a new batch of points, using the current delta
		while (NCOL(newPar) < batchSize) {
			n <- max(100,batchSize - NCOL(newPar)) # don't bother suggesting anything too small
			proposed <- proposed + n
			if (messages) message(sprintf("proposing %i new points.",n))
			## resample from previous batch:
			k <- sample(seq_along(curWeight),n,replace=TRUE,prob=curWeight)
			canPar <- curPar[,k]
			## jitter, using a Gaussian transition kernel
			canPar <- canPar + t(mvtnorm::rmvnorm(NCOL(canPar),numeric(NROW(startPar)),Sigma))
			l <- apply(canPar,2,parAcceptable) # hard-rejection
			if (!any(l)) next
			canPar <- canPar[,l,drop=FALSE]
			canPrior <- apply(canPar,2,dprior) # prior values
			l <- (canPrior>0.0)
			if (!any(l)) next
			canPar <- canPar[,l,drop=FALSE]
			canPrior <- canPrior[l]
			canDistance <- colMeans(objectiveFunction(canPar)) # vectorized
			l <- (canDistance < delta)
			if (any(l)){
				newPar <- cbind(newPar,canPar[,l,drop=FALSE])
				newPrior <- c(newPrior,canPrior[l])
				newDistance <- c(newDistance,canDistance[l])
				accepted <- accepted + sum(l)
			}
		}
		if (messages) message(sprintf("accepted: %i",accepted))
		if (messages) message(sprintf("proposed: %i",proposed))
		## newPar is now an aggregate new batch, possibly bigger than batchSize
		## we now calculate weights for newPar:
		t_curPar <- t(curPar)
		sum_W_K <- apply(newPar,2,\(z) {sum(curWeight*mvtnorm::dmvnorm(t_curPar,z,Sigma))})
		newWeight <- newPrior/sum_W_K
		newWeight <- newWeight/sum(newWeight)
		## save the generated batch as current values
		curPar <- newPar
		curDistance <- newDistance
		curWeight <- newWeight
		## update Sigma
		WSigma <- stats::cov.wt(t(curPar))
		Sigma <- 2*WSigma$cov
		## update delta
		k <- sample(seq_along(curWeight),batchSize,prob=curWeight,replace=TRUE)
		delta <- median(curDistance[k])
		acceptanceRate <- accepted/proposed
		if (messages) message(sprintf("acceptance rate: %i %%",round(acceptanceRate*100)))
	}
	draws <- t(curPar[,k])
	colnames(draws) <- rownames(startPar)
	return(
		list(
			draws = draws,
			scores = curDistance[k],
			acceptanceRate = acceptanceRate
		)
	)
}

