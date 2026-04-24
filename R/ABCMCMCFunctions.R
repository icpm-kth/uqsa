#' Performs and Approximate Bayesian Computation Sampling of Model Parameters
#'
#' Given a set of simulation experiments (list), a model, parameter
#' boundaries, this function will draw a sample of parameters from the
#' posterior probability density of the given problem.
#'
#' Normally, ABC would produce a highly auto-correlated sample,
#' wasting lots of disk-space; `batchSize` can be used to thin out the
#' sample, recording every batchSize-th point: with `batchSize=100`,
#' we perform 100 updates to the Markov chain variable and then save
#' the state to the sample. Higher batchSize numbers improve the
#' apparent quality of the sample, but create more work per sampled
#' point.
#'
#' @export
#' @param objectiveFunction function that, given a (vectorial)
#'     parameter as input, simulates the model, and outputs the
#'     distance between experimental data and data simulated from the
#'     model with the parameter provided in input
#' @param startPar starting value for the parameter vector
#' @param nSims requested sample size
#' @param Sigma0 multivariate normal covariance of Markov chain
#'     transition kernel
#' @param delta ABC acceptance threshold
#' @param dprior a function that returns prior probability density
#'     values
#' @param batchSize number of points to produce to record one point
#' @param allow.reg allow regularization (logical), if TRUE, then
#'     Sigma will be made smaller once a very low acceptance rate is
#'     detected: one accepted update per batch
#' @return a list containing a sample matrix and a vector of scores
#'     (values of delta for each sample)
#' @examples
#' \dontrun{
#'   f <- uqsa_example("AKAR4")
#'   m <- model_from_tsv(f)
#'   o <- as_ode(m)
#'   ex <- experiments(m,o)
#'   C <- generate_code(o)
#'   c_path(o) <- write_c_code(C)
#'   so_path(o) <- shlib(o)
#'   s <- simulator.c(ex,o)
#'   objFunc <- makeObjective(ex,s)
#'   startPar <- values(m$Parameter)
#'   lowerBound <- startPar - m$Paramster$stdv # m$Parameter$min
#'   upperBound <- startPar + m$Paramster$stdv # m$Parameter$max
#'   abcSample <- ABCMCMC(
#'     objFunc,
#'     startPar=values(m$Parameter),
#'     100,
#'     cov(rprior(1000)),
#'     delta=1,
#'     dprior=dUniformPrior(lowerBound,upperBound),
#'     batchSize = 100
#'   )
#' }
ABCMCMC <- function(objectiveFunction, startPar, nSims, Sigma0, delta, dprior, batchSize = 100, parAcceptable=\(p){all(is.finite(p))},allow.reg = FALSE){
	cat("Started chain.\n")
	Sigma1 <- diag(0.25*diag(Sigma0))
	curDelta <- Inf
	np <- length(startPar)
	## current vlaues
	curPar  <- startPar
	curDelta <- max(objectiveFunction(curPar))
	stopifnot(is.finite(curDelta))

	curPrior <- dprior(curPar)
	draws <- matrix(NA, nSims, np, dimnames=list(NULL,names(startPar)))
	scores <- rep(NA, nSims)

	nRegularizations <- 0
	accRate <- 0.0
	for (i in seq(nSims)) {
		a <- 0
		message(i)
		for (j in seq(batchSize)) {
			## candidate values:
			if (runif(1)<=0.95) {
				canPar <- MASS::mvrnorm(n=1, curPar, Sigma0)
			} else {
				canPar <- MASS::mvrnorm(n=1, curPar, Sigma1)
			}
			out <- parUpdate(
				objectiveFunction,
				curPar,
				canPar,
				curDelta,
				curPrior,
				delta,
				dprior,
				parAcceptable
			)
			curPar <- out$curPar
			curDelta <- out$curDelta
			curPrior <- out$curPrior
			a <- a + out$acceptance
		}
		accRate <- accRate + a
		if (a <= 1 && as.logical(allow.reg)){
			nRegularizations <- nRegularizations + 1
			Sigma0 <- solve(
				solve(Sigma0)+solve(0.1*norm(Sigma0)*diag(1,np,np))
			)
			Sigma1 <- diag(0.25*diag(Sigma0))
		}
		draws[i,]  <- curPar
		scores[i] <- curDelta
	}
	accRate <- accRate/(nSims*batchSize)
	return(
		list(
			draws = draws,
			scores = scores,
			acceptanceRate = accRate,
			nRegularizations = nRegularizations
		)
	)
}


#' Updates Parameter Values
#'
#' under valid ABC update conditions (successful simulation) the
#' parameters are updated to new values.
#'
#' @noRd
#' @param objectiveFunction function that, given a (vectorial)
#'     parameter as input, simulates the model, and outputs the
#'     distance between experimental data and data simulated from the
#'     model with the parameter provided in input
#' @param curPar current parameter values (as ABC samples them)
#' @param canPar candidate parameter values (for MCMC)
#' @param curDelta current distance between data and simulation, if
#'     the MCMC chain has not yet reached any point where this is
#'     below the threshold (delta), this can be accepted as the new
#'     current state for the chain.
#' @param curPrior current Prior values given curPar
#' @param delta distance threshold for ABC
#' @param dprior prior probability density function
#' @param parAcceptable user-specified constraint function, must return a scalar Boolean
#' @return updated values for curPar, curDelta, and curPrior
#' @examples
#' \dontrun{
#'   parUpdate(objectiveFunction, curPar, canPar, curDelta, curPrior, delta, dprior, parAcceptable)
#' }
parUpdate <- function(objectiveFunction, curPar, canPar, curDelta, curPrior, delta, dprior, parAcceptable){
	## candidate's prior
	canPrior <- dprior(canPar)
	acceptance <- (runif(1) <= canPrior/curPrior) && parAcceptable(canPar)
	if (acceptance){
		canDelta <- max(objectiveFunction(canPar))
		if (is.na(canDelta)) {
			warning("[parUpdate] candidate's distance-value is NA. Replacing it with Inf")
			canDelta <- Inf
		}
		if (canDelta <= max(delta, curDelta)){
			curDelta <- canDelta
			curPrior <- canPrior
			curPar <- canPar
		} else {
			acceptance <- FALSE
		}
	}
	return(
		list(
			curPar=curPar,
			curDelta=curDelta,
			curPrior=curPrior,
			acceptance=acceptance
		)
	)
}


#' ABC acceptance of currently sampled values given old data (Prior)
#'
#' The prior probability density model using copulas and vines is not
#' perfect, so values sampled from an imperfect prior estimate can be
#' checked against old data.
#'
#' @export
#' @param draws matrix of sampled values (to be filtered).
#' @param objectiveFunction function that, given a (vectorial)
#'     parameter as input, simulates the model, and outputs the
#'     distance between experimental data and data simulated from the
#'     model with the parameter provided in input
#' @param delta the acceptance threshold.
#' @return a filtered subset of acceptable parameter draws
#' @examples
#' \dontrun{
#'   posterior <- ABCMCMC(...)
#'   passed <- checkFitWithPreviousExperiments(posterior$draws, objFunc, delta=0.5)
#'   posterior$draws <- passed
#' }
checkFitWithPreviousExperiments <- function(draws, objectiveFunction, delta){
	cat("\n-Checking fit with previous data\n")
	nDraws = dim(draws)[1]
	scores <- objectiveFunction(t(draws))
	if (is.matrix(scores)){
		acceptable <- apply(scores <= delta,2,all)
	} else {
		acceptable <- (as.numeric(scores) <= delta)
	}
	stopifnot(length(acceptable)==nDraws)
	if (any(acceptable)){
		draws <- draws[acceptable,]
		nPickedDraws <- sum(acceptable)
		nonFits <- nDraws - nPickedDraws;
		cat("-- ", nonFits, " samples  did not fit previous datasets")
	} else {
		print(scores)
		warning("none of the draws have been accepted.")
	}
	return(draws)
}
