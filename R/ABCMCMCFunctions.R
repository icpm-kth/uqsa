#' Performs and Approximate Bayesian Computation Sampling of Model Parameters
#'
#' ABC replaces the need for an exact likelihood function and uses a
#' distance function instead: the distance between data and
#' simulation. This distance function is very similar to a likelihood
#' function but lacks a statistical justification. Nevertheless, this
#' distance function, like the likelihood function of a deterministic
#' model, performs a simulation of the scientific model, be it fully
#' stochastic or a stochasically embedded, but deterministic in its
#' core.
#'
#' The distance of the ABC setting is compared to a threshold value
#' \eqn{\\delta}{delta}. The threshold doesn't need to be explicitly
#' provided. You can however provide a span of acceptabel values in
#' any order, the smaller value will be used as a lower bound, the
#' larger value will be used initially.
#'
#' The ABC procedure will attempt to converge first, using the initial
#' delta value, and decreasing it slowly using observed distance
#' values.
#'
#' This function always operates on a bundle of Markov chains. The
#' size of this bundle can be determined through `startPar` (a matrix
#' of column vectors). Each column will be used as the initial point
#' of a Markov chain. The chains will be resampled at each step during
#' the convergence phase, and decouple from one another once the
#' burn-in is complete.
#'
#' ABC methods (distance function, threshold delta) can be combined
#' with several other methds (like particle filters). Here we use
#' several parallle Markov chains to sample from the approximate
#' posterior.
#'
#' @export
#' @param objectiveFunction function that, given a parameter matrix as
#'     input, simulates the model, and outputs the distance between
#'     experimental data and data simulated from the model with the
#'     parameter provided in input. This has to be a closure that
#'     contains the experimental data within itself. The closure must
#'     be vectorized over the columns of its matrix argument.
#' @param startPar starting values for the parameter vector, can (and
#'     should) be a matrix with n columns, where each column is a
#'     valid parameter vector; the number of columns determines the
#'     batch size.
#' @param N requested number of batches to return, the sample will be
#'     of size `batchSize*N` (batch size is the number of parallel
#'     Markov chains).
#' @param burnIn number of batches where the transition kernel will be
#'     adjusted to achieve an acceptance rate of below 10%. These are
#'     printed as negative numbers during progress printouts:
#'     -N...0,1...N, where no adjustments to the MCMC parameters
#'     happen on positive loop-counts.
#' @param Sigma0 multivariate normal covariance of Markov chain
#'     transition kernel, defaults to the covariance of the initial
#'     parameters. If startPar is one vector, this matrix must be
#'     provided explicitly.
#' @param dprior a function that returns prior probability density
#'     values.
#' @param deltaSpan either an initial and final value for the ABC
#'     threshold delta, or a fixed value for delta that will never
#'     change.
#' @param batchSize the size of each batch, this should be a number
#'     that could be sufficient to calculate the covariance of in the
#'     given parameter space.
#' @param parAcceptable a function that can reject a parameter vector
#'     early based on user-requirements. Has to return a scalar
#'     Boolean. Use this to test for inequalities that you find
#'     difficult to encode in the prior.
#' @return a list containing a sample matrix and a vector of scores
#'     (values of delta for each sample)
#' @examples
#'   f <- uqsa_example("AKAR4")
#'   m <- model_from_tsv(f)
#'   o <- as_ode(m)
#'   ex <- experiments(m,o)
#'   C <- generate_code(o)
#'   c_path(o) <- write_c_code(C)
#'   so_path(o) <- shlib(o)
#'   s <- simulator.c(ex,o,parMap=log10ParMap)
#'   objFunc <- makeObjective(ex,s)
#'   p0 <- log10(values(m$Parameter))
#'   lowerBound <- p0 - 3
#'   upperBound <- p0 + 3
#'   dprior <- dUniformPrior(lowerBound,upperBound)
#'   rprior <- rUniformPrior(lowerBound,upperBound)
#'   X <- rprior(24)   # increase this number ...
#'   abcSample <- abc_mcmc(
#'     objFunc,
#'     startPar=t(X),
#'     N=10,           # ... and this number
#'     Sigma0=cov(X),
#'     dprior=dprior
#'   )
abc_mcmc <- function(objectiveFunction, startPar, N, burnIn=ceiling(sqrt(N)), Sigma0=cov(t(startPar)), dprior=NULL, deltaSpan=NULL,batchSize=100*NROW(Sigma0), parAcceptable=\(p){all(is.finite(p))}, verbose=TRUE){
	show_progress <- verbose && interactive()
	if (is.null(deltaSpan)){
		delta <- 2*max(objectiveFunction(startPar))
		delta_LB <- 0 # lower bound
	} else {
		delta <- max(deltaSpan)
		delta_LB <- min(deltaSpan)
	}
	np <- nrow(Sigma0)
	b <- sample.int(NCOL(startPar),size=batchSize,replace=TRUE)
	## make startPar _batch shaped_ (with batchSize columns), and add a little noise to it, in case it was exactly one vector
	curPar <- as.matrix(startPar)[,b] + matrix(rnorm(np*batchSize,0,norm(Sigma0)*0.01),np,batchSize)
	if (missing(dprior) || is.null(dprior)) { # construct something useful
		LB <- apply(curPar,1,min)
		UB <- apply(curPar,1,max)
		if (any(UB<=LB)){
			MD <- apply(curPar,1,median)
			LB <- MD - 2
			UB <- MD + 2
		}
		dprior <- dUniformPrior(LB,UB)
		warning("dprior is missing, will use uniform prior, with these bounds: ")
		print(data.frame(lower.bound=LB,upper.bound=UB))
	}
	curPrior <- dprior(t(curPar))
	curDistance <- NULL
	mu <- rowMeans(curPar)                # the posterior's mean (tracking)
	S <- cov(t(curPar)) * (batchSize - 1) # the posterior's covariance (tracking)
	I <- diag(nrow=nrow(Sigma0))
	Sigma1 <- diag(diag(Sigma0))
	batchSample <- matrix(nrow=batchSize,ncol=length(startPar))
	draws <- NULL
	distanceRecord <- NULL
	acceptanceRate <- NULL
	deltaRecord <- NULL
	n <- 1
	if (show_progress) {
		cli::cli_progress_bar(name="abc", total = N+burnIn+1)
	}
	for (i in seq(-burnIn,N)) {
		L <- chol(Sigma0)
		Z <- matrix(rnorm(batchSize*np),np,batchSize) # this shape is expected by the Objective
		canPar <- curPar + L %*% Z                    # canPar is a multivariate normal batch derived from the previous batch
		canPrior <- dprior(t(canPar))                 # a vector
		canWeight <- canPrior/curPrior                # element-wise division
		canWeight[!is.finite(canWeight)] <- 0
		l <- as.logical((runif(batchSize) <= canWeight) & apply(canPar,2,parAcceptable))   # l marks parameters that can be investigated further
		canWeight[!l] <- 0
		canDistance <- rep(Inf,batchSize)
		if (any(l)){
			canDistance[l] <- colMeans(objectiveFunction(canPar[,l,drop=FALSE])) # Obj is a matrix: nExperiment x length(l)
			canWeight <- canWeight * (canDistance <= delta) * is.finite(canDistance)
		}
		a <- mean(canWeight>0,na.rm=TRUE)                      # acceptance rate
		if (any(canWeight>0,na.rm=TRUE)){
			canWeight[is.na(canWeight)] <- 0
			if (i<=0){
				j <- sample.int(batchSize,size=batchSize,replace=TRUE,prob=canWeight)
				curDistance <- canDistance[j]
				curPrior <- canPrior[j]
				curPar <- canPar[,j]
			} else {
				j <- which(canWeight>0)
				curDistance[j] <- canDistance[j]
				curPrior[j] <- canPrior[j]
				curPar[,j] <- canPar[,j]
			}
		} else {
			a <- 0
		}
		if (show_progress){
			cli::cli_progress_update(inc=1,status=ifelse(i<=0,"burn-in","recording"))
		}
		if (i <= 0) {
			if (delta>delta_LB &&  a>0.1) { # adjust delta down
				md <- median(curDistance,na.rm=TRUE)
				delta <- max(delta_LB,md)
			}
			A <- a^2/(0.1^2 + a^2) + 0.5 # or exp(2.5 * (a - 0.10))
			w <- a/(0.1 + a)
			## calculate a new value for Sigma0, by tracking global covariance S (unscaled)
			mu_batch <- rowMeans(curPar)
			S_batch  <- cov(t(curPar)) * (batchSize - 1)
			Dmu <- (mu_batch - mu)
			S <- S + S_batch + (Dmu %*% t(Dmu)) * ((n/(n+1))*batchSize)
			mu <- mu + Dmu/(n+1)
			n <- n+1
			C <- S/(n*batchSize-1) # C is now the updated global covariance
			if (abs(det(C)) < 1e-6 + 1e-6*norm(C)) C <- I*norm(Sigma1) # reset if in danger
			Sigma0 <- A*((1-w)*Sigma0 + w*C + 1e-8*Sigma1)
		} else {
			draws  <- rbind(draws,t(curPar))
			distanceRecord <- c(distanceRecord,curDistance)
			acceptanceRate <- c(acceptanceRate,a)
			deltaRecord <- c(deltaRecord,delta)
		}
	}
	cli::cli_progress_done()
	return(
		list(
			draws = draws,
			distances = distanceRecord,
			acceptanceRate = acceptanceRate,
			delta = deltaRecord,
			Sigma = Sigma0
		)
	)
}
