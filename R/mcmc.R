#' Initialize the Markov chain
#'
#' This function must append all required attributes to the MCMC
#' varible, for the Markov chain to update correctly.
#'
#' @param parMCMC a plain starting value for the Markov chain
#' @param logLikelihood a function that maps simulations to
#'     logLikelihood values
#' @param gradLogLikelihood the gradient function of the logLikelihood
#'     (optional) -- only if the algorithm requires it
#' @param fisherInformation a function that calculates the Fisher Information matrix
#' @return the same starting parameter vector, but with attributes.
#' @export
mcmcInit <- function(parMCMC,simulate,logLikelihood,gradLogLikelihood=NULL,fisherInformation=NULL){
	simulations <- simulate(parMCMC)
	attr(parMCMC,"simulations") <- simulations
	attr(parMCMC,"logLikelihood") <- logLikelihood(parMCMC)
	if(!is.null(fisherInformation)) attr(parMCMC,"fisherInformation") <- fisherInformation(parMCMC)
	if(!is.null(gradLogLikelihood)) attr(parMCMC,"gradLogLikelihood") <- gradLogLikelihood(parMCMC)
	return(parMCMC)
}

#' Markov Chain Monte Carlo
#'
#' This function creates an MCMC function for a given set of
#' experiments.
#'
#' The algorithm is entirely determined by the update function.  Any
#' intermediate values that updates requires aside from simulation
#' results have to be attributes of the MCMC variable: parMCMC.
#'
#' The update function: update(parGiven) -> parUpdate depends only on
#' the given parameters, all other dependencies have to be either
#' implicit (as a closure) or attributes of parGiven.
#'
#' @param update and update function
#' @return M(initiPar,N), a function of initial starting values and
#'     number of Markov chain steps
#' @export
mcmc <- function(update){
	M <- function(parMCMC,N=1000,eps=1e-4){
		sample <- matrix(NA,N,length(parMCMC))
		a <- 0
		for (i in seq(N)){
			parMCMC <- update(parMCMC,eps)
			sample[i,] <- parMCMC
			a <- a + as.numeric(attr(parMCMC,"accepted"))
		}
		attr(sample,"acceptanceRate") <- a/N
		return(sample)
	}
}

#' This function proposes an MCMC candidate variable, and either accepts or rejects the candidate
#'
#' This function receives a current MCMC variable, then calculates a
#' possible successor and returns it in the case of acceptance. It
#' returns the (old) current state upon rejection of the candidate.
#'
#' The Markov chain has a current state (the MCMC variable, often x in
#' literature), but in the context of sampling the MCMC variables are
#' used as the parameters to a scientific model of some sort (and
#' these often have state variables, also x, or y). This is why we
#' call the variables parMCMC (parABC), or
#' par{Current|Given|Proposal}, and similar.
#'
#' @export
#' @param logLikelihood function, returns log(p(simulations|experiments))
#' @param dprior prior probability density function
#' @param eps a step size parameter for Markov chain moves (propotional to step size)
#' @return a function that returns possibly updated states of the Markov chain
mcmcUpdate <- function(simulate, experiments, model, logLikelihood, gradLogLikelihood, fisherInformation, fisherInformationPrior, dprior, parMap=identity, parMapJac=1){
	U <- function(parGiven, eps=1e-4){
		r <- runif(1)
		LGiven <- attr(parGiven,"logLikelihood")
		fiGiven <- attr(parGiven,"fisherInformation")
		gradLGiven <- attr(parGiven,"gradLogLikelihood")
		fi <- fiGiven+fisherInformationPrior
		n <- length(parGiven)
		## the very important step: suggest a successor to parGiven and simulate the model
		parProposal <- as.numeric(mvtnorm::rmvnorm(1,
			mean=parGiven+0.5*eps*solve(fiGiven+fisherInformationPrior,gradLGiven),
			sigma=eps*eps*(fiGiven+fisherInformationPrior)))
		#cat("parProposal: ",parProposal)
		attr(parProposal,"simulations") <- simulate(parProposal)
		## re-calculate things that depend on parProposal
		LProposal <- logLikelihood(parProposal)
		fiProposal <- fisherInformation(parProposal)
		gradLProposal <- gradLogLikelihood(parProposal)
		## in this specific case, the proposal is asymmetric, so we need the forward and backward transition densitity values
		fwdDensity <- mvtnorm::dmvnorm(as.numeric(parProposal),
			mean=as.numeric(parGiven+0.5*eps*solve(fiGiven+fisherInformationPrior,gradLGiven)),
			sigma=eps*eps*solve(fiGiven+fisherInformationPrior))

		bwdDensity <- mvtnorm::dmvnorm(as.numeric(parGiven),
			mean=as.numeric(parProposal+0.5*eps*solve(fiProposal+fisherInformationPrior,gradLProposal)),
			sigma=eps*eps*solve(fiProposal+fisherInformationPrior))
		if (r < exp(LProposal - LGiven) * (dprior(parProposal)/dprior(parGiven)) * (bwdDensity/fwdDensity)){
			attr(parProposal,"logLikelihood") <- LProposal
			attr(parProposal,"fisherInformation") <- fiProposal
			attr(parProposal,"gradLogLikelihood") <- gradLProposal
			attr(parProposal,"accepted") <- TRUE
			return(parProposal)
		} else {
			attr(parGiven,"accepted") <- FALSE
			return(parGiven)
		}
	}
	return(U)
}

#' Calculate Global Fisher Information
#'
#' Given a sample, this performs global sensitivity analysis, and then
#' squares the sensitivity.
#'
#' @param Sample an MCMC sample, or ABC sample
#' @param yf the simulations of the sample
#' @param E experiments list
fisherInformationFromGSA <- function(Sample,yf=NULL,E){
	funcDim <- dim(yf[[1]]$func)
	mcmcDim <- dim(Sample)
	nF <- funcDim[1]
	nT <- funcDim[2]
	nP <- mcmcDim[2]
	N <- mcmcDim[1] # sample size
	stopifnot(N>nP)
	n <- length(yf)
	fi  <- matrix(0.0,nP,nP)
	for (i in seq(n)){
		outF <- aperm(yf[[i]]$func,c(3,1,2))
		errF <- t(E[[i]]$errorValues)
		errF[is.na(errF)] <- Inf
		for (j in seq(nT)){
			sigma <- matrix(errF[,j],nF,nP)
			S <- globalSensitivity(Sample,outF[,,j])/sigma
			fi  <-  fi + t(S) %*% S
		}
	}
	return(fi)
}

#' Fisher Information from Sensitivity
#'
#' Given a list of simulation sensitivities, this function returns the
#' fisher information (sum over all experiments). The actual work is
#' done in the returned function that implicitly depends on the model,
#' experiments, and parameter mapping
#'
#' return value: function(par, simulations, sensitivity) ->
#' fisherInformation (matrix)
#'
#' where par refers to the model parameters
#' (possibly transformed), and simulations performed with those
#' parameters.
#'
#' @export
#' @param model list of R functions for the ODE model
#' @param experiments list of experiments, with inputs
#' @param parMap mapping between MCMC variables and ODE parameters
#' @param parMapJac the jacobian of the above map
#' @return fisher information calculating funciton
fisherInformation <- function(model, experiments, parMap=identity, parMapJac=1){
	nF <- length(model$func(0.0,model$init(),model$par()))
	l10 <- log(10)
	F <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		np <- length(parMCMC)
		fi  <- matrix(0.0,np,np)
		for (i in seq(length(experiments))){
			errF <- t(experiments[[i]]$errorValues)
			errF[is.na(errF)] <- Inf
			nt <- length(experiments[[i]]$outputTimes)
			for (j in seq(nt)){
				sigma_j <- matrix(errF[,j],nF,np)
				# output function sensitivity:
				Sh <- matrix(simulations[[i]]$funcsens[,seq(np),j],nF,np)/sigma_j
				fi  <-  fi + t(Sh) %*% Sh
			}
		}
		if (any(is.na(as.numeric(fi)))) {
			stop()
		}
		return(fi)
	}
	return(F)
}

#' Default log-likelihood function
#'
#' This returns a function f(simulations), which maps simulation
#' results to log(likelihood) values. The experiments are used
#' implicitly; simulations is a list as returned by
#' rgsl::r_gsl_odeiv2_outer().
#'
#' @param experiment will be compared tp the simulation results
#' @export
logLikelihood <- function(experiments){
	N <- length(experiments)
	llf <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		dimFunc <- dim(simulations[[1]]$func)
		n <- dimFunc[3]
		m <- dimFunc[1]*dimFunc[2]
		L <- rep(-0.5*m*N*log(2*pi),n)
		for (i in seq(N)){
			y <- t(experiments[[i]]$outputValues)
			stdv <- t(experiments[[i]]$errorValues)
			for (k in seq(n)){
				h <- simulations[[i]]$func[,,k]
				L[k] <- L[k] - 0.5*sum(((y - h)/stdv)^2,na.rm=TRUE) + sum(log(stdv))
			}
		}
		return(L)
	}
	return(llf)
}

#' LOG10 parameter mapping used by the MCMC module
#'
#' This map is used by the simulator to transform sampling variables
#' into ODE-model porameters.
#'
#' @param parMCMC the sampling variables (numeric vector)
#' @export
log10ParMap <- function(parMCMC){
	return(10^(parMCMC))
}

#' LOG10 parameter mapping, jacobian
#'
#' This map is used by the simulator to transform sampling variables
#' into ODE-model porameters. As we often calculate sensitivites, we
#' alos need the jacobian of the map, due to the chain rule of
#' differentiation.
#'
#' @param parMCMC the sampling variables (numeric vector)
#' @export
log10ParMapJac <- function(parMCMC){
	return(diag(10^(parMCMC) * log(10)))
}

#' Default log-likelihood function, gradient
#'
#' This returns a function g(x,simulations), which maps simulation
#' results and the MCMC variables x to the gradient of log(likelihood)
#' values withj respect to x. The experiments are used implicitly;
#' simulations is a list as returned by rgsl::r_gsl_odeiv2_outer().
#'
#' @param experiment will be compared tp the simulation results
#' @export
gradLogLikelihood <- function(model,experiments,parMap=identity,parMapJac=1){
	N <- length(experiments)
	gradLL <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		lMCMC <- length(parMCMC) # the dimension of the MCMC variable (parMCMC)
		gL <- rep(0,length(parMCMC))
		for (i in seq(N)){
			d <- dim(simulations[[i]]$func)
			T <- length(experiments[[i]]$outputTimes)
			y <- t(experiments[[i]]$outputValues)
			h <- simulations[[i]]$func[,,1]
			dim(h) <- head(d,2)
			stdv <- t(experiments[[i]]$errorValues)
			Sh <- simulations[[i]]$funcsens
			dS <- dim(Sh)
			for (j in seq(T)){
				Shj <- Sh[,,j]
				dim(Shj) <- head(dS,2)
				gL <- gL + as.numeric(t(((y[,j,drop=FALSE] - h[,j,drop=FALSE])/stdv[,j,drop=FALSE]^2)) %*% Shj)
			}
		}
		return(gL)
	}
	return(gradLL)
}

