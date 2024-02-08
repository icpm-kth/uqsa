#' This function can be used to specify default values
#'
#' When attributes are missing, the `base::attr()` function returns
#' NULL. In those cases this function can be used to find an
#' alternative value in one expression:
#' `attr(x,"dim") %otherwise% length(x)`
#' @param a value to check for NULL
#' @param b value to substitute
#' @return a, or b if a is NULL
#' @export
`%otherwise%` <- function(a,b){
	if (is.null(a) || any(is.na(a))) {
		return(b)
	} else {
		return(a)
	}
}

#' Initialize the Markov chain
#'
#' This function must append all required attributes to the MCMC
#' varible, for the Markov chain to update correctly.
#'
#' @param beta inverse temperature for the Markov chain (parallel tempering)
#' @param parMCMC a plain starting value for the Markov chain
#' @param logLikelihood a function that maps simulations to
#'     logLikelihood values
#' @param gradLogLikelihood the gradient function of the logLikelihood
#'     (optional) -- only if the algorithm requires it
#' @param fisherInformation a function that calculates the Fisher Information matrix
#' @return the same starting parameter vector, but with attributes.
#' @export
mcmcInit <- function(beta=1.0,parMCMC,simulate,dprior,logLikelihood,gradLogLikelihood=NULL,fisherInformation=NULL){
	simulations <- simulate(parMCMC)
	attr(parMCMC,"beta") <- beta
	attr(parMCMC,"simulations") <- simulations
	attr(parMCMC,"prior")  <- dprior(parMCMC)
	attr(parMCMC,"logLikelihood") <- beta*logLikelihood(parMCMC)
	if(!is.null(fisherInformation)) attr(parMCMC,"fisherInformation") <- beta*beta*fisherInformation(parMCMC)
	if(!is.null(gradLogLikelihood)) attr(parMCMC,"gradLogLikelihood") <- beta*gradLogLikelihood(parMCMC)
	return(parMCMC)
}

#' Swap the end-points of two Markov chains
#'
#' This is a conditional swap, according to the rules of parallel
#' tempering.
#'
#' @export
#' @param parMCMC a list of Markov chain end points
#' @return a list with some members swapped
swap_points <- function(parMCMC){
	for (j in seq(1,nChains,by=2)){
		L1 <- attr(parMCMC[[j  ]],"logLikelihood")
		L2 <- attr(parMCMC[[j+1]],"logLikelihood")
		B1 <- attr(parMCMC[[j  ]],"beta")
		B2 <- attr(parMCMC[[j+1]],"beta")
		a <- (B2-B1)*(L1-L2)
		r <- runif(1)
		if (r<a){
			X <- parMCMC[[j]]
			parMCMC[[j]] <- parMCMC[[j+1]]
			parMCMC[[j+1]] <- X
			# update all attributes
			attr(parMCMC[[j]],"beta") <- B1
			attr(parMCMC[[j+1]],"beta") <- B2
		}
	}
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
		ll <- numeric(N)
		a <- 0
		for (i in seq(N)){
			parMCMC <- update(parMCMC,eps)
			ll[[i]] <- attr(parMCMC,"logLikelihood")
			sample[i,] <- parMCMC
			a <- a + as.numeric(attr(parMCMC,"accepted"))
		}
		attr(sample,"acceptanceRate") <- a/N
		attr(sample,"logLikelihood") <- ll
		attr(sample,"lastPoint") <- parMCMC
		return(sample)
	}
}

#' SMMALA move
#'
#' The Simiplified Manifold Metropolis Adjusted Langevin Algorithm uses a
#' move instriction that uses a Gaussian kernel that is shifted away
#' from the current point
#'
#' @export
#' @param beta inverse temperature (parallel tempering)
#' @param parGiven given point
#' @return SMMALA proposal point
smmala_move <- function(beta=1.0,parGiven,fisherInformationPrior,eps=1e-2){
	fiGiven <- attr(parGiven,"fisherInformation")
	gradLGiven <- attr(parGiven,"gradLogLikelihood")
	G <- (beta^2*fiGiven)+fisherInformationPrior
	g <- solve(G,beta*gradLGiven)
	parProposal <- mvtnorm::rmvnorm(1,
		mean=parGiven+0.5*eps*g,
		sigma=eps*eps*G
	)
 return(as.numeric(parProposal))
}

#' SMMALA transition kernel density
#'
#' The Simiplified Manifold Metropolis Adjusted Langevin Algorithm uses a
#' move instriction that uses a Gaussian kernel that is shifted away
#' from the current point.
#'
#' @export
#' @param beta inverse temperature (parallel tempering)
#' @param parGiven given point
#' @return SMMALA proposal point
smmala_move_density <- function(beta=1.0,parProposal,parGiven,fisherInformationPrior,eps=1e-2){
	fiGiven <- attr(parGiven,"fisherInformation")
	gradLGiven <- attr(parGiven,"gradLogLikelihood")
	G <- (beta^2*fiGiven)+fisherInformationPrior
	g <- solve(G,beta*gradLGiven)
	return(mvtnorm::dmvnorm(as.numeric(parProposal),
		mean=parGiven+0.5*eps*g,
		sigma=eps*eps*G)
	)
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
mcmcUpdate <- function(simulate, experiments, model, logLikelihood, gradLogLikelihood, fisherInformation, fisherInformationPrior, dprior){
	if (is.null(gradLogLikelihood)){
		cat("unhandled case.\n")
		U <- NULL
	} else {
	U <- function(parGiven, eps=1e-4, beta=1.0){
		r <- runif(1)
		fp <- fisherInformationPrior
		beta <- attr(parGiven,"beta") %otherwise% 1.0
		llGiven <- attr(parGiven,"logLikelihood")
		priorGiven <- attr(parGiven,"prior")
		n <- length(parGiven)
		## the very important step: suggest a successor to parGiven and simulate the model
		parProposal <- smmala_move(beta,parGiven,fp,eps)
		attr(parProposal,"simulations") <- simulate(parProposal)
		llProposal <- logLikelihood(parProposal)
		priorProposal <- dprior(parProposal)
		attr(parProposal,"beta") <- beta
		attr(parProposal,"logLikelihood") <- llProposal
		attr(parProposal,"prior") <- priorProposal
		attr(parProposal,"fisherInformation") <- fisherInformation(parProposal)
		attr(parProposal,"gradLogLikelihood") <- gradLogLikelihood(parProposal)

		fwdDensity <- smmala_move_density(beta,parProposal,parGiven,fp,eps)
		bwdDensity <- smmala_move_density(beta,parGiven,parProposal,fp,eps)

		L <- exp(beta*(llProposal - llGiven))
		P <- priorProposal/priorGiven
		K <- bwdDensity/fwdDensity

		if (r < L * P * K){
			attr(parProposal,"accepted") <- TRUE
			return(parProposal)
		} else {
			attr(parGiven,"accepted") <- FALSE
			return(parGiven)
		}
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

