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
mcmcInit <- function(parMCMC,logLikelihood,gradLogLikelihood=NULL,fisherInformation=NULL){
	attr(parMCMC,"logLikelihood") <- logLikelihood(simulations)
	if(!is.null(fisherInformation)) attr(parMCMC,"fisherInformation") <- fisherInformation(parMCMC)
	if(!is.null(gradLogLikelihood)) attr(parMCMC,"gradLogLikelihood") <- gradLogLikelihood(parMCMC,simulations)
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
#' @param simulate a function that simulates the model, given MCMC
#'     values
#' @param experiments list of experiments
#' @param update and update function
#' @param dprior probability density of the prior distribution
#' @param model a list of model functions (R functions) for the
#'     underlying ODE model
#' @return M(initiPar,N), a function of initial starting values and
#'     number of Markov chain steps
#' @export
mcmc <- function(simulate,experiments,update,dprior,model){
	M <- function(parMCMC,N=1000){
		simulations <- simulate(parMCMC)
		sample <- matrix(NA,N,length(parMCMC))
		for (i in seq(N)){
			simulations <- simulate(parMCMC)
			parMCMC <- update(parMCMC,simulations)
			sample[i,] <- parMCMC
		}
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
mcmcUpdate <- function(logLikelihood, gradLogLikelihood, fisherInformation, dprior, model, experiments, parMap=defaultParMap, parMapJac=defaultParMapJac,eps=1e-5){
	U <- function(parGiven){
		r <- runif(1)
		LGiven <- attr(parGiven,"logLikelihood")
		fiGiven <- attr(parGiven,"fisherInformation")
		gradLGiven <- attr(parGiven,"gradLogLikelihood")
		n <- length(parGiven)
		## the very important step: suggest a successor to parGiven
		parProposal <- mvtnorm::rmvnorm(1,
			mu=parGiven+0.5*eps*eps*solve(fiGiven,gradLGiven),
			sigma=eps*eps*fiGiven)
		## re-calculate things that depend on parProposal
		LProposal <- logLikelihood(parProposal)
		fiProposal <- fisherInformation(parProposal)
		gradLProposal <- gradLogLikelihood(parProposal)
		## in this specific case, the proposal is asymmetric, so we need the forward and backward transition densitity values
		fwdDensity <- mvtnorm::dmvnorm(parProposal,
			mu=parGiven+0.5*eps*eps*solve(fiGiven,gradLGiven),
			sigma=eps*eps*fiGiven)
		bwdDensity <- mvtnorm::dmvnorm(parGiven,
			mu=parProposal+0.5*eps*eps*solve(fiProposal,gradLProposal),
			sigma=eps*eps*fiProposal)
		if (r < exp(LProposal - LGiven) * (dprior(parProposal)/dprior(parGiven)) * (bwdDensity/fwdDensity)){
			attr(parProposal,"logLikelihood") <- LProposal
			attr(parProposal,"fisherInformation") <- fiProposal
			attr(parProposal,"gradLogLikelihood") <- gradLProposal
			return(parProposal)
		} else {
			return(parGiven)
		}
	}
	return(U)
}

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
fisherInformation <- function(model, experiments, parMap=defaultParMap, parMapJac=defaultParMapJac){
	nF <- length(model$func(0.0,model$init(),model$par()))
	l10 <- log(10)
	F <- function(parMCMC, simulations, sensitivity){
		np <- length(parMCMC)
		fi  <- matrix(0.0,np,np)
		for (i in seq(length(experiments))){
			errF <- t(experiments[[i]]$errorValues)
			errF[is.na(errF)] <- Inf
			oT <- experiments[[i]]$outputTimes
			for (j in seq(length(oT))){
				u <- experiments[[i]]$input
				modelPar <- c(parMap(parMCMC),u)
				sigma_j <- matrix(errF[,j],nF,np)
				x_ij <- simulations[[i]]$state[,j,1]
				# output function sensitivity:
				Sx <- sensitivity[[i]][,seq(np),j]
				funcJac <- model$funcJac(oT[j],x_ij,modelPar)
				funcJacp <- model$funcJacp(oT[j],x_ij,modelPar)[,seq(np)]
				S <- (funcJac %*% Sx + funcJacp)/sigma_j
				if (sensitivityMap=="log10") {
					P <- matrix(modelPar[seq(np)],nF,np,byrow=TRUE)
					S <- S*P*l10;
				} else if (sensitivityMap=="log"){
					P <- matrix(modelPar[seq(np)],nF,np,byrow=TRUE)
					S <- S*P;
				}
				fi  <-  fi + t(S) %*% S
			}
		}
		return(fi)
	}
	return(F)
}


logLikelihood <- function(experiments,simulations){
	N <- length(experiments)
	llf <- function(simulations){
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

defaultParMap <- function(parMCMC){
	return(10^(parMCMC))
}

defaultParMapJac <- function(parMCMC){
	return(diag(10^(parMCMC) * log(10)))
}

gradLogLikelihood <- function(model,experiments,parMap=defaultParMap,parMapJac=defaultParMapJac){
	N <- length(experiments)
	gradLL <- function(parMCMC,simulations){
		lMCMC <- length(parMCMC) # the dimension of the MCMC variable (parMCMC)
		gL <- rep(0,length(parMCMC))
		Sh <- funcSensitivity(parMCMC,experiments,simulations,model,parMap,parMapJac)
		for (i in seq(N)){
			T <- length(simulations[[i]]$time)
			y <- t(experiments[[i]]$outputValues)
			h <- simulations[[i]]$func[,,1]
			stdv <- t(experiments[[i]]$errorValues)
			for (j in seq(T)){
				gL <- gL + ((y[,j] - h[,j])/stdv[,j]^2) %*% Sh[[i]][,,j]
			}
		}
		return(gL)
	}
	return(gradLL)
}

