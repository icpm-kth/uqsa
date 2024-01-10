mcmc <- function(simulate,experiments,update,logLikelihood,gradLogLikelihood,fisherInformation,sensApprox,dprior,model){
	M <- function(parMCMC=model$par(),N=1000){
		simulations <- simulate(parMCMC)
		sample <- matrix(NA,N,length(parMCMC))
		attr(parMCMC,"logLikelihood") <- logLikelihood(simulations)
		attr(parMCMC,"fisherInformation") <- fisherInformation(parMCMC)
		attr(parMCMC,"gradLogLikelihood") <- gradLogLikelihood(parMCMC,simulations)

		for (i in seq(N)){
			simulations <- simulate(parMCMC)
			sensitivity <- sensApprox(parMCMC,simulations)
			parMCMC <- update(parMCMC,simulations,sensitivity)
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
#' @param model
#' @param experiments
#' @param parMap
#' @param parMapJac
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

#' Equilibrium state approximation of the solution sensitivity for ODE systems
#'
#' In this context, the sensitivity S(t;x,p) is dx(t;p)/dp, where
#' x(t;p) is the parameterized solution to an initial value problem
#' for ordinary differential equations and t is the independent
#' varibale: x'=f(t,x;p), where «'» indicates the derivative with
#' respect to t. In cases where you have a proxy variable for p,
#' e.g. r=log(p), the chain rule applies. Similarly, we also have an
#' output sensitivity for the function g(x(t;p)).  The equilibrium
#' approximation is exact for state-variable values close to an
#' equilibrium point q(p) (fixed-point): f(t,q(p);p)=0.
#'
#' Typically, the sensitivity needs to be known at different
#' time-points t_k, therefore, this function returns S as a
#' 3-dimensional array S[i,j,k], where the index k corrsponds to time
#' t_k; the closer x(t_k) is to equilibrium, the better the
#' approximation; near the initial state, the sensitivity is also
#' correct (only the intermediate time-span is approximate).
#'
#' This function requires pracma::expm to work.
#'
#' @param experiments a list of simulation experiments
#' @param simulations an equivalent list of simulation results, for one parameter vector
#' @param model a list of functions for the model the experiments are applicable to
#' @param parMCMC the parameters that are used in Markov chain Monte Carlo as the MC variable
#' @param parMap a map to transform parMCMC into p, parameters the model accepts
#' @return S, the state sensitivity matrix length(x) × length(p) × length(t)
sensitivityEquilibriumApproximation <- function(experiments, simulations, model, parMCMC, parMap=identity){
	y0 <- model$init(0.0)
	n  <- length(y0)
	SEA <- function(parMCMC,simulations){
		S <- list()
		for (i in seq(length(experiments))){
			t <- experiments[[i]]$outputTimes
			t0 <- experiments[[i]]$initialTime
			p <- c(parMap(parMCMC),experiments[[1]]$input)
			S[[i]] <- array(0.0,dim=c(n,length(p),length(t)))
			for (j in seq(length(t))){
				if (abs(t[j]-t0) > 1e-15){
					u <- experiments[[i]]$input
					u_pos <- seq(length(p)-length(u)+1,length(p))
					p[u_pos] <- u
					A <- model$jac(t[j], simulations[[i]]$state[,j,1], p)
					B <- model$jacp(t[j], simulations[[i]]$state[,j,1], p)
					AB <- solve(A,B)
					C <- (pracma::expm((t[j]-t0)*A) %*% AB) - AB
					S[[i]][,,j] <- C
				}
			}
		}
		return(S)
	}
	return(SEA)
}

#' funcSensitivity transforms state sensitivity into the sensitivity of output functions
#'
#' The state sensitivity matrix:
#'
#'            d state(time[k],state, param)[i]
#' S[i,j,k] = --------------------------------  ,
#'            d param[j]
#'
#' where param are the raw model parameters.
#' This matrix is calculated as an intermediate and then transformed into:
#'
#'             d func(time[k], state, c(parMap(parMCMC),input))[i]
#' Sh[i,j,k] = --------------------------------------------------
#'             d parMCMC[j]
#'
#' where parMCMC is the Markov chain variable and usually shorter than
#' param as we typically don't sample all of the model's
#' parameters. Some model parameters may be known, some may be input
#' parameters not intrinsic to the model but related to the
#' experimental setup (that is why parMCMC and param are different).
#'
#' This transformation requires the output function jacobian (funcJac)
#' and the parameter jacobian (funcJacp) in the model variable.
#'
#' As we transform the parameters themselves, the chain rule requests
#' parMapJac[l,k] = d param[l] / d parMCMC[k]
#'
#' @param parMCMC Markov chain variables
#' @param experiments list of experiments
#' @param simulations list of simulation results
#' @param model a list of model functions, including funcJac, and funcJacp
#' @param parMap a parameter transformation function
#' @param parMapJac transformation jacobian: t(∇) %*% parMap
#' @return Sh[[i]][,,k] a list of 3d-arrays, one item per expariment
#'     (i), and a third dimension for the time variable (finite set of
#'     observations)
#' @export
funcSensitivity <- function(parMCMC,experiments,simulations,model,parMap=defaultParMap,parMapJac=defaultParMapJac){
	N <- length(experiments)
	fSens <- function(parMCMC,experiments){
		lMCMC <- length(parMCMC) # the dimension of the MCMC variable (parMCMC)
		S  <- sensitivityEquilibriumApproximation(experiments, simulations, model, parMCMC, parMap)
		Sh <- list()
		for (i in seq(N)){
			T <- length(simulations[[i]]$time)
			F <- length(experiments[[i]]$outputValues)
			y <- t(experiments[[i]]$outputValues)
			parMapped <- parMap(parMCMC)
			parODE <- c(parMapped,experiments[[i]]$input)
			Sh[[i]] <- array(0,dim=c(F,lMCMC,T))
			for (j in seq(T)){
				tj <- experiments[[i]]$outputTimes[j]
				ShODE <- model$funcJac(tj,y[,j],parODE) %*% S[[i]][,,j] + model$funcJacp(tj,y[,j],parODE)
				Sh[[i]][,,j] <- ShODE[,seq(lMCMC),drop=FALSE] %*% parMapJac(parMCMC)
			}
		}
		return(Sh)
	}
	return(fSens)
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

