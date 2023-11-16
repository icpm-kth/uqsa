#' This function proposes an MCMC candidate variable, and either accepts or rejects the candidate
#'
#' This function receives a current MCMC variable, then calculates a
#' possible successor and returns it in the case of acceptance. It
#' returns the (old) current state upon rejection of the candidate.
mcmcUpdate <- function(parGiven, logLikelihood, dprior){
	r <- runif(1)
	LGiven <- attr(parGiven,"logLikelihood")
	LProposal <- logLikelihood(parProposal)
	attr(parProposal,"logLikelihood") <- LProposal
	if (r < exp(LProposal - LGiven) * (dprior(parProposal)/dprior(parGiven))){
		return(parProposal)
	} else {
		return(parGiven)
	}
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

fisherInformation <- function(model, experiments, parMap=identity, sensitivityMap="log10"){
	nF <- length(model$func(0.0,model$init(),model$par()))
	l10 <- log(10)
	FI <- function(parMCMC, sensitivity=NULL, simulations=NULL, simulate=NULL){
		np <- length(parMCMC)
		fi  <- matrix(0.0,np,np)
		if (missing(simulate)) stopifnot(is.list(simulations))
		if (is.null(simulations)) simulations <- simulate(parMCMC)
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
	return(FI)
}

#' Equilibrium state approximation of the solution sensitivity for ODE systems
#'
#' In this context, the sensitivity S(t;x,p) is dx(t;p)/dp, where
#' x(t;p) is the parameterized solution to an initial value problem
#' for ordinary differential equations: x'=f(t,x;p). In cases where
#' you have a proxy variable for p, e.g. r=log(p), the chain rule
#' applies. Similarly, we also have an output sensitivity for the
#' function g(x(t;p)).  The equilibrium approximation is correct for
#' state-variable values close to an equilibrium point q(p)
#' (fixed-point): f(t,q(p);p)=0.
#'
#' Typically, the sensitivity needs to be known at different
#' time-points t_k, therefore, this function returns S as a
#' 3-dimensional array S[i,j,k], where the index k corrsponds to time
#' t_k.
#'
#' This function requires pracma::expm to work.
#' 
#' @param experiments a list of simulation experiments
#' @param simulations an equivalent list of simulation results, for one parameter vector
#' @param model a list of functions for the model the experiments are applicable to
#' @param parMCMC the parameters that are used in Markov chain Monte Carlo as the MC variable
#' @param parMap a map to transform parMCMC into p, parameters the model accepts
#' @return S, the state sensitivity matrix length(x) × length(p) × length(t)
sensitivityEquilibriumApproximation <- function(experiments, simulations, model, parMCMC=NULL, parMap=identity){
	S <- list()
	y0 <- model$init(0.0)
	n  <- length(y0)

	for (i in seq(length(experiments))){
		t <- experiments[[i]]$outputTimes
		t0 <- experiments[[i]]$initialTime
		if (missing(parMCMC)) {
			p <- model$par()
		} else {
			p <- c(parMap(parMCMC),experiments[[1]]$input)
		}
		S[[i]] <- array(0.0,dim=c(n,length(p),length(t)))
		for (j in seq(length(t))){
			if (abs(t[j]-t0)>1e-15){
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

logLikelihood <- function(experiments,simulations){
	N <- length(experiments)
	dimFunc <- dim(simulations[[1]]$func)
	n <- dimFunc[3]
	m <- dimFunc[1]*dimFunc[2]
	logNormalizingConstant <- -0.5*(m*log(2*pi)+sum(experiments[[i]]$errorValue^2))
	L <- logNormalizingConstant
	for (i in 1:N){
		for (k in 1:n){
			L <- L - 0.5*sum(((t(experiments[[i]]$outputValues) - simulations[[i]]$func[,,k])/t(experiments[[i]]$errorValue))^2)
		}
	}
	return(L)
}

