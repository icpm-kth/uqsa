mcmcUpdate <- function(parProposal,parGiven){
	return(NULL)
}


fisherInformationFromGSA <- function(Sample,yf,E){
	funcDim <- dim(yf[[1]]$func)
	mcmcDim <- dim(Sample)
	nF <- funcDim[1]
	nT <- funcDim[2]
	nP <- mcmcDim[2]
	N <- mcmcDim[1] # sample size
	stopifnot(N>nP)
	n <- length(yf)
	fi  <- matrix(0.0,nP,nP)
	for (i in 1:n){
		outF <- aperm(yf[[i]]$func,c(3,1,2))
		errF <- t(E[[i]]$errorValues)
		errF[is.na(errF)] <- Inf
		for (j in 1:nT){
			sigma <- matrix(errF[,j],nF,nP)
			S <- sensitivity(Sample,outF[,,j])/sigma
			fi  <-  fi + t(S) %*% S
		}
	}
	return(fi)
}


metricTensorApprox <- function(Sample,yf,E,model){
	fi <- fisherInformationFromGSA(Sample,yf,E)
	mcmcDim <- dim(Sample)
	n <- mcmcDim[2]
	G <- function(sim,p,a){
		for (i in 1:length(E)){
			pu <- c(parMap(p),E[[i]]$input)
			t <- E[[i]]$outputTimes
			for (j in 1:length(t)){
				K <- model$jacp(t[j],sim[[i]]$state[,j,1],pu)[,1:n] # this should be the output jacobian (we don't make one of those yet) and perhaps also Hessian?
				fi  <- fi + a * t(K) %*% K
			}
		}
		return(solve(fi))
	}
	return(G)
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

logLikelihoodHessian <- function(model,experiments,simulation) {
	N <- length(experiments)
	dimFunc <- dim(simulations[[1]]$func)
	n <- dimFunc[3]
	m <- dimFunc[1]*dimFunc[2]
	logNormalizingConstant <- -0.5*(m*log(2*pi)+sum(experiments[[i]]$errorValue^2))
	L <- logNormalizingConstant
	for (i in 1:N){
		for (k in 1:n){
			H[i,j] <- H[i,j] #+ theActualSecondDerivative
		}
	}
	return(L)
}
