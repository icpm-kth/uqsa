mcmcUpdate <- function(parProposal,parGiven){
	return(NULL)
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
	FI <- function(p, sensitivity=NULL, simulations=NULL, simulate=NULL){
		np <- length(p)
		fi  <- matrix(0.0,np,np)
		if (missing(simulate)) stopifnot(is.list(simulations))
		if (is.null(simulations)) simulations <- simulate(p)
		for (i in seq(length(experiments))){
			errF <- t(experiments[[i]]$errorValues)
			errF[is.na(errF)] <- Inf
			oT <- experiments[[i]]$outputTimes
			for (j in seq(length(oT))){
				u <- experiments[[i]]$input
				modelPar <- c(parMap(p),u)
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
