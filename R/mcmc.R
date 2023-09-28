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

fisherInformation <- function(experiments,simulations,model,sensitivity=NULL, par=NULL, parMap=identity){
	stopifnot(is.list(S))
	stopifnot(is.list(experiments))
	stopifnot(length(S)==length(experiments))
	nP <- dim(S[[1]])[2]
	nF <- length(model$func(0.0,model$init(),model$par()))
	fi  <- matrix(0.0,nP,nP)
	if (missing(sensitivity)){
		sensitivity <- sensitivityEquilibriumApprox(experiments, simulations, model, par, parMap)
	}
	if (missing(par)) {
		p <- model$par()
	} else {
		p <- c(parMap(par),experiments[[1]]$input)
	}
	for (i in seq(length(S))){
		errF <- t(experiments[[i]]$errorValues)
		errF[is.na(errF)] <- Inf
		nT <- length(experiments[[i]]$outputTimes)
		t <- experiments[[i]]$outputTimes
		for (j in 1:nT){
			u <- experiments[[i]]$input
			u_pos <- seq(length(p)-length(u)+1,length(p))
			p[u_pos] <- u
			sigma <- matrix(errF[,j],nF,nP)
			x <- simulations[[i]]$state[,j,1]
			S <- (model$funcJac(t[j],x,p) %*% sensitivity[[i]][,,j] + model$funcJacp(t[j],x,p))/sigma
			fi  <-  fi + t(S) %*% S
		}
	}
	return(fi)
}

sensitivityEquilibriumApprox <- function(experiments, simulations, model, par=NULL, parMap=identity){
	S <- list()
	for (i in seq(length(experiments))){
		t <- experiments[[i]]$outputTimes
		t0 <- experiments[[i]]$initialTime
		if (missing(par)) {
			p <- model$par()
		} else {
			p <- c(parMap(par),experiments[[1]]$input)
		}
		for (j in seq(length(t))){
			u <- experiments[[i]]$input
			u_pos <- seq(length(p)-length(u)+1,length(p))
			p[u_pos] <- u
			A <- model$jac(t[j],simulations[[i]]$state[,j,1],p)
			B <- model$jacp(t[j],simulations[[i]]$state[,j,1],p)
			AB <- solve(A+diag(.Machine$double.eps,n,n),B)
			S[[i]][,,j] <- (pracma::expm((t[j]-t0)*A) %*% AB) - AB
		}
	}
	return(S)
}

sea <- sensitivityEquilibriumApprox
