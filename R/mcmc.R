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

sensitivityEquilibriumApproximator <- function(experiments, simulations, model, par=NULL, parMap=identity){
	S <- list()
	y0 <- model$init(0.0)
	n  <- length(y0)

	for (i in seq(length(experiments))){
		t <- experiments[[i]]$outputTimes
		t0 <- experiments[[i]]$initialTime
		if (missing(par)) {
			p <- model$par()
		} else {
			p <- c(parMap(par),experiments[[1]]$input)
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

