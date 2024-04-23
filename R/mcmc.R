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


#' checks whether a variable has the named attributes
#'
#' @param var a variable to check for attributes
#' @param attrNames named attributes
#' @return TRUE if all attributes are present
#' @export
`%has%` <- function(var,attrNames){
	return(all(attrNames %in% names(attributes(var))))
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
mcmcInit <- function(beta,parMCMC,simulate,logLikelihood,dprior,gradLogLikelihood=NULL,gprior=NULL,fisherInformation=NULL){
	simulations <- simulate(parMCMC)
	attr(parMCMC,"beta") <- beta
	attr(parMCMC,"simulations") <- simulations
	attr(parMCMC,"prior")  <- dprior(parMCMC)
	attr(parMCMC,"logLikelihood") <- logLikelihood(parMCMC)
	if(!is.null(fisherInformation))
		attr(parMCMC,"fisherInformation") <- fisherInformation(parMCMC)
	if(!is.null(gradLogLikelihood))
		attr(parMCMC,"gradLogLikelihood") <- gradLogLikelihood(parMCMC)
	if(!is.null(gprior))
		attr(parMCMC,"gradLogPrior") <- gprior(parMCMC)
	return(parMCMC)
}

#' Should 2 Markov chains exchange their temperatures
#'
#' This function makes a Boolean choice about chnages in temperature,
#' based on the log(liklihood) values of two Markov chains in a
#' parallel tempering setting.
#'
#' This function is useful if `mpi.send()` and `mpi.recv()` are used.
#'
#' @param b1 the inverse temperature of chain 1
#' @param ll1 the log-liklihood of chain 1
#' @param b2 the inverse temperature of chain 2
#' @param ll2 the log-lilihood of chain 2
#' @return TRUE is the chains should swap their temperatures
#' @export
change_temperature <- function(b1,ll1,b2,ll2){
	a <- exp((b2-b1)*(ll1-ll2))
	r <- runif(1)
	return(r<a)
}

#' Swap the end-points of two Markov chains
#'
#' This is a conditional swap, according to the rules of parallel
#' tempering. This function is only useful if the Markov chains have
#' returned to the global scope and one process will make the decision
#' and perform the swap, i.e.: the current state of each chain is
#' locally available.
#'
#' @export
#' @param parMCMC a list of Markov chain end points, each entry
#'     annotated with a temperature attribute
#'     `attr(parMCMC[[i]],"beta")`
#' @return a list with some members swapped
swap_points_locally <- function(parMCMC){
	for (j in seq(1,nChains-1)){
		L1 <- attr(parMCMC[[j  ]],"logLikelihood")
		L2 <- attr(parMCMC[[j+1]],"logLikelihood")
		B1 <- attr(parMCMC[[j  ]],"beta")
		B2 <- attr(parMCMC[[j+1]],"beta")
		a <- exp((B2-B1)*(L1-L2))
		r <- runif(1)
		if (r<a){
			cat(sprintf("swapping %i with %i\n",j,j+1))
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
#' experiments. The Markov chains have no communication between them
#' if more than one is created using this mechanism.
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
		sample <- matrix(nrow=N,ncol=length(parMCMC))
		ll <- numeric(N)
		b <- numeric(N)
		a <- 0
		for (i in seq(N)){
			parMCMC <- update(parMCMC,eps)
			ll[[i]] <- attr(parMCMC,"logLikelihood")
			sample[i,] <- as.numeric(parMCMC)
			b[i] <- attr(parMCMC,"beta")
			a <- a + as.numeric(attr(parMCMC,"accepted"))
		}
		attr(sample,"acceptanceRate") <- a/N
		attr(sample,"logLikelihood") <- ll
		attr(sample,"lastPoint") <- parMCMC
		attr(sample,"beta") <- b
		attr(sample,"stepSize") <- eps
		return(sample)
	}
}

#' Communicate with other ranks and swap beta
#'
#' Given a current log-likelihood, temperature and step-size, this
#' decides whether to send or receive the same variables from a
#' neighboring process and swap temperatures with them.
#'
#' @param i MCMC iteration
#' @param B inverse temperature (parallel tempering)
#' @param LL log-likelihood value of current point
#' @param H algorithm's step size (often called epsilon in literature)
#' @param r MPI rank
#' @param comm MPI communicator
#' @param cs MPI comm size
#' @export
Rmpi_swap_temperatures <- function(i, B, LL, H, r, comm, cs){
	if (r %% 2 == i %% 2) { # alternating, to swap with r-1 or r+1
		Rmpi::mpi.send.Robj(LL, dest=(r+1)%%cs, tag=1, comm=comm)
		Rmpi::mpi.send.Robj(B, dest=(r+1)%%cs, tag=2, comm=comm)
		Rmpi::mpi.send.Robj(H, dest=(r+1)%%cs, tag=4, comm=comm)
	} else { # get e(x)ternal objects, from the other chain
		xLL <- Rmpi::mpi.recv.Robj(source=(r-1) %% cs, tag=1, comm=comm)
		xB <- Rmpi::mpi.recv.Robj(source=(r-1) %% cs, tag=2, comm=comm)
		xH <- Rmpi::mpi.recv.Robj(source=(r-1) %% cs, tag=4, comm=comm)
	}
	if (r %% 2 == i %% 2){
		B <- Rmpi::mpi.recv.Robj(source=(r+1)%%cs, tag=3, comm=comm)
		H <- Rmpi::mpi.recv.Robj(source=(r+1)%%cs, tag=5, comm=comm)
	} else if(change_temperature(B,LL,xB,xLL)){
		Rmpi::mpi.send.Robj(B,dest=(r-1) %% cs, tag=3, comm=comm)
		Rmpi::mpi.send.Robj(H,dest=(r-1) %% cs, tag=5, comm=comm)
		swaps <- swaps + 1
		B <- xB
		H <- xH
	} else {
		Rmpi::mpi.send.Robj(xB,dest=(r-1) %% cs, tag=3, comm=comm)
		Rmpi::mpi.send.Robj(xH,dest=(r-1) %% cs, tag=5, comm=comm)
	}
	return(list(B,LL,H))
}

#' Communicate with other ranks and swap beta
#'
#' Given a current log-likelihood, temperature and step-size, this
#' decides whether to send or receive the same variables from a
#' neighboring process and swap temperatures with them.
#'
#' Only one of the two communicating ranks is allowed to make the
#' decision to swap, because a random variable is used to make this
#' decision. Which rank will make this swap decision needs to be
#' determined somehow. Currently we alternate this reposibility based
#' on the current iteration `i`.
#'
#' @param i MCMC iteration
#' @param B inverse temperature (parallel tempering)
#' @param LL log-likelihood value of current point
#' @param H algorithm's step size (often called epsilon in literature)
#' @param r MPI rank
#' @param comm MPI communicator
#' @param cs MPI comm size
#' @export
pbdMPI_swap_temperatures <- function(i, B, LL, H, r, comm, cs){
	## 1. either send LL, B, and H to neighbor or receive them
	if (r %% 2 == i %% 2) { # alternating, to swap with r-1 or r+1
		pbdMPI::send(LL, rank.dest=(r+1)%%cs, tag=1, comm=comm)
		pbdMPI::send(B, rank.dest=(r+1)%%cs, tag=2, comm=comm)
		pbdMPI::send(H, rank.dest=(r+1)%%cs, tag=4, comm=comm)
	} else { # get e(x)ternal objects, from the other chain
		xLL <- pbdMPI::recv(rank.source=(r-1) %% cs, tag=1, comm=comm)
		xB <- pbdMPI::recv(rank.source=(r-1) %% cs, tag=2, comm=comm)
		xH <- pbdMPI::recv(rank.source=(r-1) %% cs, tag=4, comm=comm)
	}
	## 2. either make a decision on swapping or accept the other rank's decision
	if (r %% 2 == i %% 2){
		B <- pbdMPI::recv(rank.source=(r+1)%%cs, tag=3, comm=comm)
		H <- pbdMPI::recv(rank.source=(r+1)%%cs, tag=5, comm=comm)
	} else if(change_temperature(B,LL,xB,xLL)){
		pbdMPI::send(B,rank.dest=(r-1) %% cs, tag=3, comm=comm)
		pbdMPI::send(H,rank.dest=(r-1) %% cs, tag=5, comm=comm)
		B <- xB
		H <- xH
	} else {
		pbdMPI::send(xB,rank.dest=(r-1) %% cs, tag=3, comm=comm)
		pbdMPI::send(xH,rank.dest=(r-1) %% cs, tag=5, comm=comm)
	}
	return(list(B=B,LL=LL,H=H))
}

#' Broadcast to other ranks and swap temperatures with any of them
#'
#' Using this function, at most two ranks will swap.
#'
#' Given a current log-likelihood, temperature and step-size, this
#' funcion will broadcast a log-likelihood value to all other ranks
#' and they can each decide to swap temperatures with the root
#' process. Root is cycled around all ranks (round-robin).
#'
#' Each other rank is allowed to make the
#' offer to swap. The root process decides which rank to swap with.
#'
#' @param i MCMC iteration
#' @param B inverse temperature (parallel tempering)
#' @param LL log-likelihood value of current point
#' @param H algorithm's step size (often called epsilon in literature)
#' @param r MPI rank
#' @param comm MPI communicator
#' @param cs MPI comm size
#' @export
pbdMPI_bcast_reduce_temperatures <- function(i, B, LL, H, r, comm, cs){
	root <- (i %% cs) # rank of root process
	rootLL <- pbdMPI::bcast(LL,root,comm)
	rootB  <- pbdMPI::bcast(B,root,comm)
	rootH  <- pbdMPI::bcast(H,root,comm)
	if (root != r && change_temperature(B,LL,rootB,rootLL)){ # each _other_ process can propose this
		msg <- r
	} else {
		msg <- -1
	}
	xr <- pbdMPI::allreduce(msg,op="max",comm=comm) # now everyone knows which two will swap
	ret <- list(B=B,LL=LL,H=H)
	if (xr>=0 && xr < cs){ # do the swap
		##message(sprintf("On iteration %i rank %i and rank %i are swapping temperatures.",i,root,xr))
		if (r == root){
			xrB  <- pbdMPI::recv(rank.source = xr, tag = 1, comm = comm)
			xrLL <- pbdMPI::recv(rank.source = xr, tag = 2, comm = comm)
			xrH  <- pbdMPI::recv(rank.source = xr, tag = 3, comm = comm)
			ret <- list(B=xrB,LL=xrLL,H=xrH)
		} else if (r == xr){
			pbdMPI::send(B,rank.dest = root, comm = comm, tag = 1)
			pbdMPI::send(LL,rank.dest = root, comm = comm, tag = 2)
			pbdMPI::send(H,rank.dest = root, comm = comm, tag = 3)
			ret <- list(B=rootB,LL=rootLL,H=rootH)
		}
	}
	return(ret)
}


#' The MPI version of the mcmc function
#'
#' this version of the MCMC function returns a Markov chain closure
#' that assumes that it is bein run in an MPI context: R was launched
#' using `runmpi` and the Rmpi package is installed. The chains shall
#' communicate using the provided `comm`.
#'
#' @param update an update function
#' @param comm an mpi comm which this function will use for send/receive operations
#' @param swapDelay swaps will be attempted every 2*swapDelay+1 iterations
#' @return an mcmc closure m(parMCMC,N,eps) that implicitly uses the supplied update function
#' @export
mcmc_mpi <- function(update, comm, swapDelay=0, swapFunc=rmpi_swap_temperatures){
	D <- max(2*swapDelay+1,1)
	M <- function(parMCMC,N=1000,eps=1e-4){
		r <- attr(comm,"rank") # 0..n-1
		cs  <- attr(comm,"size")
		sample <- matrix(nrow=N,ncol=length(parMCMC))
		ll <- numeric(N)
		b <- numeric(N)
		a <- 0
		swaps <- 0
		for (i in seq(N)){
			parMCMC <- update(parMCMC,eps)
			sample[i,] <- parMCMC
			LL <- attr(parMCMC,"logLikelihood")
			B <- attr(parMCMC,"beta")
			a <- a + as.numeric(attr(parMCMC,"accepted"))
			ll[i] <- LL
			b[i]  <- B
			if (i %% D == 0){ # e.g. i = 3,6,9, or i = 5,10,15
				res <- swapFunc(i,B,LL,eps,r,comm,cs)
				if (res$B != B) swaps <- swaps+1
				attr(parMCMC,"beta") <- res$B
				attr(parMCMC,"logLikelihood") <- res$LL
				eps <- res$H
			}
		}
		attr(sample,"acceptanceRate") <- a/N
		attr(sample,"logLikelihood") <- ll
		attr(sample,"lastPoint") <- parMCMC
		attr(sample,"beta") <- b
		attr(sample,"swapRate") <- swaps/N
		attr(sample,"stepSize") <- eps
		return(sample)
	}
	return(M)
}

#' This function merges mpi-samples into one
#'
#' When using MPI, we save the sample immediately into a file, each
#' rank saves to its own file. This function collects all of these
#' smaller samples into one. The samples should be saved with
#' `saveRDS()`.
#'
#' @export
#' @param files the files where the individual samples are stored
#' @return one matrix where all samples are concatenated.
loadSample_mpi <- function(files){
	s <- lapply(files,readRDS)
	betaTrace <- Reduce(function(a,b) c(a,attr(b,"beta")),s,init=NULL)
	uB <- sort(unique(betaTrace),decreasing=TRUE)
	bSelection <- lapply(uB,function(b) abs(betaTrace-b)<1e-8)
	acc <- Reduce(function(a,b) c(a,attr(b,"acceptanceRate")),s,init=NULL)
	sR <- Reduce(function(a,b) c(a,attr(b,"swapRate")),s,init=NULL)
	ll <- Reduce(function(a,b) c(a,attr(b,"logLikelihood")),s,init=NULL)
	cat("loading sample files with acceptances:\n")
	print(acc)
	Sample <- Reduce(rbind,s)
	return(list(Sample=Sample,beta=betaTrace,acceptanceRate=acc,swapRate=sR,logLikelihood=ll,betaSelection=bSelection))
}

is.invertible <- function(G=NULL,abs_tol=1e-11){
	return(!is.null(G) && is.numeric(G) && is.matrix(G) && all(dim(G)>0) && !any(is.na(G)) && all(is.finite(G)) && isSymmetric(G) && rcond(G) > abs_tol)
}

dmvnorm <- function(x,mean,precision){
	mu <- mean
	stopifnot(!is.null(x) && is.numeric(x) && is.numeric(mu) && length(x) == length(mu))
	stopifnot(!is.null(precision) && is.matrix(precision))
	k <- length(mu)
	xm <- (x-mu)
	xmPxm <- (x-mu) %*% precision %*% (x-mu)
	C <- (2*pi)^(-0.5*k) * sqrt(abs(det(precision))) * exp(-0.5 * xmPxm)
	return(C)
}

rmvnorm <- function(mean,precision){
	mu <- mean
	stopifnot(is.numeric(mu) && is.matrix(precision))
	k <- length(mu)
	P <- chol(0.5*(t(precision) %*% precision))
	x <- solve(P,rnorm(k,0,1))+mu
	return(x)
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
smmala_move <- function(beta,parGiven,fisherInformationPrior,eps=1e-2){
	stopifnot(parGiven %has% c("fisherInformation","gradLogLikelihood","gradLogPrior"))
	stopifnot(!is.null(beta) && !any(is.na(beta)) && is.finite(beta))
	fiGiven <- attr(parGiven,"fisherInformation")
	gradLGiven <- attr(parGiven,"gradLogLikelihood")
	gradPGiven <- attr(parGiven,"gradLogPrior")
	G0 <- fisherInformationPrior
	G <- (beta^2*fiGiven)+G0
	stopifnot(!is.null(G) && is.matrix(G))
	#cat("is.invertible(G): ",is.invertible(G), " (rcond: ",rcond(G),").\n")
	if (isTRUE(is.invertible(G))){
		g <- solve(G,beta*gradLGiven+gradPGiven)
		#Sigma <- solve(G)
	} else {
		stopifnot(is.matrix(G0) && isSymmetric(G0))
		g <- solve(G0,gradPGiven)
		G <- G0
		#Sigma <- solve(G0)
	}
	parProposal <- rmvnorm(
		mean=parGiven+0.5*eps*g,
		precision=G/eps
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
smmala_move_density <- function(beta,parProposal,parGiven,fisherInformationPrior,eps=1e-2){
	stopifnot(parGiven %has% c("fisherInformation","gradLogLikelihood","gradLogPrior"))
	stopifnot(!is.null(beta) && !any(is.na(beta)) && is.finite(beta))
	fiGiven <- attr(parGiven,"fisherInformation")
	gradLGiven <- attr(parGiven,"gradLogLikelihood")
	gradPGiven <- attr(parGiven,"gradLogPrior")
	G0 <- fisherInformationPrior
	G <- (beta^2*fiGiven)+G0
	stopifnot(!is.null(G) && is.matrix(G))
	#cat("is.invertible(G): ",is.invertible(G), " (rcond: ",rcond(G),").\n")
	if (isTRUE(is.invertible(G))){
		g <- solve(G,beta*gradLGiven+gradPGiven)
		#Sigma <- solve(G)
	} else {
		stopifnot(is.matrix(G0) && isSymmetric(G0))
		g <- solve(G0,gradPGiven)
		G <- G0
		#Sigma <- solve(G0)
	}
	return(dmvnorm(as.numeric(parProposal),
		mean=parGiven+0.5*eps*g,
		precision=G/eps)
	)
}


metropolisUpdate <- function(simulate, experiments, model, logLikelihood, dprior){
	U <- function(parGiven, eps=1e-4){
		stopifnot(is.numeric(parGiven) && length(parGiven)>0 && all(is.finite(parGiven)))
		stopifnot(parGiven %has% c("logLikelihood","prior"))
		beta <- attr(parGiven,"beta") %otherwise% 1.0
		llGiven <- attr(parGiven,"logLikelihood")
		priorGiven <- attr(parGiven,"prior")
		stopifnot(is.numeric(llGiven) && length(llGiven)==1 && is.finite(llGiven))
		stopifnot(is.numeric(priorGiven) && length(priorGiven)==1 && is.finite(priorGiven))
		stopifnot(!is.null(eps) && is.numeric(eps) && length(eps)==1 && is.finite(eps))
		parProposal <- parGiven + rnorm(length(parGiven),0,eps)
		stopifnot(is.numeric(parProposal) && all(is.finite(priorGiven)))
		attr(parProposal,"simulations") <- simulate(parProposal)
		llProposal <- logLikelihood(parProposal)
		stopifnot(is.numeric(llProposal) && length(llProposal)==1 && is.finite(llProposal))
		attr(parProposal,"logLikelihood") <- llProposal
		priorProposal <- dprior(parProposal)
		attr(parProposal,"prior") <- priorProposal
		L <- exp(beta*(llProposal - llGiven))
		if (is.null(L)) cat("llProposal: ",llProposal," llGiven: ",llGiven," beta: ",beta," L:",L,"\n")
		P <- priorProposal/priorGiven
		if (!is.null(L) && is.finite(L) && runif(1) < L*P){
			attr(parProposal,"accepted") <- TRUE
			return(parProposal)
		} else {
			attr(parGiven,"accepted") <- FALSE
			return(parGiven)
		}
	}
	attr(U,"algorithm") <- "Metropolis-Hastings"
	return(U)
}

smmalaUpdate <- function(simulate, experiments, model, logLikelihood, dprior, gradLogLikelihood, gprior, fisherInformation, fisherInformationPrior){
	#np <- length(model$par())
	#nu <- length(experiments[[1]]$input)
	U <- function(parGiven, eps=1e-4){
		stopifnot(parGiven %has% c("logLikelihood","prior","fisherInformation","gradLogLikelihood","gradLogPrior"))
		fp <- fisherInformationPrior
		stopifnot(is.matrix(fp) && isSymmetric(fp))
		beta <- attr(parGiven,"beta") %otherwise% 1.0
		llGiven <- attr(parGiven,"logLikelihood")
		priorGiven <- attr(parGiven,"prior")
		n <- length(parGiven)
		## the very important step: suggest a successor to parGiven and simulate the model
		parProposal <- smmala_move(beta,parGiven,fp,eps)
		#if (any(is.na(parProposal))) {
		#	attr(parGiven,"accepted") <- FALSE
		#	return(parGiven)
		#}
		##cat("rank ",r,", parProposal: ",as.numeric(parProposal),"\n")
		flush.console()
		attr(parProposal,"simulations") <- simulate(parProposal)
		llProposal <- logLikelihood(parProposal)
		priorProposal <- dprior(parProposal)
		attr(parProposal,"beta") <- beta
		attr(parProposal,"logLikelihood") <- llProposal
		attr(parProposal,"prior") <- priorProposal
		attr(parProposal,"fisherInformation") <- fisherInformation(parProposal)
		##cat("rank ",r," rcond(fisherInformation): ", rcond(attr(parProposal,"fisherInformation")),".\n")
		attr(parProposal,"gradLogLikelihood") <- gradLogLikelihood(parProposal)
		##cat("rank",r,"gradient-LL: ", attr(parProposal,"gradLogLikelihood"),"\n")
		attr(parProposal,"gradLogPrior") <- gprior(parProposal)
		##cat("rank",r,"gradient-PR: ", attr(parProposal,"gradLogPrior"),"\n")
		fwdDensity <- smmala_move_density(beta,parProposal,parGiven,fp,eps)
		bwdDensity <- smmala_move_density(beta,parGiven,parProposal,fp,eps)
		L <- exp(beta*(llProposal - llGiven))
		P <- priorProposal/priorGiven
		K <- bwdDensity/fwdDensity
		##cat("rank: ",r,"; beta: ",beta, "; llProposal: ",llProposal,"; llGiven: ",llGiven,".\n")
		##cat("rank: ",r,";L ",L,";P: ",P,";K: ",K,".\n")
		flush.console()
		if (is.finite(L) && is.finite(K) && runif(1) < L * P * K){
			attr(parProposal,"accepted") <- TRUE
			return(parProposal)
		} else {
			attr(parGiven,"accepted") <- FALSE
			return(parGiven)
		}
	}
	attr(U,"algorithm") <- "smmala"
	return(U)
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
#' @param simulate a function that simulates the model for a given
#'     parMCMC
#' @param experiments the list of experiments (with simulation
#'     instructions)
#' @param model the list of model functions
#' @param logLikelihood a function that calculates log-likelihood
#'     values for given parMCMC
#' @param gradLogLikelihood a function that calculates the gradient of
#'     the log-likelihood for given parMCMC
#' @param fisherInformation a function that calculates approximates
#'     Fisher information matrices
#' @param fisherInformationPrior a constant matrix, the prior
#'     distributions fisher information
#' @param dprior prior density function
#' @param gprior gradient of the prior density
#' @return a function that returns possibly updated states of the
#'     Markov chain
mcmcUpdate <- function(simulate, experiments, model, logLikelihood, dprior, gradLogLikelihood=NULL, gprior=NULL, fisherInformation=NULL, fisherInformationPrior=NULL){
	if (is.null(gradLogLikelihood)){ # Metropolis Hastings
		return(metropolisUpdate(simulate, experiments, model, logLikelihood, dprior))
	} else {
		return(smmalaUpdate(simulate, experiments, model, logLikelihood, dprior, gradLogLikelihood, gprior, fisherInformation, fisherInformationPrior))
	}
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
fisherInformationFunc <- function(model, experiments, parMap=identity, parMapJac=function (x) {diag(1,length(x))}){
	nF <- length(model$func(0.0,model$init(),model$par()))
	l10 <- log(10)
	F <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		np <- length(parMCMC)
		fi  <- matrix(0.0,np,np)
		for (i in seq(length(experiments))){
			errF <- t(experiments[[i]]$errorValues)
			if (any(is.na(errF))){
				errF[is.na(errF)] <- Inf
			}
			nt <- length(experiments[[i]]$outputTimes)
			for (j in seq(nt)){
				sigma_j <- matrix(errF[,j],nF,np)
				if (any(is.na(sigma_j))){
					message(sprintf("sigma_j has invalid elements"));
				}
				# output function sensitivity:
				Sh <- simulations[[i]]$funcSensitivity[[1]][,seq(np),j]
				lNA <- is.na(Sh)
				if (any(lNA)) {
					Sh[lNA] <- 0.0
					##message(sprintf("Sensitivity has %i missing values. Replacing with 0.",sum(lNA)))
				}
				dim(Sh) <- c(nF,np)
				pmj <- parMapJac(as.numeric(parMCMC))
				Sh <- (Sh %*% pmj)/sigma_j
				if (any(is.finite(Sh))){
					fi  <-  fi + t(Sh) %*% Sh
				}
			}
		}
		lNAfi <- is.na(fi)
		if (any(lNAfi)) {
			##message(sprintf("Fisher-information has %i missing values. Setting them to 0.0.\n",sum(lNAfi)))
			fi[!is.finite(fi)] <- 0.0
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
logLikelihoodFunc <- function(experiments){
	N <- length(experiments)
	llf <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		dimFunc <- dim(simulations[[1]]$func)
		n <- dimFunc[3]
		m <- head(dimFunc,2)
		L <- rep(-0.5*prod(m)*N*log(2*pi),n)
		for (i in seq(N)){
			y <- t(experiments[[i]]$outputValues)
			stdv <- t(experiments[[i]]$errorValues)
			for (k in seq(n)){
				h <- simulations[[i]]$func[,,k]
				dim(h) <- m
				stopifnot(all(dim(h)==dim(y)) && all(dim(y)==dim(stdv)))
				L[k] <- L[k] - 0.5*sum(((y - h)/stdv)^2,na.rm=TRUE) - sum(log(stdv),na.rm=TRUE)
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
gradLogLikelihoodFunc <- function(model,experiments,parMap=identity,parMapJac=function(x) {diag(1,length(x))})
{
	N <- length(experiments)
	nF <- length(model$func(0.0,model$init(),model$par()))
	gradLL <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		np <- length(parMCMC) # the dimension of the MCMC variable (parMCMC)
		z <- rep(0,np)
		gL <- z
		for (i in seq(N)){
			d <- dim(simulations[[i]]$func)
			nt <- length(experiments[[i]]$outputTimes)
			y <- t(experiments[[i]]$outputValues)
			y[is.na(y)] <- 0.0
			h <- simulations[[i]]$func[,,1]
			dim(h) <- head(d,2)
			stdv2 <- t(experiments[[i]]$errorValues)^2 # sigma^2
			stdv2[is.na(stdv2)] <- Inf
			Sh <- simulations[[i]]$funcSensitivity[[1]]
			for (j in seq(nt)){
				Shj <- Sh[,seq(np),j]
				dim(Shj) <- c(nF,np)
				Shj <- Shj %*% parMapJac(as.numeric(parMCMC))
				gj <- as.numeric(t((y[,j] - h[,j])/stdv2[,j]) %*% Shj)
				if (all(is.finite(gj))) {
					gL <- gL + gj
				}
			}
		}
		return(gL)
	}
	return(gradLL)
}
