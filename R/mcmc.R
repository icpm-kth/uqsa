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

#' Determine a prefix from a character vector str of similar contents
#'
#' The result is such that `all(startsWith(str,determinePrefix(str)))` is `TRUE`.
#'
#' By default, the strings are assumed to be '-' separated words, and
#' a series of words is found to be thje prefix if all entries start
#' with that set of words.
#'
#' The normal case is c("abc-1","abc-2b","abc-2a") maps to "abc"
#' @export
#' @param str a character vector
#' @param split the token to use for strsplit instead of '-', this should be `character(0)` if you want to split letter by letter
#' @return the prefix common to all entries of str.
determinePrefix <- function(str,split="-",collapse="-"){
return(paste(
	Reduce(
		function(a,b) {
			m<-seq(min(length(a),length(b)));
			a<-a[m]; b<-b[m];
			i<-which(unlist(mapply(identical,a,b)));
			return(a[i])
		},
		strsplit(str,split)),
		collapse=collapse)
	)
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
mcmcInit <- function(beta,parMCMC,simulate,logLikelihood=ll,dprior=\(x) prod(rnorm(x)),gradLogLikelihood=NULL,gprior=NULL,fisherInformation=NULL){
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
#' @param ll1 the log-likelihood of chain 1
#' @param b2 the inverse temperature of chain 2
#' @param ll2 the log-likelihood of chain 2
#' @return TRUE is the chains should swap their temperatures
#' @export
change_temperature <- function(b1,ll1,b2,ll2){
	a <- exp((b2-b1)*(ll1-ll2))
	r <- runif(1)
	return(r<a)
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
		colnames(sample) <- names(parMCMC)
		ll <- numeric(N)
		b <- numeric(N)
		a <- logical(N)
		for (i in seq(N)){
			parMCMC <- update(parMCMC,eps)
			ll[[i]] <- attr(parMCMC,"logLikelihood")
			sample[i,] <- as.numeric(parMCMC)
			b[i] <- attr(parMCMC,"beta")
			a[i] <- attr(parMCMC,"accepted")
		}
		attr(sample,"acceptanceRate") <- mean(a)
		attr(sample,"acceptance") <- a
		attr(sample,"logLikelihood") <- ll
		attr(sample,"lastPoint") <- parMCMC
		attr(sample,"beta") <- b
		attr(sample,"stepSize") <- eps
		return(sample)
	}
	comment(M) <- "function(parMCMC,N=1000,eps=1e-4) where eps is the step-size (numeric scalar)"
	return(M)
}
#' Broadcast to other ranks and swap temperatures with any of them
#'
#' Using this function, at most two ranks will swap.
#'
#' Given a current log-likelihood, temperature and step-size, this
#' function will broadcast a log-likelihood value to all other ranks
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
#' @param swapDelay swaps will be attempted every 2*swapDelay+1 iterations [deprecated]
#' @return an mcmc closure m(parMCMC,N,eps) that implicitly uses the supplied update function
#' @export
mcmc_mpi <- function(update, comm, swapDelay=0, swapFunc=pbdMPI_bcast_reduce_temperatures){
	D <- max(2*swapDelay+1,1)
	M <- function(parMCMC,N=1000,eps=1e-4){
		r <- attr(comm,"rank") # 0..n-1
		cs  <- attr(comm,"size")
		sample <- matrix(nrow=N,ncol=length(parMCMC))
		colnames(sample) <- names(parMCMC)
		ll <- numeric(N)
		b <- numeric(N)
		a <- logical(N)
		swaps <- logical(N)
		h <- numeric(N)
		for (i in seq(N)){
			parMCMC <- update(parMCMC,eps)
			sample[i,] <- as.numeric(parMCMC)
			LL <- attr(parMCMC,"logLikelihood")
			B <- attr(parMCMC,"beta")
			a[i] <- attr(parMCMC,"accepted")
			ll[i] <- LL
			b[i]  <- B
			h[i] <- eps
			for (j in seq(0,cs-1)){
				res <- swapFunc(j,B,LL,eps,r,comm,cs)
				B <- res$B
				LL <- res$LL
				eps <- res$H
			}
			attr(parMCMC,"beta") <- B
			attr(parMCMC,"logLikelihood") <- LL
			swaps[i] <- (B != b[i])
		}
		attr(sample,"acceptanceRate") <- mean(a)
		attr(sample,"acceptance") <- a
		attr(sample,"logLikelihood") <- ll
		attr(sample,"lastPoint") <- parMCMC
		attr(sample,"beta") <- b
		attr(sample,"swapRate") <- mean(swaps)
		attr(sample,"swaps") <- swaps
		attr(sample,"stepSize") <- h
		return(sample)
	}
	comment(M) <- "function(parMCMC, N=1000, eps=1e-4) where eps is the step size, and N the sample size"
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
	bSelection <- lapply(uB,function(b) abs(betaTrace-b)<1e-10)
	acc <- Reduce(function(a,b) c(a,attr(b,"acceptanceRate")),s,init=NULL)
	sR <- Reduce(function(a,b) c(a,attr(b,"swapRate")),s,init=NULL)
	ll <- Reduce(function(a,b) c(a,attr(b,"logLikelihood")),s,init=NULL)
	cat("loading sample files with acceptances:\n")
	print(acc)
	Sample <- Reduce(rbind,s)
	return(list(Sample=Sample,beta=betaTrace,acceptanceRate=acc,swapRate=sR,logLikelihood=ll,betaSelection=bSelection,uB=uB))
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
#' @param size sub sample size, if not set, the whole sample is
#'     returned
#' @param selection integer index vector or logical vector indicating
#'     which temperatures to return: beta\[selection\] is returned, in
#'     decreasing order of beta.
#' @param mc.cores defaults to the total number of cores, but can be
#'     reduced with this option.
#' @return a list of matrices, by temperature, concatenated.
loadSubSample_mpi <- function(files,size=NA,selection=NA,mc.cores=parallel::detectCores()){
	S <- parallel::mclapply(files,function(f){
		s <- readRDS(f)
		b <- attr(s,"beta")
		cat(sprintf("acceptance rate: %02f\n",attr(s,"acceptanceRate")))
		if (!any(is.na(size))){
			j <- seq(1,NROW(s),length.out=size)
			s <- s[j,]
			b <- b[j]
		}
		attr(s,"beta") <- b
		return (s)
	},mc.cores=mc.cores)
	uB <- sort(unique(unlist(parallel::mclapply(S,function(s) {return(unique(attr(s,"beta")))}))),decreasing=TRUE)
	cat("unique temperatures:",uB,"\n")
	if (length(uB)!=length(files)) {
		warning(sprintf("number of temperatures (%i) not the same as number of files (%i).",length(uB),length(files)))
	}
	if (!any(is.na(selection))) uB <- uB[selection]
	else selection <- seq_along(uB)
	if (is.list(S) && length(S) == length(files) && is.numeric(S[[1]]) && is.matrix(S[[1]])){
		n <- NROW(S[[1]])
		m <- NCOL(S[[1]])
		l <- length(uB)
	} else {
		print(S)
		error("loading samples did not succeed")
	}
	x <- array(NA,dim=c(n,m,l))
	dimnames(x) <- list(NULL,colnames(S[[1]]),sprintf("beta_%02i",selection))
	for (i in seq_along(S)){
		for (k in seq_along(uB)){
			j <- which(attr(S[[i]],"beta") == uB[k])
			x[j,,k] <- S[[i]][j,]
		}
	}
	attr(x,"beta") <- uB
	return(x)
}

#' gatherSample collects all sample points, from all files, with the
#' given temperature
#'
#' This function assumes that each supplied RDS file contains a matrix
#' of model MCMC parameters, with an attribute called "beta" that
#' lists the temperature of each row.
#'
#' This function selects and collects all rows, from all files with
#' the same (given) temperature.
#' @export
#' @param files a list of file names
#' @param beta the inverse temperture to extract sample for
#' @param size a size the is smaller than the actual sample size, if
#'     left unchanged, all sampled points are returned
#' @return a matrix of sampled points, all with the same temperature
gatherSample <- function(files,beta=1.0,size=NA){
	x <- numeric(0)
	lL <- numeric(0)
	H <- numeric(0)
	for (f in files){
		s <- readRDS(f)
		b <- attr(s,"beta")
		l <- attr(s,"logLikelihood")
		h <- attr(s,"stepSize")
		i <- which(abs(b - beta) <= 1e-9*beta + 1e-15)
		s <- s[i,,drop=FALSE]
		l <- l[i]
		if (!is.null(h) && length(h) == length(b)) {
			h <- h[i]
		}
		x <- rbind(x,s)
		lL <- c(lL,l)
		H <- c(H,h)
		if (!any(is.na(size)) && NROW(x)>=size) break
	}
	attr(x,"beta") <- beta
	attr(x,"logLikelihood") <- lL
	attr(x,"stepSize") <- H
	return(x)
}

#' gatherReplicas collects all sample points, from all files, which
#' are assumed to be exact replicas, with different seeds (and
#' possibly sizes). This function uses mclapply to process the files,
#' which may be quicker than gatherSample.  The temperature is
#' disregarded, assuming that no parallel tempering was used.  To
#' facilitate the loading of a very big sample, this function will
#' analyse the auto-correltation within each file and returned a
#' thinned subsample of size N/2*tauint (effective sample size). There
#' is no need to further reduce the rfesult.
#'
#' For small samples, it is better to load the entire sample and
#' analyse it in full. This function is intended for samples that are
#' so big that they challenge the memory of the machine.
#'
#' This function is quicker if you have used trivial parallelism,
#' without mpi communication between the ranks (or another method of
#' obtaining several replicas, like forking or sequetial reprtition).
#'
#' This function assumes that each supplied RDS file contains a matrix
#' of model MCMC parameters. The returned value X will be the rbind of
#' all the smaller x contained in the individual files. The value X
#' will have several attributes attached to it:
#'
#' - logLikelihood: log(likelihood(X\[i,\])), one value per row of X
#' - stepSize: the MCMC step size used in each given file#'
#'
#' @export
#' @param files a list of file names
#' @return a matrix of sampled points, all with the same temperature
gatherReplicas <- function(files){
	if (requireNamespace("hadron") && requireNamespace("errors")){
		tau <- lapply(
			lapply(
				parallel::mclapply(files,readRDS),
				attr,"logLikelihood"
			),
			function(l) {
				res <- hadron::uwerr(data=l)
				tau <- errors::set_errors(res$tauint,res$dtauint)
				return(tau)
			}
		)
	} else {
		tau <- lapply(
			lapply(
				parallel::mclapply(files,readRDS),
				attr,"logLikelihood"
			),
			function(l) {
				res <- acf(l,plot=FALSE,lag.max=length(l))
				tau <- sum(res$acf)
				return(tau)
			}
		)
	}
	x <- Reduce(
		rbind,
		parallel::mcmapply(
			function(x,tau){
				return(x[seq(1,NROW(x),by=ceiling(tau)),])
			},
			parallel::mclapply(files,readRDS),
			tau
		),
		init=numeric(0)
	)
	l <- unlist(
		parallel::mcmapply(
			function(l,tau){
				return(l[seq(1,length(l),by=ceiling(tau))])
			},
			lapply(parallel::mclapply(files,readRDS),attr,"logLikelihood"),
			tau
		)
	)
	h <- sapply(parallel::mclapply(files,readRDS),attr,"stepSize")
	attr(x,"logLikelihood") <- l
	attr(x,"stepSize") <- h
	attr(x,"tau") <- as.data.frame(Reduce(c,tau))
	return(x)
}


#' checks whether a given matrix is a valid, invertible fisherInformation
#'
#' This matrix has to be symmetric and invertible. But, because the
#' matrix has a perhaps sketchy origin, it could be defective in all
#' possible ways.
#'
#' @param G a matrix
#' @param abs_tol absolute tolerance for the reciprocal condition number of G
#' @return TRUE or FALSE
is.invertible <- function(G=NULL,abs_tol=1e-11){
	return(!is.null(G) && is.numeric(G) && is.matrix(G) && all(dim(G)>0) && !any(is.na(G)) && all(is.finite(G)) && isSymmetric(G) && rcond(G) > abs_tol)
}

dmvnorm <- function(x,mean,precision){
	mu <- mean
	stopifnot(!is.null(x) && is.numeric(x) && is.numeric(mu) && length(x) == length(mu))
	stopifnot(!is.null(precision) && is.matrix(precision))
	k <- length(mu)
	stopifnot(dim(precision) == c(k,k))
	xm <- (x-mu)
	xmPxm <- (x-mu) %*% (precision %*% (x-mu))
	C <- (2*pi)^(-0.5*k) * sqrt(abs(det(precision))) * exp(-0.5 * xmPxm)
	return(C)
}

rmvnorm <- function(mean,precision){
	mu <- mean
	stopifnot(!is.null(precision) && !is.null(mu))
	stopifnot(is.numeric(mu) && is.matrix(precision))
	stopifnot(all(is.finite(mu)))
	stopifnot(all(is.finite(precision)))
	k <- length(mu)
	P <- chol(precision) # chol(0.5*(t(precision) %*% precision))
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
	fiGiven <- parGiven %@% "fisherInformation"
	gradLGiven <- parGiven %@% "gradLogLikelihood"
	gradPGiven <- parGiven %@% "gradLogPrior"
	G0 <- fisherInformationPrior
	G <- (beta^2*fiGiven)+G0
	stopifnot(!is.null(G) && is.matrix(G))
	#cat("is.invertible(G): ",is.invertible(G), " (rcond: ",rcond(G),").\n")
	if (isTRUE(is.invertible(G))){
		g <- solve(G,as.numeric(beta*gradLGiven+gradPGiven))
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
	names(parProposal) <- names(parGiven)
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
	fiGiven <- parGiven %@% "fisherInformation"
	gradLGiven <- parGiven %@% "gradLogLikelihood"
	gradPGiven <- parGiven %@% "gradLogPrior"
	G0 <- fisherInformationPrior
	G <- (beta^2*fiGiven)+G0
	stopifnot(!is.null(G) && is.matrix(G))
	if (isTRUE(is.invertible(G))){ # isTRUE has some additional safety for breakage cases we didn't consider
		g <- solve(G,beta*gradLGiven+gradPGiven)
	} else {
		stopifnot(is.matrix(G0) && isSymmetric(G0))
		g <- solve(G0,gradPGiven)
		G <- G0
	}
	return(dmvnorm(as.numeric(parProposal),
		mean=as.numeric(parGiven+0.5*eps*g),
		precision=G/eps)
	)
}

metropolisUpdate <- function(simulate, experiments, logLikelihood, dprior, Sigma=NULL, parAcceptable=\(p) {TRUE}){
	if (!is.null(Sigma)){
		cSigma <- chol(Sigma)
	}
	if (is.null(Sigma)){
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
			names(parProposal) <- names(parGiven)
			priorProposal <- dprior(parProposal)
			attr(parProposal,"prior") <- priorProposal
			if (priorProposal == 0.0 || !parAcceptable(parProposal)) {
				attr(parGiven,"accepted") <- FALSE
				return(parGiven)
			}
			attr(parProposal,"simulations") <- simulate(parProposal)
			llProposal <- logLikelihood(parProposal)
			if (!is.numeric(llProposal) || length(llProposal)!=1){
				warning(sprintf("metropolisUpdate encountered an invalid likelihood value: %f\n",llProposal[1]))
				print(llProposal)
				print(as.numeric(parGiven))
			}
			attr(parProposal,"logLikelihood") <- llProposal
			L <- exp(beta*(llProposal - llGiven))
			if (is.null(L)) cat("llProposal: ",llProposal," llGiven: ",llGiven," beta: ",beta," L:",L,"\n")
			P <- priorProposal/priorGiven
			if (!is.null(L) && !is.null(P) && is.finite(L) && is.finite(P) && runif(1) < L*P){
				attr(parProposal,"accepted") <- TRUE
				return(parProposal)
			} else {
				attr(parGiven,"accepted") <- FALSE
				return(parGiven)
			}
		}
	} else {
		U <- function(parGiven,eps=1e-4){
			stopifnot(is.numeric(parGiven) && length(parGiven)>0 && all(is.finite(parGiven)))
			stopifnot(parGiven %has% c("logLikelihood","prior"))
			beta <- attr(parGiven,"beta") %otherwise% 1.0
			llGiven <- attr(parGiven,"logLikelihood")
			priorGiven <- attr(parGiven,"prior")
			parProposal <- parGiven + solve(cSigma,rnorm(length(parGiven),0,eps))
			stopifnot(is.numeric(parProposal) && all(is.finite(priorGiven)))
			names(parProposal) <- names(parGiven)
			priorProposal <- dprior(parProposal)
			attr(parProposal,"prior") <- priorProposal
			if (priorProposal < 1e-15*priorGiven+1e-15 || !parAcceptable(parProposal)) {
				attr(parGiven,"accepted") <- FALSE
				return(parGiven)
			}
			attr(parProposal,"simulations") <- simulate(parProposal)
			llProposal <- logLikelihood(parProposal)
			if (!is.numeric(llProposal) || length(llProposal)!=1){
				warning(sprintf("metropolis update encountered an invalid likelihood value: %f\n",llProposal[1]))
				print(llProposal)
				print(as.numeric(parGiven))
			}
			attr(parProposal,"logLikelihood") <- llProposal
			L <- exp(beta*(llProposal - llGiven))
			if (is.null(L)) cat("llProposal: ",llProposal," llGiven: ",llGiven," beta: ",beta," L:",L,"\n")
			P <- priorProposal/priorGiven
			if (!is.null(L) && !is.null(P) && is.finite(L) && is.finite(P) && runif(1) < L*P){
				attr(parProposal,"accepted") <- TRUE
				return(parProposal)
			} else {
				attr(parGiven,"accepted") <- FALSE
				return(parGiven)
			}
		}
	}
	attr(U,"algorithm") <- "Metropolis-Hastings"
	return(U)
}

#' Default Log-likelihood Function
#'
#' Extracts the `logLikelihood` value from the simulations attribute
#' of the parMCMC argument, requires:
#' - parMCMC has simulations attribute
#' - simulations list includes logLikelihood values
#'
#' This function will take the ll-values claculated by the ode solver
#' in this package, and return the sum of those values over all
#' experiments. The ll-value the simulator returns is calculated with
#' the assumption of a normal distribution on measurement errors.
#'
#' @param parMCMC a numeric vector, with attributes for MCMC, specifically smmala
#' @return a scalar value: log(likelihood(data|parMCMC))
#' @export
ll <- function(parMCMC){
	y <- parMCMC %@% "simulations"
	return(sum(sapply(y,\(sim) sim$logLikelihood)))
}

#' Default Fisher Information Function
#'
#' Extracts the `FisherInformation` values from the simulations attribute
#' of the parMCMC argument, requires:
#' - parMCMC has simulations attribute
#' - simulations list includes Fisher-Information values
#'
#' This function will take the fi-values claculated by the ode solver
#' in this package, and return the sum of those values over all
#' experiments. The gll-value the simulator returns is calculated with
#' the assumption of a normal distribution on measurement errors.
#'
#' This function uses the log10ParMapJac(parMCMC) function by default,
#' which assumes that sampling takes place in logarothmic space.
#'
#' @param parMCMC a numeric vector, with attributes for MCMC, specifically smmala
#' @return a numeric vector: grad(log(likelihood(data|parMCMC)))
#' @export
gllf <- function(parMapJac=log10ParMapJac) {
	return(
		function(parMCMC){
			g <- Reduce(
				\(a,b) a+as.numeric(b$gradLogLikelihood[seq_along(parMCMC)]),
				parMCMC %@% "simulations",
				init=0.0
			)
			return(
				as.numeric(g %*% parMapJac(parMCMC))
			)
		}
	)
}

#' Default gradient-Log-likelihood Function
#'
#' Extracts the `gradLogLikelihood` values from the simulations attribute
#' of the parMCMC argument, requires:
#' - parMCMC has simulations attribute
#' - simulations list includes gradLogLikelihood values
#'
#' This function will take the gll-values claculated by the ode solver
#' in this package, and return the sum of those values over all
#' experiments. The gll-value the simulator returns is calculated with
#' the assumption of a normal distribution on measurement errors.
#'
#' This function uses the log10ParMapJac(parMCMC)
#'
#' @param parMCMC a numeric vector, with attributes for MCMC, specifically smmala
#' @return a scalar value: log(likelihood(data|parMCMC))
#' @export
fi <- function(parMapJac=log10ParMapJac){
	return(
		function(parMCMC) {
			i <- seq_along(parMCMC)
			f <- Reduce(
				\(a,b) a + b$FisherInformation[i,i,1],
				parMCMC %@% "simulations",
				init = 0.0
			)
			return(
				t(parMapJac(parMCMC)) %*% f %*% parMapJac(parMCMC)
			)
		}
	)
}

smmalaUpdate <- function(simulate, experiments, logLikelihood=ll, dprior=\(x) prod(dnorm(x)), gradLogLikelihood=gllf(log10ParMapJac), gprior=\(x) (-x), fisherInformation=fi, fisherInformationPrior=0, parAcceptable=\(p) all(is.finite(p))){
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
		priorProposal <- dprior(parProposal)
		if (priorProposal < 1e-15*priorGiven+1e-15 || !parAcceptable(parProposal)){
			attr(parGiven,"accepted") <- FALSE
			return(parGiven)
		}
		attr(parProposal,"simulations") <- simulate(parProposal)
		llProposal <- logLikelihood(parProposal)
		if (!is.numeric(llProposal) || length(llProposal)!=1){
			warning(sprintf("smmala update encountered an invalid likelihood value: %f\n",llProposal[1]))
			print(llProposal)
			print(as.numeric(parGiven))
		}
		attr(parProposal,"beta") <- beta
		attr(parProposal,"logLikelihood") <- llProposal
		attr(parProposal,"prior") <- priorProposal
		attr(parProposal,"fisherInformation") <- fisherInformation(parProposal)
		attr(parProposal,"gradLogLikelihood") <- gradLogLikelihood(parProposal)
		attr(parProposal,"gradLogPrior") <- gprior(parProposal)
		fwdDensity <- smmala_move_density(beta,parProposal,parGiven,fp,eps)
		bwdDensity <- smmala_move_density(beta,parGiven,parProposal,fp,eps)
		L <- exp(beta*(llProposal - llGiven)) %otherwise% 0.0
		P <- priorProposal/priorGiven %otherwise% 0.0
		K <- bwdDensity/fwdDensity %otherwise% 0.0
		flush.console()
		if (is.finite(L) && is.finite(K) && is.finite(P) && runif(1) < L * P * K){
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
#' @param logLikelihood a function that calculates log-likelihood
#'     values for given parMCMC
#' @param dprior prior density function
#' @param gradLogLikelihood a function that calculates the gradient of
#'     the log-likelihood for given parMCMC
#' @param gprior gradient of the prior density
#' @param fisherInformation a function that calculates approximates
#'     Fisher information matrices
#' @param fisherInformationPrior a constant matrix, the prior
#'     distributions fisher information
#' @param Sigma alternatively, `Sigma=solve(fisherInformationPrior)`,
#'     \[the inverse to fisherInformationPrior\] can be specified for
#'     the metropolis algorithm
#' @param parAcceptable a user shaped function that returns a Boolean
#'     scalar indicating whether a proposed parameter satisfies any
#'     user chosen constraints. A value of FALSE will trigger an early
#'     return and the model will not be simulated.
#' @return a function that returns possibly updated states of the
#'     Markov chain
mcmcUpdate <- function(simulate, experiments, logLikelihood, dprior, gradLogLikelihood=NULL, gprior=NULL, fisherInformation=NULL, fisherInformationPrior=NULL, Sigma=NULL, parAcceptable=\(p) {TRUE}){
	if (is.null(Sigma) && !is.null(fisherInformationPrior)) {
		Sigma <- solve(fisherInformationPrior)
	}
	if (is.null(gradLogLikelihood)) { # Metropolis Hastings
		return(metropolisUpdate(simulate, experiments, logLikelihood, dprior, Sigma=Sigma, parAcceptable=parAcceptable))
	} else {
		return(smmalaUpdate(simulate, experiments, logLikelihood, dprior, gradLogLikelihood, gprior, fisherInformation, fisherInformationPrior, parAcceptable))
	}
}

#' Default log-likelihood function
#'
#' This returns a function f(simulations), which maps simulation
#' results to log(likelihood) values. The experiments are used
#' implicitly; simulations is a list as returned by
#' rgsl::r_gsl_odeiv2_outer().
#'
#' @param experiment will be compared tp the simulation results
#' @param perExpLLF (optional) a user supplied function with the
#'     interface `perExpLLF(p,s,e)`, where `p` are the parameters, `s`
#'     are the simulations and `e` are the experiments (with
#'     data). Supply this function if some of your experiments need to
#'     be normalized by the other experiments (and other complex
#'     cases).
#' @param simpleUserLLF (optional) a user supplied function that is
#'     used instead of the default sum of ((y-h)/stdv)^2 terms. The
#'     interface is: `simpleUserLLF(y,h,stdv,name=NULL)`, where each
#'     of them is an N-M-matrix where N is the dimensionality of the
#'     model output and M the number of data time-points.  Here, `y`
#'     is `t(experiments[[i]]$data)` and may contain NA
#'     values.  This function should also accept an optional _name_
#'     argument (this is the name of the experiment this function is
#'     currently called for).
#' @return `llf(parMCMC)`, a closure (function) of the mcmc-variable:
#'     parMCMC; returns a scalar log-likelihood value. Alternatively,
#'     the user can define such a function:
#'     `parMCMC -> log(Likelihood(parMCMC))`,
#'     and use that during sampling. A test simulation of `p`:
#'     `y <- simulate(p)` will reveal which values the simulator produces.
#'     These values will be attached to p during sampling, as an
#'     attribute.  mcmcInit will attach the same values for the
#'     initial Markov chain state.  The log-likelihood function can
#'     use these attributes.
#' @export
logLikelihoodFunc <- function(experiments,perExpLLF=NULL,simpleUserLLF=NULL){
	N <- length(experiments)
	n.out <- sum(unlist(lapply(experiments,\(e) sum(!is.na(e$data))))) # total number of valid values
	message(sprintf("experiments contain %i non-missing values",n.out))
	if (!is.null(simpleUserLLF)){
		llf <- function(parMCMC){
			if (!("simulations" %in% names(attributes(parMCMC)))) {
				return(-Inf)
			} else if (any(is.na(attr(parMCMC,"simulations")))) {
				return(-Inf)
			} else {
				simulations <- attr(parMCMC,"simulations")
			}
			n <- NCOL(parMCMC)
			L <- rep(0,n)
			for (i in seq(N)){
				if (!("func" %in% names(simulations[[i]]))) {
					warning("simulations contain no 'func' values (output function values, i.e. observables).")
					return(-Inf)
				} else if (any(is.na(simulations[[i]]$func))) {
					warning("simultions contain NA values which shouldn't happen, normally.")
					return(-Inf)
				}
				dimFunc <- dim(simulations[[i]]$func)
				m <- head(dimFunc,2)
				y <- experiments[[i]]$data %otherwise% t(experiments[[i]]$measurements)
				stdv <- standard_error_matrix(y) %otherwise% t(experiments[[i]]$standardError)
				for (k in seq(n)){
					h <- simulations[[i]]$func[,,k]
					dim(h) <- m
					stopifnot(all(dim(h)==dim(y)) && all(dim(y)==dim(stdv)))
					L[k] <- L[k] + simpleUserLLF(y,h,stdv,names(experiments)[i])
				}
			}
			return(L)
		}
	} else if (!is.null(perExpLLF)){
		llf <- function(parMCMC){
			if (!("simulations" %in% names(attributes(parMCMC))) || any(is.na(attr(parMCMC,"simulations")))) {
				warning("no simulations attached to parameter vector: attr(parMCMC,'simulations') is missing.")
				return(-Inf)
			} else {
				simulations <- attr(parMCMC,"simulations")
			}
			simulations <- attr(parMCMC,"simulations")
			L <- perExpLLF(parMCMC,simulations,experiments)
			return(L)
		}
	} else {
		llf <- function(parMCMC){
			if (!(parMCMC %has% "simulations")) {
				warning("no simulations attached to parameter vector: attr(parMCMC,'simulations') is missing.")
				return(-Inf)
			} else if (any(is.na(attr(parMCMC,"simulations")))) {
				warning("simulated values contain NA, which shouldn't happen normally.")
				return(-Inf)
			} else {
				simulations <- attr(parMCMC,"simulations")
			}
			simulations <- attr(parMCMC,"simulations")
			n <- NCOL(parMCMC)
			L <- rep(-0.5*n.out*log(2*pi),n)
			for (i in seq(N)){
				if (!("func" %in% names(simulations[[i]])) || any(is.na(simulations[[i]]$func))){
					return(-Inf)
				}
				dimFunc <- dim(simulations[[i]]$func)
				m <- head(dimFunc,2)
				y <- experiments[[i]]$data
				stdv <- standard_error_matrix(y) %otherwise% t(experiments[[i]]$standardError)
				for (k in seq(n)){
					h <- simulations[[i]]$func[,,k]
					dim(h) <- m
					stopifnot(all(dim(h)==dim(y)) && all(dim(y)==dim(stdv)))
					L[k] <- L[k] - 0.5*sum(((y - h)/stdv)^2,na.rm=TRUE) - sum(log(stdv),na.rm=TRUE)
				}
			}
			return(L)
		}
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

#' LOG2 parameter mapping used by the MCMC module
#'
#' This map is used by the simulator to transform sampling variables
#' into ODE-model porameters.
#'
#' @param parMCMC the sampling variables (numeric vector)
#' @export
log2ParMap <- function(parMCMC){
	return(2^(parMCMC))
}

#' LOG2 parameter mapping, jacobian
#'
#' This map is used by the simulator to transform sampling variables
#' into ODE-model porameters. As we often calculate sensitivites, we
#' alos need the jacobian of the map, due to the chain rule of
#' differentiation.
#'
#' @param parMCMC the sampling variables (numeric vector)
#' @export
log2ParMapJac <- function(parMCMC){
	return(diag(2^(parMCMC) * log(2)))
}


#' NATURAL LOG parameter mapping used by the MCMC module
#'
#' This map is used by the simulator to transform sampling variables
#' into ODE-model porameters.
#'
#' @param parMCMC the sampling variables (numeric vector)
#' @export
logParMap <- function(parMCMC){
	return(exp(parMCMC))
}

#' NATURAL LOG parameter mapping, jacobian
#'
#' This map is used by the simulator to transform sampling variables
#' into ODE-model porameters. As we often calculate sensitivites, we
#' alos need the jacobian of the map, due to the chain rule of
#' differentiation.
#'
#' @param parMCMC the sampling variables (numeric vector)
#' @export
logParMapJac <- function(parMCMC){
	return(diag(exp(parMCMC)))
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
gradLogLikelihoodFunc <- function(experiments,parMap=identity,parMapJac=function(x) {diag(1,length(x))})
{
	N <- length(experiments)
	nF <- NROW(experiments[[1]]$data)
	gradLL <- function(parMCMC){
		simulations <- attr(parMCMC,"simulations")
		np <- NROW(parMCMC) # the dimension of the MCMC variable (parMCMC)
		z <- rep(0,np)
		gL <- z
		for (i in seq(N)){
			d <- dim(simulations[[i]]$func)
			nt <- length(experiments[[i]]$outputTimes)
			y <- experiments[[i]]$data
			y[is.na(y)] <- 0.0
			h <- simulations[[i]]$func[,,1]
			dim(h) <- head(d,2)
			stdv <- standard_error_matrix(y)
			stdv2 <- (stdv)^2 # sigma^2
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

#' SMMALA -- The default Extractor of the log-likelihood computed by the simfi solver
#'
#' This function will only extract the log-likelihood value from the
#' solution via simfi. Simfi returns the likelihhod value based on the
#' assumption that the data provided in the list of experiments has a
#' Gaussian standard error. The value includes the normalising factor 1/sqrt(2*pi*sigma^2)
#' for each measured value (the logarithm).
#'
#' @param init a base value that will be added to the log-likelihood
#' @export
#' @return a log-likelihood value for the set of experiments the
#'     solver was set up with.
simfiGaussianLogLikelihood <- function(init = 0.0){
	llf <- function(parMCMC){
		i <- seq_along(parMCMC)
		return(Reduce(\(a,b) a + b$logLikelihood[1], attr(parMCMC,"simulations")))
	}
	return(llf)
}

#' SMMALA -- Extract the gradient of the log-likelihood from the simfi solver's return value
#'
#' This function extracts the approximate gradient of the log-likelihood and
#' transforms the gradient using the Jacobian of the parameter map
#' between Markov chain variables and model parameters.
#'
#' The `simfi()` gradient is with respect to the raw model parameters
#'
#' @param ParMapJac Jacobian of the parMap function
#' @export
#' @return the gradient of the Gaussian log-likelihood with respect to
#'     the MCMC variable
simfiGaussianGradLogLikelihood <- function(ParMapJac=function (x) {diag(1,length(x))}){
	gllf <- function(parMCMC){
		i <- seq_along(parMCMC)
		J <- ParMapJac(parMCMC)
		return(as.numeric(Reduce(\(a,b) a + b$gradLogLikelihood[i,1], attr(parMCMC,"simulations"), init = 0.0) %*% J))
	}
	return(gllf)
}

#' SMMALA -- Extract the approximate Fisher infomration from the simfi results
#'
#' This function extracts the approximate Fisher information matrix `G`
#' of the log-likelihood and transforms it using the
#' Jacobian of the parameter map between Markov chain variables and
#' model parameters.
#'
#' The `simfi()` values are with respect to the raw model parameters,
#' while this function rephrases them in terms of the Markov chain's
#' position.
#'
#' @param ParMapJac Jacobian of the parMap function
#' @export
#' @return the approximate Fisher information of the Gaussian log-likelihood with respect to
#'     the MCMC variable, useful for SMMALA
simfiGaussianFILL <- function(ParMapJac=function (x) {diag(1,length(x))}){
	fi <- function(parMCMC){
		i <- seq_along(parMCMC)
		J <- ParMapJac(parMCMC)
		return(t(J) %*% Reduce(\(a,b) a + b$FisherInformation[i,i,1], attr(parMCMC,"simulations"), init = 0.0) %*% J)
	}
	return(fi)
}

#' High Level SMMALA function
#'
#' This function uses default assumption everywhere and returns a
#' function that will sample from the given model. This funciton will
#' generate code, compile the code, create an ODE solver for it, infer
#' the sampling space from the scale of the parameters, create all
#' necessary functions to move in parameter space (gradients of
#' likelihood and prior), as well as Fisher Information functions.
#'
#' @export
#' @param m the model's TSV representation read via `model_from_tsv`
#' @param o (optional) ode representation of `m`
#' @param ex experiments of `m`, with simulation instructions for `o`.
#' @return `smmala` a function of three arguments: p0, N, eps; where p0 is the starting point, N is the desired sample-size, and eps is the step size.
high_level_smmala <- function(m,o=as_ode(m,cla=TRUE),ex=experiments(m,o)){
	o <- as_ode(m,cla=TRUE)
	C <- generateCode(o)
	odeModel <- write_and_compile(C)
	message(comment(odeModel))

	p <- values(m$Parameter)
	if (all(m$Parameter$scale == "log10")){
		parMap <- log10ParMap
		parMapJac <- log10ParMapJac
	} else if (all(m$Parameter$scale == "log2") || all(m$Parameter$scale == "ld")) {
		parMap <- log2ParMap
		parMapJac <- log2ParMapJac
	} else if (all(m$Parameter$scale == "log") || all(m$Parameter$scale == "ln")){
		parMap <- logParMap
		parMapJac <- logParMapJac
	} else {
		parMap <- identity
		parMapJac <- diag(rep(1.0,length(p)))
	}

	s <- simulator.c(ex,odeModel,parMap=parMap,omit=0)

	dprior <- dNormalPrior(
		mean=m$Parameter$median %otherwise% values(m$Parameter),
		sd=m$Parameter$stdv %otherwise% m$Parameter$sd
	)
	gprior <- gNormalPrior(
		mean=m$Parameter$median %otherwise% values(m$Parameter),
		sd=m$Parameter$stdv %otherwise% m$Parameter$sd
	)
	parMCMC <- mcmcInit(
		beta=1.0,
		values(m$Parameter),
		simulate=s,
		dprior=dprior,
		gradLogLikelihood=gllf(parMapJac),
		gprior=gprior,
		fisherInformation=fi(parMapJac)
	)
	smmala <- mcmc(
		smmalaUpdate(
			simulate=s,
			experiments=ex,
			dprior=dprior,
			gprior=gprior,
			fisherInformation=fi(log10ParMapJac),
			fisherInformationPrior=diag(1.0/m$Parameter$stdv^2)
		)
	)
	attr(smmala,"init") <- parMCMC
	return(smmala)
}
