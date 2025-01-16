#' The mean value of an observable value given a parameter bin
#'
#' The model parameters were binned with a histogram function. In each
#' bin one of the parameters is almost fixed (it varies much less than
#' the other parameters [full range]). This function returns the mean
#' of the observable for each bin, as a vector.
#'
#' @param id is an integer vector that identifies the bin the
#'     parameter vector falls into to create the same row of the
#'     outputSample (the output stems from a model simulation with
#'     parameters). id has the same length as outputSample rows.
#' @param outputSample a matrix of output values, one output vector
#'     per row (different rows are results at different parameter
#'     values)
#' @return M[i,j] the mean of each observable[j] in bin[i]
observable.mean.in.bin <- function(id,outputSample){
	d<-dim(outputSample)
	stopifnot(length(id)==d[1])
	n<-max(id)
	bin.mean <- matrix(NA,n,d[2])
	for (i in 1:n){
		l <- id == i
		if (any(l)){
			bin.mean[i,] <- colMeans(outputSample[l,,drop=FALSE])
		}
	}
	return(bin.mean)
}

#' Weighted Sum of Bin-specific variances
#'
#' This function calculates the variance sum of a vector valued
#' observable.
#'
#' @usage sum.of.bin.variance(hst, binMeans, totalMean)
#' @param hst the histogram of the parameter sample
#' @param binMeans the means of the observable within each bin (rows of means)
#' @param totalMean the mean of the observable over the entire sample (vector)
#' @return The weighted sum of square differences between the binMean and the totalMean
sum.of.bin.variance  <- function(hst,binMeans,totalMean){
	B <- dim(binMeans) # binning dimensions
	stopifnot(B[2] == length(totalMean))
	return(colSums(hst$counts*(t(t(binMeans)-totalMean))^2,na.rm=TRUE)/sum(hst$counts))
}

#' Global Sensitivity Analysis
#'
#' This function performs a binning based estimation of the global
#' sensitivity of a model's output with respect to the model's
#' parameters. The output can be a prediction of the model's behaviour
#' in a scenario of interest (parameters, input, intial values,
#' boundary conditions, scheduled events etc.). The output models a
#' potentially measurable value (the "observable"). The sample-rows
#' and the output rows must correspond (they must be from the same
#' model simulation).
#'
#' @param parSample a matrix of parameter vectors (rows)
#' @param outputSample a matrix, with rows of outputs (row-index is the sample index)
#' @param nBins number of bins, if unset defaults to the default of the hist function
#' @export
#' @return sensitivity S[i,j] of output[i] with respect to parameter[j]
globalSensitivity<-function(parSample,outputSample,nBins="Sturges"){
	isNA <- apply(is.na(outputSample),1,any)
	parSample <- parSample[!isNA,,drop=FALSE]
	outputSample <- outputSample[!isNA,,drop=FALSE]
	meanOutput <- colMeans(outputSample)
	varOutput <- diag(cov(outputSample))
	outputSize <- dim(outputSample)
	hst <- vector("list",NCOL(parSample))
	id <- vector("list",NCOL(parSample))
	for (i in 1:NCOL(parSample)){
		hst[[i]] <- hist(parSample[,i],plot=FALSE,breaks=nBins)
		id[[i]] <- findInterval(parSample[,i],hst[[i]]$breaks,all.inside=TRUE)
	}
	# a list, one item per fixed parameter
	binMeans <- lapply(id,observable.mean.in.bin,outputSample=outputSample)
	S <- matrix(0.0,NCOL(outputSample),NCOL(parSample))
	for (i in 1:NCOL(parSample)){
		Vi <- sum.of.bin.variance(hst[[i]],binMeans[[i]],totalMean=meanOutput)
		S[,i] <- Vi/(1e-300+varOutput)
	}
	colnames(S) <- colnames(parSample)
	rownames(S) <- colnames(outputSample)
	return(S)
}

#' plot the sensitivity matrix
#'
#' Produce a cumulative shaded area plot for the sensitivity matrix.
#'
#' @export
#' @param u the values of the x-axis for the plot, if named, the names
#'     are put at the tick-marks
#' @param S the sensitivity matrix as returned by `globalSensitivity()`,
#'     S[i,j] is with respect to model output i and parameter j
#' @param color the list of colors to use for the shaded areas, e.g.:
#'     rainbow(24)
#' @param line.color the color of the lines drawn between the shaded
#'     areas
#' @param do.sort the parameter sensitivities are sorted according to
#'     the mean over all outputs, the parameter with the most
#'     sensitivity is plotted first, at the bottom
#' @param decreasing direction of sort, the first item in the sorted
#'     list (the parameter) will be plotted first, and thus at the
#'     bottom of the plot
#' @param title string, written above, as a title
#' @return nothing
sensitivity.graph <- function(u,S,color=hcl.colors(dim(S)[2]),line.color=hcl.colors(dim(S)[2]+1),do.sort=TRUE,decreasing=FALSE,title="Sensitivity"){
	d <- dim(S)
	n <- d[2]-1
	if (do.sort) {
		m <- colMeans(S)
		I <- order(m,decreasing=decreasing)
		S <- S[,I]
		ylabel <- "sorted cumulative sensitivity"
	} else {
		ylabel <- "cumulative sensitivity"
	}
	C <- t(apply(S,1,cumsum))
	x <- c(u,rev(u))

	plot(u,C[,1],type='l',ylim=c(0,max(C)*1.1),ylab=ylabel,xlab="output",main=title,axes=FALSE,col=line.color[1])
	axis(1,at=u,labels=names(u))
	axis(2)
	z <- c(S[,1]*0,rev(S[,1]))
	polygon(x,z,col=color[1],lty=0)
	print(n)
	for (i in 1:n){
		y <- c(C[,i],rev(C[,i+1]))
		polygon(x,y,col=color[i+1],lty=0)
		lines(u,C[,i],col=line.color[i+1],lwd=2)
	}
	legend(x="topright",fill=color[1:d[2]],legend=colnames(S),ncol=2)
}


#' Outpts the random sample on which to perform the Sobol-Homma-Saltelli global sensitivity
#' analysis
#'
#' Each parameter vector has length nPars,
#' The sample consists of two random (nSamples x nPars) matrices M1,
#' M2 and a third (nSamples x nPars x nPars) array N. N consists of
#' nPars copies of M2, except that in each M2-matrix one column has
#' been replaced by the corresponding column of M1. M1 and M2 consists
#' of random numbers from a normal distribution.
#'
#' These matrices provide prior distribution samples to be further
#' processed by the simulator, similar to this:
#'
#' ```
#' sim <- simulator.c(experiments,modelName)
#' fM1 <- t(sim(t(M1))[[1]]$state[,ti,])       # or similar
#' ```
#'
#' For details see: Halnes, Geir, et al. J. comp. neuroscience 27.3 (2009): 471.
#'
#' @export
#' @param nSamples number of rows to return
#' @param rprior a function that samples from the prior distribution
#' @return a list with the components M1, M2 (both matrices) and N (a
#'     3D-array).
shs_prior <- function(nSamples,rprior){
	nPars<-NCOL(rprior(1))
	M1 <- rprior(nSamples)
	M2 <- rprior(nSamples)
	N <- array(NA, dim=c(nSamples,nPars,nPars))
	for (i in 1:nPars){
		# Replace the i:th column in M2 by the i:th column from M1 to obtain Ni
		N[,,i] <- M2
		N[,i,i] <- M1[,i]
	}
	return(list(M1=M1,M2=M2,N=N))
}

#' Outputs the global sensitivity scores SI and SIT, calculated by the Sobol-Homma-Saltelli method
#'
#' M1, M2, and N are matrices prepared by `uqsa::shs_prior()`.  The
#' parameters (rows) from these matrices need to be simulated (using
#' any method), to obtain fM1, fM2 and fN.
#'
#' These matrices are shaped similarly to M1, M2 and N respectively,
#' but now the parameters are replaced by the effects they have on a
#' observable of interest (the output). It can be the vector
#' valued output at a specific (single) time-point or a scalar output
#' at different time-points.
#'
#' See Geir Halnes et al. (Halnes, Geir, et al. J. comp. neuroscience 27.3 (2009): 471.
#'
#' @export
#' @param fM1 output (f)unction values for M1, nSamples × nOuts
#' @param fM2 output (f)unction values for M2, nSamples × nOuts
#' @param fN output (f)unction values for N, nSamples × nOuts × nPars
#' @return a list with sensitivity indices $SI and total sensitivities $SIT
shs_gsa<- function(fM1,fM2,fN, subtractMean = TRUE){
	nSamples <- dim(fM1)[1]
	nOuts <- dim(fM1)[2]
	nPars <- dim(fN)[3]
	if(subtractMean){
		fM1 <- fM1 - matrix(colMeans(fM1), nrow=nSamples, ncol=nOuts, byrow=1)
		fM2 <- fM2 - matrix(colMeans(fM2), nrow=nSamples, ncol=nOuts, byrow=1)
		for (i in 1:nPars){
			fN[,,i] <- fN[,,i] - matrix(colMeans(fN[,,i]), nrow=nSamples, ncol=nOuts, byrow=1)
		}
	}

	EY2 <- colMeans(fM1*fM2)
	VY <- colSums(fM1*fM1)/(nSamples-1) - EY2
	VYT <- colSums(fM2*fM2)/(nSamples-1) - EY2

	SI <- matrix(0, nrow=nOuts,ncol=nPars)
	SIT <- matrix(0, nrow=nOuts,ncol=nPars)
	for(i in 1:nPars){
		SI[,i] <- (colSums(fM1*fN[,,i])/(nSamples-1) - EY2)/VY
		SIT[,i] <- 1 - (colSums(fM2*fN[,,i])/(nSamples-1) - EY2)/VYT
	}
	return(list(SI=SI,SIT=SIT))
}

