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
sensitivity<-function(parSample,outputSample,nBins="Sturges"){
	meanOutput <- colMeans(outputSample)
	varOutput <- diag(cov(outputSample))
	SampleSize <- dim(parSample)
	outputSize <- dim(outputSample)
	hst <- vector("list",SampleSize[2])
	id <- vector("list",SampleSize[2])
	for (i in 1:SampleSize[2]){
		hst[[i]] <- hist(parSample[,i],plot=FALSE,breaks=nBins)
		id[[i]] <- findInterval(parSample[,i],hst[[i]]$breaks)
	}
	# a list, one item per fixed parameter
	binMeans <- lapply(id,observable.mean.in.bin,outputSample=outputSample)
	S <- matrix(0,outputSize[2],SampleSize[2])
	for (i in 1:SampleSize[2]){
		Vi <- sum.of.bin.variance(hst[[i]],binMeans[[i]],totalMean=meanOutput)
		S[,i] <- Vi/varOutput
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
#' @param u the values of the x-axis for the plot
#' @param S the sensitivity matrix as returned by `sensitivity()`
#' @return nothing
sensitivity.graph <- function(u,S,color=rainbow(dim(S)[2]),do.sort=TRUE,title="Sensitivity"){
	d <- dim(S)
	n <- d[2]-1
	if (do.sort) {
		m <- colMeans(S)
		I <- order(m,decreasing=TRUE)
		S <- S[,I]
		ylabel <- "sorted cumulative sensitivity"
	} else {
		ylabel <- "cumulative sensitivity"
	}
	C <- t(apply(S,1,cumsum))
	x <- c(u,rev(u))

	plot(u,C[,1],type='l',ylim=c(0,max(C)*1.1),ylab=ylabel,xlab="output",main=title,axes=FALSE)
	axis(1,at=u,labels=names(u))
	axis(2)
	z <- c(S[,1]*0,rev(S[,1]))
	polygon(x,z,col=color[1],lty=0)
	print(n)
	for (i in 1:n){
		y <- c(C[,i],rev(C[,i+1]))
		polygon(x,y,col=color[i+1],lty=0)
		lines(u,C[,i],col=color[i+1],lwd=2)
	}
	legend(x="topright",fill=color[1:d[2]],legend=colnames(S),ncol=2)
}
