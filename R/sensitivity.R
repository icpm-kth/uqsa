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
	parSample <- parSample[!isNA,]
	outputSample <- outputSample[!isNA,]
	meanOutput <- colMeans(outputSample)
	varOutput <- diag(cov(outputSample))
	SampleSize <- dim(parSample)
	outputSize <- dim(outputSample)
	hst <- vector("list",SampleSize[2])
	id <- vector("list",SampleSize[2])
	for (i in 1:SampleSize[2]){
		hst[[i]] <- hist(parSample[,i],plot=FALSE,breaks=nBins)
		id[[i]] <- findInterval(parSample[,i],hst[[i]]$breaks,all.inside=TRUE)
	}
	# a list, one item per fixed parameter
	binMeans <- lapply(id,observable.mean.in.bin,outputSample=outputSample)
	S <- matrix(0.0,outputSize[2],SampleSize[2])
	for (i in 1:SampleSize[2]){
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
#' The state sensitivity matrix:
#'
#' ```
#'            d state(time[k],state, param)[i]
#' S[i,j,k] = --------------------------------  ,
#'            d param[j]
#' ```
#'
#' where param are the raw model parameters.
#' This matrix is calculated as an intermediate and then transformed into:
#'
#' ```
#'             d func(time[k], state, c(parMap(parMCMC),input))[i]
#' Sh[i,j,k] = --------------------------------------------------
#'             d parMCMC[j]
#' ```
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
#' Typically, the sensitivity needs to be known at different
#' time-points t_k. The 3-dimensional array S[i,j,k], where the index k corrsponds to time
#' t_k; the closer x(t_k) is to equilibrium, the better the
#' approximation; near the initial state, the sensitivity is also
#' correct (only the intermediate time-span is approximate).
#'
#' This function requires pracma::expm to work.
#'
#' The f
#'
#' @export
#' @param experiments a list of simulation experiments
#' @param simulations an equivalent list of simulation results, for
#'     one parameter vector
#' @param model a list of functions for the model the experiments are
#'     applicable to
#' @param parMCMC the parameters that are used in Markov chain Monte
#'     Carlo as the MC variable
#' @param parMap a map to transform parMCMC into p, parameters the
#'     model accepts
#' @return a function S(parMCMC) ->
#'     simulations_with_sensitivity, which attaches the state
#'     sensitivity matrix array length(x) × length(p) × length(t) to
#'     the simulations (solutions to the ODE).
#' @examples \dontrun{
#' y <- simulate(parMCMC)
#' S <- sensitivityEquilibriumApproximation(experiments, model, parMap, parMapJac)
#' y <- S(parMap,y)
#' }
sensitivityEquilibriumApproximation <- function(experiments, model, parMap=identity, parMapJac=1.0){
	y0 <- model$init(0.0)
	n  <- length(y0)
	m  <- ncol(experiments[[1]]$outputValues)
	u <- experiments[[1]]$input
	nu <- length(u)
	defaultPar <- c(model$par(),u)
	np <- length(defaultPar)
	d <- np - nu
	N <- length(experiments)
	SEA <- function(parMCMC,simulations){
		stopifnot(length(simulations) == N)
		for (i in seq(N)){
			p <- c(parMap(parMCMC),experiments[[i]]$input)
			tm <- experiments[[i]]$outputTimes
			t0 <- experiments[[i]]$initialTime
			simulations[[i]]$sens <- array(0.0,dim=c(n,length(parMCMC),length(tm)))
			simulations[[i]]$funcsens <- array(0.0,dim=c(m,length(parMCMC),length(tm)))
			for (j in seq(length(t))){
				if (abs(tm[j]-t0) > 1e-16 * abs(t0)){
					y <- simulations[[i]]$state[,j,1]
					A <- model$jac(tm[j], y, p)
					B <- head(model$jacp(tm[j], y, p), c(n,d))
					AB <- solve(A,B)
					# state variables:
					C <- ((pracma::expm((t[j]-t0)*A) %*% AB) - AB)
					simulations[[i]]$sens[,,j] <- C %*% parMapJac(parMCMC)
					# functions
					CF <- model$funcJac(tm[j],y,p) %*% C + head(model$funcJacp(tm[j],y,p),c(m,d))
					simulations[[i]]$funcsens[,,j] <- CF %*% parMapJac(parMCMC)
				}
			}
		}
		return(simulations)
	}
	return(SEA)
}
