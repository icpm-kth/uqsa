#' copulaPrior creates a prior probability density function
#'
#' This function accepts the return list of fitCopula() or
#' makeIndepCopula() and creates a density function from it.
#'
#' @export
#' @importFrom VineCopula RVinePDF
#' @param Copula a list, as returned by fitCopula() or makeIndepCopula
#' @return a function that maps parameters (a vector) to probability density values (scalar)
#' @examples
#' prior_pdf<-dCopulaPrior(makeIndepCopula(ll=c(-1,-1,-1),ul=c(1,1,1)))
#' prior_pdf(runif(3))
#' 0.125
dCopulaPrior <- function(Copula){
	U <- Copula$U
	Y <- Copula$Y
	Z <- Copula$Z
	copula <- Copula$copula
	np  <- ncol(U)
	priorPDF<-function(inx){
		if(all(!is.na(inx))){
		  lbU <- sapply(1:np, function(i) min(U[,i]) - 0.5*c(-1,1) %*% range(U[,i]))
		  ubU <- sapply(1:np, function(i) max(U[,i]) + 0.5*c(-1,1) %*% range(U[,i]))
		  ed <- sapply(1:np, function(i) approx(c(unique(U[,i]), lbU[i], ubU[i]), c(unique(Z[,i]), 0, 1), xout=inx[i])$y)
			mpdf <- sapply(1:np, function(i) approx(c(unique(U[,i]), lbU[i], ubU[i]), c(unique(Y[,i]), 0, 0), xout=inx[i])$y)
			if(any(is.na(ed)) || any(is.na(mpdf))){ # outside of copula defined limits
				jpdf <- 0
			} else {
				jpdf <- RVinePDF(ed, copula, verbose = TRUE)*prod(mpdf)
			}
		} else {
			jpdf <- 0
		}
		return(jpdf)
	}
	return(priorPDF)
}

#' rCopulaPrior returns a function that generates random values from the copula model
#'
#' The returned function generates n random vectors, as rows of a matrix.
#'
#' @importFrom VineCopula RVineSim
#' @export
#' @param Copula the return value of fitCopula()
#' @return a matrix of random values
#' @examples
#' rcp<-rCopulaPrior(makeIndepCopula(ll=c(0,1,2),ul=c(1,2,3)))
#' rcp(12)
#' # this returns a matrix of random values, the first column has values between 0 and 1
#' # the second column has values between 1 and 2, etc.
rCopulaPrior <- function(Copula){
  copula <- Copula$copula
  Z <- Copula$Z
  U <- Copula$U
  np <- ncol(Z)
  rprior <- function(npc){
    R <- RVineSim(npc, copula)
    prePar <- matrix(0, npc, np)
    for(i in 1:np){
      prePar[,i] = spline(unique(Z[,i]),unique(U[,i]),xout=R[,i])$y
    }
    return(prePar)
  }
  return(rprior)
}

#' dUniformPrior creates a uniform density function
#'
#' The returned denisty function takes vectors of the same size as ll
#' and ul. It returns the product of the component's one-dimensional
#' uniform distribtions.
#'
#' @export
#' @param ll lower limit of the random variables (a vector)
#' @param ul upper limit of the random variables (same size vector as ll)
#' @return a probability density function on vectors withthe same length as ll and ul.
#' @examples
#' dup<-dUniformPrior(ll=c(0,1,2),ul=c(1,2,3))
#' dup(c(0.5,1.5,2.5))
#' [1] 1
dUniformPrior <- function(ll,ul){
  dprior <- function(x){
    return(prod(dunif(x,min=ll,max=ul)))
  }
  return(dprior)
}

#' rUniformPrior returns a random vector generator
#'
#' The return value is a function that generates random vectors of the
#' same size as ll and ul from a uniform distribution within the
#' limits defined by ul and ll. The random vectors are returned as n
#' rows of a matrix, where n is the only argument of the returned
#' function.
#'
#' @export
#' @param ll lower limit of the random variables (a vector)
#' @param ul upper limit of the random variables (same size vector as
#'     ll)
#' @return a uniform random vector generating function: runiform(n),
#'     where n is the requested number of vectors (rows)
#' @examples
#' rup<-rUniformPrior(ll=c(0,1,2),ul=c(1,2,3))
#' rup(12)
#'           [,1]     [,2]     [,3]
#'  [1,] 0.7879692 1.344263 2.452547
#'  [2,] 0.7154937 1.846704 2.174703
#'  [3,] 0.6644920 1.257256 2.764870
#'  [4,] 0.3470803 1.524026 2.715229
#'  [5,] 0.7879692 1.344263 2.452547
#'  [6,] 0.7154937 1.846704 2.174703
#'  [7,] 0.6644920 1.257256 2.764870
#'  [8,] 0.3470803 1.524026 2.715229
#'  [9,] 0.7879692 1.344263 2.452547
#' [10,] 0.7154937 1.846704 2.174703
#' [11,] 0.6644920 1.257256 2.764870
#' [12,] 0.3470803 1.524026 2.715229
rUniformPrior <- function(ll,ul){
  np <- length(ll)
  stopifnot(np==length(ul))
  rprior <- function(n){
    r <- matrix(runif(n*np,min=ll,max=ul),n,np,byrow=TRUE)
    return(r)
  }
  return(rprior)
}



#' dNormalPrior creates the density function of a multivariate normal distribution with independent components
#'
#' The returned density function takes vectors of the same size as mean and sd.
#' It returns the product of the components' one-dimensional normal distribution,
#' with mean "mean" and standard deviation "sd".
#'
#' @export
#' @param mean mean of the random variables (a vector)
#' @param sd standard deviation of the random variables (same size vector as mean)
#' @return a probability density function on vectors with the same length as mean and sd.
#' @examples
#' dnp<-dNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
#' dnp(c(0.5,1.5,2.5))
#' [1] 0.008926651
dNormalPrior <- function(mean,sd){
  dprior <- function(x){
    return(prod(dnorm(x, mean=mean, sd = sd)))
  }
  return(dprior)
}

#' rNormalPrior returns a random vector generator
#'
#' The return value is a function that generates random vectors of the
#' same size as mean and sd from a multivariate normal distribution with 
#' independent components with mean "mean" and standard deviation "sd".
#' The random vectors are returned as n
#' rows of a matrix, where n is the only argument of the returned
#' function.
#'
#' @export
#' @param mean mean of the random variables (a vector)
#' @param sd standard deviation of the random variables (same size vector as
#'     mean)
#' @return an independentent multivariate normal random vector generating function: rprior(n),
#'     where n is the requested number of vectors (rows)
#' @examples
#' rnp<-rNormalPrior(mean=c(0,1,2),sd=c(1,2,3))
#' rnp(12)
#'             [,1]       [,2]       [,3]
#' [1,]  0.99368106  1.0074638  0.2411801
#' [2,] -1.39127288  0.3608953 -0.7704979
#' [3,] -0.23957518 -0.4376285  7.2133888
#' [4,]  0.55225848 -3.0625581  2.1468860
#' [5,] -1.74950066 -1.9512875 -0.1923920
#' [6,]  0.36377999 -0.2222019  4.0832047
#' [7,]  0.62605893 -0.8562764 -0.2743055
#' [8,] -0.99829859  2.8863743  7.4753572
#' [9,] -0.71302121  2.8879101  1.4377714
#' [10,]  0.04040756 -1.6638517  0.7293233
#' [11,] -0.04951403  1.9498103  2.1706586
#' [12,]  0.49326957  2.7241179  7.4083735

rNormalPrior <- function(mean,sd){
  np <- length(mean)
  stopifnot(np==length(sd))
  rprior <- function(n){
    r <- matrix(rnorm(n*np,mean=mean,sd=sd),n,np,byrow=TRUE)
    return(r)
  }
  return(rprior)
}



