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