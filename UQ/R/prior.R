#' copulaPrior creates a prior probability density function
#'
#' This function accepts the return list of fitCopula() or
#' makeIndepCopula() and creates a density function from it.
#'
#' @export
#' @param Copula a list, as returned by fitCopula() or makeIndepCopula
#' @return a function that maps parameters (a vector) to probability density values (scalar)
#' @examples
#' prior_pdf<-copulaPrior(makeIndepCopula(ll=c(-1,-1,-1),ul=c(1,1,1)))
#' prior_pdf(runif(3))
#' 0.125
copulaPrior <- function(Copula){
	U <- Copula$U
	Y <- Copula$Y
	Z <- Copula$Z
	copula <- Copula$copula
	ll <- apply(U,2,min)
	ul <- apply(U,2,max)
	np  <- length(ll)
	priorPDF<-function(inx){
		if(all(!is.na(inx))){
			ed <- sapply(1:np, function(i) approx(U[,i], Z[,i], xout=inx[i])$y)
			mpdf <- sapply(1:np, function(i) approx(U[,i], Y[,i], xout=inx[i])$y)
			if(any(is.na(ed)) || any(is.na(mpdf))){ # outside of copula defined limits
				jpdf <- 0
			}else if(!(all(inx >=ll) && all(inx<=ul))){ # outside of prior
				jpdf <- 0
			}else{
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
#' @export
#' @param Copula the return value of fitCopula()
#' @return a matrix of random values
rCopulaPrior <- function(Copula){
	copula <- Copula$copula
	Z <- Copula$Z
	U <- Copula$U
	np <- ncol(Z)
	rprior <- function(npc){
		R <- RVineSim(npc, copula)
		prePar <- matrix(0, npc, np)
		for(i in 1:np){
			prePar[,i] = spline(Z[,i],U[,i],xout=R[,i])$y
		}
		return(prePar)
	}
	return(rprior)
}
