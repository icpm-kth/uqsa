#' Kullback Leibler Divergence
#'
#' This function calculates the Kullback Leibler Divergence $D(P\|Q)$ value
#' between two distributions $P$ and $Q$, represented by their two samples, `X` and
#' `Y`. Both samples will have their density inferred.
#' @param X sample from distribution P
#' @param Y sample from distribution Q
#' @return D, a scalar value, the Kullback Leibler Divergence
#' @export
KLD <- function(X,Y,use.kde=FALSE){
	if (use.kde) { # quick and dirty KDE
		density_fn <- function(ZETA, h=NULL) {
			if (missing(h) || is.null(h)) {
				SD <- apply(ZETA,2,sd)
				iqr <- apply(ZETA,2,IQR)
				ssize <- NROW(ZETA)
				n <- NCOL(ZETA)
				h <- 0.9*pmin(SD,iqr/1.34)*exp(-1/(2*n+1)*log(ssize)) # crude rule of thumb
			}
			return(
				function(xi) { # p_Y(X)
					sapply(
						seq(NROW(xi)),
						\(i) {
							mean(mvtnorm::dmvnorm(ZETA, mean = xi[i, ], sigma = diag(h)))
						}
					)
				}
			)
		}
		P <- density_fn(X)
		Q <- density_fn(Y)
	} else {
		P <- \(Z) apply(Z,1,dCopulaPrior(fitCopula(X))) # density estimate for X
		Q <- \(Z) apply(Z,1,dCopulaPrior(fitCopula(Y))) # density estimate for Y
	}
	D <- mean(log(P(X))-log(Q(X))) # X ~ P
	return(D)
}

#' Multivariate Normal Distribution KLD
#'
#' Like the [KLD] function, this function calculates KLD values, but
#' for the specific case of multivariate normal distributions
#' $\mathcal{N}_{A}$ and $\mathcal{N}_{B}$.  The two distributions are
#' specified using $\mu$ and $\Sigma$ values (mean and covariance).
#'
#' @param muA
#' @param SigmaA
#' @param muB
#' @param SigmaB
#' @return the KLD value D(A|B)
#' @export
exact_normal_kld <- function(muA,SigmaA,muB,SigmaB){
	tr <-\(M) sum(diag(M))
	k <- length(muA)
	log_det_A <- as.numeric(determinant(SigmaA, logarithm = TRUE)$modulus)
	log_det_B <- as.numeric(determinant(SigmaB, logarithm = TRUE)$modulus)
	D <- tr(solve(SigmaB,SigmaA)) - k + t(muB-muA) %*% solve(SigmaB,muB-muA) + log_det_B - log_det_A
	return(0.5*drop(D))
}
