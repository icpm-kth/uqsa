#!/usr/bin/env Rscript

library(uqsa)
library(mvtnorm)
library(bench)
library(mclust)
library(ks)
library(errors)

muA <- c(0,0)
SigmaA <- matrix(c(1,0.2,0.2,1),2,2)
muB <- c(0.1,-0.2)
SigmaB <- matrix(c(1,0.4,0.4,1),2,2)
set.seed(42)

f <- function(method="copula"){
	N <- 1e3
	n <- 6 # repetitions
	E <- exact_normal_kld(muA,SigmaA,muB,SigmaB)
	D <- replicate(n,
		KLD(
			rmvnorm(N,muA,SigmaA),
			rmvnorm(N,muB,SigmaB),
			method
		)-E
	)
	ae <- set_errors(mean(abs(D)),sd(abs(D))/sqrt(n))
	return(ae)
}

B <- bench::mark(
	fitCopula=f("copula"),
	mclust=f("mclust"),
	ks=f("ks"),
	mvtnorm=f("mvtnorm"),
	max_iterations=3,
	min_time=Inf,
	check=\(a,b){abs(a-b)<1}
)

err <- as.data.frame(Reduce(\(a,b) c(a,b),B$result))
B <- B |> dplyr::mutate(abs_err=Reduce(\(a,b) c(a,b),B$result))
print(B[,c("expression","median","abs_err")])
