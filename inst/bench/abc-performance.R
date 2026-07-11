#!/usr/bin/env Rscript

library(uqsa)
library(errors)
library(bench)

f <- uqsa_example("AKAR4")
m <- model_from_tsv(f)
o <- as_ode(m)
c_path(o) <- write_c_code(generate_code(o))
so_path(o) <- shlib(o)

ex <- experiments(m,o)
time_out_seconds <- 2.5
integrator_step_limit <- 10000
s <- simulator.c(
	ex,             # experiments
	o,              # the model
	log10ParMap,    # reverse map, to get back from Markov chain space to model
	omit=2,         # omit Fisher-Information and Gradient calculations
	num.steps=intergrator_step_limit,
	time.out=time_out_seconds
)
Obj <- makeObjective(ex,s)

p0 <- log10(values(m$Parameter))
dprior <- dUniformPrior(p0-3,p0+3)
rprior <- rUniformPrior(p0-3,p0+3)

P <- p0 + matrix(rnorm(3*100),3,100)

B <- bench::mark(
	"ABCMCMC" = { # old algorithm
		ret <- ABCMCMC(Obj, p0, 100, Sigma0=cov(t(P))*0.1,delta=1,dprior=dprior, allow.reg.TRUE)
		ACF <- acf(ret$scores)$acf
		tau <- sum(A[A>0.2])
	}
	"abc mcmc" = { # new algorithm
		ret <- abc_mcmc(Obj,P,100,burnIn=50,Sigma0=cov(t(P))*0.1,dprior=dprior)
		ACF <- acf(ret$distances)$acf
		tau <- sum(A[A>0.2])
	},
	"abc mcmc, no batches" = { # all adaptive features turned off
		ret <- abc_mcmc(Obj,P,10000,burnIn=0,Sigma0=cov(t(P))*0.1,batchSize=1,dprior=dprior)
		ACF <- acf(ret$distances)$acf
		tau <- sum(A[A>0.2])
	},
	max_iterations=3,
	min_time=Inf,
	check=\(a,b){TRUE},
	memory = FALSE
)
##hexbin::hexplom(ret$draws)

tau <- as.data.frame(Reduce(\(a,b) c(a,b),B$result))
v <- 10000/(2*tau)
B <- B |> dplyr::mutate(effective_speed=v)
print(B[,c("expression","median","effective_speed")])
