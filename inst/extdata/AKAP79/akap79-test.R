#!/usr/bin/env Rscript

library(rgsl)
library(uqsa)

comm  <- 0
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"
sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

parMCMC <- c(0.298079363763137, 1.92414841193426, -3.46777437505845, -0.526142972399993, -0.335094425404467, 0.0154698674191168, 0.0589600659244421, -3.61722729859103, -2.81780027714603, -0.6294158469834, -0.337804704712963, -0.321516588604802, 1.32711476805281, -2.46115660220286, -2.43983426718858, 2.23812031125173, -3.72472670364874, -0.436575562361984, 2.33085374648653, -0.185827208095606, -1.06956748735169, -0.371748694335462, -0.810380219155343, 2.76171905370505, 2.09976590481578, 1.72733346358505, -2.19495503976565)
stdv <- rep(2,length(parMCMC))

## ----simulate-----------------------------------------------------------------

sim_with_sensitivity <- list()
sim_without_sensitivity <- list()

for (m in seq(0,10)){
	sim_with_sensitivity[[m+1]] <- simc(ex,modelName,log10ParMap,method=m)
	sim_without_sensitivity[[m+1]] <- simcf(ex,modelName,log10ParMap,method=m)
}

yf <- list()
yS <- list()
b <- list()
for (i in seq(11)){
	t_1 <- Sys.time()
	yS[[i]] <- sim_with_sensitivity[[i]](parMCMC)
	t_1 <- Sys.time() - t_1
	t_2 <- Sys.time()
	yS[[i]] <- sim_without_sensitivity[[i]](parMCMC)
	t_2 <- Sys.time() - t_2
	cat(sprintf("integration method: «%s»\t%f\t%f\t%f\n",nameMethod(i-1),t_1,t_2,t_1/t_2))
	b[[i]] <- rbenchmark::benchmark(
		"with sensitivity"={sim_with_sensitivity[[i]](parMCMC)},
		"without sensitivity"={sim_without_sensitivity[[i]](parMCMC)}
	)
	print(b[[i]])
	flush.console()
}

