#!/usr/bin/env Rscript
library(parallel)

f  <- dir(pattern="[.]tsv$",full.names=TRUE)
sb <- SBtabVFGEN::sbtab_from_tsv(f)
ex <- SBtabVFGEN::sbtab.data(sb)
options(mc.cores = length(ex))

modelName <- uqsa::checkModel(comment(sb),paste0("./",comment(sb),".so"))
source("CaMKIIs.R")
p <- CaMKIIs_default()
u <- ex[[1]]$input
par <- head(p,-length(u))

## plot(log10(p),type='l')
## lines(log10(c(par,u+1e-20)))

F <- c(1,1,1,1,1,1) # time dilation factor, this is to try different durations to reach steady state
Methods <- seq(0,10)
for (M in Methods[-3]){
	ex <- SBtabVFGEN::sbtab.data(sb)

	## 1. Tthis is a normal simulation, as it would be done during MCMC
	for (i in seq_along(ex)){
		t_ <- ex[[i]]$outputTimes
		ex[[i]]$_t <- t_
		ex[[i]]$events$time <- ex[[i]]$events$time*F[i]
		ex[[i]]$outputTimes <- t_*F[i]
	}
	s <- uqsa::simulator.c(ex,modelName,method = M)

	T <- Sys.time()
	y <- tryCatch(s(as.matrix(par)),error = function(e) {print(e); return(NA)})
	Success <- unlist(lapply(y,\(y) !any(is.na(y$func))))
	cpuSeconds <- unlist(lapply(y,function(l){l$cpuSeconds}))
	cat(M,rgsl::nameMethod(M),Success,difftime(Sys.time(),T,units="mins"),cpuSeconds,"\n")
	if (length(y)==length(ex)){
		dev.new()
		par(mfrow=c(2,6))
		for (i in seq_along(ex)){
			d <- ex[[i]]$events$dose
			o <- apply(is.na(ex[[i]]$outputValues),2,any)
			if (!any(is.na(y[[i]]$func))){
				plot(d,y[[i]]$func[!o,,1],type='l',xlab='dose',ylab=names(ex[[i]]$outputValues)[!o],main=names(ex)[i])
			} else {
				plot.new()
			}
		}
		mtext(rgsl::nameMethod(M), side = 3)
	}

	## 2. This is a simulation with time points inserted between the measurement times, to see what happens between them
	for (i in seq_along(ex)){
		t_ <- ex[[i]]$outputTimes
		ex[[i]]$outputTimes <- sort(c(t_,seq(min(0,min(t_)),max(t_),length.out=length(t_)*10)))
	}

	s2 <- uqsa::simulator.c(ex,modelName,method=M)
	 T <- Sys.time()
	tryCatch(y2 <- s2(as.matrix(par)),error = function(e) {print(e); return(NA)})
	cpuSeconds <- unlist(lapply(y,function(l){l$cpuSeconds}))
	#cat("max(F)=",max(F)," M=",M," time=",difftime(Sys.time(),T,units="mins"),"minutes"," cpuSeconds=[",cpuSeconds,"]\n")
	if (length(y)==length(ex)){
		for (i in seq_along(ex)){
			o <- apply(is.na(ex[[i]]$outputValues),2,any)
			if (!any(is.na(y2[[i]]$func))){
				plot(ex[[i]]$outputTimes,y2[[i]]$func[!o,,1],type='l',xlab='time',ylab=names(ex[[i]]$outputValues)[!o],main=names(ex)[i])
			}
		}
	}
}
