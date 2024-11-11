#!/usr/bin/env Rscript
f  <- dir(pattern="[.]tsv$",full.names=TRUE)
sb <- SBtabVFGEN::sbtab_from_tsv(f)
ex <- SBtabVFGEN::sbtab.data(sb)[seq(5)]
names(ex)
modelName <- uqsa::checkModel(comment(sb),paste0("./",comment(sb),".so"))

source("CaMKIIs.R")
p <- CaMKIIs_default()
u <- ex[[1]]$input
s <- uqsa::simulator.c(ex,modelName)

par <- head(p,-length(u))
plot(log10(p),type='l')
lines(log10(c(par,u+1e-20)))
print(par)
print(ex[[1]]$input)


## 1. Tthis is a normal simulation, as it would be done during MCMC
T <- Sys.time()
y <- s(as.matrix(par))
print(difftime(T,Sys.time()))

for (i in seq_along(ex)){
	dev.new()
	d <- ex[[i]]$events$dose
	o <- apply(is.na(ex[[i]]$outputValues),2,any)
	plot(d,y[[i]]$func[!o,,1],type='l',xlab='dose',ylab=names(ex[[i]]$outputValues)[!o],main=names(ex)[i])
}

## 2. This is a simulation with time points inserted between the measurement times, to see what happens between them
F <- 1 # time dilation factor, this is to try different durations to reach steady state
for (i in seq_along(ex)){
	t_ <- ex[[i]]$outputTimes
	ex[[i]]$events$time <- ex[[i]]$events$time*F
	ex[[i]]$outputTimes <- sort(c(t_,seq(0,max(t_),length.out=length(t_)*10)))*F
}

s2 <- uqsa::simulator.c(ex,modelName)
 T <- Sys.time()
y2 <- s2(as.matrix(par))
print(difftime(T,Sys.time()))

for (i in seq_along(ex)){
	dev.new()
	o <- apply(is.na(ex[[i]]$outputValues),2,any)
	plot(ex[[i]]$outputTimes,y2[[i]]$func[!o,,1],type='l',xlab='time',ylab=names(ex[[i]]$outputValues)[!o],main=names(ex)[i])
}
