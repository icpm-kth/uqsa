#!/usr/bin/env Rscript
library(uqsa)
library(SBtabVFGEN)
library(ggplot2)

files <- commandArgs(trailingOnly=TRUE)
modelName <- "CaMKIIs"
comment(modelName) <- "./CaMKIIs.so"
#sz <- 1e3
cat("Loading sample\n")
PREFIX <- paste(
	Reduce(
	function(a,b) {
		m<-seq(min(length(a),length(b)));
		a<-a[m]; b<-b[m];
		i<-which(unlist(mapply(identical,a,b)));
		return(a[i])},
	strsplit(sub("[.]RDS$",'',files),'-')),
collapse='-')
print(PREFIX)

x <- uqsa::gatherSample(files,beta=1.0)
l <- attr(x,"logLikelihood")
I <- order(l,decreasing=TRUE)[seq(round(length(l)/2))]
MLE <- I[1]
x <- x[I,]
l <- l[I]
attr(x,"logLikelihood") <- l
MLE <- 1
sb <- sbtab_from_tsv(dir(pattern="[.]tsv$",full.names=TRUE))
ex <- sbtab.data(sb)

## ---- R functions for this model:
source(sprintf("%s.R",modelName))
t0 <- -1
p0 <- CaMKIIs_default(t0)
y0 <- CaMKIIs_init(t0,p0)
C <- NCOL(x)
colnames(x) <- names(p0[seq(C)])

## ---- override time-points for plotting
for (i in seq_along(ex)){
	t_ <- ex[[i]]$outputTimes
	ex[[i]]$measurementTimes <- t_
	mx_ <- max(t_)
	mn_ <- min(t_)
	if (abs(mx_ - mn_) < 1e-3) {
		mx_ <- mx_ + 1
		mn_ <- mn_ - 1
	}
	r <- abs(mx_ - mn_)
	L <- round(500.0 + 1000.0*r*r/(100.0**2+r*r))
	cat(sprintf("Experiment %i gets %i time points for plotting.\n",i,L))
	ex[[i]]$outputTimes <- sort(unique(c(t_,seq(mn_,mx_,length.out=L))))
}

## ---- Simulations:
simulate <- simulator.c(ex,modelName,parMap=uqsa::log10ParMap)
cat(sprintf("Simulating %i x %i trajectories, this will take a while.\n",length(ex),NROW(x)))
T0 <- Sys.time()
y <- simulate(t(x))

for (i in seq_along(y)){
	dimnames(y[[i]]$state) <- list(names(y0),NULL,NULL)
	dimnames(y[[i]]$func) <- list(names(ex[[i]]$outputValues),NULL,NULL)
}

print(difftime(Sys.time(),T0))
rm(x)

## ---- plot the function values and log10-state-variables:
pngFile <- sprintf("%s-sample-simulations-lines.png",PREFIX)
png(file=pngFile,width=length(ex)*1920,height=(NROW(y[[1]]$state)+NROW(y[[1]]$func))*1080,res=100)
print(ggplotTimeSeriesStates(y,ex,type="lines",MLE=MLE))
dev.off()
