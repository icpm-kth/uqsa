#!/usr/bin/env Rscript
library(uqsa)
library(SBtabVFGEN)
library(ggplot2)

showBestHalf <- TRUE

files <- commandArgs(trailingOnly=TRUE)
modelName <- "CaMKIIs"
comment(modelName) <- "./CaMKIIs.so"
sz <- 1e4
cat(sprintf("Loading sample of size: %i, this can take a few seconds.\n",sz))
print(files)

PREFIX <- uqsa::determinePrefix(files)
print(PREFIX)

x <- gatherSample(files,beta=1.0,size=sz)
cat("dim(x): ",paste(dim(x)),"\n")
l <- attr(x,"logLikelihood")
I <- order(l,decreasing=TRUE)[seq(round(length(l)/2))]
MLE <- I[1]
if (showBestHalf) {
	x <- x[I,]
	l <- l[I]
	attr(x,"logLikelihood") <- l
	MLE <- 1
}
sb <- sbtab_from_tsv(dir("..",pattern="tsv$",full.names=TRUE))
ex <- sbtab.data(sb)

## ---- R functions for this model:
source(sprintf("../R/%s.R",modelName))
t0 <- -1
p0 <- CaMKIIs_default(t0)
y0 <- CaMKIIs_init(t0,p0)
print(names(y0))
C <- NCOL(x)
colnames(x) <- names(p0[seq(C)])

## ---- Simulations:
simulate <- simcf(ex,modelName,parMap=uqsa::log10ParMap)
cat(sprintf("Simulating %i x %i trajectories, this will take a while.\n",length(ex),sz))
T0 <- Sys.time()
y <- simulate(t(x))
for (i in seq_along(y)){
	dimnames(y[[i]]$state) <- list(names(y0),NULL,NULL)
	dimnames(y[[i]]$func) <- list(names(ex[[i]]$outputValues),NULL,NULL)
}
print(difftime(Sys.time(),T0))

## ---- plot the function values and state-variables:
T <- "boxes"
png(file=sprintf("%s-sample-simulations-%s.png",PREFIX,T),width=length(ex)*1920,height=(NROW(y[[1]]$state)+NROW(y[[1]]$func))*1080,res=150)
print(ggplotTimeSeriesStates(y,ex,type=T,ttf=seq_along,yl.func=yl.func,yl.state=yl.state,MLE=MLE))
dev.off()
