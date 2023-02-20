plotDataAndSimulations<- function(draws,experiment,parMap=identity,nCores=detectCores()){
nDraws = dim(draws)[1]
stopifnot(length(experiment)==1)
out <- runModel(experiment, modelName, t(draws), parMap, nCores)
output_yy <- out[[1]]$func
df_ <- mclapply(1:length(output_yy), function(i) as.data.frame(x = list(output_yy[i,]/0.2,experiment[[1]][['outputTimes']]), col.names = c("y","t")))

df__ <- melt(df_,id=c("t","y"))
yy_exp <- (experiment[[1]][["outputValues"]]-108.6)/(183.9-108.6)
#yy_exp <- (experiments[[expInd]][["outputValues"]]-100)/(171.67-100)
dfExpA<- data.frame(t=experiment[[1]][["outputTimes"]], y=yy_exp)
ggplot(df__,aes(x=t, y=y, group=L1))+
  geom_line(color="blue")+
  geom_point(data=dfExpA, aes(x=t, y=y), inherit.aes=FALSE)
}

# an alternative function with simpler design, but doesn't work either
plotSample<-function(draws,experiment,parMap=identity,nCores=detectCores()){
	N<-dim(draws)[1]
	out <- runModel(experiment, modelName, t(draws), parMap, nCores)
	plot.new()
	for(i in 1:N){
		lines(experiment[[1]]$outputTimes,out[[1]]$func[1,,i])
	}
	points(experiment[[1]]$outputTimes,experiment[[1]]$outputValues)
}
