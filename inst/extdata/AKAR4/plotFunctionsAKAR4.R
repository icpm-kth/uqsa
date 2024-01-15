plotDataAndSimulations<- function(draws,experiment,parMap=identity,nCores=detectCores()){
nDraws = dim(draws)[1]
stopifnot(length(experiment)==1)
<<<<<<< HEAD

simulate <- simulator.c(experiment,modelName,parMap)
output_yy <- simulate(t(draws))

# old: when we used runModel
#output_yy <- runModel(experiment, modelName, t(draws), parMap)

=======
simulate <- simulator.c(experiment,modelName,parMap)
output_yy <- simulate(t(draws))
>>>>>>> bfe23a4858fefc378f942dfe5b801def1e063668
df_ <- mclapply(1:dim(output_yy[[1]][["func"]])[3], function(i) as.data.frame(x = list(output_yy[[1]][["func"]][1,,i],experiment[[1]][['outputTimes']]), col.names = c("y","t")))

df__ <- melt(df_,id=c("t","y"))
yy_exp <- experiment[[1]][["outputValues"]]
<<<<<<< HEAD
#yy_exp <- (experiments[[expInd]][["outputValues"]]-100)/(171.67-100)
=======
>>>>>>> bfe23a4858fefc378f942dfe5b801def1e063668
dfExpA<- data.frame(t=experiment[[1]][["outputTimes"]], y=yy_exp)
names(dfExpA)[2] <- "y"
ggplot(df__,aes(x=t, y=y, group=L1))+
  geom_line(color="blue")+
  geom_point(data=dfExpA, aes(x=t, y=y), inherit.aes=FALSE)
}

# an alternative function with simpler design, but doesn't work either
plotSample<-function(draws,experiment,parMap=identity,nCores=detectCores()){
	N<-dim(draws)[1]
<<<<<<< HEAD
	out <- runModel(experiment, modelName, t(draws), parMap)
=======
	simulate <- simulator.c(experiment,modelName,parMap)
	output_yy <- simulate(t(draws))
>>>>>>> bfe23a4858fefc378f942dfe5b801def1e063668
	plot.new()
	for(i in 1:N){
		lines(experiment[[1]]$outputTimes,out[[1]]$func[1,,i])
	}
	points(experiment[[1]]$outputTimes,experiment[[1]]$outputValues)
}
