plotDataAndSimulations<- function(draws,experiment,parMap=identity,nCores=detectCores()){
nDraws = dim(draws)[1]
stopifnot(length(experiment)==1)
output_yy <- runModel(experiment, modelName, t(draws), parMap, nCores)
df_ <- mclapply(1:dim(output_yy[[1]][["func"]])[3], function(i) as.data.frame(x = list(output_yy[[1]][["func"]][1,,i]/0.2,experiment[[1]][['outputTimes']]), col.names = c("y","t")))

df__ <- melt(df_,id=c("t","y"))
yy_exp <- (experiment[[1]][["outputValues"]]-108.6)/(183.9-108.6)
#yy_exp <- (experiments[[expInd]][["outputValues"]]-100)/(171.67-100)
dfExpA<- data.frame(t=experiment[[1]][["outputTimes"]], y=yy_exp)
ggplot(df__,aes(x=t, y=y, group=L1))+
  geom_line(color="blue")+
  geom_point(data=dfExpA, aes(x=t, y=y), inherit.aes=FALSE)
}
