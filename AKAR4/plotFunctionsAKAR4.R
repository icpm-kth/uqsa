plotDataAndSimulations<- function(draws, expInd){
nDraws = dim(draws)[1]

tmp_list <- mclapply(experiments[expInd], function(x) replicate(nDraws, c(parVal,x[["input"]])),  mc.preschedule = FALSE, mc.cores = nCores)
params_inputs <- do.call(cbind, tmp_list)
params_inputs[parIdx,] <- 10^t(draws)

tmp_list <- mclapply(experiments[expInd], function(x) replicate(nDraws, x[["initialState"]]),  mc.preschedule = FALSE, mc.cores = nCores)
y0 <- do.call(cbind, tmp_list)

outputTimes_list <- list()
outputFunctions_list <- list()
for(k in 1:length(expInd)){
  outputTimes_list <- c(outputTimes_list, replicate(nDraws, list(experiments[[k]][["outputTimes"]])))
  outputFunctions_list <- c(outputFunctions_list, replicate(nDraws, list(experiments[[k]][["outputFunction"]])))
}

output_yy <- runModel(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment, nCores)
df_ <- mclapply(1:length(output_yy), function(i) as.data.frame(x = list(output_yy[[i]]/0.2,outputTimes_list[[i]]), col.names = c("y","t")))

df__ <- melt(df_,id=c("t","y"))
yy_exp <- (experiments[[expInd]][["outputValues"]]-108.6)/(183.9-108.6)
#yy_exp <- (experiments[[expInd]][["outputValues"]]-100)/(171.67-100)
dfExpA<- data.frame(t=experiments[[expInd]][["outputTimes"]], y=yy_exp)
ggplot(df__,aes(x=t, y=y, group=L1))+
  geom_line(color="blue")+
  geom_point(data=dfExpA, aes(x=t, y=y), inherit.aes=FALSE)
}
