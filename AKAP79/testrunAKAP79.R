
testRunParallelExperiments <-function(n=1){
	## test simulation, with one parameter vector
	print(experiments[[1]][['input']])
	print(parVal)
	## this will be parallel in *experiments*
	n<-length(experiments)
	## mclapply will return a list, so we remove one unnecessary layer of lists
	out <- mclapply(1:n,function(i) {unlist(runModel(experiments[i], modelName, as.matrix(parVal), parMap))})
	return(out)
}
