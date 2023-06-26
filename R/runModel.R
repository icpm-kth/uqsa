# Uncertainty Quantification: run model ordinary differential equation
# Copyright (C) 2022 Federica Milinanni (fedmil@kth.se)

#' Simulate an Experiment using the named ODE Model
#'
#' Simulation experiments consist at least of initial values for the
#' state variables, a parameter vector, and a list of times at which
#' the solution needs to be known.
#'
#' This function will use the GSL solvers, or deSolve [default].  In
#' addition, a model usually has observables: values that depend on
#' the state variables and can be measured in a real experiment. These
#' are modeled by output functions.
#'
#' We distinguish normal parameters and input parameters. Input
#' parameters are known and not subject to any estimation
#' procedure. Furthermore, they are meant to represent the
#' experimental conditions, so they are either under direct control of
#' the experimenter or very carefully measured. The inputs are
#' probably different for each simulation experiment in at least one
#' value.
#'
#' If the environment variable is set to "C", then this function will
#' attempt to compile the file modelName_gvf.c to a shared library
#' modelName.so, if it doesn't already exist.
#'
#' @export
#' @param experiments list of experiments to simulate
#' @param modelName used to find model files and functions within the
#'     file (a prefix)
#' @param parABC a matrix of column vectors; each column
#'     contains a vector of both normal parameters (e.g. kinetic
#'     params like kf and kr) and input_parameters (concatenated in
#'     that order). With N columns, N simulations will be performed.
#' @return output function values
runModel <- function(experiments, modelName,  parABC, parMap=identity){
	if (is.matrix(parABC)){
		npc <- ncol(parABC)
	} else {
		npc <- 1
	}
	numExperiments <- length(experiments)
	# transform the ABC parameters if necessary (parMap is a user supplied function),
	# then determine the effective number of parameters, after transformation
	modelPar <- parMap(parABC)
	if (is.matrix(modelPar)) {
		np <- nrow(modelPar)
		N <- ncol(modelPar)
	} else {
		np <- length(modelPar)
		N <- 1
	}
	N <- N*numExperiments

	# an experiment can have an optional input
	if ('input' %in% names(experiments[[1]])){
		nu <- length(experiments[[1]][['input']])
	} else {
		nu <- 0
	}

	if (is.null(comment(modelName))) {
		modelFile <- sprintf("%s.R",modelName)
		stopifnot(file.exists(modelFile))
	} else {
		modelFile <- comment(modelName)
	}

	if (grepl('.so$',modelFile,useBytes=TRUE) && requireNamespace("rgsl")){
		so <- modelFile
		output_yy <- rgsl::r_gsl_odeiv2_outer(modelName, experiments, matrix(modelPar,np))
	} else if (grepl('.[Rr]$', modelFile) && requireNamespace("deSolve")) {
		## densely repeat the model parameters npc times
		modelPar <- matrix(modelPar,np,npc*numExperiments)

		## create a matrix that has all experimental inputs, repeated (densely) npc times per experiment
		if (nu>0){
			V <- vapply(experiments,function(E) matrix(E[['input']],nu,npc),FUN.VALUE=matrix(0,nu,npc))
			dim(V) <- c(nu,npc*numExperiments)
			modelPar <- rbind(modelPar,V)
		}
		stopifnot(ncol(modelPar)==N)
		## create a matrix of initial states, repeated npc times per experiment, as with the inputs and parameters
		ny <- length(experiments[[1]][['initialState']])
		y0 <- vapply(experiments, function(E) matrix(E[['initialState']],ny,npc), FUN.VALUE=matrix(0,ny,npc))
		dim(y0) <- c(ny,npc*numExperiments)

		outputTimes_list <- list()
		outputFunctions_list <- list()
		for(i in 1:numExperiments){
			outputTimes_list <- c(outputTimes_list, replicate(npc, list(experiments[[i]][["outputTimes"]])))
			outputFunctions_list <- c(outputFunctions_list, replicate(npc, list(experiments[[i]][["outputFunction"]])))
		}

		stopifnot(file.exists(modelFile))
		source(modelFile)
		func <- eval(as.name(paste0(modelName,"_vf")))
		yy <- mclapply(1:N, function(i) matrix(t(deSolve::lsode(y0[,i], c(0,outputTimes_list[[i]]), func=func, parms=modelPar[,i])[-1, -1]), ncol=length(outputTimes_list[[i]])))

		out_yy <- mclapply(1:N, function(i) apply(yy[[i]],2,outputFunctions_list[[i]]))

		fun <- function(yy_, out_yy_){
		  out_state <- do.call(cbind, yy_)
		  dim(out_state) <- c(dim(yy_[[1]]), length(yy_))

		  out_func <- do.call(cbind, out_yy_)
		  dim(out_func) <- c(length(out_func)/(dim(out_state)[2]*npc), dim(out_state)[2], npc)
		  return(list(state = out_state, func = out_func))
		}
		output_yy <- mclapply(1:numExperiments, function(i) fun(yy[1:npc + npc*(i-1)], out_yy[1:npc + npc*(i-1)]))
	} else {
		stop("one of the two solvers must be installed: rgsl, deSolve")
	}
	return(output_yy)
}

simulator.c <- function(experiments, modelName,  parABC, parMap=identity){
	require(rgsl)
## simulator:
	sim <- function(parABC){
		modelPar <- parMap(parABC)
		yf <- rgsl::r_gsl_odeiv2_outer(modelName, experiments, as.matrix(modelPar))
		return(yf)
	}
	return(sim)
}

simulate.R  <- function(experiments, model,  parABC, parMap=identity){
	if (is.matrix(parABC)){
		npc <- ncol(parABC)
	} else {
		npc <- 1
	}
	numExperiments <- length(experiments)
	# transform the ABC parameters if necessary (parMap is a user supplied function),
	# then determine the effective number of parameters, after transformation

	# an experiment can have an optional input
	if ('input' %in% names(experiments[[1]])){
		nu <- length(experiments[[1]][['input']])
	} else {
		nu <- 0
	}
	sim <- function(parABC){
		modelPar <- parMap(parABC)
		if (is.matrix(modelPar)) {
			np <- nrow(modelPar)
			N <- ncol(modelPar)
		} else {
			np <- length(modelPar)
			N <- 1
		}
		N <- N*numExperiments
		modelPar <- matrix(modelPar,np,npc*numExperiments)

		## create a matrix that has all experimental inputs, repeated (densely) npc times per experiment
		if (nu>0){
			V <- vapply(experiments,function(E) matrix(E[['input']],nu,npc),FUN.VALUE=matrix(0,nu,npc))
			dim(V) <- c(nu,npc*numExperiments)
			modelPar <- rbind(modelPar,V)
		}
		stopifnot(ncol(modelPar)==N)
		## create a matrix of initial states, repeated npc times per experiment, as with the inputs and parameters
		ny <- length(experiments[[1]][['initialState']])
		y0 <- vapply(experiments, function(E) matrix(E[['initialState']],ny,npc), FUN.VALUE=matrix(0,ny,npc))
		dim(y0) <- c(ny,npc*numExperiments)

		outputTimes_list <- list()
		outputFunctions_list <- list()
		for(i in 1:numExperiments){
			outputTimes_list <- c(outputTimes_list, replicate(npc, list(experiments[[i]][["outputTimes"]])))
			outputFunctions_list <- c(outputFunctions_list, replicate(npc, list(experiments[[i]][["outputFunction"]])))
		}
		func <- eval(as.name(paste0(modelName,"_vf")))
		yy <- mclapply(1:N, function(i) matrix(t(deSolve::lsode(y0[,i], c(0,outputTimes_list[[i]]), func=model$vf, parms=modelPar[,i])[-1, -1]), ncol=length(outputTimes_list[[i]])))

		out_yy <- mclapply(1:N, function(i) apply(yy[[i]],2,outputFunctions_list[[i]]))

		fun <- function(yy_, out_yy_){
			out_state <- do.call(cbind, yy_)
			dim(out_state) <- c(dim(yy_[[1]]), length(yy_))
			out_func <- do.call(cbind, out_yy_)
			dim(out_func) <- c(length(out_func)/(dim(out_state)[2]*npc), dim(out_state)[2], npc)
			return(list(state = out_state, func = out_func))
		}
		output_yy <- mclapply(1:numExperiments, function(i) fun(yy[1:npc + npc*(i-1)], out_yy[1:npc + npc*(i-1)]))
		return(output_yy)
	}
	return(sim)
}

#' checkModel tries to establish the simulation file for a given model
#'
#' This function returns the model name, with some additional comments
#' about the file
#'
#' As an alternative to this function, it is sufficient to write
#' modelName <- "test_ode_model"             # or some other model name
#' comment(modelName) <- "test_ode_model.so" # or .R
#'
#' This function will not attempt to find a model file, other than in
#' the current directory. But, checkModel will compile a GSL
#' compatible C source file into a shared object if modelFile ends
#' with `.c` and stop if that doesn't work.
#'
#' In any case, this function stops execution if the model file
#' doesn't exist.
#'
#' @export
#' @param modelName a string
#' @param modelFile a string, if the model file is different from
#'     "modelName.R". If the file name ends in .c, the c source will be
#'     compiled to a shared library.
#' @return modelName with an additional comment about which file to use for simulations
checkModel <- function(modelName,modelFile=NULL){
	if (is.null(modelFile)) {
		modelFile <- paste0(modelName,c('.R','_gvf.c','.so'));
		modelFile <- modelFile[file.exists(modelFile)]
		stopifnot(!is.null(modelFile) && length(modelFile)>0)
		modelFile <- modelFile[1]
	}
	if (grepl('[.]c$',modelFile,useBytes=TRUE)){
		stopifnot(file.exists(modelFile))
		message('building a shared library from c source, and using GSL odeiv2 as backend (pkg-config is used here).')
		LIBS <- "`pkg-config --libs gsl`"
		CFLAGS <- "-shared -fPIC `pkg-config --cflags gsl`"
		so <- sprintf("./%s.so",modelName)
		command_args <- sprintf("%s -o '%s' '%s' %s",CFLAGS,so,modelFile,LIBS)
		message(paste("cc",command_args))
		system2("cc",command_args)
		stopifnot(file.exists(so))
		comment(modelName)<-so
	} else if (grepl('[.]so$',modelFile,useBytes=TRUE)) {
		stopifnot(file.exists(modelFile))
		message(sprintf('Will use pre-existing %s for simulations.',modelFile))
		comment(modelName) <- modelFile
	} else if (grepl('[.]R$',modelFile,useBytes=TRUE)){
		message(sprintf('will use the %s for simulations (deSolve backend)',modelFile))
		comment(modelName) <- modelFile
	}
	return(modelName)
}

#' creates Objective functions from ingredients
#'
#' the returned objective function has only one argument (the ABC
#' parameters)
#'
#' @export
#' @param experiments a list of simulation experiments
#' @param modelName and model storage file as comment
#' @param getScore a function that calculates ABC scores
#' @param parMap a function that transforms ABC variables into acceptable model parameters
#' @return an objective function
makeObjective <- function(experiments,modelName,distance,parMap=identity,simulate=NULL)
{
	if (is.null(simulate)){
		simulate <- function(par){
			return(runModel(experiments, modelName,  par, parMap))
		}
	}
	Objective <- function(parABC){
		out <- simulate(parABC)
		S <- c()
		for(i in 1:length(experiments)){
		  S <- c(S, unlist(mclapply(1:dim(out[[i]]$func)[3], function(j) distance(out[[i]]$func[,,j], experiments[[i]]$outputValues, experiments[[i]]$errorValues))))
		}
		return(S)
	}
	return(Objective)
}

#' Ths function Creates an acceptanceProbability function
#'
#' The returned closure needs only the sampling variables (parABC) as
#' inputand calculates a probability of accepting Markov chainmoves.
#'
#' @export
#' @param experiments a list of experiments
#' @param modelName an annotated string, with the model name and model file as comment
#' @param getAcceptanceProbability an R function that mape the results of a simulation and experimental data to an acceptance probability
#' @param parMap an optional mapping between sampling parameters (parABC) and model parameters (e.g. rescaling,re-ordering).
#' @return a function that calculates probabilities given only parABC as input; it implicitly uses all the argiments to this function.
makeAcceptanceProbability <- function(experiments, modelName, getAcceptanceProbability, parMap=identity)
{
	acceptanceProbability <- function(parABC){
		out <- runModel(experiments, modelName,  parABC, parMap)
		S <- c()
		for(i in 1:length(experiments)){
		  S <- c(S, unlist(mclapply(1:dim(out[[i]]$func)[3], function(j) getAcceptanceProbability(out[[i]]$func[,,j], experiments[[i]]$outputValues, experiments[[i]]$errorValues))))
		}
		return(min(S))
	}
	return(acceptanceProbability)
}
