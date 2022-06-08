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
#' @param y0 Initial values: y'=f(t,y,p), y(t=0)=y0
#' @param modelName used to find model files and functions within the
#'     file (a prefix)
#' @param params_inputs a matrix of column vectors; each column
#'     contains a vector of both normal parameters (e.g. kinetic
#'     params like kf and kr) and input_parameters (concatenated in
#'     that order). With N columns, N simulations will be performed.
#' @param outputTimes_list a list of output time vectors, one per
#'     simulation experiment.
#' @param outputFunctions_list a list of vector valued output functions
#' @param environment "C" selects GSL solvers, "R" (default) selects deSolve as backend
#' @param mc.cores number of cores to use (defaults to 8)
#' @return output function values
runModel <- function(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, mc.cores = 8){
	if (is.null(comment(modelName))) {
		modelFile <- sprintf("%s.R",modelName))
	} else {
		modelFile <- comment(modelName)
	}
	if (grepl('.so$',modelFile)){
		so <- modelFile
		N <- dim(params_inputs)[2]
		outputTimes <- outputTimes_list[[1]] ##TO CHANGE WHEN WE WILL HAVE A MORE GENERAL r_gsl_odeiv2 THAT ACCEPTS DIFFERENT OUTPUTTIMES FOR EACH PARAMETER/INPUT SET
		yy_gsl<-r_gsl_odeiv2(modelName, as.double(outputTimes), y0, params_inputs)
		yy_gsl_as_list <- mclapply(seq(dim(yy_gsl)[3]), function(x) yy_gsl[ , , x], mc.preschedule = FALSE, mc.cores = mc.cores)
		output_yy <- mclapply(1:N, function(i)	apply(yy_gsl_as_list[[i]],2,outputFunctions_list[[i]]), mc.preschedule = FALSE, mc.cores = mc.cores)
	} else if (grepl('.[Rr]$',modelFile)) {
		stopifnot(file.exists(modelFile))
		source(modelFile)
		func <- eval(as.name(modelName))
		numSimulations <- dim(params_inputs)[2]
		yy <- mclapply(1:numSimulations, function(i) matrix(t(lsode(y0[,i], c(0,outputTimes_list[[i]]), func=func, parms=params_inputs[,i])[-1, -1]),ncol=length(outputTimes_list[[i]])), mc.preschedule = FALSE, mc.cores = mc.cores)
		output_yy <- mclapply(1:numSimulations, function(i)	apply(yy[[i]],2,outputFunctions_list[[i]]), mc.preschedule = FALSE, mc.cores = mc.cores)
	}
	return(output_yy)
}

#' checkModel tries to establish the simulation file for a given model
#'
#' This function returns the model name, with some additional comments about the file
#'
#' As an alternative to this function, it is sufficient to write
#' modelName <- "test_ode_model"             # or some other model name
#' comment(modelName) <- "test_ode_model.so" # or .R
#'
#' But, this function will compile a c source into a sharde object if modelFile ends with `.c`.
#' It also checks whether the model file exists.
#'
#' @export
#' @param modelName a string
#' @param modelFile a string, if the model file is different from "modelName.R"
#' @return modelName with an additional comment about which file to use for simulations
checkModel <- function(modelName,modelFile=NULL){
	if (modelFile==NULL) {
		modelFile <- paste0(modelName,'.R');
		stopifnot(file.exists(modelFile))
		message(sprintf('will use the file %s %s for simulations (deSolve backaned)',getwd(),modelFile))
		comment(modelName) <- modelFile
	} else if (grepl('.c$',modelFile)){
		message('building a shared library from c source, and using GSL odeiv2 as backend.')
		LIBS <- "-lgsl -lgslcblas -lm"
		CFLAGS <- "-shared -fPIC -Wall -O2"
		so <- sprintf("%s.so",modelName)
		system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,modelName,LIBS))
		stopifnot(file.exists(so))
		comment(modelName)<-so
	} else if (grepl('.so$',modelFile)) {
		stopifnot(file.exists(modelFile))
		message(sprintf('Will use pre-existing %s for simulations.',modelFile))
		comment(modelName) <- modelFile
	}
	return(modelName)
}
