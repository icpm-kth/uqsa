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
runModel <- function(y0, modelName, params_inputs, outputTimes_list, outputFunctions_list, environment="R", mc.cores = 8){
  if(environment=="C")
  {
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",modelName)
    if (!file.exists(so)){
      system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,modelName,LIBS))
    }
    N <- dim(params_inputs)[2]
    outputTimes <- outputTimes_list[[1]] ##TO CHANGE WHEN WE WILL HAVE A MORE GENERAL r_gsl_odeiv2 THAT ACCEPTS DIFFERENT OUTPUTTIMES FOR EACH PARAMETER/INPUT SET
    yy_gsl<-r_gsl_odeiv2(modelName, as.double(outputTimes), y0, params_inputs)

    yy_gsl_as_list <- mclapply(seq(dim(yy_gsl)[3]), function(x) yy_gsl[ , , x], mc.preschedule = FALSE, mc.cores = mc.cores)
    output_yy <- mclapply(1:N, function(i)  apply(yy_gsl_as_list[[i]],2,outputFunctions_list[[i]]), mc.preschedule = FALSE, mc.cores = mc.cores)
  }
  else
  {
    if(environment!="R")
    {
      str <- sprintf("%s is not admitted as value for the variable 'environment'. The default value 'R' will be used instead.", environment)
      warning(str)
    }
    source(sprintf("%s.R",modelName))
    func <- eval(as.name(modelName))
    numSimulations <- dim(params_inputs)[2]
    yy <- mclapply(1:numSimulations, function(i) matrix(t(lsode(y0[,i], c(0,outputTimes_list[[i]]), func=func, parms=params_inputs[,i])[-1, -1]),ncol=length(outputTimes_list[[i]])), mc.preschedule = FALSE, mc.cores = mc.cores)
    output_yy <- mclapply(1:numSimulations, function(i)  apply(yy[[i]],2,outputFunctions_list[[i]]), mc.preschedule = FALSE, mc.cores = mc.cores)
  }
  return(output_yy)
}


