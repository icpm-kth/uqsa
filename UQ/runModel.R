remotes::install_github("a-kramer/rgsl", ref="OpenMP")

runModel <- function(y0, modelFunctionName, params_inputs, outputTimes_list, outputFunctions_list, environment="R", mc.cores = 8){
    
  if(environment=="C")
  {
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",modelFunctionName)
    if (!file.exists(so)){
      system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,modelFunctionName,LIBS))
    }
    
    if (require("rgsl") && require("parallel")){
      N <- dim(params_inputs)[2]
      outputTimes <- outputTimes_list[[1]] ##TO CHANGE WHEN WE WILL HAVE A MORE GENERAL r_gsl_odeiv2 THAT ACCEPTS DIFFERENT OUTPUTTIMES FOR EACH PARAMETER/INPUT SET
      yy_gsl<-r_gsl_odeiv2(modelFunctionName, as.double(outputTimes), y0, params_inputs)
      
      yy_gsl_as_list <- mclapply(seq(dim(yy_gsl)[3]), function(x) yy_gsl[ , , x], mc.preschedule = FALSE, mc.cores = mc.cores)
      output_yy <- mclapply(1:N, function(i)  apply(yy_gsl_as_list[[i]],2,outputFunctions_list[[i]]), mc.preschedule = FALSE, mc.cores = mc.cores)
    }
  }
  else
  {
    if(environment!="R")
    {
      str <- sprintf("%s is not admitted as value for the variable 'environment'. The default value 'R' will be used instead.", environment)
      warning(str)
    }
    source(sprintf("%s.R",modelFunctionName,modelFunctionName))
    func <- eval(as.name(modelFunctionName))
    if(require("deSolve") && require("parallel")){
      
      numSimulations <- dim(params_inputs)[2]
      yy <- mclapply(1:numSimulations, function(i) matrix(t(lsode(y0[,i], c(0,outputTimes_list[[i]]), func=func, parms=params_inputs[,i])[-1, -1]),ncol=length(outputTimes_list[[i]])), mc.preschedule = FALSE, mc.cores = mc.cores)
      output_yy <- mclapply(1:numSimulations, function(i)  apply(yy[[i]],2,outputFunctions_list[[i]]), mc.preschedule = FALSE, mc.cores = mc.cores)
    }
  }
  return(output_yy)
}


