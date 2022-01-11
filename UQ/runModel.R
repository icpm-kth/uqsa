remotes::install_github("a-kramer/rgsl", ref="OpenMP")

runModel <- function(y0, modelFunctionName, params, inputs, outputTimes, outputFunction, environment="R", mc.cores = 8){

  if(environment=="C")
  {
    
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",modelFunctionName)
    if (!file.exists(so)){
      system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,modelFunctionName,LIBS))
    }
    
    n_inputs <- length(inputs)
    param_mat <- matrix(params,nrow=length(params),ncol=n_inputs)
    input_mat <- do.call(cbind,inputs)
    p <- rbind(param_mat,input_mat)
    
    if (require("rgsl") && require("parallel")){
      yy_gsl<-r_gsl_odeiv2(modelFunctionName,outputTimes,y0,p)
      yy_gsl_as_list <- mclapply(seq(dim(yy_gsl)[3]), function(x) yy_gsl[ , , x], mc.preschedule = FALSE, mc.cores = mc.cores)
      output_yy <- mclapply(yy_gsl_as_list, outputFunction, mc.preschedule = FALSE, mc.cores = mc.cores)
    }
  }
  else
  {
    if(environment!="R")
    {
      str <- sprintf("%s is not admitted as value for the variable 'environment'. The default value 'R' will be used instead.",environment)
      warning(str)
    }
    source(sprintf("%s.R",modelFunctionName))
    func <- eval(as.name(modelFunctionName))
    if(require("deSolve") && require("parallel")){
      yy <- mclapply(inputs, function(inp) matrix(t(lsode(y0, outputTimes, func=func, parms=c(params,inp))[, -1]),ncol=length(outputTimes)), mc.preschedule = FALSE, mc.cores = mc.cores)
      output_yy <- mclapply(yy, outputFunction, mc.preschedule = FALSE, mc.cores = mc.cores)
    }
  }
  #return(matrix(unlist(yy),ncol=length(yy)))
  return(output_yy)
}


