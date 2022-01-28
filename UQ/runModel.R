remotes::install_github("a-kramer/rgsl") #, ref="OpenMP")

runModel <- function(y0, modelFunctionName, params_inputs, outputTimes, outputFunction, environment="R", mc.cores = 8){
    
  if(environment=="C")
  {
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",modelFunctionName)
    if (!file.exists(so)){
      system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,modelFunctionName,LIBS))
    }
    
    if (require("rgsl") && require("parallel")){
      N <- ncol(params_inputs)
      y0 <- matrix(rep(y0,N),ncol=N)
      yy_gsl<-r_gsl_odeiv2(modelFunctionName,outputTimes,y0,params_inputs)
      yy_gsl_as_list <- mclapply(seq(dim(yy_gsl)[3]), function(x) yy_gsl[ , , x], mc.preschedule = FALSE, mc.cores = mc.cores)
      output_yy <- mclapply(yy_gsl_as_list, function(yy_) apply(yy_,2,outputFunction), mc.preschedule = FALSE, mc.cores = mc.cores)
    }
  }
  else
  {
    if(environment!="R")
    {
      str <- sprintf("%s is not admitted as value for the variable 'environment'. The default value 'R' will be used instead.", environment)
      warning(str)
    }
    source(sprintf("%s.R",modelFunctionName))
    func <- eval(as.name(modelFunctionName))
    if(require("deSolve") && require("parallel")){
      
      params_inputs_as_list <- mclapply(seq(dim(params_inputs)[2]), function(x) params_inputs[,x], mc.preschedule = FALSE, mc.cores = mc.cores)
      yy <- mclapply(params_inputs_as_list, function(p) matrix(t(lsode(y0, outputTimes, func=func, parms=p)[, -1]),ncol=length(outputTimes)), mc.preschedule = FALSE, mc.cores = mc.cores)
      output_yy <- mclapply(yy, function(yy_) apply(yy_,2,outputFunction), mc.preschedule = FALSE, mc.cores = mc.cores)
    }
  }
  return(output_yy)
}


