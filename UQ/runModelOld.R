# Uncertainty Quantification: Model simulation
# Copyright (C) 2018 
# Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Olivia Eriksson
# Andrei Kramer
# Anu G Nair

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
remotes::install_github("a-kramer/rgsl", ref="OpenMP")
library(rgsl)
runModel <- function(y0, modelFunction, params, inputs, outputTimes, outputFunction, envir="R", mc.cores = 8){
  
  if(envir=="R")
  {
    yy <- mclapply(inputs, function(inp) matrix(t(lsode(y0, outputTimes, func=modelFunction, parms=c(params,inp))[, -1]),ncol=length(outputTimes)), mc.preschedule = FALSE, mc.cores = mc.cores)
    output_yy <- mclapply(yy, outputFunction, mc.preschedule = FALSE, mc.cores = mc.cores)
    
    ## MAKE THE OUTPUT INTO A LIST OF VECTORS (EACH ELEMENT OF THE LIST CORRESPONDS TO A DIFFERENT SET OF INPUTS)
    #output_yy <- do.call(cbind, output_yy)
    
    ## NON PARALLELIZED
    #yy <- lapply(inputs, function(inp) matrix(t(lsode(y0, outputTimes, func=modelFunction, parms=c(params,inp))[-1, -1]),ncol=length(outputTimes)-1))
    #output_yy <- lapply(yy)
  }
  else if(envir=="C")
  {
    n_inputs <- length(inputs)
    param_mat <- matrix(params,nrow=length(params),ncol=n_inputs)
    input_mat <- do.call(cbind,inputs)
    p <- rbind(param_mat,input_mat)
    LIBS <- "-lgsl -lgslcblas -lm"
    CFLAGS <- "-shared -fPIC -Wall -O2"
    so <- sprintf("%s.so",modelFunction)
    if (!file.exists(so)){
      system2("gcc",sprintf("%s -o %s %s_gvf.c %s",CFLAGS,so,modelFunction,LIBS))
    }
    yy_gsl<-r_gsl_odeiv2(modelFunction,outputTimes,y0,p)
    
    ## NON PARALLELIZED
    #output_yy_gsl <- mclapply(yy_gsl, outputFunction, mc.preschedule = FALSE, mc.cores = mc.cores)
    
    yy_gsl_as_list <- mclapply(seq(dim(yy_gsl)[3]), function(x) yy_gsl[ , , x], mc.preschedule = FALSE, mc.cores = mc.cores)
    output_yy_gsl <- mclapply(yy_gsl_as_list, outputFunction, mc.preschedule = FALSE, mc.cores = mc.cores)
  }
  #return(matrix(unlist(yy),ncol=length(yy)))
  return(yy_gsl)
}


runModelOld <- function(tpar, parIdx, input, rInd){
  
  parDefVal <- c(0.0021524, 0.0088968, 0.028311, 0.041046, 0.046, 0.046, 0.046, 
              0.046, 0.046, 0.0016667, 0.016667, 0.008125, 0.092857, 15836, 
              22959, 1124, 4646, 0.028, 70, 800, 60, 600, 0.021, 0.04, 0.02, 
              0.002, 0.0044, 0.021, 0.021, 0.021, 0.021, 52, 2000, 4000, 5000, 
              2272.7, 0.04, 0.011, 0.011, 0.011, 0.011, 0.02, 0.002, 0.011, 
              0.0044, 2000, 4000, 5000, 2272.7, 0.73, 50, 0.05, 10, 10)
  
  pars <- parDefVal
  pars[parIdx] <- 10^tpar
  yy <- numeric(length(input))
  
  if(rInd=='A'){
    yy <- SSsolutionFigA(pars, input) #param, Ca
  }else if(rInd=='B'){
    yy <- SSsolutionFigB(pars, input, 100) #param, CaM, PP2B
  }else if (rInd=='C1'){
    yy<- SSsolutionFigC(pars, input, 30,3) #param, Ca, totalCaM, totalPP2B
  }
  xx = log10(input)
  return(list(xx=xx, yy=yy))
}
