# Simulate a model imported from SBtab and plot simulations and experimental data, given a parameter vector in input
# Copyright (C) 2023 Federica Milinanni (fedmil@kth.se)

#' Simulate and plot Data and Simulation
#'
#' This function imports model and experimental data saved as SBtab files,
#' it simualtes the model with the same initial conditions and input as in
#' the corresponding experiments, and it plots the simulations together with
#' the corresponding experimental data.
#' It currently works with one dimensional output.
#'
#' @export
#' @param SBtabDir the directory that contains `.tsv` files (with SBtab content)
#' @param paramVal parameter vector with which the model has to be simulated
#' @param plotDir the directory where the plots will be saved as .pdf and .RData variables
#' @param width width of the plot window
#' @param heigth heigth of the plot window
#' @return a vector of R plot objects
plotSimualtionsFromSBtab <- function(SBtabDir, paramVal, plotDir = NULL, width = 15, heigth = 10){
  require(SBtabVFGEN)
  require(rgsl)
  require(uqsa)
  require(deSolve)
  require(reshape2)
  require(ggplot2)

  model = import_from_SBtab(SBtabDir)
  #setwd(SBtabDir)
  rFile<-paste0(SBtabDir,'/',comment(model),'.R')
  cFile<-paste0(SBtabDir,'/',comment(model),'_gvf.c')
  modelName <- checkModel(comment(model),cFile) # uncomment to use gsl (C) as backend
  #modelName <- checkModel(comment(model),rFile)  # uncomment to use deSolve as backend
  #source(rFile)
  experiments <- import_experiments(modelName, SBtabDir)

  paramNames <- model[["Parameter"]][["!Name"]]
  names(paramVal)<-paramNames

  # This function maps parameter variables from logarithmic to linear scale
  parMap <- function(parABC) {
    return(10^parABC)
  }

  # Simualte the model
  output_yy <- runModel(experiments, modelName, paramVal, parMap = parMap)

  # Create folder for plots (if not provided in input)
  if(is.null(plotDir)){
    plotDir <- paste0(SBtabDir,"/../plots",modelName)
    if(!dir.exists(plotDir))
      dir.create(plotDir)
  }

  plots <- list()
  for(i in 1:length(experiments)){

    # Save the output of the simulation as a data frame with time and output value
    df <- as.data.frame(x = list(output_yy[[i]][["func"]][1,,1], experiments[[1]][["outputTimes"]]), col.names = c("y","t"))

    # Save also the experimental data as a dataframe
    yy_exp <- experiments[[i]][["outputValues"]]
    dfExp <- data.frame(t = experiments[[i]][["outputTimes"]],y = yy_exp)

    pl <- ggplot(df, aes(x = t, y = y)) + geom_line(color = "blue") + geom_point(data = na.omit(dfExp), aes(x = t, y = y)) + ggtitle(model$Experiments[["!Name"]])
    print(pl)
    ggsave(paste0(plotDir,"/",modelName,"_PlotExperimentAndSimulation_",model$Experiments[["!Name"]][i],".pdf"), plot = pl, width = width, height = heigth, units = "cm")

    plots[[i]] <- pl

    # lengthState <- dim(output_yy[[i]][["state"]])[1]
    # nrows <- ceil(sqrt(lengthState))

    # gtable()
    # par(mfrow=c(nrows,nrows))
    # print(pl)
    # for(j in 1:lengthState){
    #   df <- as.data.frame(x = list(output_yy[[i]][["state"]][j,,1], experiments[[1]][["outputTimes"]]), col.names = c("y","t"))
    #   plState <- ggplot(df, aes(x = t, y = y)) + geom_line(color = "blue") + ggtitle(model$Compound[["!Name"]][j])
    #   print(plState)
    # }
  }

  save(plots, file = paste0(plotDir,"/",modelName,"_PlotsExperimentsAndSimulations.RData"))

  return(plots)
}



# Example
# paramVal <- c( 5.090291e-01,  1.484507e-01, -2.615819e+00, -2.317419e-01, -7.612798e-01,
# -3.306151e+00, -2.324251e+00, -3.556083e+00, -2.872133e+00, -7.930016e-05,
# 1.103114e+00, -3.955432e-01, -7.729774e-01, -4.968397e-01,  7.941385e-01,
# 1.180435e+00, -4.999689e+00, -6.042806e-01, -1.305900e+00, -4.115583e-01,
# 2.965258e+00, -1.669173e+00, -9.119168e-01,  1.018882e+00,  1.933725e+00,
# -6.643307e-02,  2.061177e-02)
#
# SBtabDir <- # path to the directory with the SBtab files
#
# plots <- plotSimualtionsFromSBtab(SBtabDir, paramVal)
#
