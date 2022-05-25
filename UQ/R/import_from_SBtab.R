#' Import Systems Biology Models given in Tabular Form 
#'
#' This function uses SBtabVFGEN to read a series of tsv files, each
#' containing one systems biology table that together create a model
#' (Reactions, Parameters, etc.). The model is converted to a vfgen
#' compatible file. This file is processed by vfgen through a system
#' call to create source code for R (deSolve) and C (GSL solvers).
#'
#' This requires vfgen to be installed
#' (https://github.com/WarrenWeckesser/vfgen) - not an R package.
#' 
#' SBtab is a particular convention on how to structure the tables
#' (sbtab.net)
#' 
#' @export
#' @param SBtabDir the directory that contains `.tsv` files (with SBtab content)
#' @return a model as a list of data.frames, one per tsv file, the
#'     Document title is attached as a comment attribute:
#'     comment(model) = Document Title
#' @examples
#'  model <- import_from_SBtab("./model")
#'  comment(model)
#'  source("model.R")
import_from_SBtab <- function(SBtabDir){
  tsvList <- dir(path = SBtabDir, pattern = ".*[.]tsv$")
  sbtab_model <- SBtabVFGEN::sbtab_from_tsv(paste(SBtabDir,"/",tsvList,sep=""))
  
  if(length(dir(path = SBtabDir, pattern = ".*[.]vf$")) == 0){
    SBtabVFGEN::sbtab_to_vfgen(sbtab_model, cla = FALSE)		#this function creates a .vf file (located in the tsvDirectory) from the SBtab model
    vfFileName <- dir(path = SBtabDir, pattern = ".*[.]vf$")
    system(paste("cd", SBtabDir))
    system(paste("vfgen r:func=yes", vfFileName))
    system(paste("vfgen gsl", vfFileName))
  }
  return(sbtab_model)
}

#' Reads the Data and Model Contained in an SBtab Document (tsv)
#'
#' An SBtab Document is a set of tables that represent reactions,
#' compounds, parameters, and measured data that correspond to
#' simulations of the model under certain input conditions and initial
#' values.
#'
#' This function assumes that this information is stored in a series
#' of tsv files. The content is imported using the SBtabVFGEN package.
#'
#' The data contents are reorganized into a list of simulation
#' experiments (initial values, measurement time points, etc.)
#'
#' @export
#' @param modelName (string) the functions of the model have this prefix
#' @param SBtabDir (string) a local directory that contains tsv files (with SBtab content)
#' @return list of simulation experiments (and the data corresponding to that simulation)
import_experiments <- function(modelName=NULL, SBtabDir){

  SBtab <- import_from_SBtab(SBtabDir)
  if (is.null(modelName)) {
		modelName <- comment(SBtab)
	}
  compoundNames <- SBtab[["Compound"]][["!Name"]]
  compoundId <- SBtab[["Compound"]][["!ID"]]
  default_y0 <- SBtab[["Compound"]][["!InitialValue"]]
  
  inputNames <- SBtab[["Input"]][["!Name"]]
  inputId <- SBtab[["Input"]][["!ID"]]
  defaultInputs <- SBtab[["Input"]][["!DefaultValue"]]
  
  outputNames <- SBtab[["Output"]][["!Name"]]
  outputId <- SBtab[["Output"]][["!ID"]]
  source(paste(SBtabDir, "/", modelName, ".R", sep = ""))
  vectorialOutputFunction <- eval(as.name(paste(modelName,"_", outputNames[1], sep="")))
  
  n_experiments <- dim(SBtab[["Experiments"]])[1]
  experiments <- vector("list", length = n_experiments)
  
  inputs_and_initState_ids <- colnames(SBtab[["Experiments"]])
  inputs_and_initState_ids <- inputs_and_initState_ids[startsWith(inputs_and_initState_ids, '>')]
  inputs_and_initState_ids <- substring(inputs_and_initState_ids, 2)
  match_input <-  match(inputs_and_initState_ids, inputId)
  match_initState <- match(inputs_and_initState_ids, compoundId)
  
  experiments_names <- SBtab[["Experiments"]][[1]]
  
  for(i in 1:n_experiments){
    inputs_and_initState_vals <- unlist(SBtab[["Experiments"]][i, paste(">",inputs_and_initState_ids,sep="")])
    
    experiments[[i]][["initialState"]] <- default_y0
    experiments[[i]][["initialState"]][match_initState[!is.na(match_initState)]] <- inputs_and_initState_vals[!is.na(match_initState)]
    
    experiments[[i]][["input"]] <- defaultInputs
    experiments[[i]][["input"]][match_input[!is.na(match_input)]]  <- inputs_and_initState_vals[!is.na(match_input)]
    
    experiment_table <- SBtab[[experiments_names[i]]]
    experiments[[i]][["outputTimes"]] <- experiment_table[["!Time"]]
    tabColnames <- colnames(experiment_table)
    tabOutputId <- tabColnames[startsWith(tabColnames,'>')]
    tabOutputId <- substring(tabOutputId, 2)
    match_output <- match(tabOutputId, outputId)
    match_output <- match_output[!is.na(match_output)]
    experiments[[i]][["outputNames"]] <- outputId[match_output]
    experiments[[i]][["outputId"]] <- outputId[match_output]
    experiments[[i]][["outputValues"]] <- as.matrix(experiment_table[paste(">", outputId[match_output], sep = "")])
    notNAidx <- apply(!is.na(experiments[[i]][["outputValues"]]) & experiments[[i]][["outputTimes"]]>0, 1, prod)
    experiments[[i]][["outputValues"]] <- experiments[[i]][["outputValues"]][as.logical(notNAidx),]
    experiments[[i]][["outputTimes"]] <- experiments[[i]][["outputTimes"]][as.logical(notNAidx)]
    experiments[[i]][["outputFunction"]] <- function(yy) {vectorialOutputFunction(0.0, yy, 0)[match_output]}
  }
  
  return(experiments)
}
