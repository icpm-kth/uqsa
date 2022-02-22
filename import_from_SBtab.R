remotes::install_github("a-kramer/SBtabVFGEN")

import_from_SBtab <- function(SBtabDir){
  tsvList <- dir(pattern = ".*[.]tsv$")
  sbtab_model <- SBtabVFGEN::sbtab_from_tsv(tsvList)
  
  if(length(dir(pattern = ".*[.]vf$")) == 0){
    SBtabVFGEN::sbtab_to_vfgen(sbtab_model, cla = FALSE)		#this function creates a .vf file (located in the tsvDirectory) from the SBtab model
    vfFileName <- dir(pattern = ".*[.]vf$")
    system(paste("cd", SBtabDir))
    system(paste("vfgen r:func=yes", vfFileName))
    system(paste("vfgen gsl", vfFileName))
  }
  #Create an R function that given a state, returns a vector of outputs
  vfFileName <- dir(pattern = ".*[.]vf$")
  modelName <- substr(vfFileName,1,nchar(vfFileName)-3)
  
  system(paste("../many_outputs_to_one.sh ", modelName, ".R > ", modelName, "_out.R", sep=""))
  return(sbtab_model)
}


import_experiments <- function(modelName, SBtabDir){

  SBtab <- import_from_SBtab(SBtabDir)
  
  compoundNames <- SBtab[["Compound"]][["!Name"]]
  compoundId <- SBtab[["Compound"]][["!ID"]]
  default_y0 <- SBtab[["Compound"]][["!InitialValue"]]
  
  inputNames <- SBtab[["Input"]][["!Name"]]
  inputId <- SBtab[["Input"]][["!ID"]]
  defaultInputs <- SBtab[["Input"]][["!DefaultValue"]]
  
  outputNames <- SBtab[["Output"]][["!Name"]]
  outputId <- SBtab[["Output"]][["!ID"]]
  source(paste(SBtabDir, "/", modelName, "_out.R", sep = ""))
  vectorialOutputFunction <- eval(as.name(paste(modelName, "_out", sep = "")))
  
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
