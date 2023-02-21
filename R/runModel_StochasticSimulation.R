# Functions to construct and run the stochastic simulation using GillespieSSA2 package

importReactionsSSA <- function(model){
  
  num_reactions <- length(model$Reaction[["!ID"]])
  num_reversible_reactions <- sum(model$Reaction[["!IsReversible"]]==TRUE)
  reactions <- vector("list", len=num_reactions + num_reversible_reactions)
  compound_names <- model$Compound[["!Name"]]
  k <- 1
  for(i in 1:num_reactions){
    kinetic_law <- model$Reaction[["!KineticLaw"]][i]
    propensity <- sub("-.*", "", kinetic_law)
    effect <- c()
    formula <- model$Reaction[["!ReactionFormula"]][i]
    reactants <- sub("<=>.*", "", formula)
    products <- sub(".*<=>", "", formula)
    for(j in 1:length(compound_names)){
      if(length(grep(paste0("\\b",compound_names[j],"\\b"),products))){
        mult <- +1
        effect <- c(effect, mult)
        names(effect)[length(effect)] <- compound_names[j]
      }
      if(length(grep(paste0("\\b",compound_names[j],"\\b"),reactants))){
        mult <- -1
        effect <- c(effect, mult)
        names(effect)[length(effect)] <- compound_names[j]
      }
    }
    # HOW TO RESCALE THE REACTION RATES INTO STOCHASTIC REACTION RATES:
    # A + B -> C => c = k/Phi
    # A -> B => c = k
    if(sum(effect<0)==2){
      propensity <- paste0(propensity," / Phi")
    }
    reactions[[k]] <- GillespieSSA2::reaction(propensity = propensity, effect = effect) #, name = model$Reaction[["!Name"]][i])
    k <- k + 1
    if(model$Reaction[["!IsReversible"]][i]){
      #also add the backward reaction 
      propensity <- sub(".*-", "", kinetic_law)
      effect <- -effect
      if(sum(effect<0)==2){
        propensity <- paste0(propensity," / Phi")
      }
      reactions[[k]] <- GillespieSSA2::reaction(propensity = propensity, effect = effect) #, name = paste0(model$Reaction[["!Name"]][i],"_backward"))
      k <- k + 1
    }
  }
  
  if(k-1 != length(reactions)){
    error("Length of reactions list doesn't match")
  }
  return(reactions)
}


#possible input: compiled_reactions, Phi
#Can maybe remove: modelName
makeObjectiveSSA <- function(experiments, parNames, distance, parMap=identity, Phi, reactions, mc.cores=detectCores(), nStochSim = 1){

  
  objectiveFunction <- function(parABC){
    
    simulateAndComputeDistance <- function(e, param){
      
      avgOutput <- rep(0, length(e[["outputTimes"]]))
      for(i in 1:nStochSim){
        out_ssa <- GillespieSSA2::ssa(
          initial_state = ceil(e[["initialState"]]*Phi),
          reactions = reactions,
          params = c(parMap(param), Phi=Phi),
          final_time = max(e[["outputTimes"]]),
          method = ssa_exact(), 
          verbose = FALSE,
          log_propensity = TRUE,
          log_firings = TRUE,
          census_interval = 0.001,
          sim_name = modelName)
        
        # out$state is a matrix of dimension (time points)x(num compounds)
        output <- apply(out_ssa$state, 1, e[["outputFunction"]])
        interpOutput <- approx(out_ssa$time, output, e[["outputTimes"]])
        interpOutput$y[is.na(interpOutput$y)] <- tail(output,1)
        avgOutput <- avgOutput + interpOutput$y
      }
      avgOutput <- avgOutput/nStochSim
      return(distance(avgOutput/Phi,e[["outputValues"]],e[["errorValues"]]))
    }
    
    if (is.matrix(parABC)) {
      rownames(parABC) <- parNames
      npc <- ncol(parABC)
      S <- mclapply(1:npc, function(i) sapply(experiments, function(e) simulateAndComputeDistance(e, parABC[,i])), mc.cores = mc.cores)
      return(unlist(S))
    }
    else {
      names(parABC) <- parNames
      S <- mclapply(experiments, function(e) simulateAndComputeDistance(e, parABC), mc.cores = mc.cores)
      return(unlist(S))
    }
  }
  return(objectiveFunction)
}

