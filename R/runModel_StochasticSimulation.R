ftsplit <- function(str,s){
	return(trimws(unlist(strsplit(str,s,fixed=TRUE))))
}

parse.formula <- function(reactionFormula){
	LR <- ftsplit(reactionFormula,"<=>")
	reactants <- ftsplit(LR[1],"+")
	products <- ftsplit(LR[2],"+")
	return(list(reactants=reactants,products=products))
}

parse.kinetic <- function(reactionKinetic){
	if (grepl("^[^-]*-[^-]*$",reactionKinetic)){
		rates<-ftsplit(reactionKinetic,"-")
	} else {
		error("The kinetic formula «%s» should match the pattern: 'A - B', where A is taken the forward rate and B the backward rate. Otherwise it's hard to determine the two.",reactionKinetic)
	}
	return(rates)
}

makeGillespieModel <- function(SBtab){
	stopifnot("Reaction" %in% names(SBtab))
	stopifnot("Compound" %in% names(SBtab))
	dR <- dim(SBtab$Reaction)
	compoundNames <- SBtab[["Compound"]][["!Name"]]
	if ("Expression" %in% names(SBtab)){
		expressionNames <- SBtab[["Expression"]][["!Name"]]
		expressionFormula <- SBtab[["Expression"]][["!Formula"]]
		names(expressionFormula)<-expressionNames
	} else {
		expressionNames <- NULL
		expressionFormula <- NULL
	}
	for (i in 1:dR[1]){
		reactionFormula <- model$Reaction[["!ReactionFormula"]][i]
		reactionKinetic <- model$Reaction[["!ReactionFormula"]][i]
		r <- parse.formula(reactionFormula)
		rVarNames <- unique(c(r$reactants,r$products))
		rCompoundNames <- rVarNames[rVarNames %in% compoundNames]
		rExpressions <- rVarNames[rVarNames %in% expressionNames]
		effect <- numeric(length(rCompoundNames))
		names(effect) <- rCompoundNames
		effect[r$re] <- effect[r$re]-1
		effect[r$pr] <- effect[r$pr]+1
		ktc <- parse.kinetic(reactionKinetic)
		if (!is.null(rExpressions)){
			propensity<-paste0(sprintf("%s = %s;",rExpressions,expressionFormula[rExpressions]),ktc[1],collapse=" ")
		} else {
			propensity<-ktc[1]
		}
	}
}


#' Functions to construct and run the stochastic simulation using GillespieSSA2 package
#'
#' This translates the Reaction network into the specific form required by GillespieSSA2
#'
#' @param model the model, represented by a list of data.frames with SBtab content
#' @return a list of GillespieSSA2::reaction items
#' @export
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


#' Function that creates the objective function
#'
#' Given a parameter set, this function computes the distance between
#' experimental data and simulated data (coresponding to the parameter
#' in input).
#'
#' @param experiments a list of experiments
#' @param parNames the names of the (biological) parameters of the
#'     model
#' @param distance a user supplied function that calculates a distance
#'     between simulation and data with an interface of
#'     distance(simulation, data, errVal), where errVal is an estimate
#'     of the measuremnet noise (e.g. standard deviation), if needed
#'     by the function.
#' @param parMap a function that translates ABC variables (parABC)
#'     into something the model will accept.
#' @param Phi Volume
#' @param reactions a list that encodes the reactions for
#'     GillespieSSA2
#' @param mc.cores same as for parallel::mclapply()
#' @param nStochSim number of stochastic simulations to average over
#' @return a closure for the objective function that implicitly
#'     depends on all of the arguments to this function but explicitly
#'     only on the ABC parameters parABC.
#' @export
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

