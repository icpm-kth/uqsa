#' Splits a formula into a left and right side
#'
#' This function splits a reaction formulka apart into its parts,
#' removing whitespace on each side: `"A + 2 B <=> AB2"` will be split
#' into a list with two entries
#'
#' ```
#' list$reactants == c("A","2*B")
#' list$products == c("AB2")
#' ```
#'
#' @param string - reactionFormula
#' @return a named list with a forward component and a backward
#'     component, each entry contains a character vector
parse.formula <- function(reactionFormula){
	LR <- ftsplit(reactionFormula,"<=>")
	reactants <- ftsplit(LR[1],"+")
	products <- ftsplit(LR[2],"+")
	return(list(reactants=reactants,products=products))
}

#' find the coefficients in a formula
#'
#' A reaction formula has reactants and products, separated by <=>,
#' with reactants on the left and products on the right (by convention).
#' Each of those is a plus separated list of reacting compounds and
#' modifiers, with optional coefficients, e.g.: `A + 2 B <=> AB2`
#'
#' Once the formula is split into left and right side, this function
#' determines the coefficients. For the above example, this function
#' returns `c(1,2)` for the left side and 1 for the right side.
#'
#' @return coefficients, as a vector
#' @param chrv a character vector as returned by parse.formula
#' @examples lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.coefficients)
match.coefficients <- function(chrv){
	cf <- numeric(length(chrv))
	l <- grepl("^[0-9]+",chrv) # leading numbers exist
	cf[l] <- as.numeric(sub("^([0-9]+)\\b.*","\\1",chrv[l]))
	cf[!l] <- 1
	return(cf)
}

#' Find the variable names in a formula
#'
#' A reaction formula has reactants and products, separated by `<=>`,
#' with reactants on the left and products on the right (by convention).
#' Each of those is a plus separated list of reacting compounds and
#' modifiers, with optional coefficients, e.g.: `A + 2 B <=> AB2`
#'
#' Once the formula is split into left and right side, this function
#' determines the names. For the above example, this function
#' returns `c("A","B")` for the left side and `"AB2"` for the right side.
#'
#' @return coefficients, as a vector
#' @param chrv a character vector as returned by parse.formula
#' @examples
#' lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.names)
#' lapply(uqsa:::parse.formula("A + 2*B <=> AB2"),uqsa:::match.names)
match.names <- function(chrv){
	vn <- character(length(chrv))
	l <- grepl("^[0-9]+",chrv) # leading numbers exist
	vn[l] <- sub("\\b([[:alpha:]]\\w*)$","\\1",chrv[l])
	vn[!l] <- chrv[!l]
	return(vn)
}

#' Find forward and backward component in a reaction kinetic
#'
#' a reaction kinetic can be almost any function, and in general it is
#' not possible to tell apart which part of a kinetic law is which.
#'
#' But for mass action kinetics, and positive reaction rate
#' coefficients, the expressions mostly look like this:
#'
#' ```
#' kf*prod(reactants.concentration) - kb*prod(product.concentrations)
#' ```
#'
#' this functions splits at `-`, and if none is present, then the
#' reaction is assumed to be irreversible.
#'
#' In a more general setting (where the `-` split is wrong),
#' the splitting has to be done by hand or more complex rules.
#' @param reactionKinetic a string with the kinetic law for a reaction
#' @return a character vector of components named 'forward' and 'backward'
#' @examples
#'  uqsa:::parse.kinetic("kf*A*B-kb*C")
#'  uqsa:::parse.kinetic("kf*A*B")
#'  uqsa:::parse.kinetic("kf*A/(Km+A)")
parse.kinetic <- function(reactionKinetic){
	if (grepl("^[^-]*-[^-]*$",reactionKinetic)){
		rates<-ftsplit(reactionKinetic,"-")
	} else {
		rates<-c(reactionKinetic,"0")
	}
	names(rates)<-c("forward","backward")
	return(rates)
}


#' Attempt to find multiplicative reaction rate coefficients
#'
#' This function assumes that the kinetic Law is mass action
#' kinetics. This function helps in converting the units of an ODE
#' into units that work in stochastic simulations. Converting the
#' units of a general formula (Michaelis Menten, Hill kinetics, etc.)
#' is difficult and dubious in stochastic simulations.
#'
#' For these reasons, this functions assumes: `kf*A*B*[...]` with the
#' first word `kf` representing the reaction rate coefficient.
#'
#' Given an SBtab document, this function finds the value (if any) and
#' unit of this coefficient. The coefficient itself can be defined as
#' a fixed constant, a parameter, or an algebraic expression in the
#' document. The most important attribute here is the unit.
#'
#' @param kineticLaw a string with a mathematical formula
#' @param tab an SBtab document, as returned by `SBtabVFGEN::sbtab_from_tsv()`
#' @return a parameter with value, and unit as attribute
parameter.from.kinetic.law <- function(kineticLaw,tab){
	if (kineticLaw == "0") return(0)
	cat("kineticLaw: ", kineticLaw,"\n")
	kName <- ftsplit(gsub("[-+*/()]"," ",kineticLaw)," ")[1]
	cat("kName: ",kName,"\n")
	if (kName %in% rownames(tab$Parameter)){
		id <- rownames(tab$Parameter)
		kValue <- tab$Parameter[["!DefaultValue"]][id==kName]
		kUnit <- as.character(tab$Parameter[["!Unit"]][id==kName])
	} else if ("Expression" %in% names(tab) && kName %in% rownames(tab$Expression)){
		id <- rownames(tab$Expression)
		kValue <- NA
		kUnit  <- as.character(tab$Expression[["!Unit"]][id==kName])
	} else if ("Constant" %in% names(tab) && kName %in% rownames(tab$Constant)){
		id <- rownames(tab$Constant)
		kValue <- tab$Constant[["!Value"]][id==kName]
		kUnit  <- as.character(tab$Constant[["!Unit"]][id==kName])
	} else {
		kValue <- 0.0
		kUnit <- "1"
	}
	k <- kValue
	attr(k,'unit') <- kUnit
	comment(k) <- kName
	return(k)
}

#' propensity creates a propensity formula
#'
#' given the custom math expressions needed to calculate a propensity,
#' the propensity coefficient and the kinetic law of the reaction,
#' this function makes a string that can be used with GillespieSSA2.
#'
#' The propensity coefficient translates between
#'
#' @param conv.coeff propensity conversion coefficient: `conv.coeff*kinetic.law = propensity function`
#' @param kinetic.law the kinetic law of this reaction (as used with ODEs)
#' @param rExpressions named math expressions that appear in the kinetic.law of this reaction
#' @return a string representation of the propensity function
propensity <- function(conv.coeff,kinetic.law,rExpressions){
	if (!is.null(rExpressions)){
		ex <- paste0(sprintf("%s = %s; ",names(rExpressions),rExpressions),collapse='')
	}
	p <- sprintf("%s %.10g*%s",ex,conv.coeff,kinetic.law)
	cat(sprintf("propensity: «%s»\n",p))
	return(p)
}


#' Functions to construct and run the stochastic simulation using GillespieSSA2 package
#'
#' This translates the Reaction network into the specific form required by GillespieSSA2
#'
#' @param model the model, represented by a list of data.frames with SBtab content
#' @param compile boolean value: if TRUE (default value), the GillespieSSA2 reactions are compiled (recommended for faster simulations)
#' @return a list of GillespieSSA2::reaction items (if compile==TRUE, the reactions are compiled)
#' @export
importReactionsSSA <- function(model.tab, compile = TRUE){
  num_reactions <- length(row.names(model.tab$Reaction))
  num_reversible_reactions <- sum(model.tab$Reaction[["!IsReversible"]]==TRUE)
  reactions <- vector("list", len=num_reactions + num_reversible_reactions)
  compound_names <- rownames(model.tab$Compound)
  k <- 1
  for(i in 1:num_reactions){
    kinetic_law <- model.tab$Reaction[["!KineticLaw"]][i]
    propensity <- sub("-.*", "", kinetic_law)
    effect <- c()
    formula <- model.tab$Reaction[["!ReactionFormula"]][i]
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
    reactions[[k]] <- GillespieSSA2::reaction(propensity = propensity, effect = effect)
    k <- k + 1
    if(model.tab$Reaction[["!IsReversible"]][i]){
      #also add the backward reaction
      propensity <- sub(".*-", "", kinetic_law)
      effect <- -effect
      if(sum(effect<0)==2){
        propensity <- paste0(propensity," / Phi")
      }
      reactions[[k]] <- GillespieSSA2::reaction(propensity = propensity, effect = effect)
      k <- k + 1
    }
  }

  if(k-1 != length(reactions)){
    error("Length of reactions list doesn't match")
  }
  if(compile){
    parVal <- model.tab$Parameter[["!DefaultValue"]]
    names(parVal) <- rownames(model.tab$Parameter)

    parameters_from_expressions <- parameters_from_expressions_func(model.tab)

    AvoNum <- 6.022e23          #Avogadro constant
    vol <- 4e-16  #volume where the reactions take place
    unit <- 1e-6  #unit of measure in meters (1e-6 corresponds to micrometer)
    Phi <- AvoNum * vol * unit #constant that will be used to simulate the random reactions

    reactions <- GillespieSSA2::compile_reactions(
      reactions = reactions,
      state_ids = rownames(model.tab$Compound),
      params = c(parVal, parameters_from_expressions(parVal), Phi=Phi)
    )
  }
  return(reactions)
}

#' Function that creates a function that reads the SBtab expression table and computes the expressions (which are functions of parameters) given the parameters passed in input
#'
#' @param model.tab an SBtab file imported in R, as it is returned by the function SBtabVFGEN::sbtab_from_tsv
#' @return a function that given a vector of (sampling) parameters returns a vector of parameters corresponding to the evaluations of the expressions in the SBtab expression table (which are functions of the parameters)
#' @export
parameters_from_expressions_func <- function(model.tab){
  param_from_expr <- function(parVal){
    # Return a named vector of parameters which are given by expressions
    if ("Input" %in% names(model.tab)){
     for(i in seq(NROW(model.tab$Input))){
      assign(rownames(model.tab$Input)[i],model.tab$Input[["!DefaultValue"]][i])
     }
    }
    if(is.null(parVal)){
      parVal <- model.tab$Parameter[["!DefaultValue"]]
    }

    parNames <- rownames(model.tab$Parameter)
    for(i in seq_along(parVal)){
      assign(parNames[i],parVal[i])
    }

    expr_formulas <- model.tab$Expression[["!Formula"]]
    expr <- c()

    for(i in seq_along(expr_formulas)){
      expr_value <- eval(str2expression(model.tab$Expression[["!Formula"]][i]))
      assign(rownames(model.tab$Expression)[i],expr_value)
      expr[i] <- expr_value
    }
    names(expr) <- rownames(model.tab$Expression)
    return(expr)
  }
  return(param_from_expr)
}



#' Function that creates a closure that simulates a stochastic trajectory with the Gillespie algorithm
#' given certain experimental conditions and a parameter vector, and computes the
#' distance between the simulation and the experimental data
#'
#' @param experiments a list of experiment.
#' @param model.tab an SBtab model.
#' @param reactions a list of GillespieSSA2 reactions.
#' @param parMap a function that translates ABC variables (parABC)
#'     into something the model will accept.
#' @param outputFunction a function that, given (t,state,param), outputs the quantity that is measured in the experiments (e.g., the outputFunction may output the concentration of one of the compounds in the system).
#' @param vol Volume in which the reactions take place.
#' @param unit unit of measure for the volume (it converts the volume into liters).
#' @param nStochSim number of stochastic simulations to average over.
#' @param distance (optional) function that computes the distance between experimental data and simulated data. Default value is NULL.
#' @return a function that, given a parameter, returns a simulated trajectory obtained via the Gillespie algorithm. Specifically, experimental time points, output value measured at these time points, and (if distance is not NULL) the distance between the simulated trajectories and the experimental data.
#' @export
simulator.stoch <- function(experiments, model.tab = model.tab, reactions = NULL, parMap = identity, outputFunction = function(t,state,param){state}, vol = 4e-16, unit = 1e-6, nStochSim = 3, distance = NULL){

    AvoNum <- 6.022e23          #Avogadro constant
    Phi <- AvoNum * vol * unit #constant that will be used to simulate the random reactions

    parameters_from_expressions <- parameters_from_expressions_func(model.tab)

    if(is.null(reactions)){
      reactions <- importReactionsSSA(model.tab, compile = TRUE)
    }

    simulate_one_experiment <- function(param, experiment){

      avgOutput <- rep(0, length(experiment[["outputTimes"]]))
      names(param) <- rownames(model.tab$Parameter)
      SSAparam <- c(parMap(param), Phi = Phi)
      if(!is.null(parameters_from_expressions)){
        SSAparam <- c(SSAparam, parameters_from_expressions(parMap(param)))
      }
      iv <- ceiling(experiment$initialState*Phi)
      names(iv) <- names(experiment$initialState)
      for(i in 1:nStochSim){
        out_ssa <- GillespieSSA2::ssa(
          initial_state = iv,
          reactions = reactions,
          params = SSAparam,
          final_time = max(experiment[["outputTimes"]]),
          method = ssa_exact(),
          verbose = FALSE,
          log_propensity = TRUE,
          log_firings = TRUE,
          census_interval = 5)

        # out$state is a matrix of dimension (time points)x(num compounds)
        output <- apply(out_ssa$state/Phi, 1, function(state){outputFunction(0, state, param=parMap(param))})
        if(sum(!is.na(out_ssa$time)) > 2){
          interpOutput <- approx(out_ssa$time, output, experiment[["outputTimes"]])
          interpOutput$y[experiment[["outputTimes"]]<min(out_ssa$time)] <- output[which.min(out_ssa$time)]
          interpOutput$y[is.na(interpOutput$y)] <- tail(output,1)
          avgOutput <- avgOutput + interpOutput$y
        } else {
          avgOutput <- Inf*avgOutput
        }
      }
      avgOutput <- avgOutput/nStochSim

      if(!is.null(distance)){
        dist <- distance(avgOutput,t(experiment[["outputValues"]]),t(experiment[["errorValues"]]))
      } else {
        dist <- NULL
      }

      return(list(time = experiment[["outputTimes"]], output = avgOutput, distance_data_simulation = dist))
    }

    sim <- function(param){
      simulators <- mclapply(experiments, function(e){simulate_one_experiment(param = param, experiment = e)})
      return(simulators)
    }
  return(sim)
}



#' Function that creates the objective function
#'
#' Given a parameter set, this function computes the distance between
#' experimental data and simulated data (coresponding to the parameter
#' in input).
#'
#' @param experiments a list of experiments.
#' @param model.tab an SBtab model.
#' @param parNames the names of the (biological) parameters of the
#'     model.
#' @param distance a user supplied function that calculates a distance
#'     between simulation and data with an interface of
#'     distance(simulation, data, errVal), where errVal is an estimate
#'     of the measurement noise (e.g. standard deviation), if needed
#'     by the function.
#' @param parMap a function that translates ABC variables (parABC)
#'     into something the model will accept.
#' @param outputFunction a function that, given (t,state,param), outputs the quantity that is measured in the experiments (e.g., the outputFunction may output the concentration of one of the compounds in the system).
#' @param vol Volume in which the reactions take place.
#' @param unit unit of measure for the volume (it converts the volume into liters).
#' @param reactions a list of GillespieSSA2 reactions
#' @param nStochSim number of stochastic simulations to average over.
#' @return a closure for the objective function that implicitly
#'     depends on all of the arguments to this function but explicitly
#'     only on the ABC parameters parABC.
#' @export
makeObjectiveSSA <- function(experiments, model.tab, parNames, distance, parMap=identity, outputFunction = identity, vol = 4e-16, unit = 1e-6, reactions, nStochSim = 1){

  simulate_experiments <- simulator.stoch(experiments = experiments,
                                          model.tab = model.tab,
                                          reactions = reactions,
                                          parMap = parMap,
                                          outputFunction = outputFunction,
                                          vol = vol,
                                          unit = unit,
                                          nStochSim = nStochSim,
                                          distance = distance)
  n_experiments <- length(experiments)
  objectiveFunction <- function(parABC){

    sim_all_experiments_one_parameter <- function(par){
      simulations <- simulate_experiments(par)
      return(sapply(1:n_experiments, function(i) simulations[[i]]$distance_data_simulation))
    }

    if (is.matrix(parABC)) {
      rownames(parABC) <- parNames
      npc <- ncol(parABC)
      S <- do.call(cbind,mclapply(1:npc, function(i) sim_all_experiments_one_parameter(parABC[,i])))
    }
    else {
      names(parABC) <- parNames
      S <- sim_all_experiments_one_parameter(parABC)
    }
    return(S)
  }
  return(objectiveFunction)
}


