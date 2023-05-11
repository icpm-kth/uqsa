#' Splits a formula into a left and right side
#'
#' This function splits a reaction formulka apart into its parts,
#' removing whitespace on each side: "A + 2*B <=> AB2" will be split
#' into a list with two entries
#' list$reactants == c("A","2*B")
#' list$products == c("AB2")
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
#' modifiers, with optional coefficients, e.g.: A + 2*B <=> AB2
#'
#' Once the formula is split into left and right side, this function
#' determines the coefficients. For the above example, this function
#' returns c(1,2) for the left side and 1 for the right side.
#'
#' @return coefficients, as a vector
#' @param chrv a character vector as returned by parse.formula
#' @examples lapply(parse.formula("A + 2*B <=> AB2"),match.coefficients)
#' $reactants
#' [1] 1 2
#' $products
#' [1] 1
match.coefficients <- function(chrv){
	cf <- numeric(length(chrv))
	l <- grepl("^[0-9]+",chrv) # leading numbers exist
	cf[l] <- as.numeric(sub("^([0-9]+)\\b","\\1",chrv[l]))
	cf[!l] <- 1
	return(cf)
}

#' Find the variable names in a formula
#'
#' A reaction formula has reactants and products, separated by <=>,
#' with reactants on the left and products on the right (by convention).
#' Each of those is a plus separated list of reacting compounds and
#' modifiers, with optional coefficients, e.g.: A + 2*B <=> AB2
#'
#' Once the formula is split into left and right side, this function
#' determines the names. For the above example, this function
#' returns c("A","B") for the left side and "AB2" for the right side.
#'
#' @return coefficients, as a vector
#' @param chrv a character vector as returned by parse.formula
#' @examples lapply(parse.formula("A + 2*B <=> AB2"),match.names)
#' lapply(parse.formula("A + 2*B <=> AB2"),match.names)
#' $reactants
#' [1] "A" "B"
#' $products
#' [1] "AB2"
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
#' kf*prod(reactants.concentration) - kb*prod(product.concentrations)
#'
#' this functions splits at '-', and if none is present, then the
#' reaction is assumed to be irreversible.
#'
#' In a more general setting (where the '-' split is wrong),
#' the splitting has to be done by hand or more complex rules.
#' @param reactionKinetic a string with the kinetic law for a reaction
#' @return a character vector of components named 'forward' and 'backward'
#' @examples
#' > parse.kinetic("kf*A*B-kb*C")
#'  forward backward
#' "kf*A*B"   "kb*C"
#'
#' > parse.kinetic("kf*A*B")
#'  forward backward
#' "kf*A*B"      "0"
#'
#' > parse.kinetic("kf*A/(Km+A)")
#'       forward      backward
#' "kf*A/(Km+A)"           "0"
parse.kinetic <- function(reactionKinetic){
	if (grepl("^[^-]*-[^-]*$",reactionKinetic)){
		rates<-ftsplit(reactionKinetic,"-")
	} else {
		rates<-c(reactionKinetic,"0")
	}
	names(rates)<-c("forward","backward")
	return(rates)
}

#' convert ODE parameter to Gillespie parameter
#'
#' ODE parameters usually have a different unit of measurement than
#' the parameters we need for stochastic simulators.  ODEs have
#' fluxes, which are in multiples of M/s (M is mol/liter), same unit
#' as the first derivative of the state variables.
#'
#' The reaction rate coefficients of mass action kinetics, kf and kb
#' have units that are compatible with the flux units, depending on
#' the order of the reaction (the order is related to the reaction's
#' stoichiometry).
#'
#' @param k the ODE reaction rate coefficient (mandatory)
#' @param n multiplicity of each reactant, if any (order > 0); omit for zero-order
#' @param LV L*V -- product of _Avogadro's number_ and _volume_ [defaults to 6.02214076e+8]
#' @return rescaled parameter for stochastic simulation with a comment of how to re-scale it
#' @examples reaction: "2 A + B -> C"
#' k <- 1.0
#' attr(k,'unit') <- "µM/s"
#' n <- c(2,1)
#' reactants <- c('A','B')
#' convert.parameter(k,n)
convert.parameter <- function(k, n=0, LV=6.02214076e+8){
	order <- sum(n)
	if (!is.null(attr(k,'unit'))){
		unit <- SBtabVFGEN::unit.from.string(attr(k,'unit'))
		unit.conversion <- 10^sum(unit$scale*unit$exponent)
	} else {
		unit.conversion <- 1
	}
	## zero-order reaction c = LVk
	## first-order reaction c = k
	## etc.
	## switch is 1-based, not 0-based:
	## switch(0,'a','b') is nothing,
	## switch(1,'a','b') is 'a'
	order.conversion <- switch(order+1,LV,1.0,1.0/LV)
	c <- unit.conversion*order.conversion
	if (order==2 && length(n)==1) c<-2*c
	propensity.coefficient <- k*c
	attr(propensity.coefficient,'conversion') <- c
	return(propensity.coefficient)
}

parameter.from.kinetic.law <- function(kineticLaw,parValue,parNames,parUnits){
	i <- parNames %in% ftsplit(kineticLaw,'*')
	if (any(i)){
		k<-parValue[i]
		comment(k)<-parNames[i]
		attr(k,'unit')<-parUnits[i]
	} else {
		k <- 0
		comment(k) <- ''
		attr(k,'unit') <- '1'
	}
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
#' @param conv.coeff propensity conversion coefficient: conv.coeff*kinetic.law = propensity function
#' @param kinetic.law the kinetic law of this reaction (as used with ODEs)
#' @param rExpressions named math expressions that appear in the kinetic.law of this reaction
#' @return a string representation of the propensity function
propensity <- function(conv.coeff,kinetic.law,rExpressions){
	ex <- ''
	if (!is.null(rExpressions)){
		ex <- paste0(sprintf("%s = %s; ",names(rExpressions),rExpressions),collapse='')
	}
	p <- sprintf("%s %.10g*%s",ex,conv.coeff,kinetic.law)
	return(p)
}

#' Create a list of reactions for GillespieSSA2
#'
#' This function takes a series of SBtab tables, as returned by
#' SBtabVFGEN::sbtab_from_tsv() and creates GillespieSSA2 reactions
#' from them. Reactions arfe made pairwise, as forward and backward
#' reaction pairs. If a backward reaction doesn't exist, the list item
#' is NULL. A valid set of reactions can be obtained with `!is.null(reactions)`
#'
#' @export
#' @import GillespieSSA2
#' @param SBtab a series of tables as returned by `sbtab_from_tsv()`
#' @param LV is the product of Avogadro's constant L and the system's
#'     volume V in litres; if unspecified this information is retrieved
#'     from the SBtab files, if missing we assume 1µm³ of volume (the
#'     approximate sizes of bacteria or synapses)
#' @return a list of reactions
#' @examples
#'  model.tsv<-dir(pattern="[.]tsv$")
#'  model.sbtab<-SBtabVFGEN::sbtab_from_tsv(model.tsv)
#'  reactions<-makeGillespieModel(model.sbtab)
#'  l<-is.null(reactions)
#'  model.ssa2<-reactions[!l]
makeGillespieModel <- function(SBtab,LV=NULL){
	stopifnot("Reaction" %in% names(SBtab))
	stopifnot("Compound" %in% names(SBtab))
	dR <- dim(SBtab$Reaction)
	n <- dR[1]
	if (is.null(LV) && 'Compartment' %in% names(SBtab) && "!Size" %in% names(SBtab[["Compartment"]])){
		L <- 6.02214076e+23
		V <- SBtab[["Compartment"]][["!Size"]][1] # litres?
		LV <- L*V
	} else {
		## L <- 6.02214076e+23
		## V <- 1e-15 # litres
		LV <- 6.02214076e+8 # L*V
	}
	parValue <- SBtab[["Parameter"]][["!DefaultValue"]]
	parNames <- SBtab[["Parameter"]][["!Name"]]
	parUnits <- SBtab[["Parameter"]][["!Unit"]]
	compoundNames <- SBtab[["Compound"]][["!Name"]]
	reactionNames <- SBtab[["Reaction"]][["!Name"]]
	isReversible <- as.logical(SBtab[["Reaction"]][["!IsReversible"]])
	SSA2reactions <- vector(mode="list",length=2*n)
	if ("Expression" %in% names(SBtab)){
		expressionNames <- SBtab[["Expression"]][["!Name"]]
		expressionFormula <- SBtab[["Expression"]][["!Formula"]]
		names(expressionFormula)<-expressionNames
	} else {
		expressionNames <- NULL
		expressionFormula <- NULL
	}
	for (i in 1:dR[1]){
		reactionFormula <- SBtab$Reaction[["!ReactionFormula"]][i]
		reactionKinetic <- SBtab$Reaction[["!KineticLaw"]][i]
		r <- parse.formula(reactionFormula)
		rVarNames <- unique(c(match.names(r$reactants),match.names(r$products)))
		rCompoundNames <- rVarNames[rVarNames %in% compoundNames]
		rExpressionNames <- rVarNames[rVarNames %in% expressionNames]
		rExpressions <- expressionFormula[rExpressionNames]
		names(rExpressions) <- rExpressionNames
		effect <- numeric(length(rCompoundNames))
		names(effect) <- rCompoundNames
		cf<-list(reactants=match.coefficients(r$re),products=match.coefficients(r$pr))
		effect[r$re] <- effect[r$re]-cf$re
		effect[r$pr] <- effect[r$pr]+cf$pr
		ktc <- parse.kinetic(reactionKinetic)
		kf <- parameter.from.kinetic.law(ktc[1],parValue,parNames,parUnits)
		kb <- parameter.from.kinetic.law(ktc[2],parValue,parNames,parUnits)
		pf <- convert.parameter(kf,cf$re,LV=LV)
		pb <- convert.parameter(kb,cf$pr,LV=LV)
		plain.reaction.name <- sprintf("%s_%s",sub("[^a-zA-Z0-9_]","_",reactionNames[i]),c('forward','backward'))
		propFormula <- list(f=propensity(attr(pf,'conversion'),ktc[1],rExpressions),b=propensity(attr(pb,'conversion'),ktc[2],rExpressions))
		SSA2reactions[[2*(i-1)+1]] <-  GillespieSSA2::reaction(propFormula$f,effect,plain.reaction.name[1])
		if (isReversible[i]) {
			SSA2reactions[[2*(i-1)+2]] <-  GillespieSSA2::reaction(propFormula$b,-1*effect,plain.reaction.name[2])
		}
	}
	return(SSA2reactions)
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

