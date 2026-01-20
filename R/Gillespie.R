#' Returns a list of reaction coefficients
#'
#' This function maps c("A","2 B") to c(1,2)
#'
#' This is a C function because it is much easier to write in C.  C
#' has the strtod() function which expects a leading number and stops
#' when the numbers end. as.character() returns NA if the input
#' contains any dirt.
#'
#' The reaction formula is as tring like this: "A + 2 B <=> C", when
#' split at `<=>` and then later at `+`, we get the strings that must
#' be parsed: "A" and "2 B" for the left side and "C" for the right
#' side. The numbers are the stoichiometric constants, or coefficients.
#'
#' @useDynLib uqsa, lstrtod
#' @param formulaList a list of character vectors, derived from the
#'     left or right side of a reaction formula:
#' @return a list of numeric coefficient vectors
#' @examples
#' as.character("12 A")
#' onlyCoefficients("12 A")
#' @export
onlyCoefficients <- function(formulaList){
	C <- lapply(formulaList,\(v) .Call(lstrtod,as.character(v)))
	return(C)
}

#' Returns only the names in a reaction formula
#'
#' This is the companion function to onlyCoefficients. It returns the
#' names of reactants, without the stoichiometry.
#'
#' @export
#' @param formulaList a list of strings like: "2 B" or "45 X"
#' @return a list of name vectors
onlyNames <- function(formulaList){
	return(lapply(formulaList,\(x) sub("^[0-9]+\\b","",x)))
}

#' This function returns a list of named stoichiometric vectors
#'
#' Given an already split list of entries such as c("3 A","B"), this
#' function returns a numeric vector c(3,1) with names c("A","B").
#' @export
#' @param formulaList a list of split reaction formula tokens
#' @return named numeric vector of stoichiometric coefficients
stoichiometry <- function(formulaList){
	sc <- onlyCoefficients(formulaList)
	nm <- onlyNames(formulaList)
	sc <- mapply(\(a,b) {names(a) <- trimws(b); return(a)}, sc, nm) # make all the vectors in the list named
	return(sc)
}

#' Returns information about parameter conversion
#'
#' Given parameters of a reaction kinetic coefficient model with
#' concentrations as state variables, this function returns the
#' parameters of the stochastic version of that model
#'
#' For this function, it is important that "2 A -> B" is not written
#' as "A + A -> B" (this may be fixed later).
#'
#' @param value parameter values (kinetic rate coefficients), a vector
#' @param unit of the kinetic parameter, a character vector
#' @param stoichiometry a list of the stoichiometric constants of each reaction, a list of integer vectors
#' @param LV Avogadro's constant multiplied with the system's volume (the default is a reasonable value for a small volume)
#' @return a vector of parameter conversion factors
#' @export
parameterConversion <- function(value, unit, stoichiometry){
	order <- unlist(lapply(stoichiometry,sum))
	l <- lapply(stoichiometry,length)
	u <- lapply(unit,SBtabVFGEN::unit.from.string) # data.frames
	# 2 A -> B reactions
	aa <- (l==1) & (order==2)
	# unit conversion factor
	f <- unlist(lapply(u,\(x) 10^sum(x$scale*x$exponent)))
	# taking reaction order into account:
	CF <- data.frame(order=order, lvpower=(1-order), factor=f*(1+aa),row.names=names(value))
	return(CF)
}

parameterReaction <- function(reactionKinetic,parNames){
	RC <- ifelse(
		grepl("-",reactionKinetic,fixed=TRUE),
		reactionKinetic,
		paste0(reactionKinetic," - 0")
	)
	RC <- trimws(unlist(strsplit(RC,"-",fixed=TRUE)))
	dim(RC) <- c(2,length(reactionKinetic))
	RC <- t(RC)
	# forward and backward reaction kinetics
	k_fwd <- lapply(parNames, \(x) grep(paste0("\\b",x,"\\b"),RC[,1]))
	names(k_fwd) <- parNames
	k_bwd <- lapply(parNames, \(x) grep(paste0("\\b",x,"\\b"),RC[,2]))
	names(k_bwd) <- parNames
	return(list(fwd=unlist(k_fwd),bwd=unlist(k_bwd)))
}

initialCount <- function(sb){
	v <- sb$Compound[["!InitialValue"]]
	names(v) <- rownames(sb$Compound)
	unit <- sb$Compound[["!Unit"]]
	u <- lapply(unit,SBtabVFGEN::unit.from.string)
	f <- unlist(lapply(u,\(x) 10^sum(x$scale*x$exponent)))
	return(v*f)
}

#' Concentration to Count Conversion
#'
#' This function returns the factor f for converting a concentration
#' to particle numbers (counts) based on the units specified in the
#' SBtab file and the volume of the system encoded as LV (Avogadro's
#' constant × Volume).
#'
#' The implicit assumption is that the LV constant the user passes to
#' this function has a compatible unit of volume to the SBtab unit of
#' concentration. LV should be in litres and the concentrations
#' also in litres up to SI prefixes. All of these are OK: M, nmol/l,
#' mol/dl, kilomol/megalitre.
#'
#' If the concentrations are 'something per cubic metre', then LV has to be in m³
#' as well. The unit of LV is not checked in any way. But the unit of
#' concentration will be parsed an converted to mol/l.
#'
#' @param sb SBtab list of data.frames with one of the members called
#'     'Compound', with a mandatory "!Unit" column.
#' @param LV Avogadro's constant multiplied by system volume
#' @return f a conversion factor: n <- f*x, where n is a particle
#'     count and x a concentration in the specified volume.
#' @export
ccc <- function(sb){
	u <- lapply(sb$Compound[["!Unit"]],SBtabVFGEN::unit.from.string)
	f <- unlist(lapply(u,\(x) 10^sum(x$scale*x$exponent)))
	names(f) <- rownames(sb$Compound)
	return(f)
}

simpleConversion <- function(v,unit){
	u <- lapply(unit,SBtabVFGEN::unit.from.string)
	x <- unlist(lapply(u,\(u) u$exponent[pmatch("mole",u$kind)]))
	x[is.na(x)] <- 0.0
	f <- unlist(lapply(u,\(u) 10^sum(u$scale*u$exponent)))
	cf <- data.frame(factor=f, lvpower=x)
	rownames(cf) <- names(v)
	return(cf)
}


#' makeGillespieModel interprets the provided SBtab file as a stochastic model
#'
#' The SBtab file is assumed to describe a reaction network.  The
#' systems biology information is assumed to be concentrations and
#' rate coefficients.
#'
#' With the information provided with the rate coefficient units and a
#' volume, this function tries to convert everything to Gillespie rate
#' constants.
#' @param sb list of data.frames
#' @param LV Avogadro's constant L multiplied by the system's volume V.
#' @return a list containing the interpreted model.
#' @export
makeGillespieModel <- function(sb){
	f <- sb$Reaction[["!ReactionFormula"]]
	rev <- as.logical(sb$Reaction[["!IsReversible"]])
	F <- sub("*"," ",trimws(unlist(strsplit(f,"<=>",fixed=TRUE))),fixed=TRUE)
	dim(F) <- c(2,length(f))
	F <- t(F)
	cl <- stoichiometry(strsplit(F[,1],"+",fixed=TRUE))
	names(cl) <- paste0(rownames(sb$Reaction),"_fwd")
	cr <- stoichiometry(strsplit(F[,2],"+",fixed=TRUE))
	names(cr) <- paste0(rownames(sb$Reaction),"_bwd")

	k <- c(sb$Parameter[["!DefaultValue"]],sb$Input[["!DefaultValue"]])
	names(k) <- c(rownames(sb$Parameter),rownames(sb$Input))
	u <- c(sb$Parameter[["!Unit"]],sb$Input[["!Unit"]])
	names(u) <- c(rownames(sb$Parameter),rownames(sb$Input))

	pr <- parameterReaction(sb$Reaction[["!KineticLaw"]],rownames(sb$Parameter))
	convFactorFwd <- parameterConversion(k[names(pr$fwd)], u[names(pr$fwd)], cl[pr$fwd])
	convFactorBwd <- parameterConversion(k[names(pr$bwd)], u[names(pr$bwd)], cr[pr$bwd])
	l <- !(names(k) %in% rownames(convFactorFwd) | names(k) %in% rownames(convFactorBwd))
	scv <- simpleConversion(k[l],u[l])
	cUnit <- sb$Compound[["!Unit"]]
	names(cUnit) <- rownames(sb$Compound)
	parUnit <- c(sb$Parameter[["!Unit"]],sb$Input[["!Unit"]])
	names(parUnit) <- c(rownames(sb$Parameter),rownames(sb$Input))
	return(list(left=cl,right=cr,scf=convFactorFwd,scb=convFactorBwd,scv=scv,initialCount=initialCount(sb),par=k,compoundUnit=cUnit, parameterUnit=parUnit))
}

reactionEffect <- function(sm){
	return(
		unlist(
			c(
				mapply(
				\(x,y,n) c(sprintf("\tcase _%s:",n),sprintf("\t\tx[_%s] -= %i;",names(x),x),sprintf("\t\tx[_%s] += %i;",names(y),y),sprintf("\t\tbreak;")),
				sm$left,
				sm$right,
				names(sm$left)
				),
				mapply(
				\(x,y,n) c(sprintf("\tcase _%s:",n),sprintf("\t\tx[_%s] -= %i;",names(x),x),sprintf("\t\tx[_%s] += %i;",names(y),y),sprintf("\t\tbreak;")),
				sm$right,
				sm$left,
				names(sm$right)
				)
			)
		)
	)
}

lengths <- function(v){
	return(sapply(as.character(v),length))
}

padding <- function(v,upper.bound=40){
	return(upper.bound - lengths(v))
}

#' Returns a string that contains code to scale the parameter
#'
#' The goal is to convertthe kinetic parameter k to the stochastic
#' parameter c.
#'
#' The input is a pre-calculated factor f, determined from the
#' parameter's unit and reaction order. LV is an importantfactor for
#' conversion between molecule numbers and concentrations, the
#' appropriate exponent for LV is passed as x: c = f * LV^x * k
#'
#' @param f the value of the scaling factor (a numeric vector), derived from the unit of the knietic parameter and stoichiometry.
#' @param x the reaction order based exponent of LV
#' @param name the name of the parameter
#' @return string with scaling instructions (C)
scaleParameter <- function(f,x,name,arg="par"){
	return(
		paste0(
		sprintf("\t%s[_%s] = ",arg,name),
		ifelse(f<1,sprintf("(%s[_%s]/%g)",arg,name,1.0/f),sprintf("(%s[_%s]*%g)",arg,name,f)),
		mapply(\(x,name)
			switch(as.character(x),
			'1'="*LV;",
			'0'=";",
			'-1'="/LV;",
			'-2'="/(LV*LV);",
			sprintf("/pow(LV,%i);",x)),x,name)
		)
	)
}

#' Generate C Code to solve a model stochastically
#'
#' This function tries to generate code for a stochastic solver, but
#' assumming that the SBtab model is written in terms of
#' concentrations and rate coefficients.
#'
#' The model is interpreted by `stochasticModel()`, parameters found
#' in reaction kinetics are converted to stochastic parameters.
#' Input parameters are not converted.
#'
#' The default system size is 1 femtolitre.
#'
#' @param sb SBtab list of data-frames, as returned by
#'     SBtabVFGEN::sbtab_from_*() functions
#' @param LV Avogadro's Constant * volume [in litres]
#' @return character vector with code
#' @export
generateGillespieCode <- function(sb,LV=6.02214076e+8){
	sm <- makeGillespieModel(sb)
	KL <- sb$Reaction[["!KineticLaw"]]
	RC <- ifelse(
		grepl("-",KL,fixed=TRUE),
		KL,
		paste0(KL," - 0.0")
	)
	RC <- trimws(unlist(strsplit(RC,"-",fixed=TRUE)))
	dim(RC) <- c(2,length(KL))
	RC <- t(RC)
	Definitions <- c(
	"\t/* constants */",
	sprintf("\t double %s = %s; /* %s */",rownames(sb$Constant),sb$Constant[["!Value"]],sb$Constant[["!Unit"]]),
	"\t/* forward parameters */",
	sprintf("\tdouble %s = c[_%s]; /* per second */",rownames(sm$scf),rownames(sm$scf)),
	"\t/* backward parameters */",
	sprintf("\tdouble %s = c[_%s]; /* per second */",rownames(sm$scb),rownames(sm$scb)),
	"\t/* parameters that do not appear in kinetric laws */",
	sprintf("\tdouble %s = c[_%s]; /* per second */",rownames(sm$scv),rownames(sm$scv)),
	"\t/*state variables */"
	)
	Counts <- c(
		sprintf("\tdouble %s = x[_%s];",rownames(sb$Compound),rownames(sb$Compound)),
		sprintf("\tdouble %s = %s;",rownames(sb$Expression),sb$Expression[["!Formula"]])
	)
	FuncMolarities <- c(
	sprintf("\tdouble %s = ((double) x[_%s])/(LV * %g); /* %s */",names(sm$initialCount),names(sm$initialCount),ccc(sb),sm$compoundUnit),
	sprintf("\tdouble %s = %s;",rownames(sb$Expression),sb$Expression[["!Formula"]])
	)
	C <- c(
	"#include <stdlib.h>",
	"#include <math.h>","",
	"/* System's volume at creation time: */",
	sprintf("static const double sys_volume = %16g;",LV/6.02214076E23),
	"/* Product of Avogadro's number L (6.02214076E23) and volume V: */",
	sprintf("static const double LV = %16g;",LV),
	"/* This model was created from an SBtab file,         */",
	"/* which always describes a systems biology model of  */",
	"/* concentrations and kinetic laws. This code         */",
	"/* re-interprets it as particle counts and stochastic */",
	"/* reaction propensities. This interpretation may be  */",
	"/*       >>INCORRECT<<       ... as it is automatic.  */",
	"/* The plan is to use the kinetic law literally,      */",
	"/* as specified in the SBtab file, but to convert the */",
	"/* rate coefficients.                                 */",
	"/* We also convert initial values to particle counts. */","",
	"/* these enums make it possible to address vector elements by name, and automatically creates lengths for these vectors*/",
	paste0("enum state {",paste0("_",rownames(sb$Compound),collapse=", "),", numStateVariables};"),
	paste0("enum parameter {",paste0("_",names(sm$par),collapse=", "),", numParameters};"),
	paste0("enum reaction {",paste0("_",rownames(sb$Reaction),"_fwd",collapse=", "),", ",paste0("_",rownames(sb$Reaction),"_bwd",collapse=", "),", numReactions};"),
	paste0("enum outputFunctions {",paste0("_",rownames(sb$Output),collapse=", "),", numFunctions};"),
	"",
	"int model_effects(double t, int *x, int j){",
	sprintf("\tif (!x) return numStateVariables;"),
	sprintf("\tswitch(j){"),
	reactionEffect(sm),
	"\t}",
	"\treturn 0;",
	"}","",
	"/* Here, we want to use the kinetic law literally, with no variable substitutions.*/",
	"/* So, 'k*A*B' will be used as written, but the stochastic propensity is c*A*B    */",
	"/* where c is different from k in value and unit. We pre-calculate c/k.           */",
	"/* We convert k to have the right value (and unit, implicitly)                    */",
	"/* and use the converted value under the name k directly so as not to change the  */",
	"/* rate. This may look misleading to the reader compared to a stochastic model    */",
	"/* made from scratch.                                                             */",
	"int model_propensities(double t, int *x, double *c, double *a){",
	"\tif (!x || !a) return numReactions;",
	Definitions,
	Counts,
	unlist(mapply(\(x,n,rf) {sprintf("\ta[_%s_fwd] = %s; /* %s */",n,x,rf)},RC[,1],rownames(sb$Reaction),sub("<=>","->",sb$Reaction[["!ReactionFormula"]],fixed=TRUE))),
	unlist(mapply(\(x,n,rf) {sprintf("\ta[_%s_bwd] = %s; /* %s */",n,x,rf)},RC[,2],rownames(sb$Reaction),sub("<=>","<-",sb$Reaction[["!ReactionFormula"]],fixed=TRUE))),
	"\treturn 0;",
	"}","",
	"/* usually, these are stochastic propensity coefficients,                                      */",
	"/* but we derive our model from a concentration based kinetic form.                            */",
	"/* So, these are not directly usable, but we'll write the propensities with conversion factors */",
	"int model_reaction_coefficients(double *c){",
	"\tif (!c) return numParameters;",
	sprintf("\tc[_%s] = %s; %*s /* %s */",names(sm$par),sm$par,40-unlist(lapply(names(sm$par),length)),"",sm$parameterUnit),
	"\treturn 0;",
	"}","",
	"int model_initial_counts(int *x){",
	"\tif(!x) return numStateVariables;",
	sprintf("\tx[_%s] = lround(%g * LV); /* %g %s */",names(sm$initialCount),sm$initialCount,sm$initialCount,sm$compoundUnit),
	"\treturn 0;",
	"}","",
	"/* This function will try to convert the model's state to concentrations */",
	"/* We assume that the model was originally phrased as concentrations and kinetic laws. */",
	"int model_func(double t, int *x, double *c, double *func){",
	"\tif(!func) return numFunctions;",
	Definitions,
	FuncMolarities,
	sprintf("\tfunc[_%s] = %s;",rownames(sb$Output),sb$Output[["!Formula"]]),
	"\treturn 0;",
	"}","",
	"/* given kinetic parameters, calculate the stochastic parameters needed for propensities */",
	"/* The changes are made in place */",
	"int model_stochastic_parameters(double t, double *par){",
	"\tif (!par) return numParameters;",
	scaleParameter(sm$scf$factor,sm$scf$lvpower,rownames(sm$scf)),
	scaleParameter(sm$scb$factor,sm$scb$lvpower,rownames(sm$scb)),
	scaleParameter(sm$scv$factor,sm$scv$lvpower,rownames(sm$scv)),
	"\treturn 0;",
	"}","",
	"/* The list of experiments we'll receive from R will cointain */",
	"/* initial concentrations rather than particle counts.        */",
	"/* We need to convert them.                                   */",
	"int model_particle_count(double t, double *molarity, int *x){",
	"\tif(!molarity || !x) return numStateVariables;",
	sprintf("\tx[_%s] = lround(%g * molarity[_%s] * LV);",rownames(sb$Compound),ccc(sb),rownames(sb$Compound)),
	"\treturn 0;",
	"}"
	)
	return(C)
}


#' Simulate stochastic model
#'
#' Simulate a stochastic model generated with
#' `uqsa::generateGillespieModel()`, using the solver in this package.
#'
#' This will simulate all experimental conditions included in the list of experiments, including applying the inputs:
#' `u <- experiments[[i]]$input` - the input will be copied to the end of the model's internal parameter vector.
#'
#' Like for deterministic models, we assume that there is a vector of
#' unknown parameter (a Markov chain variable, a vector of
#' optimization variables) and also known parameters (aka the input
#' parameters). The model itself does not distinguish between the two,
#' but one is the same between the experiments and one is different
#' between different experiments: `modelParam <- c(mcmcParam, inputParam)`
#'
#' @param experiments list of expperiments, same as for the deterministic solvers.
#' @param model.so compiled C code for the model, this has to be a path with at least one slash in it, e.g.: ./model.so
#' @param parameters a numeric vector of appropriate size
#' @export
#' @return simulation result list
#' @useDynLib uqsa gillespie
simstoch <- function(experiments, model.so, parMap=identity){
	if (!file.exists(model.so)) {
		warning(sprintf("model.so «%s» not found.",model.so))
		return(NULL)
	}
	for (i in seq_along(experiments)){
		if ("input" %in% names(experiments)){
			cat(sprintf("Experiment %i, input has length: %i\n",i,length(experiments[[i]]$input)))
		} else {
			warning("experiments have no input parameters, defaults to numeric(0).")
			experiments[[i]]$input <- numeric(0)
		}
	}
	return(function(parMCMC){
		p <- as.matrix(parMap(parMCMC))
		y <- .Call(gillespie, model.so, experiments, p)
		names(y) <- names(experiments)
		for (i in seq_along(y)){
			rownames(y[[i]]$state) <- names(experiments[[i]]$initialState)
			rownames(y[[i]]$func) <- names(experiments[[i]]$outputValues)
		}
		return(y)
	})
}
