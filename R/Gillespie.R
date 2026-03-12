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
	return(ifelse(nchar(formulaList)>0,C,0))
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

is_empty <- function(b){
	return(b == "NULL" | b == "_" | b == "Ø" | b == "\U00002205" | b=="")
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
	sc <- mapply(
		function(a,b) {
			if (length(a)==0) {
				return(NULL)
			} else if (length(a) == 1 && is_empty(b)){
				return(NULL)
			} else {
				names(a) <- trimws(b)
				return(a)
			}
		},
		sc,
		nm,
		SIMPLIFY=FALSE) # make all the vectors in the list named
	names(sc) <- names(formulaList)
	return(sc)
}

#' Find the Reaction a parameter appears in
#'
#' Usually Parameters are a list of names and values, Reactions have
#' reaction kinetics, which can contain any number of parameters. To
#' make matters worse, a parameter may influence two or more
#' reactions. The goal is to decide on the parameter's conversion
#' factor. We shall assume that if the parameter appears in more than
#' one reaction, it is the same kind of context, and thus it must be
#' converted exactly the sam eway.  This function finds the first
#' reaction a given parameter appears in.
#'
#' Note that the entire conversion may not work at all for complicated
#' kinetics, we use very simple rules here.
#'
#' @param reactionKinetic a character vector with all of the kinetic
#'     law expressions
#' @param parNames parameter names
#' @return list with two named components: fwd, bwd; each is a named
#'     vector of integers, th enames are taken from parNames, the
#'     integer indicates a reaction.
parameterReaction <- function(reactionKinetic,parNames){
	N <- length(reactionKinetic)
	RC <- ifelse(
		grepl("-",reactionKinetic,fixed=TRUE),
		reactionKinetic,
		paste0(reactionKinetic," - 0")
	)
	RC <- trimws(unlist(strsplit(RC,"-",fixed=TRUE)))
	dim(RC) <- c(2,N)
	RC <- t(RC)
	# forward and backward reaction kinetics
	k_fwd <- lapply(
		parNames,
		function(x){
			return(
				head(
					grep(
						paste0("\\b",x,"\\b"),RC[,1]
					),
					1
				)
			)
		}
	)
	names(k_fwd) <- parNames
	k_bwd <- lapply(
		parNames,
		function(x) {
			return(
				head(
					grep(
						paste0("\\b",x,"\\b"),RC[,2]
					),
					1
				)
			)
		}
	)
	names(k_bwd) <- parNames
	return(list(fwd=unlist(k_fwd),bwd=unlist(k_bwd)))
}

initialCount <- function(m){
	v <- values(m$Compound)
	unit <- m$Compound$unit
	u <- lapply(unit,unit.from.string)
	f <- unlist(lapply(u,\(x) 10^sum(x$scale*x$exponent)))
	if (is.numeric(v)){
		IC <- v*f
	} else {
		IC <- sprintf("%s * %s",v,f)
		names(IC) <- names(v)
	}
	attr(IC,"unit") <- unit
	return(IC)
}

reaction_formula_from_stoichiometry <- function(reactants,products){
	return(
		paste(
			sapply(reactants,FUN=\(r) paste(sprintf("%i %s",r,names(r)),collapse=" + ")),
			sapply(products,FUN=\(p) paste(sprintf("%i %s",p,names(p)),collapse=" + ")),
			sep=" -> "
		)
	)
}

#' Creates a Matrix with reaction kinetic entries
#'
#' Given a data.frame that describes reaction, this function extracts
#' the forward flux and backward flux from the kinetic.law column.
#'
#' @param Reaction a data.frame with a column named 'kinetic.law'
flux_matrix <- function(Reaction){
	kl <- ifelse(
		grepl("-",KL,fixed=TRUE),
		Reaction$kinetic.law,
		paste0(KL," - 0.0")
	)
	names(kl) <- rownames(Reaction)
	return(
		t(matrix(
			trimws(
				unlist(
					strsplit(kl,"-",fixed=TRUE)
				)
			),
			nrow=2,
			dimnames=list(c('fwd','bwd'),rownames(Reaction))
		))
	)
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
#' @param m list of data.frames, obtained via `model_from_tsv()`
#' @param LV Avogadro's constant L multiplied by the system's volume V.
#' @return a list containing the interpreted model.
#' @export
makeGillespieModel <- function(m){
	rev <- as.logical(m$Reaction[["is.reversible"]])
	F <- matrix(
		c(
			sub("*"," ",trimws(m$Reaction$reactants),fixed=TRUE),
			sub("*"," ",trimws(m$Reaction$reactants),fixed=TRUE)
		),
		ncol=2
	)
	cl <- stoichiometry(strsplit(F[,1],"+",fixed=TRUE))
	names(cl) <- paste0(rownames(m$Reaction),"_fwd")
	cr <- stoichiometry(strsplit(F[,2],"+",fixed=TRUE))
	names(cr) <- paste0(rownames(m$Reaction),"_bwd")
	k <- c(values(m$Parameter),values(m$Input))
	u <- c(units_from_table(m$Parameter),units_from_table(m$Input))
	attr(k,"unit") <- u
	kl <- flux_matrix(m$Reaction)
	return(
		list(
			left=cl,
			right=cr,
			parConversion=parameterConversion(u,cl,cr,kl),
			initialCount=initialCount(m),
			par=k,
			expression=formulae(m$Expression),
			reactionMultiplicityFWD=sapply(cl,\(x) 1+sum(x-1)),
			reactionMultiplicityBWD=sapply(cr,\(x) 1+sum(x-1)),
			const=values(m$Constant),
			kinetic.law=kl,
			output=formulae(m$Output)
		)
	)
}

reactionEffect <- function(sm){
	return(
		unlist(
			c(
				mapply(
				\(x,y,n) c(sprintf("\tcase _%s:",n),sprintf("\t\tx[_%s] -= %i;",names(x),x),sprintf("\t\tx[_%s] += %i;",names(y),y),sprintf("\t\tbreak;")),
				sm$left,
				sm$right,
				names(sm$left),
				SIMPLIFY=FALSE
				),
				mapply(
				\(x,y,n) c(sprintf("\tcase _%s:",n),sprintf("\t\tx[_%s] -= %i;",names(x),x),sprintf("\t\tx[_%s] += %i;",names(y),y),sprintf("\t\tbreak;")),
				sm$right,
				sm$left,
				names(sm$right),
				SIMPLIFY=FALSE
				)
			)
		)
	)
}

lengths <- function(v){
	if (is.null(v) || length(v)==0) return(0)
	else return(sapply(as.character(v),nchar))
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
#' @param f the value of the scaling factor (a numeric vector),
#'     derived from the unit of the kinetic parameter and
#'     stoichiometry.
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

spacing <- function(v,max.width=40){
	if (is.null(v)) return(max.width)
	if (!is.null(names(v))) return(max.width - lengths(v) - lengths(names(v)))
	else return(max.width - lengths(v))
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
#' @param sb a stochatsic Gillespie model obtained via `makeGillespieModel`
#' @param LV Avogadro's Constant * volume (in litres)
#' @return character vector with code
#' @export
generateGillespieCode <- function(sm,LV=6.02214076e+8){
	RC <- trimws(unlist(strsplit(sm$kinetic.law,"-",fixed=TRUE)))
	dim(RC) <- c(2,length(KL))
	RC <- t(RC)
	Definitions <- c(
		"\t/* constants */",
		sprintf(
			"\t double %s = %s; %*s /* %s */",
			names(sm$const),
			sm$const,
			43 - lengths(sm$const) - lengths(names(sm$const))," ",
			sm$const %@% "unit"
		),
		"\t/* parameters */",
		sprintf(
			"\tdouble %s = c[_%s]; %*s /* %s */",
			names(sm$par),
			names(sm$par),
			40 - 2*lengths(names(sm$par))," ",
			sm$parConversion$effectively
		),
		"\t/*state variables */"
	)
	Counts <- c(
		sprintf("\tdouble %s = x[_%s];",names(sm$initialCount),names(sm$initialCount)),
		sprintf("\tdouble %s = %s;",names(sm$expression),replace_powers(sm$expression))
	)
	conc_to_count <- moleCountConversion(sm$initialCount %@% "unit")
	FuncMolarities <- c(
		sprintf(
			"\tdouble %s = ((double) x[_%s])/(LV * %g); %*s /* %s */",
			names(sm$initialCount),
			names(sm$initialCount),
			conc_to_count$factor,
			40 - 2*lengths(names(sm$initialCount))-nchar(conc_to_count$factor),"",
			sm$initialCount %@% "unit" # original unit
		),
		sprintf(
			"\tdouble %s = %s;",
			names(sm$expression),
			replace_powers(formulae(sm$expression))
		)
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
	paste0("enum state {",paste0("_",names(sm$initialCount),collapse=", "),", numStateVariables};"),
	paste0("enum parameter {",paste0("_",names(sm$par),collapse=", "),", numParameters};"),
	paste0("enum reaction {",paste0("_",names(sm$kinetic.law),"_fwd",collapse=", "),", ",paste0("_",rownames(sb$Reaction),"_bwd",collapse=", "),", numReactions};"),
	paste0("enum outputFunctions {",paste0("_",names(sm$output),collapse=", "),", numFunctions};"),
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
	unlist(
		mapply(
			\(x,n,f,rf) {
				sprintf(
					"\ta[_%s_fwd] = %i*%s; %*s /* %s */",
					n,f,x,
					40-lengths(x)-lengths(n)," ",
					rf
				)
			},
			RC[,1],
			rownames(sm$kinetic.law),
			sm$reactionMultiplicityFWD,
			reaction_formula_from_stoichiometry(sm$left,sm$right)
		)
	),
	unlist(
		mapply(
			\(x,n,f,rf) {
				sprintf(
					"\ta[_%s_bwd] = %i*%s; %*s /* %s */",
					n,f,x,
					40-lengths(x)-lengths(n)," ",
					rf
				)
			},
			RC[,2],
			rownames(sb$Reaction),
			sm$reactionMultiplicityBWD,
			reaction_formula_from_stoichiometry(sm$right,sm$left)
		)
	),
	"\treturn 0;",
	"}","",
	"/* usually, these are stochastic propensity coefficients,                                      */",
	"/* but we derive our model from a concentration based kinetic form.                            */",
	"/* So, these are not directly usable, but we'll write the propensities with conversion factors */",
	"int model_reaction_coefficients(double *c){",
	"\tif (!c) return numParameters;",
	sprintf("\tc[_%s] = %s; %*s /* %s */",names(sm$par),sm$par,40-lengths(names(sm$par))-lengths(sm$par),"",sm$par %@% "unit"),
	"\treturn 0;",
	"}","",
	"int model_initial_counts(int *x, double *c){",
	"\tif(!x) return numStateVariables;",
	Definitions,
	sprintf(
		"\tx[_%s] = lround(%s * LV); %*s /* %s %s */",
		names(sm$initialCount),
		as.character(sm$initialCount),
		40-nchar(as.character(sm$initialCount))-nchar(names(sm$initialCount))," ",
		as.character(sm$initialCount),
		sm$initialCount %@% "unit"
	),
	"\treturn 0;",
	"}","",
	"/* This function will try to convert the model's state to concentrations */",
	"/* We assume that the model was originally phrased as concentrations and kinetic laws. */",
	"int model_func(double t, int *x, double *c, double *func){",
	"\tif(!func) return numFunctions;",
	Definitions,
	FuncMolarities,
	sprintf("\tfunc[_%s] = %s;",names(sm$output),sm$output),
	"\treturn 0;",
	"}","",
	"/* given kinetic parameters, calculate the stochastic parameters needed for propensities */",
	"/* The changes are made in place */",
	"int model_stochastic_parameters(double t, double *par){",
	"\tif (!par) return numParameters;",
	scaleParameter(
		sm$parConversion$factor,
		sm$parConversion$lvpower,
		rownames(sm$parConversion)
	),
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
		if ("input" %in% names(experiments[[i]])){
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
