#' Represents the kinetic laws of the model as a two column matrix
#'
#' This function performs a very simplified split of a kinetic law
#' into a forward part and a backward part, if it isn't pre-split in
#' the file.
#'
#' If the data.frame contains separate forward and backward rates these will
#' be returned instead. Instead of using this function,
#' `m$Reaction[,c("fwd","bwd")]` would accomplish a very similar
#' thing.
#'
#' @param r the reaction table (data.frame)
#' @export
#' @return a character matrix with a forward and backward column
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAP79")) # not pre-split
#' print(colnames(m$Reaction))
#' k <- kinetic_law_matrix(m$Reaction)
#' print(k)
kinetic_law_matrix <- function(r){
	if ("kinetic.law" %in% names(r)){
		F <- ifelse(
			grepl("-",r$kinetic.law,fixed=TRUE),
			r$kinetic.law,
			paste0(r$kinetic.law," - 0")
		)
		F <- matrix(
			trimws(
				unlist(
					strsplit(F,"-",fixed=TRUE)
				)
			),
			ncol=2,
			byrow=TRUE,
			dimnames=list(rownames(r),c("fwd","bwd"))
		)
	} else {
		direction <- list(c("fwd","forward"),c("bwd","backward","rev","reverse"))
		F <- matrix("",NROW(r),2,dimnames=list(rownames(r),c("fwd","bwd")))
		for (d in direction){
			j <- na.omit(pmatch(d,colnames(r)))
			if (length(j)>0){
				F[,head(d,1)] <- r[[j[1]]]
			} else {
				warning(sprintf("no %s rates found in reaction table.",tail(d,1)))
			}
		}

	}
	return(F)
}


#' Returns a list of reaction coefficients
#'
#' This function maps c("A","2 B") to c(1,2)
#'
#' This is a C function because it is much easier to write in C.  C
#' has the strtod() function which expects a leading number and stops
#' when the numbers end. as.character() returns NA if the input
#' contains any dirt.
#'
#' The reaction formula is as string like this: "A + 2 B <=> C", when
#' split at `<=>` and then later at `+`, we get the strings that must
#' be parsed: "A" and "2 B" for the left side and "C" for the right
#' side. The numbers are the stoichiometric constants, or coefficients.
#' @keywords internal
#' @param formulaList a list of character vectors, derived from the
#'     left or right side of a reaction formula:
#' @return a list of numeric coefficient vectors
#' @export
#' @examples
#' print(onlyCoefficients("12 A"))
onlyCoefficients <- function(formulaList){
	C <- lapply(formulaList,\(v) .Call("lstrtod",as.character(v)))
	return(ifelse(nchar(formulaList)>0,C,0))
}

#' Returns only the names in a reaction formula
#'
#' This is the companion function to onlyCoefficients. It returns the
#' names of reactants, without the stoichiometry.
#' @keywords internal
#' @export
#' @param formulaList a list of strings like: "2 B" or "45 X"
#' @return a list of name vectors
#' @examples
#' print(onlyNames("12 A"))
onlyNames <- function(formulaList){
	return(lapply(formulaList,\(x) sub("^[^[:alpha:]]*","",x)))
}

## the best thing is to keep empty lists of reactants or products just
## empty (type nothing). But, this catches many ways of writing this.
is_empty <- function(b){
	if (nzchar(b)){
		return(trimws(b) %in% c("NULL","Nil","[]","{}","_","\u2205","\u2300","\u00D8","\u00F8"," "))
	} else {
		return(TRUE)
	}
}

#' Finds the stoichiometry for a given parameter name
#'
#' Given a list representation of stoichiometry (list of named integer
#' vectors), and a character vector of kinetic laws, this function
#' returns the correct stoichiometry for the given parameter.
#'
#' @noRd
#' @param p parameter name
#' @param reactants stoichiometry of the reaction's left side
#' @param products stoichiometry of the reaction's right side
#' @param kinetic.law a character matrix of fluxes, column 1 for
#'     forward reaction, column2 for backward reactions
#' @return the stoichiometry entry that belongs to the given parameter
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAP79"))
#' r <- stoichiometry(strsplit(m$Reaction$reactants,"+",fixed=TRUE))
#' p <- stoichiometry(strsplit(m$Reaction$products ,"+",fixed=TRUE))
#' k <- kinetic_law_matrix(m$Reaction)
#' j <- parameter_stoichiometry("k5_1",r,p,k)
#' print(j)
parameter_stoichiometry <- function(p,reactants,products,kinetic.law){
	i <- grepl(paste0("\\b",p,"\\b"),kinetic.law[,1])
	j <- grepl(paste0("\\b",p,"\\b"),kinetic.law[,2])
	if (any(i)){
		return(reactants[[which(i)[1]]])
	}
	if (any(j)){
		return(products[[which(j)[1]]])
	} else {
		return(NULL)
	}
}

#' Returns a list of named stoichiometric vectors
#'
#' Given an already split list of entries such as c("3 A","B"), this
#' function returns a numeric vector c(3,1) with names c("A","B").
#' @export
#' @param formulaList reaction formulae, either as a pre-split list or character vector
#' @return named numeric vector of stoichiometric coefficients
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAP79"))
#' r <- stoichiometry(m$Reaction$reactants)
#' print(head(r))
#' print(tail(r))
stoichiometry <- function(formulaList){
	if (is.character(formulaList)){
		nm <- names(formulaList)
		formulaList <- lapply(strsplit(formulaList,"+",fixed=TRUE),trimws)
		names(formulaList) <- nm
	}
	sc <- onlyCoefficients(formulaList)
	nm <- lapply(onlyNames(formulaList),trimws)
	sc <- mapply(
		function(COUNTs,NAMEs) {
			if (length(COUNTs)==0) {
				return(NULL)
			} else if (length(COUNTs) == 1 && is_empty(NAMEs)){
				return(NULL)
			} else {
				return(tapply(COUNTs,NAMEs,sum))
			}
		},
		sc,
		nm,
		SIMPLIFY=FALSE) # make all the vectors in the list named
	names(sc) <- names(formulaList)
	return(sc)
}


#' Initial count of reacting Compounds
#'
#' This function returns the molecule count apart from LV.  This
#' returned number must be multiplied by Avogadro's constant and
#' volume.
#'
#' @noRd
#' @param m the model data.frames obtained from [model_from_tsv]
#' @param LV product of Avogadro's number and volume
#' @return numeric vector or character vector, depending on how the
#'     initial concentration was provided
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAP79"))
#' print(initialCount(m)
initialCount <- function(m,LV){
	v <- values(m$Compound)
	unit <- units_from_table(m$Compound)
	u <- lapply(unit,unit.from.string)
	## f is the conversion-factor between the given substance units and mol/liter
	f <- unlist(lapply(u,\(x) 10^sum(x$scale*x$exponent)))
	if (is.numeric(v)){
		IC <- round(v*f*LV); # LV is in: particles/mol * liter
	} else {
		## the initial value is a formula of some sort,
		## an expression of parameters, constants, and other expressions.
		IC <- sprintf("(%s) * %g * LV",v,f)
		names(IC) <- names(v)
	}
	attr(IC,"unit") <- rep("1",length(IC))
	return(IC)
}

#' Determines a printable reaction Formula from the Stoichiometry
#'
#' Given a list of reactants and products, this function will print a
#' reaction fomrula for each list item.
#'
#' The stoichiometry (sparse representation), consists of two lists,
#' giving the names and quantities for each reaction.
#'
#' @param reactants a list of named numeric vector
#' @param products a list of named numeric vectors
#' @return a character vector with reaction formulas
#' @noRd
#' @examples
#' rf <- reaction_formula_from_stoichiometry(
#'   list(c(A=1,B=2),c(C=3)),
#'   list(c(C=3),c(A=1,B=2))
#' )
#' cat(rf,sep="\n")
reaction_formula_from_stoichiometry <- function(reactants,products){
	stopifnot(length(reactants)==length(products))
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
#' @return a matrix with two columns, forward fluxes and backward
#'     fluxes, and as many rows as there are reactions.
#' @noRd
flux_matrix <- function(Reaction){
	kl <- ifelse(
		grepl("-",Reaction$kinetic.law,fixed=TRUE),
		Reaction$kinetic.law,
		paste0(Reaction$kinetic.law," - 0.0")
	)
	names(kl) <- rownames(Reaction)
	return(
		t(
			matrix(
				trimws(
					unlist(strsplit(kl,"-",fixed=TRUE))
				),
				nrow=2,
				dimnames=list(c('fwd','bwd'),rownames(Reaction))
			)
		)
	)
}

#' Interprets the provided biological reaction model as stochastic
#'
#' The chemical master equation can be simulated as a Markov jump
#' process (or continuous time Markov chain). One of the stochastic
#' solver algorithms is the Gillespie algorithm. This function return
#' sa data structure that can be used to generate code for the
#' Gillespie solver in this package.
#'
#' This function interprets the continuous model `m` as a discrete
#' state model with molecule counts and propensities. For this reason,
#' we need to specify a volume for the simulations to take place in.
#'
#' The model `m` is assumed to describe a reaction network, as a list
#' of data.frames (as retuned by [model_from_tsv]).  The systems
#' biology information in the file is assumed to be concentrations and
#' rate coefficients, regardless of the interpretation this function
#' will derive from it. This is to make the model format of the TSV
#' file fairly uniform and independent of how we want to solve the
#' derived equations, be it ODE or CME.
#'
#' Like the ode object, the returned object can also store the paths
#' of files we create for this model, with: [c_path<-], and [so_path<-]
#'
#' With the information provided with the rate coefficient units and a
#' volume, this function tries to convert everything to Gillespie rate
#' constants.
#' @param m list of data.frames, obtained via `model_from_tsv()`
#' @param LV product of Avogadro's number (Josef Loschmidt) and system volume
#' @return a list containing the interpreted model.
#' @export
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' stochasticModel <- as_cme(m)
#' print(stochasticModel)
as_cme <- function(m,LV=6.02214076e+8){
	rev <- as.logical(m$Reaction[["is.reversible"]])
	F <- matrix(
		c(
			sub("*"," ",trimws(m$Reaction$reactants),fixed=TRUE),
			sub("*"," ",trimws(m$Reaction$products),fixed=TRUE)
		),
		ncol=2
	)
	cl <- stoichiometry(strsplit(F[,1],"+",fixed=TRUE))
	names(cl) <- paste0(rownames(m$Reaction),"_fwd")
	cr <- stoichiometry(strsplit(F[,2],"+",fixed=TRUE))
	names(cr) <- paste0(rownames(m$Reaction),"_bwd")
	k <- c(
		linear_scale(m$Parameter),
		linear_scale(m$Input)
	)
	u <- c(
		units_from_table(m$Parameter),
		units_from_table(m$Input)
	)
	attr(k,"unit") <- u
	kl <- flux_matrix(m$Reaction)
	ic <- linear_scale(m$Compound) # will have units, as linear_scale returns them
	mcc <- moleCountConversion(attr(ic,"unit"))
	pcc <- moleCountConversion(attr(k,"unit"))
	stochastic_par <- convert_units(k,LV)
	constants <- linear_scale(m$Constant) # preserves unit attribute
	CME <- list(
		left=cl,                              # reactants
		right=cr,                             # products
		k=k,                                  # original (kinetic) parameters
		initialConcentration=ic,
		mcc=mcc,                              # mole to count conversion structure
		initialCount=initialCount(m,LV),      # as number of particles
		par=stochastic_par,                   # Gillespie parameters
		pcc=pcc,                              # parameter conversion to counts and seconds
		expression=formulae(m$Expression),
		reactionMultiplicityFWD=sapply(cl,\(x) 1+sum(x-1)),
		reactionMultiplicityBWD=sapply(cr,\(x) 1+sum(x-1)),
		const=constants,
		kinetic.law=kl,
		reaction.formula=cbind(
			fwd=paste(m$Reaction$reactants,m$Reaction$products,sep=" -> "),
			bwd=paste(m$Reaction$reactants,m$Reaction$products,sep=" <- ")
		),
		output=formulae(m$Output),
		c_path=NULL,
		c.date=NULL,
		so_path=NULL,
		so.date=NULL,
		LV=LV,
		name=comment(m)
	)
	comment(CME) <- comment(m)
	class(CME) <- "cme"
	return(CME)
}

reactionEffect <- function(sm){
	return(
		unlist(
			mapply(
				\(x,y,n,rf) {
					c(
						sprintf("\tcase _%s: /* %s */",n,rf),
						sprintf("\t\tx[_%s] -= %i;",names(x),x),
						sprintf("\t\tx[_%s] += %i;",names(y),y),
						sprintf("\t\tbreak;")
					)
				},
				c(sm$left,sm$right),
				c(sm$right,sm$left),
				c(names(sm$left),names(sm$right)),
				sm$reaction.formula,
				SIMPLIFY=FALSE
			)
		)
	)
}

lengths <- function(v){
	if (is.null(v) || length(v)==0) return(0)
	else return(nchar(as.character(v)))
}

padding <- function(...,upper.bound=40){
	return(pmax(0,upper.bound - lengths(paste(...))))
}

#' Returns a string that contains code to scale the parameter
#'
#' The goal is to convert the kinetic parameter k to the stochastic
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
#' @noRd
scaleParameter <- function(f,x,name,arg="par"){
	return(
		paste0(
		sprintf("\t%s[_%s] = ",arg,name),
		ifelse(f<1,sprintf("(%s[_%s]/%g)",arg,name,1/f),sprintf("(%s[_%s]*%g)",arg,name,f)),
		mapply(\(x,name)
			switch(as.character(x),
			'1'="*LV;",
			'0'=";",
			'-1'="/LV;",
			'-2'="/(LV*LV);",
			sprintf("*pow(LV,%i);",x)),x,name)
		)
	)
}

spacing <- function(v,max.width=40){
	if (is.null(v)) return(max.width)
	if (!is.null(names(v))) return(max.width - lengths(v) - lengths(names(v)))
	else return(max.width - lengths(v))
}

prod_1 <- function(str, n) {
	if (n <= 1) return(str)
	return(
		paste(
			c(str, sprintf("(%s-%i)", str, seq_len(n-1))),
			collapse="*"
		)
	)
}

rate_to_propensity <- function(reaction.rate,verbose=getOption("uqsa.verbose",TRUE)){
	if (any(grepl("[-+/(]",reaction.rate))){
		cli::cli_alert_info("reaction rate '{reaction.rate}' is not a mass action reaction.")
		return(reaction.rate)
	}
	extract_factors <- function(e) {
		if (is.call(e) && identical(e[[1]], as.name("*"))) {
			c(extract_factors(e[[2]]), extract_factors(e[[3]]))
		} else if (is.call(e) && identical(e[[1]], as.name("^")) && is.numeric(e[[3]])) {
			rep(list(e[[2]]), e[[3]])
		} else {
			list(e)
		}
	}
	## factors of mass action reaction
	f <- vapply(extract_factors(str2lang(reaction.rate)),deparse1,character(1))
	u <- unique(f)
	n <- setNames(tabulate(match(f, u)), u)
	prop <- mapply(prod_1,u,n)
	return(paste(prop,collapse="*"))
}

## # alternative implementation for rate to propensity for mass action reactions with up to two reactants:
## x <- re_substitutes(
## 	x,
## 	rex(boundary,capture(alpha,any_of(alpha,digit,"_"),name="A"),"*",capture_group("A"),boundary),
## 	sprintf("%s*(%s-1)",cg(1),cg(1))
## )
## x <- re_substitutes(
## 	x,
## 	rex(boundary,capture(alpha,any_of(alpha,digit,"_"),name="A"),"^2",boundary),
## 	sprintf("%s*(%s-1)",cg(1),cg(1))
## )
## x <- sub("\\b([[:alnum:]_]+)^2\\b","\\1*(\\1-1)",x)
## fmt <- "\ta[_%s%s] = %s; %*s /* 0.5*%i*%s (%s)*/"


#' Generates C code for Propensities
#'
#' Propensities can be slightly different depending on whether the
#' reaction is combining two different molecules or two of the same.
#'
#' `A` and `B` shall be two different compounds and the number of
#' molecules for each, while `a` and `b` are the concentrations for
#' the two substances (e.g. in nM).
#'
#' Combinatorics:
#'
#' - in `A + B -> C` the number of combinations is `A*B`
#' - in `A + A -> C` the number of uniqe combinations is `0.5*A*(A-1)`
#'
#' Propensity: consists of a parameter (factor) `c`, and the number of combinations, e.g.: `c*A*B`
#'
#' We are provided with a string from the kinetic formulation of the
#' model, e.g.: kf*a*b, using some unit of cencentration in kf and the concentrations, e.g. nM
#' (nanomoles per litre): kf is in 1/(nM*s), a & b are in nM.
#'
#' So, the two expressions are very similar.
#'
#' Given a function `f(kf,a,b) = kf*a*b`, we can plug in counts into
#' the slots for `a` and `b` and convert the unit of `k` to remain
#' compatible, we can get the propensity like this: `f(kf*nM,A,B)`
#'
#' For heterodimerization, and similar processes, it is a bit different:
#'
#' kinetic law: kf*A*A
#' propensity: 0.5*c*a*(a-1)
#'
#' So, given the user supplied function `g(kf,a)` we can plug in `g(0.5*c,A)`
#' the appropriate `c` is `2*kf*nM` (to get rid of the unit).
#'
#' The 2 and 0.5 cancel, so we can plug in directly: `g(kf*nM,A)`, we
#' can also try to do string substitution in `g` (it is given as a
#' user supplied string).
#'
#' @param multiplicity reaction multiplicity, 1 for heterodimerization, 2 for homodimerization
#' @param flux values of the fluxes as character vector
#' @param suffix e.g. "_fwd"
#' @param formula reaction formula (will become a comment)
#' @param do.substitution replace A*A with A*(A-1), and A^2 with A*(A-1)
#' @noRd
#' @return character vector of c code
propensity_in_c <- function(multiplicity,flux,suffix,formula,do.substitution=TRUE){
	cg <- \(n) sprintf(r"(\%i)",n) # captured groups
	propensity <- unlist(
		mapply(
			function(n,f,x,rf) {
				if (f>1 && isTRUE(do.substitution)){
					p <- rate_to_propensity(x)
					fmt <- "\ta[_%s%s] = %s; %*s /* 0.5*%i*%s (%s) made from flux %s */"
				} else {
					p <- x
					fmt <- "\ta[_%s%s] = %s; %*s /* %.0s%s (%s) %.0s*/"
				}
				return(
					sprintf(
						fmt,
						n,suffix,
						p,
						padding(n,suffix,x)," ",
						f,p,
						rf,x
					)
				)
			},
			names(flux),
			multiplicity,
			flux,
			formula
		)
	)
	return(propensity)
}

#' Generates C Code to solve a model stochastically
#'
#' This function tries to generate code for a stochastic solver, but
#' assumming that the systems biology model is written in terms of
#' concentrations and rate coefficients. The stochastic model `sm` is
#' of the type returned by `as_cme`.
#'
#' The arameters found in reaction kinetics are converted to
#' stochastic parameters.  Input parameters are not converted.
#'
#' The default system volume is 1 femtolitre.  As volumes are quite
#' small in cellular systems, we specify L * V, where L is Avogadro's
#' constant, sometimes written as $N_\text{A}$. The volume V needs to
#' be given in a unit that is compatible with the concentration unit
#' (disregarding SI-prefixes). If the concentrations are given in nM
#' (nanomole/L), then the volume needs to be in litres. Whenever the
#' volume is a very small number (and Avogadro's constant just is a
#' very big number), as a product, they often become a reasonable in
#' scale.
#'
#' @param sm a stochatsic Gillespie model obtained via `makeGillespieModel`
#' @return character vector with code
#' @noRd
generateGillespieCode <- function(sm){
	LV <- sm$LV
	C_Constants <- sprintf(
			"\t double %s = %s; %*s /* %s */",
			names(sm$const),
			sm$const,
			padding(sm$const,names(sm$const))," ",
			attr(sm$const,"unit")
	)
	C_Expressions <- sprintf("\tdouble %s = %s;",names(sm$expression),replace_powers(sm$expression))
	Definitions <- c(
		"\t/* constants */",
		C_Constants,
		"\t/* parameters */",
		sprintf(
			"\tdouble %s = c[_%s]; %*s /* %s; originally: %s = %.2g %s */",
			names(sm$par),
			names(sm$par),
			padding(names(sm$par),names(sm$par))," ",
			attr(sm$par,"unit"),
			names(sm$k),
			sm$k,
			attr(sm$k,"unit")
		),
		"\t/*state variables */"
	)
	C_Counts <- c(
		sprintf("\tdouble %s = x[_%s];",names(sm$initialCount),names(sm$initialCount)),
		C_Expressions
	)
	FuncMolarities <- c(
		sprintf(
			"\tdouble %s = c[_%s]/(%g * pow(LV,%i)); %*s /* %s */",
			names(sm$k),
			names(sm$par),
			sm$pcc$factor,
			sm$pcc$lvpower,
			padding(names(sm$k),names(sm$par),sm$pcc$factor,sm$pcc$lvpower)," ",
			attr(sm$k,"unit")
		),
		sprintf(
			"\tdouble %s = ((double) x[_%s])/(LV * %g); %*s /* %s */",
			names(sm$initialCount),
			names(sm$initialCount),
			sm$mcc$factor,
			padding(names(sm$initialCount),names(sm$initialCount),sm$mcc$factor),"",
			attr(sm$initialConcentration,"unit") # original unit
		),
		C_Expressions
	)
	if (is.numeric(sm$initialCount)){
		C_initialCount <- sprintf(
			"\tx[_%s] = %i; %*s /* originally: %.2g %s */",
			names(sm$initialCount),
			sm$initialCount,
			padding(sm$initialCount,names(sm$initialCount))," ",
			sm$initialConcentration,
			attr(sm$initialConcentration,"unit")
		)
	} else if (is.character(sm$initialCount)){
		C_initialCount <- c(
			Definitions,
			sprintf(
				"\tx[_%s] = lround(%s); %*s /* originally: %s %s */",
				names(sm$initialCount),
				sm$initialCount,
				padding(sm$initialCount,names(sm$initialCount))," ",
				sm$initialConcentration,
				attr(sm$initialConcentration,"unit")
			)
		)
	} else {
		stop("initial count is broken.")
	}
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
	paste0("enum reaction {",paste0("_",rownames(sm$kinetic.law),"_fwd",collapse=", "),", ",paste0("_",rownames(sm$kinetic.law),"_bwd",collapse=", "),", numReactions};"),
	paste0("enum outputFunctions {",paste0("_",names(sm$output),collapse=", "),", numFunctions};"),
	"",
	"int model_effects(double t, int *x, double *c, double* a, int j){",
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
	C_Counts,
	propensity_in_c(
		sm$reactionMultiplicityFWD,
		sm$kinetic.law[,'fwd'],
		"_fwd",
		reaction_formula_from_stoichiometry(sm$left,sm$right)
	),
	propensity_in_c(
		sm$reactionMultiplicityBWD,
		sm$kinetic.law[,'bwd'],
		"_bwd",
		reaction_formula_from_stoichiometry(sm$right,sm$left)
	),
	"\treturn 0;",
	"}","",
	"/* usually, these are stochastic propensity coefficients,                                      */",
	"/* but we derive our model from a concentration based kinetic form.                            */",
	"/* So, these are not directly usable, but we'll write the propensities with conversion factors */",
	"int model_reaction_coefficients(double *c){",
	"\tif (!c) return numParameters;",
	sprintf(
		"\tc[_%s] = %.6g; %*s /* %s */",
		names(sm$par),
		sm$par,
		padding(names(sm$par),sprintf("%.6g",sm$par)),"",
		attr(sm$par,"unit")
	),
	"\treturn 0;",
	"}","",
	"int model_initial_counts(int *x, double *c){",
	"\tif(!x) return numStateVariables;",
	C_initialCount,
	"\treturn 0;",
	"}","",
	"/* This function will try to convert the model's state to concentrations */",
	"/* We assume that the model was originally phrased as concentrations and kinetic laws. */",
	"int model_func(double t, int *x, double *c, double *func){",
	"\tif(!func) return numFunctions;",
	C_Constants,
	FuncMolarities,
	sprintf("\tfunc[_%s] = %s;",names(sm$output),sm$output),
	"\treturn 0;",
	"}","",
	"/* given kinetic parameters, calculate the stochastic parameters needed for propensities */",
	"/* The changes are made in place */",
	"int model_stochastic_parameters(double t, double *par){",
	"\tif (!par) return numParameters;",
	scaleParameter(
		sm$pcc$factor,
		sm$pcc$lvpower,
		rownames(sm$pcc)
	),
	"\treturn 0;",
	"}","",
	"/* The list of experiments we'll receive from R will contain  */",
	"/* initial concentrations rather than particle counts.        */",
	"/* We need to convert them.                                   */",
	"int model_particle_count(double t, double *molarity, int *x){",
	"\tif(!molarity || !x) return numStateVariables;",
	sprintf(
		"\tx[_%s] = lround(%g * molarity[_%s] * LV); %*s /* from: %s */",
		names(sm$initialCount),
		sm$mcc$factor,
		names(sm$initialConcentration),
		padding(names(sm$initialCount),sm$mcc$factor,names(sm$initialConcentration))," ",
		attr(sm$initialConcentration,"unit")
	),
	"\treturn 0;",
	"}"
	)
	return(C)
}


#' Simulates a Stochastic Model
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
#' The path to the shared library, is required to contain at least one
#' slash in it, e.g.: "./model.so", "/tmp/Rsdkljhskjdhf/model.so" But,
#' not just "model.so", otherwise the shared library is interpreted as
#' a system library by `dlopen()` (it will not be found).
#'
#' The number of reactions can be limited by `nstep` (micro time-steps
#' forward). The maximum can be set by performing a good simulation
#' and reading out the number of steps taken in that reference
#' simulation: `y[[i]]$numSteps`.
#'
#' @param ex list of experiments, same as for the deterministic
#'     solvers.
#' @param cmeModel Either the cmeModel from [as_cme], with a shared
#'     library path stored inside, or the path to the so file
#' @param parMap map from MCMC variable (or ABC variable) to
#'     model-parameters.
#' @param time.out in seconds
#' @param nstep number of reactions for early exit, defaults to
#'     unlimited (0)
#' @export
#' @return a closure that simulates the model in `model.so`
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' cme <- as_cme(m)
#' C <- generate_code(cme)
#' c_path(cme) <- write_c_code(C)
#' so_path(cme) <- shlib(cme)
#' ex <- experiments(m)
#' p0 <- values(m$Parameter)
#' s <- simstoch(ex,cme)
#' res <- s(p0)
#' require(errors)
#' plot(as.errors(ex[[1]]$outputTimes),ex[[1]]$data,xlab="time",ylab="AKAR4p",main=names(ex)[1])
#' lines(ex[[1]]$outputTimes,res[[1]]$func,type="s",lwd=2,col="red3")
simstoch <- function(ex, cmeModel, parMap=identity, time.out=1, nstep=0){
	if (is.character(cmeModel)) {
		model.so <- cmeModel
	} else if (is(cmeModel,"cme")) {
		model.so <- so_path(cmeModel)
	}
	if (!grepl("/",model.so)) {
		stop(sprintf("The shared library path must have at least one slash: %s",model.so))
	}
	if (!file.exists(model.so)) {
		warning(sprintf("The file \u00ab%s\u00bb does not exist.",model.so))
		return(NULL)
	}
	for (i in seq_along(ex)){
		if (!("input" %in% names(ex[[i]]))){
			ex[[i]]$input <- numeric(0)
		}
	}
	return(
		function(parMCMC){
			p <- as.matrix(parMap(parMCMC))
			y <- .Call("gillespie", model.so, ex, p, as.double(time.out), nstep)
			names(y) <- names(ex)
			stopifnot(length(y) == length(ex))
			for (i in seq_along(y)){
				rownames(y[[i]]$state) <- names(ex[[i]]$initialState)
				rownames(y[[i]]$func) <- rownames(ex[[i]]$data)
			}
			class(y) <- "simulation"
			return(y)
		}
	)
}
