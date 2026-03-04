## Here we define functions that are meant to replace SBtab, by making
## it as simple as possible, with no need to remember anything.

#' model_from_tsv loads the content from a series of tsv files
#'
#' The argument can be either a series of tsv file-names, or a
#' directory with tsv files. If it is a directory, all tsv files
#' therein wil be used.
#'
#' @export
#' @return a list of data.frames, one per file, named like the files.
#' @param src either a vector of files, or a directory with tsv files
model_from_tsv <- function(src="."){
	stopifnot(is.character(src))
	if (length(src)==1 && dir.exists(src)){
		modelName <- basename(normalizePath(src))
		src <- dir(src,pattern="tsv$")
	} else {
		modelName <- basename(dirname(normalizePath(src[1])))
	}
	stopifnot(all(file.exists(src)))
	if (is.null(names(src))){
		names(src) <- sub("[.]tsv$","",basename(src))
	}
	m <- lapply(src,read.delim,row.names=1)
	comment(m) <- modelName
	names(m) <- names(src)
	return(m)
}


#' Find values in a data.frame that is derived from a tsv file or similar
#'
#' given a data.frame, this function will look for a column that
#' contains some kind of value and retrieve it. The returned numeric
#' vector will be named.
#'
#' If the values contain a standard error, the returned value is of
#' class _errors_.
#'
#' @param df a data frame with a "value" column
values <- function(df){
	if ("value" %in% colnames(df)){
		v <- df$value
	} else {
		v <- df[[grep("[Vv]alue",colnames(df))[1]]]
	}
	if (is.null(v)) stop("data.frame contains no values column.")
	if (is.character(v)){
		v <- parse_concise(v)
	}
	names(v) <- rownames(df)
	return(v)
}

#' Find a column that contains some kind of mathematic expression in a data.frame
#'
#' Given a data.frame that should contain a column that assigns a math
#' expression to a name (in row names), this function returns a named
#' character vector with the expressions. The formula column shoul dbe
#' named "formula" (if it exists, only this column will be used). But
#' some other spellings will also work as fallback. As a fallback
#' "value" is acceptable as well, because it makes sense to say
#' "the value of x is 'y/2+1'", even though it is not an atomic value
#' (but an expression).
#'
#' @param df a data.frame with a "formula" column
#' @return character vector with names taken from the row names of df
formulae <- function(df){
	if ("formula" %in% colnames(df)){
		f <- df$formula
	} else {
		f <- df[[grep("[Ff]ormula|[Vv]alue|[Kk]inetic[.][Ll]aw",colnames(df))[1]]]
	}
	if (is.null(f)) stop("data.frame contains no formula column")
	names(f) <- rownames(df)
	return(f)
}

#' Modifies a value
#'
#' This function does the same as `x <- x + v`, but without repeating x
#' modify(x) <- 1 will add one to it
#'
#' This is meant to add entries in a matrix 
#' @param x
#' @param v
#' @return `x+v`
`modify<-` <- function(x,i,j,sign=+1,value){
	x[i,j] <- x[i,j] + sign*value
	return(x)
}

#' Convert to linear space
#'
#' A number given in some logarithmic space can be transformed back to linear space
#' A call like `base(x) <- 10` means that x was provided in decadic logarithm form.
#' This will adjust `x` so that it is now in linear space.
#'
#' If `x` was provided in logarithmic space, then it is an exponent to
#' the given base (value).
#' 
#' @param x a numeric vector
#' @param i a subset of values in x, defaults to all values of x
#' @param value the base of the logarithm x was provided in
#' @return x will be changed to be in linear space
#' @example x <- 1; base(x) <- 10; print(x)
`base<-` <- function(x,i=seq_along(x),value){
	x[i] <- exp(log(value)*x[i])
	return(x)
}

stoichiometric_matrix <- function(m) {
	nu <- matrix(0,NROW(m$Compound),NROW(m$Reaction),dimnames=list(rownames(m$Compound),rownames(m$Reaction)))
	reactants <- stoichiometry(lapply(strsplit(m$Reaction$reactants,"+",fixed=TRUE),trimws))
	products <- stoichiometry(lapply(strsplit(m$Reaction$products,"+",fixed=TRUE),trimws))
	f <- m$Reaction$kinetic.law
	for (j in seq_along(reactants)){
		r <- reactants[[j]]
		p <- products[[j]]
		modify(nu,names(r),j) <- -r
		modify(nu,names(p),j) <- p
	}
	attr(nu,"reactants") <- reactants
	attr(nu,"products") <- products
	return(nu)
}

#' Fetch an Attribute
#'
#' This function differs from `rlang::%@%` in that it stops if the
#' attribute doesn't exist.
#' @param x an R object (variable with attributes)
#' @param a the name of an attribute
#' @return the value of the attribute: `attr(x,a)`
`%@%` <- function(x,a){
	stopifnot(x %has% a)
	return(attr(x,a))
}

#' Interprets a character vector as names of logarithms
#'
#' The values in `x` are possibly given in a logarithmic space. The
#' parameter `str_scale` gives this logarithmic scale (provided in a
#' language agnostic form), by a human. An empty string causes no
#' transformations. Similarly, providing no scale at all causes no
#' transformations.
#'
#' The words in `str_scale` name a logarithm, e.g. "log10". Currently understood scales:
#' - log10
#' - log2, ld
#' - ln, log
#' @param x values
#' @param str_scale character vector
#' @return a copy of `x`,  transformed in to linear space
linear_scale <- function(x,str_scale=""){
	if (is.character(str_scale)){
		# general case: logXX
		l <- grepl("^log[0-9]+$",str_scale)
		base(x,l) <- as.numeric(sub("^log([0-9]+)$","\\1",str_scale[l]))
		# some special cases:
		base(x,grepl("^log$|^ln$",str_scale)) <- exp(1)
		base(x,grepl("^ld$",str_scale)) <- 2
	}
	return(x)
}


#' Add a Reaction to an ODE
#'
#' This function operates on an ODE vector field. It makes it possible
#' to construct the right hand side one reaction at a time:
#' `reaction(vf,r,p) <- flux` will add a flux term to the ODE's right
#' hand side (vector field), adding the given `flux` (reaction
#' kinetic) in the right places, informed by the reactants and
#' products.
#'
#' The reactants, e.g.: `c(A=2,B=1)` name which state variables are affected by the flux (negatively, for reactants)
#' The products, e.g.: `c(C=1)` name which state variables are affected positively by the flux.
#' @param vf a named character vector of the right length (number of state variables)
#' @param r a named vector of stoichiometric coefficients for the reactants
#' @param p a named vector of stoichiometric coefficients for the products
#' @return updated vf
#' @example
#' Reaction <- "A + B <=> C"
#' r <- c(A=1,B=1)
#' p <- c(C=1)
#' vf <- c(A="",B="",C="") # empty
#' reaction(vf,r,p) <- "A*B-C"
`reaction<-` <- function(vf,r,p,value){
	stopifnot(r %has% "names")
	stopifnot(p %has% "names")
	stopifnot(vf %has% "names")
	if (value %has% "names") {
		f <- names(value)
	} else {
		f <- value
	}
	if (all(r==1)){
		vf[names(r)] <- sprintf("%s-%s",vf[names(r)],f)
	} else {
		vf[names(r)] <- sprintf("%s-%i*%s",vf[names(r)],r,f)
	}
	if (all(p==1)){
		vf[names(p)] <- sprintf("%s+%s",vf[names(p)],f)
	} else {
		vf[names(p)] <- sprintf("%s+%i*%s",vf[names(p)],p,f)
	}
	return(vf)
}

#' Reduce the size of the system
#'
#' Given a stoichiometric matrix, this function performs model
#' reduction via linear algebra operations, with [pracma].
#' @param nu stoichiometric matrix
#' @useDynLib uqsa, lstrtod
#' @return a list of conservation laws
conservation_law_analysis <- function(nu,iv,verbose=FALSE) {
	C <- pracma::rref(t(pracma::nullspace(t(nu))))
	stopifnot(norm(C %*% nu)<1e-6)
	nm <- rownames(nu)
	colnames(C) <- rownames(nu)
	K <- numeric(NROW(C))
	allText <- character(NROW(C))
	if (is.character(iv)) {
		warning("initial values are not numeric.")
		iv <- .Call(lstrtod,iv) # guarantees no NA values
	}
	for (i in seq(NROW(C))){
		l <- which(abs(C[i,])>1e-3)
		for (k in l){
			if (all(k!=K)){
				K[i] <- k
				break
			}
		}
		j <- l[-k]
		if (C[i,k]!=1){
			Text <- paste0(
				sprintf("(%s_ConservedConst - (",nm[k]),
				paste0(
					sprintf("%+i*%s",round(C[i,j]),nm[j]),
					collapse=""
				),
				sprintf("))*(%i)",1.0/C[i,k])
			)
		} else {
			Text <- paste0(
				sprintf("%s_ConservedConst - (",nm[k]),
				paste0(
					sprintf("%+i*%s",round(C[i,j]),nm[j]),
					collapse=""
				),
				sprintf(")")
			)
		}
		if (as.logical(verbose)){
			cat("d/dt(",sprintf("%+i*%s",round(C[i,l]),nm[l]),") == 0\n")
			cat(Text,"\n")
		}
		allText[i] <- Text
	}
	if (any(diff(c(1,2,4,11,14))==0)) stop("conservation law analysis is trying to replace the same compound several times.")
	CL <- data.frame(
		Eliminates=rev(K),
		value=rev(C %*% iv),
		Constant=rev(C %*% iv),
		ConstantName=rev(sprintf("%s_ConservedConst",nm[K])),
		Formula=rev(gsub("+1*","+",allText,fixed=TRUE)),
		row.names=rev(nm[K])
	)
	attr(CL,"lawMatrix") <- t(C)
	return(CL)
}

#' Interpret a model as an ODE
#'
#' This function accepts a list generated from a collection of TSV
#' files (or a similar format) and interprets the contents as an
#' ordinary differential equation (ODE).
#'
#' The argument `m` can be obtained via `model_from_tsv()`.
#' It has the components:
#' - `m$Constant`
#' - `m$Parameter`
#' - `m$Input`
#' - `m$Expression`
#' - `m$Compound`
#' - `m$Reaction`
#' - `m$Experiment`
#'
#' There can be additional components describing measured data for
#' this model.
#' @param m a list of data.frames, each corresponding to a TSV file or
#'     sheet in a spreadsheet.
#' @param cla a Boolean value indicating whether conservation law
#'     analysis should be performed.
#' @return a list that contains a summary of this model interpreted as
#'     an ODE, crucially, the list contains the element `vf`, the
#'     right-hand-side (vector field) of the ODE, this is the main
#'     result of this function.
#' @export
as_ode <- function(m,cla=requireNamespace("pracma")){
	iv <- values(m$Compound)
	pv <- values(m$Parameter)
	xp <- c(formulae(m$Expression),formulae(m$Reaction))
	inp <- values(m$Input)
	par <- c(linear_scale(pv,m$Parameter$scale),linear_scale(inp,m$Input$scale))
	out <- formulae(m$Output)
	nu <- stoichiometric_matrix(m)
	flux <- formulae(m$Reaction)
	vf <- character(NROW(iv))
	names(vf) <- names(iv)
	r <- nu %@% 'reactants'
	p <- nu %@% 'products'
	for (i in seq_along(flux)){
		reaction(vf,r[[i]],p[[i]]) <- flux[i]
	}
	if (as.logical(cla)){
		CL <- conservation_law_analysis(nu,iv)
		iv <- iv[-CL$Eliminates]
		cq <- CL$Formula
		names(cq) <- rownames(CL)
		xp <- c(formulae(m$Expression),cq,formulae(m$Reaction))
		C <- CL$Constant
		names(C) <- CL$ConstantName
		par <- c(par,C)
		vf <- vf[-CL$Eliminates]
	} else {
		CL <- NULL
		xp <- c(formulae(m$Expression),formulae(m$Reaction))
	}
	ode <- list(vf=vf,par=par,var=iv,exp=xp,func=out,stoichiometric_matrix=nu,conservationLaws=CL)
	comment(ode) <- comment(m)
	return(ode)
}

time_series_experiment <- function(t,measurements){

}

#'
#'
#' Updates the named values of vector `v` with values mentioned in data.frame d
#'
update_values <- function(v,d){
	if (is.null(names(v))) stop("[update_values] v must be named")
	ret <- matrix(v,length(v),NROW(d),dimnames=list(names(v),rownames(d)))
	for (i in seq(NCOL(ret))){
		j <- names(v) %in% colnames(d)
		if (any(j)){
			ret[j,i] <- d[i,names(v)[j]]
		}
	}
	return(ret)
}

#' Extract Measured Data and Simulation Experiment Instructions
#'
#' This function accepts the model obtained via `model_fromt_tsv` or a
#' similar function. It finds the data tables for this model (if any
#' are present), and finds the simulation instructions to reproduce
#' these data sets using the model.
#'
#' This function requires that the files the model is stored as
#' contains measurements (data) that can be interpreted fairly
#' easily. Each data file needs columns that are named like the
#' observable quantities listed in the Output table.
#'
#' If the data is very indirectly related to the model, then we don't
#' interpret the data files themselves and the user needs to write a
#' specialized likelihood function to relate the raw data in the files
#' with something that the model does. In such cases, don't use this function.
#'
#' The instructions must be organised in a table called Experiment(s).
#' @param m the model (with data), as obtained via `model_from_tsv()`, or similar.
#' @param CL conservation laws, obtained through conservation law analysis.
#' @return a list of simulation instructions
data_with_instructions <- function(m,o){
	if (!is.list(m)) {
		stop("the first argument needs to be a list of file contents.")
	}
	# The Experiment table is orchestrating the whole simulation procedure
	if (all(is.finite(pmatch('Experiment',names(m))))){
		E <- m$Experiment
	} else {
		cat(names(m),"\n",sep=", ")
		stop("argument must contain an item named 'Experiment(s)'.")
	}
	# If there isn't an output table, then all state variables must be measurable
	if (all(is.finite(pmatch('Output',names(m))))){
		out <- rownames(m$Output)
		print(out)
	} else {
		out <- rownames(m$Compound)
	}
	if (all(is.finite(pmatch('Input',names(m))))){
		C <- o$conservationLaws$Constant
		names(C) <- o$conservationLaws$ConstantName
		input <- update_values(
			c(values(m$Input),C),
			m$Experiment
		)
	} else {
		input <- NULL
	}
	iv <- update_values(o$var,m$Experiment)
	D <- vector("list",NROW(E))
	for (i in seq(NROW(E))){
		d <- m[[rownames(E)[i]]]
		if (is.null(d)) {
			message("no data provided for experiment ",i)
		} else {
			DATA <-parse_concise(
				t(d[,colnames(d) %in% out, drop=FALSE]),
				use.errors=TRUE
			)  # this is a matrix
			rownames(DATA) <- out
		}
		D[[i]] <- list(
			measurements=cbind(time=d$time,as.data.frame(t(DATA))),
			data=DATA,
			input=input[,i],
			initialTime=E$t0[i] %otherwise% min(d$time),
			initialState=iv[,i],
			outputTimes=d$time
		)
	}
	names(D) <- rownames(E)
	return(D)
}
