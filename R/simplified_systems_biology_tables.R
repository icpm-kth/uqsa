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

#' Interpret a model as an ODE
as_ode <- function(m){
	iv <- values(m$Compound)
	pv <- values(m$Parameter)
	xp <- c(formulae(m$Expression),formulae(m$Reaction))
	inp <- values(m$Input)
	par <- c(linear_scale(pv,m$Parameter$scale),linear_scale(inp,m$Input$scale))
	out <- formulae(m$Output)
	nu <- stoichiometric_matrix(m)
	flux <- formulae(m$Reaction)
	vf <- character(NROW(iv))
	names(vf) <- rownames(m$Compound)
	r <- nu %@% 'reactants'
	p <- nu %@% 'products'
	for (i in seq_along(flux)){
		reaction(vf,r[[i]],p[[i]]) <- flux[i]
	}
	ode <- list(vf=vf,par=par,var=iv,exp=xp,func=out,stoichiometric_matrix=nu)
	return(ode)
}
