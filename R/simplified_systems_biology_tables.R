## Here we define functions that are meant to replace SBtab, by making
## it as simple as possible, with no need to remember anything.

#' Standard Error Matrix from an errors object
#'
#' If a matrix has an `errors` attribute, it is usually a vector.
#' THis function returns the values of this attribute as a matrix (it
#' preserves the dimensions of the host matrix).
#'
#' @param M a matrix with errors (uncertainties)
#' @export
#' @return A matrix similar to E, with standard error values
#' @examples
#' M <- matrix(seq(12),3,4,dimnames=liest(letters[seq(3)],LETTERS[seq(4)]))
#' errors(M) <- abs(M*0.1 + 0.1)
#' E <- standard_error_matrix(M)
#' print(E)
standard_error_matrix <- function(M){
	if (is.matrix(M) && is(M,"errors")){
		d <- dim(M)
		E <- errors::errors(M)
		dim(E) <- d
		dimnames(E) <- dimnames(M)
	} else {
		E <- NULL
	}
	return(E)
}

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

#' Get j-th column with names
#'
#' When indexing a matrix, rownames are lost. This function will
#' return a column of a matrix, as a vector (dropping rank so to
#' speak), but the vector will retain the rownames of the matrix
#' @param m a matrix
#' @param j a column index
#' @return a named vector
column <- function(m,j=1){
	v <- m[,j]
	names(v) <- rownames(m)
	return(v)
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
#' @return a named numeric vector
#' @export
values <- function(df){
	if (is.null(df)) return(NULL)
	else stopifnot(is.data.frame(df))
	##
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
	attr(v,"unit") <- units_from_table(df)
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
	attr(f,"unit") <- units_from_table(df)
	return(f)
}

#' Modifies a value
#'
#' This function does the same as `x <- x + sign*value`, but without
#' repeating `x`. The expression `modify(x) <- rnorm(length(x),0,1)`
#' will add Gaussian noise to it. This is meant as a replacement for
#' the `x += 1` syntax of C, it exists only for aesthetic reasons.
#'
#' Specifically, this function should work for matrices, and it is
#' possible to supply row and column index vectors: x[i,j] will be
#' modified.
#'
#' This function is quite useful if `x` has a very long name,
#' e.g. `experiments[[1]]$func`.
#'
#' @param x a numeric value to be modified
#' @param i row-indices of `x` to be modified
#' @param j column-indices of `x` to be modified
#' @param sgn modification
#' @return The value of x is modified additively ihn place: `x <- x+value`
#' @export
#' @examples
#'  x <- matrix(seq(12),3,4)
#'  modify(x,seq(2),seq(2)) <- 10
#'  print(x)
`modify<-` <- function(x,i=seq_along(NROW(x)),j=seq_along(NCOL(x)),sgn=+1,value){
	if (is.matrix(x)){
		x[i,j] <- x[i,j] + sgn*value
	} else {
		x[i] <- x[i] + sgn*value
	}
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
#' @examples
#'  x <- 1
#'  base(x) <- 10
#'  print(x)
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
#'
#' This function tries to find a similarly named attribute
#' disregarding capitalization and using partila matching.
#'
#' The only way from this function to return NULL is when `x` is null
#' (the object that supposedly has the attribute). For the purposes of
#' this function , NULL objects are treated as optional things, and
#' thus their attributes do not matter. Non-NULL objects that should
#' have an attribute, but don't are considered erroneous.
#'
#' @param x an R object (variable with attributes)
#' @param a the name of an attribute
#' @return the value of the attribute: `attr(x,a)`
#' @export
`%@%` <- function(x,a){
	if (is.null(x)) return(NULL)
	stopifnot(is.character(a))
	if (x %has% a){
		return(attr(x,a))
	} else if (any(is.finite(pmatch(tolower(a),tolower(names(attributes(x))))))){
		return(attributes(x)[[pmatch(tolower(a),tolower(names(attributes(x))))]])
	} else {
		stop(sprintf("Attribute «%s» not found in the list of objects's attributes: %s.\n",a,paste0(names(attributes(x)),collapse=", ")))
	}
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
#' @examples
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
#' reduction via linear algebra operations, with [pracma::null].
#' @param nu stoichiometric matrix
#' @useDynLib uqsa, lstrtod
#' @return a list of conservation laws
#' @export
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
	if (is.finite(pmatch("Transformation",names(m)))){
		tf <- character(length(iv)+length(par))
		names(tf) <- c(names(par),names(iv))
		tf <- update_values(tf,m$Transformation,"character")
		attr(tf,"type") <- c(rep("par",length(par)),rep("var",length(iv)))
	} else {
		tf <- NULL
	}
	ode <- list(vf=vf,par=par,var=iv,exp=xp,func=out,stoichiometric_matrix=nu,conservationLaws=CL,tf=tf)
	comment(ode) <- comment(m)
	return(ode)
}

#' Updates the named values of vector `v` with values mentioned in data.frame d
#'
#' Given a named vector, this function will create a mtarix, with
#' several copies of the vector, with different values, dependeing on
#' some context. The context is each row of a `data.frame`, with
#' columns referring to the names of `v`, and values for the entries
#' in `v` for that row's context. Example: the rows can be different
#' experinemnts, with each experiment assigning new values to some
#' members of `v` (but not necessarily all).
#' @param v a named vector
#' @param d data.frame with column names that correspond to those of `v`
#' @param as_type a character scalar indicating a type ('character','numeric','logical',etc.)
#' @return a matrix of dimension length(v) × NROW(d)
update_values <- function(v,d,as_type="numeric"){
	if (is.matrix(v) && NCOL(v)==1) {
		v <- column(v,1)
	}
	if (is.null(names(v))) stop("[update_values] v must be named")
	if (is.null(d) || prod(dim(d))==0) stop("non-empty data.frame d is mandatory")
	ret <- matrix(v,length(v),NROW(d),dimnames=list(names(v),rownames(d)))
	for (i in seq(NCOL(ret))){
		j <- names(v) %in% colnames(d)
		if (any(j)){
			ret[j,i] <- as(d[i,names(v)[j]],as_type)
		}
	}
	return(ret)
}

time_series_experiments <- function(m,E,iv,input,out){
	if (is.null(E) || NROW(E)==0) return(NULL)
	D <- vector("list",NROW(E))
	eventSchedule <- E$event
	tr <- m$Transformation
	if (is.null(tr) && !is.null(eventSchedule)){
		stop("When an 'event' column is present in the table of experiments (Experiment.tsv), ",
			 "then a transformation table must exist (named «Transformation.tsv») as well.")
	}
	for (i in seq(NROW(E))){
		d <- m[[rownames(E)[i]]]
		ev <- m[[eventSchedule[i]]]
		if (is.null(d)) {
			message("no data provided for experiment ",rownames(E)[i])
			DATA <- NULL
		} else {
			DATA <- parse_concise(
				t(d[,colnames(d) %in% out, drop=FALSE]),
				use.errors=TRUE
			)  # this is a matrix
			rownames(DATA) <- out
		}
		if (is.null(m[[eventSchedule[i]]])){
			event_list <- NULL
		} else {
			stopifnot(!is.null(tr))
			event_list <- list(
				time=as.double(ev$time),
				label=as.integer(match(ev$transformation,rownames(tr))-1),
				dose=as.double(ev$dose)
			)
		}
		D[[i]] <- list(
			measurements=cbind(time=d$time,as.data.frame(t(DATA))),
			data=DATA,
			input=as.double(input[,i]),
			initialTime=as.double(E$t0[i] %otherwise% min(d$time)),
			initialState=iv[,i],
			outputTimes=as.double(d$time),
			events=event_list
		)
	}
	names(D) <- rownames(E)
	return(D)
}

## dose_response_experiments(m,E[l,,drop=FALSE],iv[,l,drop=FALSE],input[,l,drop=FALSE],out)
dose_response_experiments <- function(m,E,iv,input,out){
	if (is.null(E) || NROW(E)==0) return(NULL)
	TS <- list() # list of time series experiments
	t0 <- E$t0
	tf <- E$tf
	if (!is.null(E$event) && any(nzchar(E$event))){
		print(E[,c("type","event")])
		warning("Dose response experiments are not (yet) compatible with events.")
	}
	for (i in seq(NROW(E))){
		d <- m[[rownames(E)[i]]] # data table
		ts <- vector("list",length=NROW(d))
		DATA <-parse_concise(
			t(d[,colnames(d) %in% out, drop=FALSE]),
			use.errors=TRUE
		)  # this is a matrix
		rownames(DATA) <- out
		for (j in seq_along(ts)){
			u <- input[,i]
			names(u) <- rownames(input)
			inputMatrix <- update_values(u,d)
			x <- iv[,i]
			names(x) <- rownames(iv)
			initialStateMatrix <- update_values(x,d)
			ts[[j]] <- list(
				outputTimes=as.double(tf[i] %otherwise% d$time[j]),
				measurements=d[j,,drop=FALSE],
				data=DATA[,j,drop=FALSE],
				input=as.double(column(inputMatrix,j)),
				initialState=as.double(column(initialStateMatrix,j)),
				initialTime=as.double(E$t0[i])
			) # several time series experiments per 1 dose response table
		}
		names(ts) <- sprintf("%s_dose_%i",rownames(E)[i],seq_along(ts))
		TS <- c(TS,ts)
	}
	return(TS)
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
#' @param o the ode derived from `m`
#' @return a list of simulation instructions
#' @export
experiments <- function(m,o){
	if (!is.list(m)) {
		stop("the first argument needs to be a list of file contents.")
	}
	if (is.finite(pmatch("Transformation",names(m))) && !is.null(o$conservationLaws)){
		warning(
			"CONFLICT: This model seems to have event based transformations and conservation laws. ",
			"These two concepts clash with one another if a compound is conserved, ",
			"but also changed by scheduled events."
		)
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
	} else {
		out <- rownames(m$Compound)
	}
	if (all(is.finite(pmatch('Input',names(m))))){
		if (is.null(o$conservationLaws)){
			conservedConstants <- NULL
		} else {
			C <- o$conservationLaws$Constant
			names(C) <- rownames(o$conservationLaws)
			conservedConstants <- update_values(C,E)
			rownames(conservedConstants) <- o$conservationLaws$ConstantName
		}
		input <- rbind(
			update_values(values(m$Input),E),
			conservedConstants
		)
	} else {
		input <- NULL
	}
	iv <- update_values(o$var,E)
	if ("type" %in% colnames(E)){
		l <- grepl("[Dd]ose[- ]?[Rr]esponse",E$type)
	} else {
		l <- logical(NROW(E))
	}
	D <- c(
		time_series_experiments(m,E[!l,,drop=FALSE],iv[,!l,drop=FALSE],input[,!l,drop=FALSE],out),
		dose_response_experiments(m,E[l,,drop=FALSE],iv[,l,drop=FALSE],input[,l,drop=FALSE],out)
	)
	return(D)
}
