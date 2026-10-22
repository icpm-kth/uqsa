boxes <- function(txt,header=names(txt),footer=NULL,...){
	if (length(txt) == 0) return(invisible(txt))
	stopifnot(is(txt,"list"))
	content_width <- max(sapply(txt,cli::ansi_nchar))
	if (missing(footer)){
		footer <- character(length(txt))
		for (i in seq_along(footer)) {
			footer[i] <- digest::digest(txt[[i]])
		}
	} else if (is.null(footer)){
		footer <- character(length(txt))
	} else if (length(footer)<length(txt)) {
		footer <- rep_len(as.character(footer), length(txt))
	}
	B <- character(length(txt))
	for (i in seq_along(txt)){
		B[i] <- cli::boxx(
			cli::ansi_align(txt[[i]], width=content_width),
			padding=c(0,1,0,1),
			header=header[i],
			footer=footer[i],
			border_style='round',
			...
		)
	}
	padding_width <- 6
	N <- max(1,floor(cli::console_width()/(content_width+padding_width)))
	for (j in seq(1,length(B),by=N)){
		b <- B[seq(j,min(j+N-1,length(B)))]
		BB <- do.call(
			paste,
			cli::ansi_strsplit(b,'\n',fixed=TRUE)
		)
		cat(BB,sep='\n')
	}
	return(invisible(txt))
}

boxes_of_reactions <- function(x,...){
	if (is.null(x) || !is(x,"list") || length(x)==0) return(invisible(x))
	if (is.null(x$Reaction) || NROW(x$Reaction) == 0) return(invisible(x))
	MR <- max(nchar(x$Reaction$reactants))
	MP <- max(nchar(x$Reaction$products))
	K <- kinetic_law_matrix(x$Reaction)
	rev <- x$Reaction$is.reversible
	if (is.null(rev)) rev <- !logical(NROW(x$Reaction))
	leftrightharpoons <- " \u21CC "
	rightharpoon <- " \u21C0 "
	W <- cli::ansi_nchar(leftrightharpoons)
	arr <- ifelse(rev,leftrightharpoons,rightharpoon)
	rev_str <- ifelse(rev,"reversible","irreversible")
	txt <- list()
	for (i in seq(NROW(x$Reaction))){
		fwd <- paste(cli::col_blue("fwd:"),cli::col_grey(K[i,'fwd']))
		bwd <- paste(cli::col_yellow("bwd:"),cli::col_grey(K[i,'bwd']))
		# take care of possible synonyms of EMPTYSET
		r <- gsub("\\b(NULL|Nil|\u00D8|\u2205x)\\b","NONE",x$Reaction$reactants[i])
		p <- gsub("\\b(NULL|Nil|\u00D8|\u2205x)\\b","NONE",x$Reaction$products[i])
		r_or_none <- ifelse(nchar(r)>0,r,"NONE")
		p_or_none <- ifelse(nchar(p)>0,p,"NONE")
		rs <- sprintf("%*s%s%s",MR,r_or_none,arr[i],p_or_none)
		rs <- gsub("\\bNONE\\b",cli::col_grey("NONE"),rs)
		txt[[i]] <- c(
			cli::ansi_align(rs,width=MR+MP+W),
			fwd,
			bwd,
			...
		)
	}
	names(txt) <- rownames(x$Reaction)
	boxes(txt,footer=rev_str,border_col='green')
	return(invisible(x))
}

boxes_of_experiments <- function(x,...){
	if (is.null(x) || !is(x,"list") || length(x)==0) return(invisible(x))
	if (is.null(x$Experiment) || NROW(x$Experiment) == 0) return(invisible(x))
	type <- x$Experiment$type
	if ("t0" %in% names(x$Experiment)) {
		t0 <- x$Experiment$t0
	} else {
		t0 <- numeric(NROW(x$Experiment))
	}
	if ("event" %in% names(x$Experiment)) {
		ev <- x$Experiment$event
	} else {
		ev <- character(NROW(x$Experiment))
	}
	txt <- list()
	type_str <- rep("time series",NROW(x$Experiment))
	for (i in seq(NROW(x$Experiment))){
		l <- colnames(x[[rownames(x$Experiment)[i]]]) %in% rownames(x$Output)
		if (any(l)){
			measured <- cli::col_grey(
				paste(
					colnames(x[[rownames(x$Experiment)[i]]])[l],
					collapse=", "
				)
			)
		} else {
			measured <- cli::col_red("NA")
		}
		DATA <- x[[rownames(x$Experiment)[i]]]
		if (is.null(DATA)) {
			d <- cli::col_red("none")
		} else {
			d <- cli::style_bold(cli::col_green(NROW(DATA)))
		}
		txt[[i]] <- c(
			paste(cli::col_magenta("data rows:"),d),
			paste(cli::col_cyan("t0:"),t0[i]),
			paste(cli::col_blue("events:"),ev[i]),
			paste(cli::col_yellow("measured:"),measured),
			...
		)
		if (!is.null(type) && isTRUE(grepl("^ *[Dd]ose[ .]?([Rr]esponse)? *$",type[i]))){
			type_str[i] <- "dose response"
		}
	}
	names(txt) <- rownames(x$Experiment)
	boxes(txt,footer=type_str,border_col='yellow')
	return(invisible(x))
}

#' Prints an overview of the model
#'
#' Print an overview of the model, arranged as boxes on screen.  The
#' primary information is the reactions within the model.  When
#' verbose output is required, boxes representing each simulation
#' experiment and data-set that corresponds to it are displayed as
#' well.
#'
#' @export
#' @param x the model (class 'sbmodel')
#' @param verbose whether to print information about the experiments and data-sets
#' @param ... ignored
#' @return the model, invisibly, called for side-effect
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' print(m)
print.sbmodel <- function(x,verbose=getOption("uqsa.verbose",TRUE),...){
	if (is.null(x) || !is(x,"list") || length(x)==0) return(invisible(x))
	d <- cli::cli_div(
		theme=list(
			rule=list("line-type" = "double", "color" = "cyan")
		)
	)
	cli::cli_rule(cli::style_bold(cli::col_br_white(comment(x))))
	cli::cli_end(d)
	n <- c(
		"constant"=NROW(x$Constant),
		"input parameters"=NROW(x$Input),
		"parameters"=NROW(x$Parameter),
		"reacting compounds"=NROW(x$Compound),
		"reactions"=NROW(x$Reaction),
		"experiments"=NROW(x$Experiment)
	)
	cat(sprintf("%*s : %i",max(nchar(names(n))),names(n),n),sep='\n')
	cli::cli_h1("Reactions")
	boxes_of_reactions(x)
	if (as.logical(verbose)) {
		cli::cli_h1("Experiments")
		boxes_of_experiments(x)
	}
	return(invisible(x))
}

#' Prints an interpretation string of a unit
#'
#' The unit object is a tagged data frame, with these columns:
#' - multiplier
#' - kind
#' - scale
#' - exponent
#'
#' The interpretation is the same as in SBML units. This function also
#' prints an inferred unit id: a string that has no special characters
#' in it and can be used in places where such characters are not
#' allowed (e.g. SBML unit `id` attribute).
#'
#' The original string that a unit was derived from is attached to the
#' unit object as a [comment].
#'
#' Units are produced by the function [unit.from.string].
#'
#' @keywords internal
#' @param x an object of type 'unit_of_measurement'
#' @param ... required by the generic print function.
#' @export
#' @return called for the side-effect; no value.
#' @examples
#' lapply(lapply(c("km/h","s^-2","1/s"),unit.from.string),print)
print.unit_of_measurement <- function(x,...){
	unit <- x
	cat(
		sprintf("\u00ab%s\u00bb has been interpreted as",comment(unit)),
		"the product of: ",
		sprintf("%30s",attr(unit,"id")),
		"==============================",
		sprintf(
			"(%g \u00D7 %s \u00D7 10^(%i))^(%i)",
			unit$multiplier,
			unit$kind,
			unit$scale,
			unit$exponent
		),
		sep='\n'
	)
	invisible(x)
}

#' Prints the simulation experiments
#'
#' The experiments, if accidentally printed, are difficult to read.
#' This function prevents these accidental prints. It summarizes the
#' data and simulation experiments instead.
#'
#' @keywords internal
#' @param x simulation experiments with data
#' @param ... ignored.
#' @export
#' @return called for side-effect (printout); no value.
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- as_ode(m)
#' ex <- experiments(m,o)
#' print(ex)
print.experiments <- function(x,...){
	ex <- x
	cat(sprintf("number of simulation experiments: %i\n",length(ex)))
	for (i in seq_along(ex)){
		cat(sprintf("%42s",names(ex)[i]),"\n")
		cat(paste0(rep("-",42),collapse=""),"\n")
		for (j in seq_along(ex[[i]])){
			x <- ex[[i]][[j]]
			if (is.array(x)){
				cat(sprintf("%24s: %s (dim)\n",names(ex[[i]])[j],paste(dim(x),collapse=", ")))
			} else if (is.data.frame(x)){
				cat(
					sprintf(
						"%24s: %i columns (%s)\n",
						names(ex[[i]])[j],
						NCOL(x),
						paste(class(x),collapse=", ")
					)
				)
			} else if (is.matrix(x)){
				cat(sprintf("%24s: %i\u00D7%i (dim)\n",names(ex[[i]])[j],NROW(x),NCOL(x)))
			} else if (is.numeric(x) && length(x)==1) {
				cat(sprintf("%24s: %g\n",names(ex[[i]])[j],x))
			} else if (is.numeric(x)) {
				cat(sprintf("%24s: %i (length)\n",names(ex[[i]])[j],length(x)))
			} else {
				cat(
					sprintf(
						"%24s: %s (class), %s (type)\n",
						names(ex[[i]])[j],
						paste(class(x),collapse=", "),
						typeof(x)
					)
				)
			}
		}
		cat("\n")
	}
	cat("experiments: ",paste(names(ex),collapse=", "),"\n")
}

#' Prints the simulation results
#'
#' The results, if accidentally printed, are difficult to read.  This
#' function prevents these accidental prints. It summarizes the
#' results instead.
#'
#' @keywords internal
#' @param x simulation results
#' @param ... requirement of print generic, not used.
#' @export
#' @return called for side-effect (printout); no value.
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- as_ode(m)
#' ex <- experiments(m,o)
#' c_path(o) <- write_c_code(generate_code(o))
#' so_path(o) <- shlib(o)
#' s <- simfi(ex,o)
#' y <- s(values(m$Parameter))
#' print(y)
print.simulation <- function(x,...){
	y <- x
	cat(sprintf("number of simulation experiments: %i\n",length(y)))
	for (i in seq_along(y)){
		cat(sprintf("%42s",names(y)[i]),"\n")
		cat(paste0(rep("-",42),collapse=""),"\n")
		for (j in seq_along(y[[i]])){
			x <- y[[i]][[j]]
			if (is.array(x)){
				cat(sprintf("%24s: %s (dim)\n",names(y[[i]])[j],paste(dim(x),collapse=", ")))
			} else if (is.numeric(x) && length(x)==1) {
				cat(sprintf("%24s: %g\n",names(y[[i]])[j],x))
			} else {
				cat(
					sprintf(
						"%24s: %s (class), %s (type)\n",
						names(y[[i]])[j],
						paste(class(x),collapse=", "),
						typeof(x)
					)
				)
			}
		}
		cat("\n")
	}
	cat("experiments: ",paste(names(y),collapse=", "),"\n")
}


#' Prints a summary about the ode
#'
#' An ODE model was created by `as_ode` can be summarized here,
#' including information about the compiled version of the model.
#'
#' The ode model is for the most part a list of named vectors and
#' matrices which together encode the mathematical structure of the
#' ode.
#'
#' @keywords internal
#' @param x the ode
#' @param ... requirement of print generic, not used.
#' @return NULL
#' @export
#' @return called for side-effect (printout); no value.
#' @examples
#' f <- uqsa_example("AKAR4")
#' m <- model_from_tsv(f)
#' o <- as_ode(m)
#' print(o)
print.ode <- function(x,...){
	o <- x
	cat(
		sprintf("%26s : %s","Model name",o$name),
		sprintf("%26s : %s [%s]","C file",o$c_path,o$c.date),
		sprintf("%26s : %s [%s]","shared library",o$so_path,o$c.date),
		sprintf("%26s : %i","Number of state variables",length(o$var)),
		sprintf("%26s : %i","Number of parameters",length(o$par)),
		sprintf("%26s : %i","Number of outputs",length(o$func)),
		sprintf("%26s : %i","Conservation laws",NROW(o$conservationLaws)),
		sprintf("%26s : %s","Transformations",ifelse(is.null(o$tf),"no","yes")),
		sep="\n"
	)
}

#' Prints information about the mcmc variable
#'
#' Some mcmc variables have many attributes, which clutter the screen
#' when accidentally printed. This function prevents these long
#' printouts.
#'
#' @keywords internal
#' @param x the variable
#' @param ... requirement of print generic, not used.
#' @export
#' @return called for side-effect (printout); no value.
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- as_ode(m)
#' c_path(o) <- write_c_code(generate_code(o))
#' so_path(o) <- shlib(o)
#' ex <- experiments(m,o)
#' dprior <- dNormalPrior(values(m$Parameter),m$Parameter$stdv)
#' s <- simfi(ex,o)
#' p <- mcmc_init(1.0,values(m$Parameter),s,dprior=dprior)
#' print(p)
print.mcmcVariable <- function(x,...){
	v <- x
	print(v[seq_along(v)])
	A <- c("simulations", "logLikelihood", "prior", "gradLogLikelihood","gradLogPrior","fisherInformation")
	for (a in A){
		if (v %has% a){
			if (is.list(attr(v,a))){
				cat(sprintf("%24s: %i (length)\n",a,length(attr(v,a))))
			} else if (is.numeric(attr(v,a)) && length(attr(v,a))==1){
				cat(sprintf("%24s: %g\n",a,attr(v,a)))
			} else if (is.numeric(attr(v,a)) && length(attr(v,a))>1){
				cat(sprintf("%24s: %i (length)\n",a,length(attr(v,a))))
			} else if (is.matrix(attr(v,a))) {
				cat(sprintf("%24s: %s (dim)\n",a,paste(dim(attr(v,a)),collapse="x")))
			} else {
				cat(
					sprintf(
						"%24s: %s (class) %s (type)\n",
						a,
							paste(class(attr(v,a)),collapse=", "),
							typeof(attr(v,a))
					)
				)
			}
		}
	}
}

#' Prints a Summary about the CME model
#'
#' This information printed on screen omits the details about the
#' interactions, only the lengths of the vectors included in the data
#' structure `CME`.
#'
#' @keywords internal
#' @param x a model created by [as_cme]
#' @param ... requirement of print generic, not used.
#' @return Nil
#' @export
#' @examples
#' f <- uqsa_example("AKAR4")
#' m <- model_from_tsv(f)
#' cmeModel <- as_cme(m)
#' print(cmeModel)
print.cme <- function(x,...){
	cmeModel <- x
	cat(
		sprintf("%26s : %s","Name",cmeModel$name),
		sprintf("%26s : %s [%s]","C file",cmeModel$c_path,cmeModel$c.date),
		sprintf("%26s : %s [%s]","shared library",cmeModel$so_path,cmeModel$so.date),
		sprintf("%26s : %i","Number of state variables",length(cmeModel$initialCount)),
		sprintf("%26s : %i","Number of parameters",length(cmeModel$par)),
		sprintf("%26s : %i","Number of outputs",length(cmeModel$output)),
		sprintf("%26s : %i","Number of constants",length(cmeModel$const)),
		sep="\n"
	)
}
