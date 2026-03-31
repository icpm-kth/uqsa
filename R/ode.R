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
#' @examples
#' f <- uqsa_example("AKAR4")
#' m <- model_from_tsv(f)
#' o <- as_ode(m)
#' print(names(o))
#' print(o$vf)
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
		IV <- update_values(iv,m$Experiment)
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
	ode <- list(
		vf=vf,
		par=par,
		var=iv,
		exp=xp,
		func=out,
		stoichiometric_matrix=nu,
		conservationLaws=CL,
		tf=tf,
		name=comment(m),
		c.path=NULL,  # path of the .c file
		c.date=NULL,  # date and time of writing the c file
		so.path=NULL, # path to the so.file
		so.date=NULL  # date and time of writing the so file
	)
	class(ode) <- "ode"
	comment(ode) <- comment(m)
	return(ode)
}

#' Print a summary about the ode
#'
#' An ODE model was crteated by `as_ode` can be summarized here,
#' including information about the compiled version of the model.
#'
#' The ode model is for the most part a list of named vectors and
#' matrices which together encode the mathematical structure of the
#' ode.
#'
#' @param o the ode
#' @return NULL
#' @export
#' @examples
#' \donttest{
#'   f <- uqsa_example("AKAR4")
#'   m <- model_from_tsv(f)
#'   o <- as_ode(m)
#'   print(o)
#' }
print.ode <- function(o){
	cat(
		sprintf("%26s : %s","Model name",o$name),
		sprintf("%26s : %s [%s]","C file",o$c.path,o$c.date),
		sprintf("%26s : %s [%s]","shared library",o$so.path,o$c.date),
		sprintf("%26s : %i","Number of state variables",length(o$var)),
		sprintf("%26s : %i","Number of parameters",length(o$par)),
		sprintf("%26s : %i","Number of outputs",length(o$func)),
		sprintf("%26s : %i","Conservation laws",NROW(o$conservationLaws)),
		sprintf("%26s : %s","Transformations",ifelse(is.null(o$tf),"no","yes")),
		sep="\n"
	)
}

#' Add information about compiled code
#'
#' Adds the path of the shared library (.so file) to the ODE model.
#'
#' @param o the ode (list of named arrays and matrices)
#' @param value the path to the compiled model
#' @return modified o, with information about compiled code
`so.path<-` <- function(o,value){
	if (!is.character(value) || length(value)>1) stop("Value must be a character scalar")
	if (!file.exists(value)) warning(sprintf("File %s does not exist",value))
	o$so.path <- value
	o$so.date <- file.info(value)$ctime
	return(o)
}

#' Add information about the model's C code
#'
#' Adds the location of the model's C code (a file).
#'
#' @param o the ode (list of named arrays and matrices)
#' @param value the path to the compiled model
#' @return modified o, with information about compiled code
`c.path<-` <- function(o,value){
	if (!is.character(value) || length(value)>1) stop("Value must be a character scalar")
	if (!file.exists(value)) warning(sprintf("File %s does not exist",value))
	o$c.path <- value
	o$c.date <- file.info(value)$ctime
	return(o)
}

#' Add information about compiled code
#'
#' Adds the path of the shared library (.so file) to the ODE model.
#'
#' @param o the ode (list of named arrays and matrices)
#' @return modified o, with information about compiled code
so.path <- function(o){
	f <- o$so.path
	if (!file.exists(f)) warning(sprintf("File %s does not exist anymore.",f))
	return(f)
}

#' Retrieve information about the model's C code
#'
#' Returns the location of the model's C code (a file).
#'
#' @param o the ode (list of named arrays and matrices)
#' @return the path where the c code is stored
c.path <- function(o){
	f <- o$c.path
	if (!file.exists(f)) warning(sprintf("File %s does not exist anymore",f))
	return(f)
}

#' Compile C code to shared library
#'
#' Calls `R CMD SHLIB` to create the model's shared library.
#'
#' @param file the c file that is to be compiled, OR an ode object
#'     with a c.file defined and recorded in it.
shlib <- function(file){
	if (is(file,"ode")) {
		so_name <- file$name
		file <- c.path(file)
	} else {
		so_name <- sub("[.]c$","",basename(file))
	}
	if (is.null(file)) stop("No C file specified")
	if (!file.exists(file)) stop(sprintf("%s does not exist.",file))

	so <- file.path(dirname(file), paste0(so_name, .Platform$dynlib.ext))
	cflags <- system2("pkg-config", c("--cflags", "gsl"), stdout = TRUE)
	libs <- system2("pkg-config", c("--libs", "gsl"), stdout = TRUE)
	compile_env <- c(
		if (nzchar(cflags)) sprintf("PKG_CPPFLAGS='%s'",cflags),
		if (nzchar(libs)) sprintf("PKG_LIBS='%s'", libs)
	)
	print(compile_env)
	status <- system2(
		command = file.path(R.home("bin"), "R"),
		args = c("CMD", "SHLIB",file),
		env = compile_env,
		stdout = TRUE,
		stderr = TRUE,
		wait = TRUE
	)
	if (!file.exists(so)) {
		cat(status)
		warning(sprintf("Building %s failed.",so))
	}
	return(so)
}


#' Write the C code to a file
#'
#' This function does not compile the code, it only writes it to a
#' file in a temporary location (tempdir). By default, the name of the
#' file will contain the hash of the entire code.
#'
#' @param C the code to write, as a character array.
#' @param model.name a string with no special characters, will be used in the file name
#' @param file override the default file name (based on hashing)
#' @return the path of the written file
write_c_code <- function(C, model.name=comment(C), file=file.path(tempdir(),digest::digest(C,"xxh3_64"),paste0(model.name,".c"))){
	cat(sprintf("Writing file: %s\n",file))
	if (!dir.exists(dirname(file))){
		dir.create(dirname(file),recursive=TRUE)
	}
	cat(C,sep="\n",file=file)
	return(file)
}

#' yacasMath converts math to Ryacas compatible math
#'
#' Given a string like `"exp(2*x)"` this function returns a string that
#' yacas can process: `"Exp(2*x)"`
#'
#' @param v a chcarcter vector with math expressions
#' @param reverse do the reverse operation
#' @return a character vector with yacas math
#' @noRd
yacasMath <- function(v,reverse=FALSE) {
	r <- c("exp","sin","cos","tan","log")
	y <- c("Exp","Sin","Cos","Tan","Ln")
	if (reverse){
		for (i in seq_along(r)){
			v <- gsub(sprintf("\\b%s\\(",y[i]),sprintf("%s(",r[i]),v)
		}
		v <- gsub("UNDERSCORE","_",v)
	} else {
		for (i in seq_along(r)){
			v <- gsub(sprintf("\\b%s\\(",r[i]),sprintf("%s(",y[i]),v)
		}
		v <- gsub("_","UNDERSCORE",v)
	}
	return(v)
}

#' Jacobian of string-math
#'
#' Given a named character array of math expressions and a vector of
#' independent variables, this function calculates the Jacobian matrix
#' of the math expressions with respect to the variables.
#'
#' @param f a character vector of length n
#' @param x a character vector of length m
#' @return a character matrix (n×m) with derivatives `df[i]/dx[j]`
#' @export
#' @examples
#' f <- c("2*x*y","exp(-k*x)")
#' x <- c("x","y")
#' J <- yJacobian(f,x)
#' print(J)
yJacobian <- function(f,x){
	f <- yacasMath(f)
	x <- yacasMath(x)
	F <- paste(f,collapse=',')
	X <- paste(x,collapse=',')
	J <- Ryacas::yac_str(Ryacas::y_fn(sprintf("{%s},{%s}",F,X),"JacobianMatrix"))
	J <- gsub("[{}]","",J)
	J <- matrix(unlist(strsplit(J,",")),nrow=length(f),ncol=length(x),byrow=TRUE)
	return(yacasMath(J,rev=TRUE))
}

#' replace_powers
#'
#' @export
#' @param v a character vector
#' @useDynLib uqsa, replace_pow=replace_pow
#' @examples
#' print(replace_powers(c("2^3.1","10^-6","x^2")))
replace_powers <- function(v){
	w<-.Call(replace_pow,as.character(v))
	return(w)
}

writeComment <- \(x) {
	n <- pmax(1,1+max(nchar(x))-nchar(x))
	return(sprintf("/* %s %*s */",x,n," "))
}

ch <- \(d) {
	w <- as.character(d[[2]])
	names(w) <- d[[1]]
	return(w)
}

cOffset <- \(d) {seq(0,NROW(d)-1)}

#' Write C Function
#'
#' This function writes a C function according to a typical template.
#'
#' @param prefix function's prefix (name of model)
#' @param fNname function's name
#' @param defs optional expressions, parameters, etc.
#' @param retValue name of return value (attr: memset?)
#' @param values a character vecotor with values or mathematic formulae to return.
#' @param body a character vector of additional content for the body of the function
#' @param otherArgs additional arguments
#' @noRd
#' @return a character array that includes the function
writeCFunction <- function(prefix, fName, defArgs=c("double t","const double y_[]"), retValue=NULL, defs=NULL, values=NULL, body=NULL, otherArgs="void *par", init0 = TRUE){
	#if (is.null(values)) return(character(0))

	if (missing(otherArgs)) {
		pDefinition <- sprintf("\tdouble *p_=par;")
	} else {
		pDefinition <- NULL
	}
	named <- \(x) "names" %in% names(attributes(x))
	if (length(retValue)>0) {
		ret <- paste("double *",retValue,sep="", collapse=", ",recycle0=TRUE)
	} else {
		ret <- NULL
	}
	arguments <- paste(c(defArgs,ret,otherArgs),sep=", ",collapse=", ")
	if (!is.null(retValue) && grepl("\\by_",arguments)) {
		errValue <- sprintf("\tif (!y_ || !%s) return %i;",retValue[1],length(values))
	} else {
		errValue <- switch(
			fName,
			default=sprintf("\tif (!p_) return numParam;"),
			event=sprintf("\tif (!y_ || EventLabel<0) return numEvents;"),
			NULL
		)
	}
	if (init0) {
		zeroBuffer<-sprintf("\tmemset(%s,0,sizeof(double)*%i); /* initialize with 0.0 */",retValue[1],length(values))
	} else {
		zeroBuffer <- NULL
	}
	# values to write into return buffer:
	z <- grepl('^0$',values)
	if ("stride" %in% names(attributes(values))){
		k <- cOffset(values)
		stride <- attr(values,"stride")
		i <- k %/% stride
		j <- k %% stride
		v <- sprintf("\t/*[%2i,%2i]*/  %s[%i] = %s;",i[!z],j[!z],retValue[1],k[!z],values[!z])
	} else if (named(values)){
		v <- sprintf("\t%s[_%s] = %s;",retValue[1],names(values)[!z],values[!z])
	} else {
		v <- sprintf("\t%s[%i] = %s;",retValue[1],cOffset(values)[!z],values[!z])
	}
	return(c(
		sprintf("int %s_%s(%s){",prefix,fName,arguments),
		pDefinition,# double *p_=par;
		errValue,   # early return
		defs,       # assign local variables
		zeroBuffer, # initialize return buffer
		v,          # the values
		body,       # additional body content
		sprintf("\treturn GSL_SUCCESS;\n}")
	))
}

#' blank is TRUE where a character array is the empty string
#'
#' This function avoid double negatives in code, because nzchar stands
#' for non-zero-char, and thus not-non-zero
#'
#' ```
#' if (!nzchar(x))
#' ```
#' becomes
#' ```
#' if (blank(x))
#' ```
#'
#' @param ch a character vector
#' @return logcal array indicating where the array is equal to `""`
#' @noRd
#' @examples
#' print(blank(c("a","","","","c")))
blank <- function(ch){
	return(!nzchar(ch))
}

eventCode <-function(odeModel){
	if ("tf" %in% names(odeModel)){
		C <- c(C,"",writeComment(c("Scheduled Events","EventLabel specifies which transformation to apply.","The scalar dose variable can be used in the transformation (on the right)")))
		affectedVector <- sapply(attr(odeModel$tf,"type"),
			\(x) switch(x,var="y",par="p")
		)
		effect <- lapply(seq(NCOL(odeModel$tf)),\(i){
			ev <- odeModel$tf[,i]
			trivial <- (names(ev) == ev) | blank(ev)
			return(
				c(
					sprintf("\tcase %s:",colnames(odeModel$tf)[i]),
					sprintf("\t\t%s_[_%s] = %s;",affectedVector[!trivial],names(ev)[!trivial],ev[!trivial]),
					sprintf("\tbreak;")
				)
			)
		})
		evsw <- c(
		 sprintf("\tswitch(EventLabel){"),
		 unlist(effect),
		 "\t}"
		)
	} else {
		evsw <- NULL
	}
	return(evsw)
}

#' Write C code
#'
#' This function expects a list of character vectors, as returned by
#' `as_ode()`. This list describes an ODE model (initial values,
#' default parameters, transformation events, output functions).  This
#' function uses this information, calculates Jacobians via Ryacas and
#' returns a character vector with C source code for the solvers in
#' the GNU Scientific Library (GSL).
#'
#' The value can be written to a file:
#' `cat(generateCode(odeModel),sep="\n",file=...)`. This file can be
#' compiled into a shared library.
#'
#' @param odeModel a list that represents an ODE
#' @export
#' @return a character vector with the generated code, one
#'     vector-element is one line of code.
#' @examples
#' f <- uqsa_example("AKAR4")
#' m <- model_from_tsv(f)
#' C <- generateCode(as_ode(m))
#' cat(head(C,12),sep='\n')
generateCode <- function(odeModel){
	warning("This function will start a background yacas process via Ryacas.\nThere is currently no working way to reset/restart that process.\n It is therefore not advisable to generate the code for two different models in the same R session.\nThe definitions for the two models will be mixed up. ")
	modelName <- comment(odeModel) %otherwise% "model"
	# simplify a data.frame with two columns to a named character vector, assuming it's name/value pairs
	makeEnum <- \(varNames, enumName, lastEntry, prefix="_") {
		if (is.null(varNames)) return(c())
		else return(sprintf("enum %s { %s, %s };",enumName,paste(prefix,varNames,sep='', collapse=', '),lastEntry))
	}
	h <- c('stdlib','math','string','gsl/gsl_errno','gsl/gsl_odeiv2','gsl/gsl_math')
	C <- c(sprintf("#include <%s.h>",h),"")
	# C expressions:
	x <- odeModel$exp
	for (i in seq_along(x)){
		CMD <- sprintf("%s := %s",yacasMath(names(x)[i]),yacasMath(x[i]))
		#cat("CMD: ",CMD,"\n")
		Ryacas::yac_silent(yacasMath(CMD))
	}
	n <- pmax(5,50 - nchar(names(odeModel$par))*2)
	m <- pmax(5,50 - nchar(names(odeModel$var))*2)
	cc <- c(writeComment("\tconstants"),
		sprintf("\tdouble %s = %s;",names(odeModel$const),as.character(odeModel$const)))
	cp <- c(writeComment("\tparameter values"),
		sprintf("\tdouble %s = p_[_%s]; %*s /* [%3i] */",names(odeModel$par),names(odeModel$par),n," ",cOffset(odeModel$par)))
	cy <- c(
		writeComment("\tstate variables"),
		sprintf("\tdouble %s = y_[_%s]; %*s /* [%3i] */",names(odeModel$var),names(odeModel$var),m," ",cOffset(odeModel$var)))
	cx <- c(
		writeComment("\texpressions"),
		sprintf("\tdouble %s = %s;",names(odeModel$exp), replace_powers(odeModel$exp)))
	C <- c(C,
	writeComment("Enums will be used for indexing purposes."),
	makeEnum(names(odeModel$var),'stateVariable','numStateVar'),
	makeEnum(names(odeModel$par),'param','numParam'),
	makeEnum(names(odeModel$func),'func','numFunc'),
	makeEnum(colnames(odeModel$tf),'eventLabel','numEvents',prefix=NULL),"", # ens with a blank line
	writeComment(c("The error codes indicate how many values a function returns.",
		"Each function expects the output buffer to be allocated with at least that many values")))
	# Jacobian
	J <- replace_powers(t(yJacobian(odeModel$vf,names(odeModel$var))))
	attr(J,"stride") <- length(odeModel$var)
	# Parameter Jacobian
	Jp <- replace_powers(t(yJacobian(odeModel$vf,names(odeModel$par))))
	attr(Jp,"stride") <- length(odeModel$par)
	# Output-Function Jacobian
	fJ <- replace_powers(t(yJacobian(odeModel$func,names(odeModel$var))))
	attr(fJ,"stride") <- length(odeModel$var)
	# Output-Function Parameter Jacobian
	fJp <- replace_powers(t(yJacobian(odeModel$func,names(odeModel$par))))
	attr(fJp,"stride") <- length(odeModel$par)

	C <- c(C,"",
		writeComment("ODE vector field: y' = f(t,y;p)"),
		writeCFunction(modelName,"vf",
			retValue="f_",
			defs=c(cc,cp,cy,cx),
			values=replace_powers(odeModel$vf),
			init0=TRUE),"",
		writeComment("ODE Jacobian: df(t,y;p)/dy"),
		writeCFunction(modelName,"jac",
			retValue=c("jac_","dfdt_"),
			defs=c(cc,cp,cy,cx),
			values=J,
			init0=TRUE),"",
		writeComment("ODE parameter Jacobian: df(t,y;p)/dp"),
		writeCFunction(modelName,"jacp",
			retValue=c("jacp_","dfdt_"),
			defs=c(cc,cp,cy,cx),
			values=Jp,
			init0=TRUE),"",
		writeComment("Output Function (Observables)"),
		writeCFunction(modelName,"func",
			retValue="func_",
			defs=c(cc,cp,cy,cx),
			values=odeModel$func,
			init0=FALSE),"",
		writeComment("Output function Jacobian: dF(t,y;p)/dx"),
		writeCFunction(modelName,"funcJac",
			retValue=c("funcJac_"),
			defs=c(cc,cp,cy,cx),
			values=fJ,
			init0=TRUE),"",
		writeComment("Output function parameter Jacobian: dF(t,y;p)/dp"),
		writeCFunction(modelName,"funcJacp",
			retValue=c("funcJacp_"),
			defs=c(cc,cp,cy,cx),
			values=fJp,
			init0=TRUE),"",
		writeCFunction(modelName,"default",
			defArgs="double t",
			retValue="p_",
			defs=cc,
			values=odeModel$par,
			otherArgs=NULL),"",
		writeCFunction(modelName,"init",
			defArgs="double t",
			retValue="y_",
			defs=c(cc,cp),
			values=odeModel$var)
		)
	if ("tf" %in% names(odeModel) && !is.null(odeModel$tf)){
		C <- c(C,"",
			writeCFunction(modelName,"event",
				defArgs=c("double t", "double *y_"),
				defs=c(cc,cp,cy,cx),
				body=eventCode(odeModel),
				otherArgs=c("double *p_","int EventLabel","double dose")
			)
		)
	}
	# event function
	comment(C) <- comment(odeModel) %otherwise% modelName
	return(C)
}

writeRFunction <- function(prefix,fName,ret,value,arguments=c('t','state','parameters'),defs=NULL){
	z <- grepl("^0$",as.character(value))
	if (is.matrix(value)){
		n <- NROW(value)
		m <- NCOL(value)
		i <- matrix(seq(n),n,m,byrow=FALSE)
		j <- matrix(seq(m),n,m,byrow=TRUE)
		v <- c(
			sprintf("\t%s <- matrix(0.0,%i,%i)",ret,n,m),
			sprintf("\t%s[%i,%i] <- %s",ret,i[!z],j[!z],as.character(value)[!z])
		)
	} else {
		v <- c(
			sprintf("\t%s <- numeric(%i)",ret,length(value)),
			sprintf("\t%s[%i] <- %s",ret,seq_along(value)[!z],value[!z])
		)
	}
	f <- c(
		sprintf("%s_%s <- function(%s) {",prefix,fName,paste(arguments,collapse=', ')),
		defs,
		v,
		ifelse(is.null(names(value)),"## no names for return value",sprintf("\tnames(%s) <- c(%s)",ret,paste0('"',names(value),'"',collapse=', '))),
		switch(fName,vf=sprintf("\treturn(list(%s))",ret),sprintf("\treturn(%s)",ret)),
		"}"
	)
}

#' Write R code
#'
#' This function expects a list of character vectors, as returned by
#' `as_ode`. This list describes an ODE model
#' (initial values, default parameters, transformation events, output
#' functions).  This function uses this information, calculates
#' Jacobians via Ryacas and returns a character vector with R source
#' code for the deSolve package.
#'
#' The value can be written to a file:
#' `cat(generateRCode(odeModel),sep="\n",file=...)`. This file can be
#' sourced later.
#'
#' @param odeModel a list of named character vectors with math
#'     expressions (which work as R code)
#' @export
#' @return a character vector with the generated code, one
#'     vector-element is one line of code.
#' @examples
#' f <- uqsa_example("AKAR4")
#' m <- model_from_tsv(f)
#' generateRCode(as_ode(m))
generateRCode <- function(odeModel){
	modelName <- comment(odeModel) %otherwise% "model"
	# simplify a data.frame with two columns to a named character vector, assuming it's name/value pairs
	x <- odeModel$exp
	for (i in seq_along(x)){
		CMD <- sprintf("%s := %s",names(x)[i],x[i])
		Ryacas::yac_silent(yacasMath(CMD))
	}
	n <- pmax(5,50 - nchar(names(odeModel$par)))
	m <- pmax(5,50 - nchar(names(odeModel$var)))
	cc <- c(sprintf("##\tconstants"),
		sprintf("\t%s = %s;",names(odeModel$const),as.character(odeModel$const)))
	cp <- c(sprintf("##\tparameter values"),
		sprintf("\t%s = parameters[%i];",names(odeModel$par),seq_along(odeModel$par)))
	cy <- c(
		sprintf("##\tstate variables"),
		sprintf("\t%s = state[%i];",names(odeModel$var),seq_along(odeModel$var)))
	cx <- c(
		sprintf("##\texpressions"),
		sprintf("\t%s = %s;",names(odeModel$exp), odeModel$exp))
	mc <- c('vf','jac','jacp','default','init','func')
	# Jacobian
	J <- yJacobian(odeModel$vf,names(odeModel$var))
	# parameter Jacobian
	Jp <- yJacobian(odeModel$vf,names(odeModel$par))
	# code string:
	C <- c(
		sprintf("# ODE vector field: y' = f(t,y;p)"),
		writeRFunction(modelName,"vf",
			ret="f_",
			defs=c(cc,cp,cy,cx),
			value=odeModel$vf,
		),"",
		sprintf("# ODE Jacobian: df(t,y;p)/dy"),
		writeRFunction(modelName,"jac",
			ret="jac_",
			defs=c(cc,cp,cy,cx),
			value=J,
		),"",
		sprintf("# ODE parameter Jacobian: df(t,y;p)/dp"),
		writeRFunction(modelName,"jacp",
			ret="jacp_",
			defs=c(cc,cp,cy,cx),
			value=Jp,
		),"",
		sprintf("# Output Function (Observables)"),
		writeRFunction(modelName,"func",
			ret="func_",
			defs=c(cc,cp,cy,cx),
			value=odeModel$func,
		),"",
		writeRFunction(modelName,"default",
			ret="parameters",
			defs=cc,
			value=odeModel$par,
		arguments="t"
		),"",
		writeRFunction(modelName,"init",
			ret="state",
			defs=c(cc,cp),
			value=odeModel$var,
			arguments=c('t','parameters')
		),"",
		"# a variable that collects all functions into one list:",
		sprintf("model <- list(%s,name='%s')",paste0(mc,'=',modelName,"_",mc,collapse=', '),modelName)
	)
	# event function
	return(C)
}

#' Writes code to file and compiles
#'
#' This function accepts the code that was written by [generateCode],
#' possibly changed by the user. It writes the contents to a c file
#' named 'modelName_gvf.c'. This file is compiled to './modelName.so'.
#'
#' This entire function can be replaced with a call to `cat()` and
#' then compiling the written file in the system's shell (or via
#' [checkModel])
#' @param C character vector with the code of this model
#' @return the model's name with annotation about file names.
#' @export
#' @examples
#' \dontrun{
#'   cCode <- generateCode(m,o)             # a character vector
#'   modelName <- write_and_compile(cCode)  # commented name
#'   print(comment(modelName))
#' }
write_and_compile <- function(C){
	c.file <- sprintf("./%s_gvf.c",comment(C))
	cat(C,sep="\n",file=c.file)
	if (file.exists(c.file)){
		message(sprintf("'%s' was created.",c.file))
	} else {
		stop("writing c-file to current working directory failed.")
	}
	return(checkModel(comment(C),c.file))
}
