#' yacasMath converts math to Ryacas compatible math
#'
#' Given a string like "exp(2*x)" this function returns a string that
#' yacas can process: "Exp(2*x)"
#'
#' @param v a chcarcter vector with math expressions
#' @param reverse do the reverse operation
#' @return a character vector with yacas math
#' @export
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
#' @param f a character vector of length n
#' @param x a character vector of length m
#' @return a character matrix (n×m) with derivatives df[i]/dx[j]
#' @export
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

#' Write C code
#'
#' This function generates C code, expecting the information about the
#' ODE to come from the files that SBtabVFGEN::sbtab_to_vfgen()
#' writes. It exists for compatibility with the scripts in the
#' RPN-derivative package.
#'
#' sbtab_to_vfgen() also returns a variable. This variable (a list)
#' has the same information as the files and can be used directly (by a different function).
#'
#' @param ode a list of data.frames with the text representation of the ordinary differential equation
#'
#' @return a character vector with the generated code, one element is one line.
generateCodeFromFile <- function(fileList){
	ode <- loadODE(fileList)
	modelName <- comment(ode) %otherwise% "model"

	# simplify a data.frame with two columns to a named character vector, assuming it's name/value pairs
	makeEnum <- \(var,enumName, lastEntry) {
		if (is.null(var)) return(c())
		else return(sprintf("enum %s { %s, %s };",enumName,paste('_',var,sep='', collapse=', '),lastEntry))
	}
	initVal <- \(y) {
		arguments <- "double t, const double y_[], void *par"
		n <- pmax(5,50 - nchar(y) - nchar(names(y)))
		m <- pmax(5,50 - 2*nchar(ode$par[[1]]))
		return(c(
			sprintf("int %s_init(%s){",modelName,arguments),
			sprintf("\tif (!y_) return %i;",length(y)),
			sprintf("\tdouble *p_=par;"),
			sprintf("\tdouble %s = %g;",ode$const[[1]],ode$const[[2]]),
			writeComment("\tparameter values"),
			sprintf("\tdouble %s = p_[_%s]; %*s /* [%3i] */",ode$par[[1]],ode$par[[1]],m," ",cOffset(ode$par)),
			sprintf("\ty[_%s] = %s; %*s /*[%2i]*/",names(y),y,n," ",cOffset(y)),
			sprintf("\treturn GSL_SUCCESS;\n}")
		))
	}
	defaultPar <- \(p) {
		arguments <- "double t, void *par"
		n <- pmax(5,50 - nchar(p) - nchar(names(p)))
		return(c(
			sprintf("int %s_default(%s){",modelName,arguments),
			sprintf("\tif (!par) return %i;",length(p)),
			sprintf("\tdouble *p_=par;"),
			sprintf("\tdouble %s = %g;",ode$const[[1]],ode$const[[2]]),
			sprintf("\tp[_%s] = %s; %*s /*[%2i]*/",names(p),p,n," ",cOffset(p)),
			sprintf("\treturn GSL_SUCCESS;\n}")
		))
	}

	writeFunc <- \(fName,retValue,body,len=length(body),Rest=NULL){
		if (length(retValue)>0) ret <- paste("double *",retValue,sep="", collapse=", ",recycle0=TRUE)
		else ret <- NULL
		arguments <- paste(c("double t","const double y_[]",ret,"void *par",Rest),sep=", ",collapse=", ")
		n <- pmax(5,50 - nchar(ode$par[[1]])*2)
		m <- pmax(5,50 - nchar(ode$var[[1]])*2)
		return(c(
			sprintf("int %s_%s(%s){",modelName,fName,arguments),
			sprintf("\tdouble *p_=par;"),
			sprintf("\tif (!y_ || !%s_) return %i;",retValue[1],len),
			writeComment("\tconstants"),
			sprintf("\tdouble %s = %g;",ode$const[[1]],ode$const[[2]]),
			writeComment("\tparameter values"),
			sprintf("\tdouble %s = p_[_%s]; %*s /* [%3i] */",ode$par[[1]],ode$par[[1]],n," ",cOffset(ode$par)),
			writeComment("\tstate variables"),
			sprintf("\tdouble %s = y_[_%s]; %*s /* [%3i] */",ode$var[[1]],ode$var[[1]],m," ",cOffset(ode$var)),
			sprintf("\tdouble %s = %s;",ode$exp[[1]],replace_powers(ode$exp[[2]])),
			sprintf("\tmemset(%s,0,sizeof(double)*%i); ",retValue[1],len),
			body,
			sprintf("\treturn GSL_SUCCESS;\n}")
		))
	}

	h <- c('stdlib','math','string','gsl/gsl_errno','gsl/gsl_odeiv2','gsl/gsl_math.h')
	C <- c(sprintf("#include <%s.h>",h),"",
	writeComment("Enums will be used for indexing purposes."),
	makeEnum(ode$vf[[1]],'stateVariable','numStateVar'),
	makeEnum(ode$par[[1]],'param','numParam'),
	makeEnum(ode$func[[1]],'func','numFunc'),
	makeEnum(ode$events[[1]],'eventLabel','numEvents'),"") # ens with a blank line

	C <- c(C,writeComment(
		c("The error codes indicate how many values a function returns.",
			"Each function expects the output buffer to be allocated with at least that many values")))
	# vf function
	C <- c(C,writeComment("ODE vector field: y'f(t,y;p)"))
	C <- c(C,writeFunc("vf","f_",sprintf("\tf_[_%s] = %s;",ode$vf[[1]],replace_powers(ode$vf[[2]]))))
	x <- ch(ode$exp)
	for (i in seq_along(x)){
		CMD <- sprintf("%s := %s",yacasMath(names(x)[i]),x[i])
		Ryacas::yac_silent(yacasMath(CMD))
	}
	# Jacobian
	C <- c(C,"",writeComment("ODE Jacobian: df(t,y;p)/dy"))
	J <- replace_powers(t(yJacobian(ode$vf[[2]],ode$var[[1]])))
	z <- grepl("^0$",J)
	i <- cOffset(J) %/% NROW(ode$var)
	j <- cOffset(J) %% NROW(ode$var)
	C <- c(C,writeFunc("jac",c("jac_","dfdt_"),sprintf("\t/*[%2i,%2i]*/  jac_[%i] = %s; ",i[!z],j[!z],cOffset(J)[!z],J[!z]),length(J)))
	# parameter Jacobian
	C <- c(C,"",writeComment("ODE parameter Jacobian: df(t,y;p)/dp"))
	Jp <- replace_powers(t(yJacobian(ode$vf[[2]],ode$par[[1]])))
	z <- grepl("^0$",Jp)
	i <- cOffset(Jp) %/% NROW(ode$par)
	j <- cOffset(Jp) %% NROW(ode$par)
	C <- c(C,writeFunc("jacp",c("jacp_","dfdt_"),sprintf("\t/*[%2i,%2i]*/  jacp_[%i] = %s; ",i[!z],j[!z],cOffset(Jp)[!z],Jp[!z]),length(Jp)))
	# output function
	C <- c(C,"",writeComment("Output Function (Observables)"))
	C <- c(C,writeFunc("func","func_",sprintf("\tfunc_[_%s] = %s;",ode$func[[1]],ode$func[[2]])))
	# initial values and default parameters
	y <- ch(ode$var)
	p <- ch(ode$par)
	C <- c(C,initVal(y),defaultPar(p))
	# event function
	if ("tf" %in% names(ode)){
		C <- c(C,"",writeComment(c("Scheduled Events","EventLabel specifies which transformation to apply.","The scalar dose variable can be used in the transformation (on the right)")))
		v <- ode$tf[[2]]=="var"
		Effect <- \(ev){
			l <- ode$tf[[1]]==ev
			y <- ch(ode$tf[l & v,c(3,4)])
			p <- ch(ode$tf[l & !v,c(3,4)])
			return(
				c(sprintf("\tcase %s:",ev),
					sprintf("\t\ty_[_%s] = %s",names(y),y),
					sprintf("\t\tp_[_%s] = %s",names(p),p),
					sprintf("\tbreak;")))
		}
		for (ev in unique(ode$tf[[1]])){
		}
		evsw <- c(
		 sprintf("\tswitch(eventLabel){"),
		 unlist(lapply(unique(ode$tf[[1]]),Effect)),
		 "\t}"
		)
		C <- c(C,writeFunc("event",NULL,evsw,Rest="int eventLabel, double dose"))
	}
	return(C)
}

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
writeCFunction <- function(prefix, fName, defArgs=c("double t","const double y_[]"), retValue=NULL, defs=NULL, values=NULL, body=NULL, otherArgs="void *par", init0 = TRUE){
	if (is.null(values)) return(character(0))

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
		errValue <- switch(fName,default=sprintf("\tif (!p_) return numParam;"),event=sprintf("\tif (!y_ || EventLabel<0) return numEvents;"),NULL)
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

eventCode <-function(odeModel){
	if ("tf" %in% names(odeModel)){
		C <- c(C,"",writeComment(c("Scheduled Events","EventLabel specifies which transformation to apply.","The scalar dose variable can be used in the transformation (on the right)")))
		affectedVector <- sapply(attr(odeModel$tf,"type"),
			\(x) switch(x,var="y",par="p")
		)
		effect <- lapply(seq(NCOL(odeModel$tf)),\(i){
			ev <- odeModel$tf[,i]
			trivial <- names(ev) == ev
			return(c(
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
#' SBtabVFGEN::sbtab_to_vfgen. This list describes an ODE model
#' (initial values, default parameters, transformation events, output
#' functions).  This function uses this information, calculates
#' Jacobians via Ryacas and returns a character vector with C source
#' code for the solvers in the GNU Scientific Library (GSL).
#'
#' The value can be written to a file:
#' `cat(generateCode(odeModel),sep="\n",file=...)`. This file can be
#' compiled into a shared library and used by the solvers in the
#' icpm-kth/rgsl package.
#'
#' @param odeModel a list, as returned from SBtabVFGEN::sbtab_to_vfgen()
#' @export
#' @return a character vector with the generated code, one vector-element is one line of code.
generateCode <- function(odeModel){
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
	# parameter Jacobian
	Jp <- replace_powers(t(yJacobian(odeModel$vf,names(odeModel$par))))
	attr(Jp,"stride") <- length(odeModel$par)

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
			values=odeModel$var),"",
		writeCFunction(modelName,"event",
			defArgs=c("double t", "double *y_"),
			defs=c(cc,cp,cy,cx),
			body=eventCode(odeModel),
			otherArgs=c("double *p_","int EventLabel","double dose"))
	)
	# event function
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
#' SBtabVFGEN::sbtab_to_vfgen. This list describes an ODE model
#' (initial values, default parameters, transformation events, output
#' functions).  This function uses this information, calculates
#' Jacobians via Ryacas and returns a character vector with R source
#' code for the deSolve package.
#'
#' The value can be written to a file:
#' `cat(generateRCode(odeModel),sep="\n",file=...)`. This file can be
#' sourced later.
#'
#' @param odeModel a list, as returned from SBtabVFGEN::sbtab_to_vfgen()
#' @export
#' @return a character vector with the generated code, one vector-element is one line of code.
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


#' Load an ODE model from a file
#'
#' The input can be either a list of files, or one compressed archive
#' file like "model.tar.gz". This function determines the name of the
#' model from the name of the archive, or from the parent directory of
#' uncompressed tsv files. This auto-determined name is attached as a
#' comment of the return value. To override this choice, replace the
#' comment. This function exists for compatibility with the scripts in
#' the RPN-derivative repository and loads the zip and tar.gz files
#' that SBtabVFGEN::sbtab_to_vfgen writes to disk.
#'
#' @param fileList list of file paths
#' @export
#' @return a list of data.frames describing the model
loadODE <- function(fileList){
	modelName <- basename(dirname(fileList[1]))
	if (length(fileList)==1 && endsWith(basename(fileList),".tar.gz")) {
		modelName <- sub("[.].*$","",basename(fileList[1]))
		tartf <- untar(fileList,list=TRUE)
		if (file.exists("/dev/shm")) {
			TMP <- "/dev/shm/ode_gen"
		} else {
			TMP <- "/tmp/ode_gen"
		}
		untar(fileList,exdir=TMP)
		fileList <- paste0(TMP,'/',tartf)
	}
	rf <- \(f){
		read.delim(file=f,header=FALSE)
	}
	ode <- lapply(fileList,rf)
	Name <- tolower(sub("[.].*$","",basename(fileList)))
	rename <- \(x) {
		switch(x,statevariables="var",parameters="par",ode="vf",expressions="exp",outputfunctions="func",transformations="tf")
	}
	names(ode) <- lapply(Name,rename)
	comment(ode) <- modelName
	return(ode)
}
