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
#' @param ode a list of data.frames with the text representation of the ordinary differential equation
#' @export
#' @return a character vector with the generated code, one element is one line.
generateCode <- function(ode){
	modelName <- comment(ode) %otherwise% "model"

	# simplify a data.frame with two columns to a named character vector, assuming it's name/value pairs
	makeEnum <- \(var,enumName, lastEntry) {
		sprintf("enum %s { %s, %s };",enumName,paste('_',var,sep='', collapse=', '),lastEntry)
	}
	writeFunc <- \(fName,retValue,body,len=length(body),Rest=NULL){
		if (length(retValue)>0) ret <- paste("double *",retValue,sep="", collapse=", ",recycle0=TRUE)
		else ret <- NULL
		arguments <- paste(c("double t","const double y_[]",ret,"void *par",Rest),sep=", ",collapse=", ")
		n <- pmax(5,50 - nchar(ode$par[[1]])*2)
		m <- pmax(5,50 - nchar(ode$var[[1]])*2)
		return(c(
			sprintf("int %s_%s(%s)",modelName,fName,arguments),
			sprintf("\tdouble *p_=par;"),
			sprintf("\tif (!y_ || !%s_) return %i;",retValue[1],len),
			sprintf("\tdouble %s = %g;",ode$const[[1]],ode$const[[2]]),
			sprintf("\tdouble %s = p_[_%s]; %*s /* [%3i] */",ode$par[[1]],ode$par[[1]],n," ",cOffset(ode$par)),
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
		CMD <- sprintf("%s := %s",names(x)[i],x[i])
		Ryacas::yac_silent(yacasMath(CMD))
	}
	# Jacobian
	C <- c(C,"",writeComment("ODE Jacobian: df(t,y;p)/dy"))
	J <- replace_powers(t(yJacobian(ode$vf[[2]],ode$var[[1]])))
	z <- grepl("^0$",J)
	C <- c(C,writeFunc("jac",c("jac_","dfdt_"),sprintf("\tjac_[%i] = %s;",cOffset(J)[!z],J[!z]),length(J)))
	# parameter Jacobian
	C <- c(C,"",writeComment("ODE parameter Jacobian: df(t,y;p)/dp"))
	Jp <- replace_powers(t(yJacobian(ode$vf[[2]],ode$par[[1]])))
	z <- grepl("^0$",Jp)
	C <- c(C,writeFunc("jacp",c("jacp_","dfdt_"),sprintf("\tjacp_[%i] = %s;",cOffset(Jp)[!z],Jp[!z]),length(Jp)))
	# output function
	C <- c(C,"",writeComment("Output Function (Observables)"))
	C <- c(C,writeFunc("func","func_",sprintf("\tfunc_[_%s] = %s;",ode$func[[1]],ode$func[[2]])))
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

#' Load an ODE model from a file
#'
#' The input can be either a list of files, or one compressed archive
#' file like "model.tar.gz". This function determines the name of the
#' model from the name of the archive, or from the parent directory of
#' uncompressed tsv files. This auto-determined name is attached as a
#' comment of the return value. To override this choice, replace the
#' comment.
#'
#' @param fileList list of file paths
#' @export
#' @return a list of data.frames describing the model
loadODE <- function(fileList){
	modelName <- basename(dirname(fileList[1]))
	if (length(fileList)==1 && endsWith(basename(fileList),".tar.gz")) {
		modelName <- sub("[.].*$","",basename(fileList[1]))
		tartf <- untar(fileList,list=TRUE)
		TMP="/dev/shm/ode_gen"
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
