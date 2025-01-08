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
	J <- yac_str(y_fn(sprintf("{%s},{%s}",F,X),"JacobianMatrix"))
	J <- gsub("[{}]","",J)
	J <- matrix(unlist(strsplit(J,",")),nrow=length(f),ncol=length(x),byrow=TRUE)
	return(yacasMath(J,rev=TRUE))
}

#' Make math C compatible
#'
#' Given a character array with math expressions, return a similar
#' vector, with c-compatible math. The most important change is
#' related to powers.
#'
#' @param v a character vector with math expressions
#' @param gsl use power functions from the GNU scientific library
#' @return a vector with c-compatible math
#' @export
cMath <- function(v,gsl=TRUE){
	# various integer powers
	v <- gsub("\\(([^()]+)\\)\\^\\(([[:digit:]]+[.][eE]?[-+]?[[:digit:]]+)\\)","pow(\\1,\\2)",v)
	v <- gsub("\\<([[:alnum:]]+)\\>\\^\\(([[:digit:]]+[.][eE]?[-+]?[[:digit:]]+)\\)","pow(\\1,\\2)",v)
	v <- gsub("\\<([[:alnum:]]+)\\>\\^([[:digit:]]+[.][eE]?[-+]?[[:digit:]]+)\\>","pow(\\1,\\2)",v)
	if (gsl){
		v <- gsub("\\(([^()]+)\\)\\^([2-9])","gsl_pow_\\2(\\1)",v)
		v <- gsub("\\b(\\w+)\\b\\^([2-9])","gsl_pow_\\2(\\1)",v)
		v <- gsub("\\b(\\d+)\\b\\^([2-9])","gsl_pow_\\2(\\1)",v)
		v <- gsub("\\(([^()]+)\\)\\^([0-9]+)","gsl_pow_int(\\1,\\2)",v)
		v <- gsub("\\b(\\w+)\\b\\^([0-9]+)","gsl_pow_int(\\1,\\2)",v)
		v <- gsub("\\b(\\d+)\\b\\^([0-9]+)","gsl_pow_int(\\1,\\2)",v)
	}
	v <- gsub("\\(([^()]+)\\)\\^\\(([^()]+)\\)","pow(\\1,\\2)",v)
	v <- gsub("\\b(\\w+)\\b\\^\\(([^()]+)\\)","pow(\\1,\\2)",v)
	v <- gsub("\\b(\\d+)\\b\\^\\(([^()]+)\\)","pow(\\1,\\2)",v)
	v <- gsub("\\(([[:alnum:]]+)\\)\\^\b([[:alnum:]]+)\b","pow(\\1,\\2)",v)
	return(v)
}

#' replace_powers
#'
#' @export
#' @param v a character vector
#' @useDynLib uqsa, replace_pow=replace_pow
replace_powers <- function(v){
	w<-.Call(replace_pow,v)
	return(w)
}
