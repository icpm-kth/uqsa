#' Read Concise Error Notation
#'
#' Convert a vector of strings of the form: c("1.2(3)E-4","1.2(3)E-2") to a matrix with
#' two rows: values, and uncertainties.
#'
#' Cincise error notation means that a floating point number is
#' followed by an integer in parentheses which indicates the
#' uncertainty of the last digits of the value:
#' 1.2345(12) = 1.2345 ± 0.0012.
#'
#' If the errors package is installed, then it will be used to
#' represent the return value.
#'
#' @export
#' @param v a character vector of numbers in concise error notation
#' @param use.errors if TRUE, the errors package will be used to
#'     retrun an object of type "errors" (from that package).
#' @return a matrix of values and uncertainties (2 rows), or errors object
#' @useDynLib uqsa, concise
parse_concise <- function(v,use.errors=requireNamespace("errors")){
	w<-.Call(concise,as.character(v))
	if (use.errors){
		return(errors::set_errors(w[1,],w[2,]))
	} else {
		return(w)
	}
}
