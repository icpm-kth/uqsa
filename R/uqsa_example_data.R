#' Load an example model for this package
#'
#' This function finds the path to an example model, given by name.
#' In the SBtab format, model and data travel together (in different
#' tables, but the same documents).
#'
#' By default this function returns the names of the tsv files
#' belonging to the named model. If no modelName is provided it
#' returns possible names (contents of the top-level example
#' directory).
#'
#' @param modelName name of model, e.g.: "AKAR4", "AKAP79", "CaMKII";
#'     if empty, this function lists all available examples.
#' @param full.names return full paths to files - defaults to TRUE
#' @param pat pattern to find specific files; if NULL, this function
#'     returns the directory of the example
#' @param f alternatively, pat can be set to files ending in f [this
#'     value]
#' @return The location of the examples in the current environment if
#'     called with no arguments, the paths to the model files if a
#'     modelName was provided or the full path to the example if the
#'     file pattern _pat_ is unset
#' @examples
#' uqsa_example()
#' uqsa_example("AKAR4")
#' uqsa_example("AKAR4",pat=NULL)
#' @export
uqsa_example<-function(modelName=NULL,full.names=TRUE,pat='[.]tsv$',f=NULL) {
	if (!is.null(f)){
		pat <- sprintf("[.]%s$",f)
	}
	if (is.null(modelName)) {
		return(dir(system.file("extdata", package = "uqsa")))
	} else if (!is.null(pat) && is.character(pat)) {
		return(dir(system.file("extdata", modelName, package = "uqsa", mustWork = TRUE),pattern=pat,full.names=full.names))
	} else {
		return(system.file("extdata", modelName, package = "uqsa", mustWork = TRUE))
	}
}
