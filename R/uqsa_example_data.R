#' Load an example model for this package
#'
#' This function finds the path to an example model, given by name.
#' In the SBtab format, model and data travel together (in different
#' tables, but the same documents).
#'
#' @param modelName name of model, e.g.: "AKAR4", "AKAP79", "CaMKII";
#'     if empty, this function lists all available examples.
#' @return The location of the example's directory (where the model is
#'     stored).
#' @export
uqsa_example<-function(modelName=NULL,full.names=TRUE,pat='[.]tsv$') {
	if (is.null(modelName)) {
		return(dir(system.file("extdata", package = "uqsa")))
	} else if (!is.null(pat) && is.character(pat)) {
		return(dir(system.file("extdata", modelName, package = "uqsa", mustWork = TRUE),pattern=pat,full.names=full.names))
	} else {
		return(system.file("extdata", modelName, package = "uqsa", mustWork = TRUE))
	}
}