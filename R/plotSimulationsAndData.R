#' Plot function for experiments
#'
#' This function uses plot.errors and the base plot functions like
#' [matplot].
#' @keywords internal
#' @param x experiment setup (a list)
#' @param y simulation results (a list)
#' @param func.names restrict the shown output functions to these
#'     named outputs
#' @param verbose prints additional information about name matching
#'     between experiments and simulations
#' @param ... forwarded to the more specific plot function
#'     [plot.errors]
#' @return plot object
#' @export
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- write_and_compile(as_ode(m))
#' ex <- experiments(m,o)
#' s <- simulator.c(ex,o)
#' p0 <- values(m$Parameter)
#' y <- s(p0)
#' plot(ex,y)
plot.experiments <- function(x, y, func.names=NULL, verbose=getOption("uqsa.verbose",FALSE),...){
	if (missing(y)) {
		simulations <- NULL
	} else {
		stopifnot(length(x)==length(y))
		simulations <- y
		nf <- dim(simulations[[1]]$func)[1]
	}
	if (missing(func.names)){
		out <- seq(NROW(x[[1]]$data))
	} else if (is.character(func.names)) {
		out <- which(rownames(x[[1]]$data) %in% func.names)
	} else if (is.numeric(func.names)){
		out <- func.names
	}
	p <- list()
	for (i in seq_along(x)){
		time <- x[[i]]$outputTimes
		if (!is.null(simulations)){
			n <- dim(simulations[[i]]$func)
		} else {
			n <- c(dim(x[[i]]$data),1)
		}
		errors(ex[[i]]$data)[!is.finite(errors(ex[[i]]$data))] <- NA
		ex[[i]]$data[!is.finite(errors(ex[[i]]$data))] <- NA
		oNames <- rownames(x[[i]]$data)
		for (j in out){
			a <- ceiling(255*exp(-0.01*n[3]))
			## only plot if the output was actually measured
			if (verbose) cli::cli_alert_info("{oNames[j]} %in% {paste(colnames(x[[i]]$measurements),collapse=' ')}")
			if (oNames[j] %in% colnames(x[[i]]$measurements)){
				d <- x[[i]]$data[j,]
				p <- c(p,
					plot(
						errors::as.errors(time),
						d,
						main=names(x)[i],
						xlab="time",
						ylab=oNames[j],
						lty=1,
						...
					)
				)
				if (!is.null(simulations)){
					matplot(
						time,
						simulations[[i]]$func[j,,],
						col=rgb(0,0,255,a,maxColorValue=255),
						lwd=2,
						add=TRUE,
						type="l",
						lty=1
					)
				}
			}
		}
	}
	return(invisible(p))
}
