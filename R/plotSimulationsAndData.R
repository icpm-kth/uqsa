#' plot function for experiments
#'
#' This function uses plot.errors and the base plot functions like
#' [matplot].
#'
#' @param x simulation results (a list)
#' @param y experiment setup (a list)
#' @param ... forwarded to the more specific plot function [plot.errors]
#' @return plot object
#' @export
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- write_and_compile(as_ode(m))
#' ex <- experiments(m,o)
#' s <- simulator.c(ex,o)
#' p0 <- values(m$Parameter)
#' y <- s(p0)
#' plot(y,ex)
plot.simulation <- function(x, y, ...){
	simulations <- x
	if (missing(y)) experiments <- NULL
	else experiments <- y
	nf <- dim(simulations[[1]]$func)[1]
	for (i in seq_along(simulations)){
		time <- experiments[[i]]$outputTimes
		n <- dim(simulations[[i]]$func)
		oNames <- rownames(simulations[[i]]$func)
		for (j in seq(n[1])){
			a <- ceiling(255*exp(-0.01*n[3]))
			if (!is.null(experiments) && !is.null(experiments[[i]]$data) && oNames[j] %in% rownames(experiments[[i]]$data)){
				k <- match(oNames[j],rownames(experiments[[i]]$data))
				d <- experiments[[i]]$data[k,]
				p <- plot(
					errors::as.errors(time),
					d,
					main=names(experiments)[i],
					xlab="time",
					ylab=oNames[j],
					lty=1,
					...
				)
				matplot(
					time,
					simulations[[i]]$func[j,,],
					col=rgb(0,0,255,a,maxColorValue=255),
					lwd=2,
					add=TRUE,
					type="l",
					lty=1
				)
			} else {
				matplot(
					time,
					simulations[[i]]$func[j,,],
					col=rgb(0,0,255,a,maxColorValue=255),
					lwd=2,
					type="l",
					lty=1
				)
			}
		}
	}
	return(p)
}

