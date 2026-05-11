
#' plot function for experiments
#'
#' This function does not use ggplot2, only base plot functions like
#' lines() and arrows().
#'
#' @param simulations simulation results (a list)
#' @param experiments experiment setup (a list)
#' @param by skip this many sampled lines in between plotted lines
#' @param ylimit y-axis limits
#' @return plot object
#' @export
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- write_and_compile(as_ode(m))
#' ex <- experiments(m,o)
#' s <- simulator.c(ex,o)
#' p0 <- values(m$Parameter)
#' y <- s(p0)
#' plot_time_series_base(y,ex)
plot_time_series_base <- function(simulations, experiments, by=1, ylimit=NULL){
	nf <- dim(simulations[[1]]$func)[1]
	par(mfrow=c(nf,length(experiments)))
	for (i in seq_along(experiments)){
		time <- experiments[[i]]$outputTimes
		n <- dim(simulations[[i]]$func)
		oNames <- names(experiments[[i]]$measurements)
		for (j in seq(n[1])){
			d <- experiments[[i]]$data[j,] %otherwise% experiments[[i]]$measurements[[j]]
			if (is.null(ylimit) || !all(is.finite(ylimit[[j]]))){
				yl <- c(min(d - 1e-1,na.rm=TRUE),max(d + 1e-1,na.rm=TRUE))
			} else {
				yl <- ylimit[[j]]
			}
			p <- plot(
				errors::as.errors(time),
				d,
				ylim=yl,
				main=names(experiments)[i],
				xlab="time",
				ylab=oNames[j]
			)
			a <- ceiling(255*exp(-0.01*n[3]))
			matplot(
				time,
				simulations[[i]]$func[j,,],
				col=rgb(0,0,255,a,maxColorValue=255),
				lwd=2,
				add=TRUE,
				type="l"
			)
		}
	}
	return(p)
}


#' Plot time series simulations with experimental data
#'
#' This function plots simulations of time series experiments and plots them
#' against experimental data. The input in the provided experiments must differ
#' only in one vector component.
#'
#' @noRd
#' @param simulations list of simualtions as output from the simulator
#' @param experiments list of experiments
#' @param show.plot boolean variable. Set `show.plot=TRUE` to display plots
#'      when running the funcion, FALSE otherwise
#' @return list of plots with simulations and experimental data
#' @examples
#' m <- model_from_tsv(uqsa_example("AKAR4"))
#' o <- write_and_compile(as_ode(m))
#' ex <- experiments(m,o)
#' s <- simulator.c(ex,o,log10ParMap)
#' p0 <- values(m$Parameter)
#' y <- s(p0)
#' plotTimeSeries(y,ex)
plotTimeSeries <- function(simulations, experiments, show.plot = TRUE){
  num.experiments <- length(simulations)
  num.simulations <- dim(simulations[[1]]$func)[3]
  stopifnot(num.experiments == length(experiments))
  num.out.funcs <- dim(experiments[[1]]$measurements)[2]
  p <- list()
  for(i in 1:num.experiments){
    experiment <- experiments[[i]]
    for(output.idx in 1:num.out.funcs){
      df.experiments <- data.frame(t=experiment[["outputTimes"]], y=experiment[["measurements"]][[output.idx]])
      y <- c(simulations[[i]]$func[output.idx,,])
      df.simulations <- data.frame(t=rep(experiment[["outputTimes"]],num.simulations), y=y, sim=rep(1:num.simulations,each=length(experiment[["outputTimes"]])))
      p[[i]] <-
        ggplot2::ggplot(df.simulations,ggplot2::aes(x=t, y=y, group=sim)) +
        ggplot2::geom_line(color="blue", alpha = 0.1, size=1.5) +
        ggplot2::geom_point(data=df.experiments, ggplot2::aes(x=t, y=y), inherit.aes=FALSE) +
        ggplot2::ggtitle(paste0("Experiment: ", names(experiments[i]), "\nOutput: ", names(experiment[["measurements"]])[output.idx]))
    }
  }
  if (show.plot) show(do.call(gridExtra::grid.arrange,p))
  return(p)
}

#' Get the values of the input for a series of dose response experiments
#'
#' This function finds the vector of inputs that  varies among a series of experiments
#' that are part of a dose response experiment
#'
#' @noRd
#' @param experiments list of experimental data from the same dose response experiment
#' @return vector of inputs (i.e. dose) that varies among the experiments provided to the function.
#'         The name of the input is saved in "comment(dose)".
getDose <- function(experiments){
  all.inputs <- sapply(experiments, function(e) e$input)
  unique.inputs <- apply(all.inputs,1,unique)
  dose.idx <- which(sapply(unique.inputs, function(i) length(i)>1))
  stopifnot(length(dose.idx)==1)
  dose <- all.inputs[dose.idx,]
  comment(dose) <- names(dose.idx)
  return(dose)
}

#' Plot dose response simulations with experimental data
#'
#' This function plots simulations of one dose response experiment and plots them
#' against experimental data.
#'
#' @noRd
#' @param simulations list of simualtions as output from the simulator
#' @param experiments list of experimental data from the same dose response experiment
#' @param dose vector of dose values to plot on the x axis
#' @param show.plot boolean variable. Set `show.plot=TRUE` to display plots
#'      when running the funcion, `FALSE` otherwise
#' @return plot with simulations and experimental data
ggplotDoseResponse <- function(simulations, experiments, dose, show.plot = TRUE){
  num.experiments <- length(simulations)
  stopifnot(num.experiments == length(experiments))
  num.simulations <- dim(simulations[[1]]$func)[3]
  num.out.funcs <- dim(experiments[[1]]$measurements)[2]
  p <- list()
  for(j in 1:num.out.funcs){
    y.exp <- sapply(experiments, function(e) e$measurements[[j]])
    E.data.frame.exp <- data.frame(x=dose,y=y.exp)
    y.sim <- sapply(simulations, function(o) o$func[j,1,])
    E.data.frame.sim <- data.frame(x=(rep(dose,each=num.simulations)),y=c(y.sim))
    p[[j]] <- ggplot2::ggplot(E.data.frame.exp, ggplot2::aes(x=x,y=y)) + geom_point(color="red")
    p[[j]] <-  p[[j]] +
      ggplot2::geom_point(data=E.data.frame.sim, ggplot2::aes(x=x,y=y), inherit.aes = FALSE, color="blue", alpha = 0.4)
    p[[j]] <-  p[[j]] +
      ggplot2::geom_point(data=E.data.frame.exp, ggplot2::aes(x=x,y=y)) +
      ggplot2::geom_point(color="red") +
      ggplot2::labs(x = comment(dose), y = names(experiments[[1]][["measurements"]])[j])
  }
  return(p)
}

