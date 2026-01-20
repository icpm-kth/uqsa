#' Plot time series simulations with experimental data
#'
#' This function plots simulations of time series experiments and plots them
#' against experimental data. The input in the provided experiments must differ
#' only in one vector component.
#'
#'
#' @export
#' @param simulations list of simualtions as output from the simulator
#' @param experiments list of experiments
#' @param show.plot boolean variable. Set show.plot=TRUE to display plots
#'      when running the funcion, FALSE otherwise
#' @return list of plots with simulations and experimental data
ggplotTimeSeries <- function(simulations, experiments, nrow=NULL, ncol=NULL, plot.state=FALSE){
	num.experiments <- length(simulations)

	stopifnot(length(simulations) == length(experiments))
	num.out.funcs <- NCOL(experiments[[1]]$outputValues)
	p <- list()
	for(i in seq(length(experiments))){
		N <- dim(simulations[[i]]$func)[3]
		oNames <- names(experiments[[i]][["outputValues"]])
		if ("measurementTimes" %in% names(experiments[[i]])){
			t_ <- experiments[[i]][["measurementTimes"]]
		} else {
			t_ <- experiments[[i]][["outputTimes"]]
		}
		for(j in seq(num.out.funcs)){
			df.experiments <- data.frame(t=t_, y=experiments[[i]][["outputValues"]][[j]])
			y <- as.numeric(simulations[[i]]$func[j,,])
			df.simulations <- data.frame(t=rep(experiments[[i]][["outputTimes"]],N), y=y, sim=rep(seq(N),each=length(experiments[[i]][["outputTimes"]])))
			p[[i]] <-
				ggplot2::ggplot(df.simulations,aes(x=t, y=y, group=sim))+
				ggplot2::geom_line(color="blue", alpha = 0.1, size=1.5)+
				ggplot2::geom_point(data=df.experiments, aes(x=t, y=y), inherit.aes=FALSE)+
				ggplot2::ggtitle(paste0("Experiment: ", names(experiments[i])))+labs(y=oNames[j])
		}
	}
	if (is.null(nrow)){
		nrow <- num.of.funcs
	}
	if (is.null(ncol)){
		ncol <- num.experiments
	}
	gridExtra::marrangeGrob(nrow=nrow,ncol=ncol,grobs=p)
}

#' Plot time series simulation with state variables
#'
#' This function plots simulations of time series experiments and plots them
#' against experimental data. The input in the provided experiments must differ
#' only in one vector component.
#'
#'
#' @export
#' @param simulations list of simualtions as output from the simulator
#' @param experiments list of experiments
#' @param var.names override the rownames of the simulation results
#' @param type 'boxes' or 'lines'
#' @param plot.states TRUE (or FALSE) - whether to plot the state
#'     variables or only the functions
#' @param ttf time transformation function - the plot will be against
#'     ttf(t), where `t` is a vector of the experiment's output times,
#'     ttf can adjust the time vector if it is very uneven or requires
#'     other modification only when plotting, e.g. `seq_along`.
#' @param xl x-axis label (time usually)
#' @param yl.func y-axis-limits of function plots, can be a list of
#'     ggplot2::ylim() objects, with NULL elements for automatic mode
#'     (the neutral element), NA elements will trigger tight bounds
#'     based on the maximum likelihood estimate and data. a simple
#'     numeric vector will be interpreted as quantiles for the
#'     quantiles function, the first and last quantile of the
#'     simulations will be used as ylim()
#' @param yl.state y-axis-limits for state variable plots, with
#'     similar rules as for yl.func
#' @return list of plots with simulations and experimental data
ggplotTimeSeriesStates <- function(simulations, experiments, var.names=NULL, type="boxes", plot.states=TRUE, ttf=identity, xl="t", yl.func=NULL, yl.state=NULL, MLE=1){
	num.experiments <- length(experiments)
	stopifnot(num.experiments == length(simulations))
	num.of.funcs <- NCOL(experiments[[1]]$outputValues)
	num.of.vars <- NROW(simulations[[1]]$state)
	if (missing(ttf)){
		t_txt <- xl
	} else {
		t_txt <- sprintf("%s(%s)",as.character(quote(ttf)),xl)
	}
	p <- list()
	M <- num.of.funcs + plot.states*num.of.vars
	T1 <- ggplot2::theme(plot.title=element_text(size=rel(2)),
	            axis.text=element_text(size=rel(1.5)),
	            axis.title=element_text(size=rel(1.6)))
	if (type == "boxes") {
		g <- ggplot2::geom_boxplot(ggplot2::aes(x=t,y=y,group=t),outlier.shape=NA) #outlier.size=0.3,outlier.color="green")
	} else {
		g <- ggplot2::geom_line(ggplot2::aes(x=t, y=y, group=sim),color="blue", alpha = 0.07, linewidth=1)
	}

	for(i in seq(length(experiments))){
		N <- dim(simulations[[i]]$func)[3]
		oNames <- names(experiments[[i]][["outputValues"]])
		if (is.null(var.names) && "initialState" %in% names(experiments[[i]])){
			xNames <- names(experiments[[i]][["initialState"]])
		} else if (is.null(var.names)){
			xNames <- rownames(simulations[[i]]$state)
		} else {
			xNames <- var.names
		}
		for(j in seq(num.of.funcs)){
			z <- experiments[[i]][["outputValues"]][[j]]
			dz <- experiments[[i]][["errorValues"]][[j]]
			if ("measurementTimes" %in% names(experiments[[i]])){
				t_ <- ttf(experiments[[i]][["measurementTimes"]])
			} else {
				t_ <- ttf(experiments[[i]][["outputTimes"]])
			}
			df.experiments <- data.frame(t=t_, y=z, upper=z+dz, lower=z-dz)
			tf <- ttf(experiments[[i]][["outputTimes"]])
			f <- as.numeric(simulations[[i]]$func[j,,])
			if (length(tf)==1){
				mle <- ggplot2::geom_point(data=data.frame(t=tf,y=simulations[[i]]$func[j,,MLE]), aes(x=t,y=y), color="magenta", size=2, inherit.aes=FALSE)
			} else {
				mle <- ggplot2::geom_line(data=data.frame(t=tf,y=simulations[[i]]$func[j,,MLE]), aes(x=t,y=y), color="magenta", linewidth=2, inherit.aes=FALSE)
			}
			if (is.null(yl.func)){
				YLIMIT <- NULL
			} else if (is.logical(yl.func[[j]]) && is.na(yl.func[[j]])){
				## the is.logical is important: is.numeric(NA) is FALSE, is.logical(NA) is TRUE. is.na(ylim()) leads to an error
				R <- range(c(z+dz,z-dz,simulations[[i]]$func[j,,MLE]),na.rm=TRUE); print(R)
				YLIMIT <- ggplot2::ylim(ifelse(diff(R)<1e-6,0.0,R[1]),R[2])
			} else if (is.numeric(yl.func[[j]])){
				Q <- quantile(f,probs=yl.func[[j]],na.rm=TRUE)
				print(Q)
				YLIMIT <- ggplot2::ylim(head(Q,1),tail(Q,1))
			} else {
				YLIMIT <- yl.func[[j]]
			}

			#cat(sprintf("[experiment %i, output function %i] N = %i; length(tf): %i, length(f): %i; tf; f\n",i,j,N,length(tf), length(f)))
			#print(tf)
			#print(f)
			df.simulations <- data.frame(t=rep(tf,N), y=f, sim=rep(seq(N),each=length(tf)))
			p[[(i-1)*M+j]] <- ggplot2::ggplot(df.simulations)+g+
				ggplot2::geom_errorbar(data=df.experiments, ggplot2::aes(x=t, y=y, ymin = lower, ymax = upper), color="red", inherit.aes=FALSE)+
				ggplot2::ggtitle(names(experiments[i]))+T1+
				ggplot2::labs(y=oNames[j],x=t_txt)+YLIMIT+mle
		}
		if (plot.states){
			for(j in seq(num.of.vars)){
				y <- as.numeric(simulations[[i]]$state[j,,])
				ty <- ttf(experiments[[i]][["outputTimes"]])
				if (length(ty)==1){
					mle <- ggplot2::geom_point(data=data.frame(t=ty,y=simulations[[i]]$state[j,,MLE]), aes(x=t,y=y), color="blue", size=2, inherit.aes=FALSE)
				} else {
					mle <- ggplot2::geom_line(data=data.frame(t=ty,y=simulations[[i]]$state[j,,MLE]), aes(x=t,y=y), color="blue", linewidth=2, inherit.aes=FALSE)
				}
				if (is.null(yl.state)){
					YLIMIT <- NULL
				} else if (is.logical(yl.state[[j]]) && is.na(yl.state[[j]])){
					R <- range(simulations[[i]]$state[j,,MLE],na.rm=TRUE)
					YLIMIT <- ggplot2::ylim(ifelse(diff(R)<1e-6,0.0,R[1]),R[2])
				} else if (is.numeric(yl.state[[j]]) && length(yl.state[[j]])==2){
					Q <- quantile(y,probs=yl.state[[j]],na.rm=TRUE)
					YLIMIT <- ggplot2::ylim(head(Q,1),tail(Q,1))
				} else {
					YLIMIT <- yl.state[[j]]
				}
				df.simulations <- data.frame(t=rep(ty,N), y=y, sim=rep(seq(N), each=length(ty)))
				p[[(i-1)*M+j+num.of.funcs]] <- ggplot2::ggplot(df.simulations,ggplot2::aes(x=t, y=y, group=sim))+g+
					ggplot2::ggtitle(names(experiments[i]))+T1+
					ggplot2::labs(y=xNames[j],x=t_txt)+YLIMIT+mle
			}
		}
	}
	n <- M
	m <- num.experiments
	return(gridExtra::marrangeGrob(grobs=p,ncol=m,nrow=n))
}


#' plot function for experiments
#'
#' This function does not use ggplot2, only base plot functions like
#' lines() and arrows().
#'
#' @param simulations simulation results (a list)
#' @param experiments experiment setup (a list)
#' @param nmax maximum index of lines to plot from sample, between 1 and size of sample
#' @param by skip this many sampled lines in between plotted lines
#' @param ylimit y-axis limits
#' @return plot object
#' @export
plotTimeSeriesBase <- function(simulations, experiments, nmax=NULL, by=1, ylimit=NULL){
	nf <- dim(simulations[[1]]$func)[1]
	par(mfrow=c(nf,length(experiments)))
	for (i in seq_along(experiments)){
		time <- experiments[[i]]$outputTimes
		n <- dim(simulations[[i]]$func)
		if (is.null(nmax)) nmax=n[3]
		oNames <- names(experiments[[i]]$outputValues)
		for (j in seq(n[1])){
			o <- experiments[[i]]$data[j,] %otherwise% experiments[[i]]$outputValues[[j]]
			e <- experiments[[i]]$stdv[j,] %otherwise% experiments[[i]]$errorValues[[j]]
			if (is.null(ylimit) || length(ylimit)<j || any(is.na(ylimit[[j]])) || !all(is.finite(ylimit[[j]]))){
				yl <- c(min(o-e-1e-1,na.rm=TRUE),max(o+e+1e-1,na.rm=TRUE))
			} else {
				yl <- ylimit[[j]]
			}
			p <- plot(time,o,type='p',ylim=yl,main=names(experiments)[i],ylab=oNames[j])
			arrows(time,o,time,o+e,angle=90,length=0.01)
			arrows(time,o,time,o-e,angle=90,length=0.01)
			for (k in seq(1,nmax,by=by)){
				y<-as.numeric(simulations[[i]]$func[j,,k])
				lines(time,y,col=rgb(0.3,0.6,0.9,0.4));
			}
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
#' @export
#' @param simulations list of simualtions as output from the simulator
#' @param experiments list of experiments
#' @param show.plot boolean variable. Set `show.plot=TRUE` to display plots
#'      when running the funcion, FALSE otherwise
#' @return list of plots with simulations and experimental data
plotTimeSeries <- function(simulations, experiments, show.plot = TRUE){
  num.experiments <- length(simulations)
  num.simulations <- dim(simulations[[1]]$func)[3]
  stopifnot(num.experiments == length(experiments))
  num.out.funcs <- dim(experiments[[1]]$outputValues)[2]

  p <- list()

  for(i in 1:num.experiments){
    experiment <- experiments[[i]]
    for(output.idx in 1:num.out.funcs){
      df.experiments <- data.frame(t=experiment[["outputTimes"]], y=experiment[["outputValues"]][[output.idx]])
      y <- c(simulations[[i]]$func[output.idx,,])
      df.simulations <- data.frame(t=rep(experiment[["outputTimes"]],num.simulations), y=y, sim=rep(1:num.simulations,each=length(experiment[["outputTimes"]])))
      p[[i]] <-
        ggplot(df.simulations,aes(x=t, y=y, group=sim)) +
        geom_line(color="blue", alpha = 0.1, size=1.5) +
        geom_point(data=df.experiments, aes(x=t, y=y), inherit.aes=FALSE) +
        ggtitle(paste0("Experiment: ", names(experiments[i]), "\nOutput: ", names(experiment[["outputValues"]])[output.idx]))
    }
  }
  if(show.plot)  show(do.call(gridExtra::grid.arrange,p))
  return(p)
}




#' Get the values of the input for a series of dose response experiments
#'
#' This function finds the vector of inputs that  varies among a series of experiments
#' that are part of a dose response experiment
#'
#' @export
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
#' @export
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
  num.out.funcs <- dim(experiments[[1]]$outputValues)[2]
  p <- list()
  for(j in 1:num.out.funcs){
    y.exp <- sapply(experiments, function(e) e$outputValues[[j]])
    E.data.frame.exp <- data.frame(x=dose,y=y.exp)
    y.sim <- sapply(simulations, function(o) o$func[j,1,])
    E.data.frame.sim <- data.frame(x=(rep(dose,each=num.simulations)),y=c(y.sim))
    p[[j]] <- ggplot(E.data.frame.exp, aes(x=x,y=y)) + geom_point(color="red")
    p[[j]] <-  p[[j]] +
      geom_point(data=E.data.frame.sim, aes(x=x,y=y), inherit.aes = FALSE, color="blue", alpha = 0.4)
    p[[j]] <-  p[[j]] +
      geom_point(data=E.data.frame.exp, aes(x=x,y=y)) +
      geom_point(color="red") +
      labs(x = comment(dose), y = names(experiments[[1]][["outputValues"]])[j])
  }
  return(p)
}

