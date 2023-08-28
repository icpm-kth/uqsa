require(ggplot2)

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
        ggplot(df.simulations,aes(x=t, y=y, group=sim))+
        geom_line(color="blue", alpha = 0.1, size=1.5)+
        geom_point(data=df.experiments, aes(x=t, y=y), inherit.aes=FALSE)+
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
#'
#' @export
#' @param experiments list of experimental data from the same dose response experiment
#' @return vector of inputs (i.e. dose) that varies among the experiments provided to the function.
#'         The name of the input is saved in "comment(dose)".
#' 

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
#'
#' @export
#' @param simulations list of simualtions as output from the simulator
#' @param experiments list of experimental data from the same dose response experiment
#' @param dose vector of dose values to plot on the x axis
#' @param show.plot boolean variable. Set show.plot=TRUE to display plots
#'      when running the funcion, FALSE otherwise
#' @return plot with simulations and experimental data
#' 
 
plotDoseResponse <- function(simulations, experiments, dose, show.plot = TRUE){
  
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
    p[[j]] <-  p[[j]] + geom_point(data=E.data.frame.sim, aes(x=x,y=y), inherit.aes = FALSE, color="blue", alpha = 0.4)
    p[[j]] <-  p[[j]] + geom_point(data=E.data.frame.exp, aes(x=x,y=y)) + geom_point(color="red") + labs(x = comment(dose), y = names(experiments[[1]][["outputValues"]])[j])
    
  }
  return(p)
}

