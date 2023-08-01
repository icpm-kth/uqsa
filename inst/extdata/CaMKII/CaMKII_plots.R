require(ggplot2)

#' Plot simulations of CaMKII model using parameter draws
#'
#' This function simulates the CaMKII model with the parameters given 
#' in input (e.g. draws obtained through ABCMCMC). The simulated data are
#' plotted together with the experimental data from the different 
#' experimental datasets.
#'
#' @export
#' @param draws a matrix containing the CaMKII parameters as rows.
#'      this name will be used to find the file and
#' @param num.sub.samples number of parameters to sample from the 
#'      rows of matrix draws to use in the simulations
#' @param show.plot boolean variable. Set show.plot=TRUE to display plots
#'      when running the funcion, FALSE otherwise
#' @return list of plots, one for each experimental dataset

plotCaMKIISimulations <- function(draws, num.sub.samples = 100, show.plot = TRUE){
  
  if(dim(draws)[1] < num.sub.samples) num.sub.samples <- dim(draws)[1]
  
  model.tsv <- uqsa_example("CaMKII",full.names=TRUE)
  model.tab <- SBtabVFGEN::sbtab_from_tsv(model.tsv)
  source(uqsa_example("CaMKII",pat="^CaMKII.*R$",full.names=TRUE))
  modelName <- checkModel(comment(model.tab), uqsa_example("CaMKII",pat="_gvf[.]c$"))
  experiments <- SBtabVFGEN::sbtab.data(model.tab)
  parMap <- function(parABC){
    return(10^parABC)
  }
  experimentsIndices <- list(
    which(startsWith(names(experiments),"E0")),
    which(startsWith(names(experiments),"E1")),
    which(startsWith(names(experiments),"E2")),
    which(startsWith(names(experiments),"E3")),
    which(startsWith(names(experiments),"E4")),
    which(startsWith(names(experiments),"E5"))
  )
  
  outputIdx <- c(2,1,3,4,4,1)
  inputIdx <- c(1,1,4,1,1,1)
  
  outputNames <- row.names(model.tab$Output)
  inputNames <- row.names(model.tab$Input)
  
  p <- list()
  for(i in 1:6){
    outputIndex <- outputIdx[i]
    inputIndex <- inputIdx[i]
    
    x <- sapply(experiments[experimentsIndices[[i]]], function(e) e$input[inputIndex])
    y_exp <- sapply(experiments[experimentsIndices[[i]]], function(e) e$outputValues[[outputIndex]])
    E.data.frame.exp <- data.frame(x=x,y=y_exp)
    
    p[[i]] <- ggplot(E.data.frame.exp, aes(x=x,y=y)) + geom_point(color="red")
    
    num.sub.samples <- 100
    simulate <- simulator.c(experiments[experimentsIndices[[i]]], modelName, parMap)
    #objectiveFunction <- makeObjective(experiments[experimentsIndices[[i]]], modelName, distanceMeasure, parMap, simulate)
    #apply(objectiveFunction(t(draws[sample(1:dim(draws)[1],num.sub.samples),])),2,max)
    
    output_yy <- simulate(t(draws[sample(1:dim(draws)[1],num.sub.samples),]))
    y <- sapply(output_yy, function(o) o$func[outputIndex,1,])
    E.data.frame.sim <- data.frame(x=(rep(x,each=num.sub.samples)),y=c(y))
    p[[i]] <-  p[[i]] + geom_point(data=E.data.frame.sim, aes(x=x,y=y), inherit.aes = FALSE, color="blue", alpha = 0.4)
    p[[i]] <-  p[[i]] + geom_point(data=E.data.frame.exp, aes(x=x,y=y)) + geom_point(color="red") + labs(x = inputNames[inputIndex], y = outputNames[outputIndex])
    
  }
  if(show.plot)  show(do.call(gridExtra::grid.arrange,p))
  return(p)
}
