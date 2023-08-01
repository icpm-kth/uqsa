require(ggplot2)

#' Plot simulations of AKAP79 model using parameter draws
#'
#' This function simulates the AKAP79 model with the parameters given 
#' in input (e.g. draws obtained through ABCMCMC). The simulated data are
#' plotted together with the experimental data from the different 
#' experimental datasets.
#'
#' @export
#' @param draws a matrix containing the AKAP79 parameters as rows.
#'      this name will be used to find the file and
#' @param num.sub.samples number of parameters to sample from the 
#'      rows of matrix draws to use in the simulations
#' @param show.plot boolean variable. Set show.plot=TRUE to display plots
#'      when running the funcion, FALSE otherwise
#' @return list of plots, one for each experimental dataset


plotAKAP79Simulations <- function(draws, num.sub.samples = 1000, show.plot = TRUE){
  
  if(dim(draws)[1] < num.sub.samples) num.sub.samples <- dim(draws)[1]
  
  model.tsv <- uqsa_example("AKAP79")
  model.tab <- sbtab_from_tsv(model.tsv)
  source(uqsa_example("AKAP79",pat="^AKAP79[.]R$"))
  modelName <- checkModel(comment(model.tab),uqsa_example("AKAP79",pat="_gvf.c"))
  experiments <- sbtab.data(model.tab)
  parMap <- function(parABC){
    return(10^parABC)
  }
  
  p<-list()
  for(i in 1:12){
    simulate <- simulator.c(experiments[i], modelName, parMap)
    #output_yy <- simulate(t(ABCMCMCoutput$draws[(ABCMCMCoutput$scores %in% sort(ABCMCMCoutput$scores)[1:300]),]))
    output_yy <- simulate(t(draws[sample(1:dim(draws)[1],num.sub.samples),]))
    
    experiment <- experiments[[i]]
    df.experiments <- data.frame(t=experiment[["outputTimes"]], y=experiment[["outputValues"]][[1]])
    
    y <- c(output_yy[[1]]$func[1,,])
    df.simulations <- data.frame(t=rep(experiment[["outputTimes"]],num.sub.samples), y=y, sim=rep(1:num.sub.samples,each=length(experiment[["outputTimes"]])))
    
    
    # df <- mclapply(1:dim(output_yy[[1]][["func"]])[3], function(i) as.data.frame(x = list(output_yy[[1]][["func"]][1,,i],experiment[['outputTimes']]), col.names = c("y","t")))
    # 
    # df <- reshape2::melt(df,id=c("t","y"))
    # 
    p[[i]] <- 
      ggplot(df.simulations,aes(x=t, y=y, group=sim))+
      geom_line(color="blue", alpha = 0.1)+
      geom_point(data=df.experiments, aes(x=t, y=y), inherit.aes=FALSE)
  }
  if(show.plot)  show(do.call(gridExtra::grid.arrange,p))
  return(p)
}
