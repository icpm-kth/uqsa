runJoao <-  function(i){
  setwd("/Users/fedmil/Documents/AKAP79_PKA/R")
  
  source('ScoringFunction_LU.R')
  source('LoadTargets_AKAP79.R')
  source('runModel_AKAP79.R')
  source('copulaFunctions_LU.R')
  source('ABCMCMCFunctions_LU.R')
  source('PreCalibration_LU.R')
  source('PKAModel.R')
  library(deSolve)
  library(ks)
  library(VineCopula)
  library(MASS)
  library(R.utils)
  library(R.matlab) #closeAllConnections()
  library(stringr)
  library(doMC)
  registerDoMC(19) # example with 20 cores throughout
  getDoParWorkers()
  
  
  
  # Load Parameters and Rii Amounts
  matdata_Parms <- read.table("pkaParms_Restri.txt", header = TRUE);
  parNames = matdata_Parms$Name.;
  parVal = matdata_Parms$Value.;
  R_Amount = 6.93 #OE: Is this used? If not remove
  
  #Experiments to run
  exp_idxs<-c(7, 13, 19, 22, 8, 14, 20, 23, 9, 15, 21, 24)
  exp_idx <- exp_idxs[i]
  # 7: 0 No_CaN
  # 8: 0 With_CaN
  # 9: 0 CaN_AKAP
  # 10: 0.1 No_CaN
  # 11: 0.1 With_CaN
  # 12: 0.1 CaN_AKAP
  # 13: 0.2 No_CaN *
  # 14: 0.2 With_CaN *
  # 15: 0.2 CaN_AKAP *
  # 16: 0.5  No_CaN
  # 17: 0.5 With_CaN
  # 18: 0.5 CaN_AKAP
  # 19:1  No_CaN *
  # 20:1 With_CaN *
  # 21:1 CaN_AKAP *
  # 22:2  No_CaN *
  # 23:2 With_CaN *
  # 24:2 CaN_AKAP *
  
  #Input values for the different experiments
  xAll <- c(-1,-1, -1,-1,-1,-1, 0, 0, 0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 1, 1, 1, 2, 2, 2)
  
  #Exp type
  exp_types<-c("","","","","","","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP","SS_No_CaN", "SS_with_CaN", "SS_CaN_AKAP")
  #OE: remove SS in exp_types
  
  
  
  
  #Initial conditions
  species_iconc_TS <- c(Rii = 6.3, cAMP = 0, RiiP = 0, Rii_C = 0.63,RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                        Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 1.5, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0,Total_C = 0.63, cycle_Rii = 0,
                        Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0);#, p_AKAP = 0); #Total C should be removed?
  
  #OE: Remove the ones below because they are not used
  species_iconc_CC <- c(Rii =0, cAMP = 0, RiiP = 0, Rii_C = 0, RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                        Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                        Cycle_RiiP = 0,Thr_unphos_Rii  = 0,Thr_phos_Rii = 0, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);
  
  species_iconc_MutA <- c(Rii = 6.3, cAMP = 1, RiiP = 0, Rii_C = 0.63,RiiP_cAMP = 0,RiiP_C = 0,RiiP_C_cAMP = 0,C = 0,
                          Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                          Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);
  
  species_iconc_MutE <- c(Rii = 0, cAMP = 1, RiiP = 6.3, Rii_C = 0,RiiP_cAMP = 0,RiiP_C = 0.63,RiiP_C_cAMP = 0,C = 0,
                          Rii_cAMP = 0,Rii_C_cAMP = 0, Total_RII = R_Amount, CaN = 0, RiiP_CaN = 0, RiiP_cAMP_CaN  = 0, Total_C = 0.63, cycle_Rii = 0,
                          Cycle_RiiP = 0,Thr_unphos_Rii  = 1,Thr_phos_Rii = 1, b_AKAP = 0, AKAR4 = 0.2, AKAR4_C = 0, AKAR4p = 0, p_AKAP = 0);
  
  
  # Set ll and up for standard parameters
  scale <- 1000
  
  #Clean up below and change to correct values
  ll = c( parVal[1:6]/scale, parVal[7]/scale, parVal[8:9]/scale, parVal[10]/scale, parVal[11:20]/scale, parVal[21]/1.9, parVal[22:24]/scale,  parVal[25:26]/1.25, parVal[27]/1.25,  parVal[28:29]/1.5,  parVal[30]/2);
  ul = c( parVal[1:6]*scale, parVal[7]*scale, parVal[8:9]*scale, parVal[10]*scale, parVal[11:20]*scale, parVal[21]*1.9, parVal[22:24]*scale,  parVal[25:26]*1.25, parVal[27]*1.25,  parVal[28:29]*1.5,  parVal[30]*2);
  
  
  
  ll = log10(ll);
  ul = log10(ul);
  
  
  parIdx <- c(1:30); 
  
  # load targets
  out <- loadTargets_TSCC()
  xtarget <- out$xtarget
  ytarget <- out$ytarget
  rm(out)
  
  ytarget_min<-0 #Zero phosphorylation
  ytarget_max<-0.2 #Full phosphorylation
  yy_min<-0 #Zero phosphorylation
  yy_max<-0.2  #Full phosphorylation
  
  
  # no of samples
  
  ns <- 10 # no of samples required from each ABC-MCMC chain #WAS 1000
  npc <- 50 # pre-calibration  WAS 50.000
  
  
  # settings
  p <- 0.01 # was 0.01 #OE: Why was this changed?
  #nChains <- 1 #was 20
  #delta <- 0.01 #was 0.1
  delta <- c(-1,-1, -1,-1,-1,-1, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
  delta <- delta
  nChains <- 19 #was 20
  
  
  set.seed(7619201)
  
  
  cat(sprintf("#####Starting run for Exp. Dataset %i #####\n", exp_idx))
  
  
  
  cat(sprintf("-Fitting independent Copula \n"))
  out <- makeIndepCopula_JA(ll, ul)
  
  copula <- out$copula
  U <- out$U
  Z <- out$Z
  Y <- out$Y
  
  ### From the preCalibration function in file 'PreCalibration_LU.R
  input <- xAll[[exp_idx]]
  parDef <- parVal
  ytarget_ <- ytarget[[exp_idx]]
  rInd <- exp_types[exp_idx]
  amounts <- species_iconc_TS
  datatype <- "T_S"
  
  
  np <- length(parIdx) #nr of parameters
  R <- RVineSim(npc, copula) #gets different params values from hyper structure 'copula'
  prePar <- matrix(0, npc, np)
  for(i in 1:np){
    prePar[,i] = spline(Z[,i],U[,i],xout=R[,i])$y #'fusion' of different params, i.e., join information from normal distribution & CDF
  }
  
  tpar <- prePar[1,]
  out <- runModel(tpar, parIdx, parDef, input, rInd, amounts, datatype)
  output_Joao<-out$yy
  
  return(list(outputJoao = output_Joao, tpar = tpar, amounts = amounts, matdata_Parms = matdata_Parms))
}

runNew <- function(i, tpar, matdata_Parms){
  setwd("~/Documents/uqsa/AKAP79/")
  #remotes::install_github("a-kramer/rgsl")
  #remotes::install_github("a-kramer/SBtabVFGEN")
  library(rgsl)
  library(SBtabVFGEN)
  library(UQ)
  source("../UQ/R/ABCMCMCFunctions.R")
  source("../UQ/R/copulaFunctions.R")
  source("../UQ/R/import_from_SBtab.R")
  source("../UQ/R/PreCalibration.R")
  source("../UQ/R/prior.R")
  source("../UQ/R/runModel.R")
  
  SBtabDir <- getwd()
  model = import_from_SBtab(SBtabDir)
  #modelName <- checkModel(comment(model),paste0(comment(model),'.R'))
  modelName <- checkModel(comment(model),paste0(comment(model),'_gvf.c'))
  
  #source(paste(SBtabDir,"/",modelName,".R",sep=""))
  
  parVal <- model[["Parameter"]][["!DefaultValue"]]
  parNames <- model[["Parameter"]][["!Name"]]
  
  # load experiments
  experiments <- import_experiments(modelName, SBtabDir)
  
  # scale to determine prior values
  defRange <- 1000
  
  # Define Lower and Upper Limits for logUniform prior distribution for the parameters
  ll <- c(parVal[1:19]/defRange, parVal[20]/1.9, parVal[21]/defRange, parVal[22:24]/1.25, parVal[25:26]/1.5, parVal[27]/2)
  ul <- c(parVal[1:19]*defRange, parVal[20]*1.9, parVal[21]*defRange, parVal[22:24]*1.25, parVal[25:26]*1.5, parVal[27]*2)
  ll = log10(ll) # log10-scale
  ul = log10(ul) # log10-scale
  
  
  # Define the experiments that have to be considered in each iteration of the for loop to compare simulations with experimental data
  experimentsIndices <- c(3, 12, 18, 9, 2, 11, 17, 8, 1, 10, 16, 7)

  expInd <- experimentsIndices[i]
  
  # Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
  npc <- 500 # pre-calibration
  
  # Define the number of Cores for the parallelization
  nCores <- parallel::detectCores() - 1 #%/% 2
  
  paramNames_new_AKAP79 <- c("kf_Rii_C__RiiP_C",
                             "kf_RiiP_CxcAMP__RiiP_C_cAMP",
                             "kf_RiiP_cAMPxC__RiiP_C_cAMP",
                             "kb_RiiP_cAMPxC__RiiP_C_cAMP",
                             "kb_RiiPXcAMP__RiiP_cAMP",
                             "kf_RiiPXcAMP__RiiP_cAMP",
                             "kf_RiiPxC__RiiP_C",
                             "kb_RiiPxC__RiiP_C",
                             "kf_cAMPxRii__Rii_cAMP",
                             "kb_cAMPxRii__Rii_cAMP",
                             "kf_Rii_CxcAMP__Rii_C_cAMP",
                             "kb_Rii_CxcAMP__Rii_C_cAMP",
                             "kf_RiixC__Rii_C",
                             "kf_Rii_cAMPxC__Rii_C_cAMP",
                             "kb_Rii_cAMPxC__Rii_C_cAMP",
                             "kf_Rii_C_cAMP__RiiP_C_cAMP",
                             "kb_RiixC__Rii_C",
                             "AKAPoff_1",
                             "AKAPoff_3",
                             "AKAPon_1",
                             "AKAPon_3",
                             "kf_C_AKAR4",
                             "kb_C_AKAR4",
                             "kcat_AKARp",
                             "kmOFF",
                             "kmON",
                             "KD_T")
  idx_paramAKAP_in_param <- match(paramNames_new_AKAP79,matdata_Parms$Name.)
  idx_paramAKAP_in_param[18] <- 18
  idx_paramAKAP_in_param[19] <- 20
  idx_paramAKAP_in_param[20] <- 21
  idx_paramAKAP_in_param[21] <- 23
  idx_paramAKAP_in_param[22] <- 25
  idx_paramAKAP_in_param[23] <- 26
  idx_paramAKAP_in_param[24] <- 27
  
  parMap <- function(parABC){
    parABC <- parABC[idx_paramAKAP_in_param]
    return(10^parABC)
  }
  
  output_yy_new <- runModel(experiments[expInd], modelName,  tpar, parMap)
  return(list(output_yy_new = output_yy_new, initialState = experiments[expInd][[1]]$initialState))
}

source("/Users/fedmil/Documents/uqsa/UQ/R/import_from_SBtab.R")
SBtabDir <- "~/Documents/uqsa/AKAP79"
modelName <- "AKAP79"
experiments <- import_experiments(modelName, SBtabDir)
experimentsIndices <- c(3, 12, 18, 9, 2, 11, 17, 8, 1, 10, 16, 7)
for(i in 1:length(experimentsIndices)){
  out <- runJoao(i)
  
  outputJoao <- out$outputJoao
  tpar <- out$tpar
  amounts <- out$amounts
  matdata_Parms <- out$matdata_Parms
  plot(outputJoao)
  
  out <- runNew(i, tpar, matdata_Parms)
  out_new <- out$output_yy_new
  initialState_new <- out$initialState
  lines(out_new[[1]],col='red')
}
