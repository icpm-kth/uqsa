Experiments.sim <-function() {
    library(SBtabVFGEN)
    library(rgsl)
    sb <- sbtab_from_tsv(dir(pattern="[.]tsv$"))
    # vectors
    k <- sbtab_quantity(sb$Parameter)
    u <- sbtab_quantity(sb$Input)
    y0 <- sbtab_quantity(sb$Compound)
    # matrices
    Y0<- update_from_table(y0,sb$Experiments)
    K <- update_from_table(k,sb$Experiments)
    U <- update_from_table(u,sb$Experiments)
    # conservation law adjustments
    if (file.exists("ConservationLaws.h5")){
        h5f <- H5File$new("ConservationLaws.h5",mode="r")
        # conservation laws must apply for _this_ model
        DocName <- h5f[['Document']]$read()
        stopifnot(DocName==comment(sb))
        # conserved constants are treated as additionl inputs
        CL <- h5f[['ConservationLaws']]$read()
        ConservedConst <- CL %*% Y0
        U <- rbind(U,ConservedConst)
        # reduce the dim of the state space
        eliminate <- -1-h5f[['EliminatedCompounds']]$read()
        Y0<-Y0[eliminate,]
        h5f$close_all()
    }
    P <- rbind(K,U)
    print(P)
    print(Y0)
    ny <- nrow(Y0)
    N <- rownames(Y0)
    t <- seq(-15,600,by=5)
    cpu.t <- Sys.time()
    Y <- r_gsl_odeiv2("AKAP79",t=t,y0=Y0,p=P)
    cpu.t <- Sys.time() - cpu.t
    print(cpu.t)
    
    for (i in 1:length(sb$Experiments)){
        plot(t,Y[ny,,i],ylab=N[i],xlab='time')
    }
}
