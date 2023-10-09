require("deSolve")

# ode vector field: y'=f(t,y;p)
AKAP79_vf <- function(t, state, parameters)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	RiiP <- state[1]
	RiiP_cAMP <- state[2]
	RiiP_C <- state[3]
	RiiP_C_cAMP <- state[4]
	C <- state[5]
	Rii_cAMP <- state[6]
	Rii_C_cAMP <- state[7]
	RiiP_CaN <- state[8]
	RiiP_cAMP_CaN <- state[9]
	AKAR4_C <- state[10]
	AKAR4p <- state[11]
	AKAR4 <- AKAR4_ConservedConst - (AKAR4_C+AKAR4p)
	CaN <- CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN)
	Rii_C <- Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	cAMP <- cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)
	Rii <- Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN <- b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP <- kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	reaction_51 <- kf_Rii_C__RiiP_C*Rii_C
	reaction_14 <- kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	reaction_12 <- kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_43 <- kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	reaction_23 <- kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_78 <- kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	reaction_56 <- kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	reaction_76 <- kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	reaction_62 <- kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	reaction_58 <- kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	reaction_44 <- kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	reaction_33 <- kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	reaction_48 <- kf_RiiP_CaN__RiixCaN*RiiP_CaN
	reaction_37 <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	f_<-vector(mode='numeric',len=11)
	f_[1] <- -reaction_14-reaction_43-reaction_44
	f_[2] <- +reaction_43-reaction_23-reaction_33
	f_[3] <- +reaction_51+reaction_14-reaction_12
	f_[4] <- +reaction_12+reaction_23+reaction_62
	f_[5] <- -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2
	f_[6] <- +reaction_78-reaction_76+reaction_37
	f_[7] <- +reaction_56+reaction_76-reaction_62
	f_[8] <- +reaction_44-reaction_48
	f_[9] <- +reaction_33-reaction_37
	f_[10] <- +reaction_1-reaction_2
	f_[11] <- +reaction_2
	names(f_) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
## for some weird reason deSolve wants this to be a list:
	return(list(f_))
}
# ode Jacobian df(t,y;p)/dy
AKAP79_jac<-function(t, state, parameters)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	RiiP <- state[1]
	RiiP_cAMP <- state[2]
	RiiP_C <- state[3]
	RiiP_C_cAMP <- state[4]
	C <- state[5]
	Rii_cAMP <- state[6]
	Rii_C_cAMP <- state[7]
	RiiP_CaN <- state[8]
	RiiP_cAMP_CaN <- state[9]
	AKAR4_C <- state[10]
	AKAR4p <- state[11]
	AKAR4 <- AKAR4_ConservedConst - (AKAR4_C+AKAR4p)
	CaN <- CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN)
	Rii_C <- Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	cAMP <- cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)
	Rii <- Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN <- b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP <- kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	reaction_51 <- kf_Rii_C__RiiP_C*Rii_C
	reaction_14 <- kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	reaction_12 <- kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_43 <- kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	reaction_23 <- kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_78 <- kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	reaction_56 <- kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	reaction_76 <- kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	reaction_62 <- kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	reaction_58 <- kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	reaction_44 <- kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	reaction_33 <- kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	reaction_48 <- kf_RiiP_CaN__RiixCaN*RiiP_CaN
	reaction_37 <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	jac_ <- matrix(NA,11,11)
# column 1 (df/dy_0)
	jac_[1,1] <- (-((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-C*kf_RiiPxC__RiiP_C-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiPXcAMP__RiiP_cAMP
	jac_[2,1] <- (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiPXcAMP__RiiP_cAMP
	jac_[3,1] <- C*kf_RiiPxC__RiiP_C
	jac_[4,1] <- 0
	jac_[5,1] <- C*kf_RiixC__Rii_C-C*kf_RiiPxC__RiiP_C
	jac_[6,1] <- -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP
	jac_[7,1] <- 0
	jac_[8,1] <- ((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[9,1] <- 0
	jac_[10,1] <- 0
	jac_[11,1] <- 0
# column 2 (df/dy_1)
	jac_[1,2] <- RiiP*kf_RiiPXcAMP__RiiP_cAMP+kb_RiiPXcAMP__RiiP_cAMP
	jac_[2,2] <- (-((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-C*kf_RiiP_cAMPxC__RiiP_C_cAMP-RiiP*kf_RiiPXcAMP__RiiP_cAMP-kb_RiiPXcAMP__RiiP_cAMP
	jac_[3,2] <- RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[4,2] <- C*kf_RiiP_cAMPxC__RiiP_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[5,2] <- C*kf_RiixC__Rii_C-C*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[6,2] <- (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP)-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP
	jac_[7,2] <- -((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,2] <- 0
	jac_[9,2] <- ((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[10,2] <- 0
	jac_[11,2] <- 0
# column 3 (df/dy_2)
	jac_[1,3] <- kb_RiiPxC__RiiP_C
	jac_[2,3] <- 0
	jac_[3,3] <- (-kf_Rii_C__RiiP_C)-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiP_CxcAMP__RiiP_C_cAMP-kb_RiiPxC__RiiP_C
	jac_[4,3] <- (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[5,3] <- kb_RiiPxC__RiiP_C-kb_RiixC__Rii_C
	jac_[6,3] <- 0
	jac_[7,3] <- -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,3] <- 0
	jac_[9,3] <- 0
	jac_[10,3] <- 0
	jac_[11,3] <- 0
# column 4 (df/dy_3)
	jac_[1,4] <- RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[2,4] <- kb_RiiP_cAMPxC__RiiP_C_cAMP-RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[3,4] <- (-kf_Rii_C__RiiP_C)+RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP+KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[4,4] <- (-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP)-KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP-kb_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[5,4] <- kb_RiiP_cAMPxC__RiiP_C_cAMP-kb_RiixC__Rii_C
	jac_[6,4] <- -((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP
	jac_[7,4] <- (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP)-((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,4] <- 0
	jac_[9,4] <- 0
	jac_[10,4] <- 0
	jac_[11,4] <- 0
# column 5 (df/dy_4)
	jac_[1,5] <- -RiiP*kf_RiiPxC__RiiP_C
	jac_[2,5] <- -RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[3,5] <- RiiP*kf_RiiPxC__RiiP_C-kf_Rii_C__RiiP_C
	jac_[4,5] <- RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[5,5] <- (-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_RiixC__Rii_C)-C*kf_RiixC__Rii_C-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-RiiP*kf_RiiPxC__RiiP_C-RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP-((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*kf_C_AKAR4-kb_RiixC__Rii_C
	jac_[6,5] <- (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP
	jac_[7,5] <- Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,5] <- 0
	jac_[9,5] <- 0
	jac_[10,5] <- ((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*kf_C_AKAR4
	jac_[11,5] <- 0
# column 6 (df/dy_5)
	jac_[1,6] <- RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[2,6] <- -RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[3,6] <- RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[4,6] <- -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[5,6] <- C*kf_RiixC__Rii_C-C*kf_Rii_cAMPxC__Rii_C_cAMP
	jac_[6,6] <- (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP)-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP-C*kf_Rii_cAMPxC__Rii_C_cAMP-kb_cAMPxRii__Rii_cAMP
	jac_[7,6] <- C*kf_Rii_cAMPxC__Rii_C_cAMP-((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,6] <- 0
	jac_[9,6] <- 0
	jac_[10,6] <- 0
	jac_[11,6] <- 0
# column 7 (df/dy_6)
	jac_[1,7] <- RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[2,7] <- -RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[3,7] <- RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-kf_Rii_C__RiiP_C
	jac_[4,7] <- kf_Rii_C_cAMP__RiiP_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[5,7] <- kb_Rii_cAMPxC__Rii_C_cAMP-kb_RiixC__Rii_C
	jac_[6,7] <- kb_Rii_cAMPxC__Rii_C_cAMP-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP
	jac_[7,7] <- (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP)-((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP-kf_Rii_C_cAMP__RiiP_C_cAMP-kb_Rii_cAMPxC__Rii_C_cAMP-kb_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,7] <- 0
	jac_[9,7] <- 0
	jac_[10,7] <- 0
	jac_[11,7] <- 0
# column 8 (df/dy_7)
	jac_[1,8] <- RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)+AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP)
	jac_[2,8] <- RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[3,8] <- 0
	jac_[4,8] <- 0
	jac_[5,8] <- C*kf_RiixC__Rii_C
	jac_[6,8] <- -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP
	jac_[7,8] <- 0
	jac_[8,8] <- (-RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-AKAPon_3*b_AKAP-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP)
	jac_[9,8] <- -RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[10,8] <- 0
	jac_[11,8] <- 0
# column 9 (df/dy_8)
	jac_[1,9] <- RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)+RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[2,9] <- RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)-RiiP*kf_RiiPXcAMP__RiiP_cAMP+AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP)
	jac_[3,9] <- RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[4,9] <- -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[5,9] <- C*kf_RiixC__Rii_C
	jac_[6,9] <- (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP)-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP+AKAPon_1*b_AKAP+AKAPoff_1*(1-b_AKAP)
	jac_[7,9] <- -((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,9] <- -RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jac_[9,9] <- (-RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-AKAPon_3*b_AKAP-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP)
	jac_[10,9] <- 0
	jac_[11,9] <- 0
# column 10 (df/dy_9)
	jac_[1,10] <- 0
	jac_[2,10] <- 0
	jac_[3,10] <- -kf_Rii_C__RiiP_C
	jac_[4,10] <- 0
	jac_[5,10] <- (-C*kf_RiixC__Rii_C)+C*kf_C_AKAR4+kcat_AKARp-kb_RiixC__Rii_C+kb_C_AKAR4
	jac_[6,10] <- (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP
	jac_[7,10] <- -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[8,10] <- 0
	jac_[9,10] <- 0
	jac_[10,10] <- (-C*kf_C_AKAR4)-kcat_AKARp-kb_C_AKAR4
	jac_[11,10] <- kcat_AKARp
# column 11 (df/dy_10)
	jac_[1,11] <- 0
	jac_[2,11] <- 0
	jac_[3,11] <- 0
	jac_[4,11] <- 0
	jac_[5,11] <- C*kf_C_AKAR4
	jac_[6,11] <- 0
	jac_[7,11] <- 0
	jac_[8,11] <- 0
	jac_[9,11] <- 0
	jac_[10,11] <- -C*kf_C_AKAR4
	jac_[11,11] <- 0
	rownames(jac_) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
	colnames(jac_) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
	return(jac_)
}
# ode parameter Jacobian df(t,y;p)/dp
AKAP79_jacp<-function(t, state, parameters)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	RiiP <- state[1]
	RiiP_cAMP <- state[2]
	RiiP_C <- state[3]
	RiiP_C_cAMP <- state[4]
	C <- state[5]
	Rii_cAMP <- state[6]
	Rii_C_cAMP <- state[7]
	RiiP_CaN <- state[8]
	RiiP_cAMP_CaN <- state[9]
	AKAR4_C <- state[10]
	AKAR4p <- state[11]
	AKAR4<-AKAR4_ConservedConst - (AKAR4_C+AKAR4p)
	CaN<-CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN)
	Rii_C<-Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	cAMP<-cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)
	Rii<-Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN<-b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN<-b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN<-b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN<-b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP<-kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	reaction_51<-kf_Rii_C__RiiP_C*Rii_C
	reaction_14<-kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	reaction_12<-kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_43<-kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	reaction_23<-kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_78<-kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	reaction_56<-kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	reaction_76<-kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	reaction_62<-kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	reaction_58<-kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	reaction_44<-kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	reaction_33<-kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	reaction_48<-kf_RiiP_CaN__RiixCaN*RiiP_CaN
	reaction_37<-kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	jacp_<-matrix(NA,11,33)
# column 1 (df/dp_1)
	jacp_[1,1] <- 0
	jacp_[2,1] <- 0
	jacp_[3,1] <- (-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C
	jacp_[4,1] <- 0
	jacp_[5,1] <- 0
	jacp_[6,1] <- 0
	jacp_[7,1] <- 0
	jacp_[8,1] <- 0
	jacp_[9,1] <- 0
	jacp_[10,1] <- 0
	jacp_[11,1] <- 0
# column 2 (df/dp_2)
	jacp_[1,2] <- 0
	jacp_[2,2] <- 0
	jacp_[3,2] <- KD_T*RiiP_C_cAMP-RiiP_C*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)
	jacp_[4,2] <- RiiP_C*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)-KD_T*RiiP_C_cAMP
	jacp_[5,2] <- 0
	jacp_[6,2] <- 0
	jacp_[7,2] <- 0
	jacp_[8,2] <- 0
	jacp_[9,2] <- 0
	jacp_[10,2] <- 0
	jacp_[11,2] <- 0
# column 3 (df/dp_3)
	jacp_[1,3] <- 0
	jacp_[2,3] <- -C*RiiP_cAMP
	jacp_[3,3] <- 0
	jacp_[4,3] <- C*RiiP_cAMP
	jacp_[5,3] <- -C*RiiP_cAMP
	jacp_[6,3] <- 0
	jacp_[7,3] <- 0
	jacp_[8,3] <- 0
	jacp_[9,3] <- 0
	jacp_[10,3] <- 0
	jacp_[11,3] <- 0
# column 4 (df/dp_4)
	jacp_[1,4] <- 0
	jacp_[2,4] <- RiiP_C_cAMP
	jacp_[3,4] <- 0
	jacp_[4,4] <- -RiiP_C_cAMP
	jacp_[5,4] <- RiiP_C_cAMP
	jacp_[6,4] <- 0
	jacp_[7,4] <- 0
	jacp_[8,4] <- 0
	jacp_[9,4] <- 0
	jacp_[10,4] <- 0
	jacp_[11,4] <- 0
# column 5 (df/dp_5)
	jacp_[1,5] <- RiiP_cAMP
	jacp_[2,5] <- -RiiP_cAMP
	jacp_[3,5] <- 0
	jacp_[4,5] <- 0
	jacp_[5,5] <- 0
	jacp_[6,5] <- 0
	jacp_[7,5] <- 0
	jacp_[8,5] <- 0
	jacp_[9,5] <- 0
	jacp_[10,5] <- 0
	jacp_[11,5] <- 0
# column 6 (df/dp_6)
	jacp_[1,6] <- -RiiP*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)
	jacp_[2,6] <- RiiP*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)
	jacp_[3,6] <- 0
	jacp_[4,6] <- 0
	jacp_[5,6] <- 0
	jacp_[6,6] <- 0
	jacp_[7,6] <- 0
	jacp_[8,6] <- 0
	jacp_[9,6] <- 0
	jacp_[10,6] <- 0
	jacp_[11,6] <- 0
# column 7 (df/dp_7)
	jacp_[1,7] <- -C*RiiP
	jacp_[2,7] <- 0
	jacp_[3,7] <- C*RiiP
	jacp_[4,7] <- 0
	jacp_[5,7] <- -C*RiiP
	jacp_[6,7] <- 0
	jacp_[7,7] <- 0
	jacp_[8,7] <- 0
	jacp_[9,7] <- 0
	jacp_[10,7] <- 0
	jacp_[11,7] <- 0
# column 8 (df/dp_8)
	jacp_[1,8] <- RiiP_C
	jacp_[2,8] <- 0
	jacp_[3,8] <- -RiiP_C
	jacp_[4,8] <- 0
	jacp_[5,8] <- RiiP_C
	jacp_[6,8] <- 0
	jacp_[7,8] <- 0
	jacp_[8,8] <- 0
	jacp_[9,8] <- 0
	jacp_[10,8] <- 0
	jacp_[11,8] <- 0
# column 9 (df/dp_9)
	jacp_[1,9] <- 0
	jacp_[2,9] <- 0
	jacp_[3,9] <- 0
	jacp_[4,9] <- 0
	jacp_[5,9] <- 0
	jacp_[6,9] <- ((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)
	jacp_[7,9] <- 0
	jacp_[8,9] <- 0
	jacp_[9,9] <- 0
	jacp_[10,9] <- 0
	jacp_[11,9] <- 0
# column 10 (df/dp_10)
	jacp_[1,10] <- 0
	jacp_[2,10] <- 0
	jacp_[3,10] <- 0
	jacp_[4,10] <- 0
	jacp_[5,10] <- 0
	jacp_[6,10] <- -Rii_cAMP
	jacp_[7,10] <- 0
	jacp_[8,10] <- 0
	jacp_[9,10] <- 0
	jacp_[10,10] <- 0
	jacp_[11,10] <- 0
# column 11 (df/dp_11)
	jacp_[1,11] <- 0
	jacp_[2,11] <- 0
	jacp_[3,11] <- 0
	jacp_[4,11] <- 0
	jacp_[5,11] <- 0
	jacp_[6,11] <- 0
	jacp_[7,11] <- ((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)
	jacp_[8,11] <- 0
	jacp_[9,11] <- 0
	jacp_[10,11] <- 0
	jacp_[11,11] <- 0
# column 12 (df/dp_12)
	jacp_[1,12] <- 0
	jacp_[2,12] <- 0
	jacp_[3,12] <- 0
	jacp_[4,12] <- 0
	jacp_[5,12] <- 0
	jacp_[6,12] <- 0
	jacp_[7,12] <- -Rii_C_cAMP
	jacp_[8,12] <- 0
	jacp_[9,12] <- 0
	jacp_[10,12] <- 0
	jacp_[11,12] <- 0
# column 13 (df/dp_13)
	jacp_[1,13] <- 0
	jacp_[2,13] <- 0
	jacp_[3,13] <- 0
	jacp_[4,13] <- 0
	jacp_[5,13] <- -C*((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)
	jacp_[6,13] <- 0
	jacp_[7,13] <- 0
	jacp_[8,13] <- 0
	jacp_[9,13] <- 0
	jacp_[10,13] <- 0
	jacp_[11,13] <- 0
# column 14 (df/dp_14)
	jacp_[1,14] <- 0
	jacp_[2,14] <- 0
	jacp_[3,14] <- 0
	jacp_[4,14] <- 0
	jacp_[5,14] <- -C*Rii_cAMP
	jacp_[6,14] <- -C*Rii_cAMP
	jacp_[7,14] <- C*Rii_cAMP
	jacp_[8,14] <- 0
	jacp_[9,14] <- 0
	jacp_[10,14] <- 0
	jacp_[11,14] <- 0
# column 15 (df/dp_15)
	jacp_[1,15] <- 0
	jacp_[2,15] <- 0
	jacp_[3,15] <- 0
	jacp_[4,15] <- 0
	jacp_[5,15] <- Rii_C_cAMP
	jacp_[6,15] <- Rii_C_cAMP
	jacp_[7,15] <- -Rii_C_cAMP
	jacp_[8,15] <- 0
	jacp_[9,15] <- 0
	jacp_[10,15] <- 0
	jacp_[11,15] <- 0
# column 16 (df/dp_16)
	jacp_[1,16] <- 0
	jacp_[2,16] <- 0
	jacp_[3,16] <- 0
	jacp_[4,16] <- Rii_C_cAMP
	jacp_[5,16] <- 0
	jacp_[6,16] <- 0
	jacp_[7,16] <- -Rii_C_cAMP
	jacp_[8,16] <- 0
	jacp_[9,16] <- 0
	jacp_[10,16] <- 0
	jacp_[11,16] <- 0
# column 17 (df/dp_17)
	jacp_[1,17] <- 0
	jacp_[2,17] <- 0
	jacp_[3,17] <- 0
	jacp_[4,17] <- 0
	jacp_[5,17] <- (-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C
	jacp_[6,17] <- 0
	jacp_[7,17] <- 0
	jacp_[8,17] <- 0
	jacp_[9,17] <- 0
	jacp_[10,17] <- 0
	jacp_[11,17] <- 0
# column 18 (df/dp_18)
	jacp_[1,18] <- -RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[2,18] <- -RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[3,18] <- 0
	jacp_[4,18] <- 0
	jacp_[5,18] <- 0
	jacp_[6,18] <- RiiP_cAMP_CaN*(1-b_AKAP)
	jacp_[7,18] <- 0
	jacp_[8,18] <- RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_CaN*(1-b_AKAP)
	jacp_[9,18] <- RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP)
	jacp_[10,18] <- 0
	jacp_[11,18] <- 0
# column 19 (df/dp_19)
	jacp_[1,19] <- RiiP_CaN*(1-b_AKAP)-RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[2,19] <- RiiP_cAMP_CaN*(1-b_AKAP)-RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)
	jacp_[3,19] <- 0
	jacp_[4,19] <- 0
	jacp_[5,19] <- 0
	jacp_[6,19] <- 0
	jacp_[7,19] <- 0
	jacp_[8,19] <- RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_CaN*(1-b_AKAP)
	jacp_[9,19] <- RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+(1-b_AKAP)^2/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP)
	jacp_[10,19] <- 0
	jacp_[11,19] <- 0
# column 20 (df/dp_20)
	jacp_[1,20] <- -RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[2,20] <- -RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[3,20] <- 0
	jacp_[4,20] <- 0
	jacp_[5,20] <- 0
	jacp_[6,20] <- RiiP_cAMP_CaN*b_AKAP
	jacp_[7,20] <- 0
	jacp_[8,20] <- RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP
	jacp_[9,20] <- RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP
	jacp_[10,20] <- 0
	jacp_[11,20] <- 0
# column 21 (df/dp_21)
	jacp_[1,21] <- RiiP_CaN*b_AKAP-RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[2,21] <- RiiP_cAMP_CaN*b_AKAP-RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)
	jacp_[3,21] <- 0
	jacp_[4,21] <- 0
	jacp_[5,21] <- 0
	jacp_[6,21] <- 0
	jacp_[7,21] <- 0
	jacp_[8,21] <- RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP
	jacp_[9,21] <- RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(b_AKAP^2/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP
	jacp_[10,21] <- 0
	jacp_[11,21] <- 0
# column 22 (df/dp_22)
	jacp_[1,22] <- 0
	jacp_[2,22] <- 0
	jacp_[3,22] <- 0
	jacp_[4,22] <- 0
	jacp_[5,22] <- -((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*C
	jacp_[6,22] <- 0
	jacp_[7,22] <- 0
	jacp_[8,22] <- 0
	jacp_[9,22] <- 0
	jacp_[10,22] <- ((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*C
	jacp_[11,22] <- 0
# column 23 (df/dp_23)
	jacp_[1,23] <- 0
	jacp_[2,23] <- 0
	jacp_[3,23] <- 0
	jacp_[4,23] <- 0
	jacp_[5,23] <- AKAR4_C
	jacp_[6,23] <- 0
	jacp_[7,23] <- 0
	jacp_[8,23] <- 0
	jacp_[9,23] <- 0
	jacp_[10,23] <- -AKAR4_C
	jacp_[11,23] <- 0
# column 24 (df/dp_24)
	jacp_[1,24] <- 0
	jacp_[2,24] <- 0
	jacp_[3,24] <- 0
	jacp_[4,24] <- 0
	jacp_[5,24] <- AKAR4_C
	jacp_[6,24] <- 0
	jacp_[7,24] <- 0
	jacp_[8,24] <- 0
	jacp_[9,24] <- 0
	jacp_[10,24] <- -AKAR4_C
	jacp_[11,24] <- AKAR4_C
# column 25 (df/dp_25)
	jacp_[1,25] <- (RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[2,25] <- (RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[3,25] <- 0
	jacp_[4,25] <- 0
	jacp_[5,25] <- 0
	jacp_[6,25] <- 0
	jacp_[7,25] <- 0
	jacp_[8,25] <- -(RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[9,25] <- -(RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF^2
	jacp_[10,25] <- 0
	jacp_[11,25] <- 0
# column 26 (df/dp_26)
	jacp_[1,26] <- (RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[2,26] <- (RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[3,26] <- 0
	jacp_[4,26] <- 0
	jacp_[5,26] <- 0
	jacp_[6,26] <- 0
	jacp_[7,26] <- 0
	jacp_[8,26] <- -(RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[9,26] <- -(RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON^2
	jacp_[10,26] <- 0
	jacp_[11,26] <- 0
# column 27 (df/dp_27)
	jacp_[1,27] <- 0
	jacp_[2,27] <- 0
	jacp_[3,27] <- RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jacp_[4,27] <- -RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jacp_[5,27] <- 0
	jacp_[6,27] <- 0
	jacp_[7,27] <- 0
	jacp_[8,27] <- 0
	jacp_[9,27] <- 0
	jacp_[10,27] <- 0
	jacp_[11,27] <- 0
# column 28 (df/dp_28)
	jacp_[1,28] <- (AKAPon_3-AKAPoff_3)*RiiP_CaN-RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)
	jacp_[2,28] <- (AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)
	jacp_[3,28] <- 0
	jacp_[4,28] <- 0
	jacp_[5,28] <- 0
	jacp_[6,28] <- (AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN
	jacp_[7,28] <- 0
	jacp_[8,28] <- RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_CaN
	jacp_[9,28] <- RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN
	jacp_[10,28] <- 0
	jacp_[11,28] <- 0
# column 29 (df/dp_29)
	jacp_[1,29] <- 0
	jacp_[2,29] <- 0
	jacp_[3,29] <- 0
	jacp_[4,29] <- 0
	jacp_[5,29] <- -C*kf_C_AKAR4
	jacp_[6,29] <- 0
	jacp_[7,29] <- 0
	jacp_[8,29] <- 0
	jacp_[9,29] <- 0
	jacp_[10,29] <- C*kf_C_AKAR4
	jacp_[11,29] <- 0
# column 30 (df/dp_30)
	jacp_[1,30] <- -RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jacp_[2,30] <- -RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jacp_[3,30] <- 0
	jacp_[4,30] <- 0
	jacp_[5,30] <- 0
	jacp_[6,30] <- 0
	jacp_[7,30] <- 0
	jacp_[8,30] <- RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jacp_[9,30] <- RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)
	jacp_[10,30] <- 0
	jacp_[11,30] <- 0
# column 31 (df/dp_31)
	jacp_[1,31] <- 0
	jacp_[2,31] <- 0
	jacp_[3,31] <- kf_Rii_C__RiiP_C
	jacp_[4,31] <- 0
	jacp_[5,31] <- kb_RiixC__Rii_C
	jacp_[6,31] <- 0
	jacp_[7,31] <- (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP
	jacp_[8,31] <- 0
	jacp_[9,31] <- 0
	jacp_[10,31] <- 0
	jacp_[11,31] <- 0
# column 32 (df/dp_32)
	jacp_[1,32] <- -RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jacp_[2,32] <- RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jacp_[3,32] <- -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jacp_[4,32] <- RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jacp_[5,32] <- 0
	jacp_[6,32] <- ((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP
	jacp_[7,32] <- ((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP
	jacp_[8,32] <- 0
	jacp_[9,32] <- 0
	jacp_[10,32] <- 0
	jacp_[11,32] <- 0
# column 33 (df/dp_33)
	jacp_[1,33] <- 0
	jacp_[2,33] <- 0
	jacp_[3,33] <- 0
	jacp_[4,33] <- 0
	jacp_[5,33] <- -C*kf_RiixC__Rii_C
	jacp_[6,33] <- (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP
	jacp_[7,33] <- 0
	jacp_[8,33] <- 0
	jacp_[9,33] <- 0
	jacp_[10,33] <- 0
	jacp_[11,33] <- 0
	rownames(jacp_) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
	colnames(jacp_) <- c("kf_Rii_C__RiiP_C", "kf_RiiP_CxcAMP__RiiP_C_cAMP", "kf_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiPXcAMP__RiiP_cAMP", "kf_RiiPXcAMP__RiiP_cAMP", "kf_RiiPxC__RiiP_C", "kb_RiiPxC__RiiP_C", "kf_cAMPxRii__Rii_cAMP", "kb_cAMPxRii__Rii_cAMP", "kf_Rii_CxcAMP__Rii_C_cAMP", "kb_Rii_CxcAMP__Rii_C_cAMP", "kf_RiixC__Rii_C", "kf_Rii_cAMPxC__Rii_C_cAMP", "kb_Rii_cAMPxC__Rii_C_cAMP", "kf_Rii_C_cAMP__RiiP_C_cAMP", "kb_RiixC__Rii_C", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP", "AKAR4_ConservedConst", "CaN_ConservedConst", "Rii_C_ConservedConst", "cAMP_ConservedConst", "Rii_ConservedConst")
	return(jacp_)
}
# ode Functions F(t,y;p)
AKAP79_func<-function(t, state, parameters)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	RiiP <- state[1]
	RiiP_cAMP <- state[2]
	RiiP_C <- state[3]
	RiiP_C_cAMP <- state[4]
	C <- state[5]
	Rii_cAMP <- state[6]
	Rii_C_cAMP <- state[7]
	RiiP_CaN <- state[8]
	RiiP_cAMP_CaN <- state[9]
	AKAR4_C <- state[10]
	AKAR4p <- state[11]
	AKAR4 <- AKAR4_ConservedConst - (AKAR4_C+AKAR4p)
	CaN <- CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN)
	Rii_C <- Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	cAMP <- cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)
	Rii <- Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN <- b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP <- kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	reaction_51 <- kf_Rii_C__RiiP_C*Rii_C
	reaction_14 <- kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	reaction_12 <- kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_43 <- kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	reaction_23 <- kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_78 <- kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	reaction_56 <- kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	reaction_76 <- kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	reaction_62 <- kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	reaction_58 <- kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	reaction_44 <- kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	reaction_33 <- kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	reaction_48 <- kf_RiiP_CaN__RiixCaN*RiiP_CaN
	reaction_37 <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	func_ <- vector(mode='numeric',len=1)
	func_[1] <- (AKAR4p*5)*71.67+100 # AKAR4pOUT
	names(func_) <- c("AKAR4pOUT")
	return(func_)
}
# output function Jacobian dF(t,y;p)/dp
AKAP79_funcJac<-function(t, state, parameters)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	RiiP <- state[1]
	RiiP_cAMP <- state[2]
	RiiP_C <- state[3]
	RiiP_C_cAMP <- state[4]
	C <- state[5]
	Rii_cAMP <- state[6]
	Rii_C_cAMP <- state[7]
	RiiP_CaN <- state[8]
	RiiP_cAMP_CaN <- state[9]
	AKAR4_C <- state[10]
	AKAR4p <- state[11]
	AKAR4<-AKAR4_ConservedConst - (AKAR4_C+AKAR4p)
	CaN<-CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN)
	Rii_C<-Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	cAMP<-cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)
	Rii<-Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN<-b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN<-b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN<-b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN<-b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP<-kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	reaction_51<-kf_Rii_C__RiiP_C*Rii_C
	reaction_14<-kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	reaction_12<-kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_43<-kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	reaction_23<-kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_78<-kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	reaction_56<-kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	reaction_76<-kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	reaction_62<-kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	reaction_58<-kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	reaction_44<-kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	reaction_33<-kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	reaction_48<-kf_RiiP_CaN__RiixCaN*RiiP_CaN
	reaction_37<-kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	jac_<-matrix(NA,1,11)
# column 1 (dF/dp_1)
	jac_[1,1] <- 0
# column 2 (dF/dp_2)
	jac_[1,2] <- 0
# column 3 (dF/dp_3)
	jac_[1,3] <- 0
# column 4 (dF/dp_4)
	jac_[1,4] <- 0
# column 5 (dF/dp_5)
	jac_[1,5] <- 0
# column 6 (dF/dp_6)
	jac_[1,6] <- 0
# column 7 (dF/dp_7)
	jac_[1,7] <- 0
# column 8 (dF/dp_8)
	jac_[1,8] <- 0
# column 9 (dF/dp_9)
	jac_[1,9] <- 0
# column 10 (dF/dp_10)
	jac_[1,10] <- 0
# column 11 (dF/dp_11)
	jac_[1,11] <- 358.35
	rownames(jac_) <- c("AKAR4pOUT")
	colnames(jac_) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p", )
	return(jac_)
}
# output function parameter Jacobian dF(t,y;p)/dp
AKAP79_funcJacp<-function(t, state, parameters)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	RiiP <- state[1]
	RiiP_cAMP <- state[2]
	RiiP_C <- state[3]
	RiiP_C_cAMP <- state[4]
	C <- state[5]
	Rii_cAMP <- state[6]
	Rii_C_cAMP <- state[7]
	RiiP_CaN <- state[8]
	RiiP_cAMP_CaN <- state[9]
	AKAR4_C <- state[10]
	AKAR4p <- state[11]
	AKAR4<-AKAR4_ConservedConst - (AKAR4_C+AKAR4p)
	CaN<-CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN)
	Rii_C<-Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C)
	cAMP<-cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN)
	Rii<-Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN<-b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN<-b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN<-b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN<-b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN<-b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP<-kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	reaction_51<-kf_Rii_C__RiiP_C*Rii_C
	reaction_14<-kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	reaction_12<-kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_43<-kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	reaction_23<-kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	reaction_78<-kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	reaction_56<-kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	reaction_76<-kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	reaction_62<-kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	reaction_58<-kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	reaction_44<-kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	reaction_33<-kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	reaction_48<-kf_RiiP_CaN__RiixCaN*RiiP_CaN
	reaction_37<-kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	jacp_<-matrix(NA,1,33)
# column 1 (dF/dp_1)
	jacp_[1,1] <- 0
# column 2 (dF/dp_2)
	jacp_[1,2] <- 0
# column 3 (dF/dp_3)
	jacp_[1,3] <- 0
# column 4 (dF/dp_4)
	jacp_[1,4] <- 0
# column 5 (dF/dp_5)
	jacp_[1,5] <- 0
# column 6 (dF/dp_6)
	jacp_[1,6] <- 0
# column 7 (dF/dp_7)
	jacp_[1,7] <- 0
# column 8 (dF/dp_8)
	jacp_[1,8] <- 0
# column 9 (dF/dp_9)
	jacp_[1,9] <- 0
# column 10 (dF/dp_10)
	jacp_[1,10] <- 0
# column 11 (dF/dp_11)
	jacp_[1,11] <- 0
# column 12 (dF/dp_12)
	jacp_[1,12] <- 0
# column 13 (dF/dp_13)
	jacp_[1,13] <- 0
# column 14 (dF/dp_14)
	jacp_[1,14] <- 0
# column 15 (dF/dp_15)
	jacp_[1,15] <- 0
# column 16 (dF/dp_16)
	jacp_[1,16] <- 0
# column 17 (dF/dp_17)
	jacp_[1,17] <- 0
# column 18 (dF/dp_18)
	jacp_[1,18] <- 0
# column 19 (dF/dp_19)
	jacp_[1,19] <- 0
# column 20 (dF/dp_20)
	jacp_[1,20] <- 0
# column 21 (dF/dp_21)
	jacp_[1,21] <- 0
# column 22 (dF/dp_22)
	jacp_[1,22] <- 0
# column 23 (dF/dp_23)
	jacp_[1,23] <- 0
# column 24 (dF/dp_24)
	jacp_[1,24] <- 0
# column 25 (dF/dp_25)
	jacp_[1,25] <- 0
# column 26 (dF/dp_26)
	jacp_[1,26] <- 0
# column 27 (dF/dp_27)
	jacp_[1,27] <- 0
# column 28 (dF/dp_28)
	jacp_[1,28] <- 0
# column 29 (dF/dp_29)
	jacp_[1,29] <- 0
# column 30 (dF/dp_30)
	jacp_[1,30] <- 0
# column 31 (dF/dp_31)
	jacp_[1,31] <- 0
# column 32 (dF/dp_32)
	jacp_[1,32] <- 0
# column 33 (dF/dp_33)
	jacp_[1,33] <- 0
	rownames(jacp_) <- c("AKAR4pOUT")
	colnames(jacp_) <- c("kf_Rii_C__RiiP_C", "kf_RiiP_CxcAMP__RiiP_C_cAMP", "kf_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiPXcAMP__RiiP_cAMP", "kf_RiiPXcAMP__RiiP_cAMP", "kf_RiiPxC__RiiP_C", "kb_RiiPxC__RiiP_C", "kf_cAMPxRii__Rii_cAMP", "kb_cAMPxRii__Rii_cAMP", "kf_Rii_CxcAMP__Rii_C_cAMP", "kb_Rii_CxcAMP__Rii_C_cAMP", "kf_RiixC__Rii_C", "kf_Rii_cAMPxC__Rii_C_cAMP", "kb_Rii_cAMPxC__Rii_C_cAMP", "kf_Rii_C_cAMP__RiiP_C_cAMP", "kb_RiixC__Rii_C", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP", "AKAR4_ConservedConst", "CaN_ConservedConst", "Rii_C_ConservedConst", "cAMP_ConservedConst", "Rii_ConservedConst")
	return(jacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
AKAP79_default<-function(t=0.0)
{
	parameters <- vector(mode='numeric',len=33)
	parameters[1] <- 33
	parameters[2] <- 0.496
	parameters[3] <- 0.00545
	parameters[4] <- 0.0156
	parameters[5] <- 0.0016
	parameters[6] <- 0.015
	parameters[7] <- 0.038
	parameters[8] <- 0.0026
	parameters[9] <- 0.015
	parameters[10] <- 0.0016
	parameters[11] <- 0.496
	parameters[12] <- 1.413
	parameters[13] <- 2.1
	parameters[14] <- 0.2984
	parameters[15] <- 0.018
	parameters[16] <- 33
	parameters[17] <- 0.0003
	parameters[18] <- 2.6
	parameters[19] <- 20
	parameters[20] <- 0.45
	parameters[21] <- 2
	parameters[22] <- 0.018
	parameters[23] <- 0.106
	parameters[24] <- 10.2
	parameters[25] <- 100
	parameters[26] <- 1
	parameters[27] <- 0.7
	parameters[28] <- 0
	parameters[29] <- 0.200000
	parameters[30] <- 1.500000
	parameters[31] <- 0.630000
	parameters[32] <- 0.000000
	parameters[33] <- 6.300000
	names(parameters) <- c("kf_Rii_C__RiiP_C", "kf_RiiP_CxcAMP__RiiP_C_cAMP", "kf_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiPXcAMP__RiiP_cAMP", "kf_RiiPXcAMP__RiiP_cAMP", "kf_RiiPxC__RiiP_C", "kb_RiiPxC__RiiP_C", "kf_cAMPxRii__Rii_cAMP", "kb_cAMPxRii__Rii_cAMP", "kf_Rii_CxcAMP__Rii_C_cAMP", "kb_Rii_CxcAMP__Rii_C_cAMP", "kf_RiixC__Rii_C", "kf_Rii_cAMPxC__Rii_C_cAMP", "kb_Rii_cAMPxC__Rii_C_cAMP", "kf_Rii_C_cAMP__RiiP_C_cAMP", "kb_RiixC__Rii_C", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP", "AKAR4_ConservedConst", "CaN_ConservedConst", "Rii_C_ConservedConst", "cAMP_ConservedConst", "Rii_ConservedConst")
	return(parameters);
}
# ode initial values
AKAP79_init<-function(t=0.0, parameters=NA)
{
	kf_Rii_C__RiiP_C <- parameters[1]
	kf_RiiP_CxcAMP__RiiP_C_cAMP <- parameters[2]
	kf_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[3]
	kb_RiiP_cAMPxC__RiiP_C_cAMP <- parameters[4]
	kb_RiiPXcAMP__RiiP_cAMP <- parameters[5]
	kf_RiiPXcAMP__RiiP_cAMP <- parameters[6]
	kf_RiiPxC__RiiP_C <- parameters[7]
	kb_RiiPxC__RiiP_C <- parameters[8]
	kf_cAMPxRii__Rii_cAMP <- parameters[9]
	kb_cAMPxRii__Rii_cAMP <- parameters[10]
	kf_Rii_CxcAMP__Rii_C_cAMP <- parameters[11]
	kb_Rii_CxcAMP__Rii_C_cAMP <- parameters[12]
	kf_RiixC__Rii_C <- parameters[13]
	kf_Rii_cAMPxC__Rii_C_cAMP <- parameters[14]
	kb_Rii_cAMPxC__Rii_C_cAMP <- parameters[15]
	kf_Rii_C_cAMP__RiiP_C_cAMP <- parameters[16]
	kb_RiixC__Rii_C <- parameters[17]
	AKAPoff_1 <- parameters[18]
	AKAPoff_3 <- parameters[19]
	AKAPon_1 <- parameters[20]
	AKAPon_3 <- parameters[21]
	kf_C_AKAR4 <- parameters[22]
	kb_C_AKAR4 <- parameters[23]
	kcat_AKARp <- parameters[24]
	kmOFF <- parameters[25]
	kmON <- parameters[26]
	KD_T <- parameters[27]
	b_AKAP <- parameters[28]
	AKAR4_ConservedConst <- parameters[29]
	CaN_ConservedConst <- parameters[30]
	Rii_C_ConservedConst <- parameters[31]
	cAMP_ConservedConst <- parameters[32]
	Rii_ConservedConst <- parameters[33]
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=11)
	state[1] <- 0
	state[2] <- 0
	state[3] <- 0
	state[4] <- 0
	state[5] <- 0
	state[6] <- 0
	state[7] <- 0
	state[8] <- 0
	state[9] <- 0
	state[10] <- 0
	state[11] <- 0
	names(state) <- c("RiiP", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4_C", "AKAR4p")
	return(state)
}
model<-list(vf=AKAP79_vf, jac=AKAP79_jac, jacp=AKAP79_jacp, func=AKAP79_func, funcJac=AKAP79_funcJac, funcJacp=AKAP79_funcJacp, init=AKAP79_init, par=AKAP79_default, name="AKAP79")
