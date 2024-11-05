# require("deSolve")

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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
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
	f_<-vector(mode='numeric',len=16)
	f_[1] <- Rii	-reaction_78-reaction_58+reaction_48
	f_[2] <- cAMP	-reaction_12-reaction_43-reaction_78-reaction_56
	f_[3] <- RiiP	-reaction_14-reaction_43-reaction_44
	f_[4] <- Rii_C	-reaction_51-reaction_56+reaction_58
	f_[5] <- RiiP_cAMP	+reaction_43-reaction_23-reaction_33
	f_[6] <- RiiP_C	+reaction_51+reaction_14-reaction_12
	f_[7] <- RiiP_C_cAMP	+reaction_12+reaction_23+reaction_62
	f_[8] <- C	-reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2
	f_[9] <- Rii_cAMP	+reaction_78-reaction_76+reaction_37
	f_[10] <- Rii_C_cAMP	+reaction_56+reaction_76-reaction_62
	f_[11] <- CaN	-reaction_44-reaction_33+reaction_48+reaction_37
	f_[12] <- RiiP_CaN	+reaction_44-reaction_48
	f_[13] <- RiiP_cAMP_CaN	+reaction_33-reaction_37
	f_[14] <- AKAR4	-reaction_1
	f_[15] <- AKAR4_C	+reaction_1-reaction_2
	f_[16] <- AKAR4p	+reaction_2
	names(f_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
## for some weird reason deSolve wants this to be a list:
	return(list(f_))
}
AKAP79_netflux <- function(t, state, parameters){
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	netFlux <- numeric(16)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN <- b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP <- kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	netFlux[1] <- kf_Rii_C__RiiP_C*Rii_C
	netFlux[2] <- kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	netFlux[3] <- kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	netFlux[4] <- kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	netFlux[5] <- kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	netFlux[6] <- kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	netFlux[7] <- kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	netFlux[8] <- kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	netFlux[9] <- kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	netFlux[10] <- kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	netFlux[11] <- kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	netFlux[12] <- kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	netFlux[13] <- kf_RiiP_CaN__RiixCaN*RiiP_CaN
	netFlux[14] <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	netFlux[15] <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	netFlux[16] <- kcat_AKARp*AKAR4_C
	return(netFlux)
}

AKAP79_fwdflux <- function(t, state, parameters){
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	fFlux <- numeric(16)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN <- b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP <- kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	# kf_Rii_C__RiiP_C*Rii_C
	fFlux[1] <- kf_Rii_C__RiiP_C*Rii_C
	# kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	fFlux[2] <- kf_RiiPxC__RiiP_C*RiiP*C  - 0
	# kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	fFlux[3] <- kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP  - 0
	# kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	fFlux[4] <- kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP  - 0
	# kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	fFlux[5] <- kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C  - 0
	# kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	fFlux[6] <- kf_cAMPxRii__Rii_cAMP*cAMP*Rii  - 0
	# kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	fFlux[7] <- kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP  - 0
	# kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	fFlux[8] <- kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C  - 0
	# kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	fFlux[9] <- kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	# kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	fFlux[10] <- kf_RiixC__Rii_C*Rii*C  - 0
	# kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	fFlux[11] <- kf_RiiPxCaN__RiiP_CaN*RiiP*CaN  - 0
	# kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	fFlux[12] <- kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP  - 0
	# kf_RiiP_CaN__RiixCaN*RiiP_CaN
	fFlux[13] <- kf_RiiP_CaN__RiixCaN*RiiP_CaN
	# kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	fFlux[14] <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	# kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	fFlux[15] <- kf_C_AKAR4*C*AKAR4  - 0
	# kcat_AKARp*AKAR4_C
	fFlux[16] <- kcat_AKARp*AKAR4_C
	return(fFlux)
}

AKAP79_bwdflux <- function(t, state, parameters){
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
	bFlux <- numeric(16)
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_RiiPxCaN__RiiP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiP_CaN__RiixCaN <- b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3
	kf_RiiPxCaN__RiiP_CaN <- b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN <- b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF
	kb_RiiP_CxcAMP__RiiP_C_cAMP <- kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T
	bFlux[1] <- 0# kf_Rii_C__RiiP_C*Rii_C
	bFlux[2] <-  kb_RiiPxC__RiiP_C*RiiP_C# kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	bFlux[3] <-  kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP# kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	bFlux[4] <-  kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP# kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	bFlux[5] <-  kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP# kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	bFlux[6] <-  kb_cAMPxRii__Rii_cAMP*Rii_cAMP# kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	bFlux[7] <-  kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP# kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	bFlux[8] <-  kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP# kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	bFlux[9] <- 0# kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	bFlux[10] <-  kb_RiixC__Rii_C*Rii_C# kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	bFlux[11] <-  kb_RiiPxCaN__RiiP_CaN*RiiP_CaN# kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	bFlux[12] <-  kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN# kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	bFlux[13] <- 0# kf_RiiP_CaN__RiixCaN*RiiP_CaN
	bFlux[14] <- 0# kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	bFlux[15] <-  kb_C_AKAR4*AKAR4_C# kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	bFlux[16] <- 0# kcat_AKARp*AKAR4_C
	return(bFlux)
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
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
	jac_ <- matrix(NA,16,16)
# column 1
	jac_[1,1] <- (-cAMP*kf_cAMPxRii__Rii_cAMP)-C*kf_RiixC__Rii_C
	jac_[2,1] <- -cAMP*kf_cAMPxRii__Rii_cAMP
	jac_[3,1] <- 0
	jac_[4,1] <- C*kf_RiixC__Rii_C
	jac_[5,1] <- 0
	jac_[6,1] <- 0
	jac_[7,1] <- 0
	jac_[8,1] <- -C*kf_RiixC__Rii_C
	jac_[9,1] <- cAMP*kf_cAMPxRii__Rii_cAMP
	jac_[10,1] <- 0
	jac_[11,1] <- 0
	jac_[12,1] <- 0
	jac_[13,1] <- 0
	jac_[14,1] <- 0
	jac_[15,1] <- 0
	jac_[16,1] <- 0
# column 2
	jac_[1,2] <- -Rii*kf_cAMPxRii__Rii_cAMP
	jac_[2,2] <- (-Rii*kf_cAMPxRii__Rii_cAMP)-Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[3,2] <- -RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[4,2] <- -Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[5,2] <- RiiP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[6,2] <- -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[7,2] <- RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[8,2] <- 0
	jac_[9,2] <- Rii*kf_cAMPxRii__Rii_cAMP
	jac_[10,2] <- Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[11,2] <- 0
	jac_[12,2] <- 0
	jac_[13,2] <- 0
	jac_[14,2] <- 0
	jac_[15,2] <- 0
	jac_[16,2] <- 0
# column 3
	jac_[1,3] <- 0
	jac_[2,3] <- -cAMP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[3,3] <- (-CaN*kf_RiiPxCaN__RiiP_CaN)-C*kf_RiiPxC__RiiP_C-cAMP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[4,3] <- 0
	jac_[5,3] <- cAMP*kf_RiiPXcAMP__RiiP_cAMP
	jac_[6,3] <- C*kf_RiiPxC__RiiP_C
	jac_[7,3] <- 0
	jac_[8,3] <- -C*kf_RiiPxC__RiiP_C
	jac_[9,3] <- 0
	jac_[10,3] <- 0
	jac_[11,3] <- -CaN*kf_RiiPxCaN__RiiP_CaN
	jac_[12,3] <- CaN*kf_RiiPxCaN__RiiP_CaN
	jac_[13,3] <- 0
	jac_[14,3] <- 0
	jac_[15,3] <- 0
	jac_[16,3] <- 0
# column 4
	jac_[1,4] <- kb_RiixC__Rii_C
	jac_[2,4] <- -cAMP*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[3,4] <- 0
	jac_[4,4] <- (-cAMP*kf_Rii_CxcAMP__Rii_C_cAMP)-kf_Rii_C__RiiP_C-kb_RiixC__Rii_C
	jac_[5,4] <- 0
	jac_[6,4] <- kf_Rii_C__RiiP_C
	jac_[7,4] <- 0
	jac_[8,4] <- kb_RiixC__Rii_C
	jac_[9,4] <- 0
	jac_[10,4] <- cAMP*kf_Rii_CxcAMP__Rii_C_cAMP
	jac_[11,4] <- 0
	jac_[12,4] <- 0
	jac_[13,4] <- 0
	jac_[14,4] <- 0
	jac_[15,4] <- 0
	jac_[16,4] <- 0
# column 5
	jac_[1,5] <- 0
	jac_[2,5] <- kb_RiiPXcAMP__RiiP_cAMP
	jac_[3,5] <- kb_RiiPXcAMP__RiiP_cAMP
	jac_[4,5] <- 0
	jac_[5,5] <- (-C*kf_RiiP_cAMPxC__RiiP_C_cAMP)-CaN*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN-kb_RiiPXcAMP__RiiP_cAMP
	jac_[6,5] <- 0
	jac_[7,5] <- C*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[8,5] <- -C*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[9,5] <- 0
	jac_[10,5] <- 0
	jac_[11,5] <- -CaN*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[12,5] <- 0
	jac_[13,5] <- CaN*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[14,5] <- 0
	jac_[15,5] <- 0
	jac_[16,5] <- 0
# column 6
	jac_[1,6] <- 0
	jac_[2,6] <- -cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[3,6] <- kb_RiiPxC__RiiP_C
	jac_[4,6] <- 0
	jac_[5,6] <- 0
	jac_[6,6] <- (-cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP)-kb_RiiPxC__RiiP_C
	jac_[7,6] <- cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[8,6] <- kb_RiiPxC__RiiP_C
	jac_[9,6] <- 0
	jac_[10,6] <- 0
	jac_[11,6] <- 0
	jac_[12,6] <- 0
	jac_[13,6] <- 0
	jac_[14,6] <- 0
	jac_[15,6] <- 0
	jac_[16,6] <- 0
# column 7
	jac_[1,7] <- 0
	jac_[2,7] <- kb_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[3,7] <- 0
	jac_[4,7] <- 0
	jac_[5,7] <- kb_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[6,7] <- kb_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[7,7] <- (-kb_RiiP_cAMPxC__RiiP_C_cAMP)-kb_RiiP_CxcAMP__RiiP_C_cAMP
	jac_[8,7] <- kb_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[9,7] <- 0
	jac_[10,7] <- 0
	jac_[11,7] <- 0
	jac_[12,7] <- 0
	jac_[13,7] <- 0
	jac_[14,7] <- 0
	jac_[15,7] <- 0
	jac_[16,7] <- 0
# column 8
	jac_[1,8] <- -Rii*kf_RiixC__Rii_C
	jac_[2,8] <- 0
	jac_[3,8] <- -RiiP*kf_RiiPxC__RiiP_C
	jac_[4,8] <- Rii*kf_RiixC__Rii_C
	jac_[5,8] <- -RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[6,8] <- RiiP*kf_RiiPxC__RiiP_C
	jac_[7,8] <- RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP
	jac_[8,8] <- (-Rii*kf_RiixC__Rii_C)-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-RiiP*kf_RiiPxC__RiiP_C-RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP-AKAR4*kf_C_AKAR4
	jac_[9,8] <- -Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP
	jac_[10,8] <- Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP
	jac_[11,8] <- 0
	jac_[12,8] <- 0
	jac_[13,8] <- 0
	jac_[14,8] <- -AKAR4*kf_C_AKAR4
	jac_[15,8] <- AKAR4*kf_C_AKAR4
	jac_[16,8] <- 0
# column 9
	jac_[1,9] <- kb_cAMPxRii__Rii_cAMP
	jac_[2,9] <- kb_cAMPxRii__Rii_cAMP
	jac_[3,9] <- 0
	jac_[4,9] <- 0
	jac_[5,9] <- 0
	jac_[6,9] <- 0
	jac_[7,9] <- 0
	jac_[8,9] <- -C*kf_Rii_cAMPxC__Rii_C_cAMP
	jac_[9,9] <- (-C*kf_Rii_cAMPxC__Rii_C_cAMP)-kb_cAMPxRii__Rii_cAMP
	jac_[10,9] <- C*kf_Rii_cAMPxC__Rii_C_cAMP
	jac_[11,9] <- 0
	jac_[12,9] <- 0
	jac_[13,9] <- 0
	jac_[14,9] <- 0
	jac_[15,9] <- 0
	jac_[16,9] <- 0
# column 10
	jac_[1,10] <- 0
	jac_[2,10] <- kb_Rii_CxcAMP__Rii_C_cAMP
	jac_[3,10] <- 0
	jac_[4,10] <- kb_Rii_CxcAMP__Rii_C_cAMP
	jac_[5,10] <- 0
	jac_[6,10] <- 0
	jac_[7,10] <- kf_Rii_C_cAMP__RiiP_C_cAMP
	jac_[8,10] <- kb_Rii_cAMPxC__Rii_C_cAMP
	jac_[9,10] <- kb_Rii_cAMPxC__Rii_C_cAMP
	jac_[10,10] <- (-kf_Rii_C_cAMP__RiiP_C_cAMP)-kb_Rii_cAMPxC__Rii_C_cAMP-kb_Rii_CxcAMP__Rii_C_cAMP
	jac_[11,10] <- 0
	jac_[12,10] <- 0
	jac_[13,10] <- 0
	jac_[14,10] <- 0
	jac_[15,10] <- 0
	jac_[16,10] <- 0
# column 11
	jac_[1,11] <- 0
	jac_[2,11] <- 0
	jac_[3,11] <- -RiiP*kf_RiiPxCaN__RiiP_CaN
	jac_[4,11] <- 0
	jac_[5,11] <- -RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[6,11] <- 0
	jac_[7,11] <- 0
	jac_[8,11] <- 0
	jac_[9,11] <- 0
	jac_[10,11] <- 0
	jac_[11,11] <- (-RiiP*kf_RiiPxCaN__RiiP_CaN)-RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[12,11] <- RiiP*kf_RiiPxCaN__RiiP_CaN
	jac_[13,11] <- RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[14,11] <- 0
	jac_[15,11] <- 0
	jac_[16,11] <- 0
# column 12
	jac_[1,12] <- kf_RiiP_CaN__RiixCaN
	jac_[2,12] <- 0
	jac_[3,12] <- kb_RiiPxCaN__RiiP_CaN
	jac_[4,12] <- 0
	jac_[5,12] <- 0
	jac_[6,12] <- 0
	jac_[7,12] <- 0
	jac_[8,12] <- 0
	jac_[9,12] <- 0
	jac_[10,12] <- 0
	jac_[11,12] <- kf_RiiP_CaN__RiixCaN+kb_RiiPxCaN__RiiP_CaN
	jac_[12,12] <- (-kf_RiiP_CaN__RiixCaN)-kb_RiiPxCaN__RiiP_CaN
	jac_[13,12] <- 0
	jac_[14,12] <- 0
	jac_[15,12] <- 0
	jac_[16,12] <- 0
# column 13
	jac_[1,13] <- 0
	jac_[2,13] <- 0
	jac_[3,13] <- 0
	jac_[4,13] <- 0
	jac_[5,13] <- kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[6,13] <- 0
	jac_[7,13] <- 0
	jac_[8,13] <- 0
	jac_[9,13] <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP
	jac_[10,13] <- 0
	jac_[11,13] <- kf_RiiP_cAMP_CaN__CaNXRii_cAMP+kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[12,13] <- 0
	jac_[13,13] <- (-kf_RiiP_cAMP_CaN__CaNXRii_cAMP)-kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	jac_[14,13] <- 0
	jac_[15,13] <- 0
	jac_[16,13] <- 0
# column 14
	jac_[1,14] <- 0
	jac_[2,14] <- 0
	jac_[3,14] <- 0
	jac_[4,14] <- 0
	jac_[5,14] <- 0
	jac_[6,14] <- 0
	jac_[7,14] <- 0
	jac_[8,14] <- -C*kf_C_AKAR4
	jac_[9,14] <- 0
	jac_[10,14] <- 0
	jac_[11,14] <- 0
	jac_[12,14] <- 0
	jac_[13,14] <- 0
	jac_[14,14] <- -C*kf_C_AKAR4
	jac_[15,14] <- C*kf_C_AKAR4
	jac_[16,14] <- 0
# column 15
	jac_[1,15] <- 0
	jac_[2,15] <- 0
	jac_[3,15] <- 0
	jac_[4,15] <- 0
	jac_[5,15] <- 0
	jac_[6,15] <- 0
	jac_[7,15] <- 0
	jac_[8,15] <- kcat_AKARp+kb_C_AKAR4
	jac_[9,15] <- 0
	jac_[10,15] <- 0
	jac_[11,15] <- 0
	jac_[12,15] <- 0
	jac_[13,15] <- 0
	jac_[14,15] <- kb_C_AKAR4
	jac_[15,15] <- (-kcat_AKARp)-kb_C_AKAR4
	jac_[16,15] <- kcat_AKARp
# column 16
	jac_[1,16] <- 0
	jac_[2,16] <- 0
	jac_[3,16] <- 0
	jac_[4,16] <- 0
	jac_[5,16] <- 0
	jac_[6,16] <- 0
	jac_[7,16] <- 0
	jac_[8,16] <- 0
	jac_[9,16] <- 0
	jac_[10,16] <- 0
	jac_[11,16] <- 0
	jac_[12,16] <- 0
	jac_[13,16] <- 0
	jac_[14,16] <- 0
	jac_[15,16] <- 0
	jac_[16,16] <- 0
	rownames(jac_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	colnames(jac_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
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
	jacp_ <- matrix(NA,16,28)
# column 1
	jacp_[1,1] <- 0
	jacp_[2,1] <- 0
	jacp_[3,1] <- 0
	jacp_[4,1] <- -Rii_C
	jacp_[5,1] <- 0
	jacp_[6,1] <- Rii_C
	jacp_[7,1] <- 0
	jacp_[8,1] <- 0
	jacp_[9,1] <- 0
	jacp_[10,1] <- 0
	jacp_[11,1] <- 0
	jacp_[12,1] <- 0
	jacp_[13,1] <- 0
	jacp_[14,1] <- 0
	jacp_[15,1] <- 0
	jacp_[16,1] <- 0
# column 2
	jacp_[1,2] <- 0
	jacp_[2,2] <- -RiiP_C*cAMP
	jacp_[3,2] <- 0
	jacp_[4,2] <- 0
	jacp_[5,2] <- 0
	jacp_[6,2] <- -RiiP_C*cAMP
	jacp_[7,2] <- RiiP_C*cAMP
	jacp_[8,2] <- 0
	jacp_[9,2] <- 0
	jacp_[10,2] <- 0
	jacp_[11,2] <- 0
	jacp_[12,2] <- 0
	jacp_[13,2] <- 0
	jacp_[14,2] <- 0
	jacp_[15,2] <- 0
	jacp_[16,2] <- 0
# column 3
	jacp_[1,3] <- 0
	jacp_[2,3] <- 0
	jacp_[3,3] <- 0
	jacp_[4,3] <- 0
	jacp_[5,3] <- -C*RiiP_cAMP
	jacp_[6,3] <- 0
	jacp_[7,3] <- C*RiiP_cAMP
	jacp_[8,3] <- -C*RiiP_cAMP
	jacp_[9,3] <- 0
	jacp_[10,3] <- 0
	jacp_[11,3] <- 0
	jacp_[12,3] <- 0
	jacp_[13,3] <- 0
	jacp_[14,3] <- 0
	jacp_[15,3] <- 0
	jacp_[16,3] <- 0
# column 4
	jacp_[1,4] <- 0
	jacp_[2,4] <- 0
	jacp_[3,4] <- 0
	jacp_[4,4] <- 0
	jacp_[5,4] <- RiiP_C_cAMP
	jacp_[6,4] <- 0
	jacp_[7,4] <- -RiiP_C_cAMP
	jacp_[8,4] <- RiiP_C_cAMP
	jacp_[9,4] <- 0
	jacp_[10,4] <- 0
	jacp_[11,4] <- 0
	jacp_[12,4] <- 0
	jacp_[13,4] <- 0
	jacp_[14,4] <- 0
	jacp_[15,4] <- 0
	jacp_[16,4] <- 0
# column 5
	jacp_[1,5] <- 0
	jacp_[2,5] <- RiiP_cAMP
	jacp_[3,5] <- RiiP_cAMP
	jacp_[4,5] <- 0
	jacp_[5,5] <- -RiiP_cAMP
	jacp_[6,5] <- 0
	jacp_[7,5] <- 0
	jacp_[8,5] <- 0
	jacp_[9,5] <- 0
	jacp_[10,5] <- 0
	jacp_[11,5] <- 0
	jacp_[12,5] <- 0
	jacp_[13,5] <- 0
	jacp_[14,5] <- 0
	jacp_[15,5] <- 0
	jacp_[16,5] <- 0
# column 6
	jacp_[1,6] <- 0
	jacp_[2,6] <- -RiiP*cAMP
	jacp_[3,6] <- -RiiP*cAMP
	jacp_[4,6] <- 0
	jacp_[5,6] <- RiiP*cAMP
	jacp_[6,6] <- 0
	jacp_[7,6] <- 0
	jacp_[8,6] <- 0
	jacp_[9,6] <- 0
	jacp_[10,6] <- 0
	jacp_[11,6] <- 0
	jacp_[12,6] <- 0
	jacp_[13,6] <- 0
	jacp_[14,6] <- 0
	jacp_[15,6] <- 0
	jacp_[16,6] <- 0
# column 7
	jacp_[1,7] <- 0
	jacp_[2,7] <- 0
	jacp_[3,7] <- -C*RiiP
	jacp_[4,7] <- 0
	jacp_[5,7] <- 0
	jacp_[6,7] <- C*RiiP
	jacp_[7,7] <- 0
	jacp_[8,7] <- -C*RiiP
	jacp_[9,7] <- 0
	jacp_[10,7] <- 0
	jacp_[11,7] <- 0
	jacp_[12,7] <- 0
	jacp_[13,7] <- 0
	jacp_[14,7] <- 0
	jacp_[15,7] <- 0
	jacp_[16,7] <- 0
# column 8
	jacp_[1,8] <- 0
	jacp_[2,8] <- 0
	jacp_[3,8] <- RiiP_C
	jacp_[4,8] <- 0
	jacp_[5,8] <- 0
	jacp_[6,8] <- -RiiP_C
	jacp_[7,8] <- 0
	jacp_[8,8] <- RiiP_C
	jacp_[9,8] <- 0
	jacp_[10,8] <- 0
	jacp_[11,8] <- 0
	jacp_[12,8] <- 0
	jacp_[13,8] <- 0
	jacp_[14,8] <- 0
	jacp_[15,8] <- 0
	jacp_[16,8] <- 0
# column 9
	jacp_[1,9] <- -Rii*cAMP
	jacp_[2,9] <- -Rii*cAMP
	jacp_[3,9] <- 0
	jacp_[4,9] <- 0
	jacp_[5,9] <- 0
	jacp_[6,9] <- 0
	jacp_[7,9] <- 0
	jacp_[8,9] <- 0
	jacp_[9,9] <- Rii*cAMP
	jacp_[10,9] <- 0
	jacp_[11,9] <- 0
	jacp_[12,9] <- 0
	jacp_[13,9] <- 0
	jacp_[14,9] <- 0
	jacp_[15,9] <- 0
	jacp_[16,9] <- 0
# column 10
	jacp_[1,10] <- Rii_cAMP
	jacp_[2,10] <- Rii_cAMP
	jacp_[3,10] <- 0
	jacp_[4,10] <- 0
	jacp_[5,10] <- 0
	jacp_[6,10] <- 0
	jacp_[7,10] <- 0
	jacp_[8,10] <- 0
	jacp_[9,10] <- -Rii_cAMP
	jacp_[10,10] <- 0
	jacp_[11,10] <- 0
	jacp_[12,10] <- 0
	jacp_[13,10] <- 0
	jacp_[14,10] <- 0
	jacp_[15,10] <- 0
	jacp_[16,10] <- 0
# column 11
	jacp_[1,11] <- 0
	jacp_[2,11] <- -Rii_C*cAMP
	jacp_[3,11] <- 0
	jacp_[4,11] <- -Rii_C*cAMP
	jacp_[5,11] <- 0
	jacp_[6,11] <- 0
	jacp_[7,11] <- 0
	jacp_[8,11] <- 0
	jacp_[9,11] <- 0
	jacp_[10,11] <- Rii_C*cAMP
	jacp_[11,11] <- 0
	jacp_[12,11] <- 0
	jacp_[13,11] <- 0
	jacp_[14,11] <- 0
	jacp_[15,11] <- 0
	jacp_[16,11] <- 0
# column 12
	jacp_[1,12] <- 0
	jacp_[2,12] <- Rii_C_cAMP
	jacp_[3,12] <- 0
	jacp_[4,12] <- Rii_C_cAMP
	jacp_[5,12] <- 0
	jacp_[6,12] <- 0
	jacp_[7,12] <- 0
	jacp_[8,12] <- 0
	jacp_[9,12] <- 0
	jacp_[10,12] <- -Rii_C_cAMP
	jacp_[11,12] <- 0
	jacp_[12,12] <- 0
	jacp_[13,12] <- 0
	jacp_[14,12] <- 0
	jacp_[15,12] <- 0
	jacp_[16,12] <- 0
# column 13
	jacp_[1,13] <- -C*Rii
	jacp_[2,13] <- 0
	jacp_[3,13] <- 0
	jacp_[4,13] <- C*Rii
	jacp_[5,13] <- 0
	jacp_[6,13] <- 0
	jacp_[7,13] <- 0
	jacp_[8,13] <- -C*Rii
	jacp_[9,13] <- 0
	jacp_[10,13] <- 0
	jacp_[11,13] <- 0
	jacp_[12,13] <- 0
	jacp_[13,13] <- 0
	jacp_[14,13] <- 0
	jacp_[15,13] <- 0
	jacp_[16,13] <- 0
# column 14
	jacp_[1,14] <- 0
	jacp_[2,14] <- 0
	jacp_[3,14] <- 0
	jacp_[4,14] <- 0
	jacp_[5,14] <- 0
	jacp_[6,14] <- 0
	jacp_[7,14] <- 0
	jacp_[8,14] <- -C*Rii_cAMP
	jacp_[9,14] <- -C*Rii_cAMP
	jacp_[10,14] <- C*Rii_cAMP
	jacp_[11,14] <- 0
	jacp_[12,14] <- 0
	jacp_[13,14] <- 0
	jacp_[14,14] <- 0
	jacp_[15,14] <- 0
	jacp_[16,14] <- 0
# column 15
	jacp_[1,15] <- 0
	jacp_[2,15] <- 0
	jacp_[3,15] <- 0
	jacp_[4,15] <- 0
	jacp_[5,15] <- 0
	jacp_[6,15] <- 0
	jacp_[7,15] <- 0
	jacp_[8,15] <- Rii_C_cAMP
	jacp_[9,15] <- Rii_C_cAMP
	jacp_[10,15] <- -Rii_C_cAMP
	jacp_[11,15] <- 0
	jacp_[12,15] <- 0
	jacp_[13,15] <- 0
	jacp_[14,15] <- 0
	jacp_[15,15] <- 0
	jacp_[16,15] <- 0
# column 16
	jacp_[1,16] <- 0
	jacp_[2,16] <- 0
	jacp_[3,16] <- 0
	jacp_[4,16] <- 0
	jacp_[5,16] <- 0
	jacp_[6,16] <- 0
	jacp_[7,16] <- Rii_C_cAMP
	jacp_[8,16] <- 0
	jacp_[9,16] <- 0
	jacp_[10,16] <- -Rii_C_cAMP
	jacp_[11,16] <- 0
	jacp_[12,16] <- 0
	jacp_[13,16] <- 0
	jacp_[14,16] <- 0
	jacp_[15,16] <- 0
	jacp_[16,16] <- 0
# column 17
	jacp_[1,17] <- Rii_C
	jacp_[2,17] <- 0
	jacp_[3,17] <- 0
	jacp_[4,17] <- -Rii_C
	jacp_[5,17] <- 0
	jacp_[6,17] <- 0
	jacp_[7,17] <- 0
	jacp_[8,17] <- Rii_C
	jacp_[9,17] <- 0
	jacp_[10,17] <- 0
	jacp_[11,17] <- 0
	jacp_[12,17] <- 0
	jacp_[13,17] <- 0
	jacp_[14,17] <- 0
	jacp_[15,17] <- 0
	jacp_[16,17] <- 0
# column 18
	jacp_[1,18] <- 0
	jacp_[2,18] <- 0
	jacp_[3,18] <- 0
	jacp_[4,18] <- 0
	jacp_[5,18] <- 0
	jacp_[6,18] <- 0
	jacp_[7,18] <- 0
	jacp_[8,18] <- 0
	jacp_[9,18] <- 0
	jacp_[10,18] <- 0
	jacp_[11,18] <- 0
	jacp_[12,18] <- 0
	jacp_[13,18] <- 0
	jacp_[14,18] <- 0
	jacp_[15,18] <- 0
	jacp_[16,18] <- 0
# column 19
	jacp_[1,19] <- 0
	jacp_[2,19] <- 0
	jacp_[3,19] <- 0
	jacp_[4,19] <- 0
	jacp_[5,19] <- 0
	jacp_[6,19] <- 0
	jacp_[7,19] <- 0
	jacp_[8,19] <- 0
	jacp_[9,19] <- 0
	jacp_[10,19] <- 0
	jacp_[11,19] <- 0
	jacp_[12,19] <- 0
	jacp_[13,19] <- 0
	jacp_[14,19] <- 0
	jacp_[15,19] <- 0
	jacp_[16,19] <- 0
# column 20
	jacp_[1,20] <- 0
	jacp_[2,20] <- 0
	jacp_[3,20] <- 0
	jacp_[4,20] <- 0
	jacp_[5,20] <- 0
	jacp_[6,20] <- 0
	jacp_[7,20] <- 0
	jacp_[8,20] <- 0
	jacp_[9,20] <- 0
	jacp_[10,20] <- 0
	jacp_[11,20] <- 0
	jacp_[12,20] <- 0
	jacp_[13,20] <- 0
	jacp_[14,20] <- 0
	jacp_[15,20] <- 0
	jacp_[16,20] <- 0
# column 21
	jacp_[1,21] <- 0
	jacp_[2,21] <- 0
	jacp_[3,21] <- 0
	jacp_[4,21] <- 0
	jacp_[5,21] <- 0
	jacp_[6,21] <- 0
	jacp_[7,21] <- 0
	jacp_[8,21] <- 0
	jacp_[9,21] <- 0
	jacp_[10,21] <- 0
	jacp_[11,21] <- 0
	jacp_[12,21] <- 0
	jacp_[13,21] <- 0
	jacp_[14,21] <- 0
	jacp_[15,21] <- 0
	jacp_[16,21] <- 0
# column 22
	jacp_[1,22] <- 0
	jacp_[2,22] <- 0
	jacp_[3,22] <- 0
	jacp_[4,22] <- 0
	jacp_[5,22] <- 0
	jacp_[6,22] <- 0
	jacp_[7,22] <- 0
	jacp_[8,22] <- -AKAR4*C
	jacp_[9,22] <- 0
	jacp_[10,22] <- 0
	jacp_[11,22] <- 0
	jacp_[12,22] <- 0
	jacp_[13,22] <- 0
	jacp_[14,22] <- -AKAR4*C
	jacp_[15,22] <- AKAR4*C
	jacp_[16,22] <- 0
# column 23
	jacp_[1,23] <- 0
	jacp_[2,23] <- 0
	jacp_[3,23] <- 0
	jacp_[4,23] <- 0
	jacp_[5,23] <- 0
	jacp_[6,23] <- 0
	jacp_[7,23] <- 0
	jacp_[8,23] <- AKAR4_C
	jacp_[9,23] <- 0
	jacp_[10,23] <- 0
	jacp_[11,23] <- 0
	jacp_[12,23] <- 0
	jacp_[13,23] <- 0
	jacp_[14,23] <- AKAR4_C
	jacp_[15,23] <- -AKAR4_C
	jacp_[16,23] <- 0
# column 24
	jacp_[1,24] <- 0
	jacp_[2,24] <- 0
	jacp_[3,24] <- 0
	jacp_[4,24] <- 0
	jacp_[5,24] <- 0
	jacp_[6,24] <- 0
	jacp_[7,24] <- 0
	jacp_[8,24] <- AKAR4_C
	jacp_[9,24] <- 0
	jacp_[10,24] <- 0
	jacp_[11,24] <- 0
	jacp_[12,24] <- 0
	jacp_[13,24] <- 0
	jacp_[14,24] <- 0
	jacp_[15,24] <- -AKAR4_C
	jacp_[16,24] <- AKAR4_C
# column 25
	jacp_[1,25] <- 0
	jacp_[2,25] <- 0
	jacp_[3,25] <- 0
	jacp_[4,25] <- 0
	jacp_[5,25] <- 0
	jacp_[6,25] <- 0
	jacp_[7,25] <- 0
	jacp_[8,25] <- 0
	jacp_[9,25] <- 0
	jacp_[10,25] <- 0
	jacp_[11,25] <- 0
	jacp_[12,25] <- 0
	jacp_[13,25] <- 0
	jacp_[14,25] <- 0
	jacp_[15,25] <- 0
	jacp_[16,25] <- 0
# column 26
	jacp_[1,26] <- 0
	jacp_[2,26] <- 0
	jacp_[3,26] <- 0
	jacp_[4,26] <- 0
	jacp_[5,26] <- 0
	jacp_[6,26] <- 0
	jacp_[7,26] <- 0
	jacp_[8,26] <- 0
	jacp_[9,26] <- 0
	jacp_[10,26] <- 0
	jacp_[11,26] <- 0
	jacp_[12,26] <- 0
	jacp_[13,26] <- 0
	jacp_[14,26] <- 0
	jacp_[15,26] <- 0
	jacp_[16,26] <- 0
# column 27
	jacp_[1,27] <- 0
	jacp_[2,27] <- 0
	jacp_[3,27] <- 0
	jacp_[4,27] <- 0
	jacp_[5,27] <- 0
	jacp_[6,27] <- 0
	jacp_[7,27] <- 0
	jacp_[8,27] <- 0
	jacp_[9,27] <- 0
	jacp_[10,27] <- 0
	jacp_[11,27] <- 0
	jacp_[12,27] <- 0
	jacp_[13,27] <- 0
	jacp_[14,27] <- 0
	jacp_[15,27] <- 0
	jacp_[16,27] <- 0
# column 28
	jacp_[1,28] <- 0
	jacp_[2,28] <- 0
	jacp_[3,28] <- 0
	jacp_[4,28] <- 0
	jacp_[5,28] <- 0
	jacp_[6,28] <- 0
	jacp_[7,28] <- 0
	jacp_[8,28] <- 0
	jacp_[9,28] <- 0
	jacp_[10,28] <- 0
	jacp_[11,28] <- 0
	jacp_[12,28] <- 0
	jacp_[13,28] <- 0
	jacp_[14,28] <- 0
	jacp_[15,28] <- 0
	jacp_[16,28] <- 0
	rownames(jacp_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	colnames(jacp_) <- c("kf_Rii_C__RiiP_C", "kf_RiiP_CxcAMP__RiiP_C_cAMP", "kf_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiPXcAMP__RiiP_cAMP", "kf_RiiPXcAMP__RiiP_cAMP", "kf_RiiPxC__RiiP_C", "kb_RiiPxC__RiiP_C", "kf_cAMPxRii__Rii_cAMP", "kb_cAMPxRii__Rii_cAMP", "kf_Rii_CxcAMP__Rii_C_cAMP", "kb_Rii_CxcAMP__Rii_C_cAMP", "kf_RiixC__Rii_C", "kf_Rii_cAMPxC__Rii_C_cAMP", "kb_Rii_cAMPxC__Rii_C_cAMP", "kf_Rii_C_cAMP__RiiP_C_cAMP", "kb_RiixC__Rii_C", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
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
	fjac_ <- matrix(NA,1,16)
# column 1
	fjac_[1,1] <- 0
# column 2
	fjac_[1,2] <- 0
# column 3
	fjac_[1,3] <- 0
# column 4
	fjac_[1,4] <- 0
# column 5
	fjac_[1,5] <- 0
# column 6
	fjac_[1,6] <- 0
# column 7
	fjac_[1,7] <- 0
# column 8
	fjac_[1,8] <- 0
# column 9
	fjac_[1,9] <- 0
# column 10
	fjac_[1,10] <- 0
# column 11
	fjac_[1,11] <- 0
# column 12
	fjac_[1,12] <- 0
# column 13
	fjac_[1,13] <- 0
# column 14
	fjac_[1,14] <- 0
# column 15
	fjac_[1,15] <- 0
# column 16
	fjac_[1,16] <- 358.35
	rownames(fjac_) <- c("AKAR4pOUT")
	colnames(fjac_) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	return(fjac_)
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
	Rii <- state[1]
	cAMP <- state[2]
	RiiP <- state[3]
	Rii_C <- state[4]
	RiiP_cAMP <- state[5]
	RiiP_C <- state[6]
	RiiP_C_cAMP <- state[7]
	C <- state[8]
	Rii_cAMP <- state[9]
	Rii_C_cAMP <- state[10]
	CaN <- state[11]
	RiiP_CaN <- state[12]
	RiiP_cAMP_CaN <- state[13]
	AKAR4 <- state[14]
	AKAR4_C <- state[15]
	AKAR4p <- state[16]
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
	fjacp_ <- matrix(NA,1,28)
# column 1
	fjacp_[1,1] <- 0
# column 2
	fjacp_[1,2] <- 0
# column 3
	fjacp_[1,3] <- 0
# column 4
	fjacp_[1,4] <- 0
# column 5
	fjacp_[1,5] <- 0
# column 6
	fjacp_[1,6] <- 0
# column 7
	fjacp_[1,7] <- 0
# column 8
	fjacp_[1,8] <- 0
# column 9
	fjacp_[1,9] <- 0
# column 10
	fjacp_[1,10] <- 0
# column 11
	fjacp_[1,11] <- 0
# column 12
	fjacp_[1,12] <- 0
# column 13
	fjacp_[1,13] <- 0
# column 14
	fjacp_[1,14] <- 0
# column 15
	fjacp_[1,15] <- 0
# column 16
	fjacp_[1,16] <- 0
# column 17
	fjacp_[1,17] <- 0
# column 18
	fjacp_[1,18] <- 0
# column 19
	fjacp_[1,19] <- 0
# column 20
	fjacp_[1,20] <- 0
# column 21
	fjacp_[1,21] <- 0
# column 22
	fjacp_[1,22] <- 0
# column 23
	fjacp_[1,23] <- 0
# column 24
	fjacp_[1,24] <- 0
# column 25
	fjacp_[1,25] <- 0
# column 26
	fjacp_[1,26] <- 0
# column 27
	fjacp_[1,27] <- 0
# column 28
	fjacp_[1,28] <- 0
	rownames(fjacp_) <- c("AKAR4pOUT")
	colnames(fjacp_) <- c("kf_Rii_C__RiiP_C", "kf_RiiP_CxcAMP__RiiP_C_cAMP", "kf_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiPXcAMP__RiiP_cAMP", "kf_RiiPXcAMP__RiiP_cAMP", "kf_RiiPxC__RiiP_C", "kb_RiiPxC__RiiP_C", "kf_cAMPxRii__Rii_cAMP", "kb_cAMPxRii__Rii_cAMP", "kf_Rii_CxcAMP__Rii_C_cAMP", "kb_Rii_CxcAMP__Rii_C_cAMP", "kf_RiixC__Rii_C", "kf_Rii_cAMPxC__Rii_C_cAMP", "kb_Rii_cAMPxC__Rii_C_cAMP", "kf_Rii_C_cAMP__RiiP_C_cAMP", "kb_RiixC__Rii_C", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
	return(fjacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
AKAP79_default<-function(t=0.0)
{
	parameters <- vector(mode='numeric',len=28)
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
	parameters[17] <- 3e-04
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
	names(parameters) <- c("kf_Rii_C__RiiP_C", "kf_RiiP_CxcAMP__RiiP_C_cAMP", "kf_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiP_cAMPxC__RiiP_C_cAMP", "kb_RiiPXcAMP__RiiP_cAMP", "kf_RiiPXcAMP__RiiP_cAMP", "kf_RiiPxC__RiiP_C", "kb_RiiPxC__RiiP_C", "kf_cAMPxRii__Rii_cAMP", "kb_cAMPxRii__Rii_cAMP", "kf_Rii_CxcAMP__Rii_C_cAMP", "kb_Rii_CxcAMP__Rii_C_cAMP", "kf_RiixC__Rii_C", "kf_Rii_cAMPxC__Rii_C_cAMP", "kb_Rii_cAMPxC__Rii_C_cAMP", "kf_Rii_C_cAMP__RiiP_C_cAMP", "kb_RiixC__Rii_C", "AKAPoff_1", "AKAPoff_3", "AKAPon_1", "AKAPon_3", "kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "kmOFF", "kmON", "KD_T", "b_AKAP")
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
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=16)
	state[1] <- 6.3
	state[2] <- 0
	state[3] <- 0
	state[4] <- 0.63
	state[5] <- 0
	state[6] <- 0
	state[7] <- 0
	state[8] <- 0
	state[9] <- 0
	state[10] <- 0
	state[11] <- 1.5
	state[12] <- 0
	state[13] <- 0
	state[14] <- 0.2
	state[15] <- 0
	state[16] <- 0
	names(state) <- c("Rii", "cAMP", "RiiP", "Rii_C", "RiiP_cAMP", "RiiP_C", "RiiP_C_cAMP", "C", "Rii_cAMP", "Rii_C_cAMP", "CaN", "RiiP_CaN", "RiiP_cAMP_CaN", "AKAR4", "AKAR4_C", "AKAR4p")
	return(state)
}
model<-list(vf=AKAP79_vf, jac=AKAP79_jac, jacp=AKAP79_jacp, func=AKAP79_func, funcJac=AKAP79_funcJac, funcJacp=AKAP79_funcJacp, init=AKAP79_init, par=AKAP79_default, name="AKAP79")
