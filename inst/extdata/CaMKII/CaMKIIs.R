# require("deSolve")

# ode vector field: y'=f(t,y;p)
CaMKIIs_vf <- function(t, state, parameters)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	CaM_Ca1 <- state[1]
	CaM_Ca2 <- state[2]
	CaM_Ca3 <- state[3]
	CaM_Ca4 <- state[4]
	PP2B_CaM <- state[5]
	PP2B_CaM_Ca1 <- state[6]
	PP2B_CaM_Ca2 <- state[7]
	PP2B_CaM_Ca3 <- state[8]
	PP2B_CaM_Ca4 <- state[9]
	CaMKII_CaM <- state[10]
	CaMKII_CaM_Ca1 <- state[11]
	CaMKII_CaM_Ca2 <- state[12]
	CaMKII_CaM_Ca3 <- state[13]
	CaMKII_CaM_Ca4 <- state[14]
	pCaMKII_CaM_Ca4 <- state[15]
	pCaMKIIa <- state[16]
	pCaMKII_CaM_Ca3 <- state[17]
	pCaMKII_CaM_Ca2 <- state[18]
	pCaMKII_CaM_Ca1 <- state[19]
	pCaMKII_CaM <- state[20]
	PP1__pCaMKIIa <- state[21]
	caa <- state[22]
	cab <- state[23]
	KD__CaM_Ca3__PP2B <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca
	KD__CaM_Ca2__PP2B <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca
	KD__CaM_Ca1__PP2B <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca
	KD__CaM__PP2B <- (KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca
	KD__CaM_Ca3__CaMKII <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__CaMKII <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__CaMKII <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca
	KD__CaM__CaMKII <- (KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca
	KD__CaM_Ca3__pCaMKIIa <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__pCaMKIIa <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__pCaMKIIa <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca
	KD__CaM__pCaMKIIa <- (KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca
	kr__CaM_Ca3__Ca <- KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca
	kr__CaM_Ca2__Ca <- KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca
	kr__CaM_Ca1__Ca <- KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca
	kr__CaM__Ca <- KD__CaM__Ca * kf__CaM__Ca
	kr__CaM_Ca4__PP2B <- KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B
	kr__CaM_Ca3__PP2B <- KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B
	kr__CaM_Ca2__PP2B <- KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B
	kr__CaM_Ca1__PP2B <- KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B
	kr__CaM__PP2B <- KD__CaM__PP2B * kf__CaM__PP2B
	kr__PP2B_CaM_Ca3__Ca <- KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca
	kr__PP2B_CaM_Ca2__Ca <- KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca
	kr__PP2B_CaM_Ca1__Ca <- KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca
	kr__PP2B_CaM__Ca <- KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca
	kr__CaM_Ca4__CaMKII <- KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII
	kr__CaM_Ca3__CaMKII <- KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII
	kr__CaM_Ca2__CaMKII <- KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII
	kr__CaM_Ca1__CaMKII <- KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII
	kr__CaM__CaMKII <- KD__CaM__CaMKII * kf__CaM__CaMKII
	kr__CaMKII_CaM_Ca3__Ca <- KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca
	kr__CaMKII_CaM_Ca2__Ca <- KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca
	kr__CaMKII_CaM_Ca1__Ca <- KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca
	kr__CaMKII_CaM__Ca <- KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca
	kr__pCaMKII_CaM_Ca3__Ca <- KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca
	kr__pCaMKII_CaM_Ca2__Ca <- KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca
	kr__pCaMKII_CaM_Ca1__Ca <- KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca
	kr__pCaMKII_CaM__Ca <- KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca
	kr__CaM_Ca4__pCaMKIIa <- KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa
	kr__CaM_Ca3__pCaMKIIa <- KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa
	kr__CaM_Ca2__pCaMKIIa <- KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa
	kr__CaM_Ca1__pCaMKIIa <- KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa
	kr__CaM__pCaMKIIa <- KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa
	Ca <- isOn*Ca_set+caa
	PP1 <- PP1_0-PP1__pCaMKIIa
	CaMKII <- CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1)
	CaM <- CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa)
	PP2B <- PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4)
	totalCaMKII <- (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM)
	ActiveCaMKII <- pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4
	r <- (CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII)
	pairedCaMKIIc <- a*(r*r)/(1+b*r)
	CaM_bound_Ca <- (1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4)
	PP2B_bound_Ca <- (1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4)
	CaMKII_bound_Ca <- (1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4))
	BoundCa <- (CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca)
	Total_CaM_CaX <- CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4
	Total_PP2B_CaM_CaX <- PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4
	Total_CaMKII_CaM_CaX <- CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4
	Total_pCaMKII_CaM_CaX <- pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4
	Total_CaM <- Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX
	Total_pCaMKII <- pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM
	ReactionFlux1 <- (kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1)
	ReactionFlux2 <- (kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2)
	ReactionFlux3 <- (kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3)
	ReactionFlux4 <- (kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4)
	ReactionFlux5 <- (kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM)
	ReactionFlux6 <- (kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1)
	ReactionFlux7 <- (kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2)
	ReactionFlux8 <- (kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3)
	ReactionFlux9 <- (kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4)
	ReactionFlux10 <- (kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1)
	ReactionFlux11 <- (kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2)
	ReactionFlux12 <- (kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3)
	ReactionFlux13 <- (kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4)
	ReactionFlux14 <- (kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM)
	ReactionFlux15 <- (kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1)
	ReactionFlux16 <- (kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2)
	ReactionFlux17 <- (kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3)
	ReactionFlux18 <- (kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4)
	ReactionFlux19 <- (kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1)
	ReactionFlux20 <- (kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2)
	ReactionFlux21 <- (kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3)
	ReactionFlux22 <- (kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4)
	ReactionFlux23 <- (kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)
	ReactionFlux24 <- (kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3)
	ReactionFlux25 <- (kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2)
	ReactionFlux26 <- (kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1)
	ReactionFlux27 <- (kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM)
	ReactionFlux28 <- (kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1)
	ReactionFlux29 <- (kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2)
	ReactionFlux30 <- (kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3)
	ReactionFlux31 <- (kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4)
	ReactionFlux32 <- (kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4)
	ReactionFlux33 <- (kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa)
	ReactionFlux34 <- (kcat__PP1__pCaMKIIa*PP1__pCaMKIIa)
	SpikeFlux1 <- cab
	SpikeFlux2 <- 0-(kca1*kca2*caa+(kca1+kca2)*cab)
	f_<-vector(mode='numeric',len=23)
	f_[1] <- +ReactionFlux1-ReactionFlux2-ReactionFlux6-ReactionFlux15-ReactionFlux26# CaM_Ca1
	f_[2] <- +ReactionFlux2-ReactionFlux3-ReactionFlux7-ReactionFlux16-ReactionFlux25# CaM_Ca2
	f_[3] <- +ReactionFlux3-ReactionFlux4-ReactionFlux8-ReactionFlux17-ReactionFlux24# CaM_Ca3
	f_[4] <- +ReactionFlux4-ReactionFlux9-ReactionFlux18-ReactionFlux23# CaM_Ca4
	f_[5] <- +ReactionFlux5-ReactionFlux10# PP2B_CaM
	f_[6] <- +ReactionFlux6+ReactionFlux10-ReactionFlux11# PP2B_CaM_Ca1
	f_[7] <- +ReactionFlux7+ReactionFlux11-ReactionFlux12# PP2B_CaM_Ca2
	f_[8] <- +ReactionFlux8+ReactionFlux12-ReactionFlux13# PP2B_CaM_Ca3
	f_[9] <- +ReactionFlux9+ReactionFlux13# PP2B_CaM_Ca4
	f_[10] <- +ReactionFlux14-ReactionFlux19# CaMKII_CaM
	f_[11] <- +ReactionFlux15+ReactionFlux19-ReactionFlux20# CaMKII_CaM_Ca1
	f_[12] <- +ReactionFlux16+ReactionFlux20-ReactionFlux21# CaMKII_CaM_Ca2
	f_[13] <- +ReactionFlux17+ReactionFlux21-ReactionFlux22# CaMKII_CaM_Ca3
	f_[14] <- +ReactionFlux18+ReactionFlux22-ReactionFlux32# CaMKII_CaM_Ca4
	f_[15] <- +ReactionFlux23+ReactionFlux31+ReactionFlux32# pCaMKII_CaM_Ca4
	f_[16] <- -ReactionFlux23-ReactionFlux24-ReactionFlux25-ReactionFlux26-ReactionFlux27-ReactionFlux33# pCaMKIIa
	f_[17] <- +ReactionFlux24+ReactionFlux30-ReactionFlux31# pCaMKII_CaM_Ca3
	f_[18] <- +ReactionFlux25+ReactionFlux29-ReactionFlux30# pCaMKII_CaM_Ca2
	f_[19] <- +ReactionFlux26+ReactionFlux28-ReactionFlux29# pCaMKII_CaM_Ca1
	f_[20] <- +ReactionFlux27-ReactionFlux28# pCaMKII_CaM
	f_[21] <- +ReactionFlux33-ReactionFlux34# PP1__pCaMKIIa
	f_[22] <- +SpikeFlux1# caa
	f_[23] <- +SpikeFlux2# cab
	names(f_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa", "caa", "cab")
## for some weird reason deSolve wants this to be a list:
	return(list(f_))
}
# ode Jacobian df(t,y;p)/dy
CaMKIIs_jac<-function(t, state, parameters)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	CaM_Ca1 <- state[1]
	CaM_Ca2 <- state[2]
	CaM_Ca3 <- state[3]
	CaM_Ca4 <- state[4]
	PP2B_CaM <- state[5]
	PP2B_CaM_Ca1 <- state[6]
	PP2B_CaM_Ca2 <- state[7]
	PP2B_CaM_Ca3 <- state[8]
	PP2B_CaM_Ca4 <- state[9]
	CaMKII_CaM <- state[10]
	CaMKII_CaM_Ca1 <- state[11]
	CaMKII_CaM_Ca2 <- state[12]
	CaMKII_CaM_Ca3 <- state[13]
	CaMKII_CaM_Ca4 <- state[14]
	pCaMKII_CaM_Ca4 <- state[15]
	pCaMKIIa <- state[16]
	pCaMKII_CaM_Ca3 <- state[17]
	pCaMKII_CaM_Ca2 <- state[18]
	pCaMKII_CaM_Ca1 <- state[19]
	pCaMKII_CaM <- state[20]
	PP1__pCaMKIIa <- state[21]
	caa <- state[22]
	cab <- state[23]
	KD__CaM_Ca3__PP2B <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca
	KD__CaM_Ca2__PP2B <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca
	KD__CaM_Ca1__PP2B <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca
	KD__CaM__PP2B <- (KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca
	KD__CaM_Ca3__CaMKII <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__CaMKII <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__CaMKII <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca
	KD__CaM__CaMKII <- (KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca
	KD__CaM_Ca3__pCaMKIIa <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__pCaMKIIa <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__pCaMKIIa <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca
	KD__CaM__pCaMKIIa <- (KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca
	kr__CaM_Ca3__Ca <- KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca
	kr__CaM_Ca2__Ca <- KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca
	kr__CaM_Ca1__Ca <- KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca
	kr__CaM__Ca <- KD__CaM__Ca * kf__CaM__Ca
	kr__CaM_Ca4__PP2B <- KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B
	kr__CaM_Ca3__PP2B <- KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B
	kr__CaM_Ca2__PP2B <- KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B
	kr__CaM_Ca1__PP2B <- KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B
	kr__CaM__PP2B <- KD__CaM__PP2B * kf__CaM__PP2B
	kr__PP2B_CaM_Ca3__Ca <- KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca
	kr__PP2B_CaM_Ca2__Ca <- KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca
	kr__PP2B_CaM_Ca1__Ca <- KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca
	kr__PP2B_CaM__Ca <- KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca
	kr__CaM_Ca4__CaMKII <- KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII
	kr__CaM_Ca3__CaMKII <- KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII
	kr__CaM_Ca2__CaMKII <- KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII
	kr__CaM_Ca1__CaMKII <- KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII
	kr__CaM__CaMKII <- KD__CaM__CaMKII * kf__CaM__CaMKII
	kr__CaMKII_CaM_Ca3__Ca <- KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca
	kr__CaMKII_CaM_Ca2__Ca <- KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca
	kr__CaMKII_CaM_Ca1__Ca <- KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca
	kr__CaMKII_CaM__Ca <- KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca
	kr__pCaMKII_CaM_Ca3__Ca <- KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca
	kr__pCaMKII_CaM_Ca2__Ca <- KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca
	kr__pCaMKII_CaM_Ca1__Ca <- KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca
	kr__pCaMKII_CaM__Ca <- KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca
	kr__CaM_Ca4__pCaMKIIa <- KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa
	kr__CaM_Ca3__pCaMKIIa <- KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa
	kr__CaM_Ca2__pCaMKIIa <- KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa
	kr__CaM_Ca1__pCaMKIIa <- KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa
	kr__CaM__pCaMKIIa <- KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa
	Ca <- isOn*Ca_set+caa
	PP1 <- PP1_0-PP1__pCaMKIIa
	CaMKII <- CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1)
	CaM <- CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa)
	PP2B <- PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4)
	totalCaMKII <- (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM)
	ActiveCaMKII <- pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4
	r <- (CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII)
	pairedCaMKIIc <- a*(r*r)/(1+b*r)
	CaM_bound_Ca <- (1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4)
	PP2B_bound_Ca <- (1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4)
	CaMKII_bound_Ca <- (1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4))
	BoundCa <- (CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca)
	Total_CaM_CaX <- CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4
	Total_PP2B_CaM_CaX <- PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4
	Total_CaMKII_CaM_CaX <- CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4
	Total_pCaMKII_CaM_CaX <- pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4
	Total_CaM <- Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX
	Total_pCaMKII <- pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM
	ReactionFlux1 <- (kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1)
	ReactionFlux2 <- (kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2)
	ReactionFlux3 <- (kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3)
	ReactionFlux4 <- (kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4)
	ReactionFlux5 <- (kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM)
	ReactionFlux6 <- (kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1)
	ReactionFlux7 <- (kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2)
	ReactionFlux8 <- (kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3)
	ReactionFlux9 <- (kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4)
	ReactionFlux10 <- (kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1)
	ReactionFlux11 <- (kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2)
	ReactionFlux12 <- (kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3)
	ReactionFlux13 <- (kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4)
	ReactionFlux14 <- (kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM)
	ReactionFlux15 <- (kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1)
	ReactionFlux16 <- (kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2)
	ReactionFlux17 <- (kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3)
	ReactionFlux18 <- (kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4)
	ReactionFlux19 <- (kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1)
	ReactionFlux20 <- (kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2)
	ReactionFlux21 <- (kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3)
	ReactionFlux22 <- (kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4)
	ReactionFlux23 <- (kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)
	ReactionFlux24 <- (kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3)
	ReactionFlux25 <- (kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2)
	ReactionFlux26 <- (kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1)
	ReactionFlux27 <- (kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM)
	ReactionFlux28 <- (kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1)
	ReactionFlux29 <- (kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2)
	ReactionFlux30 <- (kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3)
	ReactionFlux31 <- (kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4)
	ReactionFlux32 <- (kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4)
	ReactionFlux33 <- (kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa)
	ReactionFlux34 <- (kcat__PP1__pCaMKIIa*PP1__pCaMKIIa)
	SpikeFlux1 <- cab
	SpikeFlux2 <- 0-(kca1*kca2*caa+(kca1+kca2)*cab)
	jac_ <- matrix(0.0,23,23)
# column 1
	jac_[1,1] <- (-kf__CaM_Ca1__pCaMKIIa*pCaMKIIa)-kr__CaM__Ca-PP2B*kf__CaM_Ca1__PP2B-CaMKII*kf__CaM_Ca1__CaMKII-Ca*kf__CaM_Ca1__Ca
	jac_[2,1] <- Ca*kf__CaM_Ca1__Ca
	jac_[6,1] <- PP2B*kf__CaM_Ca1__PP2B
	jac_[11,1] <- CaMKII*kf__CaM_Ca1__CaMKII
	jac_[16,1] <- -kf__CaM_Ca1__pCaMKIIa*pCaMKIIa
	jac_[19,1] <- kf__CaM_Ca1__pCaMKIIa*pCaMKIIa
# column 2
	jac_[1,2] <- kr__CaM_Ca1__Ca
	jac_[2,2] <- (-kf__CaM_Ca2__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca1__Ca-PP2B*kf__CaM_Ca2__PP2B-CaMKII*kf__CaM_Ca2__CaMKII-Ca*kf__CaM_Ca2__Ca
	jac_[3,2] <- Ca*kf__CaM_Ca2__Ca
	jac_[7,2] <- PP2B*kf__CaM_Ca2__PP2B
	jac_[12,2] <- CaMKII*kf__CaM_Ca2__CaMKII
	jac_[16,2] <- -kf__CaM_Ca2__pCaMKIIa*pCaMKIIa
	jac_[18,2] <- kf__CaM_Ca2__pCaMKIIa*pCaMKIIa
# column 3
	jac_[2,3] <- kr__CaM_Ca2__Ca
	jac_[3,3] <- (-kf__CaM_Ca3__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca2__Ca-PP2B*kf__CaM_Ca3__PP2B-CaMKII*kf__CaM_Ca3__CaMKII-Ca*kf__CaM_Ca3__Ca
	jac_[4,3] <- Ca*kf__CaM_Ca3__Ca
	jac_[8,3] <- PP2B*kf__CaM_Ca3__PP2B
	jac_[13,3] <- CaMKII*kf__CaM_Ca3__CaMKII
	jac_[16,3] <- -kf__CaM_Ca3__pCaMKIIa*pCaMKIIa
	jac_[17,3] <- kf__CaM_Ca3__pCaMKIIa*pCaMKIIa
# column 4
	jac_[3,4] <- kr__CaM_Ca3__Ca
	jac_[4,4] <- (-kf__CaM_Ca4__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca3__Ca-PP2B*kf__CaM_Ca4__PP2B-CaMKII*kf__CaM_Ca4__CaMKII
	jac_[9,4] <- PP2B*kf__CaM_Ca4__PP2B
	jac_[14,4] <- CaMKII*kf__CaM_Ca4__CaMKII
	jac_[15,4] <- kf__CaM_Ca4__pCaMKIIa*pCaMKIIa
	jac_[16,4] <- -kf__CaM_Ca4__pCaMKIIa*pCaMKIIa
# column 5
	jac_[5,5] <- (-kr__CaM__PP2B)-Ca*kf__PP2B_CaM__Ca
	jac_[6,5] <- Ca*kf__PP2B_CaM__Ca
# column 6
	jac_[1,6] <- kr__CaM_Ca1__PP2B
	jac_[5,6] <- kr__PP2B_CaM__Ca
	jac_[6,6] <- (-kr__PP2B_CaM__Ca)-kr__CaM_Ca1__PP2B-Ca*kf__PP2B_CaM_Ca1__Ca
	jac_[7,6] <- Ca*kf__PP2B_CaM_Ca1__Ca
# column 7
	jac_[2,7] <- kr__CaM_Ca2__PP2B
	jac_[6,7] <- kr__PP2B_CaM_Ca1__Ca
	jac_[7,7] <- (-kr__PP2B_CaM_Ca1__Ca)-kr__CaM_Ca2__PP2B-Ca*kf__PP2B_CaM_Ca2__Ca
	jac_[8,7] <- Ca*kf__PP2B_CaM_Ca2__Ca
# column 8
	jac_[3,8] <- kr__CaM_Ca3__PP2B
	jac_[7,8] <- kr__PP2B_CaM_Ca2__Ca
	jac_[8,8] <- (-kr__PP2B_CaM_Ca2__Ca)-kr__CaM_Ca3__PP2B-Ca*kf__PP2B_CaM_Ca3__Ca
	jac_[9,8] <- Ca*kf__PP2B_CaM_Ca3__Ca
# column 9
	jac_[4,9] <- kr__CaM_Ca4__PP2B
	jac_[8,9] <- kr__PP2B_CaM_Ca3__Ca
	jac_[9,9] <- (-kr__PP2B_CaM_Ca3__Ca)-kr__CaM_Ca4__PP2B
# column 10
	jac_[10,10] <- (-kr__CaM__CaMKII)-Ca*kf__CaMKII_CaM__Ca
	jac_[11,10] <- Ca*kf__CaMKII_CaM__Ca
# column 11
	jac_[1,11] <- kr__CaM_Ca1__CaMKII
	jac_[10,11] <- kr__CaMKII_CaM__Ca
	jac_[11,11] <- (-kr__CaM_Ca1__CaMKII)-kr__CaMKII_CaM__Ca-Ca*kf__CaMKII_CaM_Ca1__Ca
	jac_[12,11] <- Ca*kf__CaMKII_CaM_Ca1__Ca
# column 12
	jac_[2,12] <- kr__CaM_Ca2__CaMKII
	jac_[11,12] <- kr__CaMKII_CaM_Ca1__Ca
	jac_[12,12] <- (-kr__CaM_Ca2__CaMKII)-kr__CaMKII_CaM_Ca1__Ca-Ca*kf__CaMKII_CaM_Ca2__Ca
	jac_[13,12] <- Ca*kf__CaMKII_CaM_Ca2__Ca
# column 13
	jac_[3,13] <- kr__CaM_Ca3__CaMKII
	jac_[12,13] <- kr__CaMKII_CaM_Ca2__Ca
	jac_[13,13] <- (-kr__CaM_Ca3__CaMKII)-kr__CaMKII_CaM_Ca2__Ca-Ca*kf__CaMKII_CaM_Ca3__Ca
	jac_[14,13] <- Ca*kf__CaMKII_CaM_Ca3__Ca
# column 14
	jac_[4,14] <- kr__CaM_Ca4__CaMKII
	jac_[13,14] <- kr__CaMKII_CaM_Ca3__Ca
	jac_[14,14] <- (-kautMax*pairedCaMKIIc)-kr__CaM_Ca4__CaMKII-kr__CaMKII_CaM_Ca3__Ca
	jac_[15,14] <- kautMax*pairedCaMKIIc
# column 15
	jac_[4,15] <- kr__CaM_Ca4__pCaMKIIa
	jac_[15,15] <- (-kr__pCaMKII_CaM_Ca3__Ca)-kr__CaM_Ca4__pCaMKIIa
	jac_[16,15] <- kr__CaM_Ca4__pCaMKIIa
	jac_[17,15] <- kr__pCaMKII_CaM_Ca3__Ca
# column 16
	jac_[1,16] <- -CaM_Ca1*kf__CaM_Ca1__pCaMKIIa
	jac_[2,16] <- -CaM_Ca2*kf__CaM_Ca2__pCaMKIIa
	jac_[3,16] <- -CaM_Ca3*kf__CaM_Ca3__pCaMKIIa
	jac_[4,16] <- -CaM_Ca4*kf__CaM_Ca4__pCaMKIIa
	jac_[15,16] <- CaM_Ca4*kf__CaM_Ca4__pCaMKIIa
	jac_[16,16] <- (-PP1*kf__PP1__pCaMKIIa)-CaM*kf__CaM__pCaMKIIa-CaM_Ca4*kf__CaM_Ca4__pCaMKIIa-CaM_Ca3*kf__CaM_Ca3__pCaMKIIa-CaM_Ca2*kf__CaM_Ca2__pCaMKIIa-CaM_Ca1*kf__CaM_Ca1__pCaMKIIa
	jac_[17,16] <- CaM_Ca3*kf__CaM_Ca3__pCaMKIIa
	jac_[18,16] <- CaM_Ca2*kf__CaM_Ca2__pCaMKIIa
	jac_[19,16] <- CaM_Ca1*kf__CaM_Ca1__pCaMKIIa
	jac_[20,16] <- CaM*kf__CaM__pCaMKIIa
	jac_[21,16] <- PP1*kf__PP1__pCaMKIIa
# column 17
	jac_[3,17] <- kr__CaM_Ca3__pCaMKIIa
	jac_[15,17] <- Ca*kf__pCaMKII_CaM_Ca3__Ca
	jac_[16,17] <- kr__CaM_Ca3__pCaMKIIa
	jac_[17,17] <- (-kr__pCaMKII_CaM_Ca2__Ca)-kr__CaM_Ca3__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca3__Ca
	jac_[18,17] <- kr__pCaMKII_CaM_Ca2__Ca
# column 18
	jac_[2,18] <- kr__CaM_Ca2__pCaMKIIa
	jac_[16,18] <- kr__CaM_Ca2__pCaMKIIa
	jac_[17,18] <- Ca*kf__pCaMKII_CaM_Ca2__Ca
	jac_[18,18] <- (-kr__pCaMKII_CaM_Ca1__Ca)-kr__CaM_Ca2__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca2__Ca
	jac_[19,18] <- kr__pCaMKII_CaM_Ca1__Ca
# column 19
	jac_[1,19] <- kr__CaM_Ca1__pCaMKIIa
	jac_[16,19] <- kr__CaM_Ca1__pCaMKIIa
	jac_[18,19] <- Ca*kf__pCaMKII_CaM_Ca1__Ca
	jac_[19,19] <- (-kr__pCaMKII_CaM__Ca)-kr__CaM_Ca1__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca1__Ca
	jac_[20,19] <- kr__pCaMKII_CaM__Ca
# column 20
	jac_[16,20] <- kr__CaM__pCaMKIIa
	jac_[19,20] <- Ca*kf__pCaMKII_CaM__Ca
	jac_[20,20] <- (-kr__CaM__pCaMKIIa)-Ca*kf__pCaMKII_CaM__Ca
# column 21
	jac_[16,21] <- kr__PP1__pCaMKIIa
	jac_[21,21] <- (-kr__PP1__pCaMKIIa)-kcat__PP1__pCaMKIIa
# column 22
	jac_[23,22] <- -kca1*kca2
# column 23
	jac_[22,23] <- 1
	jac_[23,23] <- (-kca2)-kca1
	rownames(jac_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa", "caa", "cab")
	colnames(jac_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa", "caa", "cab")
	return(jac_)
}
# ode parameter Jacobian df(t,y;p)/dp
CaMKIIs_jacp<-function(t, state, parameters)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	CaM_Ca1 <- state[1]
	CaM_Ca2 <- state[2]
	CaM_Ca3 <- state[3]
	CaM_Ca4 <- state[4]
	PP2B_CaM <- state[5]
	PP2B_CaM_Ca1 <- state[6]
	PP2B_CaM_Ca2 <- state[7]
	PP2B_CaM_Ca3 <- state[8]
	PP2B_CaM_Ca4 <- state[9]
	CaMKII_CaM <- state[10]
	CaMKII_CaM_Ca1 <- state[11]
	CaMKII_CaM_Ca2 <- state[12]
	CaMKII_CaM_Ca3 <- state[13]
	CaMKII_CaM_Ca4 <- state[14]
	pCaMKII_CaM_Ca4 <- state[15]
	pCaMKIIa <- state[16]
	pCaMKII_CaM_Ca3 <- state[17]
	pCaMKII_CaM_Ca2 <- state[18]
	pCaMKII_CaM_Ca1 <- state[19]
	pCaMKII_CaM <- state[20]
	PP1__pCaMKIIa <- state[21]
	caa <- state[22]
	cab <- state[23]
	KD__CaM_Ca3__PP2B<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca
	KD__CaM_Ca2__PP2B<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca
	KD__CaM_Ca1__PP2B<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca
	KD__CaM__PP2B<-(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca
	KD__CaM_Ca3__CaMKII<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__CaMKII<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__CaMKII<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca
	KD__CaM__CaMKII<-(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca
	KD__CaM_Ca3__pCaMKIIa<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__pCaMKIIa<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__pCaMKIIa<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca
	KD__CaM__pCaMKIIa<-(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca
	kr__CaM_Ca3__Ca<-KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca
	kr__CaM_Ca2__Ca<-KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca
	kr__CaM_Ca1__Ca<-KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca
	kr__CaM__Ca<-KD__CaM__Ca * kf__CaM__Ca
	kr__CaM_Ca4__PP2B<-KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B
	kr__CaM_Ca3__PP2B<-KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B
	kr__CaM_Ca2__PP2B<-KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B
	kr__CaM_Ca1__PP2B<-KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B
	kr__CaM__PP2B<-KD__CaM__PP2B * kf__CaM__PP2B
	kr__PP2B_CaM_Ca3__Ca<-KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca
	kr__PP2B_CaM_Ca2__Ca<-KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca
	kr__PP2B_CaM_Ca1__Ca<-KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca
	kr__PP2B_CaM__Ca<-KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca
	kr__CaM_Ca4__CaMKII<-KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII
	kr__CaM_Ca3__CaMKII<-KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII
	kr__CaM_Ca2__CaMKII<-KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII
	kr__CaM_Ca1__CaMKII<-KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII
	kr__CaM__CaMKII<-KD__CaM__CaMKII * kf__CaM__CaMKII
	kr__CaMKII_CaM_Ca3__Ca<-KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca
	kr__CaMKII_CaM_Ca2__Ca<-KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca
	kr__CaMKII_CaM_Ca1__Ca<-KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca
	kr__CaMKII_CaM__Ca<-KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca
	kr__pCaMKII_CaM_Ca3__Ca<-KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca
	kr__pCaMKII_CaM_Ca2__Ca<-KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca
	kr__pCaMKII_CaM_Ca1__Ca<-KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca
	kr__pCaMKII_CaM__Ca<-KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca
	kr__CaM_Ca4__pCaMKIIa<-KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa
	kr__CaM_Ca3__pCaMKIIa<-KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa
	kr__CaM_Ca2__pCaMKIIa<-KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa
	kr__CaM_Ca1__pCaMKIIa<-KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa
	kr__CaM__pCaMKIIa<-KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa
	Ca<-isOn*Ca_set+caa
	PP1<-PP1_0-PP1__pCaMKIIa
	CaMKII<-CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1)
	CaM<-CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa)
	PP2B<-PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4)
	totalCaMKII<-(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM)
	ActiveCaMKII<-pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4
	r<-(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII)
	pairedCaMKIIc<-a*(r*r)/(1+b*r)
	CaM_bound_Ca<-(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4)
	PP2B_bound_Ca<-(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4)
	CaMKII_bound_Ca<-(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4))
	BoundCa<-(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca)
	Total_CaM_CaX<-CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4
	Total_PP2B_CaM_CaX<-PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4
	Total_CaMKII_CaM_CaX<-CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4
	Total_pCaMKII_CaM_CaX<-pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4
	Total_CaM<-Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX
	Total_pCaMKII<-pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM
	ReactionFlux1<-(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1)
	ReactionFlux2<-(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2)
	ReactionFlux3<-(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3)
	ReactionFlux4<-(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4)
	ReactionFlux5<-(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM)
	ReactionFlux6<-(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1)
	ReactionFlux7<-(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2)
	ReactionFlux8<-(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3)
	ReactionFlux9<-(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4)
	ReactionFlux10<-(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1)
	ReactionFlux11<-(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2)
	ReactionFlux12<-(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3)
	ReactionFlux13<-(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4)
	ReactionFlux14<-(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM)
	ReactionFlux15<-(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1)
	ReactionFlux16<-(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2)
	ReactionFlux17<-(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3)
	ReactionFlux18<-(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4)
	ReactionFlux19<-(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1)
	ReactionFlux20<-(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2)
	ReactionFlux21<-(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3)
	ReactionFlux22<-(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4)
	ReactionFlux23<-(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)
	ReactionFlux24<-(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3)
	ReactionFlux25<-(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2)
	ReactionFlux26<-(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1)
	ReactionFlux27<-(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM)
	ReactionFlux28<-(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1)
	ReactionFlux29<-(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2)
	ReactionFlux30<-(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3)
	ReactionFlux31<-(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4)
	ReactionFlux32<-(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4)
	ReactionFlux33<-(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa)
	ReactionFlux34<-(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa)
	SpikeFlux1<-cab
	SpikeFlux2<-0-(kca1*kca2*caa+(kca1+kca2)*cab)
	jacp_ <- matrix(0.0,23,62)
# column 1
	jacp_[1,1] <- Ca*CaM
# column 2
	jacp_[1,2] <- -Ca*CaM_Ca1
	jacp_[2,2] <- Ca*CaM_Ca1
# column 3
	jacp_[2,3] <- -Ca*CaM_Ca2
	jacp_[3,3] <- Ca*CaM_Ca2
# column 4
	jacp_[3,4] <- -Ca*CaM_Ca3
	jacp_[4,4] <- Ca*CaM_Ca3
# column 5
	jacp_[5,5] <- CaM*PP2B
# column 6
	jacp_[1,6] <- -CaM_Ca1*PP2B
	jacp_[6,6] <- CaM_Ca1*PP2B
# column 7
	jacp_[2,7] <- -CaM_Ca2*PP2B
	jacp_[7,7] <- CaM_Ca2*PP2B
# column 8
	jacp_[3,8] <- -CaM_Ca3*PP2B
	jacp_[8,8] <- CaM_Ca3*PP2B
# column 9
	jacp_[4,9] <- -CaM_Ca4*PP2B
	jacp_[9,9] <- CaM_Ca4*PP2B
# column 10
	jacp_[5,10] <- -Ca*PP2B_CaM
	jacp_[6,10] <- Ca*PP2B_CaM
# column 11
	jacp_[6,11] <- -Ca*PP2B_CaM_Ca1
	jacp_[7,11] <- Ca*PP2B_CaM_Ca1
# column 12
	jacp_[7,12] <- -Ca*PP2B_CaM_Ca2
	jacp_[8,12] <- Ca*PP2B_CaM_Ca2
# column 13
	jacp_[8,13] <- -Ca*PP2B_CaM_Ca3
	jacp_[9,13] <- Ca*PP2B_CaM_Ca3
# column 14
# column 15
# column 16
# column 17
# column 18
# column 19
# column 20
# column 21
# column 22
# column 23
	jacp_[10,23] <- CaM*CaMKII
# column 24
	jacp_[13,24] <- -Ca*CaMKII_CaM_Ca3
	jacp_[14,24] <- Ca*CaMKII_CaM_Ca3
# column 25
	jacp_[12,25] <- -Ca*CaMKII_CaM_Ca2
	jacp_[13,25] <- Ca*CaMKII_CaM_Ca2
# column 26
	jacp_[11,26] <- -Ca*CaMKII_CaM_Ca1
	jacp_[12,26] <- Ca*CaMKII_CaM_Ca1
# column 27
	jacp_[10,27] <- -Ca*CaMKII_CaM
	jacp_[11,27] <- Ca*CaMKII_CaM
# column 28
	jacp_[1,28] <- -CaMKII*CaM_Ca1
	jacp_[11,28] <- CaMKII*CaM_Ca1
# column 29
	jacp_[2,29] <- -CaMKII*CaM_Ca2
	jacp_[12,29] <- CaMKII*CaM_Ca2
# column 30
	jacp_[3,30] <- -CaMKII*CaM_Ca3
	jacp_[13,30] <- CaMKII*CaM_Ca3
# column 31
	jacp_[4,31] <- -CaMKII*CaM_Ca4
	jacp_[14,31] <- CaMKII*CaM_Ca4
# column 32
# column 33
# column 34
# column 35
# column 36
# column 37
	jacp_[15,37] <- Ca*pCaMKII_CaM_Ca3
	jacp_[17,37] <- -Ca*pCaMKII_CaM_Ca3
# column 38
	jacp_[16,38] <- -CaM*pCaMKIIa
	jacp_[20,38] <- CaM*pCaMKIIa
# column 39
	jacp_[1,39] <- -CaM_Ca1*pCaMKIIa
	jacp_[16,39] <- -CaM_Ca1*pCaMKIIa
	jacp_[19,39] <- CaM_Ca1*pCaMKIIa
# column 40
	jacp_[2,40] <- -CaM_Ca2*pCaMKIIa
	jacp_[16,40] <- -CaM_Ca2*pCaMKIIa
	jacp_[18,40] <- CaM_Ca2*pCaMKIIa
# column 41
	jacp_[3,41] <- -CaM_Ca3*pCaMKIIa
	jacp_[16,41] <- -CaM_Ca3*pCaMKIIa
	jacp_[17,41] <- CaM_Ca3*pCaMKIIa
# column 42
	jacp_[17,42] <- Ca*pCaMKII_CaM_Ca2
	jacp_[18,42] <- -Ca*pCaMKII_CaM_Ca2
# column 43
	jacp_[18,43] <- Ca*pCaMKII_CaM_Ca1
	jacp_[19,43] <- -Ca*pCaMKII_CaM_Ca1
# column 44
	jacp_[4,44] <- -CaM_Ca4*pCaMKIIa
	jacp_[15,44] <- CaM_Ca4*pCaMKIIa
	jacp_[16,44] <- -CaM_Ca4*pCaMKIIa
# column 45
	jacp_[19,45] <- Ca*pCaMKII_CaM
	jacp_[20,45] <- -Ca*pCaMKII_CaM
# column 46
# column 47
# column 48
# column 49
# column 50
# column 51
	jacp_[14,51] <- -CaMKII_CaM_Ca4*pairedCaMKIIc
	jacp_[15,51] <- CaMKII_CaM_Ca4*pairedCaMKIIc
# column 52
	jacp_[16,52] <- -PP1*pCaMKIIa
	jacp_[21,52] <- PP1*pCaMKIIa
# column 53
	jacp_[16,53] <- PP1__pCaMKIIa
	jacp_[21,53] <- -PP1__pCaMKIIa
# column 54
	jacp_[21,54] <- -PP1__pCaMKIIa
# column 55
# column 56
# column 57
# column 58
# column 59
# column 60
	jacp_[23,60] <- (-caa*kca2)-cab
# column 61
	jacp_[23,61] <- (-caa*kca1)-cab
# column 62
	rownames(jacp_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa", "caa", "cab")
	colnames(jacp_) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "Ca_set", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0", "kca1", "kca2", "isOn")
	return(jacp_)
}
# ode Functions F(t,y;p)
CaMKIIs_func<-function(t, state, parameters)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	CaM_Ca1 <- state[1]
	CaM_Ca2 <- state[2]
	CaM_Ca3 <- state[3]
	CaM_Ca4 <- state[4]
	PP2B_CaM <- state[5]
	PP2B_CaM_Ca1 <- state[6]
	PP2B_CaM_Ca2 <- state[7]
	PP2B_CaM_Ca3 <- state[8]
	PP2B_CaM_Ca4 <- state[9]
	CaMKII_CaM <- state[10]
	CaMKII_CaM_Ca1 <- state[11]
	CaMKII_CaM_Ca2 <- state[12]
	CaMKII_CaM_Ca3 <- state[13]
	CaMKII_CaM_Ca4 <- state[14]
	pCaMKII_CaM_Ca4 <- state[15]
	pCaMKIIa <- state[16]
	pCaMKII_CaM_Ca3 <- state[17]
	pCaMKII_CaM_Ca2 <- state[18]
	pCaMKII_CaM_Ca1 <- state[19]
	pCaMKII_CaM <- state[20]
	PP1__pCaMKIIa <- state[21]
	caa <- state[22]
	cab <- state[23]
	KD__CaM_Ca3__PP2B <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca
	KD__CaM_Ca2__PP2B <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca
	KD__CaM_Ca1__PP2B <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca
	KD__CaM__PP2B <- (KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca
	KD__CaM_Ca3__CaMKII <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__CaMKII <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__CaMKII <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca
	KD__CaM__CaMKII <- (KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca
	KD__CaM_Ca3__pCaMKIIa <- (KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__pCaMKIIa <- (KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__pCaMKIIa <- (KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca
	KD__CaM__pCaMKIIa <- (KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca
	kr__CaM_Ca3__Ca <- KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca
	kr__CaM_Ca2__Ca <- KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca
	kr__CaM_Ca1__Ca <- KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca
	kr__CaM__Ca <- KD__CaM__Ca * kf__CaM__Ca
	kr__CaM_Ca4__PP2B <- KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B
	kr__CaM_Ca3__PP2B <- KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B
	kr__CaM_Ca2__PP2B <- KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B
	kr__CaM_Ca1__PP2B <- KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B
	kr__CaM__PP2B <- KD__CaM__PP2B * kf__CaM__PP2B
	kr__PP2B_CaM_Ca3__Ca <- KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca
	kr__PP2B_CaM_Ca2__Ca <- KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca
	kr__PP2B_CaM_Ca1__Ca <- KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca
	kr__PP2B_CaM__Ca <- KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca
	kr__CaM_Ca4__CaMKII <- KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII
	kr__CaM_Ca3__CaMKII <- KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII
	kr__CaM_Ca2__CaMKII <- KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII
	kr__CaM_Ca1__CaMKII <- KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII
	kr__CaM__CaMKII <- KD__CaM__CaMKII * kf__CaM__CaMKII
	kr__CaMKII_CaM_Ca3__Ca <- KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca
	kr__CaMKII_CaM_Ca2__Ca <- KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca
	kr__CaMKII_CaM_Ca1__Ca <- KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca
	kr__CaMKII_CaM__Ca <- KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca
	kr__pCaMKII_CaM_Ca3__Ca <- KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca
	kr__pCaMKII_CaM_Ca2__Ca <- KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca
	kr__pCaMKII_CaM_Ca1__Ca <- KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca
	kr__pCaMKII_CaM__Ca <- KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca
	kr__CaM_Ca4__pCaMKIIa <- KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa
	kr__CaM_Ca3__pCaMKIIa <- KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa
	kr__CaM_Ca2__pCaMKIIa <- KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa
	kr__CaM_Ca1__pCaMKIIa <- KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa
	kr__CaM__pCaMKIIa <- KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa
	Ca <- isOn*Ca_set+caa
	PP1 <- PP1_0-PP1__pCaMKIIa
	CaMKII <- CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1)
	CaM <- CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa)
	PP2B <- PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4)
	totalCaMKII <- (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM)
	ActiveCaMKII <- pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4
	r <- (CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII)
	pairedCaMKIIc <- a*(r*r)/(1+b*r)
	CaM_bound_Ca <- (1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4)
	PP2B_bound_Ca <- (1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4)
	CaMKII_bound_Ca <- (1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4))
	BoundCa <- (CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca)
	Total_CaM_CaX <- CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4
	Total_PP2B_CaM_CaX <- PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4
	Total_CaMKII_CaM_CaX <- CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4
	Total_pCaMKII_CaM_CaX <- pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4
	Total_CaM <- Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX
	Total_pCaMKII <- pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM
	ReactionFlux1 <- (kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1)
	ReactionFlux2 <- (kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2)
	ReactionFlux3 <- (kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3)
	ReactionFlux4 <- (kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4)
	ReactionFlux5 <- (kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM)
	ReactionFlux6 <- (kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1)
	ReactionFlux7 <- (kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2)
	ReactionFlux8 <- (kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3)
	ReactionFlux9 <- (kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4)
	ReactionFlux10 <- (kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1)
	ReactionFlux11 <- (kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2)
	ReactionFlux12 <- (kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3)
	ReactionFlux13 <- (kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4)
	ReactionFlux14 <- (kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM)
	ReactionFlux15 <- (kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1)
	ReactionFlux16 <- (kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2)
	ReactionFlux17 <- (kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3)
	ReactionFlux18 <- (kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4)
	ReactionFlux19 <- (kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1)
	ReactionFlux20 <- (kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2)
	ReactionFlux21 <- (kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3)
	ReactionFlux22 <- (kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4)
	ReactionFlux23 <- (kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)
	ReactionFlux24 <- (kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3)
	ReactionFlux25 <- (kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2)
	ReactionFlux26 <- (kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1)
	ReactionFlux27 <- (kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM)
	ReactionFlux28 <- (kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1)
	ReactionFlux29 <- (kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2)
	ReactionFlux30 <- (kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3)
	ReactionFlux31 <- (kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4)
	ReactionFlux32 <- (kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4)
	ReactionFlux33 <- (kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa)
	ReactionFlux34 <- (kcat__PP1__pCaMKIIa*PP1__pCaMKIIa)
	SpikeFlux1 <- cab
	SpikeFlux2 <- 0-(kca1*kca2*caa+(kca1+kca2)*cab)
	func_ <- vector(mode='numeric',len=5)
	func_[1] <- BoundCa / (s + Total_CaM) # CaPerCaM
	func_[2] <- 100*(Total_pCaMKII / (s + totalCaMKII)) # AutoCaMKII
	func_[3] <- Total_PP2B_CaM_CaX / (s + PP2B + Total_PP2B_CaM_CaX) # CaMPerPP2B
	func_[4] <- 100 * (PP2B_CaM_Ca4 / (s + PP2B + Total_PP2B_CaM_CaX)) # ActivePP2B
	func_[5] <- Ca # Camonitor
	names(func_) <- c("CaPerCaM", "AutoCaMKII", "CaMPerPP2B", "ActivePP2B", "Camonitor")
	return(func_)
}
# output function Jacobian dF(t,y;p)/dp
CaMKIIs_funcJac<-function(t, state, parameters)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	CaM_Ca1 <- state[1]
	CaM_Ca2 <- state[2]
	CaM_Ca3 <- state[3]
	CaM_Ca4 <- state[4]
	PP2B_CaM <- state[5]
	PP2B_CaM_Ca1 <- state[6]
	PP2B_CaM_Ca2 <- state[7]
	PP2B_CaM_Ca3 <- state[8]
	PP2B_CaM_Ca4 <- state[9]
	CaMKII_CaM <- state[10]
	CaMKII_CaM_Ca1 <- state[11]
	CaMKII_CaM_Ca2 <- state[12]
	CaMKII_CaM_Ca3 <- state[13]
	CaMKII_CaM_Ca4 <- state[14]
	pCaMKII_CaM_Ca4 <- state[15]
	pCaMKIIa <- state[16]
	pCaMKII_CaM_Ca3 <- state[17]
	pCaMKII_CaM_Ca2 <- state[18]
	pCaMKII_CaM_Ca1 <- state[19]
	pCaMKII_CaM <- state[20]
	PP1__pCaMKIIa <- state[21]
	caa <- state[22]
	cab <- state[23]
	KD__CaM_Ca3__PP2B<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca
	KD__CaM_Ca2__PP2B<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca
	KD__CaM_Ca1__PP2B<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca
	KD__CaM__PP2B<-(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca
	KD__CaM_Ca3__CaMKII<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__CaMKII<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__CaMKII<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca
	KD__CaM__CaMKII<-(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca
	KD__CaM_Ca3__pCaMKIIa<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__pCaMKIIa<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__pCaMKIIa<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca
	KD__CaM__pCaMKIIa<-(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca
	kr__CaM_Ca3__Ca<-KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca
	kr__CaM_Ca2__Ca<-KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca
	kr__CaM_Ca1__Ca<-KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca
	kr__CaM__Ca<-KD__CaM__Ca * kf__CaM__Ca
	kr__CaM_Ca4__PP2B<-KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B
	kr__CaM_Ca3__PP2B<-KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B
	kr__CaM_Ca2__PP2B<-KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B
	kr__CaM_Ca1__PP2B<-KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B
	kr__CaM__PP2B<-KD__CaM__PP2B * kf__CaM__PP2B
	kr__PP2B_CaM_Ca3__Ca<-KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca
	kr__PP2B_CaM_Ca2__Ca<-KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca
	kr__PP2B_CaM_Ca1__Ca<-KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca
	kr__PP2B_CaM__Ca<-KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca
	kr__CaM_Ca4__CaMKII<-KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII
	kr__CaM_Ca3__CaMKII<-KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII
	kr__CaM_Ca2__CaMKII<-KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII
	kr__CaM_Ca1__CaMKII<-KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII
	kr__CaM__CaMKII<-KD__CaM__CaMKII * kf__CaM__CaMKII
	kr__CaMKII_CaM_Ca3__Ca<-KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca
	kr__CaMKII_CaM_Ca2__Ca<-KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca
	kr__CaMKII_CaM_Ca1__Ca<-KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca
	kr__CaMKII_CaM__Ca<-KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca
	kr__pCaMKII_CaM_Ca3__Ca<-KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca
	kr__pCaMKII_CaM_Ca2__Ca<-KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca
	kr__pCaMKII_CaM_Ca1__Ca<-KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca
	kr__pCaMKII_CaM__Ca<-KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca
	kr__CaM_Ca4__pCaMKIIa<-KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa
	kr__CaM_Ca3__pCaMKIIa<-KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa
	kr__CaM_Ca2__pCaMKIIa<-KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa
	kr__CaM_Ca1__pCaMKIIa<-KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa
	kr__CaM__pCaMKIIa<-KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa
	Ca<-isOn*Ca_set+caa
	PP1<-PP1_0-PP1__pCaMKIIa
	CaMKII<-CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1)
	CaM<-CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa)
	PP2B<-PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4)
	totalCaMKII<-(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM)
	ActiveCaMKII<-pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4
	r<-(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII)
	pairedCaMKIIc<-a*(r*r)/(1+b*r)
	CaM_bound_Ca<-(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4)
	PP2B_bound_Ca<-(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4)
	CaMKII_bound_Ca<-(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4))
	BoundCa<-(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca)
	Total_CaM_CaX<-CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4
	Total_PP2B_CaM_CaX<-PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4
	Total_CaMKII_CaM_CaX<-CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4
	Total_pCaMKII_CaM_CaX<-pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4
	Total_CaM<-Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX
	Total_pCaMKII<-pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM
	ReactionFlux1<-(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1)
	ReactionFlux2<-(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2)
	ReactionFlux3<-(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3)
	ReactionFlux4<-(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4)
	ReactionFlux5<-(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM)
	ReactionFlux6<-(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1)
	ReactionFlux7<-(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2)
	ReactionFlux8<-(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3)
	ReactionFlux9<-(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4)
	ReactionFlux10<-(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1)
	ReactionFlux11<-(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2)
	ReactionFlux12<-(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3)
	ReactionFlux13<-(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4)
	ReactionFlux14<-(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM)
	ReactionFlux15<-(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1)
	ReactionFlux16<-(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2)
	ReactionFlux17<-(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3)
	ReactionFlux18<-(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4)
	ReactionFlux19<-(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1)
	ReactionFlux20<-(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2)
	ReactionFlux21<-(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3)
	ReactionFlux22<-(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4)
	ReactionFlux23<-(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)
	ReactionFlux24<-(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3)
	ReactionFlux25<-(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2)
	ReactionFlux26<-(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1)
	ReactionFlux27<-(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM)
	ReactionFlux28<-(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1)
	ReactionFlux29<-(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2)
	ReactionFlux30<-(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3)
	ReactionFlux31<-(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4)
	ReactionFlux32<-(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4)
	ReactionFlux33<-(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa)
	ReactionFlux34<-(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa)
	SpikeFlux1<-cab
	SpikeFlux2<-0-(kca1*kca2*caa+(kca1+kca2)*cab)
	fjac_ <- matrix(0.0,5,23)
# column 1
# column 2
# column 3
# column 4
# column 5
	fjac_[3,5] <- 1/(s+PP2B_0)
# column 6
	fjac_[3,6] <- 1/(s+PP2B_0)
# column 7
	fjac_[3,7] <- 1/(s+PP2B_0)
# column 8
	fjac_[3,8] <- 1/(s+PP2B_0)
# column 9
	fjac_[3,9] <- 1/(s+PP2B_0)
	fjac_[4,9] <- 100/(s+PP2B_0)
# column 10
	fjac_[2,10] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 11
	fjac_[2,11] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 12
	fjac_[2,12] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 13
	fjac_[2,13] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 14
	fjac_[2,14] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 15
	fjac_[2,15] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 16
	fjac_[2,16] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 17
	fjac_[2,17] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 18
	fjac_[2,18] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 19
	fjac_[2,19] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 20
	fjac_[2,20] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
# column 21
# column 22
	fjac_[5,22] <- 1
# column 23
	rownames(fjac_) <- c("CaPerCaM", "AutoCaMKII", "CaMPerPP2B", "ActivePP2B", "Camonitor")
	colnames(fjac_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa", "caa", "cab")
	return(fjac_)
}
# output function parameter Jacobian dF(t,y;p)/dp
CaMKIIs_funcJacp<-function(t, state, parameters)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	CaM_Ca1 <- state[1]
	CaM_Ca2 <- state[2]
	CaM_Ca3 <- state[3]
	CaM_Ca4 <- state[4]
	PP2B_CaM <- state[5]
	PP2B_CaM_Ca1 <- state[6]
	PP2B_CaM_Ca2 <- state[7]
	PP2B_CaM_Ca3 <- state[8]
	PP2B_CaM_Ca4 <- state[9]
	CaMKII_CaM <- state[10]
	CaMKII_CaM_Ca1 <- state[11]
	CaMKII_CaM_Ca2 <- state[12]
	CaMKII_CaM_Ca3 <- state[13]
	CaMKII_CaM_Ca4 <- state[14]
	pCaMKII_CaM_Ca4 <- state[15]
	pCaMKIIa <- state[16]
	pCaMKII_CaM_Ca3 <- state[17]
	pCaMKII_CaM_Ca2 <- state[18]
	pCaMKII_CaM_Ca1 <- state[19]
	pCaMKII_CaM <- state[20]
	PP1__pCaMKIIa <- state[21]
	caa <- state[22]
	cab <- state[23]
	KD__CaM_Ca3__PP2B<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca
	KD__CaM_Ca2__PP2B<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca
	KD__CaM_Ca1__PP2B<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca
	KD__CaM__PP2B<-(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca
	KD__CaM_Ca3__CaMKII<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__CaMKII<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__CaMKII<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca
	KD__CaM__CaMKII<-(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca
	KD__CaM_Ca3__pCaMKIIa<-(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca
	KD__CaM_Ca2__pCaMKIIa<-(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca
	KD__CaM_Ca1__pCaMKIIa<-(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca
	KD__CaM__pCaMKIIa<-(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca
	kr__CaM_Ca3__Ca<-KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca
	kr__CaM_Ca2__Ca<-KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca
	kr__CaM_Ca1__Ca<-KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca
	kr__CaM__Ca<-KD__CaM__Ca * kf__CaM__Ca
	kr__CaM_Ca4__PP2B<-KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B
	kr__CaM_Ca3__PP2B<-KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B
	kr__CaM_Ca2__PP2B<-KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B
	kr__CaM_Ca1__PP2B<-KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B
	kr__CaM__PP2B<-KD__CaM__PP2B * kf__CaM__PP2B
	kr__PP2B_CaM_Ca3__Ca<-KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca
	kr__PP2B_CaM_Ca2__Ca<-KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca
	kr__PP2B_CaM_Ca1__Ca<-KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca
	kr__PP2B_CaM__Ca<-KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca
	kr__CaM_Ca4__CaMKII<-KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII
	kr__CaM_Ca3__CaMKII<-KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII
	kr__CaM_Ca2__CaMKII<-KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII
	kr__CaM_Ca1__CaMKII<-KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII
	kr__CaM__CaMKII<-KD__CaM__CaMKII * kf__CaM__CaMKII
	kr__CaMKII_CaM_Ca3__Ca<-KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca
	kr__CaMKII_CaM_Ca2__Ca<-KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca
	kr__CaMKII_CaM_Ca1__Ca<-KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca
	kr__CaMKII_CaM__Ca<-KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca
	kr__pCaMKII_CaM_Ca3__Ca<-KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca
	kr__pCaMKII_CaM_Ca2__Ca<-KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca
	kr__pCaMKII_CaM_Ca1__Ca<-KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca
	kr__pCaMKII_CaM__Ca<-KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca
	kr__CaM_Ca4__pCaMKIIa<-KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa
	kr__CaM_Ca3__pCaMKIIa<-KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa
	kr__CaM_Ca2__pCaMKIIa<-KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa
	kr__CaM_Ca1__pCaMKIIa<-KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa
	kr__CaM__pCaMKIIa<-KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa
	Ca<-isOn*Ca_set+caa
	PP1<-PP1_0-PP1__pCaMKIIa
	CaMKII<-CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1)
	CaM<-CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa)
	PP2B<-PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4)
	totalCaMKII<-(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM)
	ActiveCaMKII<-pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4
	r<-(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII)
	pairedCaMKIIc<-a*(r*r)/(1+b*r)
	CaM_bound_Ca<-(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4)
	PP2B_bound_Ca<-(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4)
	CaMKII_bound_Ca<-(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4))
	BoundCa<-(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca)
	Total_CaM_CaX<-CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4
	Total_PP2B_CaM_CaX<-PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4
	Total_CaMKII_CaM_CaX<-CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4
	Total_pCaMKII_CaM_CaX<-pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4
	Total_CaM<-Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX
	Total_pCaMKII<-pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM
	ReactionFlux1<-(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1)
	ReactionFlux2<-(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2)
	ReactionFlux3<-(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3)
	ReactionFlux4<-(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4)
	ReactionFlux5<-(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM)
	ReactionFlux6<-(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1)
	ReactionFlux7<-(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2)
	ReactionFlux8<-(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3)
	ReactionFlux9<-(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4)
	ReactionFlux10<-(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1)
	ReactionFlux11<-(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2)
	ReactionFlux12<-(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3)
	ReactionFlux13<-(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4)
	ReactionFlux14<-(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM)
	ReactionFlux15<-(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1)
	ReactionFlux16<-(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2)
	ReactionFlux17<-(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3)
	ReactionFlux18<-(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4)
	ReactionFlux19<-(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1)
	ReactionFlux20<-(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2)
	ReactionFlux21<-(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3)
	ReactionFlux22<-(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4)
	ReactionFlux23<-(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)
	ReactionFlux24<-(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3)
	ReactionFlux25<-(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2)
	ReactionFlux26<-(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1)
	ReactionFlux27<-(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM)
	ReactionFlux28<-(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1)
	ReactionFlux29<-(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2)
	ReactionFlux30<-(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3)
	ReactionFlux31<-(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4)
	ReactionFlux32<-(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4)
	ReactionFlux33<-(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa)
	ReactionFlux34<-(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa)
	SpikeFlux1<-cab
	SpikeFlux2<-0-(kca1*kca2*caa+(kca1+kca2)*cab)
	fjacp_ <- matrix(0.0,5,62)
# column 1
# column 2
# column 3
# column 4
# column 5
# column 6
# column 7
# column 8
# column 9
# column 10
# column 11
# column 12
# column 13
# column 14
# column 15
# column 16
# column 17
# column 18
# column 19
# column 20
# column 21
# column 22
# column 23
# column 24
# column 25
# column 26
# column 27
# column 28
# column 29
# column 30
# column 31
# column 32
# column 33
# column 34
# column 35
# column 36
# column 37
# column 38
# column 39
# column 40
# column 41
# column 42
# column 43
# column 44
# column 45
# column 46
# column 47
# column 48
# column 49
# column 50
# column 51
# column 52
# column 53
# column 54
# column 55
	fjacp_[5,55] <- isOn
# column 56
# column 57
# column 58
# column 59
	fjacp_[3,59] <- -(PP2B_CaM_Ca4+PP2B_CaM_Ca3+PP2B_CaM_Ca2+PP2B_CaM_Ca1+PP2B_CaM)/(s+PP2B_0)^2
	fjacp_[4,59] <- -(100*PP2B_CaM_Ca4)/(s+PP2B_0)^2
# column 60
# column 61
# column 62
	fjacp_[5,62] <- Ca_set
	rownames(fjacp_) <- c("CaPerCaM", "AutoCaMKII", "CaMPerPP2B", "ActivePP2B", "Camonitor")
	colnames(fjacp_) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "Ca_set", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0", "kca1", "kca2", "isOn")
	return(fjacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
CaMKIIs_default<-function(t=0.0)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	parameters <- vector(mode='numeric',len=62)
	parameters[1] <- 0.151909999999999
	parameters[2] <- 3.42450000000002e-05
	parameters[3] <- 0.0738930000000004
	parameters[4] <- 0.0061814
	parameters[5] <- 0.0202929999999999
	parameters[6] <- 0.00452479999999998
	parameters[7] <- 0.051176
	parameters[8] <- 0.274210000000001
	parameters[9] <- 0.0833569999999997
	parameters[10] <- 0.0011578
	parameters[11] <- 0.0047884
	parameters[12] <- 0.0350789999999999
	parameters[13] <- 0.0455659999999999
	parameters[14] <- 7271.30000000002
	parameters[15] <- 37062.0000000017
	parameters[16] <- 1827.9
	parameters[17] <- 2662.30000000001
	parameters[18] <- 0.0399700000000002
	parameters[19] <- 91.5429999999998
	parameters[20] <- 916.149999999996
	parameters[21] <- 285.030000000001
	parameters[22] <- 82.8369999999998
	parameters[23] <- 0.237449999999999
	parameters[24] <- 0.0258579999999999
	parameters[25] <- 0.130860000000001
	parameters[26] <- 0.0755390000000001
	parameters[27] <- 0.000797720000000003
	parameters[28] <- 0.0558820000000003
	parameters[29] <- 0.0460280000000002
	parameters[30] <- 0.208550000000001
	parameters[31] <- 0.0226620000000001
	parameters[32] <- 8.28490000000001
	parameters[33] <- 483.480000000002
	parameters[34] <- 1143.6
	parameters[35] <- 645.069999999998
	parameters[36] <- 3081.60000000001
	parameters[37] <- 0.000829840000000001
	parameters[38] <- 0.00032583
	parameters[39] <- 0.0589279999999999
	parameters[40] <- 0.0231899999999999
	parameters[41] <- 0.0302520000000001
	parameters[42] <- 0.0384980000000001
	parameters[43] <- 0.0004565
	parameters[44] <- 0.0722370000000001
	parameters[45] <- 0.00212670000000001
	parameters[46] <- 539.41
	parameters[47] <- 1784.00000000001
	parameters[48] <- 57728.0000000003
	parameters[49] <- 1342.9
	parameters[50] <- 3.74600000000001
	parameters[51] <- 0.00375589999999998
	parameters[52] <- 0.0016604
	parameters[53] <- 0.205170000000001
	parameters[54] <- 0.302249999999999
	parameters[55] <- 2187.8
	parameters[56] <- 0
	parameters[57] <- 0
	parameters[58] <- 30
	parameters[59] <- 3
	parameters[60] <- 0.00978422
	parameters[61] <- 0.03448
	parameters[62] <- 0
	names(parameters) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "Ca_set", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0", "kca1", "kca2", "isOn")
	return(parameters);
}
# ode initial values
CaMKIIs_init<-function(t=0.0, parameters=NA)
{
	a<-3.90264
	b<-2.86972
	s<-1e-05
	kf__CaM__Ca <- parameters[1]
	kf__CaM_Ca1__Ca <- parameters[2]
	kf__CaM_Ca2__Ca <- parameters[3]
	kf__CaM_Ca3__Ca <- parameters[4]
	kf__CaM__PP2B <- parameters[5]
	kf__CaM_Ca1__PP2B <- parameters[6]
	kf__CaM_Ca2__PP2B <- parameters[7]
	kf__CaM_Ca3__PP2B <- parameters[8]
	kf__CaM_Ca4__PP2B <- parameters[9]
	kf__PP2B_CaM__Ca <- parameters[10]
	kf__PP2B_CaM_Ca1__Ca <- parameters[11]
	kf__PP2B_CaM_Ca2__Ca <- parameters[12]
	kf__PP2B_CaM_Ca3__Ca <- parameters[13]
	KD__CaM_Ca3__Ca <- parameters[14]
	KD__CaM_Ca2__Ca <- parameters[15]
	KD__CaM_Ca1__Ca <- parameters[16]
	KD__CaM__Ca <- parameters[17]
	KD__CaM_Ca4__PP2B <- parameters[18]
	KD__PP2B_CaM_Ca3__Ca <- parameters[19]
	KD__PP2B_CaM_Ca2__Ca <- parameters[20]
	KD__PP2B_CaM_Ca1__Ca <- parameters[21]
	KD__PP2B_CaM__Ca <- parameters[22]
	kf__CaM__CaMKII <- parameters[23]
	kf__CaMKII_CaM_Ca3__Ca <- parameters[24]
	kf__CaMKII_CaM_Ca2__Ca <- parameters[25]
	kf__CaMKII_CaM_Ca1__Ca <- parameters[26]
	kf__CaMKII_CaM__Ca <- parameters[27]
	kf__CaM_Ca1__CaMKII <- parameters[28]
	kf__CaM_Ca2__CaMKII <- parameters[29]
	kf__CaM_Ca3__CaMKII <- parameters[30]
	kf__CaM_Ca4__CaMKII <- parameters[31]
	KD__CaM_Ca4__CaMKII <- parameters[32]
	KD__CaMKII_CaM_Ca3__Ca <- parameters[33]
	KD__CaMKII_CaM_Ca2__Ca <- parameters[34]
	KD__CaMKII_CaM_Ca1__Ca <- parameters[35]
	KD__CaMKII_CaM__Ca <- parameters[36]
	kf__pCaMKII_CaM_Ca3__Ca <- parameters[37]
	kf__CaM__pCaMKIIa <- parameters[38]
	kf__CaM_Ca1__pCaMKIIa <- parameters[39]
	kf__CaM_Ca2__pCaMKIIa <- parameters[40]
	kf__CaM_Ca3__pCaMKIIa <- parameters[41]
	kf__pCaMKII_CaM_Ca2__Ca <- parameters[42]
	kf__pCaMKII_CaM_Ca1__Ca <- parameters[43]
	kf__CaM_Ca4__pCaMKIIa <- parameters[44]
	kf__pCaMKII_CaM__Ca <- parameters[45]
	KD__pCaMKII_CaM_Ca3__Ca <- parameters[46]
	KD__pCaMKII_CaM_Ca2__Ca <- parameters[47]
	KD__pCaMKII_CaM_Ca1__Ca <- parameters[48]
	KD__pCaMKII_CaM__Ca <- parameters[49]
	KD__CaM_Ca4__pCaMKIIa <- parameters[50]
	kautMax <- parameters[51]
	kf__PP1__pCaMKIIa <- parameters[52]
	kr__PP1__pCaMKIIa <- parameters[53]
	kcat__PP1__pCaMKIIa <- parameters[54]
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	kca1 <- parameters[60]
	kca2 <- parameters[61]
	isOn <- parameters[62]
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=23)
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
	state[12] <- 0
	state[13] <- 0
	state[14] <- 0
	state[15] <- 0
	state[16] <- 0
	state[17] <- 0
	state[18] <- 0
	state[19] <- 0
	state[20] <- 0
	state[21] <- 0
	state[22] <- 0
	state[23] <- 0
	names(state) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa", "caa", "cab")
	return(state)
}
model<-list(vf=CaMKIIs_vf, jac=CaMKIIs_jac, jacp=CaMKIIs_jacp, func=CaMKIIs_func, funcJac=CaMKIIs_funcJac, funcJacp=CaMKIIs_funcJacp, init=CaMKIIs_init, par=CaMKIIs_default, name="CaMKIIs")
