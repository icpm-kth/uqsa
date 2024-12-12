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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
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
	Ca <- isOn*CaBase
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
	f_<-vector(mode='numeric',len=21)
	f_[1] <- CaM_Ca1	+ReactionFlux1-ReactionFlux2-ReactionFlux6-ReactionFlux15-ReactionFlux26
	f_[2] <- CaM_Ca2	+ReactionFlux2-ReactionFlux3-ReactionFlux7-ReactionFlux16-ReactionFlux25
	f_[3] <- CaM_Ca3	+ReactionFlux3-ReactionFlux4-ReactionFlux8-ReactionFlux17-ReactionFlux24
	f_[4] <- CaM_Ca4	+ReactionFlux4-ReactionFlux9-ReactionFlux18-ReactionFlux23
	f_[5] <- PP2B_CaM	+ReactionFlux5-ReactionFlux10
	f_[6] <- PP2B_CaM_Ca1	+ReactionFlux6+ReactionFlux10-ReactionFlux11
	f_[7] <- PP2B_CaM_Ca2	+ReactionFlux7+ReactionFlux11-ReactionFlux12
	f_[8] <- PP2B_CaM_Ca3	+ReactionFlux8+ReactionFlux12-ReactionFlux13
	f_[9] <- PP2B_CaM_Ca4	+ReactionFlux9+ReactionFlux13
	f_[10] <- CaMKII_CaM	+ReactionFlux14-ReactionFlux19
	f_[11] <- CaMKII_CaM_Ca1	+ReactionFlux15+ReactionFlux19-ReactionFlux20
	f_[12] <- CaMKII_CaM_Ca2	+ReactionFlux16+ReactionFlux20-ReactionFlux21
	f_[13] <- CaMKII_CaM_Ca3	+ReactionFlux17+ReactionFlux21-ReactionFlux22
	f_[14] <- CaMKII_CaM_Ca4	+ReactionFlux18+ReactionFlux22-ReactionFlux32
	f_[15] <- pCaMKII_CaM_Ca4	+ReactionFlux23+ReactionFlux31+ReactionFlux32
	f_[16] <- pCaMKIIa	-ReactionFlux23-ReactionFlux24-ReactionFlux25-ReactionFlux26-ReactionFlux27-ReactionFlux33
	f_[17] <- pCaMKII_CaM_Ca3	+ReactionFlux24+ReactionFlux30-ReactionFlux31
	f_[18] <- pCaMKII_CaM_Ca2	+ReactionFlux25+ReactionFlux29-ReactionFlux30
	f_[19] <- pCaMKII_CaM_Ca1	+ReactionFlux26+ReactionFlux28-ReactionFlux29
	f_[20] <- pCaMKII_CaM	+ReactionFlux27-ReactionFlux28
	f_[21] <- PP1__pCaMKIIa	+ReactionFlux33-ReactionFlux34
	names(f_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
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
	Ca <- isOn*CaBase
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
	jac_ <- matrix(NA,21,21)
# column 1
	jac_[1,1] <- (-kf__CaM_Ca1__pCaMKIIa*pCaMKIIa)-kr__CaM__Ca-PP2B*kf__CaM_Ca1__PP2B-CaMKII*kf__CaM_Ca1__CaMKII-Ca*kf__CaM_Ca1__Ca
	jac_[2,1] <- Ca*kf__CaM_Ca1__Ca
	jac_[3,1] <- 0
	jac_[4,1] <- 0
	jac_[5,1] <- 0
	jac_[6,1] <- PP2B*kf__CaM_Ca1__PP2B
	jac_[7,1] <- 0
	jac_[8,1] <- 0
	jac_[9,1] <- 0
	jac_[10,1] <- 0
	jac_[11,1] <- CaMKII*kf__CaM_Ca1__CaMKII
	jac_[12,1] <- 0
	jac_[13,1] <- 0
	jac_[14,1] <- 0
	jac_[15,1] <- 0
	jac_[16,1] <- -kf__CaM_Ca1__pCaMKIIa*pCaMKIIa
	jac_[17,1] <- 0
	jac_[18,1] <- 0
	jac_[19,1] <- kf__CaM_Ca1__pCaMKIIa*pCaMKIIa
	jac_[20,1] <- 0
	jac_[21,1] <- 0
# column 2
	jac_[1,2] <- kr__CaM_Ca1__Ca
	jac_[2,2] <- (-kf__CaM_Ca2__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca1__Ca-PP2B*kf__CaM_Ca2__PP2B-CaMKII*kf__CaM_Ca2__CaMKII-Ca*kf__CaM_Ca2__Ca
	jac_[3,2] <- Ca*kf__CaM_Ca2__Ca
	jac_[4,2] <- 0
	jac_[5,2] <- 0
	jac_[6,2] <- 0
	jac_[7,2] <- PP2B*kf__CaM_Ca2__PP2B
	jac_[8,2] <- 0
	jac_[9,2] <- 0
	jac_[10,2] <- 0
	jac_[11,2] <- 0
	jac_[12,2] <- CaMKII*kf__CaM_Ca2__CaMKII
	jac_[13,2] <- 0
	jac_[14,2] <- 0
	jac_[15,2] <- 0
	jac_[16,2] <- -kf__CaM_Ca2__pCaMKIIa*pCaMKIIa
	jac_[17,2] <- 0
	jac_[18,2] <- kf__CaM_Ca2__pCaMKIIa*pCaMKIIa
	jac_[19,2] <- 0
	jac_[20,2] <- 0
	jac_[21,2] <- 0
# column 3
	jac_[1,3] <- 0
	jac_[2,3] <- kr__CaM_Ca2__Ca
	jac_[3,3] <- (-kf__CaM_Ca3__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca2__Ca-PP2B*kf__CaM_Ca3__PP2B-CaMKII*kf__CaM_Ca3__CaMKII-Ca*kf__CaM_Ca3__Ca
	jac_[4,3] <- Ca*kf__CaM_Ca3__Ca
	jac_[5,3] <- 0
	jac_[6,3] <- 0
	jac_[7,3] <- 0
	jac_[8,3] <- PP2B*kf__CaM_Ca3__PP2B
	jac_[9,3] <- 0
	jac_[10,3] <- 0
	jac_[11,3] <- 0
	jac_[12,3] <- 0
	jac_[13,3] <- CaMKII*kf__CaM_Ca3__CaMKII
	jac_[14,3] <- 0
	jac_[15,3] <- 0
	jac_[16,3] <- -kf__CaM_Ca3__pCaMKIIa*pCaMKIIa
	jac_[17,3] <- kf__CaM_Ca3__pCaMKIIa*pCaMKIIa
	jac_[18,3] <- 0
	jac_[19,3] <- 0
	jac_[20,3] <- 0
	jac_[21,3] <- 0
# column 4
	jac_[1,4] <- 0
	jac_[2,4] <- 0
	jac_[3,4] <- kr__CaM_Ca3__Ca
	jac_[4,4] <- (-kf__CaM_Ca4__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca3__Ca-PP2B*kf__CaM_Ca4__PP2B-CaMKII*kf__CaM_Ca4__CaMKII
	jac_[5,4] <- 0
	jac_[6,4] <- 0
	jac_[7,4] <- 0
	jac_[8,4] <- 0
	jac_[9,4] <- PP2B*kf__CaM_Ca4__PP2B
	jac_[10,4] <- 0
	jac_[11,4] <- 0
	jac_[12,4] <- 0
	jac_[13,4] <- 0
	jac_[14,4] <- CaMKII*kf__CaM_Ca4__CaMKII
	jac_[15,4] <- kf__CaM_Ca4__pCaMKIIa*pCaMKIIa
	jac_[16,4] <- -kf__CaM_Ca4__pCaMKIIa*pCaMKIIa
	jac_[17,4] <- 0
	jac_[18,4] <- 0
	jac_[19,4] <- 0
	jac_[20,4] <- 0
	jac_[21,4] <- 0
# column 5
	jac_[1,5] <- 0
	jac_[2,5] <- 0
	jac_[3,5] <- 0
	jac_[4,5] <- 0
	jac_[5,5] <- (-kr__CaM__PP2B)-Ca*kf__PP2B_CaM__Ca
	jac_[6,5] <- Ca*kf__PP2B_CaM__Ca
	jac_[7,5] <- 0
	jac_[8,5] <- 0
	jac_[9,5] <- 0
	jac_[10,5] <- 0
	jac_[11,5] <- 0
	jac_[12,5] <- 0
	jac_[13,5] <- 0
	jac_[14,5] <- 0
	jac_[15,5] <- 0
	jac_[16,5] <- 0
	jac_[17,5] <- 0
	jac_[18,5] <- 0
	jac_[19,5] <- 0
	jac_[20,5] <- 0
	jac_[21,5] <- 0
# column 6
	jac_[1,6] <- kr__CaM_Ca1__PP2B
	jac_[2,6] <- 0
	jac_[3,6] <- 0
	jac_[4,6] <- 0
	jac_[5,6] <- kr__PP2B_CaM__Ca
	jac_[6,6] <- (-kr__PP2B_CaM__Ca)-kr__CaM_Ca1__PP2B-Ca*kf__PP2B_CaM_Ca1__Ca
	jac_[7,6] <- Ca*kf__PP2B_CaM_Ca1__Ca
	jac_[8,6] <- 0
	jac_[9,6] <- 0
	jac_[10,6] <- 0
	jac_[11,6] <- 0
	jac_[12,6] <- 0
	jac_[13,6] <- 0
	jac_[14,6] <- 0
	jac_[15,6] <- 0
	jac_[16,6] <- 0
	jac_[17,6] <- 0
	jac_[18,6] <- 0
	jac_[19,6] <- 0
	jac_[20,6] <- 0
	jac_[21,6] <- 0
# column 7
	jac_[1,7] <- 0
	jac_[2,7] <- kr__CaM_Ca2__PP2B
	jac_[3,7] <- 0
	jac_[4,7] <- 0
	jac_[5,7] <- 0
	jac_[6,7] <- kr__PP2B_CaM_Ca1__Ca
	jac_[7,7] <- (-kr__PP2B_CaM_Ca1__Ca)-kr__CaM_Ca2__PP2B-Ca*kf__PP2B_CaM_Ca2__Ca
	jac_[8,7] <- Ca*kf__PP2B_CaM_Ca2__Ca
	jac_[9,7] <- 0
	jac_[10,7] <- 0
	jac_[11,7] <- 0
	jac_[12,7] <- 0
	jac_[13,7] <- 0
	jac_[14,7] <- 0
	jac_[15,7] <- 0
	jac_[16,7] <- 0
	jac_[17,7] <- 0
	jac_[18,7] <- 0
	jac_[19,7] <- 0
	jac_[20,7] <- 0
	jac_[21,7] <- 0
# column 8
	jac_[1,8] <- 0
	jac_[2,8] <- 0
	jac_[3,8] <- kr__CaM_Ca3__PP2B
	jac_[4,8] <- 0
	jac_[5,8] <- 0
	jac_[6,8] <- 0
	jac_[7,8] <- kr__PP2B_CaM_Ca2__Ca
	jac_[8,8] <- (-kr__PP2B_CaM_Ca2__Ca)-kr__CaM_Ca3__PP2B-Ca*kf__PP2B_CaM_Ca3__Ca
	jac_[9,8] <- Ca*kf__PP2B_CaM_Ca3__Ca
	jac_[10,8] <- 0
	jac_[11,8] <- 0
	jac_[12,8] <- 0
	jac_[13,8] <- 0
	jac_[14,8] <- 0
	jac_[15,8] <- 0
	jac_[16,8] <- 0
	jac_[17,8] <- 0
	jac_[18,8] <- 0
	jac_[19,8] <- 0
	jac_[20,8] <- 0
	jac_[21,8] <- 0
# column 9
	jac_[1,9] <- 0
	jac_[2,9] <- 0
	jac_[3,9] <- 0
	jac_[4,9] <- kr__CaM_Ca4__PP2B
	jac_[5,9] <- 0
	jac_[6,9] <- 0
	jac_[7,9] <- 0
	jac_[8,9] <- kr__PP2B_CaM_Ca3__Ca
	jac_[9,9] <- (-kr__PP2B_CaM_Ca3__Ca)-kr__CaM_Ca4__PP2B
	jac_[10,9] <- 0
	jac_[11,9] <- 0
	jac_[12,9] <- 0
	jac_[13,9] <- 0
	jac_[14,9] <- 0
	jac_[15,9] <- 0
	jac_[16,9] <- 0
	jac_[17,9] <- 0
	jac_[18,9] <- 0
	jac_[19,9] <- 0
	jac_[20,9] <- 0
	jac_[21,9] <- 0
# column 10
	jac_[1,10] <- 0
	jac_[2,10] <- 0
	jac_[3,10] <- 0
	jac_[4,10] <- 0
	jac_[5,10] <- 0
	jac_[6,10] <- 0
	jac_[7,10] <- 0
	jac_[8,10] <- 0
	jac_[9,10] <- 0
	jac_[10,10] <- (-kr__CaM__CaMKII)-Ca*kf__CaMKII_CaM__Ca
	jac_[11,10] <- Ca*kf__CaMKII_CaM__Ca
	jac_[12,10] <- 0
	jac_[13,10] <- 0
	jac_[14,10] <- 0
	jac_[15,10] <- 0
	jac_[16,10] <- 0
	jac_[17,10] <- 0
	jac_[18,10] <- 0
	jac_[19,10] <- 0
	jac_[20,10] <- 0
	jac_[21,10] <- 0
# column 11
	jac_[1,11] <- kr__CaM_Ca1__CaMKII
	jac_[2,11] <- 0
	jac_[3,11] <- 0
	jac_[4,11] <- 0
	jac_[5,11] <- 0
	jac_[6,11] <- 0
	jac_[7,11] <- 0
	jac_[8,11] <- 0
	jac_[9,11] <- 0
	jac_[10,11] <- kr__CaMKII_CaM__Ca
	jac_[11,11] <- (-kr__CaM_Ca1__CaMKII)-kr__CaMKII_CaM__Ca-Ca*kf__CaMKII_CaM_Ca1__Ca
	jac_[12,11] <- Ca*kf__CaMKII_CaM_Ca1__Ca
	jac_[13,11] <- 0
	jac_[14,11] <- 0
	jac_[15,11] <- 0
	jac_[16,11] <- 0
	jac_[17,11] <- 0
	jac_[18,11] <- 0
	jac_[19,11] <- 0
	jac_[20,11] <- 0
	jac_[21,11] <- 0
# column 12
	jac_[1,12] <- 0
	jac_[2,12] <- kr__CaM_Ca2__CaMKII
	jac_[3,12] <- 0
	jac_[4,12] <- 0
	jac_[5,12] <- 0
	jac_[6,12] <- 0
	jac_[7,12] <- 0
	jac_[8,12] <- 0
	jac_[9,12] <- 0
	jac_[10,12] <- 0
	jac_[11,12] <- kr__CaMKII_CaM_Ca1__Ca
	jac_[12,12] <- (-kr__CaM_Ca2__CaMKII)-kr__CaMKII_CaM_Ca1__Ca-Ca*kf__CaMKII_CaM_Ca2__Ca
	jac_[13,12] <- Ca*kf__CaMKII_CaM_Ca2__Ca
	jac_[14,12] <- 0
	jac_[15,12] <- 0
	jac_[16,12] <- 0
	jac_[17,12] <- 0
	jac_[18,12] <- 0
	jac_[19,12] <- 0
	jac_[20,12] <- 0
	jac_[21,12] <- 0
# column 13
	jac_[1,13] <- 0
	jac_[2,13] <- 0
	jac_[3,13] <- kr__CaM_Ca3__CaMKII
	jac_[4,13] <- 0
	jac_[5,13] <- 0
	jac_[6,13] <- 0
	jac_[7,13] <- 0
	jac_[8,13] <- 0
	jac_[9,13] <- 0
	jac_[10,13] <- 0
	jac_[11,13] <- 0
	jac_[12,13] <- kr__CaMKII_CaM_Ca2__Ca
	jac_[13,13] <- (-kr__CaM_Ca3__CaMKII)-kr__CaMKII_CaM_Ca2__Ca-Ca*kf__CaMKII_CaM_Ca3__Ca
	jac_[14,13] <- Ca*kf__CaMKII_CaM_Ca3__Ca
	jac_[15,13] <- 0
	jac_[16,13] <- 0
	jac_[17,13] <- 0
	jac_[18,13] <- 0
	jac_[19,13] <- 0
	jac_[20,13] <- 0
	jac_[21,13] <- 0
# column 14
	jac_[1,14] <- 0
	jac_[2,14] <- 0
	jac_[3,14] <- 0
	jac_[4,14] <- kr__CaM_Ca4__CaMKII
	jac_[5,14] <- 0
	jac_[6,14] <- 0
	jac_[7,14] <- 0
	jac_[8,14] <- 0
	jac_[9,14] <- 0
	jac_[10,14] <- 0
	jac_[11,14] <- 0
	jac_[12,14] <- 0
	jac_[13,14] <- kr__CaMKII_CaM_Ca3__Ca
	jac_[14,14] <- (-kautMax*pairedCaMKIIc)-kr__CaM_Ca4__CaMKII-kr__CaMKII_CaM_Ca3__Ca
	jac_[15,14] <- kautMax*pairedCaMKIIc
	jac_[16,14] <- 0
	jac_[17,14] <- 0
	jac_[18,14] <- 0
	jac_[19,14] <- 0
	jac_[20,14] <- 0
	jac_[21,14] <- 0
# column 15
	jac_[1,15] <- 0
	jac_[2,15] <- 0
	jac_[3,15] <- 0
	jac_[4,15] <- kr__CaM_Ca4__pCaMKIIa
	jac_[5,15] <- 0
	jac_[6,15] <- 0
	jac_[7,15] <- 0
	jac_[8,15] <- 0
	jac_[9,15] <- 0
	jac_[10,15] <- 0
	jac_[11,15] <- 0
	jac_[12,15] <- 0
	jac_[13,15] <- 0
	jac_[14,15] <- 0
	jac_[15,15] <- (-kr__pCaMKII_CaM_Ca3__Ca)-kr__CaM_Ca4__pCaMKIIa
	jac_[16,15] <- kr__CaM_Ca4__pCaMKIIa
	jac_[17,15] <- kr__pCaMKII_CaM_Ca3__Ca
	jac_[18,15] <- 0
	jac_[19,15] <- 0
	jac_[20,15] <- 0
	jac_[21,15] <- 0
# column 16
	jac_[1,16] <- -CaM_Ca1*kf__CaM_Ca1__pCaMKIIa
	jac_[2,16] <- -CaM_Ca2*kf__CaM_Ca2__pCaMKIIa
	jac_[3,16] <- -CaM_Ca3*kf__CaM_Ca3__pCaMKIIa
	jac_[4,16] <- -CaM_Ca4*kf__CaM_Ca4__pCaMKIIa
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
	jac_[15,16] <- CaM_Ca4*kf__CaM_Ca4__pCaMKIIa
	jac_[16,16] <- (-PP1*kf__PP1__pCaMKIIa)-CaM*kf__CaM__pCaMKIIa-CaM_Ca4*kf__CaM_Ca4__pCaMKIIa-CaM_Ca3*kf__CaM_Ca3__pCaMKIIa-CaM_Ca2*kf__CaM_Ca2__pCaMKIIa-CaM_Ca1*kf__CaM_Ca1__pCaMKIIa
	jac_[17,16] <- CaM_Ca3*kf__CaM_Ca3__pCaMKIIa
	jac_[18,16] <- CaM_Ca2*kf__CaM_Ca2__pCaMKIIa
	jac_[19,16] <- CaM_Ca1*kf__CaM_Ca1__pCaMKIIa
	jac_[20,16] <- CaM*kf__CaM__pCaMKIIa
	jac_[21,16] <- PP1*kf__PP1__pCaMKIIa
# column 17
	jac_[1,17] <- 0
	jac_[2,17] <- 0
	jac_[3,17] <- kr__CaM_Ca3__pCaMKIIa
	jac_[4,17] <- 0
	jac_[5,17] <- 0
	jac_[6,17] <- 0
	jac_[7,17] <- 0
	jac_[8,17] <- 0
	jac_[9,17] <- 0
	jac_[10,17] <- 0
	jac_[11,17] <- 0
	jac_[12,17] <- 0
	jac_[13,17] <- 0
	jac_[14,17] <- 0
	jac_[15,17] <- Ca*kf__pCaMKII_CaM_Ca3__Ca
	jac_[16,17] <- kr__CaM_Ca3__pCaMKIIa
	jac_[17,17] <- (-kr__pCaMKII_CaM_Ca2__Ca)-kr__CaM_Ca3__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca3__Ca
	jac_[18,17] <- kr__pCaMKII_CaM_Ca2__Ca
	jac_[19,17] <- 0
	jac_[20,17] <- 0
	jac_[21,17] <- 0
# column 18
	jac_[1,18] <- 0
	jac_[2,18] <- kr__CaM_Ca2__pCaMKIIa
	jac_[3,18] <- 0
	jac_[4,18] <- 0
	jac_[5,18] <- 0
	jac_[6,18] <- 0
	jac_[7,18] <- 0
	jac_[8,18] <- 0
	jac_[9,18] <- 0
	jac_[10,18] <- 0
	jac_[11,18] <- 0
	jac_[12,18] <- 0
	jac_[13,18] <- 0
	jac_[14,18] <- 0
	jac_[15,18] <- 0
	jac_[16,18] <- kr__CaM_Ca2__pCaMKIIa
	jac_[17,18] <- Ca*kf__pCaMKII_CaM_Ca2__Ca
	jac_[18,18] <- (-kr__pCaMKII_CaM_Ca1__Ca)-kr__CaM_Ca2__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca2__Ca
	jac_[19,18] <- kr__pCaMKII_CaM_Ca1__Ca
	jac_[20,18] <- 0
	jac_[21,18] <- 0
# column 19
	jac_[1,19] <- kr__CaM_Ca1__pCaMKIIa
	jac_[2,19] <- 0
	jac_[3,19] <- 0
	jac_[4,19] <- 0
	jac_[5,19] <- 0
	jac_[6,19] <- 0
	jac_[7,19] <- 0
	jac_[8,19] <- 0
	jac_[9,19] <- 0
	jac_[10,19] <- 0
	jac_[11,19] <- 0
	jac_[12,19] <- 0
	jac_[13,19] <- 0
	jac_[14,19] <- 0
	jac_[15,19] <- 0
	jac_[16,19] <- kr__CaM_Ca1__pCaMKIIa
	jac_[17,19] <- 0
	jac_[18,19] <- Ca*kf__pCaMKII_CaM_Ca1__Ca
	jac_[19,19] <- (-kr__pCaMKII_CaM__Ca)-kr__CaM_Ca1__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca1__Ca
	jac_[20,19] <- kr__pCaMKII_CaM__Ca
	jac_[21,19] <- 0
# column 20
	jac_[1,20] <- 0
	jac_[2,20] <- 0
	jac_[3,20] <- 0
	jac_[4,20] <- 0
	jac_[5,20] <- 0
	jac_[6,20] <- 0
	jac_[7,20] <- 0
	jac_[8,20] <- 0
	jac_[9,20] <- 0
	jac_[10,20] <- 0
	jac_[11,20] <- 0
	jac_[12,20] <- 0
	jac_[13,20] <- 0
	jac_[14,20] <- 0
	jac_[15,20] <- 0
	jac_[16,20] <- kr__CaM__pCaMKIIa
	jac_[17,20] <- 0
	jac_[18,20] <- 0
	jac_[19,20] <- Ca*kf__pCaMKII_CaM__Ca
	jac_[20,20] <- (-kr__CaM__pCaMKIIa)-Ca*kf__pCaMKII_CaM__Ca
	jac_[21,20] <- 0
# column 21
	jac_[1,21] <- 0
	jac_[2,21] <- 0
	jac_[3,21] <- 0
	jac_[4,21] <- 0
	jac_[5,21] <- 0
	jac_[6,21] <- 0
	jac_[7,21] <- 0
	jac_[8,21] <- 0
	jac_[9,21] <- 0
	jac_[10,21] <- 0
	jac_[11,21] <- 0
	jac_[12,21] <- 0
	jac_[13,21] <- 0
	jac_[14,21] <- 0
	jac_[15,21] <- 0
	jac_[16,21] <- kr__PP1__pCaMKIIa
	jac_[17,21] <- 0
	jac_[18,21] <- 0
	jac_[19,21] <- 0
	jac_[20,21] <- 0
	jac_[21,21] <- (-kr__PP1__pCaMKIIa)-kcat__PP1__pCaMKIIa
	rownames(jac_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
	colnames(jac_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
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
	Ca<-isOn*CaBase
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
	jacp_ <- matrix(NA,21,60)
# column 1
	jacp_[1,1] <- Ca*CaM
	jacp_[2,1] <- 0
	jacp_[3,1] <- 0
	jacp_[4,1] <- 0
	jacp_[5,1] <- 0
	jacp_[6,1] <- 0
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
	jacp_[17,1] <- 0
	jacp_[18,1] <- 0
	jacp_[19,1] <- 0
	jacp_[20,1] <- 0
	jacp_[21,1] <- 0
# column 2
	jacp_[1,2] <- -Ca*CaM_Ca1
	jacp_[2,2] <- Ca*CaM_Ca1
	jacp_[3,2] <- 0
	jacp_[4,2] <- 0
	jacp_[5,2] <- 0
	jacp_[6,2] <- 0
	jacp_[7,2] <- 0
	jacp_[8,2] <- 0
	jacp_[9,2] <- 0
	jacp_[10,2] <- 0
	jacp_[11,2] <- 0
	jacp_[12,2] <- 0
	jacp_[13,2] <- 0
	jacp_[14,2] <- 0
	jacp_[15,2] <- 0
	jacp_[16,2] <- 0
	jacp_[17,2] <- 0
	jacp_[18,2] <- 0
	jacp_[19,2] <- 0
	jacp_[20,2] <- 0
	jacp_[21,2] <- 0
# column 3
	jacp_[1,3] <- 0
	jacp_[2,3] <- -Ca*CaM_Ca2
	jacp_[3,3] <- Ca*CaM_Ca2
	jacp_[4,3] <- 0
	jacp_[5,3] <- 0
	jacp_[6,3] <- 0
	jacp_[7,3] <- 0
	jacp_[8,3] <- 0
	jacp_[9,3] <- 0
	jacp_[10,3] <- 0
	jacp_[11,3] <- 0
	jacp_[12,3] <- 0
	jacp_[13,3] <- 0
	jacp_[14,3] <- 0
	jacp_[15,3] <- 0
	jacp_[16,3] <- 0
	jacp_[17,3] <- 0
	jacp_[18,3] <- 0
	jacp_[19,3] <- 0
	jacp_[20,3] <- 0
	jacp_[21,3] <- 0
# column 4
	jacp_[1,4] <- 0
	jacp_[2,4] <- 0
	jacp_[3,4] <- -Ca*CaM_Ca3
	jacp_[4,4] <- Ca*CaM_Ca3
	jacp_[5,4] <- 0
	jacp_[6,4] <- 0
	jacp_[7,4] <- 0
	jacp_[8,4] <- 0
	jacp_[9,4] <- 0
	jacp_[10,4] <- 0
	jacp_[11,4] <- 0
	jacp_[12,4] <- 0
	jacp_[13,4] <- 0
	jacp_[14,4] <- 0
	jacp_[15,4] <- 0
	jacp_[16,4] <- 0
	jacp_[17,4] <- 0
	jacp_[18,4] <- 0
	jacp_[19,4] <- 0
	jacp_[20,4] <- 0
	jacp_[21,4] <- 0
# column 5
	jacp_[1,5] <- 0
	jacp_[2,5] <- 0
	jacp_[3,5] <- 0
	jacp_[4,5] <- 0
	jacp_[5,5] <- CaM*PP2B
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
	jacp_[17,5] <- 0
	jacp_[18,5] <- 0
	jacp_[19,5] <- 0
	jacp_[20,5] <- 0
	jacp_[21,5] <- 0
# column 6
	jacp_[1,6] <- -CaM_Ca1*PP2B
	jacp_[2,6] <- 0
	jacp_[3,6] <- 0
	jacp_[4,6] <- 0
	jacp_[5,6] <- 0
	jacp_[6,6] <- CaM_Ca1*PP2B
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
	jacp_[17,6] <- 0
	jacp_[18,6] <- 0
	jacp_[19,6] <- 0
	jacp_[20,6] <- 0
	jacp_[21,6] <- 0
# column 7
	jacp_[1,7] <- 0
	jacp_[2,7] <- -CaM_Ca2*PP2B
	jacp_[3,7] <- 0
	jacp_[4,7] <- 0
	jacp_[5,7] <- 0
	jacp_[6,7] <- 0
	jacp_[7,7] <- CaM_Ca2*PP2B
	jacp_[8,7] <- 0
	jacp_[9,7] <- 0
	jacp_[10,7] <- 0
	jacp_[11,7] <- 0
	jacp_[12,7] <- 0
	jacp_[13,7] <- 0
	jacp_[14,7] <- 0
	jacp_[15,7] <- 0
	jacp_[16,7] <- 0
	jacp_[17,7] <- 0
	jacp_[18,7] <- 0
	jacp_[19,7] <- 0
	jacp_[20,7] <- 0
	jacp_[21,7] <- 0
# column 8
	jacp_[1,8] <- 0
	jacp_[2,8] <- 0
	jacp_[3,8] <- -CaM_Ca3*PP2B
	jacp_[4,8] <- 0
	jacp_[5,8] <- 0
	jacp_[6,8] <- 0
	jacp_[7,8] <- 0
	jacp_[8,8] <- CaM_Ca3*PP2B
	jacp_[9,8] <- 0
	jacp_[10,8] <- 0
	jacp_[11,8] <- 0
	jacp_[12,8] <- 0
	jacp_[13,8] <- 0
	jacp_[14,8] <- 0
	jacp_[15,8] <- 0
	jacp_[16,8] <- 0
	jacp_[17,8] <- 0
	jacp_[18,8] <- 0
	jacp_[19,8] <- 0
	jacp_[20,8] <- 0
	jacp_[21,8] <- 0
# column 9
	jacp_[1,9] <- 0
	jacp_[2,9] <- 0
	jacp_[3,9] <- 0
	jacp_[4,9] <- -CaM_Ca4*PP2B
	jacp_[5,9] <- 0
	jacp_[6,9] <- 0
	jacp_[7,9] <- 0
	jacp_[8,9] <- 0
	jacp_[9,9] <- CaM_Ca4*PP2B
	jacp_[10,9] <- 0
	jacp_[11,9] <- 0
	jacp_[12,9] <- 0
	jacp_[13,9] <- 0
	jacp_[14,9] <- 0
	jacp_[15,9] <- 0
	jacp_[16,9] <- 0
	jacp_[17,9] <- 0
	jacp_[18,9] <- 0
	jacp_[19,9] <- 0
	jacp_[20,9] <- 0
	jacp_[21,9] <- 0
# column 10
	jacp_[1,10] <- 0
	jacp_[2,10] <- 0
	jacp_[3,10] <- 0
	jacp_[4,10] <- 0
	jacp_[5,10] <- -Ca*PP2B_CaM
	jacp_[6,10] <- Ca*PP2B_CaM
	jacp_[7,10] <- 0
	jacp_[8,10] <- 0
	jacp_[9,10] <- 0
	jacp_[10,10] <- 0
	jacp_[11,10] <- 0
	jacp_[12,10] <- 0
	jacp_[13,10] <- 0
	jacp_[14,10] <- 0
	jacp_[15,10] <- 0
	jacp_[16,10] <- 0
	jacp_[17,10] <- 0
	jacp_[18,10] <- 0
	jacp_[19,10] <- 0
	jacp_[20,10] <- 0
	jacp_[21,10] <- 0
# column 11
	jacp_[1,11] <- 0
	jacp_[2,11] <- 0
	jacp_[3,11] <- 0
	jacp_[4,11] <- 0
	jacp_[5,11] <- 0
	jacp_[6,11] <- -Ca*PP2B_CaM_Ca1
	jacp_[7,11] <- Ca*PP2B_CaM_Ca1
	jacp_[8,11] <- 0
	jacp_[9,11] <- 0
	jacp_[10,11] <- 0
	jacp_[11,11] <- 0
	jacp_[12,11] <- 0
	jacp_[13,11] <- 0
	jacp_[14,11] <- 0
	jacp_[15,11] <- 0
	jacp_[16,11] <- 0
	jacp_[17,11] <- 0
	jacp_[18,11] <- 0
	jacp_[19,11] <- 0
	jacp_[20,11] <- 0
	jacp_[21,11] <- 0
# column 12
	jacp_[1,12] <- 0
	jacp_[2,12] <- 0
	jacp_[3,12] <- 0
	jacp_[4,12] <- 0
	jacp_[5,12] <- 0
	jacp_[6,12] <- 0
	jacp_[7,12] <- -Ca*PP2B_CaM_Ca2
	jacp_[8,12] <- Ca*PP2B_CaM_Ca2
	jacp_[9,12] <- 0
	jacp_[10,12] <- 0
	jacp_[11,12] <- 0
	jacp_[12,12] <- 0
	jacp_[13,12] <- 0
	jacp_[14,12] <- 0
	jacp_[15,12] <- 0
	jacp_[16,12] <- 0
	jacp_[17,12] <- 0
	jacp_[18,12] <- 0
	jacp_[19,12] <- 0
	jacp_[20,12] <- 0
	jacp_[21,12] <- 0
# column 13
	jacp_[1,13] <- 0
	jacp_[2,13] <- 0
	jacp_[3,13] <- 0
	jacp_[4,13] <- 0
	jacp_[5,13] <- 0
	jacp_[6,13] <- 0
	jacp_[7,13] <- 0
	jacp_[8,13] <- -Ca*PP2B_CaM_Ca3
	jacp_[9,13] <- Ca*PP2B_CaM_Ca3
	jacp_[10,13] <- 0
	jacp_[11,13] <- 0
	jacp_[12,13] <- 0
	jacp_[13,13] <- 0
	jacp_[14,13] <- 0
	jacp_[15,13] <- 0
	jacp_[16,13] <- 0
	jacp_[17,13] <- 0
	jacp_[18,13] <- 0
	jacp_[19,13] <- 0
	jacp_[20,13] <- 0
	jacp_[21,13] <- 0
# column 14
	jacp_[1,14] <- 0
	jacp_[2,14] <- 0
	jacp_[3,14] <- 0
	jacp_[4,14] <- 0
	jacp_[5,14] <- 0
	jacp_[6,14] <- 0
	jacp_[7,14] <- 0
	jacp_[8,14] <- 0
	jacp_[9,14] <- 0
	jacp_[10,14] <- 0
	jacp_[11,14] <- 0
	jacp_[12,14] <- 0
	jacp_[13,14] <- 0
	jacp_[14,14] <- 0
	jacp_[15,14] <- 0
	jacp_[16,14] <- 0
	jacp_[17,14] <- 0
	jacp_[18,14] <- 0
	jacp_[19,14] <- 0
	jacp_[20,14] <- 0
	jacp_[21,14] <- 0
# column 15
	jacp_[1,15] <- 0
	jacp_[2,15] <- 0
	jacp_[3,15] <- 0
	jacp_[4,15] <- 0
	jacp_[5,15] <- 0
	jacp_[6,15] <- 0
	jacp_[7,15] <- 0
	jacp_[8,15] <- 0
	jacp_[9,15] <- 0
	jacp_[10,15] <- 0
	jacp_[11,15] <- 0
	jacp_[12,15] <- 0
	jacp_[13,15] <- 0
	jacp_[14,15] <- 0
	jacp_[15,15] <- 0
	jacp_[16,15] <- 0
	jacp_[17,15] <- 0
	jacp_[18,15] <- 0
	jacp_[19,15] <- 0
	jacp_[20,15] <- 0
	jacp_[21,15] <- 0
# column 16
	jacp_[1,16] <- 0
	jacp_[2,16] <- 0
	jacp_[3,16] <- 0
	jacp_[4,16] <- 0
	jacp_[5,16] <- 0
	jacp_[6,16] <- 0
	jacp_[7,16] <- 0
	jacp_[8,16] <- 0
	jacp_[9,16] <- 0
	jacp_[10,16] <- 0
	jacp_[11,16] <- 0
	jacp_[12,16] <- 0
	jacp_[13,16] <- 0
	jacp_[14,16] <- 0
	jacp_[15,16] <- 0
	jacp_[16,16] <- 0
	jacp_[17,16] <- 0
	jacp_[18,16] <- 0
	jacp_[19,16] <- 0
	jacp_[20,16] <- 0
	jacp_[21,16] <- 0
# column 17
	jacp_[1,17] <- 0
	jacp_[2,17] <- 0
	jacp_[3,17] <- 0
	jacp_[4,17] <- 0
	jacp_[5,17] <- 0
	jacp_[6,17] <- 0
	jacp_[7,17] <- 0
	jacp_[8,17] <- 0
	jacp_[9,17] <- 0
	jacp_[10,17] <- 0
	jacp_[11,17] <- 0
	jacp_[12,17] <- 0
	jacp_[13,17] <- 0
	jacp_[14,17] <- 0
	jacp_[15,17] <- 0
	jacp_[16,17] <- 0
	jacp_[17,17] <- 0
	jacp_[18,17] <- 0
	jacp_[19,17] <- 0
	jacp_[20,17] <- 0
	jacp_[21,17] <- 0
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
	jacp_[17,18] <- 0
	jacp_[18,18] <- 0
	jacp_[19,18] <- 0
	jacp_[20,18] <- 0
	jacp_[21,18] <- 0
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
	jacp_[17,19] <- 0
	jacp_[18,19] <- 0
	jacp_[19,19] <- 0
	jacp_[20,19] <- 0
	jacp_[21,19] <- 0
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
	jacp_[17,20] <- 0
	jacp_[18,20] <- 0
	jacp_[19,20] <- 0
	jacp_[20,20] <- 0
	jacp_[21,20] <- 0
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
	jacp_[17,21] <- 0
	jacp_[18,21] <- 0
	jacp_[19,21] <- 0
	jacp_[20,21] <- 0
	jacp_[21,21] <- 0
# column 22
	jacp_[1,22] <- 0
	jacp_[2,22] <- 0
	jacp_[3,22] <- 0
	jacp_[4,22] <- 0
	jacp_[5,22] <- 0
	jacp_[6,22] <- 0
	jacp_[7,22] <- 0
	jacp_[8,22] <- 0
	jacp_[9,22] <- 0
	jacp_[10,22] <- 0
	jacp_[11,22] <- 0
	jacp_[12,22] <- 0
	jacp_[13,22] <- 0
	jacp_[14,22] <- 0
	jacp_[15,22] <- 0
	jacp_[16,22] <- 0
	jacp_[17,22] <- 0
	jacp_[18,22] <- 0
	jacp_[19,22] <- 0
	jacp_[20,22] <- 0
	jacp_[21,22] <- 0
# column 23
	jacp_[1,23] <- 0
	jacp_[2,23] <- 0
	jacp_[3,23] <- 0
	jacp_[4,23] <- 0
	jacp_[5,23] <- 0
	jacp_[6,23] <- 0
	jacp_[7,23] <- 0
	jacp_[8,23] <- 0
	jacp_[9,23] <- 0
	jacp_[10,23] <- CaM*CaMKII
	jacp_[11,23] <- 0
	jacp_[12,23] <- 0
	jacp_[13,23] <- 0
	jacp_[14,23] <- 0
	jacp_[15,23] <- 0
	jacp_[16,23] <- 0
	jacp_[17,23] <- 0
	jacp_[18,23] <- 0
	jacp_[19,23] <- 0
	jacp_[20,23] <- 0
	jacp_[21,23] <- 0
# column 24
	jacp_[1,24] <- 0
	jacp_[2,24] <- 0
	jacp_[3,24] <- 0
	jacp_[4,24] <- 0
	jacp_[5,24] <- 0
	jacp_[6,24] <- 0
	jacp_[7,24] <- 0
	jacp_[8,24] <- 0
	jacp_[9,24] <- 0
	jacp_[10,24] <- 0
	jacp_[11,24] <- 0
	jacp_[12,24] <- 0
	jacp_[13,24] <- -Ca*CaMKII_CaM_Ca3
	jacp_[14,24] <- Ca*CaMKII_CaM_Ca3
	jacp_[15,24] <- 0
	jacp_[16,24] <- 0
	jacp_[17,24] <- 0
	jacp_[18,24] <- 0
	jacp_[19,24] <- 0
	jacp_[20,24] <- 0
	jacp_[21,24] <- 0
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
	jacp_[12,25] <- -Ca*CaMKII_CaM_Ca2
	jacp_[13,25] <- Ca*CaMKII_CaM_Ca2
	jacp_[14,25] <- 0
	jacp_[15,25] <- 0
	jacp_[16,25] <- 0
	jacp_[17,25] <- 0
	jacp_[18,25] <- 0
	jacp_[19,25] <- 0
	jacp_[20,25] <- 0
	jacp_[21,25] <- 0
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
	jacp_[11,26] <- -Ca*CaMKII_CaM_Ca1
	jacp_[12,26] <- Ca*CaMKII_CaM_Ca1
	jacp_[13,26] <- 0
	jacp_[14,26] <- 0
	jacp_[15,26] <- 0
	jacp_[16,26] <- 0
	jacp_[17,26] <- 0
	jacp_[18,26] <- 0
	jacp_[19,26] <- 0
	jacp_[20,26] <- 0
	jacp_[21,26] <- 0
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
	jacp_[10,27] <- -Ca*CaMKII_CaM
	jacp_[11,27] <- Ca*CaMKII_CaM
	jacp_[12,27] <- 0
	jacp_[13,27] <- 0
	jacp_[14,27] <- 0
	jacp_[15,27] <- 0
	jacp_[16,27] <- 0
	jacp_[17,27] <- 0
	jacp_[18,27] <- 0
	jacp_[19,27] <- 0
	jacp_[20,27] <- 0
	jacp_[21,27] <- 0
# column 28
	jacp_[1,28] <- -CaMKII*CaM_Ca1
	jacp_[2,28] <- 0
	jacp_[3,28] <- 0
	jacp_[4,28] <- 0
	jacp_[5,28] <- 0
	jacp_[6,28] <- 0
	jacp_[7,28] <- 0
	jacp_[8,28] <- 0
	jacp_[9,28] <- 0
	jacp_[10,28] <- 0
	jacp_[11,28] <- CaMKII*CaM_Ca1
	jacp_[12,28] <- 0
	jacp_[13,28] <- 0
	jacp_[14,28] <- 0
	jacp_[15,28] <- 0
	jacp_[16,28] <- 0
	jacp_[17,28] <- 0
	jacp_[18,28] <- 0
	jacp_[19,28] <- 0
	jacp_[20,28] <- 0
	jacp_[21,28] <- 0
# column 29
	jacp_[1,29] <- 0
	jacp_[2,29] <- -CaMKII*CaM_Ca2
	jacp_[3,29] <- 0
	jacp_[4,29] <- 0
	jacp_[5,29] <- 0
	jacp_[6,29] <- 0
	jacp_[7,29] <- 0
	jacp_[8,29] <- 0
	jacp_[9,29] <- 0
	jacp_[10,29] <- 0
	jacp_[11,29] <- 0
	jacp_[12,29] <- CaMKII*CaM_Ca2
	jacp_[13,29] <- 0
	jacp_[14,29] <- 0
	jacp_[15,29] <- 0
	jacp_[16,29] <- 0
	jacp_[17,29] <- 0
	jacp_[18,29] <- 0
	jacp_[19,29] <- 0
	jacp_[20,29] <- 0
	jacp_[21,29] <- 0
# column 30
	jacp_[1,30] <- 0
	jacp_[2,30] <- 0
	jacp_[3,30] <- -CaMKII*CaM_Ca3
	jacp_[4,30] <- 0
	jacp_[5,30] <- 0
	jacp_[6,30] <- 0
	jacp_[7,30] <- 0
	jacp_[8,30] <- 0
	jacp_[9,30] <- 0
	jacp_[10,30] <- 0
	jacp_[11,30] <- 0
	jacp_[12,30] <- 0
	jacp_[13,30] <- CaMKII*CaM_Ca3
	jacp_[14,30] <- 0
	jacp_[15,30] <- 0
	jacp_[16,30] <- 0
	jacp_[17,30] <- 0
	jacp_[18,30] <- 0
	jacp_[19,30] <- 0
	jacp_[20,30] <- 0
	jacp_[21,30] <- 0
# column 31
	jacp_[1,31] <- 0
	jacp_[2,31] <- 0
	jacp_[3,31] <- 0
	jacp_[4,31] <- -CaMKII*CaM_Ca4
	jacp_[5,31] <- 0
	jacp_[6,31] <- 0
	jacp_[7,31] <- 0
	jacp_[8,31] <- 0
	jacp_[9,31] <- 0
	jacp_[10,31] <- 0
	jacp_[11,31] <- 0
	jacp_[12,31] <- 0
	jacp_[13,31] <- 0
	jacp_[14,31] <- CaMKII*CaM_Ca4
	jacp_[15,31] <- 0
	jacp_[16,31] <- 0
	jacp_[17,31] <- 0
	jacp_[18,31] <- 0
	jacp_[19,31] <- 0
	jacp_[20,31] <- 0
	jacp_[21,31] <- 0
# column 32
	jacp_[1,32] <- 0
	jacp_[2,32] <- 0
	jacp_[3,32] <- 0
	jacp_[4,32] <- 0
	jacp_[5,32] <- 0
	jacp_[6,32] <- 0
	jacp_[7,32] <- 0
	jacp_[8,32] <- 0
	jacp_[9,32] <- 0
	jacp_[10,32] <- 0
	jacp_[11,32] <- 0
	jacp_[12,32] <- 0
	jacp_[13,32] <- 0
	jacp_[14,32] <- 0
	jacp_[15,32] <- 0
	jacp_[16,32] <- 0
	jacp_[17,32] <- 0
	jacp_[18,32] <- 0
	jacp_[19,32] <- 0
	jacp_[20,32] <- 0
	jacp_[21,32] <- 0
# column 33
	jacp_[1,33] <- 0
	jacp_[2,33] <- 0
	jacp_[3,33] <- 0
	jacp_[4,33] <- 0
	jacp_[5,33] <- 0
	jacp_[6,33] <- 0
	jacp_[7,33] <- 0
	jacp_[8,33] <- 0
	jacp_[9,33] <- 0
	jacp_[10,33] <- 0
	jacp_[11,33] <- 0
	jacp_[12,33] <- 0
	jacp_[13,33] <- 0
	jacp_[14,33] <- 0
	jacp_[15,33] <- 0
	jacp_[16,33] <- 0
	jacp_[17,33] <- 0
	jacp_[18,33] <- 0
	jacp_[19,33] <- 0
	jacp_[20,33] <- 0
	jacp_[21,33] <- 0
# column 34
	jacp_[1,34] <- 0
	jacp_[2,34] <- 0
	jacp_[3,34] <- 0
	jacp_[4,34] <- 0
	jacp_[5,34] <- 0
	jacp_[6,34] <- 0
	jacp_[7,34] <- 0
	jacp_[8,34] <- 0
	jacp_[9,34] <- 0
	jacp_[10,34] <- 0
	jacp_[11,34] <- 0
	jacp_[12,34] <- 0
	jacp_[13,34] <- 0
	jacp_[14,34] <- 0
	jacp_[15,34] <- 0
	jacp_[16,34] <- 0
	jacp_[17,34] <- 0
	jacp_[18,34] <- 0
	jacp_[19,34] <- 0
	jacp_[20,34] <- 0
	jacp_[21,34] <- 0
# column 35
	jacp_[1,35] <- 0
	jacp_[2,35] <- 0
	jacp_[3,35] <- 0
	jacp_[4,35] <- 0
	jacp_[5,35] <- 0
	jacp_[6,35] <- 0
	jacp_[7,35] <- 0
	jacp_[8,35] <- 0
	jacp_[9,35] <- 0
	jacp_[10,35] <- 0
	jacp_[11,35] <- 0
	jacp_[12,35] <- 0
	jacp_[13,35] <- 0
	jacp_[14,35] <- 0
	jacp_[15,35] <- 0
	jacp_[16,35] <- 0
	jacp_[17,35] <- 0
	jacp_[18,35] <- 0
	jacp_[19,35] <- 0
	jacp_[20,35] <- 0
	jacp_[21,35] <- 0
# column 36
	jacp_[1,36] <- 0
	jacp_[2,36] <- 0
	jacp_[3,36] <- 0
	jacp_[4,36] <- 0
	jacp_[5,36] <- 0
	jacp_[6,36] <- 0
	jacp_[7,36] <- 0
	jacp_[8,36] <- 0
	jacp_[9,36] <- 0
	jacp_[10,36] <- 0
	jacp_[11,36] <- 0
	jacp_[12,36] <- 0
	jacp_[13,36] <- 0
	jacp_[14,36] <- 0
	jacp_[15,36] <- 0
	jacp_[16,36] <- 0
	jacp_[17,36] <- 0
	jacp_[18,36] <- 0
	jacp_[19,36] <- 0
	jacp_[20,36] <- 0
	jacp_[21,36] <- 0
# column 37
	jacp_[1,37] <- 0
	jacp_[2,37] <- 0
	jacp_[3,37] <- 0
	jacp_[4,37] <- 0
	jacp_[5,37] <- 0
	jacp_[6,37] <- 0
	jacp_[7,37] <- 0
	jacp_[8,37] <- 0
	jacp_[9,37] <- 0
	jacp_[10,37] <- 0
	jacp_[11,37] <- 0
	jacp_[12,37] <- 0
	jacp_[13,37] <- 0
	jacp_[14,37] <- 0
	jacp_[15,37] <- Ca*pCaMKII_CaM_Ca3
	jacp_[16,37] <- 0
	jacp_[17,37] <- -Ca*pCaMKII_CaM_Ca3
	jacp_[18,37] <- 0
	jacp_[19,37] <- 0
	jacp_[20,37] <- 0
	jacp_[21,37] <- 0
# column 38
	jacp_[1,38] <- 0
	jacp_[2,38] <- 0
	jacp_[3,38] <- 0
	jacp_[4,38] <- 0
	jacp_[5,38] <- 0
	jacp_[6,38] <- 0
	jacp_[7,38] <- 0
	jacp_[8,38] <- 0
	jacp_[9,38] <- 0
	jacp_[10,38] <- 0
	jacp_[11,38] <- 0
	jacp_[12,38] <- 0
	jacp_[13,38] <- 0
	jacp_[14,38] <- 0
	jacp_[15,38] <- 0
	jacp_[16,38] <- -CaM*pCaMKIIa
	jacp_[17,38] <- 0
	jacp_[18,38] <- 0
	jacp_[19,38] <- 0
	jacp_[20,38] <- CaM*pCaMKIIa
	jacp_[21,38] <- 0
# column 39
	jacp_[1,39] <- -CaM_Ca1*pCaMKIIa
	jacp_[2,39] <- 0
	jacp_[3,39] <- 0
	jacp_[4,39] <- 0
	jacp_[5,39] <- 0
	jacp_[6,39] <- 0
	jacp_[7,39] <- 0
	jacp_[8,39] <- 0
	jacp_[9,39] <- 0
	jacp_[10,39] <- 0
	jacp_[11,39] <- 0
	jacp_[12,39] <- 0
	jacp_[13,39] <- 0
	jacp_[14,39] <- 0
	jacp_[15,39] <- 0
	jacp_[16,39] <- -CaM_Ca1*pCaMKIIa
	jacp_[17,39] <- 0
	jacp_[18,39] <- 0
	jacp_[19,39] <- CaM_Ca1*pCaMKIIa
	jacp_[20,39] <- 0
	jacp_[21,39] <- 0
# column 40
	jacp_[1,40] <- 0
	jacp_[2,40] <- -CaM_Ca2*pCaMKIIa
	jacp_[3,40] <- 0
	jacp_[4,40] <- 0
	jacp_[5,40] <- 0
	jacp_[6,40] <- 0
	jacp_[7,40] <- 0
	jacp_[8,40] <- 0
	jacp_[9,40] <- 0
	jacp_[10,40] <- 0
	jacp_[11,40] <- 0
	jacp_[12,40] <- 0
	jacp_[13,40] <- 0
	jacp_[14,40] <- 0
	jacp_[15,40] <- 0
	jacp_[16,40] <- -CaM_Ca2*pCaMKIIa
	jacp_[17,40] <- 0
	jacp_[18,40] <- CaM_Ca2*pCaMKIIa
	jacp_[19,40] <- 0
	jacp_[20,40] <- 0
	jacp_[21,40] <- 0
# column 41
	jacp_[1,41] <- 0
	jacp_[2,41] <- 0
	jacp_[3,41] <- -CaM_Ca3*pCaMKIIa
	jacp_[4,41] <- 0
	jacp_[5,41] <- 0
	jacp_[6,41] <- 0
	jacp_[7,41] <- 0
	jacp_[8,41] <- 0
	jacp_[9,41] <- 0
	jacp_[10,41] <- 0
	jacp_[11,41] <- 0
	jacp_[12,41] <- 0
	jacp_[13,41] <- 0
	jacp_[14,41] <- 0
	jacp_[15,41] <- 0
	jacp_[16,41] <- -CaM_Ca3*pCaMKIIa
	jacp_[17,41] <- CaM_Ca3*pCaMKIIa
	jacp_[18,41] <- 0
	jacp_[19,41] <- 0
	jacp_[20,41] <- 0
	jacp_[21,41] <- 0
# column 42
	jacp_[1,42] <- 0
	jacp_[2,42] <- 0
	jacp_[3,42] <- 0
	jacp_[4,42] <- 0
	jacp_[5,42] <- 0
	jacp_[6,42] <- 0
	jacp_[7,42] <- 0
	jacp_[8,42] <- 0
	jacp_[9,42] <- 0
	jacp_[10,42] <- 0
	jacp_[11,42] <- 0
	jacp_[12,42] <- 0
	jacp_[13,42] <- 0
	jacp_[14,42] <- 0
	jacp_[15,42] <- 0
	jacp_[16,42] <- 0
	jacp_[17,42] <- Ca*pCaMKII_CaM_Ca2
	jacp_[18,42] <- -Ca*pCaMKII_CaM_Ca2
	jacp_[19,42] <- 0
	jacp_[20,42] <- 0
	jacp_[21,42] <- 0
# column 43
	jacp_[1,43] <- 0
	jacp_[2,43] <- 0
	jacp_[3,43] <- 0
	jacp_[4,43] <- 0
	jacp_[5,43] <- 0
	jacp_[6,43] <- 0
	jacp_[7,43] <- 0
	jacp_[8,43] <- 0
	jacp_[9,43] <- 0
	jacp_[10,43] <- 0
	jacp_[11,43] <- 0
	jacp_[12,43] <- 0
	jacp_[13,43] <- 0
	jacp_[14,43] <- 0
	jacp_[15,43] <- 0
	jacp_[16,43] <- 0
	jacp_[17,43] <- 0
	jacp_[18,43] <- Ca*pCaMKII_CaM_Ca1
	jacp_[19,43] <- -Ca*pCaMKII_CaM_Ca1
	jacp_[20,43] <- 0
	jacp_[21,43] <- 0
# column 44
	jacp_[1,44] <- 0
	jacp_[2,44] <- 0
	jacp_[3,44] <- 0
	jacp_[4,44] <- -CaM_Ca4*pCaMKIIa
	jacp_[5,44] <- 0
	jacp_[6,44] <- 0
	jacp_[7,44] <- 0
	jacp_[8,44] <- 0
	jacp_[9,44] <- 0
	jacp_[10,44] <- 0
	jacp_[11,44] <- 0
	jacp_[12,44] <- 0
	jacp_[13,44] <- 0
	jacp_[14,44] <- 0
	jacp_[15,44] <- CaM_Ca4*pCaMKIIa
	jacp_[16,44] <- -CaM_Ca4*pCaMKIIa
	jacp_[17,44] <- 0
	jacp_[18,44] <- 0
	jacp_[19,44] <- 0
	jacp_[20,44] <- 0
	jacp_[21,44] <- 0
# column 45
	jacp_[1,45] <- 0
	jacp_[2,45] <- 0
	jacp_[3,45] <- 0
	jacp_[4,45] <- 0
	jacp_[5,45] <- 0
	jacp_[6,45] <- 0
	jacp_[7,45] <- 0
	jacp_[8,45] <- 0
	jacp_[9,45] <- 0
	jacp_[10,45] <- 0
	jacp_[11,45] <- 0
	jacp_[12,45] <- 0
	jacp_[13,45] <- 0
	jacp_[14,45] <- 0
	jacp_[15,45] <- 0
	jacp_[16,45] <- 0
	jacp_[17,45] <- 0
	jacp_[18,45] <- 0
	jacp_[19,45] <- Ca*pCaMKII_CaM
	jacp_[20,45] <- -Ca*pCaMKII_CaM
	jacp_[21,45] <- 0
# column 46
	jacp_[1,46] <- 0
	jacp_[2,46] <- 0
	jacp_[3,46] <- 0
	jacp_[4,46] <- 0
	jacp_[5,46] <- 0
	jacp_[6,46] <- 0
	jacp_[7,46] <- 0
	jacp_[8,46] <- 0
	jacp_[9,46] <- 0
	jacp_[10,46] <- 0
	jacp_[11,46] <- 0
	jacp_[12,46] <- 0
	jacp_[13,46] <- 0
	jacp_[14,46] <- 0
	jacp_[15,46] <- 0
	jacp_[16,46] <- 0
	jacp_[17,46] <- 0
	jacp_[18,46] <- 0
	jacp_[19,46] <- 0
	jacp_[20,46] <- 0
	jacp_[21,46] <- 0
# column 47
	jacp_[1,47] <- 0
	jacp_[2,47] <- 0
	jacp_[3,47] <- 0
	jacp_[4,47] <- 0
	jacp_[5,47] <- 0
	jacp_[6,47] <- 0
	jacp_[7,47] <- 0
	jacp_[8,47] <- 0
	jacp_[9,47] <- 0
	jacp_[10,47] <- 0
	jacp_[11,47] <- 0
	jacp_[12,47] <- 0
	jacp_[13,47] <- 0
	jacp_[14,47] <- 0
	jacp_[15,47] <- 0
	jacp_[16,47] <- 0
	jacp_[17,47] <- 0
	jacp_[18,47] <- 0
	jacp_[19,47] <- 0
	jacp_[20,47] <- 0
	jacp_[21,47] <- 0
# column 48
	jacp_[1,48] <- 0
	jacp_[2,48] <- 0
	jacp_[3,48] <- 0
	jacp_[4,48] <- 0
	jacp_[5,48] <- 0
	jacp_[6,48] <- 0
	jacp_[7,48] <- 0
	jacp_[8,48] <- 0
	jacp_[9,48] <- 0
	jacp_[10,48] <- 0
	jacp_[11,48] <- 0
	jacp_[12,48] <- 0
	jacp_[13,48] <- 0
	jacp_[14,48] <- 0
	jacp_[15,48] <- 0
	jacp_[16,48] <- 0
	jacp_[17,48] <- 0
	jacp_[18,48] <- 0
	jacp_[19,48] <- 0
	jacp_[20,48] <- 0
	jacp_[21,48] <- 0
# column 49
	jacp_[1,49] <- 0
	jacp_[2,49] <- 0
	jacp_[3,49] <- 0
	jacp_[4,49] <- 0
	jacp_[5,49] <- 0
	jacp_[6,49] <- 0
	jacp_[7,49] <- 0
	jacp_[8,49] <- 0
	jacp_[9,49] <- 0
	jacp_[10,49] <- 0
	jacp_[11,49] <- 0
	jacp_[12,49] <- 0
	jacp_[13,49] <- 0
	jacp_[14,49] <- 0
	jacp_[15,49] <- 0
	jacp_[16,49] <- 0
	jacp_[17,49] <- 0
	jacp_[18,49] <- 0
	jacp_[19,49] <- 0
	jacp_[20,49] <- 0
	jacp_[21,49] <- 0
# column 50
	jacp_[1,50] <- 0
	jacp_[2,50] <- 0
	jacp_[3,50] <- 0
	jacp_[4,50] <- 0
	jacp_[5,50] <- 0
	jacp_[6,50] <- 0
	jacp_[7,50] <- 0
	jacp_[8,50] <- 0
	jacp_[9,50] <- 0
	jacp_[10,50] <- 0
	jacp_[11,50] <- 0
	jacp_[12,50] <- 0
	jacp_[13,50] <- 0
	jacp_[14,50] <- 0
	jacp_[15,50] <- 0
	jacp_[16,50] <- 0
	jacp_[17,50] <- 0
	jacp_[18,50] <- 0
	jacp_[19,50] <- 0
	jacp_[20,50] <- 0
	jacp_[21,50] <- 0
# column 51
	jacp_[1,51] <- 0
	jacp_[2,51] <- 0
	jacp_[3,51] <- 0
	jacp_[4,51] <- 0
	jacp_[5,51] <- 0
	jacp_[6,51] <- 0
	jacp_[7,51] <- 0
	jacp_[8,51] <- 0
	jacp_[9,51] <- 0
	jacp_[10,51] <- 0
	jacp_[11,51] <- 0
	jacp_[12,51] <- 0
	jacp_[13,51] <- 0
	jacp_[14,51] <- -CaMKII_CaM_Ca4*pairedCaMKIIc
	jacp_[15,51] <- CaMKII_CaM_Ca4*pairedCaMKIIc
	jacp_[16,51] <- 0
	jacp_[17,51] <- 0
	jacp_[18,51] <- 0
	jacp_[19,51] <- 0
	jacp_[20,51] <- 0
	jacp_[21,51] <- 0
# column 52
	jacp_[1,52] <- 0
	jacp_[2,52] <- 0
	jacp_[3,52] <- 0
	jacp_[4,52] <- 0
	jacp_[5,52] <- 0
	jacp_[6,52] <- 0
	jacp_[7,52] <- 0
	jacp_[8,52] <- 0
	jacp_[9,52] <- 0
	jacp_[10,52] <- 0
	jacp_[11,52] <- 0
	jacp_[12,52] <- 0
	jacp_[13,52] <- 0
	jacp_[14,52] <- 0
	jacp_[15,52] <- 0
	jacp_[16,52] <- -PP1*pCaMKIIa
	jacp_[17,52] <- 0
	jacp_[18,52] <- 0
	jacp_[19,52] <- 0
	jacp_[20,52] <- 0
	jacp_[21,52] <- PP1*pCaMKIIa
# column 53
	jacp_[1,53] <- 0
	jacp_[2,53] <- 0
	jacp_[3,53] <- 0
	jacp_[4,53] <- 0
	jacp_[5,53] <- 0
	jacp_[6,53] <- 0
	jacp_[7,53] <- 0
	jacp_[8,53] <- 0
	jacp_[9,53] <- 0
	jacp_[10,53] <- 0
	jacp_[11,53] <- 0
	jacp_[12,53] <- 0
	jacp_[13,53] <- 0
	jacp_[14,53] <- 0
	jacp_[15,53] <- 0
	jacp_[16,53] <- PP1__pCaMKIIa
	jacp_[17,53] <- 0
	jacp_[18,53] <- 0
	jacp_[19,53] <- 0
	jacp_[20,53] <- 0
	jacp_[21,53] <- -PP1__pCaMKIIa
# column 54
	jacp_[1,54] <- 0
	jacp_[2,54] <- 0
	jacp_[3,54] <- 0
	jacp_[4,54] <- 0
	jacp_[5,54] <- 0
	jacp_[6,54] <- 0
	jacp_[7,54] <- 0
	jacp_[8,54] <- 0
	jacp_[9,54] <- 0
	jacp_[10,54] <- 0
	jacp_[11,54] <- 0
	jacp_[12,54] <- 0
	jacp_[13,54] <- 0
	jacp_[14,54] <- 0
	jacp_[15,54] <- 0
	jacp_[16,54] <- 0
	jacp_[17,54] <- 0
	jacp_[18,54] <- 0
	jacp_[19,54] <- 0
	jacp_[20,54] <- 0
	jacp_[21,54] <- -PP1__pCaMKIIa
# column 55
	jacp_[1,55] <- 0
	jacp_[2,55] <- 0
	jacp_[3,55] <- 0
	jacp_[4,55] <- 0
	jacp_[5,55] <- 0
	jacp_[6,55] <- 0
	jacp_[7,55] <- 0
	jacp_[8,55] <- 0
	jacp_[9,55] <- 0
	jacp_[10,55] <- 0
	jacp_[11,55] <- 0
	jacp_[12,55] <- 0
	jacp_[13,55] <- 0
	jacp_[14,55] <- 0
	jacp_[15,55] <- 0
	jacp_[16,55] <- 0
	jacp_[17,55] <- 0
	jacp_[18,55] <- 0
	jacp_[19,55] <- 0
	jacp_[20,55] <- 0
	jacp_[21,55] <- 0
# column 56
	jacp_[1,56] <- 0
	jacp_[2,56] <- 0
	jacp_[3,56] <- 0
	jacp_[4,56] <- 0
	jacp_[5,56] <- 0
	jacp_[6,56] <- 0
	jacp_[7,56] <- 0
	jacp_[8,56] <- 0
	jacp_[9,56] <- 0
	jacp_[10,56] <- 0
	jacp_[11,56] <- 0
	jacp_[12,56] <- 0
	jacp_[13,56] <- 0
	jacp_[14,56] <- 0
	jacp_[15,56] <- 0
	jacp_[16,56] <- 0
	jacp_[17,56] <- 0
	jacp_[18,56] <- 0
	jacp_[19,56] <- 0
	jacp_[20,56] <- 0
	jacp_[21,56] <- 0
# column 57
	jacp_[1,57] <- 0
	jacp_[2,57] <- 0
	jacp_[3,57] <- 0
	jacp_[4,57] <- 0
	jacp_[5,57] <- 0
	jacp_[6,57] <- 0
	jacp_[7,57] <- 0
	jacp_[8,57] <- 0
	jacp_[9,57] <- 0
	jacp_[10,57] <- 0
	jacp_[11,57] <- 0
	jacp_[12,57] <- 0
	jacp_[13,57] <- 0
	jacp_[14,57] <- 0
	jacp_[15,57] <- 0
	jacp_[16,57] <- 0
	jacp_[17,57] <- 0
	jacp_[18,57] <- 0
	jacp_[19,57] <- 0
	jacp_[20,57] <- 0
	jacp_[21,57] <- 0
# column 58
	jacp_[1,58] <- 0
	jacp_[2,58] <- 0
	jacp_[3,58] <- 0
	jacp_[4,58] <- 0
	jacp_[5,58] <- 0
	jacp_[6,58] <- 0
	jacp_[7,58] <- 0
	jacp_[8,58] <- 0
	jacp_[9,58] <- 0
	jacp_[10,58] <- 0
	jacp_[11,58] <- 0
	jacp_[12,58] <- 0
	jacp_[13,58] <- 0
	jacp_[14,58] <- 0
	jacp_[15,58] <- 0
	jacp_[16,58] <- 0
	jacp_[17,58] <- 0
	jacp_[18,58] <- 0
	jacp_[19,58] <- 0
	jacp_[20,58] <- 0
	jacp_[21,58] <- 0
# column 59
	jacp_[1,59] <- 0
	jacp_[2,59] <- 0
	jacp_[3,59] <- 0
	jacp_[4,59] <- 0
	jacp_[5,59] <- 0
	jacp_[6,59] <- 0
	jacp_[7,59] <- 0
	jacp_[8,59] <- 0
	jacp_[9,59] <- 0
	jacp_[10,59] <- 0
	jacp_[11,59] <- 0
	jacp_[12,59] <- 0
	jacp_[13,59] <- 0
	jacp_[14,59] <- 0
	jacp_[15,59] <- 0
	jacp_[16,59] <- 0
	jacp_[17,59] <- 0
	jacp_[18,59] <- 0
	jacp_[19,59] <- 0
	jacp_[20,59] <- 0
	jacp_[21,59] <- 0
# column 60
	jacp_[1,60] <- 0
	jacp_[2,60] <- 0
	jacp_[3,60] <- 0
	jacp_[4,60] <- 0
	jacp_[5,60] <- 0
	jacp_[6,60] <- 0
	jacp_[7,60] <- 0
	jacp_[8,60] <- 0
	jacp_[9,60] <- 0
	jacp_[10,60] <- 0
	jacp_[11,60] <- 0
	jacp_[12,60] <- 0
	jacp_[13,60] <- 0
	jacp_[14,60] <- 0
	jacp_[15,60] <- 0
	jacp_[16,60] <- 0
	jacp_[17,60] <- 0
	jacp_[18,60] <- 0
	jacp_[19,60] <- 0
	jacp_[20,60] <- 0
	jacp_[21,60] <- 0
	rownames(jacp_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
	colnames(jacp_) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "CaBase", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0", "isOn")
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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
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
	Ca <- isOn*CaBase
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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
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
	Ca<-isOn*CaBase
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
	fjac_ <- matrix(NA,5,21)
# column 1
	fjac_[1,1] <- 0
	fjac_[2,1] <- 0
	fjac_[3,1] <- 0
	fjac_[4,1] <- 0
	fjac_[5,1] <- 0
# column 2
	fjac_[1,2] <- 0
	fjac_[2,2] <- 0
	fjac_[3,2] <- 0
	fjac_[4,2] <- 0
	fjac_[5,2] <- 0
# column 3
	fjac_[1,3] <- 0
	fjac_[2,3] <- 0
	fjac_[3,3] <- 0
	fjac_[4,3] <- 0
	fjac_[5,3] <- 0
# column 4
	fjac_[1,4] <- 0
	fjac_[2,4] <- 0
	fjac_[3,4] <- 0
	fjac_[4,4] <- 0
	fjac_[5,4] <- 0
# column 5
	fjac_[1,5] <- 0
	fjac_[2,5] <- 0
	fjac_[3,5] <- 1/(s+PP2B_0)
	fjac_[4,5] <- 0
	fjac_[5,5] <- 0
# column 6
	fjac_[1,6] <- 0
	fjac_[2,6] <- 0
	fjac_[3,6] <- 1/(s+PP2B_0)
	fjac_[4,6] <- 0
	fjac_[5,6] <- 0
# column 7
	fjac_[1,7] <- 0
	fjac_[2,7] <- 0
	fjac_[3,7] <- 1/(s+PP2B_0)
	fjac_[4,7] <- 0
	fjac_[5,7] <- 0
# column 8
	fjac_[1,8] <- 0
	fjac_[2,8] <- 0
	fjac_[3,8] <- 1/(s+PP2B_0)
	fjac_[4,8] <- 0
	fjac_[5,8] <- 0
# column 9
	fjac_[1,9] <- 0
	fjac_[2,9] <- 0
	fjac_[3,9] <- 1/(s+PP2B_0)
	fjac_[4,9] <- 100/(s+PP2B_0)
	fjac_[5,9] <- 0
# column 10
	fjac_[1,10] <- 0
	fjac_[2,10] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,10] <- 0
	fjac_[4,10] <- 0
	fjac_[5,10] <- 0
# column 11
	fjac_[1,11] <- 0
	fjac_[2,11] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,11] <- 0
	fjac_[4,11] <- 0
	fjac_[5,11] <- 0
# column 12
	fjac_[1,12] <- 0
	fjac_[2,12] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,12] <- 0
	fjac_[4,12] <- 0
	fjac_[5,12] <- 0
# column 13
	fjac_[1,13] <- 0
	fjac_[2,13] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,13] <- 0
	fjac_[4,13] <- 0
	fjac_[5,13] <- 0
# column 14
	fjac_[1,14] <- 0
	fjac_[2,14] <- -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,14] <- 0
	fjac_[4,14] <- 0
	fjac_[5,14] <- 0
# column 15
	fjac_[1,15] <- 0
	fjac_[2,15] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,15] <- 0
	fjac_[4,15] <- 0
	fjac_[5,15] <- 0
# column 16
	fjac_[1,16] <- 0
	fjac_[2,16] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,16] <- 0
	fjac_[4,16] <- 0
	fjac_[5,16] <- 0
# column 17
	fjac_[1,17] <- 0
	fjac_[2,17] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,17] <- 0
	fjac_[4,17] <- 0
	fjac_[5,17] <- 0
# column 18
	fjac_[1,18] <- 0
	fjac_[2,18] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,18] <- 0
	fjac_[4,18] <- 0
	fjac_[5,18] <- 0
# column 19
	fjac_[1,19] <- 0
	fjac_[2,19] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,19] <- 0
	fjac_[4,19] <- 0
	fjac_[5,19] <- 0
# column 20
	fjac_[1,20] <- 0
	fjac_[2,20] <- 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)^2
	fjac_[3,20] <- 0
	fjac_[4,20] <- 0
	fjac_[5,20] <- 0
# column 21
	fjac_[1,21] <- 0
	fjac_[2,21] <- 0
	fjac_[3,21] <- 0
	fjac_[4,21] <- 0
	fjac_[5,21] <- 0
	rownames(fjac_) <- c("CaPerCaM", "AutoCaMKII", "CaMPerPP2B", "ActivePP2B", "Camonitor")
	colnames(fjac_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
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
	Ca<-isOn*CaBase
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
	fjacp_ <- matrix(NA,5,60)
# column 1
	fjacp_[1,1] <- 0
	fjacp_[2,1] <- 0
	fjacp_[3,1] <- 0
	fjacp_[4,1] <- 0
	fjacp_[5,1] <- 0
# column 2
	fjacp_[1,2] <- 0
	fjacp_[2,2] <- 0
	fjacp_[3,2] <- 0
	fjacp_[4,2] <- 0
	fjacp_[5,2] <- 0
# column 3
	fjacp_[1,3] <- 0
	fjacp_[2,3] <- 0
	fjacp_[3,3] <- 0
	fjacp_[4,3] <- 0
	fjacp_[5,3] <- 0
# column 4
	fjacp_[1,4] <- 0
	fjacp_[2,4] <- 0
	fjacp_[3,4] <- 0
	fjacp_[4,4] <- 0
	fjacp_[5,4] <- 0
# column 5
	fjacp_[1,5] <- 0
	fjacp_[2,5] <- 0
	fjacp_[3,5] <- 0
	fjacp_[4,5] <- 0
	fjacp_[5,5] <- 0
# column 6
	fjacp_[1,6] <- 0
	fjacp_[2,6] <- 0
	fjacp_[3,6] <- 0
	fjacp_[4,6] <- 0
	fjacp_[5,6] <- 0
# column 7
	fjacp_[1,7] <- 0
	fjacp_[2,7] <- 0
	fjacp_[3,7] <- 0
	fjacp_[4,7] <- 0
	fjacp_[5,7] <- 0
# column 8
	fjacp_[1,8] <- 0
	fjacp_[2,8] <- 0
	fjacp_[3,8] <- 0
	fjacp_[4,8] <- 0
	fjacp_[5,8] <- 0
# column 9
	fjacp_[1,9] <- 0
	fjacp_[2,9] <- 0
	fjacp_[3,9] <- 0
	fjacp_[4,9] <- 0
	fjacp_[5,9] <- 0
# column 10
	fjacp_[1,10] <- 0
	fjacp_[2,10] <- 0
	fjacp_[3,10] <- 0
	fjacp_[4,10] <- 0
	fjacp_[5,10] <- 0
# column 11
	fjacp_[1,11] <- 0
	fjacp_[2,11] <- 0
	fjacp_[3,11] <- 0
	fjacp_[4,11] <- 0
	fjacp_[5,11] <- 0
# column 12
	fjacp_[1,12] <- 0
	fjacp_[2,12] <- 0
	fjacp_[3,12] <- 0
	fjacp_[4,12] <- 0
	fjacp_[5,12] <- 0
# column 13
	fjacp_[1,13] <- 0
	fjacp_[2,13] <- 0
	fjacp_[3,13] <- 0
	fjacp_[4,13] <- 0
	fjacp_[5,13] <- 0
# column 14
	fjacp_[1,14] <- 0
	fjacp_[2,14] <- 0
	fjacp_[3,14] <- 0
	fjacp_[4,14] <- 0
	fjacp_[5,14] <- 0
# column 15
	fjacp_[1,15] <- 0
	fjacp_[2,15] <- 0
	fjacp_[3,15] <- 0
	fjacp_[4,15] <- 0
	fjacp_[5,15] <- 0
# column 16
	fjacp_[1,16] <- 0
	fjacp_[2,16] <- 0
	fjacp_[3,16] <- 0
	fjacp_[4,16] <- 0
	fjacp_[5,16] <- 0
# column 17
	fjacp_[1,17] <- 0
	fjacp_[2,17] <- 0
	fjacp_[3,17] <- 0
	fjacp_[4,17] <- 0
	fjacp_[5,17] <- 0
# column 18
	fjacp_[1,18] <- 0
	fjacp_[2,18] <- 0
	fjacp_[3,18] <- 0
	fjacp_[4,18] <- 0
	fjacp_[5,18] <- 0
# column 19
	fjacp_[1,19] <- 0
	fjacp_[2,19] <- 0
	fjacp_[3,19] <- 0
	fjacp_[4,19] <- 0
	fjacp_[5,19] <- 0
# column 20
	fjacp_[1,20] <- 0
	fjacp_[2,20] <- 0
	fjacp_[3,20] <- 0
	fjacp_[4,20] <- 0
	fjacp_[5,20] <- 0
# column 21
	fjacp_[1,21] <- 0
	fjacp_[2,21] <- 0
	fjacp_[3,21] <- 0
	fjacp_[4,21] <- 0
	fjacp_[5,21] <- 0
# column 22
	fjacp_[1,22] <- 0
	fjacp_[2,22] <- 0
	fjacp_[3,22] <- 0
	fjacp_[4,22] <- 0
	fjacp_[5,22] <- 0
# column 23
	fjacp_[1,23] <- 0
	fjacp_[2,23] <- 0
	fjacp_[3,23] <- 0
	fjacp_[4,23] <- 0
	fjacp_[5,23] <- 0
# column 24
	fjacp_[1,24] <- 0
	fjacp_[2,24] <- 0
	fjacp_[3,24] <- 0
	fjacp_[4,24] <- 0
	fjacp_[5,24] <- 0
# column 25
	fjacp_[1,25] <- 0
	fjacp_[2,25] <- 0
	fjacp_[3,25] <- 0
	fjacp_[4,25] <- 0
	fjacp_[5,25] <- 0
# column 26
	fjacp_[1,26] <- 0
	fjacp_[2,26] <- 0
	fjacp_[3,26] <- 0
	fjacp_[4,26] <- 0
	fjacp_[5,26] <- 0
# column 27
	fjacp_[1,27] <- 0
	fjacp_[2,27] <- 0
	fjacp_[3,27] <- 0
	fjacp_[4,27] <- 0
	fjacp_[5,27] <- 0
# column 28
	fjacp_[1,28] <- 0
	fjacp_[2,28] <- 0
	fjacp_[3,28] <- 0
	fjacp_[4,28] <- 0
	fjacp_[5,28] <- 0
# column 29
	fjacp_[1,29] <- 0
	fjacp_[2,29] <- 0
	fjacp_[3,29] <- 0
	fjacp_[4,29] <- 0
	fjacp_[5,29] <- 0
# column 30
	fjacp_[1,30] <- 0
	fjacp_[2,30] <- 0
	fjacp_[3,30] <- 0
	fjacp_[4,30] <- 0
	fjacp_[5,30] <- 0
# column 31
	fjacp_[1,31] <- 0
	fjacp_[2,31] <- 0
	fjacp_[3,31] <- 0
	fjacp_[4,31] <- 0
	fjacp_[5,31] <- 0
# column 32
	fjacp_[1,32] <- 0
	fjacp_[2,32] <- 0
	fjacp_[3,32] <- 0
	fjacp_[4,32] <- 0
	fjacp_[5,32] <- 0
# column 33
	fjacp_[1,33] <- 0
	fjacp_[2,33] <- 0
	fjacp_[3,33] <- 0
	fjacp_[4,33] <- 0
	fjacp_[5,33] <- 0
# column 34
	fjacp_[1,34] <- 0
	fjacp_[2,34] <- 0
	fjacp_[3,34] <- 0
	fjacp_[4,34] <- 0
	fjacp_[5,34] <- 0
# column 35
	fjacp_[1,35] <- 0
	fjacp_[2,35] <- 0
	fjacp_[3,35] <- 0
	fjacp_[4,35] <- 0
	fjacp_[5,35] <- 0
# column 36
	fjacp_[1,36] <- 0
	fjacp_[2,36] <- 0
	fjacp_[3,36] <- 0
	fjacp_[4,36] <- 0
	fjacp_[5,36] <- 0
# column 37
	fjacp_[1,37] <- 0
	fjacp_[2,37] <- 0
	fjacp_[3,37] <- 0
	fjacp_[4,37] <- 0
	fjacp_[5,37] <- 0
# column 38
	fjacp_[1,38] <- 0
	fjacp_[2,38] <- 0
	fjacp_[3,38] <- 0
	fjacp_[4,38] <- 0
	fjacp_[5,38] <- 0
# column 39
	fjacp_[1,39] <- 0
	fjacp_[2,39] <- 0
	fjacp_[3,39] <- 0
	fjacp_[4,39] <- 0
	fjacp_[5,39] <- 0
# column 40
	fjacp_[1,40] <- 0
	fjacp_[2,40] <- 0
	fjacp_[3,40] <- 0
	fjacp_[4,40] <- 0
	fjacp_[5,40] <- 0
# column 41
	fjacp_[1,41] <- 0
	fjacp_[2,41] <- 0
	fjacp_[3,41] <- 0
	fjacp_[4,41] <- 0
	fjacp_[5,41] <- 0
# column 42
	fjacp_[1,42] <- 0
	fjacp_[2,42] <- 0
	fjacp_[3,42] <- 0
	fjacp_[4,42] <- 0
	fjacp_[5,42] <- 0
# column 43
	fjacp_[1,43] <- 0
	fjacp_[2,43] <- 0
	fjacp_[3,43] <- 0
	fjacp_[4,43] <- 0
	fjacp_[5,43] <- 0
# column 44
	fjacp_[1,44] <- 0
	fjacp_[2,44] <- 0
	fjacp_[3,44] <- 0
	fjacp_[4,44] <- 0
	fjacp_[5,44] <- 0
# column 45
	fjacp_[1,45] <- 0
	fjacp_[2,45] <- 0
	fjacp_[3,45] <- 0
	fjacp_[4,45] <- 0
	fjacp_[5,45] <- 0
# column 46
	fjacp_[1,46] <- 0
	fjacp_[2,46] <- 0
	fjacp_[3,46] <- 0
	fjacp_[4,46] <- 0
	fjacp_[5,46] <- 0
# column 47
	fjacp_[1,47] <- 0
	fjacp_[2,47] <- 0
	fjacp_[3,47] <- 0
	fjacp_[4,47] <- 0
	fjacp_[5,47] <- 0
# column 48
	fjacp_[1,48] <- 0
	fjacp_[2,48] <- 0
	fjacp_[3,48] <- 0
	fjacp_[4,48] <- 0
	fjacp_[5,48] <- 0
# column 49
	fjacp_[1,49] <- 0
	fjacp_[2,49] <- 0
	fjacp_[3,49] <- 0
	fjacp_[4,49] <- 0
	fjacp_[5,49] <- 0
# column 50
	fjacp_[1,50] <- 0
	fjacp_[2,50] <- 0
	fjacp_[3,50] <- 0
	fjacp_[4,50] <- 0
	fjacp_[5,50] <- 0
# column 51
	fjacp_[1,51] <- 0
	fjacp_[2,51] <- 0
	fjacp_[3,51] <- 0
	fjacp_[4,51] <- 0
	fjacp_[5,51] <- 0
# column 52
	fjacp_[1,52] <- 0
	fjacp_[2,52] <- 0
	fjacp_[3,52] <- 0
	fjacp_[4,52] <- 0
	fjacp_[5,52] <- 0
# column 53
	fjacp_[1,53] <- 0
	fjacp_[2,53] <- 0
	fjacp_[3,53] <- 0
	fjacp_[4,53] <- 0
	fjacp_[5,53] <- 0
# column 54
	fjacp_[1,54] <- 0
	fjacp_[2,54] <- 0
	fjacp_[3,54] <- 0
	fjacp_[4,54] <- 0
	fjacp_[5,54] <- 0
# column 55
	fjacp_[1,55] <- 0
	fjacp_[2,55] <- 0
	fjacp_[3,55] <- 0
	fjacp_[4,55] <- 0
	fjacp_[5,55] <- isOn
# column 56
	fjacp_[1,56] <- 0
	fjacp_[2,56] <- 0
	fjacp_[3,56] <- 0
	fjacp_[4,56] <- 0
	fjacp_[5,56] <- 0
# column 57
	fjacp_[1,57] <- 0
	fjacp_[2,57] <- 0
	fjacp_[3,57] <- 0
	fjacp_[4,57] <- 0
	fjacp_[5,57] <- 0
# column 58
	fjacp_[1,58] <- 0
	fjacp_[2,58] <- 0
	fjacp_[3,58] <- 0
	fjacp_[4,58] <- 0
	fjacp_[5,58] <- 0
# column 59
	fjacp_[1,59] <- 0
	fjacp_[2,59] <- 0
	fjacp_[3,59] <- -(PP2B_CaM_Ca4+PP2B_CaM_Ca3+PP2B_CaM_Ca2+PP2B_CaM_Ca1+PP2B_CaM)/(s+PP2B_0)^2
	fjacp_[4,59] <- -(100*PP2B_CaM_Ca4)/(s+PP2B_0)^2
	fjacp_[5,59] <- 0
# column 60
	fjacp_[1,60] <- 0
	fjacp_[2,60] <- 0
	fjacp_[3,60] <- 0
	fjacp_[4,60] <- 0
	fjacp_[5,60] <- CaBase
	rownames(fjacp_) <- c("CaPerCaM", "AutoCaMKII", "CaMPerPP2B", "ActivePP2B", "Camonitor")
	colnames(fjacp_) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "CaBase", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0", "isOn")
	return(fjacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
CaMKIIs_default<-function(t=0.0)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	parameters <- vector(mode='numeric',len=60)
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
	parameters[55] <- 0
	parameters[56] <- 0
	parameters[57] <- 0
	parameters[58] <- 30
	parameters[59] <- 3
	parameters[60] <- 0
	names(parameters) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "CaBase", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0", "isOn")
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
	CaBase <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
	isOn <- parameters[60]
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=21)
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
	names(state) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
	return(state)
}
model<-list(vf=CaMKIIs_vf, jac=CaMKIIs_jac, jacp=CaMKIIs_jacp, func=CaMKIIs_func, funcJac=CaMKIIs_funcJac, funcJacp=CaMKIIs_funcJacp, init=CaMKIIs_init, par=CaMKIIs_default, name="CaMKIIs")
