require("deSolve")

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
	logistic <- 1.0/(1+exp(-10*t))
	Ca <- logistic*Ca_set
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
	f_[1] <- +ReactionFlux1-ReactionFlux2-ReactionFlux6-ReactionFlux15-ReactionFlux26
	f_[2] <- +ReactionFlux2-ReactionFlux3-ReactionFlux7-ReactionFlux16-ReactionFlux25
	f_[3] <- +ReactionFlux3-ReactionFlux4-ReactionFlux8-ReactionFlux17-ReactionFlux24
	f_[4] <- +ReactionFlux4-ReactionFlux9-ReactionFlux18-ReactionFlux23
	f_[5] <- +ReactionFlux5-ReactionFlux10
	f_[6] <- +ReactionFlux6+ReactionFlux10-ReactionFlux11
	f_[7] <- +ReactionFlux7+ReactionFlux11-ReactionFlux12
	f_[8] <- +ReactionFlux8+ReactionFlux12-ReactionFlux13
	f_[9] <- +ReactionFlux9+ReactionFlux13
	f_[10] <- +ReactionFlux14-ReactionFlux19
	f_[11] <- +ReactionFlux15+ReactionFlux19-ReactionFlux20
	f_[12] <- +ReactionFlux16+ReactionFlux20-ReactionFlux21
	f_[13] <- +ReactionFlux17+ReactionFlux21-ReactionFlux22
	f_[14] <- +ReactionFlux18+ReactionFlux22-ReactionFlux32
	f_[15] <- +ReactionFlux23+ReactionFlux31+ReactionFlux32
	f_[16] <- -ReactionFlux23-ReactionFlux24-ReactionFlux25-ReactionFlux26-ReactionFlux27-ReactionFlux33
	f_[17] <- +ReactionFlux24+ReactionFlux30-ReactionFlux31
	f_[18] <- +ReactionFlux25+ReactionFlux29-ReactionFlux30
	f_[19] <- +ReactionFlux26+ReactionFlux28-ReactionFlux29
	f_[20] <- +ReactionFlux27-ReactionFlux28
	f_[21] <- +ReactionFlux33-ReactionFlux34
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
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
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
	logistic <- 1.0/(1+exp(-10*t))
	Ca <- logistic*Ca_set
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
# column 1 (df/dy_0)
	jac_[1,1] <- ((0-(kf__CaM_Ca1__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,1] <- 0
# column 2 (df/dy_1)
	jac_[1,2] <- ((0-(kf__CaM_Ca2__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,2] <- 0
# column 3 (df/dy_2)
	jac_[1,3] <- ((0-(kf__CaM_Ca3__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,3] <- 0
# column 4 (df/dy_3)
	jac_[1,4] <- ((-1*(kf__CaM_Ca4__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,4] <- 0
# column 5 (df/dy_4)
	jac_[1,5] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,5] <- 0
# column 6 (df/dy_5)
	jac_[1,6] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,6] <- 0
# column 7 (df/dy_6)
	jac_[1,7] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,7] <- 0
# column 8 (df/dy_7)
	jac_[1,8] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,8] <- 0
# column 9 (df/dy_8)
	jac_[1,9] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,9] <- 0
# column 10 (df/dy_9)
	jac_[1,10] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,10] <- 0
# column 11 (df/dy_10)
	jac_[1,11] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,11] <- 0
# column 12 (df/dy_11)
	jac_[1,12] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,12] <- 0
# column 13 (df/dy_12)
	jac_[1,13] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,13] <- 0
# column 14 (df/dy_13)
	jac_[1,14] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,14] <- 0
# column 15 (df/dy_14)
	jac_[1,15] <- ((-1*(0-(KD__CaM_Ca4__pCaMKIIa*kf__CaM_Ca4__pCaMKIIa)))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,15] <- 0
# column 16 (df/dy_15)
	jac_[1,16] <- ((((((-1*(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4))-(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3))-(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2))-(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1))-(kf__CaM__pCaMKIIa*((((CaM_0+PP1_0)-CaMKII_0)-(((((((((((CaM_Ca1+CaM_Ca2)+CaM_Ca3)+CaM_Ca4)+PP2B_CaM)+PP2B_CaM_Ca1)+PP2B_CaM_Ca2)+PP2B_CaM_Ca3)+PP2B_CaM_Ca4)+(PP1_0-PP1__pCaMKIIa))-((CaMKII_0-PP1_0)-(((((((((((CaMKII_CaM+CaMKII_CaM_Ca1)+CaMKII_CaM_Ca2)+CaMKII_CaM_Ca3)+CaMKII_CaM_Ca4)+pCaMKII_CaM_Ca4)+pCaMKIIa)+pCaMKII_CaM_Ca3)+pCaMKII_CaM_Ca2)+pCaMKII_CaM_Ca1)+pCaMKII_CaM)-(PP1_0-PP1__pCaMKIIa))))-pCaMKIIa))-(0*pCaMKII_CaM))))-(kf__PP1__pCaMKIIa*(PP1_0-PP1__pCaMKIIa)))
	jac_[2,16] <- (kf__PP1__pCaMKIIa*(PP1_0-PP1__pCaMKIIa))
# column 17 (df/dy_16)
	jac_[1,17] <- ((0-(0-(kf__CaM_Ca3__pCaMKIIa*((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca))))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,17] <- 0
# column 18 (df/dy_17)
	jac_[1,18] <- ((0-(0-(KD__CaM_Ca2__Ca*(kf__CaM_Ca2__pCaMKIIa*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)))))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,18] <- 0
# column 19 (df/dy_18)
	jac_[1,19] <- ((0-(0-(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(kf__CaM_Ca1__pCaMKIIa*((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca))))))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM))))
	jac_[2,19] <- 0
# column 20 (df/dy_19)
	jac_[1,20] <- (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(KD__CaM__Ca*(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)/KD__pCaMKII_CaM__Ca)))))))
	jac_[2,20] <- 0
# column 21 (df/dy_20)
	jac_[1,21] <- (0-(((kf__PP1__pCaMKIIa*pCaMKIIa)*-1)-kr__PP1__pCaMKIIa))
	jac_[2,21] <- ((((kf__PP1__pCaMKIIa*pCaMKIIa)*-1)-kr__PP1__pCaMKIIa)-kcat__PP1__pCaMKIIa)
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
	Ca_set <- parameters[55]
	PP1_0 <- parameters[56]
	CaMKII_0 <- parameters[57]
	CaM_0 <- parameters[58]
	PP2B_0 <- parameters[59]
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
	logistic<-1.0/(1+exp(-10*t))
	Ca<-logistic*Ca_set
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
	jacp_<-matrix(NA,21,59)
# column 1 (df/dp_1)
	jacp_[1,1] <- 0
	jacp_[2,1] <- 0
# column 2 (df/dp_2)
	jacp_[1,2] <- 0
	jacp_[2,2] <- 0
# column 3 (df/dp_3)
	jacp_[1,3] <- 0
	jacp_[2,3] <- 0
# column 4 (df/dp_4)
	jacp_[1,4] <- 0
	jacp_[2,4] <- 0
# column 5 (df/dp_5)
	jacp_[1,5] <- 0
	jacp_[2,5] <- 0
# column 6 (df/dp_6)
	jacp_[1,6] <- 0
	jacp_[2,6] <- 0
# column 7 (df/dp_7)
	jacp_[1,7] <- 0
	jacp_[2,7] <- 0
# column 8 (df/dp_8)
	jacp_[1,8] <- 0
	jacp_[2,8] <- 0
# column 9 (df/dp_9)
	jacp_[1,9] <- 0
	jacp_[2,9] <- 0
# column 10 (df/dp_10)
	jacp_[1,10] <- 0
	jacp_[2,10] <- 0
# column 11 (df/dp_11)
	jacp_[1,11] <- 0
	jacp_[2,11] <- 0
# column 12 (df/dp_12)
	jacp_[1,12] <- 0
	jacp_[2,12] <- 0
# column 13 (df/dp_13)
	jacp_[1,13] <- 0
	jacp_[2,13] <- 0
# column 14 (df/dp_14)
	jacp_[1,14] <- ((((0-(0-(((KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca)*kf__CaM_Ca3__pCaMKIIa)*pCaMKII_CaM_Ca3)))-(0-((((KD__CaM_Ca2__Ca*(KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,14] <- 0
# column 15 (df/dp_15)
	jacp_[1,15] <- (((0-(0-(((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,15] <- 0
# column 16 (df/dp_16)
	jacp_[1,16] <- ((0-(0-((((KD__CaM_Ca2__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca2__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,16] <- 0
# column 17 (df/dp_17)
	jacp_[1,17] <- (0-(kf__CaM__pCaMKIIa*(0-(((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,17] <- 0
# column 18 (df/dp_18)
	jacp_[1,18] <- 0
	jacp_[2,18] <- 0
# column 19 (df/dp_19)
	jacp_[1,19] <- 0
	jacp_[2,19] <- 0
# column 20 (df/dp_20)
	jacp_[1,20] <- 0
	jacp_[2,20] <- 0
# column 21 (df/dp_21)
	jacp_[1,21] <- 0
	jacp_[2,21] <- 0
# column 22 (df/dp_22)
	jacp_[1,22] <- 0
	jacp_[2,22] <- 0
# column 23 (df/dp_23)
	jacp_[1,23] <- 0
	jacp_[2,23] <- 0
# column 24 (df/dp_24)
	jacp_[1,24] <- 0
	jacp_[2,24] <- 0
# column 25 (df/dp_25)
	jacp_[1,25] <- 0
	jacp_[2,25] <- 0
# column 26 (df/dp_26)
	jacp_[1,26] <- 0
	jacp_[2,26] <- 0
# column 27 (df/dp_27)
	jacp_[1,27] <- 0
	jacp_[2,27] <- 0
# column 28 (df/dp_28)
	jacp_[1,28] <- 0
	jacp_[2,28] <- 0
# column 29 (df/dp_29)
	jacp_[1,29] <- 0
	jacp_[2,29] <- 0
# column 30 (df/dp_30)
	jacp_[1,30] <- 0
	jacp_[2,30] <- 0
# column 31 (df/dp_31)
	jacp_[1,31] <- 0
	jacp_[2,31] <- 0
# column 32 (df/dp_32)
	jacp_[1,32] <- 0
	jacp_[2,32] <- 0
# column 33 (df/dp_33)
	jacp_[1,33] <- 0
	jacp_[2,33] <- 0
# column 34 (df/dp_34)
	jacp_[1,34] <- 0
	jacp_[2,34] <- 0
# column 35 (df/dp_35)
	jacp_[1,35] <- 0
	jacp_[2,35] <- 0
# column 36 (df/dp_36)
	jacp_[1,36] <- 0
	jacp_[2,36] <- 0
# column 37 (df/dp_37)
	jacp_[1,37] <- 0
	jacp_[2,37] <- 0
# column 38 (df/dp_38)
	jacp_[1,38] <- (0-(((((CaM_0+PP1_0)-CaMKII_0)-(((((((((((CaM_Ca1+CaM_Ca2)+CaM_Ca3)+CaM_Ca4)+PP2B_CaM)+PP2B_CaM_Ca1)+PP2B_CaM_Ca2)+PP2B_CaM_Ca3)+PP2B_CaM_Ca4)+(PP1_0-PP1__pCaMKIIa))-((CaMKII_0-PP1_0)-(((((((((((CaMKII_CaM+CaMKII_CaM_Ca1)+CaMKII_CaM_Ca2)+CaMKII_CaM_Ca3)+CaMKII_CaM_Ca4)+pCaMKII_CaM_Ca4)+pCaMKIIa)+pCaMKII_CaM_Ca3)+pCaMKII_CaM_Ca2)+pCaMKII_CaM_Ca1)+pCaMKII_CaM)-(PP1_0-PP1__pCaMKIIa))))-pCaMKIIa))*pCaMKIIa)-((KD__CaM__Ca*(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)/KD__pCaMKII_CaM__Ca))))*pCaMKII_CaM)))
	jacp_[2,38] <- 0
# column 39 (df/dp_39)
	jacp_[1,39] <- (0-((CaM_Ca1*pCaMKIIa)-((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)))*pCaMKII_CaM_Ca1)))
	jacp_[2,39] <- 0
# column 40 (df/dp_40)
	jacp_[1,40] <- (0-((CaM_Ca2*pCaMKIIa)-((KD__CaM_Ca2__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))*pCaMKII_CaM_Ca2)))
	jacp_[2,40] <- 0
# column 41 (df/dp_41)
	jacp_[1,41] <- (0-((CaM_Ca3*pCaMKIIa)-(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)*pCaMKII_CaM_Ca3)))
	jacp_[2,41] <- 0
# column 42 (df/dp_42)
	jacp_[1,42] <- 0
	jacp_[2,42] <- 0
# column 43 (df/dp_43)
	jacp_[1,43] <- 0
	jacp_[2,43] <- 0
# column 44 (df/dp_44)
	jacp_[1,44] <- (-1*((CaM_Ca4*pCaMKIIa)-(KD__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)))
	jacp_[2,44] <- 0
# column 45 (df/dp_45)
	jacp_[1,45] <- 0
	jacp_[2,45] <- 0
# column 46 (df/dp_46)
	jacp_[1,46] <- ((((0-(0-((((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca))*kf__CaM_Ca3__pCaMKIIa)*pCaMKII_CaM_Ca3)))-(0-((((KD__CaM_Ca2__Ca*((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca)))/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca)))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca)))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,46] <- 0
# column 47 (df/dp_47)
	jacp_[1,47] <- (((0-(0-((((KD__CaM_Ca2__Ca*(0-((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)))/(KD__pCaMKII_CaM_Ca2__Ca*KD__pCaMKII_CaM_Ca2__Ca))*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(0-((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)))/(KD__pCaMKII_CaM_Ca2__Ca*KD__pCaMKII_CaM_Ca2__Ca)))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(0-((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)))/(KD__pCaMKII_CaM_Ca2__Ca*KD__pCaMKII_CaM_Ca2__Ca)))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,47] <- 0
# column 48 (df/dp_48)
	jacp_[1,48] <- ((0-(0-((((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(0-(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))))/(KD__pCaMKII_CaM_Ca1__Ca*KD__pCaMKII_CaM_Ca1__Ca))*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(0-(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))))/(KD__pCaMKII_CaM_Ca1__Ca*KD__pCaMKII_CaM_Ca1__Ca)))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,48] <- 0
# column 49 (df/dp_49)
	jacp_[1,49] <- (0-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(0-((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)))))/(KD__pCaMKII_CaM__Ca*KD__pCaMKII_CaM__Ca))*pCaMKII_CaM))))
	jacp_[2,49] <- 0
# column 50 (df/dp_50)
	jacp_[1,50] <- (((((-1*(0-(kf__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)))-(0-(((KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca)*kf__CaM_Ca3__pCaMKIIa)*pCaMKII_CaM_Ca3)))-(0-((((KD__CaM_Ca2__Ca*(KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM))))
	jacp_[2,50] <- 0
# column 51 (df/dp_51)
	jacp_[1,51] <- 0
	jacp_[2,51] <- 0
# column 52 (df/dp_52)
	jacp_[1,52] <- (0-(pCaMKIIa*(PP1_0-PP1__pCaMKIIa)))
	jacp_[2,52] <- (pCaMKIIa*(PP1_0-PP1__pCaMKIIa))
# column 53 (df/dp_53)
	jacp_[1,53] <- (0-(0-PP1__pCaMKIIa))
	jacp_[2,53] <- (0-PP1__pCaMKIIa)
# column 54 (df/dp_54)
	jacp_[1,54] <- 0
	jacp_[2,54] <- (0-PP1__pCaMKIIa)
# column 55 (df/dp_55)
	jacp_[1,55] <- 0
	jacp_[2,55] <- 0
# column 56 (df/dp_56)
	jacp_[1,56] <- (0-(kf__PP1__pCaMKIIa*pCaMKIIa))
	jacp_[2,56] <- (kf__PP1__pCaMKIIa*pCaMKIIa)
# column 57 (df/dp_57)
	jacp_[1,57] <- 0
	jacp_[2,57] <- 0
# column 58 (df/dp_58)
	jacp_[1,58] <- (0-(kf__CaM__pCaMKIIa*(pCaMKIIa-(0*pCaMKII_CaM))))
	jacp_[2,58] <- 0
# column 59 (df/dp_59)
	jacp_[1,59] <- 0
	jacp_[2,59] <- 0
	rownames(jacp_) <- c("CaM_Ca1", "CaM_Ca2", "CaM_Ca3", "CaM_Ca4", "PP2B_CaM", "PP2B_CaM_Ca1", "PP2B_CaM_Ca2", "PP2B_CaM_Ca3", "PP2B_CaM_Ca4", "CaMKII_CaM", "CaMKII_CaM_Ca1", "CaMKII_CaM_Ca2", "CaMKII_CaM_Ca3", "CaMKII_CaM_Ca4", "pCaMKII_CaM_Ca4", "pCaMKIIa", "pCaMKII_CaM_Ca3", "pCaMKII_CaM_Ca2", "pCaMKII_CaM_Ca1", "pCaMKII_CaM", "PP1__pCaMKIIa")
	colnames(jacp_) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "Ca_set", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0")
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
	logistic <- 1.0/(1+exp(-10*t))
	Ca <- logistic*Ca_set
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
	func_ <- vector(mode='numeric',len=4)
	func_[1] <- BoundCa / (s + Total_CaM) # MolCaPerMolCaM
	func_[2] <- 100*(Total_pCaMKII / (s + totalCaMKII)) # AutoCamkiiPercentage
	func_[3] <- Total_PP2B_CaM_CaX / (s + PP2B + Total_PP2B_CaM_CaX) # MolCaMPerMolPP2B
	func_[4] <- 100 * (PP2B_CaM_Ca4 / (s + PP2B + Total_PP2B_CaM_CaX)) # ActivePP2BPercentage
	names(func_) <- c("MolCaPerMolCaM", "AutoCamkiiPercentage", "MolCaMPerMolPP2B", "ActivePP2BPercentage")
	return(func_)
}
# ode default parameters; can depend on constants, and time  of initialization
CaMKIIs_default<-function(t)
{
	a <- 3.90264
	b <- 2.86972
	s <- 1e-05
	parameters <- vector(mode='numeric',len=59)
	parameters[1] <- 0.15191
	parameters[2] <- 3.4245e-05
	parameters[3] <- 0.073893
	parameters[4] <- 0.0061814
	parameters[5] <- 0.020293
	parameters[6] <- 0.0045248
	parameters[7] <- 0.051176
	parameters[8] <- 0.27421
	parameters[9] <- 0.083357
	parameters[10] <- 0.0011578
	parameters[11] <- 0.0047884
	parameters[12] <- 0.035079
	parameters[13] <- 0.045566
	parameters[14] <- 7271.3
	parameters[15] <- 37062
	parameters[16] <- 1827.9
	parameters[17] <- 2662.3
	parameters[18] <- 0.03997
	parameters[19] <- 91.543
	parameters[20] <- 916.15
	parameters[21] <- 285.03
	parameters[22] <- 82.837
	parameters[23] <- 0.23745
	parameters[24] <- 0.025858
	parameters[25] <- 0.13086
	parameters[26] <- 0.075539
	parameters[27] <- 0.00079772
	parameters[28] <- 0.055882
	parameters[29] <- 0.046028
	parameters[30] <- 0.20855
	parameters[31] <- 0.022662
	parameters[32] <- 8.2849
	parameters[33] <- 483.48
	parameters[34] <- 1143.6
	parameters[35] <- 645.07
	parameters[36] <- 3081.6
	parameters[37] <- 0.00082984
	parameters[38] <- 0.00032583
	parameters[39] <- 0.058928
	parameters[40] <- 0.02319
	parameters[41] <- 0.030252
	parameters[42] <- 0.038498
	parameters[43] <- 0.0004565
	parameters[44] <- 0.072237
	parameters[45] <- 0.0021267
	parameters[46] <- 539.41
	parameters[47] <- 1784
	parameters[48] <- 57728
	parameters[49] <- 1342.9
	parameters[50] <- 3.746
	parameters[51] <- 0.0037559
	parameters[52] <- 0.0016604
	parameters[53] <- 0.20517
	parameters[54] <- 0.30225
	parameters[55] <- 2187.8
	parameters[56] <- 0
	parameters[57] <- 0
	parameters[58] <- 30
	parameters[59] <- 3
	names(parameters) <- c("kf__CaM__Ca", "kf__CaM_Ca1__Ca", "kf__CaM_Ca2__Ca", "kf__CaM_Ca3__Ca", "kf__CaM__PP2B", "kf__CaM_Ca1__PP2B", "kf__CaM_Ca2__PP2B", "kf__CaM_Ca3__PP2B", "kf__CaM_Ca4__PP2B", "kf__PP2B_CaM__Ca", "kf__PP2B_CaM_Ca1__Ca", "kf__PP2B_CaM_Ca2__Ca", "kf__PP2B_CaM_Ca3__Ca", "KD__CaM_Ca3__Ca", "KD__CaM_Ca2__Ca", "KD__CaM_Ca1__Ca", "KD__CaM__Ca", "KD__CaM_Ca4__PP2B", "KD__PP2B_CaM_Ca3__Ca", "KD__PP2B_CaM_Ca2__Ca", "KD__PP2B_CaM_Ca1__Ca", "KD__PP2B_CaM__Ca", "kf__CaM__CaMKII", "kf__CaMKII_CaM_Ca3__Ca", "kf__CaMKII_CaM_Ca2__Ca", "kf__CaMKII_CaM_Ca1__Ca", "kf__CaMKII_CaM__Ca", "kf__CaM_Ca1__CaMKII", "kf__CaM_Ca2__CaMKII", "kf__CaM_Ca3__CaMKII", "kf__CaM_Ca4__CaMKII", "KD__CaM_Ca4__CaMKII", "KD__CaMKII_CaM_Ca3__Ca", "KD__CaMKII_CaM_Ca2__Ca", "KD__CaMKII_CaM_Ca1__Ca", "KD__CaMKII_CaM__Ca", "kf__pCaMKII_CaM_Ca3__Ca", "kf__CaM__pCaMKIIa", "kf__CaM_Ca1__pCaMKIIa", "kf__CaM_Ca2__pCaMKIIa", "kf__CaM_Ca3__pCaMKIIa", "kf__pCaMKII_CaM_Ca2__Ca", "kf__pCaMKII_CaM_Ca1__Ca", "kf__CaM_Ca4__pCaMKIIa", "kf__pCaMKII_CaM__Ca", "KD__pCaMKII_CaM_Ca3__Ca", "KD__pCaMKII_CaM_Ca2__Ca", "KD__pCaMKII_CaM_Ca1__Ca", "KD__pCaMKII_CaM__Ca", "KD__CaM_Ca4__pCaMKIIa", "kautMax", "kf__PP1__pCaMKIIa", "kr__PP1__pCaMKIIa", "kcat__PP1__pCaMKIIa", "Ca_set", "PP1_0", "CaMKII_0", "CaM_0", "PP2B_0")
	return(parameters);
}
# ode initial values
CaMKIIs_init<-function(t, parameters)
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
model<-list(vf=CaMKIIs_vf, jac=CaMKIIs_jac, jacp=CaMKIIs_jacp, func=CaMKIIs_func, init=CaMKIIs_init, par=CaMKIIs_default, name="CaMKIIs")
