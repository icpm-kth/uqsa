# require("deSolve")

# ode vector field: y'=f(t,y;p)
AKAR4cl_vf <- function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	AKAR4p <- state[1]
	C <- state[2]
	AKAR4_C <- AKAR4_C_ConservedConst - (C)
	AKAR4 <- AKAR4_ConservedConst - (AKAR4p-C)
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	f_<-vector(mode='numeric',len=2)
	f_[1] <- +reaction_2
	f_[2] <- -reaction_1+reaction_2
	names(f_) <- c("AKAR4p", "C")
## for some weird reason deSolve wants this to be a list:
	return(list(f_))
}
# ode Jacobian df(t,y;p)/dy
AKAR4cl_jac<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	AKAR4p <- state[1]
	C <- state[2]
	AKAR4_C <- AKAR4_C_ConservedConst - (C)
	AKAR4 <- AKAR4_ConservedConst - (AKAR4p-C)
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	jac_ <- matrix(NA,2,2)
# column 1
	jac_[1,1] <- 0
	jac_[2,1] <- C*kf_C_AKAR4
# column 2
	jac_[1,2] <- -kcat_AKARp
	jac_[2,2] <- (-(C-AKAR4p+AKAR4_ConservedConst)*kf_C_AKAR4)-C*kf_C_AKAR4-kcat_AKARp-kb_C_AKAR4
	rownames(jac_) <- c("AKAR4p", "C")
	colnames(jac_) <- c("AKAR4p", "C")
	return(jac_)
}
# ode parameter Jacobian df(t,y;p)/dp
AKAR4cl_jacp<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	AKAR4p <- state[1]
	C <- state[2]
	AKAR4_C<-AKAR4_C_ConservedConst - (C)
	AKAR4<-AKAR4_ConservedConst - (AKAR4p-C)
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	jacp_ <- matrix(NA,2,5)
# column 1
	jacp_[1,1] <- 0
	jacp_[2,1] <- -C*(C-AKAR4p+AKAR4_ConservedConst)
# column 2
	jacp_[1,2] <- 0
	jacp_[2,2] <- AKAR4_C_ConservedConst-C
# column 3
	jacp_[1,3] <- AKAR4_C_ConservedConst-C
	jacp_[2,3] <- AKAR4_C_ConservedConst-C
# column 4
	jacp_[1,4] <- kcat_AKARp
	jacp_[2,4] <- kcat_AKARp+kb_C_AKAR4
# column 5
	jacp_[1,5] <- 0
	jacp_[2,5] <- -C*kf_C_AKAR4
	rownames(jacp_) <- c("AKAR4p", "C")
	colnames(jacp_) <- c("kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "AKAR4_C_ConservedConst", "AKAR4_ConservedConst")
	return(jacp_)
}
# ode Functions F(t,y;p)
AKAR4cl_func<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	AKAR4p <- state[1]
	C <- state[2]
	AKAR4_C <- AKAR4_C_ConservedConst - (C)
	AKAR4 <- AKAR4_ConservedConst - (AKAR4p-C)
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	func_ <- vector(mode='numeric',len=1)
	func_[1] <- 108 + 380*AKAR4p # AKAR4pOUT
	names(func_) <- c("AKAR4pOUT")
	return(func_)
}
# output function Jacobian dF(t,y;p)/dp
AKAR4cl_funcJac<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	AKAR4p <- state[1]
	C <- state[2]
	AKAR4_C<-AKAR4_C_ConservedConst - (C)
	AKAR4<-AKAR4_ConservedConst - (AKAR4p-C)
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	fjac_ <- matrix(NA,1,2)
# column 1
	fjac_[1,1] <- 380
# column 2
	fjac_[1,2] <- 0
	rownames(fjac_) <- c("AKAR4pOUT")
	colnames(fjac_) <- c("AKAR4p", "C")
	return(fjac_)
}
# output function parameter Jacobian dF(t,y;p)/dp
AKAR4cl_funcJacp<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	AKAR4p <- state[1]
	C <- state[2]
	AKAR4_C<-AKAR4_C_ConservedConst - (C)
	AKAR4<-AKAR4_ConservedConst - (AKAR4p-C)
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	fjacp_ <- matrix(NA,1,5)
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
	rownames(fjacp_) <- c("AKAR4pOUT")
	colnames(fjacp_) <- c("kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "AKAR4_C_ConservedConst", "AKAR4_ConservedConst")
	return(fjacp_)
}
# ode default parameters; can depend on constants, and time  of initialization
AKAR4cl_default<-function(t=0.0)
{
	parameters <- vector(mode='numeric',len=5)
	parameters[1] <- 0.018
	parameters[2] <- 0.106
	parameters[3] <- 10.2
	parameters[4] <- 0.000000
	parameters[5] <- 0.200000
	names(parameters) <- c("kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp", "AKAR4_C_ConservedConst", "AKAR4_ConservedConst")
	return(parameters);
}
# ode initial values
AKAR4cl_init<-function(t=0.0, parameters=NA)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4_C_ConservedConst <- parameters[4]
	AKAR4_ConservedConst <- parameters[5]
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=2)
	state[1] <- 0
	state[2] <- 0
	names(state) <- c("AKAR4p", "C")
	return(state)
}
model<-list(vf=AKAR4cl_vf, jac=AKAR4cl_jac, jacp=AKAR4cl_jacp, func=AKAR4cl_func, funcJac=AKAR4cl_funcJac, funcJacp=AKAR4cl_funcJacp, init=AKAR4cl_init, par=AKAR4cl_default, name="AKAR4cl")
