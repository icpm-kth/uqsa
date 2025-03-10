# require("deSolve")

# ode vector field: y'=f(t,y;p)
AKAR4_vf <- function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4 <- state[1]
	AKAR4_C <- state[2]
	AKAR4p <- state[3]
	C <- state[4]
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	f_<-vector(mode='numeric',len=4)
	f_[1] <- -reaction_1
	f_[2] <- +reaction_1-reaction_2
	f_[3] <- +reaction_2
	f_[4] <- -reaction_1+reaction_2
	names(f_) <- c("AKAR4", "AKAR4_C", "AKAR4p", "C")
## for some weird reason deSolve wants this to be a list:
	return(list(f_))
}
# ode Jacobian df(t,y;p)/dy
AKAR4_jac<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4 <- state[1]
	AKAR4_C <- state[2]
	AKAR4p <- state[3]
	C <- state[4]
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	jac_ <- matrix(NA,4,4)
# column 1 (df/dy_0)
	jac_[1,1] <- -C*kf_C_AKAR4
	jac_[2,1] <- C*kf_C_AKAR4
	jac_[3,1] <- 0
	jac_[4,1] <- -C*kf_C_AKAR4
# column 2 (df/dy_1)
	jac_[1,2] <- kb_C_AKAR4
	jac_[2,2] <- (-kcat_AKARp)-kb_C_AKAR4
	jac_[3,2] <- kcat_AKARp
	jac_[4,2] <- kcat_AKARp+kb_C_AKAR4
# column 3 (df/dy_2)
	jac_[1,3] <- 0
	jac_[2,3] <- 0
	jac_[3,3] <- 0
	jac_[4,3] <- 0
# column 4 (df/dy_3)
	jac_[1,4] <- -AKAR4*kf_C_AKAR4
	jac_[2,4] <- AKAR4*kf_C_AKAR4
	jac_[3,4] <- 0
	jac_[4,4] <- -AKAR4*kf_C_AKAR4
	rownames(jac_) <- c("AKAR4", "AKAR4_C", "AKAR4p", "C")
	colnames(jac_) <- c("AKAR4", "AKAR4_C", "AKAR4p", "C")
	return(jac_)
}
# ode parameter Jacobian df(t,y;p)/dp
AKAR4_jacp<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4 <- state[1]
	AKAR4_C <- state[2]
	AKAR4p <- state[3]
	C <- state[4]
	reaction_1<-kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2<-kcat_AKARp*AKAR4_C
	jacp_<-matrix(NA,4,3)
# column 1 (df/dp_1)
	jacp_[1,1] <- -AKAR4*C
	jacp_[2,1] <- AKAR4*C
	jacp_[3,1] <- 0
	jacp_[4,1] <- -AKAR4*C
# column 2 (df/dp_2)
	jacp_[1,2] <- AKAR4_C
	jacp_[2,2] <- -AKAR4_C
	jacp_[3,2] <- 0
	jacp_[4,2] <- AKAR4_C
# column 3 (df/dp_3)
	jacp_[1,3] <- 0
	jacp_[2,3] <- -AKAR4_C
	jacp_[3,3] <- AKAR4_C
	jacp_[4,3] <- AKAR4_C
	rownames(jacp_) <- c("AKAR4", "AKAR4_C", "AKAR4p", "C")
	colnames(jacp_) <- c("kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp")
	return(jacp_)
}
# ode Functions F(t,y;p)
AKAR4_func<-function(t, state, parameters)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	AKAR4 <- state[1]
	AKAR4_C <- state[2]
	AKAR4p <- state[3]
	C <- state[4]
	reaction_1 <- kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	reaction_2 <- kcat_AKARp*AKAR4_C
	func_ <- vector(mode='numeric',len=1)
	func_[1] <- 108 + 380*AKAR4p # AKAR4pOUT
	names(func_) <- c("AKAR4pOUT")
	return(func_)
}
# ode default parameters; can depend on constants, and time  of initialization
AKAR4_default<-function(t=0.0)
{
	parameters <- vector(mode='numeric',len=3)
	parameters[1] <- 0.018
	parameters[2] <- 0.106
	parameters[3] <- 10.2
	names(parameters) <- c("kf_C_AKAR4", "kb_C_AKAR4", "kcat_AKARp")
	return(parameters);
}
# ode initial values
AKAR4_init<-function(t=0.0, parameters=NA)
{
	kf_C_AKAR4 <- parameters[1]
	kb_C_AKAR4 <- parameters[2]
	kcat_AKARp <- parameters[3]
	# the initial value may depend on the parameters. 
	state<-vector(mode='numeric',len=4)
	state[1] <- 0.2
	state[2] <- 0
	state[3] <- 0
	state[4] <- 0
	names(state) <- c("AKAR4", "AKAR4_C", "AKAR4p", "C")
	return(state)
}
model<-list(vf=AKAR4_vf, jac=AKAR4_jac, jacp=AKAR4_jacp, func=AKAR4_func, init=AKAR4_init, par=AKAR4_default, name="AKAR4")
