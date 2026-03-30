#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* Enums will be used for indexing purposes.   */
enum stateVariable { _AKAR4p, _C, numStateVar };
enum param { _kf_C_AKAR4, _kb_C_AKAR4, _kcat_AKARp, _AKAR4_C_ConservedConst, _AKAR4_ConservedConst, numParam };
enum func { _AKAR4pOUT, numFunc };

/* The error codes indicate how many values a function returns.                             */
/* Each function expects the output buffer to be allocated with at least that many values   */

/* ODE vector field: y' = f(t,y;p)   */
int AKAR4_vf(double t, const double y_[], double *f_, void *par){
	double *p_=par;
	if (!y_ || !f_) return 2;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
/* 	state variables   */
	double AKAR4p = y_[_AKAR4p];                                        /* [  0] */
	double C = y_[_C];                                                  /* [  1] */
/* 	expressions   */
	double AKAR4_C = AKAR4_C_ConservedConst - (+C);
	double AKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(f_,0,sizeof(double)*2); /* initialize with 0.0 */
	f_[0] = +reaction_2;
	f_[1] = -reaction_1+reaction_2;
	return GSL_SUCCESS;
}

/* ODE Jacobian: df(t,y;p)/dy   */
int AKAR4_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par){
	double *p_=par;
	if (!y_ || !jac_) return 4;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
/* 	state variables   */
	double AKAR4p = y_[_AKAR4p];                                        /* [  0] */
	double C = y_[_C];                                                  /* [  1] */
/* 	expressions   */
	double AKAR4_C = AKAR4_C_ConservedConst - (+C);
	double AKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(jac_,0,sizeof(double)*4); /* initialize with 0.0 */
	/*[ 0, 1]*/  jac_[1] = -kcat_AKARp;
	/*[ 1, 0]*/  jac_[2] = kf_C_AKAR4*C;
	/*[ 1, 1]*/  jac_[3] = -(kcat_AKARp+kb_C_AKAR4+kf_C_AKAR4*C+kf_C_AKAR4*(AKAR4_ConservedConst-(AKAR4p-C)));
	return GSL_SUCCESS;
}

/* ODE parameter Jacobian: df(t,y;p)/dp   */
int AKAR4_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par){
	double *p_=par;
	if (!y_ || !jacp_) return 10;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
/* 	state variables   */
	double AKAR4p = y_[_AKAR4p];                                        /* [  0] */
	double C = y_[_C];                                                  /* [  1] */
/* 	expressions   */
	double AKAR4_C = AKAR4_C_ConservedConst - (+C);
	double AKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(jacp_,0,sizeof(double)*10); /* initialize with 0.0 */
	/*[ 0, 2]*/  jacp_[2] = AKAR4_C_ConservedConst-C;
	/*[ 0, 3]*/  jacp_[3] = kcat_AKARp;
	/*[ 1, 0]*/  jacp_[5] = -C*(AKAR4_ConservedConst-(AKAR4p-C));
	/*[ 1, 1]*/  jacp_[6] = AKAR4_C_ConservedConst-C;
	/*[ 1, 2]*/  jacp_[7] = AKAR4_C_ConservedConst-C;
	/*[ 1, 3]*/  jacp_[8] = kb_C_AKAR4+kcat_AKARp;
	/*[ 1, 4]*/  jacp_[9] = -kf_C_AKAR4*C;
	return GSL_SUCCESS;
}

/* Output Function (Observables)   */
int AKAR4_func(double t, const double y_[], double *func_, void *par){
	double *p_=par;
	if (!y_ || !func_) return 1;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
/* 	state variables   */
	double AKAR4p = y_[_AKAR4p];                                        /* [  0] */
	double C = y_[_C];                                                  /* [  1] */
/* 	expressions   */
	double AKAR4_C = AKAR4_C_ConservedConst - (+C);
	double AKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	func_[_AKAR4pOUT] = 108 + 380*AKAR4p;
	return GSL_SUCCESS;
}

/* Output function Jacobian: dF(t,y;p)/dx   */
int AKAR4_funcJac(double t, const double y_[], double *funcJac_, void *par){
	double *p_=par;
	if (!y_ || !funcJac_) return 2;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
/* 	state variables   */
	double AKAR4p = y_[_AKAR4p];                                        /* [  0] */
	double C = y_[_C];                                                  /* [  1] */
/* 	expressions   */
	double AKAR4_C = AKAR4_C_ConservedConst - (+C);
	double AKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(funcJac_,0,sizeof(double)*2); /* initialize with 0.0 */
	/*[ 0, 0]*/  funcJac_[0] = 380;
	return GSL_SUCCESS;
}

/* Output function parameter Jacobian: dF(t,y;p)/dp   */
int AKAR4_funcJacp(double t, const double y_[], double *funcJacp_, void *par){
	double *p_=par;
	if (!y_ || !funcJacp_) return 5;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
/* 	state variables   */
	double AKAR4p = y_[_AKAR4p];                                        /* [  0] */
	double C = y_[_C];                                                  /* [  1] */
/* 	expressions   */
	double AKAR4_C = AKAR4_C_ConservedConst - (+C);
	double AKAR4 = AKAR4_ConservedConst - (+AKAR4p-1*C);
	double reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2 = kcat_AKARp*AKAR4_C;
	memset(funcJacp_,0,sizeof(double)*5); /* initialize with 0.0 */
	return GSL_SUCCESS;
}

int AKAR4_default(double t, double *p_){
	if (!p_) return numParam;
/* 	constants   */
	memset(p_,0,sizeof(double)*5); /* initialize with 0.0 */
	p_[_kf_C_AKAR4] = 0.018;
	p_[_kb_C_AKAR4] = 0.106;
	p_[_kcat_AKARp] = 10.2;
	p_[_AKAR4_ConservedConst] = 0.2;
	return GSL_SUCCESS;
}

int AKAR4_init(double t, double *y_, void *par){
	double *p_=par;
	if (!y_ || !y_) return 2;
/* 	constants   */
/* 	parameter values   */
	double kf_C_AKAR4 = p_[_kf_C_AKAR4];                                /* [  0] */
	double kb_C_AKAR4 = p_[_kb_C_AKAR4];                                /* [  1] */
	double kcat_AKARp = p_[_kcat_AKARp];                                /* [  2] */
	double AKAR4_C_ConservedConst = p_[_AKAR4_C_ConservedConst];        /* [  3] */
	double AKAR4_ConservedConst = p_[_AKAR4_ConservedConst];            /* [  4] */
	memset(y_,0,sizeof(double)*2); /* initialize with 0.0 */
	return GSL_SUCCESS;
}
