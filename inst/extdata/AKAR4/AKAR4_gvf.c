#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

/* The error code indicates how to pre-allocate memory
 * for output values such as `f_`. The _vf function returns
 * the number of state variables, if any of the args are `NULL`.
 * evaluation errors can be indicated by negative return values.
 * GSL_SUCCESS (0) is returned when no error occurred.
 */

/* ode vector field: y'=f(t,y;p) */
int AKAR4_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 4;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4=y_[0];
	double AKAR4_C=y_[1];
	double AKAR4p=y_[2];
	double C=y_[3];
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	f_[0] = -reaction_1;
	f_[1] = +reaction_1-reaction_2;
	f_[2] = +reaction_2;
	f_[3] = -reaction_1+reaction_2;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int AKAR4_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 4*4;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4=y_[0];
	double AKAR4_C=y_[1];
	double AKAR4p=y_[2];
	double C=y_[3];
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dy_0) */
	jac_[0] = (-1*(kf_C_AKAR4*C)); /* [0, 0] */
	jac_[4] = (kf_C_AKAR4*C); /* [1, 0] */
	jac_[8] = 0; /* [2, 0] */
	jac_[12] = (-1*(kf_C_AKAR4*C)); /* [3, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = (-1*(0-kb_C_AKAR4)); /* [0, 1] */
	jac_[5] = ((0-kb_C_AKAR4)-kcat_AKARp); /* [1, 1] */
	jac_[9] = kcat_AKARp; /* [2, 1] */
	jac_[13] = ((-1*(0-kb_C_AKAR4))+kcat_AKARp); /* [3, 1] */
/* column 3 (df/dy_2) */
	jac_[2] = 0; /* [0, 2] */
	jac_[6] = 0; /* [1, 2] */
	jac_[10] = 0; /* [2, 2] */
	jac_[14] = 0; /* [3, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = (-1*(kf_C_AKAR4*AKAR4)); /* [0, 3] */
	jac_[7] = (kf_C_AKAR4*AKAR4); /* [1, 3] */
	jac_[11] = 0; /* [2, 3] */
	jac_[15] = (-1*(kf_C_AKAR4*AKAR4)); /* [3, 3] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAR4_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 4*3;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4=y_[0];
	double AKAR4_C=y_[1];
	double AKAR4p=y_[2];
	double C=y_[3];
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dp_0) */
	jacp_[0] = (-1*(C*AKAR4)); /* [0, 0] */
	jacp_[3] = (C*AKAR4); /* [1, 0] */
	jacp_[6] = 0; /* [2, 0] */
	jacp_[9] = (-1*(C*AKAR4)); /* [3, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = (-1*(0-AKAR4_C)); /* [0, 1] */
	jacp_[4] = (0-AKAR4_C); /* [1, 1] */
	jacp_[7] = 0; /* [2, 1] */
	jacp_[10] = (-1*(0-AKAR4_C)); /* [3, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 0; /* [0, 2] */
	jacp_[5] = (0-AKAR4_C); /* [1, 2] */
	jacp_[8] = AKAR4_C; /* [2, 2] */
	jacp_[11] = AKAR4_C; /* [3, 2] */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int AKAR4_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4=y_[0];
	double AKAR4_C=y_[1];
	double AKAR4p=y_[2];
	double C=y_[3];
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	func_[0] = AKAR4p; /* AKAR4pOUT */
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAR4_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 3;
	p_[0] = 0.018;
	p_[1] = 0.106;
	p_[2] = 10.2;
	return GSL_SUCCESS;
}
/* ode initial values */
int AKAR4_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return        4;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	/* the initial value of y may depend on the parameters. */
	y_[0] = 0.2;
	y_[1] = 0;
	y_[2] = 0;
	y_[3] = 0;
	return GSL_SUCCESS;
}
