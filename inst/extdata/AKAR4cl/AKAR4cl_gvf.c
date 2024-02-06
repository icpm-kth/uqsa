#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* The error code indicates how to pre-allocate memory
 * for output values such as `f_`. The _vf function returns
 * the number of state variables, if any of the args are `NULL`.
 * evaluation errors can be indicated by negative return values.
 * GSL_SUCCESS (0) is returned when no error occurred.
 */

/* ode vector field: y'=f(t,y;p) */
int AKAR4cl_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 2;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	double AKAR4p=y_[0];
	double C=y_[1];
	double AKAR4_C=AKAR4_C_ConservedConst - (C);
	double AKAR4=AKAR4_ConservedConst - (AKAR4p-C);
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	f_[0] = +reaction_2;
	f_[1] = -reaction_1+reaction_2;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int AKAR4cl_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 2*2;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	double AKAR4p=y_[0];
	double C=y_[1];
	double AKAR4_C=AKAR4_C_ConservedConst - (C);
	double AKAR4=AKAR4_ConservedConst - (AKAR4p-C);
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dy_0) */
	jac_[0] = 0; /* [0, 0] */
	jac_[2] = C*kf_C_AKAR4; /* [1, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = -kcat_AKARp; /* [0, 1] */
	jac_[3] = (-(C-AKAR4p+AKAR4_ConservedConst)*kf_C_AKAR4)-C*kf_C_AKAR4-kcat_AKARp-kb_C_AKAR4; /* [1, 1] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAR4cl_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 2*5;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	double AKAR4p=y_[0];
	double C=y_[1];
	double AKAR4_C=AKAR4_C_ConservedConst - (C);
	double AKAR4=AKAR4_ConservedConst - (AKAR4p-C);
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dp_0) */
	jacp_[0] = 0; /* [0, 0] */
	jacp_[5] = -C*(C-AKAR4p+AKAR4_ConservedConst); /* [1, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = 0; /* [0, 1] */
	jacp_[6] = AKAR4_C_ConservedConst-C; /* [1, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = AKAR4_C_ConservedConst-C; /* [0, 2] */
	jacp_[7] = AKAR4_C_ConservedConst-C; /* [1, 2] */
/* column 4 (df/dp_3) */
	jacp_[3] = kcat_AKARp; /* [0, 3] */
	jacp_[8] = kcat_AKARp+kb_C_AKAR4; /* [1, 3] */
/* column 5 (df/dp_4) */
	jacp_[4] = 0; /* [0, 4] */
	jacp_[9] = -C*kf_C_AKAR4; /* [1, 4] */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int AKAR4cl_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	double AKAR4p=y_[0];
	double C=y_[1];
	double AKAR4_C=AKAR4_C_ConservedConst - (C);
	double AKAR4=AKAR4_ConservedConst - (AKAR4p-C);
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	func_[0] = 108 + 380*AKAR4p; /* AKAR4pOUT */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int AKAR4cl_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 2;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	double AKAR4p=y_[0];
	double C=y_[1];
	double AKAR4_C=AKAR4_C_ConservedConst - (C);
	double AKAR4=AKAR4_ConservedConst - (AKAR4p-C);
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (dF/dy_0) */
	funcJac_[0] = 380; /* [0, 0] */
/* column 2 (dF/dy_1) */
	funcJac_[1] = 0; /* [0, 1] */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int AKAR4cl_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 5;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	double AKAR4p=y_[0];
	double C=y_[1];
	double AKAR4_C=AKAR4_C_ConservedConst - (C);
	double AKAR4=AKAR4_ConservedConst - (AKAR4p-C);
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (dF/dp_0) */
	funcJacp_[0] = 0; /* [0, 0] */
/* column 2 (dF/dp_1) */
	funcJacp_[1] = 0; /* [0, 1] */
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAR4cl_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 5;
	p_[0] = 0.018;
	p_[1] = 0.106;
	p_[2] = 10.2;
	p_[3] = 0.000000;
	p_[4] = 0.200000;
	return GSL_SUCCESS;
}
/* ode initial values */
int AKAR4cl_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 2;
	double kf_C_AKAR4=p_[0];
	double kb_C_AKAR4=p_[1];
	double kcat_AKARp=p_[2];
	double AKAR4_C_ConservedConst=p_[3];
	double AKAR4_ConservedConst=p_[4];
	/* the initial value of y may depend on the parameters. */
	y_[0] = 0;
	y_[1] = 0;
	return GSL_SUCCESS;
}
