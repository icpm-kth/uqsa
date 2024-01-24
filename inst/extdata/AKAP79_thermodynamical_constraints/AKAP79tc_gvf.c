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
int AKAP79tc_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 16;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double kf_RiixC__Rii_C=(kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP)/(kb_cAMPxRii__Rii_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP);
	double kf_RiiPXcAMP__RiiP_cAMP=(kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP *kb_RiiP_cAMPxC__RiiP_C_cAMP)/(kb_RiiPxC__RiiP_C*kb_RiiP_CxcAMP__RiiP_C_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP);
	double reaction_51=kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14=kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12=kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43=kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23=kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78=kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56=kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76=kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62=kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58=kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	f_[0] = -reaction_78-reaction_58+reaction_48;
	f_[1] = -reaction_12-reaction_43-reaction_78-reaction_56;
	f_[2] = -reaction_14-reaction_43-reaction_44;
	f_[3] = -reaction_51-reaction_56+reaction_58;
	f_[4] = +reaction_43-reaction_23-reaction_33;
	f_[5] = +reaction_51+reaction_14-reaction_12;
	f_[6] = +reaction_12+reaction_23+reaction_62;
	f_[7] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2;
	f_[8] = +reaction_78-reaction_76+reaction_37;
	f_[9] = +reaction_56+reaction_76-reaction_62;
	f_[10] = -reaction_44-reaction_33+reaction_48+reaction_37;
	f_[11] = +reaction_44-reaction_48;
	f_[12] = +reaction_33-reaction_37;
	f_[13] = -reaction_1;
	f_[14] = +reaction_1-reaction_2;
	f_[15] = +reaction_2;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int AKAP79tc_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 16*16;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double kf_RiixC__Rii_C=(kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP)/(kb_cAMPxRii__Rii_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP);
	double kf_RiiPXcAMP__RiiP_cAMP=(kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP *kb_RiiP_cAMPxC__RiiP_C_cAMP)/(kb_RiiPxC__RiiP_C*kb_RiiP_CxcAMP__RiiP_C_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP);
	double reaction_51=kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14=kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12=kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43=kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23=kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78=kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56=kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76=kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62=kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58=kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dy_0) */
/* column 2 (df/dy_1) */
/* column 3 (df/dy_2) */
/* column 4 (df/dy_3) */
/* column 5 (df/dy_4) */
/* column 6 (df/dy_5) */
/* column 7 (df/dy_6) */
/* column 8 (df/dy_7) */
/* column 9 (df/dy_8) */
/* column 10 (df/dy_9) */
/* column 11 (df/dy_10) */
/* column 12 (df/dy_11) */
/* column 13 (df/dy_12) */
/* column 14 (df/dy_13) */
/* column 15 (df/dy_14) */
/* column 16 (df/dy_15) */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAP79tc_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 16*26;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double kf_RiixC__Rii_C=(kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP)/(kb_cAMPxRii__Rii_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP);
	double kf_RiiPXcAMP__RiiP_cAMP=(kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP *kb_RiiP_cAMPxC__RiiP_C_cAMP)/(kb_RiiPxC__RiiP_C*kb_RiiP_CxcAMP__RiiP_C_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP);
	double reaction_51=kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14=kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12=kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43=kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23=kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78=kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56=kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76=kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62=kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58=kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dp_0) */
/* column 2 (df/dp_1) */
/* column 3 (df/dp_2) */
/* column 4 (df/dp_3) */
/* column 5 (df/dp_4) */
/* column 6 (df/dp_5) */
/* column 7 (df/dp_6) */
/* column 8 (df/dp_7) */
/* column 9 (df/dp_8) */
/* column 10 (df/dp_9) */
/* column 11 (df/dp_10) */
/* column 12 (df/dp_11) */
/* column 13 (df/dp_12) */
/* column 14 (df/dp_13) */
/* column 15 (df/dp_14) */
/* column 16 (df/dp_15) */
/* column 17 (df/dp_16) */
/* column 18 (df/dp_17) */
/* column 19 (df/dp_18) */
/* column 20 (df/dp_19) */
/* column 21 (df/dp_20) */
/* column 22 (df/dp_21) */
/* column 23 (df/dp_22) */
/* column 24 (df/dp_23) */
/* column 25 (df/dp_24) */
/* column 26 (df/dp_25) */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int AKAP79tc_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double kf_RiixC__Rii_C=(kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP)/(kb_cAMPxRii__Rii_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP);
	double kf_RiiPXcAMP__RiiP_cAMP=(kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP *kb_RiiP_cAMPxC__RiiP_C_cAMP)/(kb_RiiPxC__RiiP_C*kb_RiiP_CxcAMP__RiiP_C_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP);
	double reaction_51=kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14=kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12=kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43=kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23=kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78=kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56=kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76=kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62=kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58=kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	func_[0] = (AKAR4p*5)*71.67+100; /* AKAR4pOUT */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int AKAP79tc_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 16;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double kf_RiixC__Rii_C=(kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP)/(kb_cAMPxRii__Rii_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP);
	double kf_RiiPXcAMP__RiiP_cAMP=(kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP *kb_RiiP_cAMPxC__RiiP_C_cAMP)/(kb_RiiPxC__RiiP_C*kb_RiiP_CxcAMP__RiiP_C_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP);
	double reaction_51=kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14=kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12=kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43=kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23=kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78=kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56=kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76=kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62=kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58=kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (dF/dy_0) */
/* column 2 (dF/dy_1) */
/* column 3 (dF/dy_2) */
/* column 4 (dF/dy_3) */
/* column 5 (dF/dy_4) */
/* column 6 (dF/dy_5) */
/* column 7 (dF/dy_6) */
/* column 8 (dF/dy_7) */
/* column 9 (dF/dy_8) */
/* column 10 (dF/dy_9) */
/* column 11 (dF/dy_10) */
/* column 12 (dF/dy_11) */
/* column 13 (dF/dy_12) */
/* column 14 (dF/dy_13) */
/* column 15 (dF/dy_14) */
/* column 16 (dF/dy_15) */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int AKAP79tc_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 26;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	double Rii=y_[0];
	double cAMP=y_[1];
	double RiiP=y_[2];
	double Rii_C=y_[3];
	double RiiP_cAMP=y_[4];
	double RiiP_C=y_[5];
	double RiiP_C_cAMP=y_[6];
	double C=y_[7];
	double Rii_cAMP=y_[8];
	double Rii_C_cAMP=y_[9];
	double CaN=y_[10];
	double RiiP_CaN=y_[11];
	double RiiP_cAMP_CaN=y_[12];
	double AKAR4=y_[13];
	double AKAR4_C=y_[14];
	double AKAR4p=y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	double kf_RiixC__Rii_C=(kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kb_Rii_CxcAMP__Rii_C_cAMP)/(kb_cAMPxRii__Rii_cAMP*kb_Rii_cAMPxC__Rii_C_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP);
	double kf_RiiPXcAMP__RiiP_cAMP=(kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP *kb_RiiP_cAMPxC__RiiP_C_cAMP)/(kb_RiiPxC__RiiP_C*kb_RiiP_CxcAMP__RiiP_C_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP);
	double reaction_51=kf_Rii_C__RiiP_C*Rii_C;
	double reaction_14=kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	double reaction_12=kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_43=kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	double reaction_23=kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	double reaction_78=kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	double reaction_56=kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_76=kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	double reaction_62=kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	double reaction_58=kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (dF/dp_0) */
/* column 2 (dF/dp_1) */
/* column 3 (dF/dp_2) */
/* column 4 (dF/dp_3) */
/* column 5 (dF/dp_4) */
/* column 6 (dF/dp_5) */
/* column 7 (dF/dp_6) */
/* column 8 (dF/dp_7) */
/* column 9 (dF/dp_8) */
/* column 10 (dF/dp_9) */
/* column 11 (dF/dp_10) */
/* column 12 (dF/dp_11) */
/* column 13 (dF/dp_12) */
/* column 14 (dF/dp_13) */
/* column 15 (dF/dp_14) */
/* column 16 (dF/dp_15) */
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAP79tc_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 26;
	p_[0] = 33;
	p_[1] = 0.496;
	p_[2] = 0.00545;
	p_[3] = 0.0156;
	p_[4] = 0.0016;
	p_[5] = 0.038;
	p_[6] = 0.0026;
	p_[7] = 0.015;
	p_[8] = 0.0016;
	p_[9] = 0.496;
	p_[10] = 1.413;
	p_[11] = 0.2984;
	p_[12] = 0.018;
	p_[13] = 33;
	p_[14] = 0.0003;
	p_[15] = 2.6;
	p_[16] = 20;
	p_[17] = 0.45;
	p_[18] = 2;
	p_[19] = 0.018;
	p_[20] = 0.106;
	p_[21] = 10.2;
	p_[22] = 100;
	p_[23] = 1;
	p_[24] = 0.7;
	p_[25] = 0;
	return GSL_SUCCESS;
}
/* ode initial values */
int AKAP79tc_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return       16;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPxC__RiiP_C=p_[5];
	double kb_RiiPxC__RiiP_C=p_[6];
	double kf_cAMPxRii__Rii_cAMP=p_[7];
	double kb_cAMPxRii__Rii_cAMP=p_[8];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[9];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[11];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[12];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[13];
	double kb_RiixC__Rii_C=p_[14];
	double AKAPoff_1=p_[15];
	double AKAPoff_3=p_[16];
	double AKAPon_1=p_[17];
	double AKAPon_3=p_[18];
	double kf_C_AKAR4=p_[19];
	double kb_C_AKAR4=p_[20];
	double kcat_AKARp=p_[21];
	double kmOFF=p_[22];
	double kmON=p_[23];
	double KD_T=p_[24];
	double b_AKAP=p_[25];
	/* the initial value of y may depend on the parameters. */
	y_[0] = 6.3;
	y_[1] = 0;
	y_[2] = 0;
	y_[3] = 0.63;
	y_[4] = 0;
	y_[5] = 0;
	y_[6] = 0;
	y_[7] = 0;
	y_[8] = 0;
	y_[9] = 0;
	y_[10] = 1.5;
	y_[11] = 0;
	y_[12] = 0;
	y_[13] = 0.2;
	y_[14] = 0;
	y_[15] = 0;
	return GSL_SUCCESS;
}
