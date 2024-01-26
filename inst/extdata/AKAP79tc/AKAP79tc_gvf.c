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
<<<<<<< HEAD
	jac_[0] = (-(C*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP))-cAMP*kf_cAMPxRii__Rii_cAMP; /* [0, 0] */
	jac_[16] = -cAMP*kf_cAMPxRii__Rii_cAMP; /* [1, 0] */
	jac_[32] = 0; /* [2, 0] */
	jac_[48] = (C*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [3, 0] */
	jac_[64] = 0; /* [4, 0] */
	jac_[80] = 0; /* [5, 0] */
	jac_[96] = 0; /* [6, 0] */
	jac_[112] = -(C*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [7, 0] */
	jac_[128] = cAMP*kf_cAMPxRii__Rii_cAMP; /* [8, 0] */
	jac_[144] = 0; /* [9, 0] */
	jac_[160] = 0; /* [10, 0] */
	jac_[176] = 0; /* [11, 0] */
	jac_[192] = 0; /* [12, 0] */
	jac_[208] = 0; /* [13, 0] */
	jac_[224] = 0; /* [14, 0] */
	jac_[240] = 0; /* [15, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = -Rii*kf_cAMPxRii__Rii_cAMP; /* [0, 1] */
	jac_[17] = (-Rii*kf_cAMPxRii__Rii_cAMP)-Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP-(RiiP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP)-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [1, 1] */
	jac_[33] = -(RiiP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [2, 1] */
	jac_[49] = -Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP; /* [3, 1] */
	jac_[65] = (RiiP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [4, 1] */
	jac_[81] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [5, 1] */
	jac_[97] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [6, 1] */
	jac_[113] = 0; /* [7, 1] */
	jac_[129] = Rii*kf_cAMPxRii__Rii_cAMP; /* [8, 1] */
	jac_[145] = Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP; /* [9, 1] */
	jac_[161] = 0; /* [10, 1] */
	jac_[177] = 0; /* [11, 1] */
	jac_[193] = 0; /* [12, 1] */
	jac_[209] = 0; /* [13, 1] */
	jac_[225] = 0; /* [14, 1] */
	jac_[241] = 0; /* [15, 1] */
/* column 3 (df/dy_2) */
	jac_[2] = 0; /* [0, 2] */
	jac_[18] = -(cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [1, 2] */
	jac_[34] = (-CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-(cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP)-C*kf_RiiPxC__RiiP_C; /* [2, 2] */
	jac_[50] = 0; /* [3, 2] */
	jac_[66] = (cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [4, 2] */
	jac_[82] = C*kf_RiiPxC__RiiP_C; /* [5, 2] */
	jac_[98] = 0; /* [6, 2] */
	jac_[114] = -C*kf_RiiPxC__RiiP_C; /* [7, 2] */
	jac_[130] = 0; /* [8, 2] */
	jac_[146] = 0; /* [9, 2] */
	jac_[162] = -CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [10, 2] */
	jac_[178] = CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [11, 2] */
	jac_[194] = 0; /* [12, 2] */
	jac_[210] = 0; /* [13, 2] */
	jac_[226] = 0; /* [14, 2] */
	jac_[242] = 0; /* [15, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = kb_RiixC__Rii_C; /* [0, 3] */
	jac_[19] = -cAMP*kf_Rii_CxcAMP__Rii_C_cAMP; /* [1, 3] */
	jac_[35] = 0; /* [2, 3] */
	jac_[51] = (-cAMP*kf_Rii_CxcAMP__Rii_C_cAMP)-kf_Rii_C__RiiP_C-kb_RiixC__Rii_C; /* [3, 3] */
	jac_[67] = 0; /* [4, 3] */
	jac_[83] = kf_Rii_C__RiiP_C; /* [5, 3] */
	jac_[99] = 0; /* [6, 3] */
	jac_[115] = kb_RiixC__Rii_C; /* [7, 3] */
	jac_[131] = 0; /* [8, 3] */
	jac_[147] = cAMP*kf_Rii_CxcAMP__Rii_C_cAMP; /* [9, 3] */
	jac_[163] = 0; /* [10, 3] */
	jac_[179] = 0; /* [11, 3] */
	jac_[195] = 0; /* [12, 3] */
	jac_[211] = 0; /* [13, 3] */
	jac_[227] = 0; /* [14, 3] */
	jac_[243] = 0; /* [15, 3] */
/* column 5 (df/dy_4) */
	jac_[4] = 0; /* [0, 4] */
	jac_[20] = kb_RiiPXcAMP__RiiP_cAMP; /* [1, 4] */
	jac_[36] = kb_RiiPXcAMP__RiiP_cAMP; /* [2, 4] */
	jac_[52] = 0; /* [3, 4] */
	jac_[68] = (-CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-C*kf_RiiP_cAMPxC__RiiP_C_cAMP-kb_RiiPXcAMP__RiiP_cAMP; /* [4, 4] */
	jac_[84] = 0; /* [5, 4] */
	jac_[100] = C*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [6, 4] */
	jac_[116] = -C*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [7, 4] */
	jac_[132] = 0; /* [8, 4] */
	jac_[148] = 0; /* [9, 4] */
	jac_[164] = -CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [10, 4] */
	jac_[180] = 0; /* [11, 4] */
	jac_[196] = CaN*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [12, 4] */
	jac_[212] = 0; /* [13, 4] */
	jac_[228] = 0; /* [14, 4] */
	jac_[244] = 0; /* [15, 4] */
/* column 6 (df/dy_5) */
	jac_[5] = 0; /* [0, 5] */
	jac_[21] = -cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [1, 5] */
	jac_[37] = kb_RiiPxC__RiiP_C; /* [2, 5] */
	jac_[53] = 0; /* [3, 5] */
	jac_[69] = 0; /* [4, 5] */
	jac_[85] = (-cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP)-kb_RiiPxC__RiiP_C; /* [5, 5] */
	jac_[101] = cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [6, 5] */
	jac_[117] = kb_RiiPxC__RiiP_C; /* [7, 5] */
	jac_[133] = 0; /* [8, 5] */
	jac_[149] = 0; /* [9, 5] */
	jac_[165] = 0; /* [10, 5] */
	jac_[181] = 0; /* [11, 5] */
	jac_[197] = 0; /* [12, 5] */
	jac_[213] = 0; /* [13, 5] */
	jac_[229] = 0; /* [14, 5] */
	jac_[245] = 0; /* [15, 5] */
/* column 7 (df/dy_6) */
	jac_[6] = 0; /* [0, 6] */
	jac_[22] = KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [1, 6] */
	jac_[38] = 0; /* [2, 6] */
	jac_[54] = 0; /* [3, 6] */
	jac_[70] = kb_RiiP_cAMPxC__RiiP_C_cAMP; /* [4, 6] */
	jac_[86] = KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [5, 6] */
	jac_[102] = (-KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP)-kb_RiiP_cAMPxC__RiiP_C_cAMP; /* [6, 6] */
	jac_[118] = kb_RiiP_cAMPxC__RiiP_C_cAMP; /* [7, 6] */
	jac_[134] = 0; /* [8, 6] */
	jac_[150] = 0; /* [9, 6] */
	jac_[166] = 0; /* [10, 6] */
	jac_[182] = 0; /* [11, 6] */
	jac_[198] = 0; /* [12, 6] */
	jac_[214] = 0; /* [13, 6] */
	jac_[230] = 0; /* [14, 6] */
	jac_[246] = 0; /* [15, 6] */
/* column 8 (df/dy_7) */
	jac_[7] = -(Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [0, 7] */
	jac_[23] = 0; /* [1, 7] */
	jac_[39] = -RiiP*kf_RiiPxC__RiiP_C; /* [2, 7] */
	jac_[55] = (Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [3, 7] */
	jac_[71] = -RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [4, 7] */
	jac_[87] = RiiP*kf_RiiPxC__RiiP_C; /* [5, 7] */
	jac_[103] = RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [6, 7] */
	jac_[119] = (-(Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP))-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-RiiP*kf_RiiPxC__RiiP_C-RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP-AKAR4*kf_C_AKAR4; /* [7, 7] */
	jac_[135] = -Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP; /* [8, 7] */
	jac_[151] = Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP; /* [9, 7] */
	jac_[167] = 0; /* [10, 7] */
	jac_[183] = 0; /* [11, 7] */
	jac_[199] = 0; /* [12, 7] */
	jac_[215] = -AKAR4*kf_C_AKAR4; /* [13, 7] */
	jac_[231] = AKAR4*kf_C_AKAR4; /* [14, 7] */
	jac_[247] = 0; /* [15, 7] */
/* column 9 (df/dy_8) */
	jac_[8] = kb_cAMPxRii__Rii_cAMP; /* [0, 8] */
	jac_[24] = kb_cAMPxRii__Rii_cAMP; /* [1, 8] */
	jac_[40] = 0; /* [2, 8] */
	jac_[56] = 0; /* [3, 8] */
	jac_[72] = 0; /* [4, 8] */
	jac_[88] = 0; /* [5, 8] */
	jac_[104] = 0; /* [6, 8] */
	jac_[120] = -C*kf_Rii_cAMPxC__Rii_C_cAMP; /* [7, 8] */
	jac_[136] = (-C*kf_Rii_cAMPxC__Rii_C_cAMP)-kb_cAMPxRii__Rii_cAMP; /* [8, 8] */
	jac_[152] = C*kf_Rii_cAMPxC__Rii_C_cAMP; /* [9, 8] */
	jac_[168] = 0; /* [10, 8] */
	jac_[184] = 0; /* [11, 8] */
	jac_[200] = 0; /* [12, 8] */
	jac_[216] = 0; /* [13, 8] */
	jac_[232] = 0; /* [14, 8] */
	jac_[248] = 0; /* [15, 8] */
/* column 10 (df/dy_9) */
	jac_[9] = 0; /* [0, 9] */
	jac_[25] = kb_Rii_CxcAMP__Rii_C_cAMP; /* [1, 9] */
	jac_[41] = 0; /* [2, 9] */
	jac_[57] = kb_Rii_CxcAMP__Rii_C_cAMP; /* [3, 9] */
	jac_[73] = 0; /* [4, 9] */
	jac_[89] = 0; /* [5, 9] */
	jac_[105] = kf_Rii_C_cAMP__RiiP_C_cAMP; /* [6, 9] */
	jac_[121] = kb_Rii_cAMPxC__Rii_C_cAMP; /* [7, 9] */
	jac_[137] = kb_Rii_cAMPxC__Rii_C_cAMP; /* [8, 9] */
	jac_[153] = (-kf_Rii_C_cAMP__RiiP_C_cAMP)-kb_Rii_cAMPxC__Rii_C_cAMP-kb_Rii_CxcAMP__Rii_C_cAMP; /* [9, 9] */
	jac_[169] = 0; /* [10, 9] */
	jac_[185] = 0; /* [11, 9] */
	jac_[201] = 0; /* [12, 9] */
	jac_[217] = 0; /* [13, 9] */
	jac_[233] = 0; /* [14, 9] */
	jac_[249] = 0; /* [15, 9] */
/* column 11 (df/dy_10) */
	jac_[10] = 0; /* [0, 10] */
	jac_[26] = 0; /* [1, 10] */
	jac_[42] = -RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [2, 10] */
	jac_[58] = 0; /* [3, 10] */
	jac_[74] = -RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [4, 10] */
	jac_[90] = 0; /* [5, 10] */
	jac_[106] = 0; /* [6, 10] */
	jac_[122] = 0; /* [7, 10] */
	jac_[138] = 0; /* [8, 10] */
	jac_[154] = 0; /* [9, 10] */
	jac_[170] = (-RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [10, 10] */
	jac_[186] = RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [11, 10] */
	jac_[202] = RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [12, 10] */
	jac_[218] = 0; /* [13, 10] */
	jac_[234] = 0; /* [14, 10] */
	jac_[250] = 0; /* [15, 10] */
/* column 12 (df/dy_11) */
	jac_[11] = AKAPon_1*b_AKAP+AKAPoff_1*(1-b_AKAP); /* [0, 11] */
	jac_[27] = 0; /* [1, 11] */
	jac_[43] = AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP); /* [2, 11] */
	jac_[59] = 0; /* [3, 11] */
	jac_[75] = 0; /* [4, 11] */
	jac_[91] = 0; /* [5, 11] */
	jac_[107] = 0; /* [6, 11] */
	jac_[123] = 0; /* [7, 11] */
	jac_[139] = 0; /* [8, 11] */
	jac_[155] = 0; /* [9, 11] */
	jac_[171] = AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP); /* [10, 11] */
	jac_[187] = (-AKAPon_3*b_AKAP)-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP); /* [11, 11] */
	jac_[203] = 0; /* [12, 11] */
	jac_[219] = 0; /* [13, 11] */
	jac_[235] = 0; /* [14, 11] */
	jac_[251] = 0; /* [15, 11] */
/* column 13 (df/dy_12) */
	jac_[12] = 0; /* [0, 12] */
	jac_[28] = 0; /* [1, 12] */
	jac_[44] = 0; /* [2, 12] */
	jac_[60] = 0; /* [3, 12] */
	jac_[76] = AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP); /* [4, 12] */
	jac_[92] = 0; /* [5, 12] */
	jac_[108] = 0; /* [6, 12] */
	jac_[124] = 0; /* [7, 12] */
	jac_[140] = AKAPon_1*b_AKAP+AKAPoff_1*(1-b_AKAP); /* [8, 12] */
	jac_[156] = 0; /* [9, 12] */
	jac_[172] = AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP); /* [10, 12] */
	jac_[188] = 0; /* [11, 12] */
	jac_[204] = (-AKAPon_3*b_AKAP)-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP); /* [12, 12] */
	jac_[220] = 0; /* [13, 12] */
	jac_[236] = 0; /* [14, 12] */
	jac_[252] = 0; /* [15, 12] */
/* column 14 (df/dy_13) */
	jac_[13] = 0; /* [0, 13] */
	jac_[29] = 0; /* [1, 13] */
	jac_[45] = 0; /* [2, 13] */
	jac_[61] = 0; /* [3, 13] */
	jac_[77] = 0; /* [4, 13] */
	jac_[93] = 0; /* [5, 13] */
	jac_[109] = 0; /* [6, 13] */
	jac_[125] = -C*kf_C_AKAR4; /* [7, 13] */
	jac_[141] = 0; /* [8, 13] */
	jac_[157] = 0; /* [9, 13] */
	jac_[173] = 0; /* [10, 13] */
	jac_[189] = 0; /* [11, 13] */
	jac_[205] = 0; /* [12, 13] */
	jac_[221] = -C*kf_C_AKAR4; /* [13, 13] */
	jac_[237] = C*kf_C_AKAR4; /* [14, 13] */
	jac_[253] = 0; /* [15, 13] */
/* column 15 (df/dy_14) */
	jac_[14] = 0; /* [0, 14] */
	jac_[30] = 0; /* [1, 14] */
	jac_[46] = 0; /* [2, 14] */
	jac_[62] = 0; /* [3, 14] */
	jac_[78] = 0; /* [4, 14] */
	jac_[94] = 0; /* [5, 14] */
	jac_[110] = 0; /* [6, 14] */
	jac_[126] = kcat_AKARp+kb_C_AKAR4; /* [7, 14] */
	jac_[142] = 0; /* [8, 14] */
	jac_[158] = 0; /* [9, 14] */
	jac_[174] = 0; /* [10, 14] */
	jac_[190] = 0; /* [11, 14] */
	jac_[206] = 0; /* [12, 14] */
	jac_[222] = kb_C_AKAR4; /* [13, 14] */
	jac_[238] = (-kcat_AKARp)-kb_C_AKAR4; /* [14, 14] */
	jac_[254] = kcat_AKARp; /* [15, 14] */
/* column 16 (df/dy_15) */
	jac_[15] = 0; /* [0, 15] */
	jac_[31] = 0; /* [1, 15] */
	jac_[47] = 0; /* [2, 15] */
	jac_[63] = 0; /* [3, 15] */
	jac_[79] = 0; /* [4, 15] */
	jac_[95] = 0; /* [5, 15] */
	jac_[111] = 0; /* [6, 15] */
	jac_[127] = 0; /* [7, 15] */
	jac_[143] = 0; /* [8, 15] */
	jac_[159] = 0; /* [9, 15] */
	jac_[175] = 0; /* [10, 15] */
	jac_[191] = 0; /* [11, 15] */
	jac_[207] = 0; /* [12, 15] */
	jac_[223] = 0; /* [13, 15] */
	jac_[239] = 0; /* [14, 15] */
	jac_[255] = 0; /* [15, 15] */
=======
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
>>>>>>> d25a8b388c0762c5d23c18e4194eb1c284c0c75a
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
<<<<<<< HEAD
	jacp_[0] = 0; /* [0, 0] */
	jacp_[26] = 0; /* [1, 0] */
	jacp_[52] = 0; /* [2, 0] */
	jacp_[78] = -Rii_C; /* [3, 0] */
	jacp_[104] = 0; /* [4, 0] */
	jacp_[130] = Rii_C; /* [5, 0] */
	jacp_[156] = 0; /* [6, 0] */
	jacp_[182] = 0; /* [7, 0] */
	jacp_[208] = 0; /* [8, 0] */
	jacp_[234] = 0; /* [9, 0] */
	jacp_[260] = 0; /* [10, 0] */
	jacp_[286] = 0; /* [11, 0] */
	jacp_[312] = 0; /* [12, 0] */
	jacp_[338] = 0; /* [13, 0] */
	jacp_[364] = 0; /* [14, 0] */
	jacp_[390] = 0; /* [15, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = 0; /* [0, 1] */
	jacp_[27] = KD_T*RiiP_C_cAMP-RiiP_C*cAMP; /* [1, 1] */
	jacp_[53] = 0; /* [2, 1] */
	jacp_[79] = 0; /* [3, 1] */
	jacp_[105] = 0; /* [4, 1] */
	jacp_[131] = KD_T*RiiP_C_cAMP-RiiP_C*cAMP; /* [5, 1] */
	jacp_[157] = RiiP_C*cAMP-KD_T*RiiP_C_cAMP; /* [6, 1] */
	jacp_[183] = 0; /* [7, 1] */
	jacp_[209] = 0; /* [8, 1] */
	jacp_[235] = 0; /* [9, 1] */
	jacp_[261] = 0; /* [10, 1] */
	jacp_[287] = 0; /* [11, 1] */
	jacp_[313] = 0; /* [12, 1] */
	jacp_[339] = 0; /* [13, 1] */
	jacp_[365] = 0; /* [14, 1] */
	jacp_[391] = 0; /* [15, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 0; /* [0, 2] */
	jacp_[28] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*gsl_pow_2(kf_RiiP_cAMPxC__RiiP_C_cAMP)); /* [1, 2] */
	jacp_[54] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*gsl_pow_2(kf_RiiP_cAMPxC__RiiP_C_cAMP)); /* [2, 2] */
	jacp_[80] = 0; /* [3, 2] */
	jacp_[106] = (-(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*gsl_pow_2(kf_RiiP_cAMPxC__RiiP_C_cAMP)))-C*RiiP_cAMP; /* [4, 2] */
	jacp_[132] = 0; /* [5, 2] */
	jacp_[158] = C*RiiP_cAMP; /* [6, 2] */
	jacp_[184] = -C*RiiP_cAMP; /* [7, 2] */
	jacp_[210] = 0; /* [8, 2] */
	jacp_[236] = 0; /* [9, 2] */
	jacp_[262] = 0; /* [10, 2] */
	jacp_[288] = 0; /* [11, 2] */
	jacp_[314] = 0; /* [12, 2] */
	jacp_[340] = 0; /* [13, 2] */
	jacp_[366] = 0; /* [14, 2] */
	jacp_[392] = 0; /* [15, 2] */
/* column 4 (df/dp_3) */
	jacp_[3] = 0; /* [0, 3] */
	jacp_[29] = -(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [1, 3] */
	jacp_[55] = -(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [2, 3] */
	jacp_[81] = 0; /* [3, 3] */
	jacp_[107] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP)+RiiP_C_cAMP; /* [4, 3] */
	jacp_[133] = 0; /* [5, 3] */
	jacp_[159] = -RiiP_C_cAMP; /* [6, 3] */
	jacp_[185] = RiiP_C_cAMP; /* [7, 3] */
	jacp_[211] = 0; /* [8, 3] */
	jacp_[237] = 0; /* [9, 3] */
	jacp_[263] = 0; /* [10, 3] */
	jacp_[289] = 0; /* [11, 3] */
	jacp_[315] = 0; /* [12, 3] */
	jacp_[341] = 0; /* [13, 3] */
	jacp_[367] = 0; /* [14, 3] */
	jacp_[393] = 0; /* [15, 3] */
/* column 5 (df/dp_4) */
	jacp_[4] = 0; /* [0, 4] */
	jacp_[30] = RiiP_cAMP-(RiiP*cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [1, 4] */
	jacp_[56] = RiiP_cAMP-(RiiP*cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [2, 4] */
	jacp_[82] = 0; /* [3, 4] */
	jacp_[108] = (RiiP*cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP)-RiiP_cAMP; /* [4, 4] */
	jacp_[134] = 0; /* [5, 4] */
	jacp_[160] = 0; /* [6, 4] */
	jacp_[186] = 0; /* [7, 4] */
	jacp_[212] = 0; /* [8, 4] */
	jacp_[238] = 0; /* [9, 4] */
	jacp_[264] = 0; /* [10, 4] */
	jacp_[290] = 0; /* [11, 4] */
	jacp_[316] = 0; /* [12, 4] */
	jacp_[342] = 0; /* [13, 4] */
	jacp_[368] = 0; /* [14, 4] */
	jacp_[394] = 0; /* [15, 4] */
/* column 6 (df/dp_5) */
	jacp_[5] = 0; /* [0, 5] */
	jacp_[31] = -(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [1, 5] */
	jacp_[57] = (-(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP))-C*RiiP; /* [2, 5] */
	jacp_[83] = 0; /* [3, 5] */
	jacp_[109] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP)/(KD_T*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [4, 5] */
	jacp_[135] = C*RiiP; /* [5, 5] */
	jacp_[161] = 0; /* [6, 5] */
	jacp_[187] = -C*RiiP; /* [7, 5] */
	jacp_[213] = 0; /* [8, 5] */
	jacp_[239] = 0; /* [9, 5] */
	jacp_[265] = 0; /* [10, 5] */
	jacp_[291] = 0; /* [11, 5] */
	jacp_[317] = 0; /* [12, 5] */
	jacp_[343] = 0; /* [13, 5] */
	jacp_[369] = 0; /* [14, 5] */
	jacp_[395] = 0; /* [15, 5] */
/* column 7 (df/dp_6) */
	jacp_[6] = 0; /* [0, 6] */
	jacp_[32] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*gsl_pow_2(kb_RiiPxC__RiiP_C)*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [1, 6] */
	jacp_[58] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*gsl_pow_2(kb_RiiPxC__RiiP_C)*kf_RiiP_cAMPxC__RiiP_C_cAMP)+RiiP_C; /* [2, 6] */
	jacp_[84] = 0; /* [3, 6] */
	jacp_[110] = -(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(KD_T*gsl_pow_2(kb_RiiPxC__RiiP_C)*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [4, 6] */
	jacp_[136] = -RiiP_C; /* [5, 6] */
	jacp_[162] = 0; /* [6, 6] */
	jacp_[188] = RiiP_C; /* [7, 6] */
	jacp_[214] = 0; /* [8, 6] */
	jacp_[240] = 0; /* [9, 6] */
	jacp_[266] = 0; /* [10, 6] */
	jacp_[292] = 0; /* [11, 6] */
	jacp_[318] = 0; /* [12, 6] */
	jacp_[344] = 0; /* [13, 6] */
	jacp_[370] = 0; /* [14, 6] */
	jacp_[396] = 0; /* [15, 6] */
/* column 8 (df/dp_7) */
	jacp_[7] = (-(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP))-Rii*cAMP; /* [0, 7] */
	jacp_[33] = -Rii*cAMP; /* [1, 7] */
	jacp_[59] = 0; /* [2, 7] */
	jacp_[85] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [3, 7] */
	jacp_[111] = 0; /* [4, 7] */
	jacp_[137] = 0; /* [5, 7] */
	jacp_[163] = 0; /* [6, 7] */
	jacp_[189] = -(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [7, 7] */
	jacp_[215] = Rii*cAMP; /* [8, 7] */
	jacp_[241] = 0; /* [9, 7] */
	jacp_[267] = 0; /* [10, 7] */
	jacp_[293] = 0; /* [11, 7] */
	jacp_[319] = 0; /* [12, 7] */
	jacp_[345] = 0; /* [13, 7] */
	jacp_[371] = 0; /* [14, 7] */
	jacp_[397] = 0; /* [15, 7] */
/* column 9 (df/dp_8) */
	jacp_[8] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*gsl_pow_2(kb_cAMPxRii__Rii_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP)+Rii_cAMP; /* [0, 8] */
	jacp_[34] = Rii_cAMP; /* [1, 8] */
	jacp_[60] = 0; /* [2, 8] */
	jacp_[86] = -(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*gsl_pow_2(kb_cAMPxRii__Rii_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP); /* [3, 8] */
	jacp_[112] = 0; /* [4, 8] */
	jacp_[138] = 0; /* [5, 8] */
	jacp_[164] = 0; /* [6, 8] */
	jacp_[190] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*gsl_pow_2(kb_cAMPxRii__Rii_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP); /* [7, 8] */
	jacp_[216] = -Rii_cAMP; /* [8, 8] */
	jacp_[242] = 0; /* [9, 8] */
	jacp_[268] = 0; /* [10, 8] */
	jacp_[294] = 0; /* [11, 8] */
	jacp_[320] = 0; /* [12, 8] */
	jacp_[346] = 0; /* [13, 8] */
	jacp_[372] = 0; /* [14, 8] */
	jacp_[398] = 0; /* [15, 8] */
/* column 10 (df/dp_9) */
	jacp_[9] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*gsl_pow_2(kf_Rii_CxcAMP__Rii_C_cAMP)); /* [0, 9] */
	jacp_[35] = -Rii_C*cAMP; /* [1, 9] */
	jacp_[61] = 0; /* [2, 9] */
	jacp_[87] = (-(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*gsl_pow_2(kf_Rii_CxcAMP__Rii_C_cAMP)))-Rii_C*cAMP; /* [3, 9] */
	jacp_[113] = 0; /* [4, 9] */
	jacp_[139] = 0; /* [5, 9] */
	jacp_[165] = 0; /* [6, 9] */
	jacp_[191] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*gsl_pow_2(kf_Rii_CxcAMP__Rii_C_cAMP)); /* [7, 9] */
	jacp_[217] = 0; /* [8, 9] */
	jacp_[243] = Rii_C*cAMP; /* [9, 9] */
	jacp_[269] = 0; /* [10, 9] */
	jacp_[295] = 0; /* [11, 9] */
	jacp_[321] = 0; /* [12, 9] */
	jacp_[347] = 0; /* [13, 9] */
	jacp_[373] = 0; /* [14, 9] */
	jacp_[399] = 0; /* [15, 9] */
/* column 11 (df/dp_10) */
	jacp_[10] = -(C*Rii*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [0, 10] */
	jacp_[36] = Rii_C_cAMP; /* [1, 10] */
	jacp_[62] = 0; /* [2, 10] */
	jacp_[88] = (C*Rii*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP)+Rii_C_cAMP; /* [3, 10] */
	jacp_[114] = 0; /* [4, 10] */
	jacp_[140] = 0; /* [5, 10] */
	jacp_[166] = 0; /* [6, 10] */
	jacp_[192] = -(C*Rii*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [7, 10] */
	jacp_[218] = 0; /* [8, 10] */
	jacp_[244] = -Rii_C_cAMP; /* [9, 10] */
	jacp_[270] = 0; /* [10, 10] */
	jacp_[296] = 0; /* [11, 10] */
	jacp_[322] = 0; /* [12, 10] */
	jacp_[348] = 0; /* [13, 10] */
	jacp_[374] = 0; /* [14, 10] */
	jacp_[400] = 0; /* [15, 10] */
/* column 12 (df/dp_11) */
	jacp_[11] = -(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [0, 11] */
	jacp_[37] = 0; /* [1, 11] */
	jacp_[63] = 0; /* [2, 11] */
	jacp_[89] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [3, 11] */
	jacp_[115] = 0; /* [4, 11] */
	jacp_[141] = 0; /* [5, 11] */
	jacp_[167] = 0; /* [6, 11] */
	jacp_[193] = (-(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP))-C*Rii_cAMP; /* [7, 11] */
	jacp_[219] = -C*Rii_cAMP; /* [8, 11] */
	jacp_[245] = C*Rii_cAMP; /* [9, 11] */
	jacp_[271] = 0; /* [10, 11] */
	jacp_[297] = 0; /* [11, 11] */
	jacp_[323] = 0; /* [12, 11] */
	jacp_[349] = 0; /* [13, 11] */
	jacp_[375] = 0; /* [14, 11] */
	jacp_[401] = 0; /* [15, 11] */
/* column 13 (df/dp_12) */
	jacp_[12] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(gsl_pow_2(kb_Rii_cAMPxC__Rii_C_cAMP)*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [0, 12] */
	jacp_[38] = 0; /* [1, 12] */
	jacp_[64] = 0; /* [2, 12] */
	jacp_[90] = -(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(gsl_pow_2(kb_Rii_cAMPxC__Rii_C_cAMP)*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [3, 12] */
	jacp_[116] = 0; /* [4, 12] */
	jacp_[142] = 0; /* [5, 12] */
	jacp_[168] = 0; /* [6, 12] */
	jacp_[194] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kb_RiixC__Rii_C*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(gsl_pow_2(kb_Rii_cAMPxC__Rii_C_cAMP)*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP)+Rii_C_cAMP; /* [7, 12] */
	jacp_[220] = Rii_C_cAMP; /* [8, 12] */
	jacp_[246] = -Rii_C_cAMP; /* [9, 12] */
	jacp_[272] = 0; /* [10, 12] */
	jacp_[298] = 0; /* [11, 12] */
	jacp_[324] = 0; /* [12, 12] */
	jacp_[350] = 0; /* [13, 12] */
	jacp_[376] = 0; /* [14, 12] */
	jacp_[402] = 0; /* [15, 12] */
/* column 14 (df/dp_13) */
	jacp_[13] = 0; /* [0, 13] */
	jacp_[39] = 0; /* [1, 13] */
	jacp_[65] = 0; /* [2, 13] */
	jacp_[91] = 0; /* [3, 13] */
	jacp_[117] = 0; /* [4, 13] */
	jacp_[143] = 0; /* [5, 13] */
	jacp_[169] = Rii_C_cAMP; /* [6, 13] */
	jacp_[195] = 0; /* [7, 13] */
	jacp_[221] = 0; /* [8, 13] */
	jacp_[247] = -Rii_C_cAMP; /* [9, 13] */
	jacp_[273] = 0; /* [10, 13] */
	jacp_[299] = 0; /* [11, 13] */
	jacp_[325] = 0; /* [12, 13] */
	jacp_[351] = 0; /* [13, 13] */
	jacp_[377] = 0; /* [14, 13] */
	jacp_[403] = 0; /* [15, 13] */
/* column 15 (df/dp_14) */
	jacp_[14] = Rii_C-(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [0, 14] */
	jacp_[40] = 0; /* [1, 14] */
	jacp_[66] = 0; /* [2, 14] */
	jacp_[92] = (C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP)-Rii_C; /* [3, 14] */
	jacp_[118] = 0; /* [4, 14] */
	jacp_[144] = 0; /* [5, 14] */
	jacp_[170] = 0; /* [6, 14] */
	jacp_[196] = Rii_C-(C*Rii*kb_Rii_CxcAMP__Rii_C_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP*kf_cAMPxRii__Rii_cAMP)/(kb_Rii_cAMPxC__Rii_C_cAMP*kb_cAMPxRii__Rii_cAMP*kf_Rii_CxcAMP__Rii_C_cAMP); /* [7, 14] */
	jacp_[222] = 0; /* [8, 14] */
	jacp_[248] = 0; /* [9, 14] */
	jacp_[274] = 0; /* [10, 14] */
	jacp_[300] = 0; /* [11, 14] */
	jacp_[326] = 0; /* [12, 14] */
	jacp_[352] = 0; /* [13, 14] */
	jacp_[378] = 0; /* [14, 14] */
	jacp_[404] = 0; /* [15, 14] */
/* column 16 (df/dp_15) */
	jacp_[15] = RiiP_CaN*(1-b_AKAP); /* [0, 15] */
	jacp_[41] = 0; /* [1, 15] */
	jacp_[67] = -CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [2, 15] */
	jacp_[93] = 0; /* [3, 15] */
	jacp_[119] = -CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [4, 15] */
	jacp_[145] = 0; /* [5, 15] */
	jacp_[171] = 0; /* [6, 15] */
	jacp_[197] = 0; /* [7, 15] */
	jacp_[223] = RiiP_cAMP_CaN*(1-b_AKAP); /* [8, 15] */
	jacp_[249] = 0; /* [9, 15] */
	jacp_[275] = (-CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF))-CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)+RiiP_cAMP_CaN*(1-b_AKAP)+RiiP_CaN*(1-b_AKAP); /* [10, 15] */
	jacp_[301] = CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_CaN*(1-b_AKAP); /* [11, 15] */
	jacp_[327] = CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP); /* [12, 15] */
	jacp_[353] = 0; /* [13, 15] */
	jacp_[379] = 0; /* [14, 15] */
	jacp_[405] = 0; /* [15, 15] */
/* column 17 (df/dp_16) */
	jacp_[16] = 0; /* [0, 16] */
	jacp_[42] = 0; /* [1, 16] */
	jacp_[68] = RiiP_CaN*(1-b_AKAP)-CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [2, 16] */
	jacp_[94] = 0; /* [3, 16] */
	jacp_[120] = RiiP_cAMP_CaN*(1-b_AKAP)-CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [4, 16] */
	jacp_[146] = 0; /* [5, 16] */
	jacp_[172] = 0; /* [6, 16] */
	jacp_[198] = 0; /* [7, 16] */
	jacp_[224] = 0; /* [8, 16] */
	jacp_[250] = 0; /* [9, 16] */
	jacp_[276] = (-CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF))-CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)+RiiP_cAMP_CaN*(1-b_AKAP)+RiiP_CaN*(1-b_AKAP); /* [10, 16] */
	jacp_[302] = CaN*RiiP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_CaN*(1-b_AKAP); /* [11, 16] */
	jacp_[328] = CaN*RiiP_cAMP*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP); /* [12, 16] */
	jacp_[354] = 0; /* [13, 16] */
	jacp_[380] = 0; /* [14, 16] */
	jacp_[406] = 0; /* [15, 16] */
/* column 18 (df/dp_17) */
	jacp_[17] = RiiP_CaN*b_AKAP; /* [0, 17] */
	jacp_[43] = 0; /* [1, 17] */
	jacp_[69] = -CaN*RiiP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [2, 17] */
	jacp_[95] = 0; /* [3, 17] */
	jacp_[121] = -CaN*RiiP_cAMP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [4, 17] */
	jacp_[147] = 0; /* [5, 17] */
	jacp_[173] = 0; /* [6, 17] */
	jacp_[199] = 0; /* [7, 17] */
	jacp_[225] = RiiP_cAMP_CaN*b_AKAP; /* [8, 17] */
	jacp_[251] = 0; /* [9, 17] */
	jacp_[277] = (-CaN*RiiP_cAMP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF))-CaN*RiiP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)+RiiP_cAMP_CaN*b_AKAP+RiiP_CaN*b_AKAP; /* [10, 17] */
	jacp_[303] = CaN*RiiP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP; /* [11, 17] */
	jacp_[329] = CaN*RiiP_cAMP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP; /* [12, 17] */
	jacp_[355] = 0; /* [13, 17] */
	jacp_[381] = 0; /* [14, 17] */
	jacp_[407] = 0; /* [15, 17] */
/* column 19 (df/dp_18) */
	jacp_[18] = 0; /* [0, 18] */
	jacp_[44] = 0; /* [1, 18] */
	jacp_[70] = RiiP_CaN*b_AKAP-CaN*RiiP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [2, 18] */
	jacp_[96] = 0; /* [3, 18] */
	jacp_[122] = RiiP_cAMP_CaN*b_AKAP-CaN*RiiP_cAMP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [4, 18] */
	jacp_[148] = 0; /* [5, 18] */
	jacp_[174] = 0; /* [6, 18] */
	jacp_[200] = 0; /* [7, 18] */
	jacp_[226] = 0; /* [8, 18] */
	jacp_[252] = 0; /* [9, 18] */
	jacp_[278] = (-CaN*RiiP_cAMP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF))-CaN*RiiP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)+RiiP_cAMP_CaN*b_AKAP+RiiP_CaN*b_AKAP; /* [10, 18] */
	jacp_[304] = CaN*RiiP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP; /* [11, 18] */
	jacp_[330] = CaN*RiiP_cAMP*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP; /* [12, 18] */
	jacp_[356] = 0; /* [13, 18] */
	jacp_[382] = 0; /* [14, 18] */
	jacp_[408] = 0; /* [15, 18] */
/* column 20 (df/dp_19) */
	jacp_[19] = 0; /* [0, 19] */
	jacp_[45] = 0; /* [1, 19] */
	jacp_[71] = 0; /* [2, 19] */
	jacp_[97] = 0; /* [3, 19] */
	jacp_[123] = 0; /* [4, 19] */
	jacp_[149] = 0; /* [5, 19] */
	jacp_[175] = 0; /* [6, 19] */
	jacp_[201] = -AKAR4*C; /* [7, 19] */
	jacp_[227] = 0; /* [8, 19] */
	jacp_[253] = 0; /* [9, 19] */
	jacp_[279] = 0; /* [10, 19] */
	jacp_[305] = 0; /* [11, 19] */
	jacp_[331] = 0; /* [12, 19] */
	jacp_[357] = -AKAR4*C; /* [13, 19] */
	jacp_[383] = AKAR4*C; /* [14, 19] */
	jacp_[409] = 0; /* [15, 19] */
/* column 21 (df/dp_20) */
	jacp_[20] = 0; /* [0, 20] */
	jacp_[46] = 0; /* [1, 20] */
	jacp_[72] = 0; /* [2, 20] */
	jacp_[98] = 0; /* [3, 20] */
	jacp_[124] = 0; /* [4, 20] */
	jacp_[150] = 0; /* [5, 20] */
	jacp_[176] = 0; /* [6, 20] */
	jacp_[202] = AKAR4_C; /* [7, 20] */
	jacp_[228] = 0; /* [8, 20] */
	jacp_[254] = 0; /* [9, 20] */
	jacp_[280] = 0; /* [10, 20] */
	jacp_[306] = 0; /* [11, 20] */
	jacp_[332] = 0; /* [12, 20] */
	jacp_[358] = AKAR4_C; /* [13, 20] */
	jacp_[384] = -AKAR4_C; /* [14, 20] */
	jacp_[410] = 0; /* [15, 20] */
/* column 22 (df/dp_21) */
	jacp_[21] = 0; /* [0, 21] */
	jacp_[47] = 0; /* [1, 21] */
	jacp_[73] = 0; /* [2, 21] */
	jacp_[99] = 0; /* [3, 21] */
	jacp_[125] = 0; /* [4, 21] */
	jacp_[151] = 0; /* [5, 21] */
	jacp_[177] = 0; /* [6, 21] */
	jacp_[203] = AKAR4_C; /* [7, 21] */
	jacp_[229] = 0; /* [8, 21] */
	jacp_[255] = 0; /* [9, 21] */
	jacp_[281] = 0; /* [10, 21] */
	jacp_[307] = 0; /* [11, 21] */
	jacp_[333] = 0; /* [12, 21] */
	jacp_[359] = 0; /* [13, 21] */
	jacp_[385] = -AKAR4_C; /* [14, 21] */
	jacp_[411] = AKAR4_C; /* [15, 21] */
/* column 23 (df/dp_22) */
	jacp_[22] = 0; /* [0, 22] */
	jacp_[48] = 0; /* [1, 22] */
	jacp_[74] = (CaN*RiiP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [2, 22] */
	jacp_[100] = 0; /* [3, 22] */
	jacp_[126] = (CaN*RiiP_cAMP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [4, 22] */
	jacp_[152] = 0; /* [5, 22] */
	jacp_[178] = 0; /* [6, 22] */
	jacp_[204] = 0; /* [7, 22] */
	jacp_[230] = 0; /* [8, 22] */
	jacp_[256] = 0; /* [9, 22] */
	jacp_[282] = (CaN*RiiP_cAMP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF)+(CaN*RiiP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [10, 22] */
	jacp_[308] = -(CaN*RiiP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [11, 22] */
	jacp_[334] = -(CaN*RiiP_cAMP*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [12, 22] */
	jacp_[360] = 0; /* [13, 22] */
	jacp_[386] = 0; /* [14, 22] */
	jacp_[412] = 0; /* [15, 22] */
/* column 24 (df/dp_23) */
	jacp_[23] = 0; /* [0, 23] */
	jacp_[49] = 0; /* [1, 23] */
	jacp_[75] = (CaN*RiiP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [2, 23] */
	jacp_[101] = 0; /* [3, 23] */
	jacp_[127] = (CaN*RiiP_cAMP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [4, 23] */
	jacp_[153] = 0; /* [5, 23] */
	jacp_[179] = 0; /* [6, 23] */
	jacp_[205] = 0; /* [7, 23] */
	jacp_[231] = 0; /* [8, 23] */
	jacp_[257] = 0; /* [9, 23] */
	jacp_[283] = (CaN*RiiP_cAMP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON)+(CaN*RiiP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [10, 23] */
	jacp_[309] = -(CaN*RiiP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [11, 23] */
	jacp_[335] = -(CaN*RiiP_cAMP*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [12, 23] */
	jacp_[361] = 0; /* [13, 23] */
	jacp_[387] = 0; /* [14, 23] */
	jacp_[413] = 0; /* [15, 23] */
/* column 25 (df/dp_24) */
	jacp_[24] = 0; /* [0, 24] */
	jacp_[50] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(gsl_pow_2(KD_T)*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP)+RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [1, 24] */
	jacp_[76] = (RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(gsl_pow_2(KD_T)*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [2, 24] */
	jacp_[102] = 0; /* [3, 24] */
	jacp_[128] = -(RiiP*cAMP*kb_RiiPXcAMP__RiiP_cAMP*kb_RiiP_cAMPxC__RiiP_C_cAMP*kf_RiiPxC__RiiP_C)/(gsl_pow_2(KD_T)*kb_RiiPxC__RiiP_C*kf_RiiP_cAMPxC__RiiP_C_cAMP); /* [4, 24] */
	jacp_[154] = RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [5, 24] */
	jacp_[180] = -RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [6, 24] */
	jacp_[206] = 0; /* [7, 24] */
	jacp_[232] = 0; /* [8, 24] */
	jacp_[258] = 0; /* [9, 24] */
	jacp_[284] = 0; /* [10, 24] */
	jacp_[310] = 0; /* [11, 24] */
	jacp_[336] = 0; /* [12, 24] */
	jacp_[362] = 0; /* [13, 24] */
	jacp_[388] = 0; /* [14, 24] */
	jacp_[414] = 0; /* [15, 24] */
/* column 26 (df/dp_25) */
	jacp_[25] = (AKAPon_1-AKAPoff_1)*RiiP_CaN; /* [0, 25] */
	jacp_[51] = 0; /* [1, 25] */
	jacp_[77] = (AKAPon_3-AKAPoff_3)*RiiP_CaN-CaN*RiiP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF); /* [2, 25] */
	jacp_[103] = 0; /* [3, 25] */
	jacp_[129] = (AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-CaN*RiiP_cAMP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF); /* [4, 25] */
	jacp_[155] = 0; /* [5, 25] */
	jacp_[181] = 0; /* [6, 25] */
	jacp_[207] = 0; /* [7, 25] */
	jacp_[233] = (AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN; /* [8, 25] */
	jacp_[259] = 0; /* [9, 25] */
	jacp_[285] = (-CaN*RiiP_cAMP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF))-CaN*RiiP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)+(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN+(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN+(AKAPon_3-AKAPoff_3)*RiiP_CaN+(AKAPon_1-AKAPoff_1)*RiiP_CaN; /* [10, 25] */
	jacp_[311] = CaN*RiiP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_CaN; /* [11, 25] */
	jacp_[337] = CaN*RiiP_cAMP*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN; /* [12, 25] */
	jacp_[363] = 0; /* [13, 25] */
	jacp_[389] = 0; /* [14, 25] */
	jacp_[415] = 0; /* [15, 25] */
=======
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
>>>>>>> d25a8b388c0762c5d23c18e4194eb1c284c0c75a
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
<<<<<<< HEAD
	funcJac_[0] = 0; /* [0, 0] */
/* column 2 (dF/dy_1) */
	funcJac_[1] = 0; /* [0, 1] */
/* column 3 (dF/dy_2) */
	funcJac_[2] = 0; /* [0, 2] */
/* column 4 (dF/dy_3) */
	funcJac_[3] = 0; /* [0, 3] */
/* column 5 (dF/dy_4) */
	funcJac_[4] = 0; /* [0, 4] */
/* column 6 (dF/dy_5) */
	funcJac_[5] = 0; /* [0, 5] */
/* column 7 (dF/dy_6) */
	funcJac_[6] = 0; /* [0, 6] */
/* column 8 (dF/dy_7) */
	funcJac_[7] = 0; /* [0, 7] */
/* column 9 (dF/dy_8) */
	funcJac_[8] = 0; /* [0, 8] */
/* column 10 (dF/dy_9) */
	funcJac_[9] = 0; /* [0, 9] */
/* column 11 (dF/dy_10) */
	funcJac_[10] = 0; /* [0, 10] */
/* column 12 (dF/dy_11) */
	funcJac_[11] = 0; /* [0, 11] */
/* column 13 (dF/dy_12) */
	funcJac_[12] = 0; /* [0, 12] */
/* column 14 (dF/dy_13) */
	funcJac_[13] = 0; /* [0, 13] */
/* column 15 (dF/dy_14) */
	funcJac_[14] = 0; /* [0, 14] */
/* column 16 (dF/dy_15) */
	funcJac_[15] = 358.35; /* [0, 15] */
=======
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
>>>>>>> d25a8b388c0762c5d23c18e4194eb1c284c0c75a
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
<<<<<<< HEAD
	funcJacp_[0] = 0; /* [0, 0] */
/* column 2 (dF/dp_1) */
	funcJacp_[1] = 0; /* [0, 1] */
/* column 3 (dF/dp_2) */
	funcJacp_[2] = 0; /* [0, 2] */
/* column 4 (dF/dp_3) */
	funcJacp_[3] = 0; /* [0, 3] */
/* column 5 (dF/dp_4) */
	funcJacp_[4] = 0; /* [0, 4] */
/* column 6 (dF/dp_5) */
	funcJacp_[5] = 0; /* [0, 5] */
/* column 7 (dF/dp_6) */
	funcJacp_[6] = 0; /* [0, 6] */
/* column 8 (dF/dp_7) */
	funcJacp_[7] = 0; /* [0, 7] */
/* column 9 (dF/dp_8) */
	funcJacp_[8] = 0; /* [0, 8] */
/* column 10 (dF/dp_9) */
	funcJacp_[9] = 0; /* [0, 9] */
/* column 11 (dF/dp_10) */
	funcJacp_[10] = 0; /* [0, 10] */
/* column 12 (dF/dp_11) */
	funcJacp_[11] = 0; /* [0, 11] */
/* column 13 (dF/dp_12) */
	funcJacp_[12] = 0; /* [0, 12] */
/* column 14 (dF/dp_13) */
	funcJacp_[13] = 0; /* [0, 13] */
/* column 15 (dF/dp_14) */
	funcJacp_[14] = 0; /* [0, 14] */
/* column 16 (dF/dp_15) */
	funcJacp_[15] = 0; /* [0, 15] */
=======
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
>>>>>>> d25a8b388c0762c5d23c18e4194eb1c284c0c75a
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
<<<<<<< HEAD
	if (!y_) return 16;
=======
	if (!y_) return       16;
>>>>>>> d25a8b388c0762c5d23c18e4194eb1c284c0c75a
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
