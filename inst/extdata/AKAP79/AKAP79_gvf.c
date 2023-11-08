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
int AKAP79_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 11;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	double RiiP=y_[0];
	double RiiP_cAMP=y_[1];
	double RiiP_C=y_[2];
	double RiiP_C_cAMP=y_[3];
	double C=y_[4];
	double Rii_cAMP=y_[5];
	double Rii_C_cAMP=y_[6];
	double RiiP_CaN=y_[7];
	double RiiP_cAMP_CaN=y_[8];
	double AKAR4_C=y_[9];
	double AKAR4p=y_[10];
	double AKAR4=AKAR4_ConservedConst - (AKAR4_C+AKAR4p);
	double CaN=CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN);
	double Rii_C=Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	double cAMP=cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN);
	double Rii=Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C);
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
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
	f_[0] = -reaction_14-reaction_43-reaction_44;
	f_[1] = +reaction_43-reaction_23-reaction_33;
	f_[2] = +reaction_51+reaction_14-reaction_12;
	f_[3] = +reaction_12+reaction_23+reaction_62;
	f_[4] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2;
	f_[5] = +reaction_78-reaction_76+reaction_37;
	f_[6] = +reaction_56+reaction_76-reaction_62;
	f_[7] = +reaction_44-reaction_48;
	f_[8] = +reaction_33-reaction_37;
	f_[9] = +reaction_1-reaction_2;
	f_[10] = +reaction_2;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int AKAP79_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 11*11;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	double RiiP=y_[0];
	double RiiP_cAMP=y_[1];
	double RiiP_C=y_[2];
	double RiiP_C_cAMP=y_[3];
	double C=y_[4];
	double Rii_cAMP=y_[5];
	double Rii_C_cAMP=y_[6];
	double RiiP_CaN=y_[7];
	double RiiP_cAMP_CaN=y_[8];
	double AKAR4_C=y_[9];
	double AKAR4p=y_[10];
	double AKAR4=AKAR4_ConservedConst - (AKAR4_C+AKAR4p);
	double CaN=CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN);
	double Rii_C=Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	double cAMP=cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN);
	double Rii=Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C);
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
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
	jac_[0] = (-((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-C*kf_RiiPxC__RiiP_C-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiPXcAMP__RiiP_cAMP; /* [0, 0] */
	jac_[11] = (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 0] */
	jac_[22] = C*kf_RiiPxC__RiiP_C; /* [2, 0] */
	jac_[33] = 0; /* [3, 0] */
	jac_[44] = C*kf_RiixC__Rii_C-C*kf_RiiPxC__RiiP_C; /* [4, 0] */
	jac_[55] = -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP; /* [5, 0] */
	jac_[66] = 0; /* [6, 0] */
	jac_[77] = ((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [7, 0] */
	jac_[88] = 0; /* [8, 0] */
	jac_[99] = 0; /* [9, 0] */
	jac_[110] = 0; /* [10, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = RiiP*kf_RiiPXcAMP__RiiP_cAMP+kb_RiiPXcAMP__RiiP_cAMP; /* [0, 1] */
	jac_[12] = (-((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-C*kf_RiiP_cAMPxC__RiiP_C_cAMP-RiiP*kf_RiiPXcAMP__RiiP_cAMP-kb_RiiPXcAMP__RiiP_cAMP; /* [1, 1] */
	jac_[23] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [2, 1] */
	jac_[34] = C*kf_RiiP_cAMPxC__RiiP_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 1] */
	jac_[45] = C*kf_RiixC__Rii_C-C*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [4, 1] */
	jac_[56] = (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP)-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP; /* [5, 1] */
	jac_[67] = -((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 1] */
	jac_[78] = 0; /* [7, 1] */
	jac_[89] = ((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [8, 1] */
	jac_[100] = 0; /* [9, 1] */
	jac_[111] = 0; /* [10, 1] */
/* column 3 (df/dy_2) */
	jac_[2] = kb_RiiPxC__RiiP_C; /* [0, 2] */
	jac_[13] = 0; /* [1, 2] */
	jac_[24] = (-kf_Rii_C__RiiP_C)-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiP_CxcAMP__RiiP_C_cAMP-kb_RiiPxC__RiiP_C; /* [2, 2] */
	jac_[35] = (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 2] */
	jac_[46] = kb_RiiPxC__RiiP_C-kb_RiixC__Rii_C; /* [4, 2] */
	jac_[57] = 0; /* [5, 2] */
	jac_[68] = -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 2] */
	jac_[79] = 0; /* [7, 2] */
	jac_[90] = 0; /* [8, 2] */
	jac_[101] = 0; /* [9, 2] */
	jac_[112] = 0; /* [10, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [0, 3] */
	jac_[14] = kb_RiiP_cAMPxC__RiiP_C_cAMP-RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 3] */
	jac_[25] = (-kf_Rii_C__RiiP_C)+RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP+KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [2, 3] */
	jac_[36] = (-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP)-KD_T*kf_RiiP_CxcAMP__RiiP_C_cAMP-kb_RiiP_cAMPxC__RiiP_C_cAMP; /* [3, 3] */
	jac_[47] = kb_RiiP_cAMPxC__RiiP_C_cAMP-kb_RiixC__Rii_C; /* [4, 3] */
	jac_[58] = -((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP; /* [5, 3] */
	jac_[69] = (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP)-((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 3] */
	jac_[80] = 0; /* [7, 3] */
	jac_[91] = 0; /* [8, 3] */
	jac_[102] = 0; /* [9, 3] */
	jac_[113] = 0; /* [10, 3] */
/* column 5 (df/dy_4) */
	jac_[4] = -RiiP*kf_RiiPxC__RiiP_C; /* [0, 4] */
	jac_[15] = -RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [1, 4] */
	jac_[26] = RiiP*kf_RiiPxC__RiiP_C-kf_Rii_C__RiiP_C; /* [2, 4] */
	jac_[37] = RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [3, 4] */
	jac_[48] = (-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_RiixC__Rii_C)-C*kf_RiixC__Rii_C-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-RiiP*kf_RiiPxC__RiiP_C-RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP-((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*kf_C_AKAR4-kb_RiixC__Rii_C; /* [4, 4] */
	jac_[59] = (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP; /* [5, 4] */
	jac_[70] = Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 4] */
	jac_[81] = 0; /* [7, 4] */
	jac_[92] = 0; /* [8, 4] */
	jac_[103] = ((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*kf_C_AKAR4; /* [9, 4] */
	jac_[114] = 0; /* [10, 4] */
/* column 6 (df/dy_5) */
	jac_[5] = RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [0, 5] */
	jac_[16] = -RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 5] */
	jac_[27] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [2, 5] */
	jac_[38] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 5] */
	jac_[49] = C*kf_RiixC__Rii_C-C*kf_Rii_cAMPxC__Rii_C_cAMP; /* [4, 5] */
	jac_[60] = (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP)-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP-C*kf_Rii_cAMPxC__Rii_C_cAMP-kb_cAMPxRii__Rii_cAMP; /* [5, 5] */
	jac_[71] = C*kf_Rii_cAMPxC__Rii_C_cAMP-((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 5] */
	jac_[82] = 0; /* [7, 5] */
	jac_[93] = 0; /* [8, 5] */
	jac_[104] = 0; /* [9, 5] */
	jac_[115] = 0; /* [10, 5] */
/* column 7 (df/dy_6) */
	jac_[6] = RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [0, 6] */
	jac_[17] = -RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 6] */
	jac_[28] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-kf_Rii_C__RiiP_C; /* [2, 6] */
	jac_[39] = kf_Rii_C_cAMP__RiiP_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 6] */
	jac_[50] = kb_Rii_cAMPxC__Rii_C_cAMP-kb_RiixC__Rii_C; /* [4, 6] */
	jac_[61] = kb_Rii_cAMPxC__Rii_C_cAMP-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP; /* [5, 6] */
	jac_[72] = (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP)-((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP-kf_Rii_C_cAMP__RiiP_C_cAMP-kb_Rii_cAMPxC__Rii_C_cAMP-kb_Rii_CxcAMP__Rii_C_cAMP; /* [6, 6] */
	jac_[83] = 0; /* [7, 6] */
	jac_[94] = 0; /* [8, 6] */
	jac_[105] = 0; /* [9, 6] */
	jac_[116] = 0; /* [10, 6] */
/* column 8 (df/dy_7) */
	jac_[7] = RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)+AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP); /* [0, 7] */
	jac_[18] = RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [1, 7] */
	jac_[29] = 0; /* [2, 7] */
	jac_[40] = 0; /* [3, 7] */
	jac_[51] = C*kf_RiixC__Rii_C; /* [4, 7] */
	jac_[62] = -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP; /* [5, 7] */
	jac_[73] = 0; /* [6, 7] */
	jac_[84] = (-RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-AKAPon_3*b_AKAP-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP); /* [7, 7] */
	jac_[95] = -RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [8, 7] */
	jac_[106] = 0; /* [9, 7] */
	jac_[117] = 0; /* [10, 7] */
/* column 9 (df/dy_8) */
	jac_[8] = RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)+RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [0, 8] */
	jac_[19] = RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF)-RiiP*kf_RiiPXcAMP__RiiP_cAMP+AKAPon_3*b_AKAP+AKAPoff_3*(1-b_AKAP); /* [1, 8] */
	jac_[30] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [2, 8] */
	jac_[41] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 8] */
	jac_[52] = C*kf_RiixC__Rii_C; /* [4, 8] */
	jac_[63] = (-(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP)-((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP+AKAPon_1*b_AKAP+AKAPoff_1*(1-b_AKAP); /* [5, 8] */
	jac_[74] = -((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 8] */
	jac_[85] = -RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [7, 8] */
	jac_[96] = (-RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF))-AKAPon_3*b_AKAP-AKAPon_1*b_AKAP-AKAPoff_3*(1-b_AKAP)-AKAPoff_1*(1-b_AKAP); /* [8, 8] */
	jac_[107] = 0; /* [9, 8] */
	jac_[118] = 0; /* [10, 8] */
/* column 10 (df/dy_9) */
	jac_[9] = 0; /* [0, 9] */
	jac_[20] = 0; /* [1, 9] */
	jac_[31] = -kf_Rii_C__RiiP_C; /* [2, 9] */
	jac_[42] = 0; /* [3, 9] */
	jac_[53] = (-C*kf_RiixC__Rii_C)+C*kf_C_AKAR4+kcat_AKARp-kb_RiixC__Rii_C+kb_C_AKAR4; /* [4, 9] */
	jac_[64] = (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP; /* [5, 9] */
	jac_[75] = -(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 9] */
	jac_[86] = 0; /* [7, 9] */
	jac_[97] = 0; /* [8, 9] */
	jac_[108] = (-C*kf_C_AKAR4)-kcat_AKARp-kb_C_AKAR4; /* [9, 9] */
	jac_[119] = kcat_AKARp; /* [10, 9] */
/* column 11 (df/dy_10) */
	jac_[10] = 0; /* [0, 10] */
	jac_[21] = 0; /* [1, 10] */
	jac_[32] = 0; /* [2, 10] */
	jac_[43] = 0; /* [3, 10] */
	jac_[54] = C*kf_C_AKAR4; /* [4, 10] */
	jac_[65] = 0; /* [5, 10] */
	jac_[76] = 0; /* [6, 10] */
	jac_[87] = 0; /* [7, 10] */
	jac_[98] = 0; /* [8, 10] */
	jac_[109] = -C*kf_C_AKAR4; /* [9, 10] */
	jac_[120] = 0; /* [10, 10] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAP79_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 11*33;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	double RiiP=y_[0];
	double RiiP_cAMP=y_[1];
	double RiiP_C=y_[2];
	double RiiP_C_cAMP=y_[3];
	double C=y_[4];
	double Rii_cAMP=y_[5];
	double Rii_C_cAMP=y_[6];
	double RiiP_CaN=y_[7];
	double RiiP_cAMP_CaN=y_[8];
	double AKAR4_C=y_[9];
	double AKAR4p=y_[10];
	double AKAR4=AKAR4_ConservedConst - (AKAR4_C+AKAR4p);
	double CaN=CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN);
	double Rii_C=Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	double cAMP=cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN);
	double Rii=Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C);
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
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
	jacp_[0] = 0; /* [0, 0] */
	jacp_[33] = 0; /* [1, 0] */
	jacp_[66] = (-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C; /* [2, 0] */
	jacp_[99] = 0; /* [3, 0] */
	jacp_[132] = 0; /* [4, 0] */
	jacp_[165] = 0; /* [5, 0] */
	jacp_[198] = 0; /* [6, 0] */
	jacp_[231] = 0; /* [7, 0] */
	jacp_[264] = 0; /* [8, 0] */
	jacp_[297] = 0; /* [9, 0] */
	jacp_[330] = 0; /* [10, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = 0; /* [0, 1] */
	jacp_[34] = 0; /* [1, 1] */
	jacp_[67] = KD_T*RiiP_C_cAMP-RiiP_C*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP); /* [2, 1] */
	jacp_[100] = RiiP_C*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)-KD_T*RiiP_C_cAMP; /* [3, 1] */
	jacp_[133] = 0; /* [4, 1] */
	jacp_[166] = 0; /* [5, 1] */
	jacp_[199] = 0; /* [6, 1] */
	jacp_[232] = 0; /* [7, 1] */
	jacp_[265] = 0; /* [8, 1] */
	jacp_[298] = 0; /* [9, 1] */
	jacp_[331] = 0; /* [10, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 0; /* [0, 2] */
	jacp_[35] = -C*RiiP_cAMP; /* [1, 2] */
	jacp_[68] = 0; /* [2, 2] */
	jacp_[101] = C*RiiP_cAMP; /* [3, 2] */
	jacp_[134] = -C*RiiP_cAMP; /* [4, 2] */
	jacp_[167] = 0; /* [5, 2] */
	jacp_[200] = 0; /* [6, 2] */
	jacp_[233] = 0; /* [7, 2] */
	jacp_[266] = 0; /* [8, 2] */
	jacp_[299] = 0; /* [9, 2] */
	jacp_[332] = 0; /* [10, 2] */
/* column 4 (df/dp_3) */
	jacp_[3] = 0; /* [0, 3] */
	jacp_[36] = RiiP_C_cAMP; /* [1, 3] */
	jacp_[69] = 0; /* [2, 3] */
	jacp_[102] = -RiiP_C_cAMP; /* [3, 3] */
	jacp_[135] = RiiP_C_cAMP; /* [4, 3] */
	jacp_[168] = 0; /* [5, 3] */
	jacp_[201] = 0; /* [6, 3] */
	jacp_[234] = 0; /* [7, 3] */
	jacp_[267] = 0; /* [8, 3] */
	jacp_[300] = 0; /* [9, 3] */
	jacp_[333] = 0; /* [10, 3] */
/* column 5 (df/dp_4) */
	jacp_[4] = RiiP_cAMP; /* [0, 4] */
	jacp_[37] = -RiiP_cAMP; /* [1, 4] */
	jacp_[70] = 0; /* [2, 4] */
	jacp_[103] = 0; /* [3, 4] */
	jacp_[136] = 0; /* [4, 4] */
	jacp_[169] = 0; /* [5, 4] */
	jacp_[202] = 0; /* [6, 4] */
	jacp_[235] = 0; /* [7, 4] */
	jacp_[268] = 0; /* [8, 4] */
	jacp_[301] = 0; /* [9, 4] */
	jacp_[334] = 0; /* [10, 4] */
/* column 6 (df/dp_5) */
	jacp_[5] = -RiiP*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP); /* [0, 5] */
	jacp_[38] = RiiP*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP); /* [1, 5] */
	jacp_[71] = 0; /* [2, 5] */
	jacp_[104] = 0; /* [3, 5] */
	jacp_[137] = 0; /* [4, 5] */
	jacp_[170] = 0; /* [5, 5] */
	jacp_[203] = 0; /* [6, 5] */
	jacp_[236] = 0; /* [7, 5] */
	jacp_[269] = 0; /* [8, 5] */
	jacp_[302] = 0; /* [9, 5] */
	jacp_[335] = 0; /* [10, 5] */
/* column 7 (df/dp_6) */
	jacp_[6] = -C*RiiP; /* [0, 6] */
	jacp_[39] = 0; /* [1, 6] */
	jacp_[72] = C*RiiP; /* [2, 6] */
	jacp_[105] = 0; /* [3, 6] */
	jacp_[138] = -C*RiiP; /* [4, 6] */
	jacp_[171] = 0; /* [5, 6] */
	jacp_[204] = 0; /* [6, 6] */
	jacp_[237] = 0; /* [7, 6] */
	jacp_[270] = 0; /* [8, 6] */
	jacp_[303] = 0; /* [9, 6] */
	jacp_[336] = 0; /* [10, 6] */
/* column 8 (df/dp_7) */
	jacp_[7] = RiiP_C; /* [0, 7] */
	jacp_[40] = 0; /* [1, 7] */
	jacp_[73] = -RiiP_C; /* [2, 7] */
	jacp_[106] = 0; /* [3, 7] */
	jacp_[139] = RiiP_C; /* [4, 7] */
	jacp_[172] = 0; /* [5, 7] */
	jacp_[205] = 0; /* [6, 7] */
	jacp_[238] = 0; /* [7, 7] */
	jacp_[271] = 0; /* [8, 7] */
	jacp_[304] = 0; /* [9, 7] */
	jacp_[337] = 0; /* [10, 7] */
/* column 9 (df/dp_8) */
	jacp_[8] = 0; /* [0, 8] */
	jacp_[41] = 0; /* [1, 8] */
	jacp_[74] = 0; /* [2, 8] */
	jacp_[107] = 0; /* [3, 8] */
	jacp_[140] = 0; /* [4, 8] */
	jacp_[173] = ((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP); /* [5, 8] */
	jacp_[206] = 0; /* [6, 8] */
	jacp_[239] = 0; /* [7, 8] */
	jacp_[272] = 0; /* [8, 8] */
	jacp_[305] = 0; /* [9, 8] */
	jacp_[338] = 0; /* [10, 8] */
/* column 10 (df/dp_9) */
	jacp_[9] = 0; /* [0, 9] */
	jacp_[42] = 0; /* [1, 9] */
	jacp_[75] = 0; /* [2, 9] */
	jacp_[108] = 0; /* [3, 9] */
	jacp_[141] = 0; /* [4, 9] */
	jacp_[174] = -Rii_cAMP; /* [5, 9] */
	jacp_[207] = 0; /* [6, 9] */
	jacp_[240] = 0; /* [7, 9] */
	jacp_[273] = 0; /* [8, 9] */
	jacp_[306] = 0; /* [9, 9] */
	jacp_[339] = 0; /* [10, 9] */
/* column 11 (df/dp_10) */
	jacp_[10] = 0; /* [0, 10] */
	jacp_[43] = 0; /* [1, 10] */
	jacp_[76] = 0; /* [2, 10] */
	jacp_[109] = 0; /* [3, 10] */
	jacp_[142] = 0; /* [4, 10] */
	jacp_[175] = 0; /* [5, 10] */
	jacp_[208] = ((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*(cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP); /* [6, 10] */
	jacp_[241] = 0; /* [7, 10] */
	jacp_[274] = 0; /* [8, 10] */
	jacp_[307] = 0; /* [9, 10] */
	jacp_[340] = 0; /* [10, 10] */
/* column 12 (df/dp_11) */
	jacp_[11] = 0; /* [0, 11] */
	jacp_[44] = 0; /* [1, 11] */
	jacp_[77] = 0; /* [2, 11] */
	jacp_[110] = 0; /* [3, 11] */
	jacp_[143] = 0; /* [4, 11] */
	jacp_[176] = 0; /* [5, 11] */
	jacp_[209] = -Rii_C_cAMP; /* [6, 11] */
	jacp_[242] = 0; /* [7, 11] */
	jacp_[275] = 0; /* [8, 11] */
	jacp_[308] = 0; /* [9, 11] */
	jacp_[341] = 0; /* [10, 11] */
/* column 13 (df/dp_12) */
	jacp_[12] = 0; /* [0, 12] */
	jacp_[45] = 0; /* [1, 12] */
	jacp_[78] = 0; /* [2, 12] */
	jacp_[111] = 0; /* [3, 12] */
	jacp_[144] = -C*((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C); /* [4, 12] */
	jacp_[177] = 0; /* [5, 12] */
	jacp_[210] = 0; /* [6, 12] */
	jacp_[243] = 0; /* [7, 12] */
	jacp_[276] = 0; /* [8, 12] */
	jacp_[309] = 0; /* [9, 12] */
	jacp_[342] = 0; /* [10, 12] */
/* column 14 (df/dp_13) */
	jacp_[13] = 0; /* [0, 13] */
	jacp_[46] = 0; /* [1, 13] */
	jacp_[79] = 0; /* [2, 13] */
	jacp_[112] = 0; /* [3, 13] */
	jacp_[145] = -C*Rii_cAMP; /* [4, 13] */
	jacp_[178] = -C*Rii_cAMP; /* [5, 13] */
	jacp_[211] = C*Rii_cAMP; /* [6, 13] */
	jacp_[244] = 0; /* [7, 13] */
	jacp_[277] = 0; /* [8, 13] */
	jacp_[310] = 0; /* [9, 13] */
	jacp_[343] = 0; /* [10, 13] */
/* column 15 (df/dp_14) */
	jacp_[14] = 0; /* [0, 14] */
	jacp_[47] = 0; /* [1, 14] */
	jacp_[80] = 0; /* [2, 14] */
	jacp_[113] = 0; /* [3, 14] */
	jacp_[146] = Rii_C_cAMP; /* [4, 14] */
	jacp_[179] = Rii_C_cAMP; /* [5, 14] */
	jacp_[212] = -Rii_C_cAMP; /* [6, 14] */
	jacp_[245] = 0; /* [7, 14] */
	jacp_[278] = 0; /* [8, 14] */
	jacp_[311] = 0; /* [9, 14] */
	jacp_[344] = 0; /* [10, 14] */
/* column 16 (df/dp_15) */
	jacp_[15] = 0; /* [0, 15] */
	jacp_[48] = 0; /* [1, 15] */
	jacp_[81] = 0; /* [2, 15] */
	jacp_[114] = Rii_C_cAMP; /* [3, 15] */
	jacp_[147] = 0; /* [4, 15] */
	jacp_[180] = 0; /* [5, 15] */
	jacp_[213] = -Rii_C_cAMP; /* [6, 15] */
	jacp_[246] = 0; /* [7, 15] */
	jacp_[279] = 0; /* [8, 15] */
	jacp_[312] = 0; /* [9, 15] */
	jacp_[345] = 0; /* [10, 15] */
/* column 17 (df/dp_16) */
	jacp_[16] = 0; /* [0, 16] */
	jacp_[49] = 0; /* [1, 16] */
	jacp_[82] = 0; /* [2, 16] */
	jacp_[115] = 0; /* [3, 16] */
	jacp_[148] = (-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C; /* [4, 16] */
	jacp_[181] = 0; /* [5, 16] */
	jacp_[214] = 0; /* [6, 16] */
	jacp_[247] = 0; /* [7, 16] */
	jacp_[280] = 0; /* [8, 16] */
	jacp_[313] = 0; /* [9, 16] */
	jacp_[346] = 0; /* [10, 16] */
/* column 18 (df/dp_17) */
	jacp_[17] = -RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [0, 17] */
	jacp_[50] = -RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [1, 17] */
	jacp_[83] = 0; /* [2, 17] */
	jacp_[116] = 0; /* [3, 17] */
	jacp_[149] = 0; /* [4, 17] */
	jacp_[182] = RiiP_cAMP_CaN*(1-b_AKAP); /* [5, 17] */
	jacp_[215] = 0; /* [6, 17] */
	jacp_[248] = RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_CaN*(1-b_AKAP); /* [7, 17] */
	jacp_[281] = RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP); /* [8, 17] */
	jacp_[314] = 0; /* [9, 17] */
	jacp_[347] = 0; /* [10, 17] */
/* column 19 (df/dp_18) */
	jacp_[18] = RiiP_CaN*(1-b_AKAP)-RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [0, 18] */
	jacp_[51] = RiiP_cAMP_CaN*(1-b_AKAP)-RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF); /* [1, 18] */
	jacp_[84] = 0; /* [2, 18] */
	jacp_[117] = 0; /* [3, 18] */
	jacp_[150] = 0; /* [4, 18] */
	jacp_[183] = 0; /* [5, 18] */
	jacp_[216] = 0; /* [6, 18] */
	jacp_[249] = RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_CaN*(1-b_AKAP); /* [7, 18] */
	jacp_[282] = RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(((1-b_AKAP)*b_AKAP)/kmON+gsl_pow_2(1-b_AKAP)/kmOFF)-RiiP_cAMP_CaN*(1-b_AKAP); /* [8, 18] */
	jacp_[315] = 0; /* [9, 18] */
	jacp_[348] = 0; /* [10, 18] */
/* column 20 (df/dp_19) */
	jacp_[19] = -RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [0, 19] */
	jacp_[52] = -RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [1, 19] */
	jacp_[85] = 0; /* [2, 19] */
	jacp_[118] = 0; /* [3, 19] */
	jacp_[151] = 0; /* [4, 19] */
	jacp_[184] = RiiP_cAMP_CaN*b_AKAP; /* [5, 19] */
	jacp_[217] = 0; /* [6, 19] */
	jacp_[250] = RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP; /* [7, 19] */
	jacp_[283] = RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP; /* [8, 19] */
	jacp_[316] = 0; /* [9, 19] */
	jacp_[349] = 0; /* [10, 19] */
/* column 21 (df/dp_20) */
	jacp_[20] = RiiP_CaN*b_AKAP-RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [0, 20] */
	jacp_[53] = RiiP_cAMP_CaN*b_AKAP-RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF); /* [1, 20] */
	jacp_[86] = 0; /* [2, 20] */
	jacp_[119] = 0; /* [3, 20] */
	jacp_[152] = 0; /* [4, 20] */
	jacp_[185] = 0; /* [5, 20] */
	jacp_[218] = 0; /* [6, 20] */
	jacp_[251] = RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_CaN*b_AKAP; /* [7, 20] */
	jacp_[284] = RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(gsl_pow_2(b_AKAP)/kmON+((1-b_AKAP)*b_AKAP)/kmOFF)-RiiP_cAMP_CaN*b_AKAP; /* [8, 20] */
	jacp_[317] = 0; /* [9, 20] */
	jacp_[350] = 0; /* [10, 20] */
/* column 22 (df/dp_21) */
	jacp_[21] = 0; /* [0, 21] */
	jacp_[54] = 0; /* [1, 21] */
	jacp_[87] = 0; /* [2, 21] */
	jacp_[120] = 0; /* [3, 21] */
	jacp_[153] = -((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*C; /* [4, 21] */
	jacp_[186] = 0; /* [5, 21] */
	jacp_[219] = 0; /* [6, 21] */
	jacp_[252] = 0; /* [7, 21] */
	jacp_[285] = 0; /* [8, 21] */
	jacp_[318] = ((-AKAR4p)+AKAR4_ConservedConst-AKAR4_C)*C; /* [9, 21] */
	jacp_[351] = 0; /* [10, 21] */
/* column 23 (df/dp_22) */
	jacp_[22] = 0; /* [0, 22] */
	jacp_[55] = 0; /* [1, 22] */
	jacp_[88] = 0; /* [2, 22] */
	jacp_[121] = 0; /* [3, 22] */
	jacp_[154] = AKAR4_C; /* [4, 22] */
	jacp_[187] = 0; /* [5, 22] */
	jacp_[220] = 0; /* [6, 22] */
	jacp_[253] = 0; /* [7, 22] */
	jacp_[286] = 0; /* [8, 22] */
	jacp_[319] = -AKAR4_C; /* [9, 22] */
	jacp_[352] = 0; /* [10, 22] */
/* column 24 (df/dp_23) */
	jacp_[23] = 0; /* [0, 23] */
	jacp_[56] = 0; /* [1, 23] */
	jacp_[89] = 0; /* [2, 23] */
	jacp_[122] = 0; /* [3, 23] */
	jacp_[155] = AKAR4_C; /* [4, 23] */
	jacp_[188] = 0; /* [5, 23] */
	jacp_[221] = 0; /* [6, 23] */
	jacp_[254] = 0; /* [7, 23] */
	jacp_[287] = 0; /* [8, 23] */
	jacp_[320] = -AKAR4_C; /* [9, 23] */
	jacp_[353] = AKAR4_C; /* [10, 23] */
/* column 25 (df/dp_24) */
	jacp_[24] = (RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [0, 24] */
	jacp_[57] = (RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [1, 24] */
	jacp_[90] = 0; /* [2, 24] */
	jacp_[123] = 0; /* [3, 24] */
	jacp_[156] = 0; /* [4, 24] */
	jacp_[189] = 0; /* [5, 24] */
	jacp_[222] = 0; /* [6, 24] */
	jacp_[255] = -(RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [7, 24] */
	jacp_[288] = -(RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*(1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmOFF); /* [8, 24] */
	jacp_[321] = 0; /* [9, 24] */
	jacp_[354] = 0; /* [10, 24] */
/* column 26 (df/dp_25) */
	jacp_[25] = (RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [0, 25] */
	jacp_[58] = (RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [1, 25] */
	jacp_[91] = 0; /* [2, 25] */
	jacp_[124] = 0; /* [3, 25] */
	jacp_[157] = 0; /* [4, 25] */
	jacp_[190] = 0; /* [5, 25] */
	jacp_[223] = 0; /* [6, 25] */
	jacp_[256] = -(RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [7, 25] */
	jacp_[289] = -(RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/gsl_pow_2(kmON); /* [8, 25] */
	jacp_[322] = 0; /* [9, 25] */
	jacp_[355] = 0; /* [10, 25] */
/* column 27 (df/dp_26) */
	jacp_[26] = 0; /* [0, 26] */
	jacp_[59] = 0; /* [1, 26] */
	jacp_[92] = RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [2, 26] */
	jacp_[125] = -RiiP_C_cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 26] */
	jacp_[158] = 0; /* [4, 26] */
	jacp_[191] = 0; /* [5, 26] */
	jacp_[224] = 0; /* [6, 26] */
	jacp_[257] = 0; /* [7, 26] */
	jacp_[290] = 0; /* [8, 26] */
	jacp_[323] = 0; /* [9, 26] */
	jacp_[356] = 0; /* [10, 26] */
/* column 28 (df/dp_27) */
	jacp_[27] = (AKAPon_3-AKAPoff_3)*RiiP_CaN-RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF); /* [0, 27] */
	jacp_[60] = (AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF); /* [1, 27] */
	jacp_[93] = 0; /* [2, 27] */
	jacp_[126] = 0; /* [3, 27] */
	jacp_[159] = 0; /* [4, 27] */
	jacp_[192] = (AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN; /* [5, 27] */
	jacp_[225] = 0; /* [6, 27] */
	jacp_[258] = RiiP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_CaN; /* [7, 27] */
	jacp_[291] = RiiP_cAMP*((-RiiP_cAMP_CaN)-RiiP_CaN+CaN_ConservedConst)*((AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmON+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*b_AKAP)/kmON-(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP))/kmOFF+((AKAPon_3+AKAPon_1-AKAPoff_3-AKAPoff_1)*(1-b_AKAP))/kmOFF)-(AKAPon_3-AKAPoff_3)*RiiP_cAMP_CaN-(AKAPon_1-AKAPoff_1)*RiiP_cAMP_CaN; /* [8, 27] */
	jacp_[324] = 0; /* [9, 27] */
	jacp_[357] = 0; /* [10, 27] */
/* column 29 (df/dp_28) */
	jacp_[28] = 0; /* [0, 28] */
	jacp_[61] = 0; /* [1, 28] */
	jacp_[94] = 0; /* [2, 28] */
	jacp_[127] = 0; /* [3, 28] */
	jacp_[160] = -C*kf_C_AKAR4; /* [4, 28] */
	jacp_[193] = 0; /* [5, 28] */
	jacp_[226] = 0; /* [6, 28] */
	jacp_[259] = 0; /* [7, 28] */
	jacp_[292] = 0; /* [8, 28] */
	jacp_[325] = C*kf_C_AKAR4; /* [9, 28] */
	jacp_[358] = 0; /* [10, 28] */
/* column 30 (df/dp_29) */
	jacp_[29] = -RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [0, 29] */
	jacp_[62] = -RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [1, 29] */
	jacp_[95] = 0; /* [2, 29] */
	jacp_[128] = 0; /* [3, 29] */
	jacp_[161] = 0; /* [4, 29] */
	jacp_[194] = 0; /* [5, 29] */
	jacp_[227] = 0; /* [6, 29] */
	jacp_[260] = RiiP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [7, 29] */
	jacp_[293] = RiiP_cAMP*((b_AKAP*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmON+((1-b_AKAP)*(AKAPon_3*b_AKAP+AKAPon_1*b_AKAP+AKAPoff_3*(1-b_AKAP)+AKAPoff_1*(1-b_AKAP)))/kmOFF); /* [8, 29] */
	jacp_[326] = 0; /* [9, 29] */
	jacp_[359] = 0; /* [10, 29] */
/* column 31 (df/dp_30) */
	jacp_[30] = 0; /* [0, 30] */
	jacp_[63] = 0; /* [1, 30] */
	jacp_[96] = kf_Rii_C__RiiP_C; /* [2, 30] */
	jacp_[129] = 0; /* [3, 30] */
	jacp_[162] = kb_RiixC__Rii_C; /* [4, 30] */
	jacp_[195] = 0; /* [5, 30] */
	jacp_[228] = (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 30] */
	jacp_[261] = 0; /* [7, 30] */
	jacp_[294] = 0; /* [8, 30] */
	jacp_[327] = 0; /* [9, 30] */
	jacp_[360] = 0; /* [10, 30] */
/* column 32 (df/dp_31) */
	jacp_[31] = -RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [0, 31] */
	jacp_[64] = RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 31] */
	jacp_[97] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [2, 31] */
	jacp_[130] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [3, 31] */
	jacp_[163] = 0; /* [4, 31] */
	jacp_[196] = ((-Rii_cAMP)+Rii_ConservedConst-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_CaN-RiiP+C+AKAR4_C)*kf_cAMPxRii__Rii_cAMP; /* [5, 31] */
	jacp_[229] = ((-Rii_C_cAMP)+Rii_C_ConservedConst-RiiP_C_cAMP-RiiP_C-C-AKAR4_C)*kf_Rii_CxcAMP__Rii_C_cAMP; /* [6, 31] */
	jacp_[262] = 0; /* [7, 31] */
	jacp_[295] = 0; /* [8, 31] */
	jacp_[328] = 0; /* [9, 31] */
	jacp_[361] = 0; /* [10, 31] */
/* column 33 (df/dp_32) */
	jacp_[32] = 0; /* [0, 32] */
	jacp_[65] = 0; /* [1, 32] */
	jacp_[98] = 0; /* [2, 32] */
	jacp_[131] = 0; /* [3, 32] */
	jacp_[164] = -C*kf_RiixC__Rii_C; /* [4, 32] */
	jacp_[197] = (cAMP_ConservedConst-Rii_cAMP-Rii_C_cAMP-RiiP_cAMP_CaN-RiiP_cAMP-RiiP_C_cAMP)*kf_cAMPxRii__Rii_cAMP; /* [5, 32] */
	jacp_[230] = 0; /* [6, 32] */
	jacp_[263] = 0; /* [7, 32] */
	jacp_[296] = 0; /* [8, 32] */
	jacp_[329] = 0; /* [9, 32] */
	jacp_[362] = 0; /* [10, 32] */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int AKAP79_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 1;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	double RiiP=y_[0];
	double RiiP_cAMP=y_[1];
	double RiiP_C=y_[2];
	double RiiP_C_cAMP=y_[3];
	double C=y_[4];
	double Rii_cAMP=y_[5];
	double Rii_C_cAMP=y_[6];
	double RiiP_CaN=y_[7];
	double RiiP_cAMP_CaN=y_[8];
	double AKAR4_C=y_[9];
	double AKAR4p=y_[10];
	double AKAR4=AKAR4_ConservedConst - (AKAR4_C+AKAR4p);
	double CaN=CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN);
	double Rii_C=Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	double cAMP=cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN);
	double Rii=Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C);
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
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
int AKAP79_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 11;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	double RiiP=y_[0];
	double RiiP_cAMP=y_[1];
	double RiiP_C=y_[2];
	double RiiP_C_cAMP=y_[3];
	double C=y_[4];
	double Rii_cAMP=y_[5];
	double Rii_C_cAMP=y_[6];
	double RiiP_CaN=y_[7];
	double RiiP_cAMP_CaN=y_[8];
	double AKAR4_C=y_[9];
	double AKAR4p=y_[10];
	double AKAR4=AKAR4_ConservedConst - (AKAR4_C+AKAR4p);
	double CaN=CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN);
	double Rii_C=Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	double cAMP=cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN);
	double Rii=Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C);
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
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
	funcJac_[10] = 358.35; /* [0, 10] */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int AKAP79_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 33;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	double RiiP=y_[0];
	double RiiP_cAMP=y_[1];
	double RiiP_C=y_[2];
	double RiiP_C_cAMP=y_[3];
	double C=y_[4];
	double Rii_cAMP=y_[5];
	double Rii_C_cAMP=y_[6];
	double RiiP_CaN=y_[7];
	double RiiP_cAMP_CaN=y_[8];
	double AKAR4_C=y_[9];
	double AKAR4p=y_[10];
	double AKAR4=AKAR4_ConservedConst - (AKAR4_C+AKAR4p);
	double CaN=CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN);
	double Rii_C=Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C);
	double cAMP=cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN);
	double Rii=Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C);
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN=b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN=b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN=b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP=kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
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
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAP79_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 33;
	p_[0] = 33;
	p_[1] = 0.496;
	p_[2] = 0.00545;
	p_[3] = 0.0156;
	p_[4] = 0.0016;
	p_[5] = 0.015;
	p_[6] = 0.038;
	p_[7] = 0.0026;
	p_[8] = 0.015;
	p_[9] = 0.0016;
	p_[10] = 0.496;
	p_[11] = 1.413;
	p_[12] = 2.1;
	p_[13] = 0.2984;
	p_[14] = 0.018;
	p_[15] = 33;
	p_[16] = 0.0003;
	p_[17] = 2.6;
	p_[18] = 20;
	p_[19] = 0.45;
	p_[20] = 2;
	p_[21] = 0.018;
	p_[22] = 0.106;
	p_[23] = 10.2;
	p_[24] = 100;
	p_[25] = 1;
	p_[26] = 0.7;
	p_[27] = 0;
	p_[28] = 0.200000;
	p_[29] = 1.500000;
	p_[30] = 0.630000;
	p_[31] = 0.000000;
	p_[32] = 6.300000;
	return GSL_SUCCESS;
}
/* ode initial values */
int AKAP79_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 11;
	double kf_Rii_C__RiiP_C=p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP=p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP=p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP=p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP=p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP=p_[5];
	double kf_RiiPxC__RiiP_C=p_[6];
	double kb_RiiPxC__RiiP_C=p_[7];
	double kf_cAMPxRii__Rii_cAMP=p_[8];
	double kb_cAMPxRii__Rii_cAMP=p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP=p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP=p_[11];
	double kf_RiixC__Rii_C=p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP=p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP=p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP=p_[15];
	double kb_RiixC__Rii_C=p_[16];
	double AKAPoff_1=p_[17];
	double AKAPoff_3=p_[18];
	double AKAPon_1=p_[19];
	double AKAPon_3=p_[20];
	double kf_C_AKAR4=p_[21];
	double kb_C_AKAR4=p_[22];
	double kcat_AKARp=p_[23];
	double kmOFF=p_[24];
	double kmON=p_[25];
	double KD_T=p_[26];
	double b_AKAP=p_[27];
	double AKAR4_ConservedConst=p_[28];
	double CaN_ConservedConst=p_[29];
	double Rii_C_ConservedConst=p_[30];
	double cAMP_ConservedConst=p_[31];
	double Rii_ConservedConst=p_[32];
	/* the initial value of y may depend on the parameters. */
	y_[0] = 0;
	y_[1] = 0;
	y_[2] = 0;
	y_[3] = 0;
	y_[4] = 0;
	y_[5] = 0;
	y_[6] = 0;
	y_[7] = 0;
	y_[8] = 0;
	y_[9] = 0;
	y_[10] = 0;
	return GSL_SUCCESS;
}
