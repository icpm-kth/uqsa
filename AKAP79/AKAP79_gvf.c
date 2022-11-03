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
int AKAP79_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 16;
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
	double reaction_44_=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33_=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_4_8=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_3_7=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	f_[0] = -reaction_78-reaction_58+reaction_4_8;
	f_[1] = -reaction_12-reaction_43-reaction_78-reaction_56;
	f_[2] = -reaction_14-reaction_43-reaction_44_;
	f_[3] = -reaction_51-reaction_56+reaction_58;
	f_[4] = +reaction_43-reaction_23-reaction_33_;
	f_[5] = +reaction_51+reaction_14-reaction_12;
	f_[6] = +reaction_12+reaction_23+reaction_62;
	f_[7] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2;
	f_[8] = +reaction_78-reaction_76+reaction_3_7;
	f_[9] = +reaction_56+reaction_76-reaction_62;
	f_[10] = -reaction_44_-reaction_33_+reaction_4_8+reaction_3_7;
	f_[11] = +reaction_44_-reaction_4_8;
	f_[12] = +reaction_33_-reaction_3_7;
	f_[13] = -reaction_1;
	f_[14] = +reaction_1-reaction_2;
	f_[15] = +reaction_2;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int AKAP79_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 16*16;
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
	double reaction_44_=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33_=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_4_8=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_3_7=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dy_0) */
	jac_[0] = ((-1*(kf_cAMPxRii__Rii_cAMP*cAMP))-(kf_RiixC__Rii_C*C)); /* [0, 0] */
	jac_[16] = (kf_RiixC__Rii_C*C); /* [1, 0] */
	jac_[32] = (kf_cAMPxRii__Rii_cAMP*cAMP); /* [2, 0] */
	jac_[48] = 0; /* [3, 0] */
	jac_[64] = 0; /* [4, 0] */
	jac_[80] = 0; /* [5, 0] */
	jac_[96] = 0; /* [6, 0] */
	jac_[112] = 0; /* [7, 0] */
	jac_[128] = 0; /* [8, 0] */
	jac_[144] = 0; /* [9, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = (-1*(kf_cAMPxRii__Rii_cAMP*Rii)); /* [0, 1] */
	jac_[17] = (0-(kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C)); /* [1, 1] */
	jac_[33] = (kf_cAMPxRii__Rii_cAMP*Rii); /* [2, 1] */
	jac_[49] = (kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C); /* [3, 1] */
	jac_[65] = 0; /* [4, 1] */
	jac_[81] = 0; /* [5, 1] */
	jac_[97] = 0; /* [6, 1] */
	jac_[113] = 0; /* [7, 1] */
	jac_[129] = 0; /* [8, 1] */
	jac_[145] = 0; /* [9, 1] */
/* column 3 (df/dy_2) */
	jac_[2] = 0; /* [0, 2] */
	jac_[18] = 0; /* [1, 2] */
	jac_[34] = 0; /* [2, 2] */
	jac_[50] = 0; /* [3, 2] */
	jac_[66] = (CaN*(-1*((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))))); /* [4, 2] */
	jac_[82] = (((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF)))*CaN); /* [5, 2] */
	jac_[98] = 0; /* [6, 2] */
	jac_[114] = 0; /* [7, 2] */
	jac_[130] = 0; /* [8, 2] */
	jac_[146] = 0; /* [9, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = (0-(0-kb_RiixC__Rii_C)); /* [0, 3] */
	jac_[19] = (((-1*kf_Rii_C__RiiP_C)-(kf_Rii_CxcAMP__Rii_C_cAMP*cAMP))+(0-kb_RiixC__Rii_C)); /* [1, 3] */
	jac_[35] = 0; /* [2, 3] */
	jac_[51] = (kf_Rii_CxcAMP__Rii_C_cAMP*cAMP); /* [3, 3] */
	jac_[67] = 0; /* [4, 3] */
	jac_[83] = 0; /* [5, 3] */
	jac_[99] = 0; /* [6, 3] */
	jac_[115] = 0; /* [7, 3] */
	jac_[131] = 0; /* [8, 3] */
	jac_[147] = 0; /* [9, 3] */
/* column 5 (df/dy_4) */
	jac_[4] = 0; /* [0, 4] */
	jac_[20] = 0; /* [1, 4] */
	jac_[36] = 0; /* [2, 4] */
	jac_[52] = 0; /* [3, 4] */
	jac_[68] = (CaN*(0-((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))))); /* [4, 4] */
	jac_[84] = 0; /* [5, 4] */
	jac_[100] = (CaN*((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF)))); /* [6, 4] */
	jac_[116] = 0; /* [7, 4] */
	jac_[132] = 0; /* [8, 4] */
	jac_[148] = 0; /* [9, 4] */
/* column 6 (df/dy_5) */
	jac_[5] = 0; /* [0, 5] */
	jac_[21] = 0; /* [1, 5] */
	jac_[37] = 0; /* [2, 5] */
	jac_[53] = 0; /* [3, 5] */
	jac_[69] = 0; /* [4, 5] */
	jac_[85] = 0; /* [5, 5] */
	jac_[101] = 0; /* [6, 5] */
	jac_[117] = 0; /* [7, 5] */
	jac_[133] = 0; /* [8, 5] */
	jac_[149] = 0; /* [9, 5] */
/* column 7 (df/dy_6) */
	jac_[6] = 0; /* [0, 6] */
	jac_[22] = 0; /* [1, 6] */
	jac_[38] = 0; /* [2, 6] */
	jac_[54] = 0; /* [3, 6] */
	jac_[70] = 0; /* [4, 6] */
	jac_[86] = 0; /* [5, 6] */
	jac_[102] = 0; /* [6, 6] */
	jac_[118] = 0; /* [7, 6] */
	jac_[134] = 0; /* [8, 6] */
	jac_[150] = 0; /* [9, 6] */
/* column 8 (df/dy_7) */
	jac_[7] = (0-(kf_RiixC__Rii_C*Rii)); /* [0, 7] */
	jac_[23] = (kf_RiixC__Rii_C*Rii); /* [1, 7] */
	jac_[39] = (0-(kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP)); /* [2, 7] */
	jac_[55] = (kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP); /* [3, 7] */
	jac_[71] = 0; /* [4, 7] */
	jac_[87] = 0; /* [5, 7] */
	jac_[103] = 0; /* [6, 7] */
	jac_[119] = (-1*(kf_C_AKAR4*AKAR4)); /* [7, 7] */
	jac_[135] = (kf_C_AKAR4*AKAR4); /* [8, 7] */
	jac_[151] = 0; /* [9, 7] */
/* column 9 (df/dy_8) */
	jac_[8] = (-1*(0-kb_cAMPxRii__Rii_cAMP)); /* [0, 8] */
	jac_[24] = 0; /* [1, 8] */
	jac_[40] = ((0-kb_cAMPxRii__Rii_cAMP)-(kf_Rii_cAMPxC__Rii_C_cAMP*C)); /* [2, 8] */
	jac_[56] = (kf_Rii_cAMPxC__Rii_C_cAMP*C); /* [3, 8] */
	jac_[72] = 0; /* [4, 8] */
	jac_[88] = 0; /* [5, 8] */
	jac_[104] = 0; /* [6, 8] */
	jac_[120] = 0; /* [7, 8] */
	jac_[136] = 0; /* [8, 8] */
	jac_[152] = 0; /* [9, 8] */
/* column 10 (df/dy_9) */
	jac_[9] = 0; /* [0, 9] */
	jac_[25] = (0-(0-kb_Rii_CxcAMP__Rii_C_cAMP)); /* [1, 9] */
	jac_[41] = (0-(0-kb_Rii_cAMPxC__Rii_C_cAMP)); /* [2, 9] */
	jac_[57] = (((0-kb_Rii_CxcAMP__Rii_C_cAMP)+(0-kb_Rii_cAMPxC__Rii_C_cAMP))-kf_Rii_C_cAMP__RiiP_C_cAMP); /* [3, 9] */
	jac_[73] = 0; /* [4, 9] */
	jac_[89] = 0; /* [5, 9] */
	jac_[105] = 0; /* [6, 9] */
	jac_[121] = 0; /* [7, 9] */
	jac_[137] = 0; /* [8, 9] */
	jac_[153] = 0; /* [9, 9] */
/* column 11 (df/dy_10) */
	jac_[10] = 0; /* [0, 10] */
	jac_[26] = 0; /* [1, 10] */
	jac_[42] = 0; /* [2, 10] */
	jac_[58] = 0; /* [3, 10] */
	jac_[74] = (((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF)))*((-1*RiiP)-RiiP_cAMP)); /* [4, 10] */
	jac_[90] = (RiiP*((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF)))); /* [5, 10] */
	jac_[106] = (((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))*((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF)))*RiiP_cAMP); /* [6, 10] */
	jac_[122] = 0; /* [7, 10] */
	jac_[138] = 0; /* [8, 10] */
	jac_[154] = 0; /* [9, 10] */
/* column 12 (df/dy_11) */
	jac_[11] = ((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)); /* [0, 11] */
	jac_[27] = 0; /* [1, 11] */
	jac_[43] = 0; /* [2, 11] */
	jac_[59] = 0; /* [3, 11] */
	jac_[75] = ((-1*(0-((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))); /* [4, 11] */
	jac_[91] = ((0-((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3)))-((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))); /* [5, 11] */
	jac_[107] = 0; /* [6, 11] */
	jac_[123] = 0; /* [7, 11] */
	jac_[139] = 0; /* [8, 11] */
	jac_[155] = 0; /* [9, 11] */
/* column 13 (df/dy_12) */
	jac_[12] = 0; /* [0, 12] */
	jac_[28] = 0; /* [1, 12] */
	jac_[44] = ((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)); /* [2, 12] */
	jac_[60] = 0; /* [3, 12] */
	jac_[76] = ((0-(0-((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))); /* [4, 12] */
	jac_[92] = 0; /* [5, 12] */
	jac_[108] = ((0-((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3)))-((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))); /* [6, 12] */
	jac_[124] = 0; /* [7, 12] */
	jac_[140] = 0; /* [8, 12] */
	jac_[156] = 0; /* [9, 12] */
/* column 14 (df/dy_13) */
	jac_[13] = 0; /* [0, 13] */
	jac_[29] = 0; /* [1, 13] */
	jac_[45] = 0; /* [2, 13] */
	jac_[61] = 0; /* [3, 13] */
	jac_[77] = 0; /* [4, 13] */
	jac_[93] = 0; /* [5, 13] */
	jac_[109] = 0; /* [6, 13] */
	jac_[125] = (-1*(kf_C_AKAR4*C)); /* [7, 13] */
	jac_[141] = (kf_C_AKAR4*C); /* [8, 13] */
	jac_[157] = 0; /* [9, 13] */
/* column 15 (df/dy_14) */
	jac_[14] = 0; /* [0, 14] */
	jac_[30] = 0; /* [1, 14] */
	jac_[46] = 0; /* [2, 14] */
	jac_[62] = 0; /* [3, 14] */
	jac_[78] = 0; /* [4, 14] */
	jac_[94] = 0; /* [5, 14] */
	jac_[110] = 0; /* [6, 14] */
	jac_[126] = (-1*(0-kb_C_AKAR4)); /* [7, 14] */
	jac_[142] = ((0-kb_C_AKAR4)-kcat_AKARp); /* [8, 14] */
	jac_[158] = kcat_AKARp; /* [9, 14] */
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
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int AKAP79_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 16*28;
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
	double reaction_44_=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33_=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_4_8=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_3_7=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
/* column 1 (df/dp_0) */
	jacp_[0] = 0; /* [0, 0] */
	jacp_[28] = (-1*Rii_C); /* [1, 0] */
	jacp_[56] = 0; /* [2, 0] */
	jacp_[84] = 0; /* [3, 0] */
	jacp_[112] = 0; /* [4, 0] */
	jacp_[140] = 0; /* [5, 0] */
	jacp_[168] = 0; /* [6, 0] */
	jacp_[196] = 0; /* [7, 0] */
	jacp_[224] = 0; /* [8, 0] */
	jacp_[252] = 0; /* [9, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = 0; /* [0, 1] */
	jacp_[29] = 0; /* [1, 1] */
	jacp_[57] = 0; /* [2, 1] */
	jacp_[85] = 0; /* [3, 1] */
	jacp_[113] = 0; /* [4, 1] */
	jacp_[141] = 0; /* [5, 1] */
	jacp_[169] = 0; /* [6, 1] */
	jacp_[197] = 0; /* [7, 1] */
	jacp_[225] = 0; /* [8, 1] */
	jacp_[253] = 0; /* [9, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 0; /* [0, 2] */
	jacp_[30] = 0; /* [1, 2] */
	jacp_[58] = 0; /* [2, 2] */
	jacp_[86] = 0; /* [3, 2] */
	jacp_[114] = 0; /* [4, 2] */
	jacp_[142] = 0; /* [5, 2] */
	jacp_[170] = 0; /* [6, 2] */
	jacp_[198] = 0; /* [7, 2] */
	jacp_[226] = 0; /* [8, 2] */
	jacp_[254] = 0; /* [9, 2] */
/* column 4 (df/dp_3) */
	jacp_[3] = 0; /* [0, 3] */
	jacp_[31] = 0; /* [1, 3] */
	jacp_[59] = 0; /* [2, 3] */
	jacp_[87] = 0; /* [3, 3] */
	jacp_[115] = 0; /* [4, 3] */
	jacp_[143] = 0; /* [5, 3] */
	jacp_[171] = 0; /* [6, 3] */
	jacp_[199] = 0; /* [7, 3] */
	jacp_[227] = 0; /* [8, 3] */
	jacp_[255] = 0; /* [9, 3] */
/* column 5 (df/dp_4) */
	jacp_[4] = 0; /* [0, 4] */
	jacp_[32] = 0; /* [1, 4] */
	jacp_[60] = 0; /* [2, 4] */
	jacp_[88] = 0; /* [3, 4] */
	jacp_[116] = 0; /* [4, 4] */
	jacp_[144] = 0; /* [5, 4] */
	jacp_[172] = 0; /* [6, 4] */
	jacp_[200] = 0; /* [7, 4] */
	jacp_[228] = 0; /* [8, 4] */
	jacp_[256] = 0; /* [9, 4] */
/* column 6 (df/dp_5) */
	jacp_[5] = 0; /* [0, 5] */
	jacp_[33] = 0; /* [1, 5] */
	jacp_[61] = 0; /* [2, 5] */
	jacp_[89] = 0; /* [3, 5] */
	jacp_[117] = 0; /* [4, 5] */
	jacp_[145] = 0; /* [5, 5] */
	jacp_[173] = 0; /* [6, 5] */
	jacp_[201] = 0; /* [7, 5] */
	jacp_[229] = 0; /* [8, 5] */
	jacp_[257] = 0; /* [9, 5] */
/* column 7 (df/dp_6) */
	jacp_[6] = 0; /* [0, 6] */
	jacp_[34] = 0; /* [1, 6] */
	jacp_[62] = 0; /* [2, 6] */
	jacp_[90] = 0; /* [3, 6] */
	jacp_[118] = 0; /* [4, 6] */
	jacp_[146] = 0; /* [5, 6] */
	jacp_[174] = 0; /* [6, 6] */
	jacp_[202] = 0; /* [7, 6] */
	jacp_[230] = 0; /* [8, 6] */
	jacp_[258] = 0; /* [9, 6] */
/* column 8 (df/dp_7) */
	jacp_[7] = 0; /* [0, 7] */
	jacp_[35] = 0; /* [1, 7] */
	jacp_[63] = 0; /* [2, 7] */
	jacp_[91] = 0; /* [3, 7] */
	jacp_[119] = 0; /* [4, 7] */
	jacp_[147] = 0; /* [5, 7] */
	jacp_[175] = 0; /* [6, 7] */
	jacp_[203] = 0; /* [7, 7] */
	jacp_[231] = 0; /* [8, 7] */
	jacp_[259] = 0; /* [9, 7] */
/* column 9 (df/dp_8) */
	jacp_[8] = (-1*(cAMP*Rii)); /* [0, 8] */
	jacp_[36] = 0; /* [1, 8] */
	jacp_[64] = (cAMP*Rii); /* [2, 8] */
	jacp_[92] = 0; /* [3, 8] */
	jacp_[120] = 0; /* [4, 8] */
	jacp_[148] = 0; /* [5, 8] */
	jacp_[176] = 0; /* [6, 8] */
	jacp_[204] = 0; /* [7, 8] */
	jacp_[232] = 0; /* [8, 8] */
	jacp_[260] = 0; /* [9, 8] */
/* column 10 (df/dp_9) */
	jacp_[9] = (-1*(0-Rii_cAMP)); /* [0, 9] */
	jacp_[37] = 0; /* [1, 9] */
	jacp_[65] = (0-Rii_cAMP); /* [2, 9] */
	jacp_[93] = 0; /* [3, 9] */
	jacp_[121] = 0; /* [4, 9] */
	jacp_[149] = 0; /* [5, 9] */
	jacp_[177] = 0; /* [6, 9] */
	jacp_[205] = 0; /* [7, 9] */
	jacp_[233] = 0; /* [8, 9] */
	jacp_[261] = 0; /* [9, 9] */
/* column 11 (df/dp_10) */
	jacp_[10] = 0; /* [0, 10] */
	jacp_[38] = (0-(Rii_C*cAMP)); /* [1, 10] */
	jacp_[66] = 0; /* [2, 10] */
	jacp_[94] = (Rii_C*cAMP); /* [3, 10] */
	jacp_[122] = 0; /* [4, 10] */
	jacp_[150] = 0; /* [5, 10] */
	jacp_[178] = 0; /* [6, 10] */
	jacp_[206] = 0; /* [7, 10] */
	jacp_[234] = 0; /* [8, 10] */
	jacp_[262] = 0; /* [9, 10] */
/* column 12 (df/dp_11) */
	jacp_[11] = 0; /* [0, 11] */
	jacp_[39] = (0-(0-Rii_C_cAMP)); /* [1, 11] */
	jacp_[67] = 0; /* [2, 11] */
	jacp_[95] = (0-Rii_C_cAMP); /* [3, 11] */
	jacp_[123] = 0; /* [4, 11] */
	jacp_[151] = 0; /* [5, 11] */
	jacp_[179] = 0; /* [6, 11] */
	jacp_[207] = 0; /* [7, 11] */
	jacp_[235] = 0; /* [8, 11] */
	jacp_[263] = 0; /* [9, 11] */
/* column 13 (df/dp_12) */
	jacp_[12] = (0-(Rii*C)); /* [0, 12] */
	jacp_[40] = (Rii*C); /* [1, 12] */
	jacp_[68] = 0; /* [2, 12] */
	jacp_[96] = 0; /* [3, 12] */
	jacp_[124] = 0; /* [4, 12] */
	jacp_[152] = 0; /* [5, 12] */
	jacp_[180] = 0; /* [6, 12] */
	jacp_[208] = 0; /* [7, 12] */
	jacp_[236] = 0; /* [8, 12] */
	jacp_[264] = 0; /* [9, 12] */
/* column 14 (df/dp_13) */
	jacp_[13] = 0; /* [0, 13] */
	jacp_[41] = 0; /* [1, 13] */
	jacp_[69] = (0-(Rii_cAMP*C)); /* [2, 13] */
	jacp_[97] = (Rii_cAMP*C); /* [3, 13] */
	jacp_[125] = 0; /* [4, 13] */
	jacp_[153] = 0; /* [5, 13] */
	jacp_[181] = 0; /* [6, 13] */
	jacp_[209] = 0; /* [7, 13] */
	jacp_[237] = 0; /* [8, 13] */
	jacp_[265] = 0; /* [9, 13] */
/* column 15 (df/dp_14) */
	jacp_[14] = 0; /* [0, 14] */
	jacp_[42] = 0; /* [1, 14] */
	jacp_[70] = (0-(0-Rii_C_cAMP)); /* [2, 14] */
	jacp_[98] = (0-Rii_C_cAMP); /* [3, 14] */
	jacp_[126] = 0; /* [4, 14] */
	jacp_[154] = 0; /* [5, 14] */
	jacp_[182] = 0; /* [6, 14] */
	jacp_[210] = 0; /* [7, 14] */
	jacp_[238] = 0; /* [8, 14] */
	jacp_[266] = 0; /* [9, 14] */
/* column 16 (df/dp_15) */
	jacp_[15] = 0; /* [0, 15] */
	jacp_[43] = 0; /* [1, 15] */
	jacp_[71] = 0; /* [2, 15] */
	jacp_[99] = (0-Rii_C_cAMP); /* [3, 15] */
	jacp_[127] = 0; /* [4, 15] */
	jacp_[155] = 0; /* [5, 15] */
	jacp_[183] = 0; /* [6, 15] */
	jacp_[211] = 0; /* [7, 15] */
	jacp_[239] = 0; /* [8, 15] */
	jacp_[267] = 0; /* [9, 15] */
/* column 17 (df/dp_16) */
	jacp_[16] = (0-(0-Rii_C)); /* [0, 16] */
	jacp_[44] = (0-Rii_C); /* [1, 16] */
	jacp_[72] = 0; /* [2, 16] */
	jacp_[100] = 0; /* [3, 16] */
	jacp_[128] = 0; /* [4, 16] */
	jacp_[156] = 0; /* [5, 16] */
	jacp_[184] = 0; /* [6, 16] */
	jacp_[212] = 0; /* [7, 16] */
	jacp_[240] = 0; /* [8, 16] */
	jacp_[268] = 0; /* [9, 16] */
/* column 18 (df/dp_17) */
	jacp_[17] = ((1-b_AKAP)*RiiP_CaN); /* [0, 17] */
	jacp_[45] = 0; /* [1, 17] */
	jacp_[73] = ((1-b_AKAP)*RiiP_cAMP_CaN); /* [2, 17] */
	jacp_[101] = 0; /* [3, 17] */
	jacp_[129] = ((1-b_AKAP)*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*(CaN*((-1*RiiP)-RiiP_cAMP)))+RiiP_CaN)+RiiP_cAMP_CaN)); /* [4, 17] */
	jacp_[157] = ((1-b_AKAP)*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*RiiP)*CaN)-RiiP_CaN)); /* [5, 17] */
	jacp_[185] = ((1-b_AKAP)*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*CaN)*RiiP_cAMP)-RiiP_cAMP_CaN)); /* [6, 17] */
	jacp_[213] = 0; /* [7, 17] */
	jacp_[241] = 0; /* [8, 17] */
	jacp_[269] = 0; /* [9, 17] */
/* column 19 (df/dp_18) */
	jacp_[18] = 0; /* [0, 18] */
	jacp_[46] = 0; /* [1, 18] */
	jacp_[74] = 0; /* [2, 18] */
	jacp_[102] = 0; /* [3, 18] */
	jacp_[130] = ((1-b_AKAP)*((-1*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*RiiP)*CaN)-RiiP_CaN))-(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*CaN)*RiiP_cAMP)-RiiP_cAMP_CaN))); /* [4, 18] */
	jacp_[158] = ((1-b_AKAP)*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*RiiP)*CaN)-RiiP_CaN)); /* [5, 18] */
	jacp_[186] = ((1-b_AKAP)*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*CaN)*RiiP_cAMP)-RiiP_cAMP_CaN)); /* [6, 18] */
	jacp_[214] = 0; /* [7, 18] */
	jacp_[242] = 0; /* [8, 18] */
	jacp_[270] = 0; /* [9, 18] */
/* column 20 (df/dp_19) */
	jacp_[19] = (b_AKAP*RiiP_CaN); /* [0, 19] */
	jacp_[47] = 0; /* [1, 19] */
	jacp_[75] = (b_AKAP*RiiP_cAMP_CaN); /* [2, 19] */
	jacp_[103] = 0; /* [3, 19] */
	jacp_[131] = (b_AKAP*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*(CaN*((-1*RiiP)-RiiP_cAMP)))+RiiP_CaN)+RiiP_cAMP_CaN)); /* [4, 19] */
	jacp_[159] = (b_AKAP*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*RiiP)*CaN)-RiiP_CaN)); /* [5, 19] */
	jacp_[187] = (b_AKAP*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*CaN)*RiiP_cAMP)-RiiP_cAMP_CaN)); /* [6, 19] */
	jacp_[215] = 0; /* [7, 19] */
	jacp_[243] = 0; /* [8, 19] */
	jacp_[271] = 0; /* [9, 19] */
/* column 21 (df/dp_20) */
	jacp_[20] = 0; /* [0, 20] */
	jacp_[48] = 0; /* [1, 20] */
	jacp_[76] = 0; /* [2, 20] */
	jacp_[104] = 0; /* [3, 20] */
	jacp_[132] = (b_AKAP*((-1*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*RiiP)*CaN)-RiiP_CaN))-(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*CaN)*RiiP_cAMP)-RiiP_cAMP_CaN))); /* [4, 20] */
	jacp_[160] = (b_AKAP*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*RiiP)*CaN)-RiiP_CaN)); /* [5, 20] */
	jacp_[188] = (b_AKAP*(((((b_AKAP*(1/kmON))+((1-b_AKAP)/kmOFF))*CaN)*RiiP_cAMP)-RiiP_cAMP_CaN)); /* [6, 20] */
	jacp_[216] = 0; /* [7, 20] */
	jacp_[244] = 0; /* [8, 20] */
	jacp_[272] = 0; /* [9, 20] */
/* column 22 (df/dp_21) */
	jacp_[21] = 0; /* [0, 21] */
	jacp_[49] = 0; /* [1, 21] */
	jacp_[77] = 0; /* [2, 21] */
	jacp_[105] = 0; /* [3, 21] */
	jacp_[133] = 0; /* [4, 21] */
	jacp_[161] = 0; /* [5, 21] */
	jacp_[189] = 0; /* [6, 21] */
	jacp_[217] = (-1*(C*AKAR4)); /* [7, 21] */
	jacp_[245] = (C*AKAR4); /* [8, 21] */
	jacp_[273] = 0; /* [9, 21] */
/* column 23 (df/dp_22) */
	jacp_[22] = 0; /* [0, 22] */
	jacp_[50] = 0; /* [1, 22] */
	jacp_[78] = 0; /* [2, 22] */
	jacp_[106] = 0; /* [3, 22] */
	jacp_[134] = 0; /* [4, 22] */
	jacp_[162] = 0; /* [5, 22] */
	jacp_[190] = 0; /* [6, 22] */
	jacp_[218] = (-1*(0-AKAR4_C)); /* [7, 22] */
	jacp_[246] = (0-AKAR4_C); /* [8, 22] */
	jacp_[274] = 0; /* [9, 22] */
/* column 24 (df/dp_23) */
	jacp_[23] = 0; /* [0, 23] */
	jacp_[51] = 0; /* [1, 23] */
	jacp_[79] = 0; /* [2, 23] */
	jacp_[107] = 0; /* [3, 23] */
	jacp_[135] = 0; /* [4, 23] */
	jacp_[163] = 0; /* [5, 23] */
	jacp_[191] = 0; /* [6, 23] */
	jacp_[219] = 0; /* [7, 23] */
	jacp_[247] = (0-AKAR4_C); /* [8, 23] */
	jacp_[275] = AKAR4_C; /* [9, 23] */
/* column 25 (df/dp_24) */
	jacp_[24] = 0; /* [0, 24] */
	jacp_[52] = 0; /* [1, 24] */
	jacp_[80] = 0; /* [2, 24] */
	jacp_[108] = 0; /* [3, 24] */
	jacp_[136] = ((((1-b_AKAP)*(0-(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))))/(kmOFF*kmOFF))*(CaN*((-1*RiiP)-RiiP_cAMP))); /* [4, 24] */
	jacp_[164] = (((((1-b_AKAP)*(0-(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))))/(kmOFF*kmOFF))*RiiP)*CaN); /* [5, 24] */
	jacp_[192] = (((((1-b_AKAP)*(0-(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))))/(kmOFF*kmOFF))*CaN)*RiiP_cAMP); /* [6, 24] */
	jacp_[220] = 0; /* [7, 24] */
	jacp_[248] = 0; /* [8, 24] */
	jacp_[276] = 0; /* [9, 24] */
/* column 26 (df/dp_25) */
	jacp_[25] = 0; /* [0, 25] */
	jacp_[53] = 0; /* [1, 25] */
	jacp_[81] = 0; /* [2, 25] */
	jacp_[109] = 0; /* [3, 25] */
	jacp_[137] = ((b_AKAP*((0-(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))/(kmON*kmON)))*(CaN*((-1*RiiP)-RiiP_cAMP))); /* [4, 25] */
	jacp_[165] = (((b_AKAP*((0-(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))/(kmON*kmON)))*RiiP)*CaN); /* [5, 25] */
	jacp_[193] = (((b_AKAP*((0-(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))/(kmON*kmON)))*CaN)*RiiP_cAMP); /* [6, 25] */
	jacp_[221] = 0; /* [7, 25] */
	jacp_[249] = 0; /* [8, 25] */
	jacp_[277] = 0; /* [9, 25] */
/* column 27 (df/dp_26) */
	jacp_[26] = 0; /* [0, 26] */
	jacp_[54] = 0; /* [1, 26] */
	jacp_[82] = 0; /* [2, 26] */
	jacp_[110] = 0; /* [3, 26] */
	jacp_[138] = 0; /* [4, 26] */
	jacp_[166] = 0; /* [5, 26] */
	jacp_[194] = 0; /* [6, 26] */
	jacp_[222] = 0; /* [7, 26] */
	jacp_[250] = 0; /* [8, 26] */
	jacp_[278] = 0; /* [9, 26] */
/* column 28 (df/dp_27) */
	jacp_[27] = ((AKAPon_1+(-1*AKAPoff_1))*RiiP_CaN); /* [0, 27] */
	jacp_[55] = 0; /* [1, 27] */
	jacp_[83] = ((AKAPon_1+(-1*AKAPoff_1))*RiiP_cAMP_CaN); /* [2, 27] */
	jacp_[111] = 0; /* [3, 27] */
	jacp_[139] = ((((-1*(((((((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))/kmON)+(b_AKAP*(((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))/kmON)))+(((-1*(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))+((1-b_AKAP)*((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))))/kmOFF))*RiiP)*CaN)-((AKAPon_3+(-1*AKAPoff_3))*RiiP_CaN)))-(((((((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))/kmON)+(b_AKAP*(((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))/kmON)))+(((-1*(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))+((1-b_AKAP)*((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))))/kmOFF))*CaN)*RiiP_cAMP)-((AKAPon_3+(-1*AKAPoff_3))*RiiP_cAMP_CaN)))+((AKAPon_1+(-1*AKAPoff_1))*RiiP_CaN))+((AKAPon_1+(-1*AKAPoff_1))*RiiP_cAMP_CaN)); /* [4, 27] */
	jacp_[167] = ((((((((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))/kmON)+(b_AKAP*(((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))/kmON)))+(((-1*(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))+((1-b_AKAP)*((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))))/kmOFF))*RiiP)*CaN)-((AKAPon_3+(-1*AKAPoff_3))*RiiP_CaN))-((AKAPon_1+(-1*AKAPoff_1))*RiiP_CaN)); /* [5, 27] */
	jacp_[195] = ((((((((((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1)))/kmON)+(b_AKAP*(((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))/kmON)))+(((-1*(((b_AKAP*AKAPon_3)+((1-b_AKAP)*AKAPoff_3))+((b_AKAP*AKAPon_1)+((1-b_AKAP)*AKAPoff_1))))+((1-b_AKAP)*((AKAPon_3+(-1*AKAPoff_3))+(AKAPon_1+(-1*AKAPoff_1)))))/kmOFF))*CaN)*RiiP_cAMP)-((AKAPon_3+(-1*AKAPoff_3))*RiiP_cAMP_CaN))-((AKAPon_1+(-1*AKAPoff_1))*RiiP_cAMP_CaN)); /* [6, 27] */
	jacp_[223] = 0; /* [7, 27] */
	jacp_[251] = 0; /* [8, 27] */
	jacp_[279] = 0; /* [9, 27] */
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
	double reaction_44_=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33_=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_4_8=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_3_7=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	func_[0] = (AKAR4p*5)*71.67+100;
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAP79_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 28;
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
	p_[16] = 3e-04;
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
	return GSL_SUCCESS;
}
/* ode initial values */
int AKAP79_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 16;
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
