#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* string.h for memset() */
enum stateVariable { _Rii,_cAMP,_RiiP,_Rii_C,_RiiP_cAMP,_RiiP_C,_RiiP_C_cAMP,_C,_Rii_cAMP,_Rii_C_cAMP,_CaN,_RiiP_CaN,_RiiP_cAMP_CaN,_AKAR4,_AKAR4_C,_AKAR4p, numStateVar }; /* state variable indexes  */
enum param { _kf_Rii_C__RiiP_C,_kf_RiiP_CxcAMP__RiiP_C_cAMP,_kf_RiiP_cAMPxC__RiiP_C_cAMP,_kb_RiiP_cAMPxC__RiiP_C_cAMP,_kb_RiiPXcAMP__RiiP_cAMP,_kf_RiiPXcAMP__RiiP_cAMP,_kf_RiiPxC__RiiP_C,_kb_RiiPxC__RiiP_C,_kf_cAMPxRii__Rii_cAMP,_kb_cAMPxRii__Rii_cAMP,_kf_Rii_CxcAMP__Rii_C_cAMP,_kb_Rii_CxcAMP__Rii_C_cAMP,_kf_RiixC__Rii_C,_kf_Rii_cAMPxC__Rii_C_cAMP,_kb_Rii_cAMPxC__Rii_C_cAMP,_kf_Rii_C_cAMP__RiiP_C_cAMP,_kb_RiixC__Rii_C,_AKAPoff_1,_AKAPoff_3,_AKAPon_1,_AKAPon_3,_kf_C_AKAR4,_kb_C_AKAR4,_kcat_AKARp,_kmOFF,_kmON,_KD_T,_b_AKAP, numParam }; /* parameter indexes  */
enum func { _AKAR4pOUT, numFunc }; /* parameter indexes  */

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
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	f_[_Rii] = -reaction_78-reaction_58+reaction_48; /* Rii */
	f_[_cAMP] = -reaction_12-reaction_43-reaction_78-reaction_56; /* cAMP */
	f_[_RiiP] = -reaction_14-reaction_43-reaction_44; /* RiiP */
	f_[_Rii_C] = -reaction_51-reaction_56+reaction_58; /* Rii_C */
	f_[_RiiP_cAMP] = +reaction_43-reaction_23-reaction_33; /* RiiP_cAMP */
	f_[_RiiP_C] = +reaction_51+reaction_14-reaction_12; /* RiiP_C */
	f_[_RiiP_C_cAMP] = +reaction_12+reaction_23+reaction_62; /* RiiP_C_cAMP */
	f_[_C] = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2; /* C */
	f_[_Rii_cAMP] = +reaction_78-reaction_76+reaction_37; /* Rii_cAMP */
	f_[_Rii_C_cAMP] = +reaction_56+reaction_76-reaction_62; /* Rii_C_cAMP */
	f_[_CaN] = -reaction_44-reaction_33+reaction_48+reaction_37; /* CaN */
	f_[_RiiP_CaN] = +reaction_44-reaction_48; /* RiiP_CaN */
	f_[_RiiP_cAMP_CaN] = +reaction_33-reaction_37; /* RiiP_cAMP_CaN */
	f_[_AKAR4] = -reaction_1; /* AKAR4 */
	f_[_AKAR4_C] = +reaction_1-reaction_2; /* AKAR4_C */
	f_[_AKAR4p] = +reaction_2; /* AKAR4p */
	return GSL_SUCCESS;
}
int AKAP79_netflux(double t, double y_[], double *flux, void *par){
	double *p_=par;
	if (!y_ || !flux) return 16;
	double kf_Rii_C__RiiP_C = p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP = p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP = p_[5];
	double kf_RiiPxC__RiiP_C = p_[6];
	double kb_RiiPxC__RiiP_C = p_[7];
	double kf_cAMPxRii__Rii_cAMP = p_[8];
	double kb_cAMPxRii__Rii_cAMP = p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[11];
	double kf_RiixC__Rii_C = p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[15];
	double kb_RiixC__Rii_C = p_[16];
	double AKAPoff_1 = p_[17];
	double AKAPoff_3 = p_[18];
	double AKAPon_1 = p_[19];
	double AKAPon_3 = p_[20];
	double kf_C_AKAR4 = p_[21];
	double kb_C_AKAR4 = p_[22];
	double kcat_AKARp = p_[23];
	double kmOFF = p_[24];
	double kmON = p_[25];
	double KD_T = p_[26];
	double b_AKAP = p_[27];
	double Rii = y_[0];
	double cAMP = y_[1];
	double RiiP = y_[2];
	double Rii_C = y_[3];
	double RiiP_cAMP = y_[4];
	double RiiP_C = y_[5];
	double RiiP_C_cAMP = y_[6];
	double C = y_[7];
	double Rii_cAMP = y_[8];
	double Rii_C_cAMP = y_[9];
	double CaN = y_[10];
	double RiiP_CaN = y_[11];
	double RiiP_cAMP_CaN = y_[12];
	double AKAR4 = y_[13];
	double AKAR4_C = y_[14];
	double AKAR4p = y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	flux[0] = kf_Rii_C__RiiP_C*Rii_C;
	flux[1] = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C;
	flux[2] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP;
	flux[3] = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP;
	flux[4] = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP;
	flux[5] = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP;
	flux[6] = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP;
	flux[7] = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP;
	flux[8] = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	flux[9] = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C;
	flux[10] = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	flux[11] = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	flux[12] = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	flux[13] = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	flux[14] = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	flux[15] = kcat_AKARp*AKAR4_C;
	return GSL_SUCCESS;
}

int AKAP79_fwdflux(double t, double y_[], double *flux, void *par){
	double *p_=par;
	if (!y_ || !flux) return 16;
	double kf_Rii_C__RiiP_C = p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP = p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP = p_[5];
	double kf_RiiPxC__RiiP_C = p_[6];
	double kb_RiiPxC__RiiP_C = p_[7];
	double kf_cAMPxRii__Rii_cAMP = p_[8];
	double kb_cAMPxRii__Rii_cAMP = p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[11];
	double kf_RiixC__Rii_C = p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[15];
	double kb_RiixC__Rii_C = p_[16];
	double AKAPoff_1 = p_[17];
	double AKAPoff_3 = p_[18];
	double AKAPon_1 = p_[19];
	double AKAPon_3 = p_[20];
	double kf_C_AKAR4 = p_[21];
	double kb_C_AKAR4 = p_[22];
	double kcat_AKARp = p_[23];
	double kmOFF = p_[24];
	double kmON = p_[25];
	double KD_T = p_[26];
	double b_AKAP = p_[27];
	double Rii = y_[0];
	double cAMP = y_[1];
	double RiiP = y_[2];
	double Rii_C = y_[3];
	double RiiP_cAMP = y_[4];
	double RiiP_C = y_[5];
	double RiiP_C_cAMP = y_[6];
	double C = y_[7];
	double Rii_cAMP = y_[8];
	double Rii_C_cAMP = y_[9];
	double CaN = y_[10];
	double RiiP_CaN = y_[11];
	double RiiP_cAMP_CaN = y_[12];
	double AKAR4 = y_[13];
	double AKAR4_C = y_[14];
	double AKAR4p = y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	// kf_Rii_C__RiiP_C*Rii_C
	flux[0] = kf_Rii_C__RiiP_C*Rii_C;
	// kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	flux[1] = kf_RiiPxC__RiiP_C*RiiP*C ;
	// kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	flux[2] = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP ;
	// kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	flux[3] = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP ;
	// kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	flux[4] = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C ;
	// kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	flux[5] = kf_cAMPxRii__Rii_cAMP*cAMP*Rii ;
	// kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	flux[6] = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP ;
	// kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	flux[7] = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C ;
	// kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	flux[8] = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP;
	// kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	flux[9] = kf_RiixC__Rii_C*Rii*C ;
	// kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	flux[10] = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN ;
	// kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	flux[11] = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP ;
	// kf_RiiP_CaN__RiixCaN*RiiP_CaN
	flux[12] = kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	// kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	flux[13] = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	// kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	flux[14] = kf_C_AKAR4*C*AKAR4 ;
	// kcat_AKARp*AKAR4_C
	flux[15] = kcat_AKARp*AKAR4_C;
	return GSL_SUCCESS;
}

int AKAP79_bwdflux(double t, double y_[], double *flux, void *par){
	double *p_=par;
	if (!y_ || !flux) return 16;
	double kf_Rii_C__RiiP_C = p_[0];
	double kf_RiiP_CxcAMP__RiiP_C_cAMP = p_[1];
	double kf_RiiP_cAMPxC__RiiP_C_cAMP = p_[2];
	double kb_RiiP_cAMPxC__RiiP_C_cAMP = p_[3];
	double kb_RiiPXcAMP__RiiP_cAMP = p_[4];
	double kf_RiiPXcAMP__RiiP_cAMP = p_[5];
	double kf_RiiPxC__RiiP_C = p_[6];
	double kb_RiiPxC__RiiP_C = p_[7];
	double kf_cAMPxRii__Rii_cAMP = p_[8];
	double kb_cAMPxRii__Rii_cAMP = p_[9];
	double kf_Rii_CxcAMP__Rii_C_cAMP = p_[10];
	double kb_Rii_CxcAMP__Rii_C_cAMP = p_[11];
	double kf_RiixC__Rii_C = p_[12];
	double kf_Rii_cAMPxC__Rii_C_cAMP = p_[13];
	double kb_Rii_cAMPxC__Rii_C_cAMP = p_[14];
	double kf_Rii_C_cAMP__RiiP_C_cAMP = p_[15];
	double kb_RiixC__Rii_C = p_[16];
	double AKAPoff_1 = p_[17];
	double AKAPoff_3 = p_[18];
	double AKAPon_1 = p_[19];
	double AKAPon_3 = p_[20];
	double kf_C_AKAR4 = p_[21];
	double kb_C_AKAR4 = p_[22];
	double kcat_AKARp = p_[23];
	double kmOFF = p_[24];
	double kmON = p_[25];
	double KD_T = p_[26];
	double b_AKAP = p_[27];
	double Rii = y_[0];
	double cAMP = y_[1];
	double RiiP = y_[2];
	double Rii_C = y_[3];
	double RiiP_cAMP = y_[4];
	double RiiP_C = y_[5];
	double RiiP_C_cAMP = y_[6];
	double C = y_[7];
	double Rii_cAMP = y_[8];
	double Rii_C_cAMP = y_[9];
	double CaN = y_[10];
	double RiiP_CaN = y_[11];
	double RiiP_cAMP_CaN = y_[12];
	double AKAR4 = y_[13];
	double AKAR4_C = y_[14];
	double AKAR4p = y_[15];
	double kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1;
	double kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3;
	double kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF;
	double kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF;
	double kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T;
	flux[0] = 0.0; // kf_Rii_C__RiiP_C*Rii_C
	flux[1] =  kb_RiiPxC__RiiP_C*RiiP_C; // kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C
	flux[2] =  kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP; // kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP
	flux[3] =  kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP; // kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP
	flux[4] =  kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP; // kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP
	flux[5] =  kb_cAMPxRii__Rii_cAMP*Rii_cAMP; // kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP
	flux[6] =  kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP; // kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP
	flux[7] =  kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP; // kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP
	flux[8] = 0.0; // kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP
	flux[9] =  kb_RiixC__Rii_C*Rii_C; // kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C
	flux[10] =  kb_RiiPxCaN__RiiP_CaN*RiiP_CaN; // kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN
	flux[11] =  kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN; // kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN
	flux[12] = 0.0; // kf_RiiP_CaN__RiixCaN*RiiP_CaN
	flux[13] = 0.0; // kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN
	flux[14] =  kb_C_AKAR4*AKAR4_C; // kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C
	flux[15] = 0.0; // kcat_AKARp*AKAR4_C
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
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(jac_,0,sizeof(double)*numStateVar*numStateVar); /* 256 */
/* column 1 (df/dy_0) */
	jac_[0] = (-cAMP*kf_cAMPxRii__Rii_cAMP)-C*kf_RiixC__Rii_C; /* [0, 0] */
	jac_[16] = -cAMP*kf_cAMPxRii__Rii_cAMP; /* [1, 0] */
	jac_[48] = C*kf_RiixC__Rii_C; /* [3, 0] */
	jac_[112] = -C*kf_RiixC__Rii_C; /* [7, 0] */
	jac_[128] = cAMP*kf_cAMPxRii__Rii_cAMP; /* [8, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = -Rii*kf_cAMPxRii__Rii_cAMP; /* [0, 1] */
	jac_[17] = (-Rii*kf_cAMPxRii__Rii_cAMP)-Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP-RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP-RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 1] */
	jac_[33] = -RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [2, 1] */
	jac_[49] = -Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP; /* [3, 1] */
	jac_[65] = RiiP*kf_RiiPXcAMP__RiiP_cAMP; /* [4, 1] */
	jac_[81] = -RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [5, 1] */
	jac_[97] = RiiP_C*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [6, 1] */
	jac_[129] = Rii*kf_cAMPxRii__Rii_cAMP; /* [8, 1] */
	jac_[145] = Rii_C*kf_Rii_CxcAMP__Rii_C_cAMP; /* [9, 1] */
/* column 3 (df/dy_2) */
	jac_[18] = -cAMP*kf_RiiPXcAMP__RiiP_cAMP; /* [1, 2] */
	jac_[34] = (-CaN*kf_RiiPxCaN__RiiP_CaN)-C*kf_RiiPxC__RiiP_C-cAMP*kf_RiiPXcAMP__RiiP_cAMP; /* [2, 2] */
	jac_[66] = cAMP*kf_RiiPXcAMP__RiiP_cAMP; /* [4, 2] */
	jac_[82] = C*kf_RiiPxC__RiiP_C; /* [5, 2] */
	jac_[114] = -C*kf_RiiPxC__RiiP_C; /* [7, 2] */
	jac_[162] = -CaN*kf_RiiPxCaN__RiiP_CaN; /* [10, 2] */
	jac_[178] = CaN*kf_RiiPxCaN__RiiP_CaN; /* [11, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = kb_RiixC__Rii_C; /* [0, 3] */
	jac_[19] = -cAMP*kf_Rii_CxcAMP__Rii_C_cAMP; /* [1, 3] */
	jac_[51] = (-cAMP*kf_Rii_CxcAMP__Rii_C_cAMP)-kf_Rii_C__RiiP_C-kb_RiixC__Rii_C; /* [3, 3] */
	jac_[83] = kf_Rii_C__RiiP_C; /* [5, 3] */
	jac_[115] = kb_RiixC__Rii_C; /* [7, 3] */
	jac_[147] = cAMP*kf_Rii_CxcAMP__Rii_C_cAMP; /* [9, 3] */
/* column 5 (df/dy_4) */
	jac_[20] = kb_RiiPXcAMP__RiiP_cAMP; /* [1, 4] */
	jac_[36] = kb_RiiPXcAMP__RiiP_cAMP; /* [2, 4] */
	jac_[68] = (-C*kf_RiiP_cAMPxC__RiiP_C_cAMP)-CaN*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN-kb_RiiPXcAMP__RiiP_cAMP; /* [4, 4] */
	jac_[100] = C*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [6, 4] */
	jac_[116] = -C*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [7, 4] */
	jac_[164] = -CaN*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [10, 4] */
	jac_[196] = CaN*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [12, 4] */
/* column 6 (df/dy_5) */
	jac_[21] = -cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [1, 5] */
	jac_[37] = kb_RiiPxC__RiiP_C; /* [2, 5] */
	jac_[85] = (-cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP)-kb_RiiPxC__RiiP_C; /* [5, 5] */
	jac_[101] = cAMP*kf_RiiP_CxcAMP__RiiP_C_cAMP; /* [6, 5] */
	jac_[117] = kb_RiiPxC__RiiP_C; /* [7, 5] */
/* column 7 (df/dy_6) */
	jac_[22] = kb_RiiP_CxcAMP__RiiP_C_cAMP; /* [1, 6] */
	jac_[70] = kb_RiiP_cAMPxC__RiiP_C_cAMP; /* [4, 6] */
	jac_[86] = kb_RiiP_CxcAMP__RiiP_C_cAMP; /* [5, 6] */
	jac_[102] = (-kb_RiiP_cAMPxC__RiiP_C_cAMP)-kb_RiiP_CxcAMP__RiiP_C_cAMP; /* [6, 6] */
	jac_[118] = kb_RiiP_cAMPxC__RiiP_C_cAMP; /* [7, 6] */
/* column 8 (df/dy_7) */
	jac_[7] = -Rii*kf_RiixC__Rii_C; /* [0, 7] */
	jac_[39] = -RiiP*kf_RiiPxC__RiiP_C; /* [2, 7] */
	jac_[55] = Rii*kf_RiixC__Rii_C; /* [3, 7] */
	jac_[71] = -RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [4, 7] */
	jac_[87] = RiiP*kf_RiiPxC__RiiP_C; /* [5, 7] */
	jac_[103] = RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP; /* [6, 7] */
	jac_[119] = (-Rii*kf_RiixC__Rii_C)-Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP-RiiP*kf_RiiPxC__RiiP_C-RiiP_cAMP*kf_RiiP_cAMPxC__RiiP_C_cAMP-AKAR4*kf_C_AKAR4; /* [7, 7] */
	jac_[135] = -Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP; /* [8, 7] */
	jac_[151] = Rii_cAMP*kf_Rii_cAMPxC__Rii_C_cAMP; /* [9, 7] */
	jac_[215] = -AKAR4*kf_C_AKAR4; /* [13, 7] */
	jac_[231] = AKAR4*kf_C_AKAR4; /* [14, 7] */
/* column 9 (df/dy_8) */
	jac_[8] = kb_cAMPxRii__Rii_cAMP; /* [0, 8] */
	jac_[24] = kb_cAMPxRii__Rii_cAMP; /* [1, 8] */
	jac_[120] = -C*kf_Rii_cAMPxC__Rii_C_cAMP; /* [7, 8] */
	jac_[136] = (-C*kf_Rii_cAMPxC__Rii_C_cAMP)-kb_cAMPxRii__Rii_cAMP; /* [8, 8] */
	jac_[152] = C*kf_Rii_cAMPxC__Rii_C_cAMP; /* [9, 8] */
/* column 10 (df/dy_9) */
	jac_[25] = kb_Rii_CxcAMP__Rii_C_cAMP; /* [1, 9] */
	jac_[57] = kb_Rii_CxcAMP__Rii_C_cAMP; /* [3, 9] */
	jac_[105] = kf_Rii_C_cAMP__RiiP_C_cAMP; /* [6, 9] */
	jac_[121] = kb_Rii_cAMPxC__Rii_C_cAMP; /* [7, 9] */
	jac_[137] = kb_Rii_cAMPxC__Rii_C_cAMP; /* [8, 9] */
	jac_[153] = (-kf_Rii_C_cAMP__RiiP_C_cAMP)-kb_Rii_cAMPxC__Rii_C_cAMP-kb_Rii_CxcAMP__Rii_C_cAMP; /* [9, 9] */
/* column 11 (df/dy_10) */
	jac_[42] = -RiiP*kf_RiiPxCaN__RiiP_CaN; /* [2, 10] */
	jac_[74] = -RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [4, 10] */
	jac_[170] = (-RiiP*kf_RiiPxCaN__RiiP_CaN)-RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [10, 10] */
	jac_[186] = RiiP*kf_RiiPxCaN__RiiP_CaN; /* [11, 10] */
	jac_[202] = RiiP_cAMP*kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [12, 10] */
/* column 12 (df/dy_11) */
	jac_[11] = kf_RiiP_CaN__RiixCaN; /* [0, 11] */
	jac_[43] = kb_RiiPxCaN__RiiP_CaN; /* [2, 11] */
	jac_[171] = kf_RiiP_CaN__RiixCaN+kb_RiiPxCaN__RiiP_CaN; /* [10, 11] */
	jac_[187] = (-kf_RiiP_CaN__RiixCaN)-kb_RiiPxCaN__RiiP_CaN; /* [11, 11] */
/* column 13 (df/dy_12) */
	jac_[76] = kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [4, 12] */
	jac_[140] = kf_RiiP_cAMP_CaN__CaNXRii_cAMP; /* [8, 12] */
	jac_[172] = kf_RiiP_cAMP_CaN__CaNXRii_cAMP+kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [10, 12] */
	jac_[204] = (-kf_RiiP_cAMP_CaN__CaNXRii_cAMP)-kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN; /* [12, 12] */
/* column 14 (df/dy_13) */
	jac_[125] = -C*kf_C_AKAR4; /* [7, 13] */
	jac_[221] = -C*kf_C_AKAR4; /* [13, 13] */
	jac_[237] = C*kf_C_AKAR4; /* [14, 13] */
/* column 15 (df/dy_14) */
	jac_[126] = kcat_AKARp+kb_C_AKAR4; /* [7, 14] */
	jac_[222] = kb_C_AKAR4; /* [13, 14] */
	jac_[238] = (-kcat_AKARp)-kb_C_AKAR4; /* [14, 14] */
	jac_[254] = kcat_AKARp; /* [15, 14] */
/* column 16 (df/dy_15) */
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
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(jacp_,0,sizeof(double)*numStateVar*numParam); /* 448 */
/* column 1 (df/dp_0) */
	jacp_[84] = -Rii_C; /* [3, 0] */
	jacp_[140] = Rii_C; /* [5, 0] */
/* column 2 (df/dp_1) */
	jacp_[29] = -RiiP_C*cAMP; /* [1, 1] */
	jacp_[141] = -RiiP_C*cAMP; /* [5, 1] */
	jacp_[169] = RiiP_C*cAMP; /* [6, 1] */
/* column 3 (df/dp_2) */
	jacp_[114] = -C*RiiP_cAMP; /* [4, 2] */
	jacp_[170] = C*RiiP_cAMP; /* [6, 2] */
	jacp_[198] = -C*RiiP_cAMP; /* [7, 2] */
/* column 4 (df/dp_3) */
	jacp_[115] = RiiP_C_cAMP; /* [4, 3] */
	jacp_[171] = -RiiP_C_cAMP; /* [6, 3] */
	jacp_[199] = RiiP_C_cAMP; /* [7, 3] */
/* column 5 (df/dp_4) */
	jacp_[32] = RiiP_cAMP; /* [1, 4] */
	jacp_[60] = RiiP_cAMP; /* [2, 4] */
	jacp_[116] = -RiiP_cAMP; /* [4, 4] */
/* column 6 (df/dp_5) */
	jacp_[33] = -RiiP*cAMP; /* [1, 5] */
	jacp_[61] = -RiiP*cAMP; /* [2, 5] */
	jacp_[117] = RiiP*cAMP; /* [4, 5] */
/* column 7 (df/dp_6) */
	jacp_[62] = -C*RiiP; /* [2, 6] */
	jacp_[146] = C*RiiP; /* [5, 6] */
	jacp_[202] = -C*RiiP; /* [7, 6] */
/* column 8 (df/dp_7) */
	jacp_[63] = RiiP_C; /* [2, 7] */
	jacp_[147] = -RiiP_C; /* [5, 7] */
	jacp_[203] = RiiP_C; /* [7, 7] */
/* column 9 (df/dp_8) */
	jacp_[8] = -Rii*cAMP; /* [0, 8] */
	jacp_[36] = -Rii*cAMP; /* [1, 8] */
	jacp_[232] = Rii*cAMP; /* [8, 8] */
/* column 10 (df/dp_9) */
	jacp_[9] = Rii_cAMP; /* [0, 9] */
	jacp_[37] = Rii_cAMP; /* [1, 9] */
	jacp_[233] = -Rii_cAMP; /* [8, 9] */
/* column 11 (df/dp_10) */
	jacp_[38] = -Rii_C*cAMP; /* [1, 10] */
	jacp_[94] = -Rii_C*cAMP; /* [3, 10] */
	jacp_[262] = Rii_C*cAMP; /* [9, 10] */
/* column 12 (df/dp_11) */
	jacp_[39] = Rii_C_cAMP; /* [1, 11] */
	jacp_[95] = Rii_C_cAMP; /* [3, 11] */
	jacp_[263] = -Rii_C_cAMP; /* [9, 11] */
/* column 13 (df/dp_12) */
	jacp_[12] = -C*Rii; /* [0, 12] */
	jacp_[96] = C*Rii; /* [3, 12] */
	jacp_[208] = -C*Rii; /* [7, 12] */
/* column 14 (df/dp_13) */
	jacp_[209] = -C*Rii_cAMP; /* [7, 13] */
	jacp_[237] = -C*Rii_cAMP; /* [8, 13] */
	jacp_[265] = C*Rii_cAMP; /* [9, 13] */
/* column 15 (df/dp_14) */
	jacp_[210] = Rii_C_cAMP; /* [7, 14] */
	jacp_[238] = Rii_C_cAMP; /* [8, 14] */
	jacp_[266] = -Rii_C_cAMP; /* [9, 14] */
/* column 16 (df/dp_15) */
	jacp_[183] = Rii_C_cAMP; /* [6, 15] */
	jacp_[267] = -Rii_C_cAMP; /* [9, 15] */
/* column 17 (df/dp_16) */
	jacp_[16] = Rii_C; /* [0, 16] */
	jacp_[100] = -Rii_C; /* [3, 16] */
	jacp_[212] = Rii_C; /* [7, 16] */
/* column 18 (df/dp_17) */
/* column 19 (df/dp_18) */
/* column 20 (df/dp_19) */
/* column 21 (df/dp_20) */
/* column 22 (df/dp_21) */
	jacp_[217] = -AKAR4*C; /* [7, 21] */
	jacp_[385] = -AKAR4*C; /* [13, 21] */
	jacp_[413] = AKAR4*C; /* [14, 21] */
/* column 23 (df/dp_22) */
	jacp_[218] = AKAR4_C; /* [7, 22] */
	jacp_[386] = AKAR4_C; /* [13, 22] */
	jacp_[414] = -AKAR4_C; /* [14, 22] */
/* column 24 (df/dp_23) */
	jacp_[219] = AKAR4_C; /* [7, 23] */
	jacp_[415] = -AKAR4_C; /* [14, 23] */
	jacp_[443] = AKAR4_C; /* [15, 23] */
/* column 25 (df/dp_24) */
/* column 26 (df/dp_25) */
/* column 27 (df/dp_26) */
/* column 28 (df/dp_27) */
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
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	func_[_AKAR4pOUT] = (AKAR4p*5)*71.67+100; /* AKAR4pOUT */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int AKAP79_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 16;
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
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(funcJac_,0,sizeof(double)*numFunc*numStateVar); /* 16 */
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
	funcJac_[15] = 358.35; /* [0, 15] */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int AKAP79_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 28;
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
	double reaction_44=kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN;
	double reaction_33=kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN;
	double reaction_48=kf_RiiP_CaN__RiixCaN*RiiP_CaN;
	double reaction_37=kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN;
	double reaction_1=kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C;
	double reaction_2=kcat_AKARp*AKAR4_C;
	memset(funcJacp_,0,sizeof(double)*numFunc*numParam); /* 28 */
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
/* column 17 (dF/dp_16) */
/* column 18 (dF/dp_17) */
/* column 19 (dF/dp_18) */
/* column 20 (dF/dp_19) */
/* column 21 (dF/dp_20) */
/* column 22 (dF/dp_21) */
/* column 23 (dF/dp_22) */
/* column 24 (dF/dp_23) */
/* column 25 (dF/dp_24) */
/* column 26 (dF/dp_25) */
/* column 27 (dF/dp_26) */
/* column 28 (dF/dp_27) */
	return GSL_SUCCESS;
}
/* ode default parameters */
int AKAP79_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 28;
	memset(p_,0,sizeof(double)*numParam);
	p_[_kf_Rii_C__RiiP_C] = 33;
	p_[_kf_RiiP_CxcAMP__RiiP_C_cAMP] = 0.496;
	p_[_kf_RiiP_cAMPxC__RiiP_C_cAMP] = 0.00545;
	p_[_kb_RiiP_cAMPxC__RiiP_C_cAMP] = 0.0156;
	p_[_kb_RiiPXcAMP__RiiP_cAMP] = 0.0016;
	p_[_kf_RiiPXcAMP__RiiP_cAMP] = 0.015;
	p_[_kf_RiiPxC__RiiP_C] = 0.038;
	p_[_kb_RiiPxC__RiiP_C] = 0.0026;
	p_[_kf_cAMPxRii__Rii_cAMP] = 0.015;
	p_[_kb_cAMPxRii__Rii_cAMP] = 0.0016;
	p_[_kf_Rii_CxcAMP__Rii_C_cAMP] = 0.496;
	p_[_kb_Rii_CxcAMP__Rii_C_cAMP] = 1.413;
	p_[_kf_RiixC__Rii_C] = 2.1;
	p_[_kf_Rii_cAMPxC__Rii_C_cAMP] = 0.2984;
	p_[_kb_Rii_cAMPxC__Rii_C_cAMP] = 0.018;
	p_[_kf_Rii_C_cAMP__RiiP_C_cAMP] = 33;
	p_[_kb_RiixC__Rii_C] = 3e-04;
	p_[_AKAPoff_1] = 2.6;
	p_[_AKAPoff_3] = 20;
	p_[_AKAPon_1] = 0.45;
	p_[_AKAPon_3] = 2;
	p_[_kf_C_AKAR4] = 0.018;
	p_[_kb_C_AKAR4] = 0.106;
	p_[_kcat_AKARp] = 10.2;
	p_[_kmOFF] = 100;
	p_[_kmON] = 1;
	p_[_KD_T] = 0.7;
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
	memset(y_,0,sizeof(double)*numStateVar);
	y_[_Rii] = 6.3;
	y_[_Rii_C] = 0.63;
	y_[_CaN] = 1.5;
	y_[_AKAR4] = 0.2;
	return GSL_SUCCESS;
}
