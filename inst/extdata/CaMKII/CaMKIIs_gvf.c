#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* string.h for memset() */
enum stateVariable { _CaM_Ca1,_CaM_Ca2,_CaM_Ca3,_CaM_Ca4,_PP2B_CaM,_PP2B_CaM_Ca1,_PP2B_CaM_Ca2,_PP2B_CaM_Ca3,_PP2B_CaM_Ca4,_CaMKII_CaM,_CaMKII_CaM_Ca1,_CaMKII_CaM_Ca2,_CaMKII_CaM_Ca3,_CaMKII_CaM_Ca4,_pCaMKII_CaM_Ca4,_pCaMKIIa,_pCaMKII_CaM_Ca3,_pCaMKII_CaM_Ca2,_pCaMKII_CaM_Ca1,_pCaMKII_CaM,_PP1__pCaMKIIa,_caa,_cab, numStateVar }; /* state variable indexes  */
enum param { _kf__CaM__Ca,_kf__CaM_Ca1__Ca,_kf__CaM_Ca2__Ca,_kf__CaM_Ca3__Ca,_kf__CaM__PP2B,_kf__CaM_Ca1__PP2B,_kf__CaM_Ca2__PP2B,_kf__CaM_Ca3__PP2B,_kf__CaM_Ca4__PP2B,_kf__PP2B_CaM__Ca,_kf__PP2B_CaM_Ca1__Ca,_kf__PP2B_CaM_Ca2__Ca,_kf__PP2B_CaM_Ca3__Ca,_KD__CaM_Ca3__Ca,_KD__CaM_Ca2__Ca,_KD__CaM_Ca1__Ca,_KD__CaM__Ca,_KD__CaM_Ca4__PP2B,_KD__PP2B_CaM_Ca3__Ca,_KD__PP2B_CaM_Ca2__Ca,_KD__PP2B_CaM_Ca1__Ca,_KD__PP2B_CaM__Ca,_kf__CaM__CaMKII,_kf__CaMKII_CaM_Ca3__Ca,_kf__CaMKII_CaM_Ca2__Ca,_kf__CaMKII_CaM_Ca1__Ca,_kf__CaMKII_CaM__Ca,_kf__CaM_Ca1__CaMKII,_kf__CaM_Ca2__CaMKII,_kf__CaM_Ca3__CaMKII,_kf__CaM_Ca4__CaMKII,_KD__CaM_Ca4__CaMKII,_KD__CaMKII_CaM_Ca3__Ca,_KD__CaMKII_CaM_Ca2__Ca,_KD__CaMKII_CaM_Ca1__Ca,_KD__CaMKII_CaM__Ca,_kf__pCaMKII_CaM_Ca3__Ca,_kf__CaM__pCaMKIIa,_kf__CaM_Ca1__pCaMKIIa,_kf__CaM_Ca2__pCaMKIIa,_kf__CaM_Ca3__pCaMKIIa,_kf__pCaMKII_CaM_Ca2__Ca,_kf__pCaMKII_CaM_Ca1__Ca,_kf__CaM_Ca4__pCaMKIIa,_kf__pCaMKII_CaM__Ca,_KD__pCaMKII_CaM_Ca3__Ca,_KD__pCaMKII_CaM_Ca2__Ca,_KD__pCaMKII_CaM_Ca1__Ca,_KD__pCaMKII_CaM__Ca,_KD__CaM_Ca4__pCaMKIIa,_kautMax,_kf__PP1__pCaMKIIa,_kr__PP1__pCaMKIIa,_kcat__PP1__pCaMKIIa,_Ca_set,_PP1_0,_CaMKII_0,_CaM_0,_PP2B_0,_kca1,_kca2,_isOn, numParam }; /* parameter indexes  */
enum eventLabel { CaSpike,Activation, numEvents }; /* event name indexes */
enum func { _CaPerCaM,_AutoCaMKII,_CaMPerPP2B,_ActivePP2B,_Camonitor, numFunc }; /* parameter indexes  */

/* The error code indicates how to pre-allocate memory
 * for output values such as `f_`. The _vf function returns
 * the number of state variables, if any of the args are `NULL`.
 * evaluation errors can be indicated by negative return values.
 * GSL_SUCCESS (0) is returned when no error occurred.
 */

/* ode vector field: y'=f(t,y;p) */
int CaMKIIs_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 23;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	f_[_CaM_Ca1] = +ReactionFlux1-ReactionFlux2-ReactionFlux6-ReactionFlux15-ReactionFlux26; /* CaM_Ca1 */
	f_[_CaM_Ca2] = +ReactionFlux2-ReactionFlux3-ReactionFlux7-ReactionFlux16-ReactionFlux25; /* CaM_Ca2 */
	f_[_CaM_Ca3] = +ReactionFlux3-ReactionFlux4-ReactionFlux8-ReactionFlux17-ReactionFlux24; /* CaM_Ca3 */
	f_[_CaM_Ca4] = +ReactionFlux4-ReactionFlux9-ReactionFlux18-ReactionFlux23; /* CaM_Ca4 */
	f_[_PP2B_CaM] = +ReactionFlux5-ReactionFlux10; /* PP2B_CaM */
	f_[_PP2B_CaM_Ca1] = +ReactionFlux6+ReactionFlux10-ReactionFlux11; /* PP2B_CaM_Ca1 */
	f_[_PP2B_CaM_Ca2] = +ReactionFlux7+ReactionFlux11-ReactionFlux12; /* PP2B_CaM_Ca2 */
	f_[_PP2B_CaM_Ca3] = +ReactionFlux8+ReactionFlux12-ReactionFlux13; /* PP2B_CaM_Ca3 */
	f_[_PP2B_CaM_Ca4] = +ReactionFlux9+ReactionFlux13; /* PP2B_CaM_Ca4 */
	f_[_CaMKII_CaM] = +ReactionFlux14-ReactionFlux19; /* CaMKII_CaM */
	f_[_CaMKII_CaM_Ca1] = +ReactionFlux15+ReactionFlux19-ReactionFlux20; /* CaMKII_CaM_Ca1 */
	f_[_CaMKII_CaM_Ca2] = +ReactionFlux16+ReactionFlux20-ReactionFlux21; /* CaMKII_CaM_Ca2 */
	f_[_CaMKII_CaM_Ca3] = +ReactionFlux17+ReactionFlux21-ReactionFlux22; /* CaMKII_CaM_Ca3 */
	f_[_CaMKII_CaM_Ca4] = +ReactionFlux18+ReactionFlux22-ReactionFlux32; /* CaMKII_CaM_Ca4 */
	f_[_pCaMKII_CaM_Ca4] = +ReactionFlux23+ReactionFlux31+ReactionFlux32; /* pCaMKII_CaM_Ca4 */
	f_[_pCaMKIIa] = -ReactionFlux23-ReactionFlux24-ReactionFlux25-ReactionFlux26-ReactionFlux27-ReactionFlux33; /* pCaMKIIa */
	f_[_pCaMKII_CaM_Ca3] = +ReactionFlux24+ReactionFlux30-ReactionFlux31; /* pCaMKII_CaM_Ca3 */
	f_[_pCaMKII_CaM_Ca2] = +ReactionFlux25+ReactionFlux29-ReactionFlux30; /* pCaMKII_CaM_Ca2 */
	f_[_pCaMKII_CaM_Ca1] = +ReactionFlux26+ReactionFlux28-ReactionFlux29; /* pCaMKII_CaM_Ca1 */
	f_[_pCaMKII_CaM] = +ReactionFlux27-ReactionFlux28; /* pCaMKII_CaM */
	f_[_PP1__pCaMKIIa] = +ReactionFlux33-ReactionFlux34; /* PP1__pCaMKIIa */
	f_[_caa] = +SpikeFlux1; /* caa */
	f_[_cab] = +SpikeFlux2; /* cab */
	return GSL_SUCCESS;
}
/* Scheduled Event function,
   EventLabel specifies which of the possible transformations to apply,
   dose can specify a scalar intensity for this transformation. */
int CaMKIIs_event(double t, double y_[], void *par, int EventLabel, double dose)
{
	double *p_=par;
	if (!y_ || !par || EventLabel<0) return 2;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	switch(EventLabel){
	case CaSpike:
		y_[_cab] = 5.658469984; /* state variable transformation */
	break;
	case Activation:
		p_[_isOn] = 1; /* parameter transformation */
	break;
	}
	return GSL_SUCCESS;
}

/* ode Jacobian df(t,y;p)/dy */
int CaMKIIs_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 23*23;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	memset(jac_,0,sizeof(double)*numStateVar*numStateVar); /* 529 */
/* column 1 (df/dy_0) */
	jac_[0] = (-kf__CaM_Ca1__pCaMKIIa*pCaMKIIa)-kr__CaM__Ca-PP2B*kf__CaM_Ca1__PP2B-CaMKII*kf__CaM_Ca1__CaMKII-Ca*kf__CaM_Ca1__Ca; /* [0, 0] */
	jac_[23] = Ca*kf__CaM_Ca1__Ca; /* [1, 0] */
	jac_[115] = PP2B*kf__CaM_Ca1__PP2B; /* [5, 0] */
	jac_[230] = CaMKII*kf__CaM_Ca1__CaMKII; /* [10, 0] */
	jac_[345] = -kf__CaM_Ca1__pCaMKIIa*pCaMKIIa; /* [15, 0] */
	jac_[414] = kf__CaM_Ca1__pCaMKIIa*pCaMKIIa; /* [18, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = kr__CaM_Ca1__Ca; /* [0, 1] */
	jac_[24] = (-kf__CaM_Ca2__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca1__Ca-PP2B*kf__CaM_Ca2__PP2B-CaMKII*kf__CaM_Ca2__CaMKII-Ca*kf__CaM_Ca2__Ca; /* [1, 1] */
	jac_[47] = Ca*kf__CaM_Ca2__Ca; /* [2, 1] */
	jac_[139] = PP2B*kf__CaM_Ca2__PP2B; /* [6, 1] */
	jac_[254] = CaMKII*kf__CaM_Ca2__CaMKII; /* [11, 1] */
	jac_[346] = -kf__CaM_Ca2__pCaMKIIa*pCaMKIIa; /* [15, 1] */
	jac_[392] = kf__CaM_Ca2__pCaMKIIa*pCaMKIIa; /* [17, 1] */
/* column 3 (df/dy_2) */
	jac_[25] = kr__CaM_Ca2__Ca; /* [1, 2] */
	jac_[48] = (-kf__CaM_Ca3__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca2__Ca-PP2B*kf__CaM_Ca3__PP2B-CaMKII*kf__CaM_Ca3__CaMKII-Ca*kf__CaM_Ca3__Ca; /* [2, 2] */
	jac_[71] = Ca*kf__CaM_Ca3__Ca; /* [3, 2] */
	jac_[163] = PP2B*kf__CaM_Ca3__PP2B; /* [7, 2] */
	jac_[278] = CaMKII*kf__CaM_Ca3__CaMKII; /* [12, 2] */
	jac_[347] = -kf__CaM_Ca3__pCaMKIIa*pCaMKIIa; /* [15, 2] */
	jac_[370] = kf__CaM_Ca3__pCaMKIIa*pCaMKIIa; /* [16, 2] */
/* column 4 (df/dy_3) */
	jac_[49] = kr__CaM_Ca3__Ca; /* [2, 3] */
	jac_[72] = (-kf__CaM_Ca4__pCaMKIIa*pCaMKIIa)-kr__CaM_Ca3__Ca-PP2B*kf__CaM_Ca4__PP2B-CaMKII*kf__CaM_Ca4__CaMKII; /* [3, 3] */
	jac_[187] = PP2B*kf__CaM_Ca4__PP2B; /* [8, 3] */
	jac_[302] = CaMKII*kf__CaM_Ca4__CaMKII; /* [13, 3] */
	jac_[325] = kf__CaM_Ca4__pCaMKIIa*pCaMKIIa; /* [14, 3] */
	jac_[348] = -kf__CaM_Ca4__pCaMKIIa*pCaMKIIa; /* [15, 3] */
/* column 5 (df/dy_4) */
	jac_[96] = (-kr__CaM__PP2B)-Ca*kf__PP2B_CaM__Ca; /* [4, 4] */
	jac_[119] = Ca*kf__PP2B_CaM__Ca; /* [5, 4] */
/* column 6 (df/dy_5) */
	jac_[5] = kr__CaM_Ca1__PP2B; /* [0, 5] */
	jac_[97] = kr__PP2B_CaM__Ca; /* [4, 5] */
	jac_[120] = (-kr__PP2B_CaM__Ca)-kr__CaM_Ca1__PP2B-Ca*kf__PP2B_CaM_Ca1__Ca; /* [5, 5] */
	jac_[143] = Ca*kf__PP2B_CaM_Ca1__Ca; /* [6, 5] */
/* column 7 (df/dy_6) */
	jac_[29] = kr__CaM_Ca2__PP2B; /* [1, 6] */
	jac_[121] = kr__PP2B_CaM_Ca1__Ca; /* [5, 6] */
	jac_[144] = (-kr__PP2B_CaM_Ca1__Ca)-kr__CaM_Ca2__PP2B-Ca*kf__PP2B_CaM_Ca2__Ca; /* [6, 6] */
	jac_[167] = Ca*kf__PP2B_CaM_Ca2__Ca; /* [7, 6] */
/* column 8 (df/dy_7) */
	jac_[53] = kr__CaM_Ca3__PP2B; /* [2, 7] */
	jac_[145] = kr__PP2B_CaM_Ca2__Ca; /* [6, 7] */
	jac_[168] = (-kr__PP2B_CaM_Ca2__Ca)-kr__CaM_Ca3__PP2B-Ca*kf__PP2B_CaM_Ca3__Ca; /* [7, 7] */
	jac_[191] = Ca*kf__PP2B_CaM_Ca3__Ca; /* [8, 7] */
/* column 9 (df/dy_8) */
	jac_[77] = kr__CaM_Ca4__PP2B; /* [3, 8] */
	jac_[169] = kr__PP2B_CaM_Ca3__Ca; /* [7, 8] */
	jac_[192] = (-kr__PP2B_CaM_Ca3__Ca)-kr__CaM_Ca4__PP2B; /* [8, 8] */
/* column 10 (df/dy_9) */
	jac_[216] = (-kr__CaM__CaMKII)-Ca*kf__CaMKII_CaM__Ca; /* [9, 9] */
	jac_[239] = Ca*kf__CaMKII_CaM__Ca; /* [10, 9] */
/* column 11 (df/dy_10) */
	jac_[10] = kr__CaM_Ca1__CaMKII; /* [0, 10] */
	jac_[217] = kr__CaMKII_CaM__Ca; /* [9, 10] */
	jac_[240] = (-kr__CaM_Ca1__CaMKII)-kr__CaMKII_CaM__Ca-Ca*kf__CaMKII_CaM_Ca1__Ca; /* [10, 10] */
	jac_[263] = Ca*kf__CaMKII_CaM_Ca1__Ca; /* [11, 10] */
/* column 12 (df/dy_11) */
	jac_[34] = kr__CaM_Ca2__CaMKII; /* [1, 11] */
	jac_[241] = kr__CaMKII_CaM_Ca1__Ca; /* [10, 11] */
	jac_[264] = (-kr__CaM_Ca2__CaMKII)-kr__CaMKII_CaM_Ca1__Ca-Ca*kf__CaMKII_CaM_Ca2__Ca; /* [11, 11] */
	jac_[287] = Ca*kf__CaMKII_CaM_Ca2__Ca; /* [12, 11] */
/* column 13 (df/dy_12) */
	jac_[58] = kr__CaM_Ca3__CaMKII; /* [2, 12] */
	jac_[265] = kr__CaMKII_CaM_Ca2__Ca; /* [11, 12] */
	jac_[288] = (-kr__CaM_Ca3__CaMKII)-kr__CaMKII_CaM_Ca2__Ca-Ca*kf__CaMKII_CaM_Ca3__Ca; /* [12, 12] */
	jac_[311] = Ca*kf__CaMKII_CaM_Ca3__Ca; /* [13, 12] */
/* column 14 (df/dy_13) */
	jac_[82] = kr__CaM_Ca4__CaMKII; /* [3, 13] */
	jac_[289] = kr__CaMKII_CaM_Ca3__Ca; /* [12, 13] */
	jac_[312] = (-kautMax*pairedCaMKIIc)-kr__CaM_Ca4__CaMKII-kr__CaMKII_CaM_Ca3__Ca; /* [13, 13] */
	jac_[335] = kautMax*pairedCaMKIIc; /* [14, 13] */
/* column 15 (df/dy_14) */
	jac_[83] = kr__CaM_Ca4__pCaMKIIa; /* [3, 14] */
	jac_[336] = (-kr__pCaMKII_CaM_Ca3__Ca)-kr__CaM_Ca4__pCaMKIIa; /* [14, 14] */
	jac_[359] = kr__CaM_Ca4__pCaMKIIa; /* [15, 14] */
	jac_[382] = kr__pCaMKII_CaM_Ca3__Ca; /* [16, 14] */
/* column 16 (df/dy_15) */
	jac_[15] = -CaM_Ca1*kf__CaM_Ca1__pCaMKIIa; /* [0, 15] */
	jac_[38] = -CaM_Ca2*kf__CaM_Ca2__pCaMKIIa; /* [1, 15] */
	jac_[61] = -CaM_Ca3*kf__CaM_Ca3__pCaMKIIa; /* [2, 15] */
	jac_[84] = -CaM_Ca4*kf__CaM_Ca4__pCaMKIIa; /* [3, 15] */
	jac_[337] = CaM_Ca4*kf__CaM_Ca4__pCaMKIIa; /* [14, 15] */
	jac_[360] = (-PP1*kf__PP1__pCaMKIIa)-CaM*kf__CaM__pCaMKIIa-CaM_Ca4*kf__CaM_Ca4__pCaMKIIa-CaM_Ca3*kf__CaM_Ca3__pCaMKIIa-CaM_Ca2*kf__CaM_Ca2__pCaMKIIa-CaM_Ca1*kf__CaM_Ca1__pCaMKIIa; /* [15, 15] */
	jac_[383] = CaM_Ca3*kf__CaM_Ca3__pCaMKIIa; /* [16, 15] */
	jac_[406] = CaM_Ca2*kf__CaM_Ca2__pCaMKIIa; /* [17, 15] */
	jac_[429] = CaM_Ca1*kf__CaM_Ca1__pCaMKIIa; /* [18, 15] */
	jac_[452] = CaM*kf__CaM__pCaMKIIa; /* [19, 15] */
	jac_[475] = PP1*kf__PP1__pCaMKIIa; /* [20, 15] */
/* column 17 (df/dy_16) */
	jac_[62] = kr__CaM_Ca3__pCaMKIIa; /* [2, 16] */
	jac_[338] = Ca*kf__pCaMKII_CaM_Ca3__Ca; /* [14, 16] */
	jac_[361] = kr__CaM_Ca3__pCaMKIIa; /* [15, 16] */
	jac_[384] = (-kr__pCaMKII_CaM_Ca2__Ca)-kr__CaM_Ca3__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca3__Ca; /* [16, 16] */
	jac_[407] = kr__pCaMKII_CaM_Ca2__Ca; /* [17, 16] */
/* column 18 (df/dy_17) */
	jac_[40] = kr__CaM_Ca2__pCaMKIIa; /* [1, 17] */
	jac_[362] = kr__CaM_Ca2__pCaMKIIa; /* [15, 17] */
	jac_[385] = Ca*kf__pCaMKII_CaM_Ca2__Ca; /* [16, 17] */
	jac_[408] = (-kr__pCaMKII_CaM_Ca1__Ca)-kr__CaM_Ca2__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca2__Ca; /* [17, 17] */
	jac_[431] = kr__pCaMKII_CaM_Ca1__Ca; /* [18, 17] */
/* column 19 (df/dy_18) */
	jac_[18] = kr__CaM_Ca1__pCaMKIIa; /* [0, 18] */
	jac_[363] = kr__CaM_Ca1__pCaMKIIa; /* [15, 18] */
	jac_[409] = Ca*kf__pCaMKII_CaM_Ca1__Ca; /* [17, 18] */
	jac_[432] = (-kr__pCaMKII_CaM__Ca)-kr__CaM_Ca1__pCaMKIIa-Ca*kf__pCaMKII_CaM_Ca1__Ca; /* [18, 18] */
	jac_[455] = kr__pCaMKII_CaM__Ca; /* [19, 18] */
/* column 20 (df/dy_19) */
	jac_[364] = kr__CaM__pCaMKIIa; /* [15, 19] */
	jac_[433] = Ca*kf__pCaMKII_CaM__Ca; /* [18, 19] */
	jac_[456] = (-kr__CaM__pCaMKIIa)-Ca*kf__pCaMKII_CaM__Ca; /* [19, 19] */
/* column 21 (df/dy_20) */
	jac_[365] = kr__PP1__pCaMKIIa; /* [15, 20] */
	jac_[480] = (-kr__PP1__pCaMKIIa)-kcat__PP1__pCaMKIIa; /* [20, 20] */
/* column 22 (df/dy_21) */
	jac_[527] = -kca1*kca2; /* [22, 21] */
/* column 23 (df/dy_22) */
	jac_[505] = 1; /* [21, 22] */
	jac_[528] = (-kca2)-kca1; /* [22, 22] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int CaMKIIs_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 23*62;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	memset(jacp_,0,sizeof(double)*numStateVar*numParam); /* 1426 */
/* column 1 (df/dp_0) */
	jacp_[0] = Ca*CaM; /* [0, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = -Ca*CaM_Ca1; /* [0, 1] */
	jacp_[63] = Ca*CaM_Ca1; /* [1, 1] */
/* column 3 (df/dp_2) */
	jacp_[64] = -Ca*CaM_Ca2; /* [1, 2] */
	jacp_[126] = Ca*CaM_Ca2; /* [2, 2] */
/* column 4 (df/dp_3) */
	jacp_[127] = -Ca*CaM_Ca3; /* [2, 3] */
	jacp_[189] = Ca*CaM_Ca3; /* [3, 3] */
/* column 5 (df/dp_4) */
	jacp_[252] = CaM*PP2B; /* [4, 4] */
/* column 6 (df/dp_5) */
	jacp_[5] = -CaM_Ca1*PP2B; /* [0, 5] */
	jacp_[315] = CaM_Ca1*PP2B; /* [5, 5] */
/* column 7 (df/dp_6) */
	jacp_[68] = -CaM_Ca2*PP2B; /* [1, 6] */
	jacp_[378] = CaM_Ca2*PP2B; /* [6, 6] */
/* column 8 (df/dp_7) */
	jacp_[131] = -CaM_Ca3*PP2B; /* [2, 7] */
	jacp_[441] = CaM_Ca3*PP2B; /* [7, 7] */
/* column 9 (df/dp_8) */
	jacp_[194] = -CaM_Ca4*PP2B; /* [3, 8] */
	jacp_[504] = CaM_Ca4*PP2B; /* [8, 8] */
/* column 10 (df/dp_9) */
	jacp_[257] = -Ca*PP2B_CaM; /* [4, 9] */
	jacp_[319] = Ca*PP2B_CaM; /* [5, 9] */
/* column 11 (df/dp_10) */
	jacp_[320] = -Ca*PP2B_CaM_Ca1; /* [5, 10] */
	jacp_[382] = Ca*PP2B_CaM_Ca1; /* [6, 10] */
/* column 12 (df/dp_11) */
	jacp_[383] = -Ca*PP2B_CaM_Ca2; /* [6, 11] */
	jacp_[445] = Ca*PP2B_CaM_Ca2; /* [7, 11] */
/* column 13 (df/dp_12) */
	jacp_[446] = -Ca*PP2B_CaM_Ca3; /* [7, 12] */
	jacp_[508] = Ca*PP2B_CaM_Ca3; /* [8, 12] */
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
	jacp_[580] = CaM*CaMKII; /* [9, 22] */
/* column 24 (df/dp_23) */
	jacp_[767] = -Ca*CaMKII_CaM_Ca3; /* [12, 23] */
	jacp_[829] = Ca*CaMKII_CaM_Ca3; /* [13, 23] */
/* column 25 (df/dp_24) */
	jacp_[706] = -Ca*CaMKII_CaM_Ca2; /* [11, 24] */
	jacp_[768] = Ca*CaMKII_CaM_Ca2; /* [12, 24] */
/* column 26 (df/dp_25) */
	jacp_[645] = -Ca*CaMKII_CaM_Ca1; /* [10, 25] */
	jacp_[707] = Ca*CaMKII_CaM_Ca1; /* [11, 25] */
/* column 27 (df/dp_26) */
	jacp_[584] = -Ca*CaMKII_CaM; /* [9, 26] */
	jacp_[646] = Ca*CaMKII_CaM; /* [10, 26] */
/* column 28 (df/dp_27) */
	jacp_[27] = -CaMKII*CaM_Ca1; /* [0, 27] */
	jacp_[647] = CaMKII*CaM_Ca1; /* [10, 27] */
/* column 29 (df/dp_28) */
	jacp_[90] = -CaMKII*CaM_Ca2; /* [1, 28] */
	jacp_[710] = CaMKII*CaM_Ca2; /* [11, 28] */
/* column 30 (df/dp_29) */
	jacp_[153] = -CaMKII*CaM_Ca3; /* [2, 29] */
	jacp_[773] = CaMKII*CaM_Ca3; /* [12, 29] */
/* column 31 (df/dp_30) */
	jacp_[216] = -CaMKII*CaM_Ca4; /* [3, 30] */
	jacp_[836] = CaMKII*CaM_Ca4; /* [13, 30] */
/* column 32 (df/dp_31) */
/* column 33 (df/dp_32) */
/* column 34 (df/dp_33) */
/* column 35 (df/dp_34) */
/* column 36 (df/dp_35) */
/* column 37 (df/dp_36) */
	jacp_[904] = Ca*pCaMKII_CaM_Ca3; /* [14, 36] */
	jacp_[1028] = -Ca*pCaMKII_CaM_Ca3; /* [16, 36] */
/* column 38 (df/dp_37) */
	jacp_[967] = -CaM*pCaMKIIa; /* [15, 37] */
	jacp_[1215] = CaM*pCaMKIIa; /* [19, 37] */
/* column 39 (df/dp_38) */
	jacp_[38] = -CaM_Ca1*pCaMKIIa; /* [0, 38] */
	jacp_[968] = -CaM_Ca1*pCaMKIIa; /* [15, 38] */
	jacp_[1154] = CaM_Ca1*pCaMKIIa; /* [18, 38] */
/* column 40 (df/dp_39) */
	jacp_[101] = -CaM_Ca2*pCaMKIIa; /* [1, 39] */
	jacp_[969] = -CaM_Ca2*pCaMKIIa; /* [15, 39] */
	jacp_[1093] = CaM_Ca2*pCaMKIIa; /* [17, 39] */
/* column 41 (df/dp_40) */
	jacp_[164] = -CaM_Ca3*pCaMKIIa; /* [2, 40] */
	jacp_[970] = -CaM_Ca3*pCaMKIIa; /* [15, 40] */
	jacp_[1032] = CaM_Ca3*pCaMKIIa; /* [16, 40] */
/* column 42 (df/dp_41) */
	jacp_[1033] = Ca*pCaMKII_CaM_Ca2; /* [16, 41] */
	jacp_[1095] = -Ca*pCaMKII_CaM_Ca2; /* [17, 41] */
/* column 43 (df/dp_42) */
	jacp_[1096] = Ca*pCaMKII_CaM_Ca1; /* [17, 42] */
	jacp_[1158] = -Ca*pCaMKII_CaM_Ca1; /* [18, 42] */
/* column 44 (df/dp_43) */
	jacp_[229] = -CaM_Ca4*pCaMKIIa; /* [3, 43] */
	jacp_[911] = CaM_Ca4*pCaMKIIa; /* [14, 43] */
	jacp_[973] = -CaM_Ca4*pCaMKIIa; /* [15, 43] */
/* column 45 (df/dp_44) */
	jacp_[1160] = Ca*pCaMKII_CaM; /* [18, 44] */
	jacp_[1222] = -Ca*pCaMKII_CaM; /* [19, 44] */
/* column 46 (df/dp_45) */
/* column 47 (df/dp_46) */
/* column 48 (df/dp_47) */
/* column 49 (df/dp_48) */
/* column 50 (df/dp_49) */
/* column 51 (df/dp_50) */
	jacp_[856] = -CaMKII_CaM_Ca4*pairedCaMKIIc; /* [13, 50] */
	jacp_[918] = CaMKII_CaM_Ca4*pairedCaMKIIc; /* [14, 50] */
/* column 52 (df/dp_51) */
	jacp_[981] = -PP1*pCaMKIIa; /* [15, 51] */
	jacp_[1291] = PP1*pCaMKIIa; /* [20, 51] */
/* column 53 (df/dp_52) */
	jacp_[982] = PP1__pCaMKIIa; /* [15, 52] */
	jacp_[1292] = -PP1__pCaMKIIa; /* [20, 52] */
/* column 54 (df/dp_53) */
	jacp_[1293] = -PP1__pCaMKIIa; /* [20, 53] */
/* column 55 (df/dp_54) */
/* column 56 (df/dp_55) */
/* column 57 (df/dp_56) */
/* column 58 (df/dp_57) */
/* column 59 (df/dp_58) */
/* column 60 (df/dp_59) */
	jacp_[1423] = (-caa*kca2)-cab; /* [22, 59] */
/* column 61 (df/dp_60) */
	jacp_[1424] = (-caa*kca1)-cab; /* [22, 60] */
/* column 62 (df/dp_61) */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int CaMKIIs_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 5;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	func_[_CaPerCaM] = BoundCa / (s + Total_CaM); /* CaPerCaM */
	func_[_AutoCaMKII] = 100*(Total_pCaMKII / (s + totalCaMKII)); /* AutoCaMKII */
	func_[_CaMPerPP2B] = Total_PP2B_CaM_CaX / (s + PP2B + Total_PP2B_CaM_CaX); /* CaMPerPP2B */
	func_[_ActivePP2B] = 100 * (PP2B_CaM_Ca4 / (s + PP2B + Total_PP2B_CaM_CaX)); /* ActivePP2B */
	func_[_Camonitor] = Ca; /* Camonitor */
	return GSL_SUCCESS;
}
/* Function Jacobian dF(t,y;p)/dy */
int CaMKIIs_funcJac(double t, const double y_[], double *funcJac_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJac_) return 115;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	memset(funcJac_,0,sizeof(double)*numFunc*numStateVar); /* 115 */
/* column 1 (dF/dy_0) */
/* column 2 (dF/dy_1) */
/* column 3 (dF/dy_2) */
/* column 4 (dF/dy_3) */
/* column 5 (dF/dy_4) */
	funcJac_[50] = 1/(s+PP2B_0); /* [2, 4] */
/* column 6 (dF/dy_5) */
	funcJac_[51] = 1/(s+PP2B_0); /* [2, 5] */
/* column 7 (dF/dy_6) */
	funcJac_[52] = 1/(s+PP2B_0); /* [2, 6] */
/* column 8 (dF/dy_7) */
	funcJac_[53] = 1/(s+PP2B_0); /* [2, 7] */
/* column 9 (dF/dy_8) */
	funcJac_[54] = 1/(s+PP2B_0); /* [2, 8] */
	funcJac_[77] = 100/(s+PP2B_0); /* [3, 8] */
/* column 10 (dF/dy_9) */
	funcJac_[32] = -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 9] */
/* column 11 (dF/dy_10) */
	funcJac_[33] = -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 10] */
/* column 12 (dF/dy_11) */
	funcJac_[34] = -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 11] */
/* column 13 (dF/dy_12) */
	funcJac_[35] = -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 12] */
/* column 14 (dF/dy_13) */
	funcJac_[36] = -(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 13] */
/* column 15 (dF/dy_14) */
	funcJac_[37] = 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 14] */
/* column 16 (dF/dy_15) */
	funcJac_[38] = 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 15] */
/* column 17 (dF/dy_16) */
	funcJac_[39] = 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 16] */
/* column 18 (dF/dy_17) */
	funcJac_[40] = 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 17] */
/* column 19 (dF/dy_18) */
	funcJac_[41] = 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 18] */
/* column 20 (dF/dy_19) */
	funcJac_[42] = 100/(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII)-(100*(pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM))/gsl_pow_2(s+pCaMKIIa+pCaMKII_CaM_Ca4+pCaMKII_CaM_Ca3+pCaMKII_CaM_Ca2+pCaMKII_CaM_Ca1+pCaMKII_CaM+CaMKII_CaM_Ca4+CaMKII_CaM_Ca3+CaMKII_CaM_Ca2+CaMKII_CaM_Ca1+CaMKII_CaM+CaMKII); /* [1, 19] */
/* column 21 (dF/dy_20) */
/* column 22 (dF/dy_21) */
	funcJac_[113] = 1; /* [4, 21] */
/* column 23 (dF/dy_22) */
	return GSL_SUCCESS;
}
/* Function parameter Jacobian dF(t,y;p)/dp */
int CaMKIIs_funcJacp(double t, const double y_[], double *funcJacp_, void *par)
{
	double *p_=par;
	if (!y_ || !funcJacp_) return 310;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	double CaM_Ca1=y_[0];
	double CaM_Ca2=y_[1];
	double CaM_Ca3=y_[2];
	double CaM_Ca4=y_[3];
	double PP2B_CaM=y_[4];
	double PP2B_CaM_Ca1=y_[5];
	double PP2B_CaM_Ca2=y_[6];
	double PP2B_CaM_Ca3=y_[7];
	double PP2B_CaM_Ca4=y_[8];
	double CaMKII_CaM=y_[9];
	double CaMKII_CaM_Ca1=y_[10];
	double CaMKII_CaM_Ca2=y_[11];
	double CaMKII_CaM_Ca3=y_[12];
	double CaMKII_CaM_Ca4=y_[13];
	double pCaMKII_CaM_Ca4=y_[14];
	double pCaMKIIa=y_[15];
	double pCaMKII_CaM_Ca3=y_[16];
	double pCaMKII_CaM_Ca2=y_[17];
	double pCaMKII_CaM_Ca1=y_[18];
	double pCaMKII_CaM=y_[19];
	double PP1__pCaMKIIa=y_[20];
	double caa=y_[21];
	double cab=y_[22];
	double KD__CaM_Ca3__PP2B=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__PP2B) / KD__PP2B_CaM_Ca3__Ca;
	double KD__CaM_Ca2__PP2B=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__PP2B) / KD__PP2B_CaM_Ca2__Ca;
	double KD__CaM_Ca1__PP2B=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__PP2B) / KD__PP2B_CaM_Ca1__Ca;
	double KD__CaM__PP2B=(KD__CaM__Ca * KD__CaM_Ca1__PP2B) / KD__PP2B_CaM__Ca;
	double KD__CaM_Ca3__CaMKII=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__CaMKII) / KD__CaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__CaMKII=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__CaMKII) / KD__CaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__CaMKII=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__CaMKII) / KD__CaMKII_CaM_Ca1__Ca;
	double KD__CaM__CaMKII=(KD__CaM__Ca * KD__CaM_Ca1__CaMKII) / KD__CaMKII_CaM__Ca;
	double KD__CaM_Ca3__pCaMKIIa=(KD__CaM_Ca3__Ca * KD__CaM_Ca4__pCaMKIIa) / KD__pCaMKII_CaM_Ca3__Ca;
	double KD__CaM_Ca2__pCaMKIIa=(KD__CaM_Ca2__Ca * KD__CaM_Ca3__pCaMKIIa) / KD__pCaMKII_CaM_Ca2__Ca;
	double KD__CaM_Ca1__pCaMKIIa=(KD__CaM_Ca1__Ca * KD__CaM_Ca2__pCaMKIIa) / KD__pCaMKII_CaM_Ca1__Ca;
	double KD__CaM__pCaMKIIa=(KD__CaM__Ca * KD__CaM_Ca1__pCaMKIIa) / KD__pCaMKII_CaM__Ca;
	double kr__CaM_Ca3__Ca=KD__CaM_Ca3__Ca*kf__CaM_Ca3__Ca;
	double kr__CaM_Ca2__Ca=KD__CaM_Ca2__Ca * kf__CaM_Ca2__Ca;
	double kr__CaM_Ca1__Ca=KD__CaM_Ca1__Ca * kf__CaM_Ca1__Ca;
	double kr__CaM__Ca=KD__CaM__Ca * kf__CaM__Ca;
	double kr__CaM_Ca4__PP2B=KD__CaM_Ca4__PP2B * kf__CaM_Ca4__PP2B;
	double kr__CaM_Ca3__PP2B=KD__CaM_Ca3__PP2B * kf__CaM_Ca3__PP2B;
	double kr__CaM_Ca2__PP2B=KD__CaM_Ca2__PP2B * kf__CaM_Ca2__PP2B;
	double kr__CaM_Ca1__PP2B=KD__CaM_Ca1__PP2B * kf__CaM_Ca1__PP2B;
	double kr__CaM__PP2B=KD__CaM__PP2B * kf__CaM__PP2B;
	double kr__PP2B_CaM_Ca3__Ca=KD__PP2B_CaM_Ca3__Ca *kf__PP2B_CaM_Ca3__Ca;
	double kr__PP2B_CaM_Ca2__Ca=KD__PP2B_CaM_Ca2__Ca * kf__PP2B_CaM_Ca2__Ca;
	double kr__PP2B_CaM_Ca1__Ca=KD__PP2B_CaM_Ca1__Ca * kf__PP2B_CaM_Ca1__Ca;
	double kr__PP2B_CaM__Ca=KD__PP2B_CaM__Ca * kf__PP2B_CaM__Ca;
	double kr__CaM_Ca4__CaMKII=KD__CaM_Ca4__CaMKII * kf__CaM_Ca4__CaMKII;
	double kr__CaM_Ca3__CaMKII=KD__CaM_Ca3__CaMKII * kf__CaM_Ca3__CaMKII;
	double kr__CaM_Ca2__CaMKII=KD__CaM_Ca2__CaMKII * kf__CaM_Ca2__CaMKII;
	double kr__CaM_Ca1__CaMKII=KD__CaM_Ca1__CaMKII * kf__CaM_Ca1__CaMKII;
	double kr__CaM__CaMKII=KD__CaM__CaMKII * kf__CaM__CaMKII;
	double kr__CaMKII_CaM_Ca3__Ca=KD__CaMKII_CaM_Ca3__Ca * kf__CaMKII_CaM_Ca3__Ca;
	double kr__CaMKII_CaM_Ca2__Ca=KD__CaMKII_CaM_Ca2__Ca * kf__CaMKII_CaM_Ca2__Ca;
	double kr__CaMKII_CaM_Ca1__Ca=KD__CaMKII_CaM_Ca1__Ca * kf__CaMKII_CaM_Ca1__Ca;
	double kr__CaMKII_CaM__Ca=KD__CaMKII_CaM__Ca * kf__CaMKII_CaM__Ca;
	double kr__pCaMKII_CaM_Ca3__Ca=KD__pCaMKII_CaM_Ca3__Ca * kf__pCaMKII_CaM_Ca3__Ca;
	double kr__pCaMKII_CaM_Ca2__Ca=KD__pCaMKII_CaM_Ca2__Ca * kf__pCaMKII_CaM_Ca2__Ca;
	double kr__pCaMKII_CaM_Ca1__Ca=KD__pCaMKII_CaM_Ca1__Ca * kf__pCaMKII_CaM_Ca1__Ca;
	double kr__pCaMKII_CaM__Ca=KD__pCaMKII_CaM__Ca * kf__pCaMKII_CaM__Ca;
	double kr__CaM_Ca4__pCaMKIIa=KD__CaM_Ca4__pCaMKIIa * kf__CaM_Ca4__pCaMKIIa;
	double kr__CaM_Ca3__pCaMKIIa=KD__CaM_Ca3__pCaMKIIa * kf__CaM_Ca3__pCaMKIIa;
	double kr__CaM_Ca2__pCaMKIIa=KD__CaM_Ca2__pCaMKIIa * kf__CaM_Ca2__pCaMKIIa;
	double kr__CaM_Ca1__pCaMKIIa=KD__CaM_Ca1__pCaMKIIa * kf__CaM_Ca1__pCaMKIIa;
	double kr__CaM__pCaMKIIa=KD__CaM__pCaMKIIa * kf__CaM__pCaMKIIa;
	double Ca=isOn*Ca_set+caa;
	double PP1=PP1_0-PP1__pCaMKIIa;
	double CaMKII=CaMKII_0 - PP1_0 - (CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM - PP1);
	double CaM=CaM_0 + PP1_0 - CaMKII_0 - (CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4 + PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4 + PP1 - CaMKII - pCaMKIIa);
	double PP2B=PP2B_0-(PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4);
	double totalCaMKII=(CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM);
	double ActiveCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM + CaMKII_CaM_Ca4;
	double r=(CaMKII_CaM_Ca4+pCaMKII_CaM_Ca4)/(s+totalCaMKII);
	double pairedCaMKIIc=a*(r*r)/(1+b*r);
	double CaM_bound_Ca=(1.0 * CaM_Ca1) + (2.0 * CaM_Ca2) + (3.0 * CaM_Ca3) + (4.0 * CaM_Ca4);
	double PP2B_bound_Ca=(1.0 * PP2B_CaM_Ca1) + (2.0 * PP2B_CaM_Ca2) + (3.0 * PP2B_CaM_Ca3) + (4.0 * PP2B_CaM_Ca4);
	double CaMKII_bound_Ca=(1.0 * (CaMKII_CaM_Ca1 + pCaMKII_CaM_Ca1)) + (2.0 * (CaMKII_CaM_Ca2 + pCaMKII_CaM_Ca2)) + (3.0 * (CaMKII_CaM_Ca3 + pCaMKII_CaM_Ca3)) + (4.0 * (CaMKII_CaM_Ca4 + pCaMKII_CaM_Ca4));
	double BoundCa=(CaM_bound_Ca + PP2B_bound_Ca + CaMKII_bound_Ca);
	double Total_CaM_CaX=CaM + CaM_Ca1 + CaM_Ca2 + CaM_Ca3 + CaM_Ca4;
	double Total_PP2B_CaM_CaX=PP2B_CaM + PP2B_CaM_Ca1 + PP2B_CaM_Ca2 + PP2B_CaM_Ca3 + PP2B_CaM_Ca4;
	double Total_CaMKII_CaM_CaX=CaMKII_CaM + CaMKII_CaM_Ca1 + CaMKII_CaM_Ca2 + CaMKII_CaM_Ca3 + CaMKII_CaM_Ca4;
	double Total_pCaMKII_CaM_CaX=pCaMKII_CaM + pCaMKII_CaM_Ca1 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca4;
	double Total_CaM=Total_CaM_CaX + Total_PP2B_CaM_CaX + Total_CaMKII_CaM_CaX + Total_pCaMKII_CaM_CaX;
	double Total_pCaMKII=pCaMKII_CaM_Ca4 + pCaMKIIa + pCaMKII_CaM_Ca3 + pCaMKII_CaM_Ca2 + pCaMKII_CaM_Ca1 + pCaMKII_CaM;
	double ReactionFlux1=(kf__CaM__Ca*Ca*CaM)-(kr__CaM__Ca*CaM_Ca1);
	double ReactionFlux2=(kf__CaM_Ca1__Ca*Ca*CaM_Ca1)-(kr__CaM_Ca1__Ca*CaM_Ca2);
	double ReactionFlux3=(kf__CaM_Ca2__Ca*Ca*CaM_Ca2)-(kr__CaM_Ca2__Ca*CaM_Ca3);
	double ReactionFlux4=(kf__CaM_Ca3__Ca*CaM_Ca3*Ca)-(kr__CaM_Ca3__Ca*CaM_Ca4);
	double ReactionFlux5=(kf__CaM__PP2B*CaM*PP2B)-(kr__CaM__PP2B*PP2B_CaM);
	double ReactionFlux6=(kf__CaM_Ca1__PP2B*CaM_Ca1*PP2B)-(kr__CaM_Ca1__PP2B*PP2B_CaM_Ca1);
	double ReactionFlux7=(kf__CaM_Ca2__PP2B*CaM_Ca2*PP2B)-(kr__CaM_Ca2__PP2B*PP2B_CaM_Ca2);
	double ReactionFlux8=(kf__CaM_Ca3__PP2B*CaM_Ca3*PP2B)-(kr__CaM_Ca3__PP2B*PP2B_CaM_Ca3);
	double ReactionFlux9=(kf__CaM_Ca4__PP2B*CaM_Ca4*PP2B)-(kr__CaM_Ca4__PP2B*PP2B_CaM_Ca4);
	double ReactionFlux10=(kf__PP2B_CaM__Ca*PP2B_CaM*Ca)-(kr__PP2B_CaM__Ca*PP2B_CaM_Ca1);
	double ReactionFlux11=(kf__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca1*Ca)-(kr__PP2B_CaM_Ca1__Ca*PP2B_CaM_Ca2);
	double ReactionFlux12=(kf__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca2*Ca)-(kr__PP2B_CaM_Ca2__Ca*PP2B_CaM_Ca3);
	double ReactionFlux13=(kf__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca3*Ca)-(kr__PP2B_CaM_Ca3__Ca*PP2B_CaM_Ca4);
	double ReactionFlux14=(kf__CaM__CaMKII*CaM*CaMKII)-(kr__CaM__CaMKII*CaMKII_CaM);
	double ReactionFlux15=(kf__CaM_Ca1__CaMKII*CaM_Ca1*CaMKII)-(kr__CaM_Ca1__CaMKII*CaMKII_CaM_Ca1);
	double ReactionFlux16=(kf__CaM_Ca2__CaMKII*CaM_Ca2*CaMKII)-(kr__CaM_Ca2__CaMKII*CaMKII_CaM_Ca2);
	double ReactionFlux17=(kf__CaM_Ca3__CaMKII*CaM_Ca3*CaMKII)-(kr__CaM_Ca3__CaMKII*CaMKII_CaM_Ca3);
	double ReactionFlux18=(kf__CaM_Ca4__CaMKII*CaM_Ca4*CaMKII)-(kr__CaM_Ca4__CaMKII*CaMKII_CaM_Ca4);
	double ReactionFlux19=(kf__CaMKII_CaM__Ca*Ca*CaMKII_CaM)-(kr__CaMKII_CaM__Ca*CaMKII_CaM_Ca1);
	double ReactionFlux20=(kf__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca1*Ca)-(kr__CaMKII_CaM_Ca1__Ca*CaMKII_CaM_Ca2);
	double ReactionFlux21=(kf__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca2*Ca)-(kr__CaMKII_CaM_Ca2__Ca*CaMKII_CaM_Ca3);
	double ReactionFlux22=(kf__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca3*Ca)-(kr__CaMKII_CaM_Ca3__Ca*CaMKII_CaM_Ca4);
	double ReactionFlux23=(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4*pCaMKIIa)-(kr__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4);
	double ReactionFlux24=(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3*pCaMKIIa)-(kr__CaM_Ca3__pCaMKIIa*pCaMKII_CaM_Ca3);
	double ReactionFlux25=(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2*pCaMKIIa)-(kr__CaM_Ca2__pCaMKIIa*pCaMKII_CaM_Ca2);
	double ReactionFlux26=(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1*pCaMKIIa)-(kr__CaM_Ca1__pCaMKIIa*pCaMKII_CaM_Ca1);
	double ReactionFlux27=(kf__CaM__pCaMKIIa*CaM*pCaMKIIa)-(kr__CaM__pCaMKIIa*pCaMKII_CaM);
	double ReactionFlux28=(kf__pCaMKII_CaM__Ca*pCaMKII_CaM*Ca)-(kr__pCaMKII_CaM__Ca*pCaMKII_CaM_Ca1);
	double ReactionFlux29=(kf__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca1*Ca)-(kr__pCaMKII_CaM_Ca1__Ca*pCaMKII_CaM_Ca2);
	double ReactionFlux30=(kf__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca2*Ca)-(kr__pCaMKII_CaM_Ca2__Ca*pCaMKII_CaM_Ca3);
	double ReactionFlux31=(kf__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca3*Ca)-(kr__pCaMKII_CaM_Ca3__Ca*pCaMKII_CaM_Ca4);
	double ReactionFlux32=(kautMax*pairedCaMKIIc*CaMKII_CaM_Ca4);
	double ReactionFlux33=(kf__PP1__pCaMKIIa*pCaMKIIa*PP1)-(kr__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double ReactionFlux34=(kcat__PP1__pCaMKIIa*PP1__pCaMKIIa);
	double SpikeFlux1=cab;
	double SpikeFlux2=0-(kca1*kca2*caa+(kca1+kca2)*cab);
	memset(funcJacp_,0,sizeof(double)*numFunc*numParam); /* 310 */
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
/* column 29 (dF/dp_28) */
/* column 30 (dF/dp_29) */
/* column 31 (dF/dp_30) */
/* column 32 (dF/dp_31) */
/* column 33 (dF/dp_32) */
/* column 34 (dF/dp_33) */
/* column 35 (dF/dp_34) */
/* column 36 (dF/dp_35) */
/* column 37 (dF/dp_36) */
/* column 38 (dF/dp_37) */
/* column 39 (dF/dp_38) */
/* column 40 (dF/dp_39) */
/* column 41 (dF/dp_40) */
/* column 42 (dF/dp_41) */
/* column 43 (dF/dp_42) */
/* column 44 (dF/dp_43) */
/* column 45 (dF/dp_44) */
/* column 46 (dF/dp_45) */
/* column 47 (dF/dp_46) */
/* column 48 (dF/dp_47) */
/* column 49 (dF/dp_48) */
/* column 50 (dF/dp_49) */
/* column 51 (dF/dp_50) */
/* column 52 (dF/dp_51) */
/* column 53 (dF/dp_52) */
/* column 54 (dF/dp_53) */
/* column 55 (dF/dp_54) */
	funcJacp_[302] = isOn; /* [4, 54] */
/* column 56 (dF/dp_55) */
/* column 57 (dF/dp_56) */
/* column 58 (dF/dp_57) */
/* column 59 (dF/dp_58) */
	funcJacp_[182] = -(PP2B_CaM_Ca4+PP2B_CaM_Ca3+PP2B_CaM_Ca2+PP2B_CaM_Ca1+PP2B_CaM)/gsl_pow_2(s+PP2B_0); /* [2, 58] */
	funcJacp_[244] = -(100*PP2B_CaM_Ca4)/gsl_pow_2(s+PP2B_0); /* [3, 58] */
/* column 60 (dF/dp_59) */
/* column 61 (dF/dp_60) */
/* column 62 (dF/dp_61) */
	funcJacp_[309] = Ca_set; /* [4, 61] */
	return GSL_SUCCESS;
}
/* ode default parameters */
int CaMKIIs_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 62;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	memset(p_,0,sizeof(double)*numParam);
	p_[_kf__CaM__Ca] = 0.151909999999999;
	p_[_kf__CaM_Ca1__Ca] = 3.42450000000002e-05;
	p_[_kf__CaM_Ca2__Ca] = 0.0738930000000004;
	p_[_kf__CaM_Ca3__Ca] = 0.0061814;
	p_[_kf__CaM__PP2B] = 0.0202929999999999;
	p_[_kf__CaM_Ca1__PP2B] = 0.00452479999999998;
	p_[_kf__CaM_Ca2__PP2B] = 0.051176;
	p_[_kf__CaM_Ca3__PP2B] = 0.274210000000001;
	p_[_kf__CaM_Ca4__PP2B] = 0.0833569999999997;
	p_[_kf__PP2B_CaM__Ca] = 0.0011578;
	p_[_kf__PP2B_CaM_Ca1__Ca] = 0.0047884;
	p_[_kf__PP2B_CaM_Ca2__Ca] = 0.0350789999999999;
	p_[_kf__PP2B_CaM_Ca3__Ca] = 0.0455659999999999;
	p_[_KD__CaM_Ca3__Ca] = 7271.30000000002;
	p_[_KD__CaM_Ca2__Ca] = 37062.0000000017;
	p_[_KD__CaM_Ca1__Ca] = 1827.9;
	p_[_KD__CaM__Ca] = 2662.30000000001;
	p_[_KD__CaM_Ca4__PP2B] = 0.0399700000000002;
	p_[_KD__PP2B_CaM_Ca3__Ca] = 91.5429999999998;
	p_[_KD__PP2B_CaM_Ca2__Ca] = 916.149999999996;
	p_[_KD__PP2B_CaM_Ca1__Ca] = 285.030000000001;
	p_[_KD__PP2B_CaM__Ca] = 82.8369999999998;
	p_[_kf__CaM__CaMKII] = 0.237449999999999;
	p_[_kf__CaMKII_CaM_Ca3__Ca] = 0.0258579999999999;
	p_[_kf__CaMKII_CaM_Ca2__Ca] = 0.130860000000001;
	p_[_kf__CaMKII_CaM_Ca1__Ca] = 0.0755390000000001;
	p_[_kf__CaMKII_CaM__Ca] = 0.000797720000000003;
	p_[_kf__CaM_Ca1__CaMKII] = 0.0558820000000003;
	p_[_kf__CaM_Ca2__CaMKII] = 0.0460280000000002;
	p_[_kf__CaM_Ca3__CaMKII] = 0.208550000000001;
	p_[_kf__CaM_Ca4__CaMKII] = 0.0226620000000001;
	p_[_KD__CaM_Ca4__CaMKII] = 8.28490000000001;
	p_[_KD__CaMKII_CaM_Ca3__Ca] = 483.480000000002;
	p_[_KD__CaMKII_CaM_Ca2__Ca] = 1143.6;
	p_[_KD__CaMKII_CaM_Ca1__Ca] = 645.069999999998;
	p_[_KD__CaMKII_CaM__Ca] = 3081.60000000001;
	p_[_kf__pCaMKII_CaM_Ca3__Ca] = 0.000829840000000001;
	p_[_kf__CaM__pCaMKIIa] = 0.00032583;
	p_[_kf__CaM_Ca1__pCaMKIIa] = 0.0589279999999999;
	p_[_kf__CaM_Ca2__pCaMKIIa] = 0.0231899999999999;
	p_[_kf__CaM_Ca3__pCaMKIIa] = 0.0302520000000001;
	p_[_kf__pCaMKII_CaM_Ca2__Ca] = 0.0384980000000001;
	p_[_kf__pCaMKII_CaM_Ca1__Ca] = 0.0004565;
	p_[_kf__CaM_Ca4__pCaMKIIa] = 0.0722370000000001;
	p_[_kf__pCaMKII_CaM__Ca] = 0.00212670000000001;
	p_[_KD__pCaMKII_CaM_Ca3__Ca] = 539.41;
	p_[_KD__pCaMKII_CaM_Ca2__Ca] = 1784.00000000001;
	p_[_KD__pCaMKII_CaM_Ca1__Ca] = 57728.0000000003;
	p_[_KD__pCaMKII_CaM__Ca] = 1342.9;
	p_[_KD__CaM_Ca4__pCaMKIIa] = 3.74600000000001;
	p_[_kautMax] = 0.00375589999999998;
	p_[_kf__PP1__pCaMKIIa] = 0.0016604;
	p_[_kr__PP1__pCaMKIIa] = 0.205170000000001;
	p_[_kcat__PP1__pCaMKIIa] = 0.302249999999999;
	p_[_Ca_set] = 2187.8;
	p_[_CaM_0] = 30;
	p_[_PP2B_0] = 3;
	p_[_kca1] = 0.00978422;
	p_[_kca2] = 0.03448;
	return GSL_SUCCESS;
}
/* ode initial values */
int CaMKIIs_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 23;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	double kf__CaM__Ca=p_[0];
	double kf__CaM_Ca1__Ca=p_[1];
	double kf__CaM_Ca2__Ca=p_[2];
	double kf__CaM_Ca3__Ca=p_[3];
	double kf__CaM__PP2B=p_[4];
	double kf__CaM_Ca1__PP2B=p_[5];
	double kf__CaM_Ca2__PP2B=p_[6];
	double kf__CaM_Ca3__PP2B=p_[7];
	double kf__CaM_Ca4__PP2B=p_[8];
	double kf__PP2B_CaM__Ca=p_[9];
	double kf__PP2B_CaM_Ca1__Ca=p_[10];
	double kf__PP2B_CaM_Ca2__Ca=p_[11];
	double kf__PP2B_CaM_Ca3__Ca=p_[12];
	double KD__CaM_Ca3__Ca=p_[13];
	double KD__CaM_Ca2__Ca=p_[14];
	double KD__CaM_Ca1__Ca=p_[15];
	double KD__CaM__Ca=p_[16];
	double KD__CaM_Ca4__PP2B=p_[17];
	double KD__PP2B_CaM_Ca3__Ca=p_[18];
	double KD__PP2B_CaM_Ca2__Ca=p_[19];
	double KD__PP2B_CaM_Ca1__Ca=p_[20];
	double KD__PP2B_CaM__Ca=p_[21];
	double kf__CaM__CaMKII=p_[22];
	double kf__CaMKII_CaM_Ca3__Ca=p_[23];
	double kf__CaMKII_CaM_Ca2__Ca=p_[24];
	double kf__CaMKII_CaM_Ca1__Ca=p_[25];
	double kf__CaMKII_CaM__Ca=p_[26];
	double kf__CaM_Ca1__CaMKII=p_[27];
	double kf__CaM_Ca2__CaMKII=p_[28];
	double kf__CaM_Ca3__CaMKII=p_[29];
	double kf__CaM_Ca4__CaMKII=p_[30];
	double KD__CaM_Ca4__CaMKII=p_[31];
	double KD__CaMKII_CaM_Ca3__Ca=p_[32];
	double KD__CaMKII_CaM_Ca2__Ca=p_[33];
	double KD__CaMKII_CaM_Ca1__Ca=p_[34];
	double KD__CaMKII_CaM__Ca=p_[35];
	double kf__pCaMKII_CaM_Ca3__Ca=p_[36];
	double kf__CaM__pCaMKIIa=p_[37];
	double kf__CaM_Ca1__pCaMKIIa=p_[38];
	double kf__CaM_Ca2__pCaMKIIa=p_[39];
	double kf__CaM_Ca3__pCaMKIIa=p_[40];
	double kf__pCaMKII_CaM_Ca2__Ca=p_[41];
	double kf__pCaMKII_CaM_Ca1__Ca=p_[42];
	double kf__CaM_Ca4__pCaMKIIa=p_[43];
	double kf__pCaMKII_CaM__Ca=p_[44];
	double KD__pCaMKII_CaM_Ca3__Ca=p_[45];
	double KD__pCaMKII_CaM_Ca2__Ca=p_[46];
	double KD__pCaMKII_CaM_Ca1__Ca=p_[47];
	double KD__pCaMKII_CaM__Ca=p_[48];
	double KD__CaM_Ca4__pCaMKIIa=p_[49];
	double kautMax=p_[50];
	double kf__PP1__pCaMKIIa=p_[51];
	double kr__PP1__pCaMKIIa=p_[52];
	double kcat__PP1__pCaMKIIa=p_[53];
	double Ca_set=p_[54];
	double PP1_0=p_[55];
	double CaMKII_0=p_[56];
	double CaM_0=p_[57];
	double PP2B_0=p_[58];
	double kca1=p_[59];
	double kca2=p_[60];
	double isOn=p_[61];
	/* the initial value of y may depend on the parameters. */
	memset(y_,0,sizeof(double)*numStateVar);
	return GSL_SUCCESS;
}
