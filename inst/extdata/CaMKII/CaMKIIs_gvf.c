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
int CaMKIIs_vf(double t, const double y_[], double f_[], void *par)
{
	double *p_=par;
	if (!y_ || !f_) return 21;
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
	double logistic=1.0/(1+exp(-10*t));
	double Ca=logistic*Ca_set;
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
	f_[0] = +ReactionFlux1-ReactionFlux2-ReactionFlux6-ReactionFlux15-ReactionFlux26;
	f_[1] = +ReactionFlux2-ReactionFlux3-ReactionFlux7-ReactionFlux16-ReactionFlux25;
	f_[2] = +ReactionFlux3-ReactionFlux4-ReactionFlux8-ReactionFlux17-ReactionFlux24;
	f_[3] = +ReactionFlux4-ReactionFlux9-ReactionFlux18-ReactionFlux23;
	f_[4] = +ReactionFlux5-ReactionFlux10;
	f_[5] = +ReactionFlux6+ReactionFlux10-ReactionFlux11;
	f_[6] = +ReactionFlux7+ReactionFlux11-ReactionFlux12;
	f_[7] = +ReactionFlux8+ReactionFlux12-ReactionFlux13;
	f_[8] = +ReactionFlux9+ReactionFlux13;
	f_[9] = +ReactionFlux14-ReactionFlux19;
	f_[10] = +ReactionFlux15+ReactionFlux19-ReactionFlux20;
	f_[11] = +ReactionFlux16+ReactionFlux20-ReactionFlux21;
	f_[12] = +ReactionFlux17+ReactionFlux21-ReactionFlux22;
	f_[13] = +ReactionFlux18+ReactionFlux22-ReactionFlux32;
	f_[14] = +ReactionFlux23+ReactionFlux31+ReactionFlux32;
	f_[15] = -ReactionFlux23-ReactionFlux24-ReactionFlux25-ReactionFlux26-ReactionFlux27-ReactionFlux33;
	f_[16] = +ReactionFlux24+ReactionFlux30-ReactionFlux31;
	f_[17] = +ReactionFlux25+ReactionFlux29-ReactionFlux30;
	f_[18] = +ReactionFlux26+ReactionFlux28-ReactionFlux29;
	f_[19] = +ReactionFlux27-ReactionFlux28;
	f_[20] = +ReactionFlux33-ReactionFlux34;
	return GSL_SUCCESS;
}
/* ode Jacobian df(t,y;p)/dy */
int CaMKIIs_jac(double t, const double y_[], double *jac_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jac_) return 21*21;
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
	double logistic=1.0/(1+exp(-10*t));
	double Ca=logistic*Ca_set;
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
/* column 1 (df/dy_0) */
	jac_[0] = ((0-(kf__CaM_Ca1__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 0] */
	jac_[21] = 0; /* [1, 0] */
/* column 2 (df/dy_1) */
	jac_[1] = ((0-(kf__CaM_Ca2__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 1] */
	jac_[22] = 0; /* [1, 1] */
/* column 3 (df/dy_2) */
	jac_[2] = ((0-(kf__CaM_Ca3__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 2] */
	jac_[23] = 0; /* [1, 2] */
/* column 4 (df/dy_3) */
	jac_[3] = ((-1*(kf__CaM_Ca4__pCaMKIIa*pCaMKIIa))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 3] */
	jac_[24] = 0; /* [1, 3] */
/* column 5 (df/dy_4) */
	jac_[4] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 4] */
	jac_[25] = 0; /* [1, 4] */
/* column 6 (df/dy_5) */
	jac_[5] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 5] */
	jac_[26] = 0; /* [1, 5] */
/* column 7 (df/dy_6) */
	jac_[6] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 6] */
	jac_[27] = 0; /* [1, 6] */
/* column 8 (df/dy_7) */
	jac_[7] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 7] */
	jac_[28] = 0; /* [1, 7] */
/* column 9 (df/dy_8) */
	jac_[8] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 8] */
	jac_[29] = 0; /* [1, 8] */
/* column 10 (df/dy_9) */
	jac_[9] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 9] */
	jac_[30] = 0; /* [1, 9] */
/* column 11 (df/dy_10) */
	jac_[10] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 10] */
	jac_[31] = 0; /* [1, 10] */
/* column 12 (df/dy_11) */
	jac_[11] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 11] */
	jac_[32] = 0; /* [1, 11] */
/* column 13 (df/dy_12) */
	jac_[12] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 12] */
	jac_[33] = 0; /* [1, 12] */
/* column 14 (df/dy_13) */
	jac_[13] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 13] */
	jac_[34] = 0; /* [1, 13] */
/* column 15 (df/dy_14) */
	jac_[14] = ((-1*(0-(KD__CaM_Ca4__pCaMKIIa*kf__CaM_Ca4__pCaMKIIa)))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 14] */
	jac_[35] = 0; /* [1, 14] */
/* column 16 (df/dy_15) */
	jac_[15] = ((((((-1*(kf__CaM_Ca4__pCaMKIIa*CaM_Ca4))-(kf__CaM_Ca3__pCaMKIIa*CaM_Ca3))-(kf__CaM_Ca2__pCaMKIIa*CaM_Ca2))-(kf__CaM_Ca1__pCaMKIIa*CaM_Ca1))-(kf__CaM__pCaMKIIa*((((CaM_0+PP1_0)-CaMKII_0)-(((((((((((CaM_Ca1+CaM_Ca2)+CaM_Ca3)+CaM_Ca4)+PP2B_CaM)+PP2B_CaM_Ca1)+PP2B_CaM_Ca2)+PP2B_CaM_Ca3)+PP2B_CaM_Ca4)+(PP1_0-PP1__pCaMKIIa))-((CaMKII_0-PP1_0)-(((((((((((CaMKII_CaM+CaMKII_CaM_Ca1)+CaMKII_CaM_Ca2)+CaMKII_CaM_Ca3)+CaMKII_CaM_Ca4)+pCaMKII_CaM_Ca4)+pCaMKIIa)+pCaMKII_CaM_Ca3)+pCaMKII_CaM_Ca2)+pCaMKII_CaM_Ca1)+pCaMKII_CaM)-(PP1_0-PP1__pCaMKIIa))))-pCaMKIIa))-(0*pCaMKII_CaM))))-(kf__PP1__pCaMKIIa*(PP1_0-PP1__pCaMKIIa))); /* [0, 15] */
	jac_[36] = (kf__PP1__pCaMKIIa*(PP1_0-PP1__pCaMKIIa)); /* [1, 15] */
/* column 17 (df/dy_16) */
	jac_[16] = ((0-(0-(kf__CaM_Ca3__pCaMKIIa*((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca))))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 16] */
	jac_[37] = 0; /* [1, 16] */
/* column 18 (df/dy_17) */
	jac_[17] = ((0-(0-(KD__CaM_Ca2__Ca*(kf__CaM_Ca2__pCaMKIIa*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)))))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 17] */
	jac_[38] = 0; /* [1, 17] */
/* column 19 (df/dy_18) */
	jac_[18] = ((0-(0-(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(kf__CaM_Ca1__pCaMKIIa*((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca))))))-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(0*pCaMKII_CaM)))); /* [0, 18] */
	jac_[39] = 0; /* [1, 18] */
/* column 20 (df/dy_19) */
	jac_[19] = (0-(kf__CaM__pCaMKIIa*((-1*pCaMKIIa)-(KD__CaM__Ca*(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)/KD__pCaMKII_CaM__Ca))))))); /* [0, 19] */
	jac_[40] = 0; /* [1, 19] */
/* column 21 (df/dy_20) */
	jac_[20] = (0-(((kf__PP1__pCaMKIIa*pCaMKIIa)*-1)-kr__PP1__pCaMKIIa)); /* [0, 20] */
	jac_[41] = ((((kf__PP1__pCaMKIIa*pCaMKIIa)*-1)-kr__PP1__pCaMKIIa)-kcat__PP1__pCaMKIIa); /* [1, 20] */
	return GSL_SUCCESS;
}
/* ode parameter Jacobian df(t,y;p)/dp */
int CaMKIIs_jacp(double t, const double y_[], double *jacp_, double *dfdt_, void *par)
{
	double *p_=par;
	if (!y_ || !jacp_) return 21*59;
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
	double logistic=1.0/(1+exp(-10*t));
	double Ca=logistic*Ca_set;
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
/* column 1 (df/dp_0) */
	jacp_[0] = 0; /* [0, 0] */
	jacp_[59] = 0; /* [1, 0] */
/* column 2 (df/dp_1) */
	jacp_[1] = 0; /* [0, 1] */
	jacp_[60] = 0; /* [1, 1] */
/* column 3 (df/dp_2) */
	jacp_[2] = 0; /* [0, 2] */
	jacp_[61] = 0; /* [1, 2] */
/* column 4 (df/dp_3) */
	jacp_[3] = 0; /* [0, 3] */
	jacp_[62] = 0; /* [1, 3] */
/* column 5 (df/dp_4) */
	jacp_[4] = 0; /* [0, 4] */
	jacp_[63] = 0; /* [1, 4] */
/* column 6 (df/dp_5) */
	jacp_[5] = 0; /* [0, 5] */
	jacp_[64] = 0; /* [1, 5] */
/* column 7 (df/dp_6) */
	jacp_[6] = 0; /* [0, 6] */
	jacp_[65] = 0; /* [1, 6] */
/* column 8 (df/dp_7) */
	jacp_[7] = 0; /* [0, 7] */
	jacp_[66] = 0; /* [1, 7] */
/* column 9 (df/dp_8) */
	jacp_[8] = 0; /* [0, 8] */
	jacp_[67] = 0; /* [1, 8] */
/* column 10 (df/dp_9) */
	jacp_[9] = 0; /* [0, 9] */
	jacp_[68] = 0; /* [1, 9] */
/* column 11 (df/dp_10) */
	jacp_[10] = 0; /* [0, 10] */
	jacp_[69] = 0; /* [1, 10] */
/* column 12 (df/dp_11) */
	jacp_[11] = 0; /* [0, 11] */
	jacp_[70] = 0; /* [1, 11] */
/* column 13 (df/dp_12) */
	jacp_[12] = 0; /* [0, 12] */
	jacp_[71] = 0; /* [1, 12] */
/* column 14 (df/dp_13) */
	jacp_[13] = ((((0-(0-(((KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca)*kf__CaM_Ca3__pCaMKIIa)*pCaMKII_CaM_Ca3)))-(0-((((KD__CaM_Ca2__Ca*(KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca4__pCaMKIIa/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 13] */
	jacp_[72] = 0; /* [1, 13] */
/* column 15 (df/dp_14) */
	jacp_[14] = (((0-(0-(((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 14] */
	jacp_[73] = 0; /* [1, 14] */
/* column 16 (df/dp_15) */
	jacp_[15] = ((0-(0-((((KD__CaM_Ca2__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca2__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 15] */
	jacp_[74] = 0; /* [1, 15] */
/* column 17 (df/dp_16) */
	jacp_[16] = (0-(kf__CaM__pCaMKIIa*(0-(((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 16] */
	jacp_[75] = 0; /* [1, 16] */
/* column 18 (df/dp_17) */
	jacp_[17] = 0; /* [0, 17] */
	jacp_[76] = 0; /* [1, 17] */
/* column 19 (df/dp_18) */
	jacp_[18] = 0; /* [0, 18] */
	jacp_[77] = 0; /* [1, 18] */
/* column 20 (df/dp_19) */
	jacp_[19] = 0; /* [0, 19] */
	jacp_[78] = 0; /* [1, 19] */
/* column 21 (df/dp_20) */
	jacp_[20] = 0; /* [0, 20] */
	jacp_[79] = 0; /* [1, 20] */
/* column 22 (df/dp_21) */
	jacp_[21] = 0; /* [0, 21] */
	jacp_[80] = 0; /* [1, 21] */
/* column 23 (df/dp_22) */
	jacp_[22] = 0; /* [0, 22] */
	jacp_[81] = 0; /* [1, 22] */
/* column 24 (df/dp_23) */
	jacp_[23] = 0; /* [0, 23] */
	jacp_[82] = 0; /* [1, 23] */
/* column 25 (df/dp_24) */
	jacp_[24] = 0; /* [0, 24] */
	jacp_[83] = 0; /* [1, 24] */
/* column 26 (df/dp_25) */
	jacp_[25] = 0; /* [0, 25] */
	jacp_[84] = 0; /* [1, 25] */
/* column 27 (df/dp_26) */
	jacp_[26] = 0; /* [0, 26] */
	jacp_[85] = 0; /* [1, 26] */
/* column 28 (df/dp_27) */
	jacp_[27] = 0; /* [0, 27] */
	jacp_[86] = 0; /* [1, 27] */
/* column 29 (df/dp_28) */
	jacp_[28] = 0; /* [0, 28] */
	jacp_[87] = 0; /* [1, 28] */
/* column 30 (df/dp_29) */
	jacp_[29] = 0; /* [0, 29] */
	jacp_[88] = 0; /* [1, 29] */
/* column 31 (df/dp_30) */
	jacp_[30] = 0; /* [0, 30] */
	jacp_[89] = 0; /* [1, 30] */
/* column 32 (df/dp_31) */
	jacp_[31] = 0; /* [0, 31] */
	jacp_[90] = 0; /* [1, 31] */
/* column 33 (df/dp_32) */
	jacp_[32] = 0; /* [0, 32] */
	jacp_[91] = 0; /* [1, 32] */
/* column 34 (df/dp_33) */
	jacp_[33] = 0; /* [0, 33] */
	jacp_[92] = 0; /* [1, 33] */
/* column 35 (df/dp_34) */
	jacp_[34] = 0; /* [0, 34] */
	jacp_[93] = 0; /* [1, 34] */
/* column 36 (df/dp_35) */
	jacp_[35] = 0; /* [0, 35] */
	jacp_[94] = 0; /* [1, 35] */
/* column 37 (df/dp_36) */
	jacp_[36] = 0; /* [0, 36] */
	jacp_[95] = 0; /* [1, 36] */
/* column 38 (df/dp_37) */
	jacp_[37] = (0-(((((CaM_0+PP1_0)-CaMKII_0)-(((((((((((CaM_Ca1+CaM_Ca2)+CaM_Ca3)+CaM_Ca4)+PP2B_CaM)+PP2B_CaM_Ca1)+PP2B_CaM_Ca2)+PP2B_CaM_Ca3)+PP2B_CaM_Ca4)+(PP1_0-PP1__pCaMKIIa))-((CaMKII_0-PP1_0)-(((((((((((CaMKII_CaM+CaMKII_CaM_Ca1)+CaMKII_CaM_Ca2)+CaMKII_CaM_Ca3)+CaMKII_CaM_Ca4)+pCaMKII_CaM_Ca4)+pCaMKIIa)+pCaMKII_CaM_Ca3)+pCaMKII_CaM_Ca2)+pCaMKII_CaM_Ca1)+pCaMKII_CaM)-(PP1_0-PP1__pCaMKIIa))))-pCaMKIIa))*pCaMKIIa)-((KD__CaM__Ca*(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)/KD__pCaMKII_CaM__Ca))))*pCaMKII_CaM))); /* [0, 37] */
	jacp_[96] = 0; /* [1, 37] */
/* column 39 (df/dp_38) */
	jacp_[38] = (0-((CaM_Ca1*pCaMKIIa)-((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)))*pCaMKII_CaM_Ca1))); /* [0, 38] */
	jacp_[97] = 0; /* [1, 38] */
/* column 40 (df/dp_39) */
	jacp_[39] = (0-((CaM_Ca2*pCaMKIIa)-((KD__CaM_Ca2__Ca*(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))*pCaMKII_CaM_Ca2))); /* [0, 39] */
	jacp_[98] = 0; /* [1, 39] */
/* column 41 (df/dp_40) */
	jacp_[40] = (0-((CaM_Ca3*pCaMKIIa)-(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)*pCaMKII_CaM_Ca3))); /* [0, 40] */
	jacp_[99] = 0; /* [1, 40] */
/* column 42 (df/dp_41) */
	jacp_[41] = 0; /* [0, 41] */
	jacp_[100] = 0; /* [1, 41] */
/* column 43 (df/dp_42) */
	jacp_[42] = 0; /* [0, 42] */
	jacp_[101] = 0; /* [1, 42] */
/* column 44 (df/dp_43) */
	jacp_[43] = (-1*((CaM_Ca4*pCaMKIIa)-(KD__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4))); /* [0, 43] */
	jacp_[102] = 0; /* [1, 43] */
/* column 45 (df/dp_44) */
	jacp_[44] = 0; /* [0, 44] */
	jacp_[103] = 0; /* [1, 44] */
/* column 46 (df/dp_45) */
	jacp_[45] = ((((0-(0-((((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca))*kf__CaM_Ca3__pCaMKIIa)*pCaMKII_CaM_Ca3)))-(0-((((KD__CaM_Ca2__Ca*((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca)))/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca)))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*((0-(KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa))/(KD__pCaMKII_CaM_Ca3__Ca*KD__pCaMKII_CaM_Ca3__Ca)))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 45] */
	jacp_[104] = 0; /* [1, 45] */
/* column 47 (df/dp_46) */
	jacp_[46] = (((0-(0-((((KD__CaM_Ca2__Ca*(0-((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)))/(KD__pCaMKII_CaM_Ca2__Ca*KD__pCaMKII_CaM_Ca2__Ca))*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(0-((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)))/(KD__pCaMKII_CaM_Ca2__Ca*KD__pCaMKII_CaM_Ca2__Ca)))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(0-((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)))/(KD__pCaMKII_CaM_Ca2__Ca*KD__pCaMKII_CaM_Ca2__Ca)))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 46] */
	jacp_[105] = 0; /* [1, 46] */
/* column 48 (df/dp_47) */
	jacp_[47] = ((0-(0-((((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(0-(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))))/(KD__pCaMKII_CaM_Ca1__Ca*KD__pCaMKII_CaM_Ca1__Ca))*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(0-(((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca))))/(KD__pCaMKII_CaM_Ca1__Ca*KD__pCaMKII_CaM_Ca1__Ca)))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 47] */
	jacp_[106] = 0; /* [1, 47] */
/* column 49 (df/dp_48) */
	jacp_[48] = (0-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*(KD__CaM_Ca1__Ca*(KD__CaM_Ca2__Ca*(0-((((KD__CaM_Ca3__Ca*KD__CaM_Ca4__pCaMKIIa)/KD__pCaMKII_CaM_Ca3__Ca)/KD__pCaMKII_CaM_Ca2__Ca)/KD__pCaMKII_CaM_Ca1__Ca)))))/(KD__pCaMKII_CaM__Ca*KD__pCaMKII_CaM__Ca))*pCaMKII_CaM)))); /* [0, 48] */
	jacp_[107] = 0; /* [1, 48] */
/* column 50 (df/dp_49) */
	jacp_[49] = (((((-1*(0-(kf__CaM_Ca4__pCaMKIIa*pCaMKII_CaM_Ca4)))-(0-(((KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca)*kf__CaM_Ca3__pCaMKIIa)*pCaMKII_CaM_Ca3)))-(0-((((KD__CaM_Ca2__Ca*(KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca)*kf__CaM_Ca2__pCaMKIIa)*pCaMKII_CaM_Ca2)))-(0-((((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca)*kf__CaM_Ca1__pCaMKIIa)*pCaMKII_CaM_Ca1)))-(kf__CaM__pCaMKIIa*(0-(((KD__CaM__Ca*((KD__CaM_Ca1__Ca*((KD__CaM_Ca2__Ca*(KD__CaM_Ca3__Ca/KD__pCaMKII_CaM_Ca3__Ca))/KD__pCaMKII_CaM_Ca2__Ca))/KD__pCaMKII_CaM_Ca1__Ca))/KD__pCaMKII_CaM__Ca)*pCaMKII_CaM)))); /* [0, 49] */
	jacp_[108] = 0; /* [1, 49] */
/* column 51 (df/dp_50) */
	jacp_[50] = 0; /* [0, 50] */
	jacp_[109] = 0; /* [1, 50] */
/* column 52 (df/dp_51) */
	jacp_[51] = (0-(pCaMKIIa*(PP1_0-PP1__pCaMKIIa))); /* [0, 51] */
	jacp_[110] = (pCaMKIIa*(PP1_0-PP1__pCaMKIIa)); /* [1, 51] */
/* column 53 (df/dp_52) */
	jacp_[52] = (0-(0-PP1__pCaMKIIa)); /* [0, 52] */
	jacp_[111] = (0-PP1__pCaMKIIa); /* [1, 52] */
/* column 54 (df/dp_53) */
	jacp_[53] = 0; /* [0, 53] */
	jacp_[112] = (0-PP1__pCaMKIIa); /* [1, 53] */
/* column 55 (df/dp_54) */
	jacp_[54] = 0; /* [0, 54] */
	jacp_[113] = 0; /* [1, 54] */
/* column 56 (df/dp_55) */
	jacp_[55] = (0-(kf__PP1__pCaMKIIa*pCaMKIIa)); /* [0, 55] */
	jacp_[114] = (kf__PP1__pCaMKIIa*pCaMKIIa); /* [1, 55] */
/* column 57 (df/dp_56) */
	jacp_[56] = 0; /* [0, 56] */
	jacp_[115] = 0; /* [1, 56] */
/* column 58 (df/dp_57) */
	jacp_[57] = (0-(kf__CaM__pCaMKIIa*(pCaMKIIa-(0*pCaMKII_CaM)))); /* [0, 57] */
	jacp_[116] = 0; /* [1, 57] */
/* column 59 (df/dp_58) */
	jacp_[58] = 0; /* [0, 58] */
	jacp_[117] = 0; /* [1, 58] */
	return GSL_SUCCESS;
}
/* ode Functions F(t,y;p) */
int CaMKIIs_func(double t, const double y_[], double *func_, void *par)
{
	double *p_=par;
	if (!y_ || !func_) return 4;
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
	double logistic=1.0/(1+exp(-10*t));
	double Ca=logistic*Ca_set;
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
	func_[0] = BoundCa / (s + Total_CaM); /* MolCaPerMolCaM */
	func_[1] = 100*(Total_pCaMKII / (s + totalCaMKII)); /* AutoCamkiiPercentage */
	func_[2] = Total_PP2B_CaM_CaX / (s + PP2B + Total_PP2B_CaM_CaX); /* MolCaMPerMolPP2B */
	func_[3] = 100 * (PP2B_CaM_Ca4 / (s + PP2B + Total_PP2B_CaM_CaX)); /* ActivePP2BPercentage */
	return GSL_SUCCESS;
}
/* ode default parameters */
int CaMKIIs_default(double t, void *par)
{
	double *p_=par;
	if (!p_) return 59;
	double a=3.90264;
	double b=2.86972;
	double s=1e-05;
	p_[0] = 0.15191;
	p_[1] = 3.4245e-05;
	p_[2] = 0.073893;
	p_[3] = 0.0061814;
	p_[4] = 0.020293;
	p_[5] = 0.0045248;
	p_[6] = 0.051176;
	p_[7] = 0.27421;
	p_[8] = 0.083357;
	p_[9] = 0.0011578;
	p_[10] = 0.0047884;
	p_[11] = 0.035079;
	p_[12] = 0.045566;
	p_[13] = 7271.3;
	p_[14] = 37062;
	p_[15] = 1827.9;
	p_[16] = 2662.3;
	p_[17] = 0.03997;
	p_[18] = 91.543;
	p_[19] = 916.15;
	p_[20] = 285.03;
	p_[21] = 82.837;
	p_[22] = 0.23745;
	p_[23] = 0.025858;
	p_[24] = 0.13086;
	p_[25] = 0.075539;
	p_[26] = 0.00079772;
	p_[27] = 0.055882;
	p_[28] = 0.046028;
	p_[29] = 0.20855;
	p_[30] = 0.022662;
	p_[31] = 8.2849;
	p_[32] = 483.48;
	p_[33] = 1143.6;
	p_[34] = 645.07;
	p_[35] = 3081.6;
	p_[36] = 0.00082984;
	p_[37] = 0.00032583;
	p_[38] = 0.058928;
	p_[39] = 0.02319;
	p_[40] = 0.030252;
	p_[41] = 0.038498;
	p_[42] = 0.0004565;
	p_[43] = 0.072237;
	p_[44] = 0.0021267;
	p_[45] = 539.41;
	p_[46] = 1784;
	p_[47] = 57728;
	p_[48] = 1342.9;
	p_[49] = 3.746;
	p_[50] = 0.0037559;
	p_[51] = 0.0016604;
	p_[52] = 0.20517;
	p_[53] = 0.30225;
	p_[54] = 2187.8;
	p_[55] = 0;
	p_[56] = 0;
	p_[57] = 30;
	p_[58] = 3;
	return GSL_SUCCESS;
}
/* ode initial values */
int CaMKIIs_init(double t, double *y_, void *par)
{
	double *p_=par;
	if (!y_) return 21;
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
	y_[11] = 0;
	y_[12] = 0;
	y_[13] = 0;
	y_[14] = 0;
	y_[15] = 0;
	y_[16] = 0;
	y_[17] = 0;
	y_[18] = 0;
	y_[19] = 0;
	y_[20] = 0;
	return GSL_SUCCESS;
}
