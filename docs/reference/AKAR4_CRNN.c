#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* This is a general mass action law model,                         */
/* with a fixed number of reactions, where                          */
/* a reaction:                                                      */
/*     n[0] X[0] + n[1] X[1] -> n[2] X[2],                          */
/* has the stoichiometry                                            */
/*     nu = {-n[0]; -n[1]; +n[2]},                                  */
/* and reaction flux:                                               */
/*     flux = exp(log[k] + n[0]*log(X[0]) + n[1]*log(X[1])),        */
/* which is the same as:                                            */
/*     flux = exp(l - nu[0]*log(X[0]) - nu[1]*log(X[1])),           */
/* but expressed in terms of stoichiometry.                         */
/*                                                                  */
/* Many quantities are expressed in logarithmic space, e.g.:        */
/*     l = log(k).                                                  */
/* Specifically, reaction fluxes have to be sums of logarithms.     */
/* For this reason, a reversible reaction has to be split up        */
/* into two reactions:                                              */
/*     log(kf*A*B) = log(kf) + log(A) + log(B), whereas             */
/* log(kf*A*B-kr*C) does not split up so neatly.                    */
/* The stoichiometry nu is a matrix in column major order, i.e.     */
/* nu[a,b] = nu[a+A*b], where a and are offsets 0,...,A-1;          */
/* a enumerates the state variables. and b the reactions.           */
/*                                                                  */
/* The stoichiometry of reverse reations nu_b = -nu is implied.     */
/* We re-use nu for backward reaction fluxes, without copying it.   */
/*                                                                  */
/* Reactions can involve a modifier: a reactant that                */
/* isn't consumed by the reaction (enzymes). These                  */
/* wouldn't normally appear in the stoichiometry                    */
/* because their number is conserved. But, they do affect           */
/* the reaction flux.                                               */
/*  For these reasons, there is a second matrix of modifiers:       */
/*     m[i+numStV*j] = {0,1,2,...};                                 */
/* The values of m are only used in flux calculations.              */
/*                                                                  */
/* This log function returns -∞ for 0   */
/* i.e.: LOG(0) = GSL_NEGINF (-∞)       */
/* log(0) returns NaN.                  */
inline static double LOG(double y){
	return y>0?log(y):GSL_NEGINF;
}

enum stateVars {_AKAR4, _AKAR4_C, _AKAR4p, _C, numStV};
enum outputFun {_AKAR4pOUT, numFun};
/* const int numStV = 4; */ /* number of state variables, see above */
const int numPar = 4; /* number of reaction rate coefficients */
const int numRct = 2; /* number of reactions */
struct par {
	double *l;  /* log(k), rate coefficients */
	double *nu; /* stoichiometry (4×4)*/
	double *m;  /* modifiers (enzymes); not consumed or produced */
};

int CRNN_flux(double t, double *y, double *fwdFlux, double *bwdFlux, struct par *p){
	double *nu=p->nu;
	double *l=p->l;
	double *m=p->m;
	int i,j;
/* 	forward rections   */
	for (i = 0; i < numRct; i++) {
		fwdFlux[i] = l[i];
		for (j = 0; j < numStV; j++) {
			if (nu[j+numStV*i]<0) {
				fwdFlux[i] -= nu[j+numStV*i]*LOG(y[j]); /* a reactant */
			}
			if (m[j+numStV*i]>0){
				fwdFlux[i] +=  m[j+numStV*i]*LOG(y[j]); /* a modifier */
			}
		}
		fwdFlux[i] = exp(fwdFlux[i]);
	}
/* 	backward rections   */
	for (i = 0; i < numRct; i++) {
		bwdFlux[i] = l[i+numRct];
		for (j = 0; j < numStV; j++) {
			if (nu[j+numStV*i]>0) {
				bwdFlux[i] += nu[j+numStV*i]*LOG(y[j]); /* a reactant */
			}
			if (m[j+numStV*i]>0){
				bwdFlux[i] +=  m[j+numStV*i]*LOG(y[j]); /* a modifier */
			}
		}
		bwdFlux[i] = exp(bwdFlux[i]);
	}
	return GSL_SUCCESS;
}
int CRNN_vf(double t, double *y, double *f, void *par){
	if (!y || !f) return(numStV);
	int i,j;
	struct par *p = par;
	double *nu=p->nu;
	double fwdFlux[numRct];
	double bwdFlux[numRct];
	CRNN_flux(t,y,fwdFlux,bwdFlux,p);
	for (i = 0; i < numStV; i++){
		f[i] = 0.0;
/*		if (y[i]<1e-10) fprintf(stderr,"[%s] warning: y[%i] is too low: %g\n",__func__,i,y[i]);*/
		for (j = 0; j < numRct; j++) {
/* 			backward reactions have the inverse stoichiometry: -nu   */
			f[i] += nu[i+numStV*j]*(fwdFlux[j] - bwdFlux[j]);
		}
	}
	return GSL_SUCCESS;
}

int CRNN_func(double t, double *y, double *func, void *par){
	if (!y || !func) return numFun;
	struct par *p=par;
	double AKAR4 = y[_AKAR4];
	double AKAR4_C = y[_AKAR4_C];
	double AKAR4p = y[_AKAR4p];
	double C = y[_C];
	func[_AKAR4pOUT] = 108 + 380*AKAR4p;
	return GSL_SUCCESS;
}

int CRNN_jac(double t, double *y, double *jac, double dfdt[], void *par){
	if (!y || !jac) return numStV*numStV;
	int i,j,k;
	struct par *p=par;
	double *l=p->l;
	double *nu=p->nu;
	double *m=p->m;
	double fwdFlux[numRct];
	double bwdFlux[numRct];
	CRNN_flux(t,y,fwdFlux,bwdFlux,p);
	for (i = 0; i < numStV; i++){
		for (j = 0; j < numStV; j++) {
			jac[i*numStV+j]=0;
			for (k = 0; k < numRct; k++) {
				if (y[j]>0) {
					jac[i*numStV+j] -= nu[i+numStV*k]*(fwdFlux[k]*(nu[j+numStV*k]<0)+bwdFlux[k]*(nu[j+numStV*k]>0))*nu[j+numStV*k]/y[j];
				}
			}
		}
	}
	return GSL_SUCCESS;
}

/* These are the default initial conditions for y.   */
/* They may depend on the parameters,                */
/* and time of initialization t.                     */
int CRNN_init(double t, double *y, void *par){
	if(!y || !par) return(numStV);
	struct par *p=par;
	y[_AKAR4] = 0.2;
	y[_AKAR4_C] = 0;
	y[_AKAR4p] = 0;
	y[_C] = 0;
	return GSL_SUCCESS;
}

/* These are default values for the parameters.    */
/* They may depend on the initialization time t.   */
int CRNN_default(double t, void *par){
	if(!par) return numRct;
	double *p=par;
/*	It is a bit unclear what to do here for CRNNs, as they have more complex parameters. */
/*	memset(p,0,lpar*sizeof(double)); */
	return GSL_SUCCESS;
}
