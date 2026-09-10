#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

/* This is a general mass action law model,                         */
/* with a free number of reactions, where                           */
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

struct par {
	double *l;        /* log(k), rate coefficients                     */
	double *nu;       /* stoichiometry                                 */
	double *mu;       /* modifiers (enzymes); not consumed or produced */
	int numStV;       /* number of state variables                     */
	int numPar;       /* number of parameters                          */
	int numRct;       /* number of reactions                           */
	size_t dim_l[2];  /* dimensions of l                               */
	size_t dim_nu[2]; /* dimensions of nu                              */
	size_t dim_m[2];  /* dimensions of m                               */
};

int CRNN_vf(double t, double *y, double *f, void *par){
	if (!y || !f) return(numStV);
	int i,j,k;
	struct par *p = par;
	double *nu=p->nu;
	const int numRct=p->numRct;
	const int numStV=p->numStV;
	double fwdFlux;
	double bwdFlux;
	double netFlux;
	CRNN_flux(t,y,fwdFlux,bwdFlux,p);
	memset(f,0,numStV*sizeof(double));
	for (j = 0; j < numRct; j++) {
		for (i = 0; i < numStV; i++) {
			if (nu[i+numStV*j] < 0.0) {
				fwdFlux -= nu[i+numStV*j]*log(y[i]); /* a reactant */
			}
			if (m[i+numStV*j] > 0.0){
				fwdFlux += mu[i+numStV*j]*log(y[i]); /* a modifier */
				bwdFlux += mu[i+numStV*j]*log(y[i]); /* a modifier */
			}
			if (nu[i+numStV*j] > 0.0) {
				bwdFlux += nu[i+numStV*j]*log(y[i]); /* a product */
			}
		}
		fwdFlux = exp(fwdFlux);
		bwdFlux = exp(bwdFlux);
		netFlux = fwdFlux - bwdFlux;
		for (i = 0; i < numStV; i++) {
			f[i] += nu[i+numStV*j]*netFlux;
		}
	}
	return GSL_SUCCESS;
}

int CRNN_jac(double t, double *y, double *jac, double dfdt[], void *par){
	if (!y || !jac) return numStV*numStV;
	int i,j,k;
	struct par *p=par;
	double *l=p->l;
	double *nu=p->nu;
	double *m=p->m;
	const int numRct=p->numRct;
	const int numStV=p->numStV;
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
	const int numStV=p->numStV;
	memset(y,0,numStV*sizeof(double));
	return GSL_SUCCESS;
}

/* These are default values for the parameters.    */
/* They may depend on the initialization time t.   */
int CRNN_default(double t, void *par){
	if(!par) return numRct;
	double *p=par;
	const int numStV=p->numStV;
	const int mumRct=p->numRct;
	const int mumPar=p->numPar;
	memset(p->l,0,(p->dim_l[0])*(p->dim_l[1])*sizeof(double));
	memset(p->nu,0,(p->dim_nu[0])*(p->dim_nu[1])*sizeof(double));
	return GSL_SUCCESS;
}
