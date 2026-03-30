#include <stdlib.h>
#include <math.h>

/* System's volume at creation time: */
static const double sys_volume =            1e-15;
/* Product of Avogadro's number L (6.02214076E23) and volume V: */
static const double LV =      6.02214e+08;
/* This model was created from an SBtab file,         */
/* which always describes a systems biology model of  */
/* concentrations and kinetic laws. This code         */
/* re-interprets it as particle counts and stochastic */
/* reaction propensities. This interpretation may be  */
/*       >>INCORRECT<<       ... as it is automatic.  */
/* The plan is to use the kinetic law literally,      */
/* as specified in the SBtab file, but to convert the */
/* rate coefficients.                                 */
/* We also convert initial values to particle counts. */

/* these enums make it possible to address vector elements by name, and automatically creates lengths for these vectors*/
enum state {_AKAR4, _AKAR4_C, _AKAR4p, _C, numStateVariables};
enum parameter {_kf_C_AKAR4, _kb_C_AKAR4, _kcat_AKARp, numParameters};
enum reaction {_reaction_1_fwd, _reaction_2_fwd, _reaction_1_bwd, _reaction_2_bwd, numReactions};
enum outputFunctions {_AKAR4pOUT, numFunctions};

int model_effects(double t, int *x, int j){
	if (!x) return numStateVariables;
	switch(j){
	case _reaction_1_fwd:
		x[_C] -= 1;
		x[_AKAR4] -= 1;
		x[_AKAR4_C] += 1;
		break;
	case _reaction_2_fwd:
		x[_AKAR4_C] -= 1;
		x[_AKAR4p] += 1;
		x[_C] += 1;
		break;
	case _reaction_1_bwd:
		x[_AKAR4_C] -= 1;
		x[_C] += 1;
		x[_AKAR4] += 1;
		break;
	case _reaction_2_bwd:
		x[_AKAR4p] -= 1;
		x[_C] -= 1;
		x[_AKAR4_C] += 1;
		break;
	}
	return 0;
}

/* Here, we want to use the kinetic law literally, with no variable substitutions.*/
/* So, 'k*A*B' will be used as written, but the stochastic propensity is c*A*B    */
/* where c is different from k in value and unit. We pre-calculate c/k.           */
/* We convert k to have the right value (and unit, implicitly)                    */
/* and use the converted value under the name k directly so as not to change the  */
/* rate. This may look misleading to the reader compared to a stochastic model    */
/* made from scratch.                                                             */
int model_propensities(double t, int *x, double *c, double *a){
	if (!x || !a) return numReactions;
	/* constants */
	/* parameters */
	double kf_C_AKAR4 = c[_kf_C_AKAR4];                      /* originally: 1/uM*s */
	double kb_C_AKAR4 = c[_kb_C_AKAR4];                      /* originally: 1/s */
	double kcat_AKARp = c[_kcat_AKARp];                      /* originally: 1/s */
	/*state variables */
	double AKAR4 = x[_AKAR4];
	double AKAR4_C = x[_AKAR4_C];
	double AKAR4p = x[_AKAR4p];
	double C = x[_C];
	a[_reaction_1_fwd] = 1*kf_C_AKAR4*C*AKAR4;              /* 1 C + 1 AKAR4 -> 1 AKAR4_C */
	a[_reaction_2_fwd] = 1*kcat_AKARp*AKAR4_C;              /* 1 AKAR4_C -> 1 AKAR4p + 1 C */
	a[_reaction_1_bwd] = 1*kb_C_AKAR4*AKAR4_C;              /* 1 AKAR4_C -> 1 C + 1 AKAR4 */
	a[_reaction_2_bwd] = 1*0.0;                             /* 1 AKAR4p + 1 C -> 1 AKAR4_C */
	return 0;
}

/* usually, these are stochastic propensity coefficients,                                      */
/* but we derive our model from a concentration based kinetic form.                            */
/* So, these are not directly usable, but we'll write the propensities with conversion factors */
int model_reaction_coefficients(double *c){
	if (!c) return numParameters;
	c[_kf_C_AKAR4] = 0.018;                           /* 1/uM*s */
	c[_kb_C_AKAR4] = 0.106;                           /* 1/s */
	c[_kcat_AKARp] = 10.2;                            /* 1/s */
	return 0;
}

int model_initial_counts(int *x, double *c){
	if(!x) return numStateVariables;
	/* constants */
	/* parameters */
	double kf_C_AKAR4 = c[_kf_C_AKAR4];                      /* originally: 1/uM*s */
	double kb_C_AKAR4 = c[_kb_C_AKAR4];                      /* originally: 1/s */
	double kcat_AKARp = c[_kcat_AKARp];                      /* originally: 1/s */
	/*state variables */
	x[_AKAR4] = lround(2e-07 * LV);                                /* originally in µM */
	x[_AKAR4_C] = lround(0 * LV);                                  /* originally in µM */
	x[_AKAR4p] = lround(0 * LV);                                   /* originally in µM */
	x[_C] = lround(0 * LV);                                        /* originally in µM */
	return 0;
}

/* This function will try to convert the model's state to concentrations */
/* We assume that the model was originally phrased as concentrations and kinetic laws. */
int model_func(double t, int *x, double *c, double *func){
	if(!func) return numFunctions;
	/* constants */
	/* parameters */
	double kf_C_AKAR4 = c[_kf_C_AKAR4];                      /* originally: 1/uM*s */
	double kb_C_AKAR4 = c[_kb_C_AKAR4];                      /* originally: 1/s */
	double kcat_AKARp = c[_kcat_AKARp];                      /* originally: 1/s */
	/*state variables */
	double AKAR4 = ((double) x[_AKAR4])/(LV * 1e-06);                           /* µM */
	double AKAR4_C = ((double) x[_AKAR4_C])/(LV * 1e-06);                       /* µM */
	double AKAR4p = ((double) x[_AKAR4p])/(LV * 1e-06);                         /* µM */
	double C = ((double) x[_C])/(LV * 1e-06);                                   /* µM */
	func[_AKAR4pOUT] = 108 + 380*AKAR4p;
	return 0;
}

/* given kinetic parameters, calculate the stochastic parameters needed for propensities */
/* The changes are made in place */
int model_stochastic_parameters(double t, double *par){
	if (!par) return numParameters;
	par[_kf_C_AKAR4] = (par[_kf_C_AKAR4]*1e+06)/LV;
	par[_kb_C_AKAR4] = (par[_kb_C_AKAR4]*1);
	par[_kcat_AKARp] = (par[_kcat_AKARp]*1);
	return 0;
}

/* The list of experiments we'll receive from R will cointain */
/* initial concentrations rather than particle counts.        */
/* We need to convert them.                                   */
int model_particle_count(double t, double *molarity, int *x){
	if(!molarity || !x) return numStateVariables;
	x[_AKAR4] = lround(1e-06 * molarity[_AKAR4] * LV);
	x[_AKAR4_C] = lround(1e-06 * molarity[_AKAR4_C] * LV);
	x[_AKAR4p] = lround(1e-06 * molarity[_AKAR4p] * LV);
	x[_C] = lround(1e-06 * molarity[_C] * LV);
	return 0;
}
