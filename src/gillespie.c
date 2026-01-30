#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <dlfcn.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

typedef SEXP Rdata;

int (*model_effects)(double t, int *x, int j);
int (*model_propensities)(double t, int *x, double *c, double*a);
int (*model_reaction_coefficients)(double *c);
int (*model_initial_counts)(int *x);
int (*model_func)(double t, int *x, double *c, double *f);
int (*model_particle_count)(double t, double *molarity, int *x);
int (*model_particle_count)(double t, double *molarity, int *x);
int (*model_stochastic_parameters)(double t, double *par);

void* load_model(const char *model_so){
	void *lib=dlopen(model_so,RTLD_LAZY);
	char *model_so_2=NULL;
	if (!lib) {
		fprintf(stderr,"[%s] %s.\n",__func__,dlerror());
		model_so_2 = malloc(strlen(model_so)+3);
		stpcpy(stpcpy(model_so_2,"./"),model_so);
		fprintf(stderr,"[%s] retrying with: «%s»\n",__func__,model_so_2);
		lib = dlopen(model_so_2,RTLD_LAZY);
	}
	if (lib){
		if ((model_effects=dlsym(lib,"model_effects"))==NULL){
			fprintf(stderr,"[%s] loading model_effects has failed.\n",__func__);
		}
		if ((model_propensities=dlsym(lib,"model_propensities"))==NULL){
			fprintf(stderr,"[%s] loading model_propensities has failed.\n",__func__);
		}
		if ((model_reaction_coefficients=dlsym(lib,"model_reaction_coefficients"))==NULL){
			fprintf(stderr,"[%s] loading model_reaction_coefficients has failed.\n",__func__);
		}
		if ((model_initial_counts=dlsym(lib,"model_initial_counts"))==NULL){
			fprintf(stderr,"[%s] loading model_initial_counts has failed.\n",__func__);
		}
		if ((model_func=dlsym(lib,"model_func"))==NULL){
			fprintf(stderr,"[%s] loading model_func has failed.\n",__func__);
		}
		if ((model_particle_count=dlsym(lib,"model_particle_count"))==NULL){
			fprintf(stderr,"[%s] loading model_particle_count has failed.\n",__func__);
		}
		if ((model_stochastic_parameters=dlsym(lib,"model_stochastic_parameters"))==NULL){
			fprintf(stderr,"[%s] loading model_stochastic_parameters has failed.\n",__func__);
		}
	}
	if (model_so_2) free(model_so_2);
	return lib;
}

int pick_reaction(gsl_vector *a, double r_sum_a){
	int j;
	size_t m=a->size;
	double psum = 0.0;
	//gsl_sort_vector(a);
	for (j=0; j<m; j++){
		psum += gsl_vector_get(a,j);
		if (psum > r_sum_a) break;
	}
	return j;
}

void print_counts(double t, gsl_vector_int *x){
	printf("%10f ",t);
	size_t n = x->size;
	int i;
	for (i=0;i<n;i++){
		printf("%10i ",x->data[i]);
	}
	printf("\n");
}

void advance_in_time(double *t, gsl_vector_int *x, gsl_vector *c, gsl_vector *a, gsl_rng *RNG){
	double tau;
	double r;
	double sum_a=0.0;
	int j;
	int i;
	model_propensities(*t,x->data,c->data,a->data);
	// 1. calculate sum(a)
	sum_a = gsl_vector_sum(a);
	// 2. pick a time step forward (tau)
	tau = gsl_ran_exponential(RNG,1.0/sum_a); /* -log(gsl_rng_uniform(RNG))/sum_a; */
	// 3. pick a reaction to perform:
	r = gsl_rng_uniform(RNG);
	j = pick_reaction(a,r*sum_a);
	model_effects(*t,x->data,j);
	*t += tau;
}

int in_list(Rdata List, const char *name);
Rdata from_list(Rdata List, const char *name);
void set_names(Rdata list, const char *names[]);

int cat_parameters(gsl_vector *c, double *p, size_t np, Rdata input){
	memcpy(c->data,p,sizeof(double)*np);
	if (input && input != R_NilValue) {
		memcpy(c->data+(c->size - length(input)),REAL(AS_NUMERIC(input)),sizeof(double)*length(input));
	}
	return 0;
}

Rdata gillespie(Rdata model_so, Rdata experiments, Rdata parameters){
	gsl_set_error_handler_off();
	int i,j,k;
	//int l;
	void *handle = load_model(CHAR(STRING_ELT(model_so,0)));
	if (!handle) return R_NilValue;
	size_t M = ncols(parameters);
	size_t nrp = nrows(parameters);
	size_t n = model_initial_counts(NULL);
	size_t m = model_reaction_coefficients(NULL);
	size_t na = model_propensities(0,NULL,NULL,NULL);
	size_t nf = model_func(0,NULL,NULL,NULL);
	/* printf("[%s] sizes: %li stateVars, %li parameters (%li columns), %li propensities, and %li functions.\n",__func__,n,m,M,na,nf); */
	fflush(stdout);
	if (m < nrp){
		fprintf(stderr,"[%s] too many parameters supplied (%i) for model (%li).\n",__func__,length(parameters),m);
	}
	gsl_vector_int *x = gsl_vector_int_alloc(n);
	gsl_vector *c = gsl_vector_alloc(m);
	gsl_vector *krc = gsl_vector_alloc(m);
	gsl_vector *a = gsl_vector_alloc(na);
	size_t nu;
	double t0 = 0;
	double tf = 10;
	double t = t0;
	double *p;
	Rdata time, input, E, initialState, initialTime, yf, y, f;
	const char *snames[] = {"state","func",NULL};
	Rdata solution = PROTECT(NEW_LIST(length(experiments)));
	SET_NAMES(solution,GET_NAMES(experiments));
	gsl_vector_int_view column;

	gsl_rng *RNG = gsl_rng_alloc(gsl_rng_ranlxs0);
	gsl_rng_set(RNG, 1337);
	model_initial_counts(x->data);
	model_reaction_coefficients(c->data);
	/* main solver loop: */
	for (i=0; i<length(experiments); i++){
		E = VECTOR_ELT(experiments,i);
		initialState = from_list(E,"initialState");
		initialTime = from_list(E,"initialTime");
		time = from_list(E,"outputTimes");
		input = from_list(E,"input");
		t0 = *REAL(AS_NUMERIC(initialTime));
		yf = PROTECT(NEW_LIST(2));
		y = PROTECT(alloc3DArray(INTSXP,n,length(time),M));
		f = PROTECT(alloc3DArray(REALSXP,nf,length(time),M));
		for (k=0;k<M;k++){
			t = t0;
			p = REAL(parameters)+(k*nrp);
			cat_parameters(c,p,nrp,input); /* c <- cat(p,input)*/
			gsl_vector_memcpy(krc,c);
			if (model_stochastic_parameters){
				model_stochastic_parameters(t0,c->data);
			}
			if (model_particle_count){
				model_particle_count(t0,REAL(initialState),x->data);
			}
			for (j=0; j<length(time); j++){
				tf = REAL(AS_NUMERIC(time))[j];
				while (t < tf){
					advance_in_time(&t,x,c,a,RNG);
				}
				column = gsl_vector_int_view_array(INTEGER(y)+(n*length(time)*k+n*j),n);
				model_func(t,x->data,krc->data,REAL(f)+(nf*length(time)*k+nf*j));
				gsl_vector_int_memcpy(&(column.vector),x);
			}
		}
		SET_VECTOR_ELT(yf,0,y);
		SET_VECTOR_ELT(yf,1,f);
		set_names(yf,snames);
		UNPROTECT(2);
		SET_VECTOR_ELT(solution,i,yf);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	gsl_rng_free(RNG);
	return solution;
}
