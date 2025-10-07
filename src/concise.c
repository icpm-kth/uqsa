#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rdefines.h>

#define Rdata SEXP

int numeric(char c){
	if ('0' <= c && c <= '9') return 1;
	else if (c=='.' || c=='+' || c=='-') return 1;
	else return 0;
}

/* numbe a = v ± u */
struct num {
	union{ /* the estimate of something */
		double v;
		double value;
		double expectedValue;
		double mean;
	};
	union { /* the uncertainty of the estimate */
		double u;
		double error;
		double stdDev;
		double uncertainty;
	};
};

/* this will only read: 12.34(56)E-7, like this, with an E*/
struct num read_concise(const char *line){
	struct num a;
	char *ptr_u=line, *ptr_e=line;
	int vscale=0;
	double v=strtod(line,&ptr_u);
	int u=0;
	//int ulen=2;
	char *decimal=strchr(line,'.');
	int digits=ptr_u-(decimal?decimal:line)-1;
	while (ptr_u && *ptr_u && !numeric(*ptr_u)) ptr_u++;
	if (numeric(*ptr_u)) {
		u=strtol(ptr_u,&ptr_e,10);
	}
	while (ptr_e && *ptr_e && !numeric(*ptr_e)) ptr_e++;
	if (numeric(*ptr_e)) vscale=strtol(ptr_e,NULL,10);
	a.v=v*pow(10,vscale);
	a.u=u*pow(10,-digits+vscale);
	return a;
}

struct num read_number(const char *line){
	struct num a;
	char *p=line;
	char *q;
	a.value = strtod(line,&p);
	char *decimal = strchr(line,'.');
	int digits = decimal?(p-decimal-1):0;
	while (p && *p && !numeric(*p)) p++;
	if (p && *p){
		a.error = strtod(p,&q);
		if (p==q) a.error = INFINITY;
	} else {
		a.error = pow(10,-digits);
	}
	return a;
}

Rdata concise(Rdata charVector){
	int j;
	struct num a;
	const char *line;
	Rdata out = PROTECT(allocMatrix(REALSXP,2,length(charVector)));
	double *d = REAL(out);
	for (j=0;j<length(charVector);j++){
		line = CHAR(STRING_ELT(charVector,j));
		if (strchr(line,'(')){
			a=read_concise(line);
		} else {
			a=read_number(line);
		}
		memcpy(d,&a,sizeof(a));
		d+=2;
	}
	UNPROTECT(1);
	return out;
}
