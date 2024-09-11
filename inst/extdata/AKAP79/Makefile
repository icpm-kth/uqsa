CC = gcc
CFLAGS = -O2 `pkg-config --cflags gsl`
LIBS = `pkg-config --libs gsl`

CC = gcc
ODE = ~/RPN-derivative/sh/ode.sh
VF = ~/SBtabVFGEN/sbtab_to_vfgen

.PHONY: all

all: $(addprefix AKAP79, _gvf.c .R)


%.so: %_gvf.c
	$(CC) -shared -fPIC $(CFLAGS) $^ $(LIBS) -o $@

%.R: %.vf
	$(ODE) --maxima -R $^ > $@

%_gvf.c: %.vf
	$(ODE) --maxima -C $^ > $@

%.vf: *.tsv
	$(VF) --cla $^
