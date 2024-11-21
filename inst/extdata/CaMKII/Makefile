CC = gcc
CFLAGS = -O3 -march=native
model = CaMKII

.PHONY: all


all: $(model)_gvf.c $(model).R $(model).so


%.R: %.tar.gz
	ode -R --maxima $^ > $@

%_gvf.c: %.tar.gz
	ode -C --maxima $^ > $@

%.tar.gz: *.tsv
	sbtab_to_vfgen $^

%.so: %_gvf.c
	$(CC) -shared -fPIC $(CFLAGS) `pkg-config --cflags gsl`-o $@ $^ `pkg-config --libs gsl`
