MODEL = AKAP79
CC = cc
CFLAGS = -O2 `pkg-config --cflags gsl`
LIBS = `pkg-config --libs gsl`


.PHONY: all

all: $(MODEL).so


$(MODEL).so: $(MODEL)_gvf.c
	$(CC) -shared -fPIC $(CFLAGS) $^ $(LIBS) -o $@

$(MODEL)_gvf.c: $(MODEL).vf
	ode $^ > $@


