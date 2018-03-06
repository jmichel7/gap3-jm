CC=cc
CFLAGS=-O -DSYS_UNIX
SMALL=-DMAXN=WORDSIZE
MEDIUM=-DMAXN=8192
LARGE=-DMAXN=8192

dreadnaut: dreadnaut.o naututil.o nauty.o nautil.o nautinv.o
	${CC} ${LDFLAGS} -o ../bin/dreadnaut dreadnaut.o naututil.o \
		nauty.o nautil.o nautinv.o
dreadnaut.o: nauty.h naututil.h dreadnaut.c
	${CC} ${CFLAGS} -c ${MEDIUM} dreadnaut.c
naututil.o: nauty.h naututil.h naututil.c
	${CC} ${CFLAGS} -c ${MEDIUM} naututil.c
nautil.o: nauty.h nautil.c
	${CC} ${CFLAGS} -c ${MEDIUM} nautil.c
nauty.o: nauty.h nauty.c
	${CC} ${CFLAGS} -c ${MEDIUM} nauty.c
nautinv.o: nauty.h naututil.h nautinv.c
	${CC} ${CFLAGS} -c ${MEDIUM} nautinv.c
nautaux.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} ${CFLAGS} -c ${MEDIUM} nautaux.c

dreadnaut1: dreadnaut1.o naututil1.o nauty1.o nautil1.o nautinv1.o
	${CC} ${LDFLAGS} -o dreadnaut1 dreadnaut1.o naututil1.o nauty1.o \
		nautil1.o nautinv1.o
dreadnaut1.o: nauty.h naututil.h dreadnaut.c
	${CC} ${CFLAGS} -c ${SMALL} -o dreadnaut1.o dreadnaut.c
naututil1.o: nauty.h naututil.h naututil.c
	${CC} ${CFLAGS} -c ${SMALL} -o naututil1.o naututil.c
nautil1.o: nauty.h nautil.c
	${CC} ${CFLAGS} -c ${SMALL} -o nautil1.o nautil.c
nauty1.o: nauty.h nauty.c
	${CC} ${CFLAGS} -c ${SMALL} -o nauty1.o nauty.c
nautinv1.o: nauty.h naututil.h nautinv.c
	${CC} ${CFLAGS} -c ${SMALL} -o nautinv1.o nautinv.c
nautaux1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} ${CFLAGS} -c ${SMALL} -o nautaux1.o nautaux.c

dreadnautx: dreadnautx.o naututilx.o nautyx.o nautilx.o nautinvx.o
	${CC} ${LDFLAGS} -o dreadnautx dreadnautx.o naututilx.o nautyx.o \
		nautilx.o nautinvx.o
dreadnautx.o: nauty.h naututil.h dreadnaut.c
	${CC} ${CFLAGS} -c ${LARGE} -o dreadnautx.o dreadnaut.c
naututilx.o: nauty.h naututil.h naututil.c
	${CC} ${CFLAGS} -c ${LARGE} -o naututilx.o naututil.c
nautilx.o: nauty.h nautil.c
	${CC} ${CFLAGS} -c ${LARGE} -o nautilx.o nautil.c
nautyx.o: nauty.h nauty.c
	${CC} ${CFLAGS} -c ${LARGE} -o nautyx.o nauty.c
nautinvx.o: nauty.h naututil.h nautinv.c
	${CC} ${CFLAGS} -c ${LARGE} -o nautinvx.o nautinv.c
nautauxx.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} ${CFLAGS} -c ${LARGE} -o nautauxx.o nautaux.c

rancubic.o : nauty.h rancubic.c
	${CC} ${CFLAGS} -c rancubic.c

cif3.o : nauty.h naututil.h cif3.c
	${CC} ${CFLAGS} -c cif3.c
cif3 : cif3.o rancubic.o naututil.o nautilx.o
	${CC} ${LDFLAGS} -o cif3 cif3.o rancubic.o naututilx.o nautilx.o

clean:
	rm -f *.o ../bin/dreadnaut
