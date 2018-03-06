# src makefile
SHELL   = /bin/sh
SYSTEM = SYSTEMV
CC      = gcc
RM = /bin/rm
LIBPATH = ../lib
BIN = ../bin
FSALIB = $(LIBPATH)/fsalib.a

COPTS = -O
CFLAGS = $(COPTS) -I$(LIBPATH)
forgap: $(BIN)/kbprog $(BIN)/ppgap $(BIN)/gpmakefsa $(BIN)/gpaxioms $(BIN)/gpminkb

all:    $(BIN)/kbprog \
	$(BIN)/fsamin $(BIN)/fsabfs $(BIN)/fsafilter \
	$(BIN)/fsaand $(BIN)/fsaor $(BIN)/fsanot $(BIN)/fsaexists \
        $(BIN)/fsacount $(BIN)/fsaenumerate\
	$(BIN)/gpwa $(BIN)/gpwab $(BIN)/gpgenmult $(BIN)/gpcheckmult \
	$(BIN)/gpmult $(BIN)/gpgenmult2 $(BIN)/gpmult2 \
	$(BIN)/gpgenmult3 $(BIN)/gpcomp $(BIN)/gpaxioms $(BIN)/gpminkb \
	$(BIN)/gpmakefsa $(BIN)/wordreduce $(BIN)/ppgap

.c.o:
	${CC} -c $(CFLAGS) $*.c

${BIN}/kbprog: kbprog.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/kbprog kbprog.o $(FSALIB)

kbprog.o: kbprog.c
	$(CC) -c $(CFLAGS) -D$(SYSTEM) kbprog.c

${BIN}/fsamin: fsamin.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsamin  fsamin.o $(FSALIB)

${BIN}/fsabfs: fsabfs.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsabfs  fsabfs.o $(FSALIB)

${BIN}/fsafilter: fsafilter.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsafilter  fsafilter.o $(FSALIB)

${BIN}/fsaand: fsaand.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsaand  fsaand.o $(FSALIB)

${BIN}/fsaor: fsaor.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsaor  fsaor.o $(FSALIB)

${BIN}/fsanot: fsanot.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsanot  fsanot.o $(FSALIB)

${BIN}/fsaexists: fsaexists.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsaexists  fsaexists.o $(FSALIB)

${BIN}/fsacount: fsacount.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsacount  fsacount.o $(FSALIB)

${BIN}/fsaenumerate: fsaenumerate.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/fsaenumerate  fsaenumerate.o $(FSALIB)

${BIN}/gpwa: gpwa.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpwa  gpwa.o $(FSALIB) 

${BIN}/gpwab: gpwab.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpwab  gpwab.o $(FSALIB)

${BIN}/gpgenmult: gpgenmult.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpgenmult  gpgenmult.o $(FSALIB)

${BIN}/gpcheckmult: gpcheckmult.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpcheckmult  gpcheckmult.o $(FSALIB)

${BIN}/gpmult: gpmult.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpmult  gpmult.o $(FSALIB)

${BIN}/gpgenmult2: gpgenmult2.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpgenmult2  gpgenmult2.o $(FSALIB)

${BIN}/gpmult2: gpmult2.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpmult2  gpmult2.o $(FSALIB)

${BIN}/gpgenmult3: gpgenmult3.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpgenmult3  gpgenmult3.o $(FSALIB)

${BIN}/gpcomp: gpcomp.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpcomp  gpcomp.o $(FSALIB)

${BIN}/gpaxioms: gpaxioms.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpaxioms  gpaxioms.o $(FSALIB)

${BIN}/gpminkb: gpminkb.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpminkb  gpminkb.o $(FSALIB)

${BIN}/gpmakefsa: gpmakefsa.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/gpmakefsa  gpmakefsa.o $(FSALIB)

${BIN}/wordreduce: wordreduce.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/wordreduce  wordreduce.o $(FSALIB)

${BIN}/ppgap: ppgap.o $(FSALIB)
	$(CC) $(CFLAGS) -o ${BIN}/ppgap  ppgap.o $(FSALIB)

clean:
	$(RM) -f *.o core a.out *\~
