CC      = gcc
CFLAGS  = -O2
DEFINES = -DINT_SIZE=32
SRCDIR  = ../../src/leon/src

GAP_PATH=../..
PKG_PATH=..
SRCDISTFILE=guava

all:
	( test -d bin || mkdir bin; \
	test -d bin/`hostname` || mkdir bin/`hostname`; cd bin/`hostname`; \
	$(MAKE) -f ../../Makefile all2 CC="$(CC)" CFLAGS="$(CFLAGS)" )

all2:	  leonconv desauto wtdist

# rules to make the executables, just link them together
leonconv: leonconv.o
	$(CC) $(CFLAGS) -o leonconv leonconv.o
desauto:  desauto.o addsgen.o bitmanp.o cdesauto.o chbase.o cmatauto.o \
          code.o compcrep.o compsg.o copy.o cparstab.o cstborb.o cstrbas.o \
          errmesg.o essentia.o factor.o field.o inform.o matrix.o new.o \
          oldcopy.o optsvec.o orbit.o orbrefn.o partn.o permgrp.o permut.o \
          primes.o ptstbref.o randgrp.o randschr.o readdes.o readgrp.o \
          readper.o rprique.o storage.o token.o util.o
	$(CC) $(CFLAGS) -o desauto \
	  desauto.o addsgen.o bitmanp.o cdesauto.o chbase.o cmatauto.o \
	  code.o compcrep.o compsg.o copy.o cparstab.o cstborb.o cstrbas.o \
	  errmesg.o essentia.o factor.o field.o inform.o matrix.o new.o \
	  oldcopy.o optsvec.o orbit.o orbrefn.o partn.o permgrp.o permut.o \
	  primes.o ptstbref.o randgrp.o randschr.o readdes.o readgrp.o \
	  readper.o rprique.o storage.o token.o util.o
wtdist:   wtdist.o bitmanp.o code.o copy.o errmesg.o essentia.o factor.o \
          field.o new.o partn.o permgrp.o permut.o primes.o readdes.o \
          storage.o token.o util.o
	$(CC) $(CFLAGS) -o wtdist \
	  wtdist.o bitmanp.o code.o copy.o errmesg.o essentia.o factor.o \
	  field.o new.o partn.o permgrp.o permut.o primes.o readdes.o \
	  storage.o token.o util.o

# rules to make the .o files, just compile the .c file
# cannot use implicit rule, because .c files are in a different directory
leonconv.o: ../../src/leonconv.c
	$(CC) -c $(CFLAGS) -o leonconv.o -c ../../src/leonconv.c
addsgen.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/essentia.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/cstborb.h $(SRCDIR)/addsgen.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/addsgen.c 
bitmanp.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/bitmanp.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/bitmanp.c 
ccent.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/cparstab.h $(SRCDIR)/errmesg.h $(SRCDIR)/inform.h $(SRCDIR)/new.h $(SRCDIR)/orbrefn.h $(SRCDIR)/permut.h $(SRCDIR)/randgrp.h $(SRCDIR)/storage.h $(SRCDIR)/ccent.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/ccent.c 
ccommut.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/addsgen.h $(SRCDIR)/copy.h $(SRCDIR)/chbase.h $(SRCDIR)/new.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/storage.h $(SRCDIR)/ccommut.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/ccommut.c 
cdesauto.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/code.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/errmesg.h $(SRCDIR)/matrix.h $(SRCDIR)/storage.h $(SRCDIR)/cdesauto.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cdesauto.c 
cent.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/ccent.h $(SRCDIR)/errmesg.h $(SRCDIR)/permgrp.h $(SRCDIR)/readgrp.h $(SRCDIR)/readper.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/cent.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cent.c 
chbase.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/addsgen.h $(SRCDIR)/cstborb.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/factor.h $(SRCDIR)/new.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/randgrp.h $(SRCDIR)/storage.h $(SRCDIR)/repinimg.h $(SRCDIR)/settoinv.h $(SRCDIR)/chbase.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/chbase.c 
cinter.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/errmesg.h $(SRCDIR)/orbrefn.h $(SRCDIR)/cinter.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cinter.c 
cjrndper.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/code.h $(SRCDIR)/copy.h $(SRCDIR)/errmesg.h $(SRCDIR)/matrix.h $(SRCDIR)/new.h $(SRCDIR)/permut.h $(SRCDIR)/readdes.h $(SRCDIR)/randgrp.h $(SRCDIR)/readgrp.h $(SRCDIR)/readpar.h $(SRCDIR)/readper.h $(SRCDIR)/readpts.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/cjrndper.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cjrndper.c 
cmatauto.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/code.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/errmesg.h $(SRCDIR)/matrix.h $(SRCDIR)/storage.h $(SRCDIR)/cmatauto.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cmatauto.c 
code.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/storage.h $(SRCDIR)/code.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/code.c 
commut.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/ccommut.h $(SRCDIR)/errmesg.h $(SRCDIR)/factor.h $(SRCDIR)/permgrp.h $(SRCDIR)/readgrp.h $(SRCDIR)/readper.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/commut.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/commut.c 
compcrep.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/cputime.h $(SRCDIR)/chbase.h $(SRCDIR)/cstrbas.h $(SRCDIR)/errmesg.h $(SRCDIR)/inform.h $(SRCDIR)/new.h $(SRCDIR)/optsvec.h $(SRCDIR)/orbit.h $(SRCDIR)/orbrefn.h $(SRCDIR)/partn.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/rprique.h $(SRCDIR)/storage.h $(SRCDIR)/compcrep.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/compcrep.c 
compgrp.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/errmesg.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/readgrp.h $(SRCDIR)/util.h $(SRCDIR)/compgrp.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/compgrp.c 
compsg.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/cputime.h $(SRCDIR)/addsgen.h $(SRCDIR)/chbase.h $(SRCDIR)/copy.h $(SRCDIR)/cstrbas.h $(SRCDIR)/errmesg.h $(SRCDIR)/inform.h $(SRCDIR)/new.h $(SRCDIR)/optsvec.h $(SRCDIR)/orbit.h $(SRCDIR)/orbrefn.h $(SRCDIR)/partn.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/ptstbref.h $(SRCDIR)/rprique.h $(SRCDIR)/storage.h $(SRCDIR)/compsg.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/compsg.c 
copy.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/essentia.h $(SRCDIR)/storage.h $(SRCDIR)/copy.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/copy.c 
cparstab.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/errmesg.h $(SRCDIR)/orbrefn.h $(SRCDIR)/cparstab.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cparstab.c 
csetstab.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/errmesg.h $(SRCDIR)/orbrefn.h $(SRCDIR)/csetstab.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/csetstab.c 
cstborb.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/factor.h $(SRCDIR)/storage.h $(SRCDIR)/cstborb.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cstborb.c 
cstrbas.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/chbase.h $(SRCDIR)/cstborb.h $(SRCDIR)/inform.h $(SRCDIR)/new.h $(SRCDIR)/orbrefn.h $(SRCDIR)/permgrp.h $(SRCDIR)/ptstbref.h $(SRCDIR)/optsvec.h $(SRCDIR)/storage.h $(SRCDIR)/cstrbas.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cstrbas.c 
cuprstab.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/compcrep.h $(SRCDIR)/compsg.h $(SRCDIR)/errmesg.h $(SRCDIR)/orbrefn.h $(SRCDIR)/cuprstab.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/cuprstab.c 
desauto.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/cdesauto.h $(SRCDIR)/errmesg.h $(SRCDIR)/permgrp.h $(SRCDIR)/readdes.h $(SRCDIR)/readgrp.h $(SRCDIR)/readper.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/desauto.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/desauto.c 
errmesg.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/errmesg.c 
essentia.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/essentia.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/essentia.c 
factor.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/factor.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/factor.c 
field.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/storage.h $(SRCDIR)/field.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/field.c 
fndelt.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/oldcopy.h $(SRCDIR)/permut.h $(SRCDIR)/readgrp.h $(SRCDIR)/readper.h $(SRCDIR)/randgrp.h $(SRCDIR)/util.h $(SRCDIR)/fndelt.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/fndelt.c 
generate.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/enum.h $(SRCDIR)/storage.h $(SRCDIR)/cputime.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/readgrp.h $(SRCDIR)/randschr.h $(SRCDIR)/stcs.h $(SRCDIR)/util.h $(SRCDIR)/generate.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/generate.c 
inform.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/cputime.h $(SRCDIR)/readgrp.h $(SRCDIR)/inform.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/inform.c 
inter.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/cinter.h $(SRCDIR)/errmesg.h $(SRCDIR)/permgrp.h $(SRCDIR)/readgrp.h $(SRCDIR)/readper.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/inter.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/inter.c 
matrix.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/storage.h $(SRCDIR)/matrix.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/matrix.c 
new.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/partn.h $(SRCDIR)/storage.h $(SRCDIR)/new.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/new.c 
oldcopy.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/storage.h $(SRCDIR)/oldcopy.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/oldcopy.c 
optsvec.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/cstborb.h $(SRCDIR)/essentia.h $(SRCDIR)/new.h $(SRCDIR)/permut.h $(SRCDIR)/permgrp.h $(SRCDIR)/storage.h $(SRCDIR)/optsvec.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/optsvec.c 
orbdes.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/chbase.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/oldcopy.h $(SRCDIR)/permut.h $(SRCDIR)/readdes.h $(SRCDIR)/readgrp.h $(SRCDIR)/storage.h $(SRCDIR)/util.h $(SRCDIR)/orbdes.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/orbdes.c 
orbit.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/cstborb.h $(SRCDIR)/orbit.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/orbit.c 
orblist.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/addsgen.h $(SRCDIR)/chbase.h $(SRCDIR)/cstborb.h $(SRCDIR)/errmesg.h $(SRCDIR)/factor.h $(SRCDIR)/new.h $(SRCDIR)/randgrp.h $(SRCDIR)/readgrp.h $(SRCDIR)/readpar.h $(SRCDIR)/readpts.h $(SRCDIR)/storage.h $(SRCDIR)/util.h $(SRCDIR)/orblist.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/orblist.c 
orbrefn.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/partn.h $(SRCDIR)/storage.h $(SRCDIR)/orbrefn.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/orbrefn.c 
partn.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/storage.h $(SRCDIR)/permgrp.h $(SRCDIR)/partn.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/partn.c 
permgrp.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/copy.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/new.h $(SRCDIR)/permut.h $(SRCDIR)/storage.h $(SRCDIR)/permgrp.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/permgrp.c 
permut.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/factor.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/storage.h $(SRCDIR)/repimg.h $(SRCDIR)/repinimg.h $(SRCDIR)/settoinv.h $(SRCDIR)/enum.h $(SRCDIR)/permut.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/permut.c 
primes.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/primes.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/primes.c 
ptstbref.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/partn.h $(SRCDIR)/cstrbas.h $(SRCDIR)/errmesg.h $(SRCDIR)/ptstbref.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/ptstbref.c 
randgrp.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/new.h $(SRCDIR)/permut.h $(SRCDIR)/randgrp.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/randgrp.c 
randobj.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/errmesg.h $(SRCDIR)/randgrp.h $(SRCDIR)/readpar.h $(SRCDIR)/readpts.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/randobj.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/randobj.c 
randschr.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/addsgen.h $(SRCDIR)/cstborb.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/factor.h $(SRCDIR)/new.h $(SRCDIR)/oldcopy.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/randgrp.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/randschr.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/randschr.c 
readdes.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/code.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/token.h $(SRCDIR)/readdes.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/readdes.c 
readgrp.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/chbase.h $(SRCDIR)/cstborb.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/factor.h $(SRCDIR)/permut.h $(SRCDIR)/permgrp.h $(SRCDIR)/randschr.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/readgrp.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/readgrp.c 
readpar.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/errmesg.h $(SRCDIR)/readpar.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/readpar.c 
readper.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/readgrp.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/readper.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/readper.c 
readpts.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/errmesg.h $(SRCDIR)/readpts.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/readpts.c 
relator.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/enum.h $(SRCDIR)/errmesg.h $(SRCDIR)/new.h $(SRCDIR)/permut.h $(SRCDIR)/stcs.h $(SRCDIR)/storage.h $(SRCDIR)/relator.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/relator.c 
rprique.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/storage.h $(SRCDIR)/rprique.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/rprique.c 
setstab.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/cparstab.h $(SRCDIR)/csetstab.h $(SRCDIR)/cuprstab.h $(SRCDIR)/errmesg.h $(SRCDIR)/permgrp.h $(SRCDIR)/readgrp.h $(SRCDIR)/readpar.h $(SRCDIR)/readper.h $(SRCDIR)/readpts.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/setstab.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/setstab.c 
stcs.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/enum.h $(SRCDIR)/repimg.h $(SRCDIR)/addsgen.h $(SRCDIR)/cstborb.h $(SRCDIR)/errmesg.h $(SRCDIR)/essentia.h $(SRCDIR)/factor.h $(SRCDIR)/new.h $(SRCDIR)/oldcopy.h $(SRCDIR)/permgrp.h $(SRCDIR)/permut.h $(SRCDIR)/randschr.h $(SRCDIR)/relator.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/stcs.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/stcs.c 
storage.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/errmesg.h $(SRCDIR)/storage.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/storage.c 
token.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/errmesg.h $(SRCDIR)/token.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/token.c 
util.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/util.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/util.c 
wtdist.o: $(SRCDIR)/group.h $(SRCDIR)/extname.h $(SRCDIR)/groupio.h $(SRCDIR)/errmesg.h $(SRCDIR)/field.h $(SRCDIR)/readdes.h $(SRCDIR)/storage.h $(SRCDIR)/token.h $(SRCDIR)/util.h $(SRCDIR)/wt.h $(SRCDIR)/swt.h $(SRCDIR)/wtdist.c
	$(CC) -c $(CFLAGS) -I $(SRCDIR) $(DEFINES) $(SRCDIR)/wtdist.c 

# pseudo targets
clean:
	( cd bin/`hostname`; rm -f *.o )

distclean:
	( rm -rf bin/`hostname` )

# for GAP distribution
src_dist:
	@(cmp ${PKG_PATH}/guava/doc/guava.tex \
         ${GAP_PATH}/doc/guava.tex \
         || echo \
         "*** WARNING: current 'guava.tex' and 'doc/guava.tex' differ ***")
	@zoo ah ${SRCDISTFILE}.zoo \
                ${PKG_PATH}/guava/Makefile \
                ${PKG_PATH}/guava/doc/guava.tex \
                ${PKG_PATH}/guava/init.g \
                `find ${PKG_PATH}/guava/lib -name "*.g" -print` \
                `find ${PKG_PATH}/guava/tbl -name "*.g" -print` \
                `find ${PKG_PATH}/guava/src -print`
	@zoo PE ${SRCDISTFILE}.zoo

