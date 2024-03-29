# Generated automatically from Makefile.in by configure.
# 
# Makefile for SISYPHOS project
# Generated automatically from Makefile.in by configure.

# SHELL = /bin/sh

#
# global options
#
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin

GCOPTS = -Dlint -DANSI -DYY_INTERACTIVE -DHAVE_CONFIG_H 
SOPTS  = -O2 -fomit-frame-pointer -DUNIX -m32

CC  = gcc
CFLAGS = -O2 -Wall -fomit-frame-pointer -DUNIX -m32
CPP	= gcc -E
YACC = bison -y
LEX = flex
#LIBS =  -lfl
LIBS = 
LOPTS = -m32

generic-unix:
	rm -f compdate.h
	@$(MAKE) sisgram.tab.o CFLAGS="$(SOPTS) $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CFLAGS="$(SOPTS) $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) sis           CFLAGS="$(CFLAGS) $(GCOPTS)" LOPTS="$(LOPTS)"

no-flex:
	rm -f compdate.h
	cp parser.c sisgram.tab.c
	cp parser.h sisgram.tab.h
	cp scanner.c lex.yy.c
	touch sisgram.tab.h
	@$(MAKE) sisgram.tab.o CFLAGS="$(SOPTS) $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CFLAGS="$(SOPTS) $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) sis           CFLAGS="$(CFLAGS) $(GCOPTS)" LOPTS="$(LOPTS)"

unknown:
	@echo "usage: 'make <target>'  where target is one of"
	@echo "'hp700-hpux-cci'     for HP 9000/7xx under HP-UX with cc - "
	@echo "                         generate version for profile dependent optimization"
	@echo "'hp700-hpux-ccp'     for HP 9000/7xx under HP-UX with cc - "
	@echo "                         relink with profile dependent optimization"
	@echo "'ibmpc-emx-gcc2'     for IBM PCs under DOS or OS/2 2.0 with emx"
	@echo "'generic-unix'       for most UNIX machines"

dev:
	rm -f compdate.h
	@$(MAKE) sisgram.tab.o CC=gcc CFLAGS="-g -ansi -Wall -pedantic \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CC=gcc CFLAGS="-g -ansi -Wall -pedantic \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)" 
	@$(MAKE) sis      CC=gcc CFLAGS="-g -Wall -ansi -Wall -pedantic \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)"

prof:
	rm -f compdate.h
	@$(MAKE) sisgram.tab.o CC=gcc CFLAGS="-g -pg -ansi -Wall -pedantic \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CC=gcc CFLAGS="-g -pg -ansi -Wall -pedantic \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)" 
	@$(MAKE) sis      CC=gcc CFLAGS="-g -pg -Wall -ansi -Wall -pedantic \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)"

hp700-hpux-cci:
	rm -f compdate.h
	@$(MAKE) sisgram.tab.o CC=cc  CFLAGS="-Aa -O -DUNIX $(GCOPTS)" \
		LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CC=cc  CFLAGS="-Aa  -DUNIX $(GCOPTS)" \
		LOPTS="$(LOPTS)"
	@$(MAKE) sis      CC=cc  CFLAGS="-Aa +I -O -DUNIX $(GCOPTS)" \
		LOPTS="+I $(LOPTS)"

hp700-hpux-ccp:
	rm -f compdate.h
	@$(MAKE) sisgram.tab.o CC=cc  CFLAGS="-Aa -O -DUNIX $(GCOPTS)" \
		LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CC=cc  CFLAGS="-Aa  -DUNIX $(GCOPTS)" \
		LOPTS="$(LOPTS)"
	@$(MAKE) sis      CC=cc  CFLAGS="-Aa +I -O -DUNIX $(GCOPTS)" \
		LOPTS="+P $(LOPTS)"

ibmpc-emx-gcc2:
	rm -f compdate.h
	@$(MAKE) sisgram.tab.o CC=emx-gcc  CFLAGS="-O2 -fomit-frame-pointer \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) lex.yy.o      CC=emx-gcc  CFLAGS="-O2 -fomit-frame-pointer \
		-DUNIX $(GCOPTS)" LOPTS="$(LOPTS)"
	@$(MAKE) sis.exe  CC=emx-gcc  CFLAGS="-O2 -Wall -fomit-frame-pointer \
		-DUNIX $(GCOPTS)"  LOPTS="$(LOPTS)"


OBFILES = pc.o parsesup.o hgroup.o matrix.o obstruct.o \
	lcontrol.o obslif.o lex.yy.o conju.o grpring.o \
	gl.o inout.o storage.o farith.o solve.o stdlib.o \
	error.o hash.o hashfunc.o symtab.o lie.o \
	space.o aut.o autgroup.o sisgram.tab.o version.o special.o aobstruc.o \
	aggroup.o cohomol.o gmodule.o sgaut.o agsupp.o \
	graut.o scansup.o dispatch.o siscode.o

sis : $(OBFILES)
	$(CC) $(LOPTS) -o sis $(OBFILES)  $(LIBS)
	mv sis ../bin/sis

sis.exe : $(OBFILES)
	$(CC) $(LOPTS) -o sis.exe $(OBFILES)
	emxbind -u $(EMX) sis.exe

sisgram.tab.c: sisgram.y
	$(YACC) -vdl sisgram.y
	mv y.tab.c sisgram.tab.c
	mv y.tab.h sisgram.tab.h
	   

sisgram.tab.h: sisgram.y
	$(YACC) -vdl sisgram.y
	mv y.tab.c sisgram.tab.c
	mv y.tab.h sisgram.tab.h

lex.yy.c: sisgram.tab.h sisscan.l 
	$(LEX) -L sisscan.l

compdate.h:
	echo "char date_string[] = \"`date`\";" >compdate.h

cleano:
	rm -f *.o *.output compdate.h lex.yy.c sisgram.tab.*

#
# dependencies
#
aggroup.o:     aggroup.c aglobals.h pc.h aggroup.h fdecla.h storage.h error.h
agsupp.o:      agsupp.c aglobals.h pc.h aggroup.h fdecla.h storage.h error.h \
               solve.h
aobstruc.o:    aobstruc.c aglobals.h fdecla.h pc.h hgroup.h grpring.h \
               storage.h solve.h
aut.o:         aut.c aglobals.h fdecla.h pc.h hgroup.h storage.h error.h \
               aut.h solve.h
autgroup.o:    autgroup.c aglobals.h fdecla.h pc.h hgroup.h storage.h \
               error.h aut.h solve.h
cohomol.o:     cohomol.c aglobals.h fdecla.h pc.h aut.h hgroup.h grpring.h \
               storage.h solve.h gmodule.h cohomol.h
conju.o:       conju.c aglobals.h fdecla.h pc.h grpring.h hgroup.h storage.h \
               solve.h
dispatch.o:    dispatch.c aglobals.h hash.h graut.h aut.h dispatch.h \
               parsesup.h grpring.h pc.h hgroup.h storage.h symtab.h error.h \
               gmodule.h cohomol.h
error.o:       error.c aglobals.h error.h
farith.o:      farith.c aglobals.h solve.h
gl.o:          gl.c aglobals.h fdecla.h pc.h storage.h solve.h
gmodule.o:     gmodule.c aglobals.h fdecla.h pc.h hgroup.h grpring.h \
               storage.h solve.h gmodule.h
graut.o:       graut.c config.h aglobals.h fdecla.h pc.h hgroup.h grpring.h \
               storage.h error.h aut.h graut.h solve.h
grpring.o:     grpring.c aglobals.h fdecla.h pc.h storage.h error.h
hash.o:        hash.c aglobals.h hash.h storage.h
hashfunc.o:    hashfunc.c debug.h hash.h
hgroup.o:      hgroup.c aglobals.h storage.h hgroup.h
inout.o:       inout.c aglobals.h fdecla.h pc.h hgroup.h storage.h aut.h
lcontrol.o:    lcontrol.c aglobals.h fdecla.h pc.h hgroup.h grpring.h \
               storage.h solve.h
lex.yy.o:      lex.yy.c config.h aglobals.h error.h dispatch.h graut.h \
               aut.h sisgram.tab.h parsesup.h
lie.o:         lie.c aglobals.h fdecla.h grpring.h pc.h storage.h error.h \
               solve.h
matrix.o:      matrix.c aglobals.h fdecla.h storage.h solve.h
obslif.o:      obslif.c aglobals.h fdecla.h pc.h hgroup.h grpring.h \
               storage.h aut.h solve.h
obstruct.o:    obstruct.c aglobals.h fdecla.h pc.h hgroup.h grpring.h \
               storage.h solve.h
parsesup.o:    parsesup.c aglobals.h fdecla.h graut.h aut.h pc.h \
               parsesup.h grpring.h gmodule.h cohomol.h storage.h siscode.h \
               error.h
pc.o:          pc.c aglobals.h graut.h aut.h pc.h fdecla.h storage.h \
               parsesup.h error.h
scansup.o:     scansup.c config.h aglobals.h error.h
sgaut.o:       sgaut.c aglobals.h pc.h aggroup.h hgroup.h fdecla.h storage.h \
               error.h solve.h aut.h
siscode.o:     siscode.c aglobals.h fdecla.h graut.h aut.h pc.h parsesup.h \
               grpring.h gmodule.h cohomol.h storage.h siscode.h error.h
sisgram.tab.o: sisgram.tab.c config.h aglobals.h fdecla.h pc.h graut.h \
               aggroup.h grpring.h hgroup.h symtab.h aut.h parsesup.h \
               storage.h error.h solve.h
solve.o:       solve.c aglobals.h fdecla.h storage.h solve.h
space.o:       space.c aglobals.h fdecla.h grpring.h pc.h storage.h solve.h
special.o:     special.c aglobals.h storage.h pc.h fdecla.h grpring.h
storage.o:     storage.c aglobals.h fdecla.h error.h storage.h
symtab.o:      symtab.c aglobals.h hash.h symtab.h graut.h
version.o:     version.c patchlev.h compdate.h
