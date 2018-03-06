#############################################################################
##
#A  Makefile			GAP source		     Martin Schoenert
##
#Y  Copyright (C)  1994,  Lehrstuhl D fuer Mathematik,	RWTH Aachen,  Germany
##
##  This file contains the definitions and actions to compile and  link	 GAP.
##
##  To compile and link GAP you simply enter
##
##	make
##
##  To specify additional compiler and linker options enter
##
##	make  <target>	COPTS=<compiler-options>  LOPTS=<linker-options>
##
#H  Revision 2014/08/22 Jean Michel
#H  remove obsolete targets. Add 64bit target.
#H
SHELL = /bin/sh

ARCH = $(shell getconf LONG_BIT)
UNAME_S = $(shell uname -s)
LINUX=$(COPTS) -DSYS_HAS_SIGNAL_PROTO -DSYS_HAS_IOCTL_PROTO -DSYS_IS_LINUX

unknown:
	@echo "###########################################################"
	@echo "#  Building GAP3 for '$(UNAME_S)' and '$(ARCH) bit'"
	@echo "###########################################################"

ifeq ($(UNAME_S),Linux)
    ifeq ($(ARCH),64)
	@$(MAKE) system.o   CC=gcc  CFLAGS="$(LINUX) -O2 -DSYS_IS_64_BIT "
	@$(MAKE) gap CC=gcc CFLAGS="$(COPTS) -O -DSYS_IS_64_BIT " LOPTS="$(LOPTS) -static"
	mv gap ../bin/gap3
    else
	@$(MAKE) system.o   CC=gcc  CFLAGS="$(LINUX) -O2"
	@$(MAKE) gap CC=gcc CFLAGS="$(COPTS) -O " LOPTS="$(LOPTS) -static"
	mv gap ../bin/gap3
    endif
endif
ifeq ($(UNAME_S),Darwin)
    ifeq ($(ARCH),64)
		@$(MAKE) system.o   CC=gcc  CFLAGS="$(COPTS) -DSYS_IS_MACOSX -DSYS_IS_64_BIT -DARCH_INCLUDE -DSYS_HAS_IOCTL_PROTO"
		@$(MAKE) gap	    CC=gcc  CFLAGS="$(COPTS) -O2 -DSYS_IS_64_BIT "  LOPTS="$(LOPTS)"
		mv gap ../bin/gap3
    else
		@$(MAKE) system.o   CC=gcc  CFLAGS="$(COPTS) -DSYS_IS_MACOSX -DARCH_INCLUDE -DSYS_HAS_IOCTL_PROTO"
		@$(MAKE) gap	    CC=gcc  CFLAGS="$(COPTS) -O2 "  LOPTS="$(LOPTS)"
		mv gap ../bin/gap3
    endif
endif

	@echo "###########################################################"
	@echo "#  Building complete, you can run GAP3 from './../bin/gap3'"
	@echo "###########################################################"

#############################################################################
##
#T  Linkage . . . . . . . . . . . . . . . . . . . . . .	 commands to link GAP
##
##  Following are the commands to link GAP.  The targets below should not  be
##  used directly, they should only be addressed by the	 recursive  calls  of
##  'make' from above.
##
##  After  linking for  DOS we  prepend the  DOS extender 'GO32.EXE' to the
##  file,  so that we get an executable that works without further run-time
##  support.
##
OBJECTS = gap.o system.o gasman.o scanner.o idents.o read.o eval.o\
          integer.o rational.o cyclotom.o unknown.o finfield.o polynom.o\
          permutat.o word.o costab.o tietze.o agcollec.o aggroup.o pcpresen.o\
          list.o plist.o set.o vector.o vecffe.o range.o blister.o string.o\
          record.o statemen.o function.o coding.o

gap	  : $(OBJECTS)
	$(CC) $(LOPTS) -o gap $(OBJECTS)

gapdjg.exe: $(OBJECTS)
	$(CC) $(LOPTS) -o gapdjg.out $(OBJECTS) -lpc
	copy /b \djgpp\bin\go32.exe+gapdjg.out gapdjg.exe

#############################################################################
##
#T  Dependencies  . . . . . . . . .  dependencies between the packages of GAP
##
##  Following are the dependencies between the various packages of GAP.	 Note
##  that the dependencies '<package>.o: <package>.c' are  not  given,  'make'
##  adds them automagically.
##
INCBASC = system.h gasman.h scanner.h idents.h eval.h integer.h
INCFLDS = rational.h cyclotom.h unknown.h finfield.h polynom.h
INCGRPS = permutat.h word.h costab.h tietze.h agcollec.h aggroup.h pcpresen.h
INCLIST = list.h plist.h set.h vector.h vecffe.h range.h blister.h string.h
INCLANG = record.h statemen.h function.h coding.h

gap.o	  :	       $(INCBASC) list.h plist.h string.h read.h $(INCLANG)
system.o  : system.h
gasman.o  : gasman.h   system.h
scanner.o : scanner.h  system.h
idents.o  : idents.h   system.h gasman.h
read.o	  : read.h     $(INCBASC)
eval.o	  : eval.h     $(INCBASC) $(INCFLDS) $(INCGRPS) $(INCLIST) $(INCLANG)
integer.o : integer.h  $(INCBASC)
rational.o: rational.h $(INCBASC)
cyclotom.o: cyclotom.h $(INCBASC) list.h
unknown.o : unknown.h  $(INCBASC)
finfield.o: finfield.h $(INCBASC)
polynom.o : polynom.h  $(INCBASC) list.h plist.h vector.h vecffe.h finfield.h
permutat.o: permutat.h $(INCBASC) list.h
word.o	  : word.h     $(INCBASC) list.h agcollec.h
costab.o  : costab.h   $(INCBASC) list.h
tietze.o  : tietze.h   $(INCBASC) list.h word.h agcollec.h
agcollec.o: agcollec.h $(INCBASC) list.h word.h aggroup.h
aggroup.o : aggroup.h  $(INCBASC) list.h plist.h word.h agcollec.h finfield.h
pcpresen.o: pcpresen.h $(INCBASC) list.h plist.h word.h agcollec.h aggroup.h
list.o	  : list.h     $(INCBASC) set.h permutat.h record.h
plist.o	  : plist.h    $(INCBASC) list.h range.h record.h
set.o	  : set.h      $(INCBASC) list.h plist.h range.h
vector.o  : vector.h   $(INCBASC) list.h plist.h range.h
vecffe.o  : vecffe.h   $(INCBASC) list.h plist.h range.h
range.o	  : range.h    $(INCBASC) list.h plist.h
blister.o : blister.h  $(INCBASC) list.h plist.h range.h set.h
string.o  : string.h   $(INCBASC) list.h plist.h range.h
record.o  : record.h   $(INCBASC) plist.h string.h
statemen.o: statemen.h $(INCBASC) list.h
function.o: function.h $(INCBASC) plist.h record.h statemen.h
coding.o  : coding.h   $(INCBASC) list.h plist.h blister.h finfield.h vecffe.h

#############################################################################
##
#T  Miscellaneous  . . . . . . . . . . . . . . . . . targets to clean up, etc.
##
clean:
	rm -f *.o gap gapdjg.out gapdjg.exe
