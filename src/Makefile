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

#############################################################################
##
#T  Targets . . . . .  targets of the form <vendor>-<cpu>-<system>-<compiler>
##
##  Following are the targets that the	user can specify.  Each target	calls
##  'make' recursively twice.
##
##  The first time it calls
##
##	@$(MAKE) <system-file>.o  CFLAGS="$(COPTS) <copts>"
##
##  where <system-file> is either 'system' or the name of another system file
##  and where <copts> are the options to  pass to  pass	 to  the  C compiler,
##  e.g., 'SYS_IS_BSD' or 'SYS_IS_WINDOWS'.
##
##  The second time it calls
##
##	@$(MAKE) gap  CC=<c-compiler> CFLAGS="$(COPTS) <opts>" \
##		      LOPTS="$(LOPTS)"
##
##  where <system-file> is the same file as above, <c-compiler> is  the	 name
##  of the C compiler to use ($(CC) is used if this is not present),  <copts>
##  are	 the options to pass ot	 the  C compiler, e.g., '-O2'.
##
unknown:
	@echo "usage: 'make <target>' where target is one of"
	@echo "'linux'   for linux with gcc (64 bit)"
	@echo "'linux32' for linux with gcc in 32bit mode"
	@echo "'macosx'  for mac os x with gcc in 32 bit mode"
	@echo "'macosx64'  for mac os x with gcc (64 bit)"
	@echo "'bsd'     for bsd with gcc"
	@echo "'solaris CC=gcc' for Solaris with GNU cc"
	@echo "'solaris CC=cc'  for Solaris with cc"
	@echo "'sunos CC=gcc'   for Sun OS with GNU cc"
	@echo "'sunos CC=cc'    for Sun OS with cc"
	@echo "	"
	@echo "	  additional C compiler and linker flags can be passed with"
	@echo "	  'make <target> COPTS=<compiler-opts> LOPTS=<linker-opts>',"
	@echo "	  for instance, 'make linux COPTS=-g LOPTS=-g for debug'"

LINUX=$(COPTS) -DSYS_HAS_SIGNAL_PROTO -DSYS_HAS_IOCTL_PROTO -DSYS_IS_LINUX

linux32:
	@$(MAKE) system.o CC=gcc CFLAGS="$(LINUX) -O2 -m32"
	@$(MAKE) gap CC=gcc CFLAGS="$(COPTS) -O2 -m32" LOPTS="$(LOPTS) -m32 -static"
	strip gap
	mv gap ../bin/gap.linux32

linux:
	@$(MAKE) system.o CC=gcc  CFLAGS="$(LINUX) -O2 -DSYS_IS_64_BIT -Wall"
	@$(MAKE) gap  CC=gcc  CFLAGS="$(LINUX) -O2 -DSYS_IS_64_BIT -Wall" LOPTS="$(LOPTS)  -static"
	strip gap
	mv gap ../bin/gap.linux

linuxd:
	@$(MAKE) system.o CC=gcc  CFLAGS="$(LINUX) -DSYS_IS_64_BIT -g"
	@$(MAKE) gap CC=gcc CFLAGS="$(COPTS) -DSYS_IS_64_BIT -g" LOPTS="$(LOPTS)"

# sbrk doesn't work in macosx so we need to use vm_allocate
macosx:
	@$(MAKE) system.o   CC=gcc  CFLAGS="$(COPTS) -m32 -DSYS_IS_MACOSX -DARCH_INCLUDE -DSYS_HAS_IOCTL_PROTO"
	@$(MAKE) gap  CC=gcc  CFLAGS="$(COPTS) -m32 -O"  LOPTS="$(LOPTS) -m32"
	strip gap
	mv gap ../bin/gap.mac

# sbrk doesn't work in macosx so we need to use vm_allocate
macosx64:
	@$(MAKE) system.o   CC=gcc  CFLAGS="$(COPTS) -O -DSYS_IS_64BIT -DSYS_IS_MACOSX -DARCH_INCLUDE -DSYS_HAS_IOCTL_PROTO"
	@$(MAKE) gap CC=gcc CFLAGS="$(COPTS) -DSYS_IS_64BIT -O" LOPTS="$(LOPTS)"
	strip gap
	mv gap ../bin/gap.mac

x86-dos-djgpp:
	@$(MAKE) system.o   CC=gcc  CFLAGS="$(COPTS) -DSYS_IS_WINDOWS"
	@$(MAKE) gapdjg.exe CC=gcc  CFLAGS="$(COPTS) -O1"  LOPTS="$(LOPTS)"

x86-dos-djgppcross:
	@$(MAKE) system.o   CC=dos-gcc	CFLAGS="$(COPTS) -DSYS_IS_WINDOWS"
	@$(MAKE) gapdjg.exe CC=dos-gcc	CFLAGS="$(COPTS) -O1"  LOPTS="$(LOPTS)"

win64:
	@$(MAKE) system.o   CC=x86_64-w64-mingw32-gcc-posix CFLAGS="$(COPTS) -DSYS_IS_64BIT -DSYS_IS_WINDOWS -DSYS_HAS_NO_GETRUSAGE"
	@$(MAKE) gap CC=x86_64-w64-mingw32-gcc CFLAGS="$(COPTS) -O -DSYS_IS_64BIT" LOPTS="$(LOPTS) -static"

solaris:
	@$(MAKE) system.o CC=$(CC)  CFLAGS="$(COPTS) -O -DSYS_IS_SOLARIS -DSOLARIS2"
	@$(MAKE) gap	  CC=$(CC) CFLAGS="$(COPTS) -O2" LOPTS="$(LOPTS)"

# 'sys/times.h' claims 'times' returns 'clock_t' (how shall it return -1?)
sunos:
	@$(MAKE) system.o   CC=$(CC)  CFLAGS="$(COPTS) -DSYS_IS_SUNOS -DSYS_HAS_SIGNAL_PROTO"
	@$(MAKE) gap	    CC=$(CC)  CFLAGS="$(COPTS) -O2"  LOPTS="$(LOPTS)"

bsd:
	@$(MAKE) system.o   CC=gcc  CFLAGS="$(COPTS) -DSYS_IS_BSD"
	@$(MAKE) gap	    CC=gcc  CFLAGS="$(COPTS) -O"  LOPTS="$(LOPTS)"

###########################################################################
##  Following are the commands to link GAP. The targets below should not be
##  used  directly, they should only be addressed by the recursive calls of
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

gap: $(OBJECTS)
	$(CC) $(LOPTS) -o gap $(OBJECTS) -lrt

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
