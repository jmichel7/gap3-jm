##  shell to use to execute the commands.
SHELL     = /bin/sh

##  global options, used only for packing the GAP distribution
SRCDISTFILE = pkg

##  set this to point to the include directory containing 'gmp.h'.
INCDIRGMP = /usr/local/include

##  set this to point to the library directory containing 'libgmp.a'.
LIBDIRGMP = /usr/local/lib

#############################################################################
##
#T  Targets . . . . .  targets of the form <vendor>-<cpu>-<system>-<compiler>
##
##  Following are the targets that the user can specify.
##
unknown:
	@echo "usage: 'make <target>' where target is one of"
	@echo "'ibm-i386-386bsd-gcc'  for IBM PC under 386BSD with gcc"
	@echo "'ibm-i386-linux-gcc'   for IBM PC under Linux with gcc"
	@echo "'ibm-i386-msdos-djgpp'  for IBM PC under MSDOS with DJ GNU cc"
	@echo "'sun-sparc-sunos-gcc'  for SUN 4 under SunOs with gcc"
	@echo "'sun-sparc-sunos-cc'    for SUN 4 under SunOS with cc"
	@echo "'usg'                   for others under System V UNIX with cc"
	@echo "   additional C compiler and linker flags can be passed with"
	@echo "   'make <target> CFLAGS=<compiler-opts> LOPTS=<linker-opts>',"
	@echo "   i.e., 'make sun-sparc-sunos-gcc CFLAGS=-g LOPTS=-g.'"

gcc     :
	$(MAKE) CC="gcc" CFLAGS="-O2 -finline-functions -m32" all

cc      :
	$(MAKE) CC="cc"  CFLAGS="-O" all

dos 	:
	$(MAKE) CC="cross-gcc" CFLAGS="-O2 -finline-functions" all

all     : me qme

cleano  :
	( cd src ; rm -f int/*.o ve/*.o rat/*.o )
clean   : cleano
	( cd bin ; rm -f *)

me      :
	( cd src ; $(MAKE) PROGNAME="me" PROGDEF="-DME" OBJSDIR="ve" \
	OBJECTX="scrut.o" INCLUDEX="" MEFLAGS="$(MEFLAGS)" \
	CC="$(CC)" CFLAGS="$(CFLAGS)" INCDIRS="" \
	LD="$(CC)" LFLAGS="$(LFLAGS)" LIBDIRS="" LIBFILES="-lm" \
	BISON="bison" -f ../Makefile ../bin/me )

qme     :
	( cd src ; $(MAKE) PROGNAME="qme" PROGDEF="-DQME" OBJSDIR="rat" \
	OBJECTX="scrut.o" INCLUDEX="" MEFLAGS="$(MEFLAGS)" \
	CC="$(CC)" CFLAGS="$(CFLAGS)" INCDIRS="-I$(INCDIRGMP)" \
	LD="$(CC)" LFLAGS="$(LFLAGS)" LIBDIRS="-L$(LIBDIRGMP)" \
        LIBFILES="-lm -lgmp" \
	BISON="bison" -f ../Makefile ../bin/qme )

zme     :
	( cd src ; $(MAKE) PROGNAME="zme" PROGDEF="-DZME" OBJSDIR="int" \
	OBJECTX="lattice.o" INCLUDEX="../latt.h" MEFLAGS="$(MEFLAGS)" \
	CC="$(CC)" CFLAGS="$(CFLAGS)" INCDIRS="-I$(INCDIRGMP)" \
	LD="$(CC)" LFLAGS="$(LFLAGS)" LIBDIRS="-L$(LIBDIRGMP)" \
        LIBFILES="-lm -lgmp" \
	BISON="bison" -f ../Makefile ../bin/zme )

ibm-i386-linux-gcc:
	@$(MAKE) me  CC=gcc  CFLAGS="$(CFLAGS) -m32" MEFLAGS="-DSYS_IS_USG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SIG_T=void"
	@$(MAKE) qme CC=gcc  CFLAGS="$(CFLAGS) -m32" MEFLAGS="-DSYS_IS_USG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SIG_T=void"

# 'stdio.h' claims 'fgets' has second argument 'size_t' (ANSI says 'int')
# 'sys/ioctl.h' claims 'ioctl' has third argument '...' (UNIX says 'char*')
# 'sys/signal.h' claims 'kill' has first argument 'pid_t' (UNIX says 'int')
ibm-i386-386bsd-gcc:
	@$(MAKE) me  CC=gcc  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_BSD -DSYS_HAS_STDIO_PROTO -DSYS_HAS_IOCTL_PROTO -DSYS_HAS_SIGNAL_PROTO"
	@$(MAKE) qme CC=gcc  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_BSD -DSYS_HAS_STDIO_PROTO -DSYS_HAS_IOCTL_PROTO -DSYS_HAS_SIGNAL_PROTO"

# 'std.h' claims 'strncat' has third argument 'int' (ANSI says 'size_t')
# 'stdio.h' claims 'setbuf' has return type 'int' (UNIX says 'void')
# 'std.h' claims 'sbrk' has return type 'void*' (UNIX says 'char*')
ibm-i386-msdos-djgpp:
	@$(MAKE) me  CC=gcc  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_MSDOS_DJGPP -DSYS_HAS_STRING_PROTO -DSYS_HAS_STDIO_PROTO -DSYS_HAS_MISC_PROTO"
	@$(MAKE) qme CC=gcc  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_MSDOS_DJGPP -DSYS_HAS_STRING_PROTO -DSYS_HAS_STDIO_PROTO -DSYS_HAS_MISC_PROTO"

# 'sys/times.h' claims 'times' returns 'clock_t' (how shall it return -1?)
sun-sparc-sunos-gcc:
	@$(MAKE) me  CC=gcc  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_USG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SIGNAL_PROTO"
	@$(MAKE) qme CC=gcc  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_USG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SIGNAL_PROTO"

# 'sys/times.h' claims 'times' returns 'clock_t' (how shall it return -1?)
sun-sparc-sunos-cc:
	@$(MAKE) me  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_USG -DSYS_HAS_TIME_PROTO"
	@$(MAKE) qme CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_USG -DSYS_HAS_TIME_PROTO"

usg:
	@$(MAKE) me  CC=$(CC)  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_USG"
	@$(MAKE) qme CC=$(CC)  CFLAGS="$(CFLAGS)" MEFLAGS="-DSYS_IS_USG"

#############################################################################
##
##  Must be called with the following arguments:
##  PROGNAME=<name of program to build>
##  PROGDEF=<define for this program, e.g., -DME>
##  OBJECTX=<additional object>
##  INCLUDEX=<additional dependency>
##  CC=<C compiler to use>
##  CFLAGS=<options for the compiler, e.g., -O>
##  INCDIRS=<additional directories to search for include files -I...>
##  LD=<linker to use>
##  LFLAGS=<options for the linker>
##  LIBDIRS=<additional directories to search for library files -L...>
##  LIBFILES=<additional libraries to link -l...>
##

#############################################################################
##
##  Command to link the executable.
##
OBJECTS1 = $(OBJSDIR)/myalloc.o $(OBJSDIR)/global.o
OBJECTS2 = $(OBJSDIR)/pack.o    $(OBJSDIR)/allocs.o $(OBJSDIR)/out.o \
           $(OBJSDIR)/push.o    $(OBJSDIR)/comline.o $(OBJSDIR)/qinput.o \
           $(OBJSDIR)/memain.o  $(OBJSDIR)/me.o $(OBJSDIR)/coin.o \
           $(OBJSDIR)/vector.o
OBJECTS3 = $(OBJSDIR)/scanner.o $(OBJSDIR)/inputfun.o $(OBJSDIR)/input.tab.o
OBJECTS4 = $(OBJSDIR)/$(OBJECTX)
OBJECTS  = $(OBJECTS1) $(OBJECTS2) $(OBJECTS3) $(OBJECTS4)

../bin/$(PROGNAME) : $(OBJECTS)
	$(LD) $(LFLAGS) -m32 $(LIBDIRS) -o ../bin/$(PROGNAME) \
	      $(OBJECTS) $(LIBFILES)

#############################################################################
##
##  Commands to compile the sources and dependencies.
##
INCLUDE1 = global.h myalloc.h
INCLUDE2 = meint.h
INCLUDE3 = input.h input.tab.h

$(OBJSDIR)/myalloc.o   : myalloc.c   $(INCLUDE1)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c myalloc.c   -o $(OBJSDIR)/myalloc.o
$(OBJSDIR)/global.o    : global.c    $(INCLUDE1)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c global.c    -o $(OBJSDIR)/global.o

$(OBJSDIR)/pack.o      : pack.c      $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c pack.c      -o $(OBJSDIR)/pack.o
$(OBJSDIR)/allocs.o    : allocs.c    $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c allocs.c    -o $(OBJSDIR)/allocs.o
$(OBJSDIR)/out.o       : out.c       $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c out.c       -o $(OBJSDIR)/out.o
$(OBJSDIR)/push.o      : push.c      $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c push.c      -o $(OBJSDIR)/push.o
$(OBJSDIR)/comline.o   : comline.c   $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c comline.c   -o $(OBJSDIR)/comline.o
$(OBJSDIR)/qinput.o    : qinput.c    $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c qinput.c    -o $(OBJSDIR)/qinput.o
$(OBJSDIR)/memain.o    : memain.c    $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c memain.c    -o $(OBJSDIR)/memain.o
$(OBJSDIR)/me.o        : me.c        $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(MEFLAGS) $(PROGDEF) -I. $(INCDIRS) -c me.c        -o $(OBJSDIR)/me.o
$(OBJSDIR)/coin.o      : coin.c      $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c coin.c      -o $(OBJSDIR)/coin.o
$(OBJSDIR)/vector.o    : vector.c    $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c vector.c    -o $(OBJSDIR)/vector.o

$(OBJSDIR)/scanner.o   : scanner.c   $(INCLUDE1) $(INCLUDE2) $(INCLUDEX) $(INCLUDE3)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c scanner.c   -o $(OBJSDIR)/scanner.o
$(OBJSDIR)/inputfun.o     : inputfun.c     $(INCLUDE1) $(INCLUDE2) $(INCLUDEX) $(INCLUDE3)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c inputfun.c -o $(OBJSDIR)/inputfun.o

$(OBJSDIR)/input.tab.o : input.tab.c input.tab.h
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c input.tab.c -o $(OBJSDIR)/input.tab.o
input.tab.c input.tab.h : input.y
	$(BISON) -d input.y -o input.tab.c

$(OBJSDIR)/scrut.o     : scrut.c     $(INCLUDE1) $(INCLUDE2) $(INCLUDEX)
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c scrut.c     -o $(OBJSDIR)/scrut.o
$(OBJSDIR)/lattice.o   : lattice.c   $(INCLUDE1) $(INCLUDE2) $(INCLUDEX) ilatt.h
	$(CC) $(CFLAGS) $(PROGDEF) -I. $(INCDIRS) -c lattice.c   -o $(OBJSDIR)/lattice.o
