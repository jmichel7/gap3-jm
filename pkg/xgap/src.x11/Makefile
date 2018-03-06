#############################################################################
##
#A  Makefile                	XGAP source                      Frank Celler
##
#H  @(#)$Id: Makefile,v 1.1.1.1 1996/12/11 12:40:06 werner Exp $
##
#Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: Makefile,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:06  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.18  1995/09/25  19:50:09  fceller
#H  removed dependency of Autofile on configure
#H
#H  Revision 1.17  1995/08/08  17:55:41  fceller
#H  added missing 'hp-hppa1.1-hpux-cc' target
#H
#H  Revision 1.16  1995/08/08  17:17:44  fceller
#H  changed makefile to use autoconf
#H
#H  Revision 1.15  1995/05/03  09:33:50  fceller
#H  added 'SYS_HAS_BROKEN_STDLIB'
#H
#H  Revision 1.14  1994/06/11  13:51:44  fceller
#H  added first (incomplete) SunOS 5.3 support
#H
#H  Revision 1.13  1994/06/11  12:44:48  fceller
#H  added "xgap.rdb" to source distribution
#H  added IRIX support
#H
#H  Revision 1.12  1994/06/06  08:55:19  fceller
#H  added FreeBSD
#H
#H  Revision 1.11  1993/12/23  08:43:54  fceller
#H  added 'ibm-i386-linux-gcc2'
#H
#H  Revision 1.11  1993/12/23  08:43:54  fceller
#H  added 'ibm-i386-linux-gcc2'
#H
#H  Revision 1.10  1993/10/21  17:02:27  fceller
#H  added first hp support
#H  
#H  Revision 1.9  1993/10/21  08:48:50  fceller
#H  fixed missing extension in 'bin_dist'
#H  
#H  Revision 1.7  1993/10/20  12:52:02  fceller
#H  added sun with gcc
#H
#H  Revision 1.6  1993/10/20  12:40:42  fceller
#H  changed exe names
#H
#H  Revision 1.5  1993/10/19  14:24:33  fceller
#H  fixed PKG_PATH
#H
#H  Revision 1.4  1993/10/18  15:45:34  fceller
#H  added 'xgap.sh' to backup
#H
#H  Revision 1.3  1993/10/18  11:13:11  fceller
#H  added 'src_dist', 'bin_dist'
#H
#H  Revision 1.2  1993/10/18  11:04:20  fceller
#H  added 'backup'
#H
#H  Revision 1.1  1993/08/18  11:11:55  fceller
#H  Initial revision
##


## <SHELL> should point to sensible shell ###################################
SHELL=/bin/sh


## <X11INC> points to the directory containing the "X11/X.h" files ##########
X11INC=/usr/include


## <X11LIB> points to the directory containing the X11 libraries ############
X11LIB=/usr/lib


## ALL: print a short help message ##########################################
all:
	@echo ; cat HELP ; echo


## MACHINES: a list of all supported machines ###############################
machines:
	@echo ; cat TARGETS ; echo


## final linkage script #####################################################
OBJS=xcmds.o utils.o gapgraph.o gaptext.o pty.o popdial.o xgap.o selfile.o
xgap: $(OBJS)
	$(CC) -o xgap $(OBJS) $(LDFLAGS) -lXaw -lXmu -lXt -lXext -lX11


## clean all "*.o" files, cores and executable ##############################
clean:
	rm -rf $(OBJS) Autofile config.* xgap core core.xgap *~


## AUTO: auto configure the system using GNU's autoconf package #############
auto: Autofile
	make -f Autofile SHELL=$(SHELL)

Autofile: Autofile.in # configure
	$(SHELL) ./configure

configure: configure.in
	autoconf
	@mv configure configure.auto
	sed -e 's|%\\\[%|%[[]%|' < configure.auto > configure
	@rm -f configure.auto config.cache
	@chmod 755 configure


## DEC ALPHA running OSF1, CC and X11R5 #####################################
dec-alpha-osf1-cc: .touchme
	make xgap SHELL="$(SHELL)" CC=cc \
                CFLAGS="$(COPTS) -I$(X11INC) -O \
		-DSYS_HAS_UNISTD -DSYS_HAS_STDLIB -D_BSD \
		-DSYS_HAS_UNION_WAIT -DSYS_HAS_STDARG" \
		LDFLAGS="$(LOPTS) -L$(X11LIB) -ldnet_stub"


.touchme:
	touch *
	touch .touchme


## Ultrix 4.2 running on a DecStation with GCC 2.5.x and X11R5 ##############
dec-mips-ultrix-gcc2: dec-mips-ultrix-gcc
dec-mips-ultrix-gcc:
	make xgap SHELL="$(SHELL)" CC=gcc \
		CFLAGS="$(COPTS) -Wall -I$(X11INC) -O2 \
		-DSYS_HAS_SIGNAL -DSYS_HAS_SYS_SIGNAL -DSYS_HAS_STDLIB \
		-DSYS_HAS_UNISTD -DSYS_HAS_UNION_WAIT -DSYS_HAS_STDARG" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## Ultrix 4.2 running on a DecStation with CC and X11R5 #####################
dec-mips-ultrix-cc:
	make xgap SHELL="$(SHELL)" CC=cc \
		CFLAGS="$(COPTS) -I$(X11INC) -O \
		-DSYS_HAS_SIGNAL -DSYS_HAS_SYS_SIGNAL -DSYS_HAS_STDLIB \
		-DSYS_HAS_UNISTD -DSYS_HAS_UNION_WAIT -DSYS_HAS_STDARG" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## HP 9000/700 running HP-UX A.09.01, CC and X11R5 ##########################
hp-hppa1.1-hpux-cc:
	make xgap SHELL="$(SHELL)" CC=cc \
		CFLAGS="$(COPTS) -I$(X11INC) -O 
		-DSYS_HAS_PID_T -DSYS_HAS_STDLIB -DSYS_HAS_SIGNAL \
		-DSYS_HAS_TCSETAW -DSYS_HAS_TERMIO -DSYS_HAS_TIME_PROTO" \
                LDFLAGS="$(LOPTS) -L$(X11LIB)"


## FreeBSD 1.1 running on a IBM-PC, GCC 2.5.x and X11R5 #####################
## FreeBSD 2.0 running on a IBM-PC, GCC 2.6.x and X11R6 #####################
ibm-i386-386bsd-cc:   ibm-i386-386bsd-gcc
ibm-i386-386bsd-gcc2: ibm-i386-386bsd-gcc
ibm-i386-386bsd-gcc:
	make xgap SHELL="$(SHELL)" CC=gcc \
		CFLAGS="$(COPTS) -Wall -I$(X11INC) -O2 \
		-DSYS_HAS_STDLIB -DSYS_HAS_UNISTD -DSYS_HAS_PID_T \
		-DSYS_HAS_STDARG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SYS_TIME \
		-DSYS_HAS_SIGNAL" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## Linux SLS 1.05 running on a IBM-PC, GCC 2.5.x and X11R5 ##################
ibm-i386-linux-gcc2: ibm-i386-linux-gcc
ibm-i386-linux-gcc:
	make xgap SHELL="$(SHELL)" CC=gcc \
		CFLAGS="$(COPTS) -Wall -I$(X11INC) -O2 \
		-DSYS_HAS_STDLIB -DSYS_HAS_UNISTD -DSYS_HAS_PID_T \
		-DSYS_HAS_STDARG -DSYS_HAS_TCSETAW -DSYS_HAS_SIGNAL \
		-DSYS_HAS_TERMIO" \
		LDFLAGS="$(LOPTS) -L$(X11LIB) -lbsd"


## NeXTSTEP 3.0 running on a M68040 with GCC 2.5.x and MouseX ###############
next-m68k-mach-gcc2: next-m68k-mach-gcc
next-m68k-mach-gcc:
	make xgap SHELL="$(SHELL)" CC=gcc \
		CFLAGS="$(COPTS) -Wall -I$(X11INC) -O2 \
		-DSYS_HAS_LIBC -DSYS_HAS_UNION_WAIT -DSYS_HAS_EXECV_CCHARPP \
		-DSYS_HAS_STDARG -DSYS_HAS_TIME_PROTO" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## NeXTSTEP 3.0 running on a M68040 with CC 1.93 and MouseX #################
## NeXTSTEP 3.3 running on a M68040 with CC 2.5.8 and MouseX ################
next-m68k-mach-cc:
	make xgap SHELL="$(SHELL)" CC=cc \
		CFLAGS="$(COPTS) -Wall -I$(X11INC) -O2 \
		-DSYS_HAS_LIBC -DSYS_HAS_UNION_WAIT -DSYS_HAS_EXECV_CCHARPP \
		-DSYS_HAS_STDARG -DSYS_HAS_TIME_PROTO" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## IRIX System V.4 running IRIX Release 5.3, GCC 2.7.x and X11R5 ############
sgi-mips-irix-gcc2: sgi-mips-irix-gcc
sgi-mips-irix-gcc:
	make xgap SHELL="$(SHELL)" CC=gcc \
		CFLAGS="$(COPTS) -I$(X11INC) -O \
		-DSYS_HAS_UNISTD -DSYS_HAS_STDLIB -DSYS_HAS_GETPTY \
		-DSYS_HAS_TERMIO -DSYS_HAS_TCSETAW -DSYS_HAS_TIME_PROTO \
		-DSYS_HAS_PID_T" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## IRIX System V.4 running IRIX Release 5.3, CC and X11R5 ###################
sgi-mips-irix-cc:
	make xgap SHELL="$(SHELL)" CC=cc \
		CFLAGS="$(COPTS) -I$(X11INC) -O \
		-DSYS_HAS_STDLIB -DSYS_HAS_GETPTY \
		-DSYS_HAS_TERMIO -DSYS_HAS_TCSETAW -DSYS_HAS_TIME_PROTO \
		-DSYS_HAS_PID_T" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## Sun running SunOS 4.1.2, GCC 2.6.x and X11R5 #############################
sun-sparc-sunos-gcc2: sun-sparc-sunos-gcc
sun-sparc-sunos-gcc:
	make xgap SHELL="$(SHELL)" CC=gcc \
                CFLAGS="$(COPTS) -I$(X11INC) -O \
		-DSYS_HAS_UNISTD -DSYS_HAS_UNION_WAIT -DSYS_HAS_SYS_TIME" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## Sun running SunOS 4.1.2, CC and X11R5 ####################################
sun-sparc-sunos-cc:
	make xgap SHELL="$(SHELL)" CC=cc \
	        CFLAGS="$(COPTS) -I$(X11INC) -O \
                -DSYS_HAS_UNISTD -DSYS_HAS_STDLIB -DSYS_HAS_UNION_WAIT \
                -DSYS_HAS_SYS_TIME" \
		LDFLAGS="$(LOPTS) -L$(X11LIB)"


## include files dependencies ###############################################
gapgraph.o: gapgraph.c utils.h gapgraph.h
gaptext.o:  gaptext.c  utils.h gaptext.h
popdial.o:  popdial.c  utils.h popdial.h
selfile.o:             utils.h

pty.o:      pty.c      utils.h pty.h xcmds.h xgap.h gaptext.h
xcmds.o:    xcmds.c    utils.h pty.h xcmds.h xgap.h gapgraph.h popdial.h
xgap.o:     xgap.c     utils.h pty.h xcmds.h xgap.h gaptext.h

utils.o:    utils.c    utils.h


