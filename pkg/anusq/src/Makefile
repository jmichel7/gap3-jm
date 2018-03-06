#############################################################################
##
#A  Makefile                    ANU SQ                      Alice C Niemeyer
##
#Y  Copyright 1994,  Department of Mathematics,   ANU,   ACT 0200,  Australia
##
VRSN="1.0, March 1994"
##
##  This file contains the definitions and dependencies in order to  compile 
##  the ANU SQ.
##  
##  Please read the file README to obtain information on which variables to
##  change and which settings are necessary.
##

VEPATH = "$(GAP_DIR)/pkg/ve/bin/"

# CFLAGS = -O -DSYS_BSD -DCHAT=1 -DVERSION='$(VRSN)' -DVEP='$(VEPATH)'
CFLAGS = -O -DSYS_BSD -DCHAT=1 -DVERSION='$(VRSN)' -DVEP='""'

sq : mem.o modmem.o presentation.o pcp.o storepc.o arith.o sq.o \
	consis.o  pcparith.o updatepres.o veinter.o eapquot.o modifypres.o \
        system.o Makefile
	gcc -o Sq $(CFLAGS)  \
	mem.o system.o presentation.o modmem.o \
	pcp.o storepc.o arith.o sq.o consis.o pcparith.o updatepres.o \
        veinter.o eapquot.o modifypres.o
	mv Sq ../bin

cleano:
	-rm *.o
clean: cleano
	-rm ../bin/Sq
system.o:       system.c
arith.o:        arith.h arith.c modmem.h pres.h sq.h 
mem.o:          mem.c pres.h
presentation.o: pres.h presentation.c
modmem.o:       modmem.c pres.h sq.h pcparith.h
sq.o:           sq.c sys.h pres.h pres.h sq.h arith.h modmem.h pcparith.h \
		veinter.h eapquot.h updatepres.h
storepc.o:      storepc.c pres.h sq.h modmem.h 
pcp.o:          pcp.c pres.h  sq.h
consis.o:       consis.c pcparith.h arith.h sq.h modmem.h
pcparith.o:     pcparith.c pres.h sq.h modmem.h  arith.h
eapquot.o:      eapquot.c pres.h pcparith.h modmem.h sq.h
veinter.o:      veinter.c modmem.h sq.h pres.h
updatepres.o:   updatepres.c pcparith.h veinter.h modmem.h pres.h sq.h
modifypres.o:   modifypres.c modmem.h sq.h pres.h
