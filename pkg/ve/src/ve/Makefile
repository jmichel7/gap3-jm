# Subdirectory makefile for Vector enumerator GF(p) version

# General control options
#CC = gcc
CFLAGS = $(COPTS)
#
# These options must allow the compiler to find the gmp header files
#  and libraries
#
CPPOPTS = -I.. -I$(GMPHEAD)
LFLAGS= -L$(GMPLIB)
#
# A pregenerated parser is supplied, so bison should not normally be needed
#
BISON=bison 
#
# Finally the ALLOCA flag should be set to alloca.o if your
# system does not supply alloca()
#
ALLOCA=

CPPFLAGS = $(VFLAG) $(CPPOPTS)



objs1= myalloc.o global.o

objs2= 	pack.o allocs.o out.o push.o comline.o qinput.o\
     memain.o me.o coin.o vector.o 

objs3=  scanner.o input.o input.tab.o


objects= $(objs1) $(objs2) $(objs3) $(XTRAOBJS) $(sysobjs)

%.c : %.y

../input.tab.c ../input.tab.h: ../input.y
	$(BISON) -d $^ -o ../$(*).tab.c

../bin/$(PROGNAME).exe : $(objects)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $(objects) -lm $(ARITHLIB)

input.tab.o : ../input.tab.c ../input.tab.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c ../input.tab.c -o $@ 

%.o : ../%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@    

$(objects) : ../global.h ../myalloc.h

$(objs2) $(objs3) $(XTRAOBJS) : ../meint.h $(XTRAGDEPS)

$(objs3) : ../input.h ../input.tab.h

$(XTRAOBJS) : $(XTRADEPS)




