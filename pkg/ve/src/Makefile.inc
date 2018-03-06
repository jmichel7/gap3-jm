# makefile includes for all versions of module enumerator
ifeq ($(ARCH),HP)
CC =c89
CFLAGS=-O
else
CC = gcc
CFLAGS = -O2  -finline-functions -DSYS_IS_USG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SIGNAL_PROT
#CFLAGS = -g
endif

VPATH=..
CPPFLAGS = $(VFLAG) -I$(VPATH) -I/local/gregory/include -Darch_$(ARCH)


LFLAGS= -L/local/gregory/lib
BISON=bison


#vpath %.a /home/maths1/sal/lib/$(ARCH)

objs1= myalloc.o global.o

objs2= 	pack.o allocs.o out.o push.o comline.o qinput.o\
     memain.o me.o coin.o vector.o 

objs3=  scanner.o input.o input.tab.o


objects= $(objs1) $(objs2) $(objs3) $(XTRAOBJS)

%.d : %.c 
	$(CC) -M $(CPPFLAGS) $< | sed 's:$*.o:& $@:g' > $@

%.c : %.y

%.s : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -S $<

$(VPATH)/%.tab.c $(VPATH)/%.tab.h: %.y
	$(BISON) -d $^ -o $(VPATH)/$(*).tab.c

$(VPATH)/bin/$(ARCH)/$(PROGNAME) : $(objects)
	$(CC) $(CFLAGS) $(LFLAGS) -o $@ $(objects) -lm $(ARITHLIB)

.PHONY  : distclean

distclean : 
	$(RM) *.d *.tab* ../*.tab* *.o

input.tab.d : $(VPATH)/input.tab.c $(VPATH)/input.tab.h
	$(CC) -M $(CPPFLAGS) $< |  sed 's:input.tab.o:& $@:g' > $@

input.tab.o : $(VPATH)/input.tab.c $(VPATH)/input.tab.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ 

include $(subst .o,.d,$(objects))

scanner.d : $(VPATH)/input.tab.h
