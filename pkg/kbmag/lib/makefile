# lib makefile
SHELL=/bin/sh
RANLIB = ranlib
RM = /bin/rm
CC= gcc
COPTS= -O
LIB = fsalib.a
LIBOBJ = fsa.o fsaio.o hash.o rwsreduce.o miscio.o rwsio.o rwsio2.o worddiff.o \
	 diffreduce.o fsalogic.o fsanosubword.o fsawa.o fsatriples.o \
	 fsacheckmult.o fsacomposite.o fsaipmin.o fsaminkb.o rabkar.o
$(LIB): $(LIBOBJ)
	ar cr $(LIB) $(LIBOBJ)
	if  which ranlib > /dev/null 2> /dev/null; then ranlib $(LIB); fi
clean:
	$(RM) -f $(LIBOBJ) $(LIB)
	$(RM) -f *\~ core a.out
.c.o:
	${CC} -c $(COPTS) $*.c 
