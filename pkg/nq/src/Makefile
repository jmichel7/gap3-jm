##############################################################################
##
##    Makefile                        NQ                       Werner Nickel
##                                         Werner.Nickel@math.rwth-aachen.de
##

GNULIB      = -lgmp
GNUINC      = /usr/local/include
COPTS       = -O
CFLAGS      = $(COPTS) -I$(GNUINC) -DVERSION='"1.1d, 18 May 1994"' -m32
TARFILE     = ../nq.tar
GAP_PATH    = ..
PKG_PATH    = ${GAP_PATH}/pkg
EXT         = .ext
SRCDISTFILE = nq

OBJ = consistency.o pc.o relations.o word.o system.o \
      addgen.o tails.o eliminate.o engel.o glimt.o \
      presentation.o mem.o pcarith.o collect.o

../bin/nq: nq.o $(OBJ)
	$(CC) $(CFLAGS) nq.o $(OBJ) $(GNULIB) -o ../bin/nq

clean:
	@rm -f $(OBJ) nq.o
