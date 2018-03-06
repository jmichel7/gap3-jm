##############################################################################
##
##    Makefile                        NQ                       Werner Nickel
##                                         Werner.Nickel@math.rwth-aachen.de
##

##    The follwoing two line define the name of the GNU multiple precision
##    library (GNULIB) and the directory that contains the corresponding
##    include files (GNUINC).
GNULIB      = /usr/lib/x86_64-linux-gnu/libgmp.a
GNUINC      = /usr/local/include
COPTS       = -O2
CFLAGS      = $(COPTS) -I$(GNUINC) -DVERSION='"1.1d, 18 May 1994"'
TARFILE     = ../nq.tar
GAP_PATH    = ..
PKG_PATH    = ${GAP_PATH}/pkg
EXT         = .ext
SRCDISTFILE = nq

OBJ = consistency.o pc.o relations.o word.o system.o \
      addgen.o tails.o eliminate.o engel.o glimt.o \
      presentation.o mem.o pcarith.o collect.o

nq : 
	@echo ""
	@echo "    Please make sure that the variables GNUINC and"
	@echo "    GNULIB in this makefile are set correctly."
	@echo ""
	( cd src; ${MAKE} -f ../Makefile nqsrc )

nqsrc : nq.o $(OBJ) $(GNULIB)
	$(CC) $(CFLAGS) nq.o $(OBJ) $(GNULIB) -o nq

clean:
	@rm -f $(OBJ) nq.o

src_dist:
	@zoo ah ${SRCDISTFILE}.zoo \
		${PKG_PATH}/nq/bin/nq \
		${PKG_PATH}/nq/bin/nq.sh \
		${PKG_PATH}/nq/init.g \
		${PKG_PATH}/nq/src/Makefile \
		${PKG_PATH}/nq/Makefile \
		`find ${PKG_PATH}/nq/examples ${PKG_PATH}/nq/gap -print`  \
		`find ${PKG_PATH}/nq/src -name "*.c" -print` \
		`find ${PKG_PATH}/nq/src -name "*.h" -print` 
	@zoo PE ${SRCDISTFILE}.zoo

bin_dist:
	@( test \! -r ${PKG_PATH}/nq/bin.orig || ( echo \
	 "*** WARNING: 'nq/bin.orig' exists ***";  exit 1 ) )
	( cd ${PKG_PATH}/nq ; mv bin bin.orig ; mkdir bin ; \
	  cp bin.orig/nq.sh bin/nq.sh ; cp bin.orig/nq-${EXT} bin/nq )
	@zoo ah ${SRCDISTFILE}-${EXT}.zoo \
		${PKG_PATH}/nq/bin/nq \
		${PKG_PATH}/nq/bin/nq.sh 
	@zoo PE ${SRCDISTFILE}-${EXT}.zoo
	( cd ${PKG_PATH}/nq ; rm -rf bin.old ; mv bin bin.old ; \
	  mv bin.orig bin )
