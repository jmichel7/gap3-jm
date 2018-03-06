## do not edit this definitions #############################################
SHELL=/bin/sh
TARFILE=../grape.tar
SRCDISTFILE=grape
EXT=".exe"
GAP_PATH=..
PKG_PATH=${GAP_PATH}/pkg

## show makefile options ####################################################
unknown:
	@echo "usage: 'make <target> EXT=<ext>' where target is one of"
	@echo "'ibm-i386-386bsd-gcc2' for IBM PCs under 386BSD with GNU cc 2"
	@echo "'ibm-i386-386bsd-cc'   for IBM PCs under 386BSD with cc (GNU)"
	@echo "'sun-sparc-sunos-cc'   for SUN SPARC under SunOS with cc"
	@echo "'bsd-gcc'              for others under Berkeley UNIX with GNU cc 2"
	@echo "'bsd'                  for others under Berkeley UNIX with cc"
	@echo "'usg-gcc'              for others under System V UNIX with GNU cc 2"
	@echo "'usg'                  for others under System V UNIX with cc"
	@echo
	@echo "where <ext> should be a sensible extension, i.e.,"
	@echo "'EXT=-sun-sparc-sunos' for SUN 4 or 'EXT=' if GRAPE only"
	@echo "runs on a single architecture"

bsd-gcc:
	( cd nauty17; make CC=gcc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O -DSYS_BSDUNIX=1" EXT=$(EXT) )
	( cd src; make CC=gcc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

bsd:
	( cd nauty17; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O -DSYS_BSDUNIX=1" EXT=$(EXT) )
	( cd src; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

usg-gcc:
	( cd nauty17; make CC=gcc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O -DSYS_UNIX=1" EXT=$(EXT) )
	( cd src; make CC=gcc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

usg:
	( cd nauty17; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O -DSYS_UNIX=1" EXT=$(EXT) )
	( cd src; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

sun-sparc-sunos-cc:
	( cd nauty17; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O -DSYS_UNIX=1" EXT=$(EXT) )
	( cd src; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

ibm-i386-386bsd-gcc2:
	( cd nauty17; make CC=gcc CFLAGS="$(COPTS) -O2 -DSYS_BSDUNIX=1" EXT=$(EXT) )
	( cd src; make CC=gcc CFLAGS="$(COPTS) -O2" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

ibm-i386-386bsd-cc:
	( cd nauty17; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O -DSYS_BSDUNIX=1" EXT=$(EXT) )
	( cd src; make CC=cc LDFLAGS="$(LOPTS)" CFLAGS="$(COPTS) -O" EXT=$(EXT) )
	@echo "*** Binaries created have extension: $(EXT)"

clean:
	( cd src; make clean EXT=$(EXT) )
	( cd nauty17;  make clean EXT=$(EXT) )

bin_dist:
	@( test \! -r ${PKG_PATH}/grape/bin.orig || ( echo \
	 "*** WARNING: 'grape/bin.orig' exists ***";  exit 1 ) )
	( cd ${PKG_PATH}/grape ; mv bin bin.orig ; mkdir bin ; \
	  cp bin.orig/drcanon3.sh bin/drcanon3.sh ; \
	  cp bin.orig/drcanon3-${EXT} bin/drcanon3 ; \
	  cp bin.orig/dreadnaut.sh bin/dreadnaut.sh ; \
	  cp bin.orig/dreadnaut-${EXT} bin/dreadnaut ; \
	  cp bin.orig/drtogap3.sh bin/drtogap3.sh ; \
	  cp bin.orig/drtogap3-${EXT} bin/drtogap3 ; \
	  cp bin.orig/enum3.sh bin/enum3.sh ; \
	  cp bin.orig/enum3.sh bin/enum3 ; \
	  cp bin.orig/enum3ca.sh bin/enum3ca.sh ; \
	  cp bin.orig/enum3ca.sh bin/enum3ca ; \
	  cp bin.orig/gap3todr.sh bin/gap3todr.sh ; \
	  cp bin.orig/gap3todr-${EXT} bin/gap3todr )
	@zoo ah ${SRCDISTFILE}-${EXT}.zoo \
		${PKG_PATH}/grape/bin/drcanon3 \
		${PKG_PATH}/grape/bin/gap3todr \
		${PKG_PATH}/grape/bin/drtogap3 \
		${PKG_PATH}/grape/bin/dreadnaut \
		${PKG_PATH}/grape/bin/enum3 \
		${PKG_PATH}/grape/bin/enum3ca \
		${PKG_PATH}/grape/bin/drcanon3.sh \
		${PKG_PATH}/grape/bin/gap3todr.sh \
		${PKG_PATH}/grape/bin/drtogap3.sh \
		${PKG_PATH}/grape/bin/dreadnaut.sh \
		${PKG_PATH}/grape/bin/enum3.sh \
		${PKG_PATH}/grape/bin/enum3ca.sh 
	@zoo PE ${SRCDISTFILE}-${EXT}.zoo
	( cd ${PKG_PATH}/grape ; rm -rf bin.old ; mv bin bin.old ; \
	  mv bin.orig bin )
