#############################################################################
##
#A  Makefile                    XGAP source                      Frank Celler
##
#H  @(#)$Id: Makefile,v 1.1.1.1 1996/12/11 12:39:45 werner Exp $
##
#Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: Makefile,v $
#H  Revision 1.1.1.1  1996/12/11 12:39:45  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.4  1995/10/30  09:06:00  fceller
#H  added missing doc.x11/Makefile
#H
#H  Revision 1.3  1995/10/25  09:59:54  fceller
#H  added zoo option
#H
#H  Revision 1.2  1995/09/25  19:42:02  fceller
#H  added german translation resource database
#H
#H  Revision 1.1  1995/09/25  11:08:48  fceller
#H  Initial revision
##
SHELL=/bin/sh
XGAP=.
GAP=gap3r4p2


## DOC: documentation for XGAP ##############################################
DOC = $(XGAP)/doc/*.tex $(XGAP)/doc/TODO $(XGAP)/doc/Makefile \
      ${XGAP}/doc/*.dvi $(XGAP)/ANNOUNCE


## GRP: a group libraray for XGAP ###########################################
GRP = $(XGAP)/grp/*.grp


## LIB: library files for XGAP ##############################################
LIB = $(XGAP)/lib/*.g


## PMG: various permutations groups for 'EpimorphismsFpGroup' ###############
PMG = $(XGAP)/pmg/*.grp


## TST: test files for XGAP #################################################
TST = $(XGAP)/tst/*.g


## BIN_X11: X11 shell script to start XGAP ##################################
BIN_X11 = $(XGAP)/bin.x11/xgap $(XGAP)/bin.x11/xgap.sh


## DOC_X11: documentation for XGAP under X11 ################################
DOC_X11 = $(XGAP)/doc.x11/*.tex $(XGAP)/doc.x11/TODO $(XGAP)/doc.x11/Makefile


## RDB_X11: X11 resource database defaults and examples #####################
RDB_X11 = $(XGAP)/rdb.x11/default $(XGAP)/rdb.x11/expl.* \
          $(XGAP)/rdb.x11/serv.* $(XGAP)/rdb.x11/german


## SRC_X11: source files for X11 version of XGAP ############################
SRC_X11 = $(XGAP)/src.x11/*.[ch] $(XGAP)/src.x11/*.bm \
	  $(XGAP)/src.x11/Makefile $(XGAP)/src.x11/HELP \
	  $(XGAP)/src.x11/TARGETS \
          $(XGAP)/src.x11/*.in $(XGAP)/src.x11/configure


## DEFAULT: give a list of possible options #################################
TAR_BCK = xgap`date | awk '{m["Jan"]=1;m["Feb"]=2;m["Mar"]=3;m["Apr"]=4;m["May"]=5;m["Jun"]=6;m["Jul"]=7;m["Aug"]=8;m["Sep"]=9;m["Oct"]=10;m["Nov"]=11;m["Dec"]=12;printf("%d%d",m[$$2],$$3);}'`.tar

default:
	@echo
	@echo "The following targets are supported:"
	@echo
	@echo " 'backup'                backup XGAP in $(TAR_BCK)"
	@echo " 'zoo'                   backup XGAP in a zoo archive"
	@echo " 'clean'                 remove all created files"
	@echo " 'x11'                   create XGAP for X-Windows X11R5"
	@echo


## BACKUP: create a tar backup ##############################################
backup:
	@echo "removing emacs ~ backup file, creating backup file $(TAR_BCK)"
	@rm -f *~ */*~
	@tar cf $(TAR_BCK)					\
		$(DOC) $(GRP) $(LIB) $(PMG) $(TST)		\
		$(BIN_X11) $(DOC_X11) $(RDB_X11) $(SRC_X11)	\
		$(XGAP)/INSTALL.X11 $(XGAP)/Makefile

## ZOO: create a zoo archive for distribution ###############################
zoo:
	@echo "removing emacs ~ backup file, creating dirs & links in /tmp"
	@rm -f *~ */*~
	@mkdir /tmp/$(GAP) /tmp/$(GAP)/pkg
	@ln -s `pwd`/../xgap /tmp/$(GAP)/pkg/xgap
	@( cd /tmp; make -f $(GAP)/pkg/xgap/Makefile \
           XGAP=./gap3r4p2/pkg/xgap zoome )
	@rm /tmp/$(GAP)/pkg/xgap
	@rmdir /tmp/$(GAP)/pkg /tmp/$(GAP)
	@mv /tmp/xgap.zoo xgap1r3.zoo
	@mv /tmp/xgap-dec-alpha-osf1-x11r5.zoo \
		xgap1r3-dec-alpha-osf1-x11r5.zoo
	@mv /tmp/xgap-dec-mips-ultrix42-x11r5.zoo \
		xgap1r3-dec-mips-ultrix42-x11r5.zoo
	@mv /tmp/xgap-hp-hppa1.1-hpux9-x11r5.zoo \
		xgap1r3-hp-hppa1.1-hpux9-x11r5.zoo
	@mv /tmp/xgap-sun-sparc-sunos412-x11r5.zoo \
		xgap1r3-sun-sparc-sunos412-x11r5.zoo

zoome:
	@zoo ah /tmp/xgap.zoo					\
		$(DOC) $(GRP) $(LIB) $(PMG) $(TST)		\
		$(BIN_X11) $(DOC_X11) $(RDB_X11) $(SRC_X11)	\
		$(XGAP)/INSTALL.X11 $(XGAP)/Makefile
	@zoo ah /tmp/xgap-dec-alpha-osf1-x11r5.zoo		\
		$(XGAP)/bin.x11/xgap-dec-alpha-osf1-x11r5
	@zoo ah /tmp/xgap-dec-mips-ultrix42-x11r5.zoo	\
		$(XGAP)/bin.x11/xgap-dec-mips-ultrix42-x11r5
	@zoo ah /tmp/xgap-hp-hppa1.1-hpux9-x11r5.zoo	\
		$(XGAP)/bin.x11/xgap-hp-hppa1.1-hpux9-x11r5
	@zoo ah /tmp/xgap-sun-sparc-sunos412-x11r5.zoo	\
		$(XGAP)/bin.x11/xgap-sun-sparc-sunos412-x11r5

src_dist:
	@zoo ah ${SRCDISTFILE}.zoo \
		$(PKG_PATH)/xgap/ANNOUNCE \
		$(PKG_PATH)/xgap/INSTALL.X11 \
		$(PKG_PATH)/xgap/Makefile \
		$(PKG_PATH)/xgap/doc/Makefile \
		$(PKG_PATH)/xgap/doc/TODO \
		$(PKG_PATH)/xgap/doc/*.tex \
		$(PKG_PATH)/xgap/grp/*.grp \
		$(PKG_PATH)/xgap/lib/*.g \
		$(PKG_PATH)/xgap/pmg/*.grp \
		$(PKG_PATH)/xgap/tst/*.g \
		$(PKG_PATH)/xgap/bin.x11/xgap \
		$(PKG_PATH)/xgap/bin.x11/xgap.sh \
		$(PKG_PATH)/xgap/doc.x11/Makefile \
		$(PKG_PATH)/xgap/doc.x11/TODO \
		$(PKG_PATH)/xgap/doc.x11/*.tex \
		$(PKG_PATH)/xgap/rdb.x11/default \
		$(PKG_PATH)/xgap/rdb.x11/german \
		$(PKG_PATH)/xgap/rdb.x11/expl.* \
		$(PKG_PATH)/xgap/rdb.x11/serv.* \
		$(PKG_PATH)/xgap/src.x11/HELP \
		$(PKG_PATH)/xgap/src.x11/Makefile \
		$(PKG_PATH)/xgap/src.x11/TARGETS \
		$(PKG_PATH)/xgap/src.x11/*.in \
		$(PKG_PATH)/xgap/src.x11/configure \
		$(PKG_PATH)/xgap/src.x11/*.[ch] \
		$(PKG_PATH)/xgap/src.x11/*.bm
	@zoo PE ${SRCDISTFILE}.zoo

## CLEAN: clean up, remove all files created during compilation #############
clean:
	( cd src.x11 ; make clean SHELL="$(SHELL)" )
	( cd doc.x11 ; make clean SHELL="$(SHELL)" )
	( cd doc     ; make clean SHELL="$(SHELL)" )


## X11: compile XGAP for X-Windows X11R5 ####################################
x11:
	( cd src.x11 ; \
          make auto SHELL="$(SHELL)" LOPTS="$(LOPTS)" COPTS="$(COPTS)" )


