# General makefile for KBMAG - all executables are put into the bin directory.
SHELL = /bin/sh
SYSTEM = SYSTEMV
RM = /bin/rm

CC = gcc
COPTS = -O

all:
	cd lib; $(MAKE) CC="$(CC)" COPTS="$(COPTS)"
	cd src; $(MAKE) CC="$(CC)" SYSTEM="$(SYSTEM)" COPTS="$(COPTS)"

clean:
	cd lib; $(MAKE) clean
	cd src; $(MAKE) clean

distclean:
	cd lib; $(MAKE) clean
	cd src; $(MAKE) clean
	$(RM) -f bin/*
	cd ag_data; $(RM) -f *.* 
	cd kb_data; $(RM) -f *.* 
	cd fsa_data; $(RM) -f *.* 
	cd doc; $(RM) -f *.aux *.toc *.log *.dvi
	cd gapdoc; $(RM) -f *.aux *.toc *.log *.dvi

test:
	cd kb_data; ../bin/kbprog nilp2; rm nilp2.*
	cd ag_data; ../bin/autgroup 235; ../bin/fsacount 235.wa; rm 235.*

# installation for GAP
GAP_PATH = ../..
PKG_PATH = $(GAP_PATH)/pkg
SRCDISTFILE = kbmag
