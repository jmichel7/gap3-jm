# Subdirectory makefile for Vector enumerator Integer version

override VFLAG=-DZME # define this when compiling
override XTRAOBJS= lattice.o #  extra object files (eg lattice.c)
override XTRADEPS= ilatt.h
override XTRAGDEPS = latt.h
override PROGNAME=zme # what to call the final program
override ARITHLIB=  -lgmp  # no extra library declarations 

include ../../Makefile.inc

