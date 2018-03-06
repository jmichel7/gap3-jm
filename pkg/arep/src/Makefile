# Makefile for the GAP/C-interface to Leon's 'desauto'
# SE, 12. Mar 1998

all:	../bin/leonin ../bin/leonout

../bin/leonin:	leonin.c
	cc -o ../bin/leonin leonin.c

../bin/leonout: leonout.c
	cc -o ../bin/leonout leonout.c
clean:
	-rm ../bin/leonin ../bin/leonout
