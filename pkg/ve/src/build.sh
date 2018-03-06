#!/bin/sh 
# Select your C compiler here
CC=gcc
# compiler flags for optimization and so on
CFLAGS='-O2 -finline-functions'
#
# flags to dodge round various system problems
# This is the sun4, sunos4, gcc version.
# for other cases look in (or use) Makefile
#
SYSFLAGS='-DSYS_IS_USG -DSYS_HAS_TIME_PROTO -DSYS_HAS_SIGNAL_PROT'

# include directories. This must ensure that the compiler finds 
# the gmp.h header file for GMP
CPPFLAGS='-I/use/slinton/include'
# link flags. These must ensure that the compiler finds libgmp.a
LFLAGS="-L/use/slinton/lib/$ARCH"
# Set this to alloca.o if your system does not provide alloca.
ALLOCA=

comp() {
$CC -I.. $CPPFLAGS $CFLAGS $SYSFLAGS $PROD -c  -o $1.o ../$1.c
}

Commoncomps() {
comp memain
comp me
comp coin
comp input.tab
comp scanner
comp input
comp myalloc
comp global
comp pack
comp allocs
comp out
comp push
comp comline
comp qinput
comp vector
# comp alloca
}

mkdir bin

# delete of comment out from here to the next comment if you do NOT
# want to build me

cd ve
PROD=-DME
Commoncomps
comp scrut
$CC $CFLAGS $LFLAGS -o ../bin/me myalloc.o global.o pack.o allocs.o out.o \
push.o comline.o qinput.o memain.o me.o coin.o vector.o scanner.o \
input.o input.tab.o scrut.o $ALLOCA -lm
cd ..

# End of the me stuff

# delete or comment out from here to the next comment if you do NOT
# want to build qme

cd rat
PROD=-DQME
Commoncomps
comp scrut
$CC $CFLAGS $LFLAGS -o ../bin/qme myalloc.o global.o pack.o allocs.o out.o \
push.o comline.o qinput.o memain.o me.o coin.o vector.o scanner.o \
input.o input.tab.o scrut.o $ALLOCA -lgmp -lm
cd ..

# End of the qme stuff

# delete or comment out from here to the next comment if you do NOT
# want to build zme

cd rat
PROD=-DZME
Commoncomps
comp lattice
$CC $CFLAGS $LFLAGS -o ../bin/zme myalloc.o global.o pack.o allocs.o out.o \
push.o comline.o qinput.o memain.o me.o coin.o vector.o scanner.o \
input.o input.tab.o lattice.o $ALLOCA -lgmp -lm
cd ..

# end of zme stuff

sh ./report.sh


