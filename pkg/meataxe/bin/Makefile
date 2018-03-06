########################################################################
#####                C MeatAxe - Makefile for UNIX                 #####
########################################################################

########################################################################
# Configuration section                                            #####
# The following definitions may need some changes                  #####
# depending on your operating system and compiler                  #####
########################################################################

# --------------------------------------------------------------------
# The source directory (relative to the directory where one step
# upstairs this Makefile resides). Normally, no change is needed here.
# --------------------------------------------------------------------
SD=../src

# --------------------------------------------------------------------
# The library directory. This is the place where arithmetic
# tables and some other input files are searched. The default
# is "../lib" - but we need the full path.
# --------------------------------------------------------------------
MTXLIB=`(cd ..; pwd)`/lib/
# Here is an example for a different file directory:
#MTXLIB=/usr/joe/mtx/lib/

# --------------------------------------------------------------------
# The name and standard options of your C compiler
# --------------------------------------------------------------------

# The default is to use GNU C:
CC=gcc	
#CC=/rem/scratch.2/gcc/bin/gcc
COMPFLAGS=-ansi -O -m32
#COMPFLAGS=-ansi -g -Wall # for debugging

# If you don"t have gcc, use the standard C compiler:
#CC=cc
#COMPFLAGS=-O

# --------------------------------------------------------------------
# Specify your operating system here. If you don't know which system
# you have, try `uname -rs'. If there is no definition for your OS,
# try DEFAULT. If this fails, try the other systems.
# If this fails, too,  edit the file 'meataxe.h' and follow the
# instructions given there...
# --------------------------------------------------------------------

# DEFAULT is known to work with:
#   OSF/1 1.2 and 2.0
#   NetBSD 0.9
#   FreeBSD 1.2
#   Linux 1.0
OS=DEFAULT

#OS=SUNOS_41   # SunOS 4.1 (Tested with 4.1 and 4.1.3)
#OS=ULTRIX     # DEC Ultrix (Tested with 4.2 and 4.3)
#OS=HPUX       # HP/UX (Tested with 9.0)
#OS=IBMVM      # IBM's VM/CMS (not tested, probably broken)
#OS=MSDOS      # MS-DOS + djgpp (not tested, probably broken)

# --------------------------------------------------------------------
# The name and options of the 'ar' command. The following is
# 'standard' on most systems. If you get errors from ar commands,
# try commenting out these lines to use the definitions which are
# (hopefully) built into your make.
# --------------------------------------------------------------------
AR=ar
ARFLAGS=rv

# --------------------------------------------------------------------
# The ranlib command. The default should work on almost any system.
# --------------------------------------------------------------------
RANLIB=`if [ -r /bin/ranlib -o -r /usr/bin/ranlib ]; \
    then echo ranlib; else echo true; fi`

# --------------------------------------------------------------------
# The name of the ZZZ module and the table file generator.
# --------------------------------------------------------------------
# small version, up to GF(256)
ZZZ=zzz
MAKETAB=maketab
# big version, up to GF(2^16)
#ZZZNAME=bigzzz
#MAKETAB=bigmktab


########################################################################
# End of configuration section                                     #####
# There should be no need to change anything below this point      #####
########################################################################

CFLAGS=$(COMPFLAGS) -DOS_$(OS)
SHELL=/bin/sh

###########################
# C MEATAXE (common part) #
###########################

# Program names
# -------------
MtxH=$(SD)/meataxe.h
MTXPRG=zad zbl zcf zcl zcp zct zcv zef zev zfr ziv zkd zmo zmu znu \
 zor zpc zpr zpt zqt zsb zsi zsm zsp zsy ztc zte ztm ztr zuk zvp \
 chop pwkond mkcycl mkinc mksub genmod mkdotl mkgraph \
 conv $(MAKETAB) zzztest checksum $(LOCALPRGS)

# Other macros
# ------------
CCOMP=$(CC) -c $(CFLAGS) -o $@
COMPILE=$(CC) -c $(CFLAGS)
BUILD=$(CC) $(CFLAGS) -o $@ $(SD)/$@.c $(LLIBS)


# Targets
# -------
all:
	$(MAKE) -f Makefile $(MTXPRG)

lattice: chop pwkond mkcycl mkinc mksub genmod mkdotl $(MAKETAB)
strip:
	strip $(MTXPRG) $(LOCALPRGS)
clean:
	rm -f *.o core
delete:
	rm -f $(MTXPRG) $(LOCALPRGS) *.o *.a core

# Auxiliary programs
# ------------------
$(MAKETAB): $(SD)/$(MAKETAB).c libmtx.a $(MtxH)
	$(BUILD) libmtx.a

conv: $(SD)/conv.c libmtx.a
	$(BUILD) libmtx.a

zzztest: $(SD)/zzztest.c libmtx.a
	$(BUILD) libmtx.a

checksum: $(SD)/checksum.c
	$(BUILD)

# MeatAxe library
# ---------------
libmtx.a:\
	libmtx.a(args.o)\
	libmtx.a(berlekmp.o)\
	libmtx.a(charpol.o)\
	libmtx.a(chbasis.o)\
	libmtx.a(fpoly.o)\
	libmtx.a(gauss.o)\
	libmtx.a(gcd.o)\
	libmtx.a(intstr.o)\
	libmtx.a(matins.o)\
	libmtx.a(message.o)\
	libmtx.a(minpol.o)\
	libmtx.a(mtxobj.o)\
	libmtx.a(os.o)\
	libmtx.a(profile.o)\
	libmtx.a(prtimes.o)\
	libmtx.a(random.o)\
	libmtx.a(sbasis.o)\
	libmtx.a(sequence.o)\
	libmtx.a(spinup.o)\
	libmtx.a(split.o)\
	libmtx.a(words.o)\
	libmtx.a(ycomp.o)\
	libmtx.a(ymatrix.o)\
	libmtx.a(yperm.o)\
	libmtx.a(ypoly.o)\
	libmtx.a(ysets.o)\
	libmtx.a(zbitstr.o)\
	libmtx.a(zfile.o)\
	libmtx.a(zgap.o)\
	libmtx.a(zgauss.o)\
	libmtx.a(zmatinv.o)\
	libmtx.a(znullsp.o)\
	libmtx.a(zpseed.o)\
	libmtx.a(zquot.o)\
	libmtx.a(zsbasis.o)\
	libmtx.a(zspin.o)\
	libmtx.a($(ZZZ).o)\
	libmtx.a(zzz2.o)
	$(RANLIB) libmtx.a

libmtx.a(args.o): $(SD)/args.c $(MtxH)
	$(COMPILE) -o args.o $(SD)/args.c
	$(AR) $(ARFLAGS) libmtx.a args.o
	rm -f args.o

libmtx.a(berlekmp.o): $(SD)/berlekmp.c $(MtxH)
	$(COMPILE) -o berlekmp.o $(SD)/berlekmp.c
	$(AR) $(ARFLAGS) libmtx.a berlekmp.o
	rm -f berlekmp.o

libmtx.a(charpol.o): $(SD)/charpol.c $(MtxH)
	$(COMPILE) -o charpol.o $(SD)/charpol.c
	$(AR) $(ARFLAGS) libmtx.a charpol.o
	rm -f charpol.o

libmtx.a(chbasis.o): $(SD)/chbasis.c $(MtxH)
	$(COMPILE) -o chbasis.o $(SD)/chbasis.c
	$(AR) $(ARFLAGS) libmtx.a chbasis.o
	rm -f chbasis.o

libmtx.a(fpoly.o): $(SD)/fpoly.c $(MtxH)
	$(COMPILE) -o fpoly.o $(SD)/fpoly.c
	$(AR) $(ARFLAGS) libmtx.a fpoly.o
	rm -f fpoly.o

libmtx.a(gauss.o): $(SD)/gauss.c $(MtxH)
	$(COMPILE) -o gauss.o $(SD)/gauss.c
	$(AR) $(ARFLAGS) libmtx.a gauss.o
	rm -f gauss.o

libmtx.a(gcd.o): $(SD)/gcd.c $(MtxH)
	$(COMPILE) -o gcd.o $(SD)/gcd.c
	$(AR) $(ARFLAGS) libmtx.a gcd.o
	rm -f gcd.o

libmtx.a(intstr.o): $(SD)/intstr.c $(MtxH)
	$(COMPILE) -o intstr.o $(SD)/intstr.c
	$(AR) $(ARFLAGS) libmtx.a intstr.o
	rm -f intstr.o

libmtx.a(matins.o): $(SD)/matins.c $(MtxH)
	$(COMPILE) -o matins.o $(SD)/matins.c
	$(AR) $(ARFLAGS) libmtx.a matins.o
	rm -f matins.o

libmtx.a(message.o): $(SD)/message.c $(MtxH)
	$(COMPILE) -o message.o $(SD)/message.c
	$(AR) $(ARFLAGS) libmtx.a message.o
	rm -f message.o

libmtx.a(minpol.o): $(SD)/minpol.c $(MtxH)
	$(COMPILE) -o minpol.o $(SD)/minpol.c
	$(AR) $(ARFLAGS) libmtx.a minpol.o
	rm -f minpol.o

libmtx.a(mtxobj.o): $(SD)/mtxobj.c $(MtxH)
	$(COMPILE) -o mtxobj.o $(SD)/mtxobj.c
	$(AR) $(ARFLAGS) libmtx.a mtxobj.o
	rm -f mtxobj.o

libmtx.a(os.o): $(SD)/os.c $(MtxH)
	$(COMPILE) $(OSFLAGS) -DMTXLIB=\"$(MTXLIB)\" -o os.o $(SD)/os.c
	$(AR) $(ARFLAGS) libmtx.a os.o
	rm -f os.o

libmtx.a(sequence.o): $(SD)/sequence.c $(MtxH)
	$(COMPILE) -o sequence.o $(SD)/sequence.c
	$(AR) $(ARFLAGS) libmtx.a sequence.o
	rm -f sequence.o

libmtx.a(words.o): $(SD)/words.c $(MtxH)
	$(COMPILE) -o words.o $(SD)/words.c
	$(AR) $(ARFLAGS) libmtx.a words.o
	rm -f words.o

libmtx.a(ycomp.o): $(SD)/ycomp.c $(MtxH)
	$(COMPILE) -o ycomp.o $(SD)/ycomp.c
	$(AR) $(ARFLAGS) libmtx.a ycomp.o
	rm -f ycomp.o

libmtx.a(ymatrix.o): $(SD)/ymatrix.c $(MtxH)
	$(COMPILE) -o ymatrix.o $(SD)/ymatrix.c
	$(AR) $(ARFLAGS) libmtx.a ymatrix.o
	rm -f ymatrix.o

libmtx.a(yperm.o): $(SD)/yperm.c $(MtxH)
	$(COMPILE) -o yperm.o $(SD)/yperm.c
	$(AR) $(ARFLAGS) libmtx.a yperm.o
	rm -f yperm.o

libmtx.a(ypoly.o): $(SD)/ypoly.c $(MtxH)
	$(COMPILE) -o ypoly.o $(SD)/ypoly.c
	$(AR) $(ARFLAGS) libmtx.a ypoly.o
	rm -f ypoly.o

libmtx.a(ysets.o): $(SD)/ysets.c $(MtxH)
	$(COMPILE) -o ysets.o $(SD)/ysets.c
	$(AR) $(ARFLAGS) libmtx.a ysets.o
	rm -f ysets.o

libmtx.a(zbitstr.o): $(SD)/zbitstr.c $(MtxH)
	$(COMPILE) -o zbitstr.o $(SD)/zbitstr.c
	$(AR) $(ARFLAGS) libmtx.a zbitstr.o
	rm -f zbitstr.o

libmtx.a(zfile.o): $(SD)/zfile.c $(MtxH)
	$(COMPILE) -o zfile.o $(SD)/zfile.c
	$(AR) $(ARFLAGS) libmtx.a zfile.o
	rm -f zfile.o

libmtx.a(zgap.o): $(SD)/zgap.c $(MtxH)
	$(COMPILE) -o zgap.o $(SD)/zgap.c
	$(AR) $(ARFLAGS) libmtx.a zgap.o
	rm -f zgap.o

libmtx.a(zgauss.o): $(SD)/zgauss.c $(MtxH)
	$(COMPILE) -o zgauss.o $(SD)/zgauss.c
	$(AR) $(ARFLAGS) libmtx.a zgauss.o
	rm -f zgauss.o

libmtx.a(zmatinv.o): $(SD)/zmatinv.c $(MtxH)
	$(COMPILE) -o zmatinv.o $(SD)/zmatinv.c
	$(AR) $(ARFLAGS) libmtx.a zmatinv.o
	rm -f zmatinv.o

libmtx.a(znullsp.o): $(SD)/znullsp.c $(MtxH)
	$(COMPILE) -o znullsp.o $(SD)/znullsp.c
	$(AR) $(ARFLAGS) libmtx.a znullsp.o
	rm -f znullsp.o

libmtx.a(zpseed.o): $(SD)/zpseed.c $(MtxH)
	$(COMPILE) -o zpseed.o $(SD)/zpseed.c
	$(AR) $(ARFLAGS) libmtx.a zpseed.o
	rm -f zpseed.o

libmtx.a(zquot.o): $(SD)/zquot.c $(MtxH)
	$(COMPILE) -o zquot.o $(SD)/zquot.c
	$(AR) $(ARFLAGS) libmtx.a zquot.o
	rm -f zquot.o

libmtx.a(zsbasis.o): $(SD)/zsbasis.c $(MtxH)
	$(COMPILE) -o zsbasis.o $(SD)/zsbasis.c
	$(AR) $(ARFLAGS) libmtx.a zsbasis.o
	rm -f zsbasis.o

libmtx.a(zspin.o): $(SD)/zspin.c $(MtxH)
	$(COMPILE) -o zspin.o $(SD)/zspin.c
	$(AR) $(ARFLAGS) libmtx.a zspin.o
	rm -f zspin.o

libmtx.a(profile.o): $(SD)/profile.c $(MtxH)
	$(COMPILE) -o profile.o $(SD)/profile.c
	$(AR) $(ARFLAGS) libmtx.a profile.o
	rm -f profile.o

libmtx.a(prtimes.o): $(SD)/prtimes.c $(MtxH)
	$(COMPILE) -o prtimes.o $(SD)/prtimes.c
	$(AR) $(ARFLAGS) libmtx.a prtimes.o
	rm -f prtimes.o

libmtx.a(random.o): $(SD)/random.c $(MtxH)
	$(COMPILE) -o random.o $(SD)/random.c
	$(AR) $(ARFLAGS) libmtx.a random.o
	rm -f random.o

libmtx.a(sbasis.o): $(SD)/sbasis.c $(MtxH)
	$(COMPILE) -o sbasis.o $(SD)/sbasis.c
	$(AR) $(ARFLAGS) libmtx.a sbasis.o
	rm -f sbasis.o

libmtx.a(spinup.o): $(SD)/spinup.c $(MtxH)
	$(COMPILE) -o spinup.o $(SD)/spinup.c
	$(AR) $(ARFLAGS) libmtx.a spinup.o
	rm -f spinup.o

libmtx.a(split.o): $(SD)/split.c $(MtxH)
	$(COMPILE) -o split.o $(SD)/split.c
	$(AR) $(ARFLAGS) libmtx.a split.o
	rm -f split.o


libmtx.a($(ZZZ).o): $(SD)/$(ZZZ).c $(MtxH)
	$(COMPILE) -DZZZ_CC='"$(CC) $(CFLAGS)"' -o $(ZZZ).o $(SD)/$(ZZZ).c
	$(AR) $(ARFLAGS) libmtx.a $(ZZZ).o
	rm -f $(ZZZ).o

libmtx.a(zzz2.o): $(SD)/zzz2.c $(MtxH)
	$(COMPILE) -o zzz2.o $(SD)/zzz2.c
	$(AR) $(ARFLAGS) libmtx.a zzz2.o
	rm -f zzz2.o

# Objects
# -------
chop.o: $(SD)/chop.c $(SD)/lattice.h $(MtxH) $(SD)/files.h
	$(CCOMP) $(SD)/chop.c

files.o: $(SD)/files.c $(SD)/lattice.h
	$(CCOMP) $(SD)/files.c

genmod.o: $(SD)/genmod.c $(SD)/lattice.h
	$(CCOMP) $(SD)/genmod.c

pwkond.o: $(SD)/pwkond.c $(MtxH) $(SD)/lattice.h
	$(CCOMP) $(SD)/pwkond.c

mkcycl.o: $(SD)/mkcycl.c $(SD)/lattice.h $(MtxH)
	$(CCOMP) $(SD)/mkcycl.c

mkdotl.o: $(SD)/mkdotl.c $(SD)/lattice.h $(MtxH)
	$(CCOMP) $(SD)/mkdotl.c

mkgraph.o: $(SD)/mkgraph.c $(SD)/lattice.h $(MtxH)
	$(CCOMP) $(SD)/mkgraph.c

mkinc.o: $(SD)/mkinc.c $(SD)/lattice.h $(MtxH)
	$(CCOMP) $(SD)/mkinc.c

mksub.o: $(SD)/mksub.c $(SD)/lattice.h $(MtxH)
	$(CCOMP) $(SD)/mksub.c

# Meat-Axe programs
# ----------------- 

zad: $(SD)/zad.c libmtx.a
	$(BUILD) libmtx.a

zbl: $(SD)/zbl.c libmtx.a
	$(BUILD) libmtx.a

zcf: $(SD)/zcf.c libmtx.a
	$(BUILD) libmtx.a

zcl: $(SD)/zcl.c libmtx.a
	$(BUILD) libmtx.a

zcp: $(SD)/zcp.c libmtx.a
	$(BUILD) libmtx.a

zct: $(SD)/zct.c libmtx.a
	$(BUILD) libmtx.a

zcv: $(SD)/zcv.c libmtx.a
	$(BUILD) libmtx.a

zef: $(SD)/zef.c libmtx.a
	$(BUILD) libmtx.a

zev: $(SD)/zev.c libmtx.a
	$(BUILD) libmtx.a

zfr: $(SD)/zfr.c libmtx.a
	$(BUILD) libmtx.a

ziv: $(SD)/ziv.c libmtx.a
	$(BUILD) libmtx.a

zkd: $(SD)/zkd.c libmtx.a
	$(BUILD) libmtx.a

zmo: $(SD)/zmo.c libmtx.a
	$(BUILD) libmtx.a

zmu: $(SD)/zmu.c libmtx.a
	$(BUILD) libmtx.a

znu: $(SD)/znu.c libmtx.a
	$(BUILD) libmtx.a

zor: $(SD)/zor.c libmtx.a
	$(BUILD) libmtx.a

zpc: $(SD)/zpc.c libmtx.a
	$(BUILD) libmtx.a

zpr: $(SD)/zpr.c libmtx.a
	$(BUILD) libmtx.a

zpt: $(SD)/zpt.c libmtx.a
	$(BUILD) libmtx.a

zqt: $(SD)/zqt.c libmtx.a
	$(BUILD) libmtx.a

zsb: $(SD)/zsb.c libmtx.a
	$(BUILD) libmtx.a

zsi: $(SD)/zsi.c libmtx.a
	$(BUILD) libmtx.a

zsm: $(SD)/zsm.c libmtx.a
	$(BUILD) libmtx.a

zsp: $(SD)/zsp.c libmtx.a
	$(BUILD) libmtx.a

zsy: $(SD)/zsy.c libmtx.a
	$(BUILD) libmtx.a

ztc: $(SD)/ztc.c libmtx.a
	$(BUILD) libmtx.a

zte: $(SD)/zte.c libmtx.a
	$(BUILD) libmtx.a

ztm: $(SD)/ztm.c libmtx.a
	$(BUILD) libmtx.a

ztr: $(SD)/ztr.c libmtx.a
	$(BUILD) libmtx.a

zuk: $(SD)/zuk.c libmtx.a
	$(BUILD) libmtx.a

zvp: $(SD)/zvp.c libmtx.a
	$(BUILD) libmtx.a


################################################################
#
#   Lattice package
#
################################################################

OBJ1=files.o chop.o libmtx.a
OBJ3=files.o pwkond.o libmtx.a
OBJ4=files.o mkcycl.o libmtx.a
OBJ5=files.o mkinc.o libmtx.a
OBJ5a=files.o mkdotl.o libmtx.a
OBJ6=mksub.o files.o libmtx.a
OBJ7=genmod.o files.o libmtx.a

# -------------------------------
# Rules for building the programs
# -------------------------------

lattice: $(LPROGS) $(MAKETAB)

chop: $(OBJ1)
	$(CC) $(CFLAGS) -o chop $(OBJ1) $(LLIBS)

pwkond: $(OBJ3)
	$(CC) $(CFLAGS) -o pwkond $(OBJ3) $(LLIBS)

mkcycl: $(OBJ4)
	$(CC) $(CFLAGS) -o mkcycl $(OBJ4) $(LLIBS)

mkdotl: $(OBJ5a)
	$(CC) $(CFLAGS) -o mkdotl $(OBJ5a) $(LLIBS)

mkinc: $(OBJ5)
	$(CC) $(CFLAGS) -o mkinc $(OBJ5) $(LLIBS)

mksub: $(OBJ6)
	$(CC) $(CFLAGS) -o mksub $(OBJ6) $(LLIBS)

genmod: $(OBJ7)
	$(CC) $(CFLAGS) -o genmod $(OBJ7) $(LLIBS)

mkgraph: mkgraph.o files.o libmtx.a
	$(CC) $(CFLAGS) -o mkgraph mkgraph.o files.o libmtx.a $(LLIBS) 

# -------------------------------------------------
