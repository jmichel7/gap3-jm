.PHONY : all qme zme me distrib report
all : report
# all : me qme zme

me : 
	(cd ve/$(ARCH); $(MAKE) -f ../GNUmakefile ../../bin/$(ARCH)/$@)

qme : 
	(cd rat/$(ARCH); $(MAKE) -f ../GNUmakefile ../../bin/$(ARCH)/$@)

zme : 
	(cd int/$(ARCH); $(MAKE) -f ../GNUmakefile ../../bin/$(ARCH)/$@)

report : me qme zme
	sh ./report.sh



distrib : nme.tar.Z

nme.tar.Z :  alloca.c allocs.c coin.c comline.c global.c input.c \
    input.tab.c lattice.c me.c memain.c myalloc.c out.c pack.c \
    push.c qinput.c scanner.c scrut.c vector.c global.h ilatt.h \
    input.h input.tab.h latt.h me.h meint.h myalloc.h input.y \
    Makefile Makefile.inc ve/Makefile rat/Makefile int/Makefile \
    build.sh docs/nme.txi README report.sh examples GNUmakefile \
    ve/GNUmakefile rat/GNUmakefile int/GNUmakefile bin/me bin/qme bin/zme
	tar cf nme.tar $^
	compress nme.tar

zoofile: nme.zoo

nme.zoo: alloca.c allocs.c coin.c comline.c global.c input.c \
    input.tab.c lattice.c me.c memain.c myalloc.c out.c pack.c \
    push.c qinput.c scanner.c scrut.c vector.c global.h ilatt.h \
    input.h input.tab.h latt.h me.h meint.h myalloc.h input.y \
    Makefile Makefile.inc ve/Makefile rat/Makefile int/Makefile \
    build.sh docs/nme.txi README report.sh examples GNUmakefile \
    ve/GNUmakefile rat/GNUmakefile int/GNUmakefile bin/me bin/qme bin/zme
	zoo ah nme.zoo $^

