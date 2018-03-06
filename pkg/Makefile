unknown:
	@echo "usage: 'make <target>' where target is one of"
	@echo "packages: compile all packages for linux/gcc"
	@echo "cleanpackageso: cleanup all non-necessary binaries built"
	@echo "                while making packages."
	@echo "cleanpackages: delete all binaries in all packages to get
	@echo "               a source-only distribution"
FORCE:
anupq: FORCE
	(cd anupq/src; $(MAKE))
cleananupqo:
	(cd anupq/src; $(MAKE) cleano)
cleananupq: cleananupqo
	(cd anupq/src; $(MAKE) clean)
anusq: FORCE
	(cd anusq/src; $(MAKE))
cleananusqo:
	(cd anusq/src; $(MAKE) cleano)
cleananusq: cleananusqo
	(cd anusq/src; $(MAKE) clean)
arep: FORCE
	(cd arep/src; $(MAKE))
cleanarep:
	(cd arep/src; $(MAKE) clean)
cohomolo: FORCE
	(cd cohomolo/src; $(MAKE))
cleancohomoloo:
	(cd cohomolo/src; $(MAKE) cleano)
cleancohomolo: cleancohomoloo
	(cd cohomolo/src; $(MAKE) clean)
grape: FORCE
	(cd grape/nauty17; $(MAKE))
	(cd grape/src; $(MAKE))
cleangrapeo:
	-rm grape/nauty17/*.o
cleangrape: cleangrapeo
	-(cd grape/bin; rm dreadnaut gap3todr drtogap3 drcanon3)
guava: FORCE
	(cd guava/src; $(MAKE))
cleanguava:
	-rm guava/bin/*
kbmag: FORCE
	(cd kbmag/lib; $(MAKE))
	(cd kbmag/src; $(MAKE))
cleankbmago:
	-(cd kbmag/lib; rm *.o *.a)
	-rm kbmag/src/*.o
cleankbmag: cleankbmago
	-(cd kbmag/bin; rm gpaxioms gpmakefsa gpminkb kbprog ppgap)
leon: FORCE
	(cd leon/src; $(MAKE) desauto)
	(cd leon/src; $(MAKE) wtdist)
cleanleono:
	-rm leon/src/*.o
cleanleon: cleanleono
	-rm leon/bin/*
meataxe: FORCE
	(cd meataxe/bin; $(MAKE))
cleanmeataxeo:
	(cd meataxe/bin; $(MAKE) clean)
cleanmeataxe: cleanmeataxeo
	(cd meataxe/bin; $(MAKE) delete)
nq: FORCE
	(cd nq/src; $(MAKE))
cleannqo:
	-rm nq/src/*.o
cleannq: cleannqo
	-rm nq/bin/*
pcqa: FORCE
	(cd pcqa/src; $(MAKE))
cleanpcqao:
	(cd pcqa/src; $(MAKE) cleano)
cleanpcqa: cleanpcqao
	-rm pcqa/bin/*
sisyphos: FORCE
	(cd sisyphos/src; $(MAKE) no-flex)
cleansisyphoso:
	(cd sisyphos/src; $(MAKE) cleano)
cleansisyphos: cleansisyphoso
	-rm sisyphos/bin/*
ve: FORCE
	(cd ve; $(MAKE) ibm-i386-linux-gcc)
cleanveo:
	(cd ve; $(MAKE) cleano)
cleanve:
	(cd ve; $(MAKE) clean)
packages: anupq anusq arep cohomolo grape guava leon nq kbmag meataxe\
  sisyphos ve
cleanpackages: cleananupq cleananusq cleanarep cleancohomolo cleangrape\
  cleanguava cleanleon cleannq cleankbmag cleanmeataxe cleanpcqa\
  cleansisyphos cleanve
cleanpackageso: cleananupqo cleananusqo cleancohomoloo cleangrapeo cleanleono\
  cleannqo cleankbmago cleanmeataxeo cleanpcqao cleansisyphoso cleanveo
