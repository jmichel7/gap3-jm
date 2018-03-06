CC = gcc
CFLAGS= -O2 -m32 -Wno-unused-result
BIN= ../bin
gap:    ${BIN}/gprun ${BIN}/egrun ${BIN}/sylrun ${BIN}/normrun ${BIN}/conrun\
	${BIN}/pcrun ${BIN}/scrun ${BIN}/selgen ${BIN}/matcalc ${BIN}/nqmrun\
	${BIN}/nqrun ${BIN}/grrun ${BIN}/extprun ${BIN}/crrun ${BIN}/readrels
all:  ${BIN}/makegp ${BIN}/gprun ${BIN}/egrun ${BIN}/sylrun ${BIN}/normrun\
	${BIN}/sylnorm ${BIN}/conrun ${BIN}/pcrun ${BIN}/scrun\
	${BIN}/selgen ${BIN}/chbrun ${BIN}/optrun ${BIN}/matcalc\
	${BIN}/matperm ${BIN}/nqmrun ${BIN}/nqrun ${BIN}/grrun\
	${BIN}/extprun ${BIN}/crrun ${BIN}/readrels ${BIN}/testmf\
	${BIN}/testpf ${BIN}/subdir ${BIN}/wreath ${BIN}/permmat\
	${BIN}/tcrun ${BIN}/tcrunb ${BIN}/pkrun ${BIN}/tens ${BIN}/matsum\
	${BIN}/cohomology ${BIN}/execcmd ${BIN}/clrgrp
.c.o:
	${CC} -c $(CFLAGS) $*.c 
${BIN}/makegp: mg.o
	${CC} $(CFLAGS) -o ${BIN}/makegp mg.o
${BIN}/gprun: gpd.o gpp.o permfns.o
	${CC} $(CFLAGS) -o ${BIN}/gprun gpd.o gpp.o permfns.o
${BIN}/egrun: egd.o egp.o permfns.o
	${CC} $(CFLAGS) -o ${BIN}/egrun egd.o egp.o permfns.o
${BIN}/sylrun: syld.o sylp.o morepermfns.o permfns.o
	${CC} $(CFLAGS) -o ${BIN}/sylrun syld.o sylp.o morepermfns.o permfns.o
${BIN}/normrun: normd.o normp1.o normp2.o chb.o permfns.o morepermfns.o
	${CC} $(CFLAGS) -o ${BIN}/normrun normd.o normp1.o normp2.o \
	chb.o permfns.o morepermfns.o
${BIN}/sylnorm: snd.o sylp2.o normp1.o normp2.o permfns.o morepermfns.o chb.o
	${CC} $(CFLAGS) -o ${BIN}/sylnorm snd.o sylp2.o normp1.o normp2.o \
	permfns.o morepermfns.o chb.o
${BIN}/conrun: conpdcd.o conpdc1.o conpdc2.o
	${CC} $(CFLAGS) -o ${BIN}/conrun conpdcd.o conpdc1.o conpdc2.o
${BIN}/pcrun: pcd.o pcp.o permfns.o pcscfns.o
	${CC} $(CFLAGS) -o ${BIN}/pcrun pcd.o pcp.o permfns.o pcscfns.o
${BIN}/scrun: scd.o scp.o permfns.o pcscfns.o
	${CC} $(CFLAGS) -o ${BIN}/scrun scd.o scp.o permfns.o pcscfns.o
${BIN}/selgen: slg.o
	${CC} $(CFLAGS) -o ${BIN}/selgen slg.o
${BIN}/chbrun: testchb.o chb.o permfns.o
	${CC} $(CFLAGS) -o ${BIN}/chbrun testchb.o chb.o permfns.o
${BIN}/optrun: optd.o optp.o permfns.o
	${CC} $(CFLAGS) -o ${BIN}/optrun optd.o optp.o permfns.o
${BIN}/matcalc: mcd.o mcp.o matfns.o
	${CC} $(CFLAGS) -o ${BIN}/matcalc mcd.o mcp.o matfns.o
${BIN}/matperm: matperm.o matfns.o
	${CC} $(CFLAGS) -o ${BIN}/matperm matperm.o matfns.o
${BIN}/nqmrun: nqmd.o
	${CC} $(CFLAGS) -o ${BIN}/nqmrun nqmd.o
${BIN}/nqrun: nqd.o nqp1.o nqp2.o nqp3.o nqfns.o matfns.o
	${CC} $(CFLAGS) -o ${BIN}/nqrun nqd.o nqp1.o nqp2.o nqp3.o nqfns.o \
	matfns.o
${BIN}/grrun: grd.o grp.o permfnsb.o
	${CC} $(CFLAGS) -o ${BIN}/grrun grd.o grp.o permfnsb.o
${BIN}/extprun: extpd.o exa.o exb.o exc.o permfns.o matfns.o
	${CC} $(CFLAGS) -o ${BIN}/extprun extpd.o exa.o exb.o exc.o \
	permfns.o matfns.o
${BIN}/crrun: crd.o crp1.o crp2.o matfns.o
	${CC} $(CFLAGS) -o ${BIN}/crrun crd.o crp1.o crp2.o matfns.o
${BIN}/readrels: readrels.o
	${CC} $(CFLAGS) -o ${BIN}/readrels readrels.o
${BIN}/testmf: testmf.o matfns.o
	${CC} $(CFLAGS) -o ${BIN}/testmf testmf.o matfns.o
${BIN}/testpf: testpf.o permfns.o
	${CC} $(CFLAGS) -o ${BIN}/testpf testpf.o permfns.o
${BIN}/subdir: subdir.o
	${CC} $(CFLAGS) -o ${BIN}/subdir subdir.o
${BIN}/wreath: wreath.o
	${CC} $(CFLAGS) -o ${BIN}/wreath wreath.o
${BIN}/permmat: permmat.o
	${CC} $(CFLAGS) -o ${BIN}/permmat permmat.o
${BIN}/tcrun: tcd.o tcp.o
	${CC} $(CFLAGS) -o ${BIN}/tcrun tcd.o tcp.o
${BIN}/tcrunb: tcdb.o tcpb.o
	${CC} $(CFLAGS) -o ${BIN}/tcrunb tcdb.o tcpb.o
${BIN}/pkrun: pkd.o pkp.o pkmatfns.o
	${CC} $(CFLAGS) -o ${BIN}/pkrun pkd.o pkp.o pkmatfns.o
${BIN}/tens: tens.o
	${CC} $(CFLAGS) -o ${BIN}/tens tens.o
${BIN}/matsum: matsum.o
	${CC} $(CFLAGS) -o ${BIN}/matsum matsum.o
cleano:
	-rm -f *.o
clean:
	-rm ../bin/*run ../bin/matcalc ../bin/selgen ../bin/readrels
