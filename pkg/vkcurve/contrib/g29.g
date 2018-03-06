ReadVK("contrib/common");
W:=ComplexReflectionGroup(29);
res:=init(W);
invars:=List(Invariants(W),f->ApplyFunc(f,vars));
res.H:=Hessian(invars[1],varnames);
#invars:=[phi1,-myDet(H)/20736,f12,f20];
common24_33(res);

#Mineurs de la matrice M
#gap> Concatenation(MinorMat(M));;
#Ensuite, calcul avec le package Groebner de Maple une base 
# de Groebner pour le degre t; on note gro29 la base. 
#gap> infty(gro29)[11];
#x^9y
#gap> MValue(infty(gro29),["x",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#  0, 0, 1381064567813536382754431904ty^11 ]
#gap> MValue(infty(gro29),["y",0]);
#[ -5130994tx^6z^3-32756157120t^2x^7z, 0, 6176tx^7z^2, 0, -2tx^8z, 0, 0, 0, 0, 
#  -16x^9z, 0, 10753689728tx^7z^3+566142649347346t^2x^8z, 
#  -24515545725115258tx^6z^4-77444314158421860480t^2x^7z^2, 0, 
#  821094700862t^2x^8z, 0, 77518tx^8z^2, 110282tx^9z, 0, 
#  3949654148671093390tx^6z^4+14813645146976315113856t^2x^7z^2, 0, 0, 
#  54373x^9z^2, 2426x^10z, 0, 0, 0, 0 ]
#gap> MValue(infty(gro29),["x",0,"y",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#  0, 0, 0 ]
#gap> MValue(infty(gro29),["x",0,"t",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#  0, 0, 0 ]
#gap> MValue(infty(gro29),["y",0,"z",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#  0, 0, 0 ]
#Donc: condition necessaire de nullite xy=0 puis
#    y<>0  =>  x=0 et t=0
#    x<>0  =>  y=0 et z=0
#On verifie qu'il y a exactement trois composantes a l'infini:
#gap> MValue(infty(gro29),["x",0,"y",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#  0, 0, 0 ]
#gap> MValue(infty(gro29),["x",0,"t",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#  0, 0, 0 ]
#gap> MValue(infty(gro29),["y",0,"z",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0 ]
#
#Conclusion: la section suivante verifie les conditions a l'infini

res.twoplane:=["x",x,"y",y,"z",x,"t",y+1];
res.curve:=Value(res.disc,res.twoplane);

#Transversalite affine: soit N la matrice M evaluee en res.twoplane;
# il s'agit de verifier que la matrice 4x6
#  N11 N12 N13 N14 1 0 
#  N21 N22 N23 N24 0 1
#  N31 N32 N33 N34 1 0 
#  N41 N42 N43 N44 0 1
#  a toujours rang 4.
# de facon equivalente, il faut verifier que [L[1]-L[3],L[2]-L[4]]
# a toujours rang 2; cela se verifie automatiquement
#Donc la section est transverse dans la partie affine
#
# VKcurve en exact donne au bout de 2h sur grobner2 (131 segments):
#
# 1: ada=dad
# 2: cac=aca
# 3: bab=aba
# 4: dbd=bdb
# 5: Dcdb=bDcd
# 6: bcbc=cbcb
# 7: cdcd=dcdc
# 8: AcdCab=bAcdCa
# 9: daDCbc=CbcdaD
# 10: cdacda=acdacd
# 11: caCbca=bAcabc
# 12: bCdabcdab=cdabcddaB
# 13: dabcBADcdabc=cdabcBADcdab
#
# b->abA gives length 122
# 
# 1: cbc=bcb
# 2: ada=dad
# 3: cac=aca
# 4: bab=aba
# 5: dCbc=bcBd
# 6: cdcd=dcdc
# 7: aDcdba=baDcdb
# 8: adbaBd=badAba
# 9: dacAda=acAdac
# 10: cdBDac=acdBDa
# 11: cbCacbac=acbacAba
# 12: aBCdbacdb=cdbCacddB
#
# a->Cac gives length 116
#
# 1: dad=ada
# 2: cbc=bcb
# 3: cac=aca
# 4: adba=dbad
# 5: dCbc=bcBd
# 6: baba=abab
# 7: cdcd=dcdc
# 8: caCdca=dcaCdc
# 9: abAcab=cabAca
# 10: dbaBdb=baBdba
# 11: dbcadbb=bcbadba
# 12: adcDACbc=bcBadcDA
#
# d->adA gives length 76
# 
# 1: bd=db
# 2: dad=ada
# 3: cdc=dcd
# 4: cbc=bcb
# 5: cac=aca
# 6: baba=abab
# 7: daDCbc=bcBAda
# 8: Cbcabc=abcBab
# 9: aCAdacda=dacAdaCd
#
# a->daD gives length 58
# 
# 1: bd=db
# 2: bcb=cbc
# 3: dad=ada
# 4: cdc=dcd
# 5: cabc=bcab
# 6: acac=caca
# 7: baba=abab
# 8: dcadca=cadcad
#
# a->caC gives length 46
# 
# 1: bd=db
# 2: ab=ba
# 3: bcb=cbc
# 4: cdc=dcd
# 5: dad=ada
# 6: acac=caca
# 7: acdacd=dacdac
#
# (abcd)^5 is central
