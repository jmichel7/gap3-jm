ReadVK("contrib/common");
W:=ComplexReflectionGroup(31);
res:=init(W);
invars:=List(Invariants(W),f->ApplyFunc(f,vars));
res.H:=Hessian(invars[1],varnames);
common24_33(res);

#Mineurs de la matrice M
#gap> Concatenation(MinorMat(M));;
#Ensuite, calcul avec le package Groebner de Maple une base 
# de Groebner pour le degre t; on note gro31 la base. 
#gap> infty(gro31);
#[ tx^6, 36tx^2y^2z^2+24xy^3z^3+2t^2x^3yz-t^3x^4-4320x^3z^4, 
#  34tx^3yz^2+24x^2y^2z^3+t^2x^4z, 17tx^4z^2+14x^3yz^3, 4x^4z^3, 
#  tx^4yz-t^2x^5+2x^3y^2z^2, tx^5z+2x^4yz^2, 2x^5z^2, x^6z, 2t^2x^5z, 
# -143208tx^3z^4-114480x^2yz^5-108t^2x^2y^2z^2-10t^3x^3yz+5t^4x^4+48y^4z^4 ]
# Donc tx=0 CN de nullite
# x<> 0  => z =0 car on a 
#gap> MValue(infty(gro31),["t",0]);
#[ 0, 24xy^3z^3-4320x^3z^4, 24x^2y^2z^3, 14x^3yz^3, 4x^4z^3, 2x^3y^2z^2, 
#  2x^4yz^2, 2x^5z^2, x^6z, 0, -114480x^2yz^5+48y^4z^4 ]
####
# t <>  => yz=0 car on a
#gap> MValue(infty(gro31),["x",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 48y^4z^4 ]
# On en deduit les trois composantes a l'infini:
#gap> MValue(infty(gro31),["x",0,"y",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
#gap> MValue(infty(gro31),["x",0,"z",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
#gap> MValue(infty(gro31),["t",0,"z",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
# Donc la section suivante est transverse a l'infini:

res.twoplane:=["x",x,"y",y,"z",y,"t",x+1];
res.curve:=Value(res.disc,res.twoplane);

#Transversalite affine: soit N la matrice M evaluee en res.twoplane;
# il s'agit de verifier que la matrice 4x6
#  N11 N12 N13 N14 1 0 
#  N21 N22 N23 N24 0 1
#  N31 N32 N33 N34 0 1 
#  N41 N42 N43 N44 1 0
#  a toujours rang 4.
# de facon equivalente, il faut verifier que [L[1]-L[4],L[2]-L[3]]
# a toujours rang 2; cela se verifie automatiquement
# Donc la section est transverse dans la partie affine

# VKcurve en exact donne au bout de 4h sur grobner2 (145 segments)
#
# 1: eb=be
# 2: ba=ab
# 3: ca=ac
# 4: cde=dec
# 5: dec=ecd
# 6: bdb=dbd
# 7: cbc=bcb
# 8: aea=eae
# 9: dad=ada
#
# (abcde)^6 is central
