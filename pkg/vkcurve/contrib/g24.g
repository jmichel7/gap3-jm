ReadVK("contrib/common"); 
W:=ComplexReflectionGroup(24);
res:=init(W);

if IsBound(VKCURVE.useChevieInvars) then
 invars:=List(Invariants(W),f->ApplyFunc(f,vars));
 res.H:=Hessian(invars[1],varnames);
 invars[2]:=-myDet(res.H)/243;
 invars[3]:=bord(invars[1],invars[2])/145152;
else
 invars:=[x^3*y+z*y^3+x*z^3];
 res.H:=Hessian(invars[1],varnames);
 invars[2]:=myDet(res.H)/54;
 invars[3]:=bord(invars[1],invars[2])/9;
fi;

common24_33(res);

#gap> M;
#[ [ 4x, 6y^2, 14z-288x^2y ], [ 6y, -1/2z, 1024xy^2-224x^4 ], 
#  [ 14z, 1024xy^3-24x^2z-224x^4y, 2296xyz-2240x^3y^2-4704y^4+1792x^6 ] ]
#gap> infty(MinorMat(M));
#[ [ -50176x^8y, 10752x^6y, -1344x^4y^2 ], [ -53760x^6y^2, 7168x^7, -896x^5y],
#  [ -1344x^4y^2, -896x^5, -36y^3 ] ]
# Donc (x,y,z)=(0,0,1) est le seul point possible a l'infini dans l'adherence
# du lieu singulier du discriminant.
# En particulier, la 2-section ci-dessous verifie les conditions
# de tranversalite a l'infini:

res.twoplane:=["x",1,"y",y,"z",x];
res.curve:=Value(res.disc,res.twoplane);

# VKcurve en exact donne au bout de 2h40 sur grobner2 (136 segments)
#
# 1: abab=baba
# 2: acac=caca
# 3: bcbc=cbcb
# 4: abcabca=bcabcab
# 5: bcabcab=cabcabc
#
# (abc)^7 is central
#
# b->abA then b->cbC gives length 36
#
# 1: bab=aba
# 2: bcb=cbc
# 3: acac=caca
# 4: abcBabca=bcBabcab <=> acabcabca=cabcabcab
#
# (abc)^7 is central
#
# b->abA gives length 34
#
# 1: bab=aba
# 2: bcbc=cbcb
# 3: acac=caca
# 4: Abacba=cAbacb <=> acbacba=cacbacb
#
# (acb)^7 is central
