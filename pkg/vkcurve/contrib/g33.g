ReadVK("contrib/common");
W:=ComplexReflectionGroup(33);
res:=init(W);
invars:=List(Invariants(W),f->ApplyFunc(f,vars));
res.H:=Hessian(invars[1],varnames);
common24_33(res);

# Lieu singulier de la partie a l'infini du discriminant:
#gap> infty(myDet(M))*65536/19683;
#-tx^7z^5+1/4ux^8z^4+2x^6yz^6
#Le lieu singulier a trois composantes (les memes que ci-dessous)
#Mineurs de la matrice M
#gap> Concatenation(MinorMat(M));;
#Ensuite, calcul avec le package Groebner de Maple une base 
# de Groebner pour le degre t; on note gro33 la base. 
#gap> infty(gro33);
#[ -4tux^7z+6ux^6yz^2+t^2x^6z^2+u^2x^8, 6tx^6z^3-3ux^7z^2, 
#  -12tux^4yz^3-1296tx^2y^2z^5+6tu^2x^6z+252ux^3y^2z^4+2592xy^3z^6+3t^2ux^5z^2+\
#24t^2x^3yz^4-9u^2x^5yz^2+4t^3x^4z^3-2u^3x^7+184320x^3z^7, 
#  -tux^5z^3-408tx^3yz^5+84ux^4yz^4+864x^2y^2z^6+4t^2x^4z^4+u^2x^6z^2, 
#  -204tx^4z^5+39ux^5z^4+504x^3yz^6, 144x^4z^6, 
#  -5tux^6z^2-12tx^4yz^4+9ux^5yz^3-4t^2x^5z^3+2u^2x^7z+72x^3y^2z^5, 
#  -12tx^5z^4-3ux^6z^3+72x^4yz^5, 72x^5z^5, 9x^6z^4, 
#  -21tx^7z^2+6ux^8z+36x^6yz^3, 9x^7z^3, -78tux^8z+144ux^7yz^2+21u^2x^9, 
#  9ux^8z^2, -72ux^5yz^4+24t^2x^5z^4+3u^2x^7z^2, 63tx^8z^2-18ux^9z, 
#  612tux^3y^2z^4+1944ux^2y^3z^5-45tu^2x^5yz^2-10tu^3x^7-6110208tx^3z^7-60t^2ux\
#^4yz^3+14653440x^2yz^8-3888t^2x^2y^2z^5+30t^2u^2x^6z+15t^3ux^5z^2+120t^3x^3yz^\
#4+20t^4x^4z^3+15552y^4z^7, 9u^2x^9z, 234tux^9z-63u^2x^10, 9u^3x^10, 
#  -756tu^2x^8yz-147t^2u^2x^9+546t^3ux^8z+126u^3x^9y ]
# donc CN de nullite : xz = 0 
#gap> MValue(infty(gro33),["z",0]);      
#[ u^2x^8, 0, -2u^3x^7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21u^2x^9, 0, 0, 0, 
#  -10tu^3x^7, 0, -63u^2x^10, 9u^3x^10, -147t^2u^2x^9+126u^3x^9y ]
#gap> MValue(infty(gro33),["x",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15552y^4z^7, 0, 0, 0, 0 ]
#Donc:
#x<>0 => u=z=0
#z<>0 => x=y=0
#Il y a trois composantes a l'infini:
#gap> MValue(infty(gro33),["x",0,"y",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
#gap> MValue(infty(gro33),["x",0,"z",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
#gap> MValue(infty(gro33),["u",0,"z",0]);
#[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
#Ainsi le deux plan ci-dessous est transverse a l'infini

res.twoplane:=["x",x,"y",y,"z",y,"t",1,"u",x];
res.curve:=Value(res.disc,res.twoplane);

#Transversalite affine: soit N la matrice M evaluee en res.twoplane;
# il s'agit de verifier que la matrice 5x7
#  N11 N12 N13 N14 N15 1 0 
#  N21 N22 N23 N24 N25 0 1
#  N31 N32 N33 N34 N35 0 1 
#  N41 N42 N43 N44 N45 0 0
#  N51 N52 N53 N54 N55 1 0
#  a toujours rang 5.
# de facon equivalente, il faut verifier que [L[1]-L[5],L[2]-L[3],L[4]]
# a toujours rang 3; cela RESTE a verifier automatiquement
######????????? Donc la section est transverse dans la partie affine

# VKcurve en exact donne au bout de 36h sur grobner1 (326 segments)
#
# 1: ad=da
# 2: bc=cb
# 3: ea=ae
# 4: ec=ce
# 5: beb=ebe
# 6: ded=ede
# 7: aca=cac
# 8: bdb=dbd
# 9: cdc=dcd
# 10: aba=bab
# 11: deDbde=bdeDbd
# 12: bdcDab=abdcDa
#
# d->Cdc then d->Ada gives length 68
# 
# 1: db=bd
# 2: cd=dc
# 3: bc=cb
# 4: ea=ae
# 5: ec=ce
# 6: ada=dad
# 7: ded=ede
# 8: beb=ebe
# 9: aca=cac
# 10: aba=bab
# 11: adEbaDABe=badEbDABd
#
# a->Dad gives length 66
# 
# 1: db=bd
# 2: cd=dc
# 3: bc=cb
# 4: ec=ce
# 5: aca=cac
# 6: aba=bab
# 7: ada=dad
# 8: ded=ede
# 9: beb=ebe
# 10: daEd=aEda
# 11: bEaBeb=EaebEa
#
# a->baB gives length 64
# 
# 1: db=bd
# 2: cd=dc
# 3: bc=cb
# 4: ec=ce
# 5: eae=aea
# 6: aca=cac
# 7: ada=dad
# 8: aba=bab
# 9: ded=ede
# 10: beb=ebe
# 11: bedabe=edabed
#
# a->daD then a->eaE gives length 62
# 
# 1: db=bd
# 2: cd=dc
# 3: bc=cb
# 4: ec=ce
# 5: ba=ab
# 6: aea=eae
# 7: aca=cac
# 8: dad=ada
# 9: ded=ede
# 10: beb=ebe
# 11: edaeAd=aedEae
# 
# (abced)^9 is central
