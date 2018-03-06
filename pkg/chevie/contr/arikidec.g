#############################################################################
##
##  arikidec.g                      Nicolas Jacon (University of Caen LMNO)
##
##
## This  file contains  functions for  computing the  canonical basis of an
## arbitrary  irreducible integrable  highest weight  representation of the
## quantum  group of the affine  special linear group U_v (\widehat{sl}_e).
## It also computes the decomposition matrix of Ariki-Koike algebras, where
## the  parameters  are  power  of  a  e-th  root  of  unity  in a field of
## characteristic  zero. The algorithm and the notations are taken from the
## paper  "An  algorithm  for  the  computation  of  the simple modules for
## Ariki-Koike algebras" (see http://arxiv.org/abs/math.RT/0409271).
##
##
## Let  l, n and e be  three positive integers and let  0 \leq v_0 \leq v_1
## \leq   ...  \leq   v_{l-1}<e  be   nonnegative  integers.  The  function
## "CrystalMatrixDecomposition(parametre,e,l,n)"   computes  the  canonical
## basis elements of rank $n$ of the irreducible $U_v
## (\widehat{sl}_e)$-module of highest weight
## $\Lambda_{v_0}+\Lambda_{v_1}+...+\Lambda_{v_{l-1}}$. Here "parametre" is
## the list [v_0,...,v_{l-1}].
##
## Example:
##  gap> CrystalMatrixDecomposition([0,1],2,2,4);
## [ [     4.,      |,      1,      0,      0,      0,      0,      0 ],
##   [     .4,      |,      0,      1,      0,      0,      0,      0 ],
##   [    3.1,      |,      v,      v,      1,      0,      0,      0 ],
##   [    1.3,      |,      v,      v,      0,      1,      0,      0 ],
##   [    .31,      |,      0,    v^2,      0,      v,      0,      0 ],
##   [    2.2,      |,    v^2,    v^2,      v,      v,      0,      0 ],
##   [    31.,      |,    v^2,      0,      v,      0,      0,      0 ],
##   [   2.11,      |,      v,    v^3,    v^2,    v^2,      0,      0 ],
##   [   1.21,      |,      0,      0,      0,      0,      1,      0 ],
##   [    .22,      |,      0,      0,      0,    v^2,      0,      0 ],
##   [   11.2,      |,    v^3,      v,    v^2,    v^2,      0,      0 ],
##   [   21.1,      |,      0,      0,      0,      0,      0,      1 ],
##   [    22.,      |,      0,      0,    v^2,      0,      0,      0 ],
##   [   211.,      |,    v^2,      0,    v^3,      0,      0,      0 ],
##   [  11.11,      |,    v^2,    v^2,    v^3,    v^3,      0,      0 ],
##   [   .211,      |,      0,    v^2,      0,    v^3,      0,      0 ],
##   [  111.1,      |,    v^3,    v^3,    v^4,      0,      0,      0 ],
##   [  1.111,      |,    v^3,    v^3,      0,    v^4,      0,      0 ],
##   [  1111.,      |,    v^4,      0,      0,      0,      0,      0 ],
##   [  .1111,      |,      0,    v^4,      0,      0,      0,      0 ] ]
##
## The  rows of the matrix corresponds to the canonical basis elements. The
## lines  are labeled by the set of  l-partitions of rank $n$ (in the above
## example,  1.111 stands for the  bipartition ((1),(1,1,1))). The function
## MatrixDecomposition(parametre,e,l,n)computes the decomposition matrix of
## the Ariki-Koike algebras with set of parameters
## {\eta_e^{v_0},\eta_e^{v_1},...,\eta_e^{v_{l-1}}}.  By  Ariki's  theorem,
## this  matrix  is  just  the  specialization  at  v=1 of the crystallized
## decomposition matrix.
##
## Example:
## gap> MatrixDecomposition([0,1],2,2,4);
##      4.   |  1.....
##      .4   |  .1....
##     3.1   |  111...
##     1.3   |  11.1..
##     .31   |  .1.1..
##     2.2   |  1111..
##     31.   |  1.1...
##    2.11   |  1111..
##    1.21   |  ....1.
##     .22   |  ...1..
##    11.2   |  1111..
##    21.1   |  .....1
##     22.   |  ..1...
##    211.   |  1.1...
##   11.11   |  1111..
##    .211   |  .1.1..
##   111.1   |  111...
##   1.111   |  11.1..
##   1111.   |  1.....
##   .1111   |  .1....
##
###############################################################################
###############################################################################

RequirePackage("chevie");

#################
#FIRST PART     #
#################
## In this part, the main function gives the set of FLOTW $d$-partitions of
## rank n with respect to the choice of parameters.

DiagramMultiPartition:=function(multipartition,l)local d,i,j,k;
  d:=[];
  for k in [1..l] do
   for i in [1..Length(multipartition[k])] do
     for j in [1..multipartition[k][i]] do
      Add(d,[i,j,k]);
     od;
   od;
  od;
  return Set(d);
end;

#EXAMPLE :
#gap> DiagramMultiPartition([[2,1],[1]],2);
#[ [ 1, 1, 1 ], [ 1, 1, 2 ], [ 1, 2, 1 ], [ 2, 1, 1 ] ]

DiagramMultiPartitionfrontiere:=function(multipartition,l)local d,i,k;
  d:=[];
    for k in [1..l] do
      for i in [1..Length(multipartition[k])] do
        Add(d,[i,multipartition[k][i] ,k]);
      od;
    od;
return Set(d);
end;

#EXAMPLE :
#gap> DiagramMultiPartitionfrontiere([[2,1],[1]],2);
#[ [ 1, 1, 2 ], [ 1, 2, 1 ], [ 2, 1, 1 ] ]

ResidueDiagram:=function(l,parametre,q,boite)
  return q^(boite[2]-boite[1]+parametre[boite[3]]);
end;

#EXAMPLE :
#gap> ResidueDiagram(2,[0,1],E(4), [ 1, 2, 1 ]);
#E(4)

Place:=function(liste,j)
 if j>Length(liste) then  return 0;else return liste[j];fi;
end;

#EXAMPLE :
#gap> Place([2,1],3);
#0
#gap> Place([2,1],1);
#2

Condition1:=function(l,parametre,multipartition,e)
 local j,i,k;
   j:=1;
    while j<l do
     i:=1;
        while i<Length(multipartition[j+1])+1 do
           if Place(multipartition[j],i)<Place(multipartition[j+1],i+
                                                 parametre[j+1]-parametre[j])
           then i:=Length(multipartition[j+1])+2;
           else  i:=i+1;
           fi;
        od;
     if i=Length(multipartition[j+1])+2 then j:=l+2;
     else j:=j+1;
     fi;
    od;
 if j=l then
   k:=1;
       while k<Length(multipartition[1])+1 do
           if Place(multipartition[j],k)<Place(multipartition[1],k+
                                               parametre[1]-parametre[j]+e)
           then   return 0;
           else  k:=k+1;
           fi;
      od;
 fi;
 if j=l+2 then  return 0;
    else return 1;
    fi;
end;

#EXAMPLE :
#gap> Condition1(2,[0,1],[[2,1],[1]],2);
#1
#gap> Condition1(2,[0,1],[[2,1,1],[1]],2);
#0

Ensembleresidus:=function(l,parametre,multipartition,q,a)local i,d,diag;
i:=1;
d:=[];
diag:=DiagramMultiPartitionfrontiere(multipartition,l);
for i in [1..Length(diag)] do
   if diag[i][2]=a then Add(d,ResidueDiagram(l,parametre,q,diag[i]));
   fi;
i:=i+1;
od;
return Set(d);
end;

#EXAMPLE :
#gap> Ensembleresidus(2,[0,1],[[2,1],[1,1]],E(4),1);
#[ 1, -E(4), E(4) ]

Grandmax:=function(l,multipartition)local i,d;
  d:=0; i:=1;
  for i in [1..l] do d:=Maximum(d,Place(multipartition[i],1)); od;
  return d;
end;

#gap> Grandmax(2,[[2,1],[5,2,3]]);
#5

Condition2:=function(l,parametre,multipartition,e)local i,j,p,q;
  q:=E(e); p:=Grandmax(l,multipartition); i:=1;
  while i in [1..p] do
    j:=Ensembleresidus(l,parametre,multipartition,q,i);
    if Length(j)=e  then i:=p+2;
    else i:=i+1;
    fi;
  od;
  if i=p+1 then return 1;
  else return 0;
  fi;
end;

#EXAMPLE :
#gap> Condition2(2,[0,1],[[2,1,1],[1]],2);
#0
#gap> Condition2(2,[0,1],[[2,1],[1]],2);
#1

# The following function gives the list of FLOTW multipartitions with respect to the list
# parametre:=[v_0,...,v_{l-1], and to q, a primitive e-root of unity.

FLOTW:=function(l,parametre,n,e)local i,ensemble,d;
  ensemble:=PartitionTuples(n,l); d:=[]; i:=1;
  for i in [1..Length(ensemble)] do
    if   Condition1(l,parametre,ensemble[i],e)=0 then i:=i+1;
    elif Condition2(l,parametre,ensemble[i],e)=0 then i:=i+1;
    else Add(d,ensemble[i]);i:=i+1;
    fi;
  od;
  return Set(d);
end;

#EXAMPLE :
#gap> FLOTW(2,[0,1],6,2);
#[ [ [  ], [ 6 ] ], [ [ 1 ], [ 4, 1 ] ], [ [ 1 ], [ 5 ] ], [ [ 2 ], [ 3, 1 ] ],
#  [ [ 2 ], [ 4 ] ], [ [ 2, 1 ], [ 3 ] ], [ [ 3 ], [ 2, 1 ] ],
#  [ [ 3, 1 ], [ 2 ] ], [ [ 4 ], [ 2 ] ], [ [ 4, 1 ], [ 1 ] ],
#  [ [ 5 ], [ 1 ] ], [ [ 6 ], [  ] ] ]

#################
#SECOND PART    #
#################
##

# The aim of this part is to give a function which computes the $a$-sequence of residues of a FLOTW
# multipartition.

Candidat:=function(l,q,multipartition)local i,j,d,p,e;
 i:=1; d:=[]; e:=[];
 p:=Grandmax(l,multipartition);
 while i<l+1 do
  j:=1;
   if Place(multipartition[i],1)<p then i:=i+1;
     else while multipartition[i][j]=Place(multipartition[i],j+1)
                  do j:=j+1;
                     Add(e,[j-1, multipartition[i][j-1],i]);
                  od;
         Add(d,[j, multipartition[i][j],i]);
           Add(e,[j, multipartition[i][j],i]);
           i:=i+1;
   fi;
 od;
return [d,e];
end;

#EXAMPLE :
#gap> Candidat(4,E(4),[[3,2],[3,1],[3,3,2],[2,1]]);
#[ [ [ 1, 3, 1 ], [ 1, 3, 2 ], [ 2, 3, 3 ] ],
#  [ [ 1, 3, 1 ], [ 1, 3, 2 ], [ 1, 3, 3 ], [ 2, 3, 3 ] ] ]

Boiteoptimale:=function(l,q,multipartition,parametre)local i,j,d,p,r,s;
 d:=Candidat(l,q,multipartition);
 p:=Length(d[1]); s:=Length(d[2]); r:=0; i:=1;
 while i<p+1 do
    j:=1;
    while j<s+1 do
      if  ResidueDiagram(l,parametre,q,d[1][i])=
          ResidueDiagram(l,parametre,q,d[2][j])*q
      then j:=s+2;
      else j:=j+1;
      fi;
    od;
    if j=s+1 then r:=d[1][i];i:=p+2;
      else i:=i+1;
    fi;
 od;
return r;
end;

#EXAMPLE
#gap> Boiteoptimale(2,E(4),[[2],[2]],[0,3]);
#[ 1, 2, 2 ]

Boitefront:=function(l,q,multipartition,a)local i,d,k;
  d:=[];
    for k in [1..l] do
      for i in [1..Length(multipartition[k])] do
         if multipartition[k][i]=a then Add(d,[i,multipartition[k][i] ,k]);
         fi;
      od;
    od;
return Set(d);
end;

#EXAMPLE :
#gap> Boitefront(3,E(4),[[3,2,2],[2,1],[6,5,2]],2);
#[ [ 1, 2, 2 ], [ 2, 2, 1 ], [ 3, 2, 1 ], [ 3, 2, 3 ] ]

Supprime:=function(l,q,multipartition,res,a,parametre)local i,d,k,e,mul;
mul:=Copy(multipartition);
d:=0;
e:=DiagramMultiPartition(multipartition,l);
k:=Boitefront(l,q,multipartition,a);
for i in [1..Length(k)] do
 if ResidueDiagram(l,parametre,q,k[i])=res then
            if
              a>Place(multipartition[k[i][3]],[k[i][1]]+1)
                                then mul[k[i][3]][k[i][1]]:=a-1;d:=d+1;
            fi;
 fi;
od;
return [mul,d];
end;

#EXEMPLE:
#gap> Supprime(2,E(4),[[2],[2,2,2]],E(4),2,[0,2]);
#[ [ [ 1 ], [ 2, 2, 1 ] ], 2 ]

asq:=function(l,q,parametre,multipartition,res)local k,p,i,mul,d,a,s;
d:=0;
mul:=Copy(multipartition);
p:=Boiteoptimale(l,q,multipartition,parametre);
a:=p[2];
while a>0 do
  k:=Boitefront(l,q,multipartition,a);
     for i in [1..Length(k)] do
      if ResidueDiagram(l,parametre,q,k[i])*q=res then return  [mul,res,d];
      fi;
     od;
s:=Supprime(l,q,mul,res,a,parametre);
mul:=s[1];d:=d+s[2];a:=a-1;od;
return [mul,res,d];
end;

#EXEMPLE :
#gap> asq(2,E(4),[0,2],[[4,3,2],[4,4,4]],E(4)^3);
#[ [ [ 3, 3, 1 ], [ 4, 4, 3 ] ], -E(4), 3 ]

# The following function computes the a-sequence of residues of a FLOTW multipartition "multipartition"
# of rank n with respect to the set of parameters "parametre" and q a primitive root of unity. The result
# is the list [[a1,b1],[a2,b2]...[an,bn]] if the a-sequence is a1,....,a1 (b1 times) a2,....,a2 (b2 times)
# ...an...an (an times).

asequence:=function(q,parametre,multipartition,n)local l,res,mul,i,j,d,p,e;
l:=Length(multipartition);mul:=Copy(multipartition);
d:=0; e:=0; j:=[];
while e<n do
 p:=Boiteoptimale(l,q,mul,parametre);
 res:=ResidueDiagram(l,parametre,q,p);
 i:=asq(l,q,parametre,mul,res);
 mul:=i[1];d:=i[3];Add(j,[res,d]);e:=e+i[3];
od;
return [j,mul];
end;

#EXEMPLE :
#gap> asequence(E(3),[0,2],[[3,2],[2,1]],8);
#[ [ [ E(3)^2, 1 ], [ 1, 2 ], [ E(3), 2 ], [ E(3)^2, 2 ], [ 1, 1 ] ],
#  [ [ 0, 0 ], [ 0, 0 ] ] ]

ordreflotw:=function(boite1,boite2,parametre)
  if boite1[2]-boite1[1]+parametre[boite1[3]]<
     boite2[2]-boite2[1]+parametre[boite2[3]]
  then return boite1;
  else if boite1[2]-boite1[1]+parametre[boite1[3]]=
          boite2[2]-boite2[1]+parametre[boite2[3]] 
       then if boite1[3]>boite2[3] then return boite1;else return boite2;fi;
       else return boite2;
       fi;
  fi;
end;

#EXEMPLE :
#gap> ordreflotw([2,1,3],[3,2,2],[0,2,3]);
#[ 3, 2, 2 ]

##################
#THIRD PART      #
##################
## The main function of this part gives the $a$-value of a multipartition "multipartition" of rank n
## with respect to the choice of parameters "parametre" and q a primitive e-th root of unity.

symbole:=function(l,multipartition,h)local i,j,mul;
  mul:=Copy(multipartition);
   for i in [1..l] do
      for j in [1..h] do
         mul[i][j]:=Place(multipartition[i],j)-j+h;
      od;
   od;
   return mul;
end;

#EXEMPLE :
#gap> symbole(2,[[2,1],[1]],2);
#[ [ 3, 1 ], [ 2, 0 ] ]

tau:=function(l,h)local i,k; i:=1; k:=0;
  while l*(h-i)+1>=2 do k:=NrCombinations([1..l*(h-i)+1],2)+k;i:=i+1;od;
  return k;
end;

afonction1:=function(l,multipartition,n,h,parametre,e)
local i,j,k,p,s,c,B,a1,a2,a3,a4,v,w,x,y,z,i1,i2,B1,afonction,poids,t;
poids:=[];
for t in [1..Length(parametre)] do poids[t]:=parametre[t]-(t-1)*e/l+e; od;
B:=symbole(l,multipartition,h);
B1:=0;
for i1 in [1..l] do for i2 in [1..h] do B1:=B1+B[i1][i2]; od; od;

a1:=0;
for i in [1..l] do a1:=a1+poids[i]; od;
a2:=0;
for j in [1..l] do
 for k in [j+1..l] do a2:=a2+Minimum(poids[k],poids[j]); od;
od;
a3:=0;
for p in [1..l] do for s in [p..l] do
    for c in [1..h] do for v in [1..h] do
         if s=p then if B[p][c]>B[s][v] then a3:=a3+
                        Minimum(B[p][c]+poids[p],B[s][v]+poids[s]);fi;
         else a3:=a3+Minimum(B[p][c]+poids[p],B[s][v]+poids[s]);fi;
      od;
    od;
  od;
od;
a4:=0;
for w in [1..l] do for x in [1..l] do
   for z in [1..h] do
    for y in [1..B[x][z]] do a4:=a4+Minimum(y+poids[x],poids[w]); od;
   od;
 od;
od;

return n*a1-tau(l,h)+B1-n-h*a2+a3-a4;
end;

# a-value of a multipartition:

afonction2:=function(l,multipartition,n,parametre,e)local i,p;
p:=0;
for i in [1..Length(multipartition)] do
   p:=Maximum(p,Length(multipartition[i]));
od;
return afonction1(l,multipartition,n,p,parametre,e) ;
end;

#EXEMPLE :
#gap> afonction2(2,[[2,1],[1]],4,[0,2],4);
#3

##################
#FOURTH PART     #
##################
## Finally,  in this part, we compute  the canonical basis elements and the
## decomposition  matrix of Ariki-Koike algebras as  it is explained in the
## introduction.

Addable:=function(l,parametre,q,mu,res)local i,j,liste,d;
liste:=[];
for i in [1..Length(mu)] do
  if ResidueDiagram(l,parametre,q,[1,Place(mu[i],1)+1,i])=res
  then Add(liste,[1,Place(mu[i],1)+1,i]);
  fi;
  j:=2;
  while j in  [2..Length(mu[i])+1] do
  if  Place(mu[i],j)=Place(mu[i],j-1) then j:=j+1;
  elif ResidueDiagram(l,parametre,q,[j,Place(mu[i],j)+1,i])=res
  then Add(liste,[j,Place(mu[i],j)+1,i]);j:=j+1;
  else j:=j+1;
  fi;od;
od;
return liste;
end;

#EXEMPLE :
#gap> Addable(2,[0,1],-1,[[2,1],[1]],1);
#[ [ 1, 3, 1 ], [ 2, 2, 1 ], [ 3, 1, 1 ], [ 1, 2, 2 ], [ 2, 1, 2 ] ]
#gap> Addable(2,[0,1],-1,[[2,1],[1]],-1);
#[  ]

Addboite:=function(l,mu,boite)local mul; mul:=Copy(mu);
  mul[boite[3]][boite[1]]:=Place(mul[boite[3]],boite[1])+1;
  return mul;
end;

#EXEMPLE :
#gap> Addboite(2,[[2,1],[1]],[1,3,1]);
#[ [ 3, 1 ], [ 1 ] ]

Removable:=function(l,parametre,q,mu,res)local i,j,liste;
liste:=[]; i:=1;
while  i<=Length(mu) do
j:=1;
 while  j<Length(mu[i])+1 do
     if mu[i][j]>Place(mu[i],j+1) then
       if ResidueDiagram(l,parametre,q,[j,Place(mu[i],j),i])=res
       then Add(liste,[j,Place(mu[i],j),i]);
       fi;
     fi;
j:=j+1;
od;
i:=i+1;
od;
return liste;
end;

#EXEMPLE
#gap> Removable(2,[0,1],-1,[[2,1],[1]],-1);
#[ [ 1, 2, 1 ], [ 2, 1, 1 ], [ 1, 1, 2 ] ]

Nib:=function(d,lambda,mu,parametre,q,res,liste)
local i,j,k,baj,bre,resu,boite,l;
k:=Length(liste);
resu:=0;
baj:=Addable(d,parametre,q,mu,res);
bre:=Removable(d,parametre,q,lambda,res);
for i in [1..k] do
 boite:=liste[i];
  for j in [1..Length(baj)] do
    if ordreflotw(boite,baj[j],parametre)=boite then resu:=resu+1;fi;od;
  for l in [1..Length(bre)] do
     if ordreflotw(boite,bre[l],parametre)=boite then resu:=resu-1;fi;od;
 od;
return resu;
end;

#EXEMPLE :
#gap> Nib(3,[[3,2],[1,1,1],[5,4,1]],[[3,2,1],[1,1,1],[5,4,1]],[1,1,2],E(3),E(3))
#3

v:=X(Rationals);v.name:="v";

Addplus:=function(l,mu,listeboite,a1)local j,nu;
  j:=a1; nu:=Copy(mu);
  while j>0 do nu:=[Addboite(l,nu[1],listeboite[j])]; j:=j-1; od;
  return nu;
end;

#EXEMPLE
#gap> Addplus(2,[[[2,1],[1]]],[[1,3,1],[2,2,1]],2);
#[ [ [ 3, 2 ], [ 1 ] ] ]

iinduction:=function(l,parametre,i1,a1,lambda,coeff,v,q)
local i,liste2,liste,li,final,j;
li:=[]; liste2:=Addable(l,parametre,q,lambda,i1);
if a1>Length(liste2) then return [];fi;
liste:=Combinations(liste2,a1);
for i in [1..Length(liste)] do
  li:=Concatenation(li,Addplus(l,[lambda],liste[i],a1));
od;
final:=[];
for j in [1..Length(li)] do
  Add(final,[coeff*v^Nib(l,lambda,li[j],parametre,q,i1,liste[j]),li[j]]);
od;
return final;
end;

#EXEMPLE
#gap> iinduction(2,[0,1],-1,2,[[1],[]],1,v,-1);
#[ [ v^0, [ [ 2 ], [ 1 ] ] ], [ v^2, [ [ 1, 1 ], [ 1 ] ] ],
#  [ v, [ [ 2, 1 ], [  ] ] ] ]

simplification:=function(liste)local i,j,p,z,lis,li2;
  lis:=[]; i:=1; li2:=Copy(liste);
  if Length(liste)=0 then return [];fi;
  while i<Length(liste)+1 do
    z:=li2[i][1]; p:=li2[i][2];
    if p=0 then i:=i+1;
    else for j in [i+1..Length(liste)] do
        if liste[j][2]=p then z:=z+liste[j][1];li2[j][2]:=0;fi;
      od;
      Add(lis,[z,p]);
      i:=i+1;
    fi;
  od;
  return lis;
end;

#

inductioneten:=function(l,parametre,i1,a1,liste,v,q)local i,n,li,resu;
li:=[]; n:=Length(liste);
for i in [1..n] do
 li:=Concatenation(li,iinduction(l,parametre,i1,a1,liste[i][2],liste[i][1],v,q));
 od;
resu:=simplification(li);
return resu;
end;

#EXEMPLE
#gap> inductioneten(3,[0,0,0],1,1,[[1,[[1],[],[]]],[v,[[],[1],[]]]],v,-1);
#[ [ v + v^(-1), [ [ 1 ], [ 1 ], [  ] ] ], [ v^0, [ [ 1 ], [  ], [ 1 ] ] ],
#  [ v, [ [  ], [ 1 ], [ 1 ] ] ] ]

aliste:=function(multiflotw,parametre,q,l,v,n)local i,aseq,provi;
  aseq:=Reversed(asequence(q,parametre,multiflotw,n)[1]);
  provi:=[[1,List([1..l],i->[])]];
  for i in [1..Length(aseq)] do
    provi:=inductioneten(l,parametre,aseq[i][1],aseq[i][2],provi,v,q);
  od;
  return provi;
end;

#gap> aliste([[2,2],[2,2,1]],[0,1],E(4),2,v,9);
#[ [ v^0, [ [ 2, 2 ], [ 2, 2, 1 ] ] ], [ v, [ [ 2, 1 ], [ 2, 2, 2 ] ] ],
#  [ v, [ [ 2, 2, 1, 1 ], [ 1, 1, 1 ] ] ],
#  [ v^2, [ [ 2, 1, 1, 1, 1 ], [ 1, 1, 1 ] ] ],
#  [ v^2, [ [ 2, 2 ], [ 1, 1, 1, 1, 1 ] ] ],
#  [ v^3, [ [ 2, 1 ], [ 1, 1, 1, 1, 1, 1 ] ] ] ]

FLOTW2:=function(l,parametre,n,e) local k,liste,i;
  liste:=FLOTW(l,parametre,n,e);
  k:=[];
  for i in [1..Length(liste)] do
    Add(k,[liste[i],afonction2(l,liste[i],n,parametre,e)]);
  od;
  SortBy(k,x->x[2]);
  return k;
end;

#EXEMPLE
#gap> FLOTW2(2,[0,1],4,2);
#[ [ [ [  ], [ 4 ] ], 0 ], [ [ [ 4 ], [  ] ], 0 ], [ [ [ 1 ], [ 3 ] ], 1 ],
#  [ [ [ 3 ], [ 1 ] ], 1 ], [ [ [ 1 ], [ 2, 1 ] ], 3 ],
#  [ [ [ 2, 1 ], [ 1 ] ], 3 ] ]

aliste2:=function(multiflotw,parametre,e,l,v,n)local al,i;
  al:=aliste(multiflotw,parametre,E(e),l,v,n);
  for i in [1..Length(al)] do
    Add(al[i],afonction2(l,al[i][2],n,parametre,e));
  od;
  return al;
end;

#EXEMPLE
#gap> aliste2([[2,1],[1]],[0,1],2,2,v,4);
#[ [ v^0, [ [ 2, 1 ], [ 1 ] ], 3 ] ]

ordonaliste2:=function(multiflotw,parametre,e,l,v,n)local k;
  k:=aliste2(multiflotw,parametre,e,l,v,n);
  SortBy(k,x->x[3]);
  return k;
end;

#gap> ordonaliste2([[2],[]],[0,1],2,2,v,2);
#[ [ v^0, [ [ 2 ], [  ] ], 0 ], [ v, [ [ 1 ], [ 1 ] ], 1 ],
#  [ v^2, [ [ 1, 1 ], [  ] ], 2 ] ]

multi:=function(liste,coeff,v)local i,li;
  li:=Copy(liste);
  for i in [1..Length(liste)] do
    li[i][1]:=coeff*(li[i][1]);
    od;
  return li;
end;

#EXEMPLE :
#gap> multi([ [ v^0, [ [ 2 ], [  ] ], 0 ], [ v, [ [ 1 ], [ 1 ] ], 1 ],
#  [ v^2, [ [ 1, 1 ], [  ] ], 2 ] ]>   [ v^2, [ [ 1, 1 ], [  ] ], 2 ] ],v,v);
#[ [ v, [ [ 2 ], [  ] ], 0 ], [ v^2, [ [ 1 ], [ 1 ] ], 1 ],
#  [ v^3, [ [ 1, 1 ], [  ] ], 2 ] ]

inversete:=function(g,v)local f,h;
f:=Copy(g); h:=Copy(g);
if f=0*v^0 then return f;fi; #peut-etre pb ici
if Degree(f)=0 then f:=f-LeadingCoefficient(f)*v^(Degree(f));fi;
  while f <> 0*v^0 and f<>0 do
    h:=h+LeadingCoefficient(f)*v^(-Degree(f));
    f:=f-LeadingCoefficient(f)*v^(Degree(f));
  od;
#h:=h+g;
return h;
end;

barinv:=function(f,v)local g;g:=f;
  while Degree(g)>0 do g:=g-LeadingCoefficient(g)*v^(Degree(g)); od;
  g:=inversete(g,v);
  return g;
end;

#gap> barinv(v^4+2*v^0+v^(-2),v);
#v^2 + 2 + v^(-2)

Esss:=function(parametre,e,l,v,n)local Ensmulti,i;
  Ensmulti:=FLOTW2(l,parametre,n,e);
  for i in [1..Length(Ensmulti)] do
    Add(Ensmulti[i],ordonaliste2(Ensmulti[i][1],parametre,e,l,v,n));
  od;
  return Ensmulti;
end;

simplification2:=function(liste)local i,j,p,z,lis,li2;
lis:=[]; i:=1; li2:=Copy(liste);
if Length(liste)=0 then return [];
else
while i<Length(liste)+1 do
   z:=li2[i][1];
   p:=li2[i][2];
   if p=0 then i:=i+1;
    else
      for j in [i+1..Length(liste)] do
        if liste[j][2]=p then z:=z+liste[j][1];li2[j][2]:=0; fi;
      od;
     Add(lis,[z,p,li2[i][3]]);
     i:=i+1;
  fi;
 od;return lis;
 fi; end;

# the canonical basis elements

Basecanonique:=function(parametre,e,l,v,n,Ensmulti)
  local i,coeff,k,loc,j,z,elem,r;
  for i in [1..Length(Ensmulti)] do
    Add(Ensmulti[i],ordonaliste2(Ensmulti[i][1],parametre,e,l,v,n));
  od;
  i:=Length(Ensmulti)-1;
  while  i>0  do
  k:=2;
  loc:=Ensmulti[i][3];
    while k<=Length(loc) do
      if Degree(loc[k][1])-EuclideanDegree(loc[k][1])<0  then
	 z:=loc[k][2];
	 coeff:=barinv(loc[k][1],v);
	 r:=i;            #on cherche la bip
	 while Ensmulti[r][1]<>z do
	   if r=Length(Ensmulti) then
	     return [Ensmulti,z,Length(Ensmulti),0,"erreur"];fi;
	   r:=r+1;
	 od;
	 elem:=Ensmulti[r][3];
	 loc:=simplification2(Concatenation(loc,multi(elem,(-coeff),v)));
	 k:=k+1;
       elif  Value(loc[k][1],0)=0 then k:=k+1;
       else z:=loc[k][2];
	 coeff:=barinv(loc[k][1],v);
	 r:=i+1;
	 while Ensmulti[r][1]<>z do
	   if r=Length(Ensmulti) then return [Ensmulti,z,r-i,0,"erreur"];fi;
	 r:=r+1;
	 od;
	 elem:=Ensmulti[r][3];
	 loc:=simplification2(Concatenation(loc,multi(elem,(-coeff),v)));
	 k:=k+1;
       fi;
    od;
    Ensmulti[i][3]:=loc;
    i:=i-1;
  od;
  return Ensmulti;
end;

#EXEMPLE:
#gap> Basecanonique([0,1],2,2,v,4,FLOTW2(2,[0,1],4,2));
#[ [ [ [  ], [ 4 ] ], 0,
#      [ [ v^0, [ [  ], [ 4 ] ], 0 ], [ v, [ [ 3 ], [ 1 ] ], 1 ],
#          [ v, [ [ 1 ], [ 3 ] ], 1 ], [ v^2, [ [ 2 ], [ 2 ] ], 2 ],
#          [ v^2, [ [  ], [ 3, 1 ] ], 2 ], [ v^3, [ [ 2 ], [ 1, 1 ] ], 3 ],
#          [ v, [ [ 1, 1 ], [ 2 ] ], 3 ], [ v^2, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ v^2, [ [  ], [ 2, 1, 1 ] ], 6 ],
#          [ v^3, [ [ 1, 1, 1 ], [ 1 ] ], 7 ],
#          [ v^3, [ [ 1 ], [ 1, 1, 1 ] ], 7 ],
#          [ v^4, [ [  ], [ 1, 1, 1, 1 ] ], 12 ] ] ],
#  [ [ [ 4 ], [  ] ], 0, [ [ v^0, [ [ 4 ], [  ] ], 0 ],
#          [ v, [ [ 1 ], [ 3 ] ], 1 ], [ v, [ [ 3 ], [ 1 ] ], 1 ],
#          [ v^2, [ [ 2 ], [ 2 ] ], 2 ], [ v^2, [ [ 3, 1 ], [  ] ], 2 ],
#          [ v^3, [ [ 1, 1 ], [ 2 ] ], 3 ], [ v, [ [ 2 ], [ 1, 1 ] ], 3 ],
#          [ v^2, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ v^2, [ [ 2, 1, 1 ], [  ] ], 6 ],
#          [ v^3, [ [ 1 ], [ 1, 1, 1 ] ], 7 ],
#          [ v^3, [ [ 1, 1, 1 ], [ 1 ] ], 7 ],
#          [ v^4, [ [ 1, 1, 1, 1 ], [  ] ], 12 ] ] ],
#  [ [ [ 1 ], [ 3 ] ], 1, [ [ v^0, [ [ 1 ], [ 3 ] ], 1 ],
#          [ v, [ [ 2 ], [ 2 ] ], 2 ], [ v, [ [  ], [ 3, 1 ] ], 2 ],
#          [ v^2, [ [ 1, 1 ], [ 2 ] ], 3 ], [ v^2, [ [ 2 ], [ 1, 1 ] ], 3 ],
#          [ v^2, [ [  ], [ 2, 2 ] ], 3 ], [ v^3, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ v^3, [ [  ], [ 2, 1, 1 ] ], 6 ],
#          [ v^4, [ [ 1 ], [ 1, 1, 1 ] ], 7 ] ] ],
#  [ [ [ 3 ], [ 1 ] ], 1, [ [ v^0, [ [ 3 ], [ 1 ] ], 1 ],
#          [ v, [ [ 2 ], [ 2 ] ], 2 ], [ v, [ [ 3, 1 ], [  ] ], 2 ],
#          [ v^2, [ [ 2 ], [ 1, 1 ] ], 3 ], [ v^2, [ [ 1, 1 ], [ 2 ] ], 3 ],
#          [ v^2, [ [ 2, 2 ], [  ] ], 3 ], [ v^3, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ v^3, [ [ 2, 1, 1 ], [  ] ], 6 ],
#          [ v^4, [ [ 1, 1, 1 ], [ 1 ] ], 7 ] ] ],
#  [ [ [ 1 ], [ 2, 1 ] ], 3, [ [ v^0, [ [ 1 ], [ 2, 1 ] ], 3 ] ] ],
#  [ [ [ 2, 1 ], [ 1 ] ], 3, [ [ v^0, [ [ 2, 1 ], [ 1 ] ], 3 ] ] ] ]

# Specialisation at v=1

decomposition:=function(parametre,e,l,n,flotw)local i,resu,v,j,liste;
i:=1; v:=X(Rationals); v.name:="v";
resu:=Basecanonique(parametre,e,l,v,n,flotw);
for i in [1..Length(resu)] do
 liste:=resu[i][3];
   for j in [1..Length(liste)] do
       resu[i][3][j]:=[Value(liste[j][1],1),liste[j][2],liste[j][3]];
       od;
 od;
return resu;
end;

#gap> decomposition([0,1],2,2,4);
#[ [ [ [  ], [ 4 ] ], 0,
#      [ [ 1, [ [  ], [ 4 ] ], 0 ], [ 1, [ [ 3 ], [ 1 ] ], 1 ],
#          [ 1, [ [ 1 ], [ 3 ] ], 1 ], [ 1, [ [ 2 ], [ 2 ] ], 2 ],
#          [ 1, [ [  ], [ 3, 1 ] ], 2 ], [ 1, [ [ 2 ], [ 1, 1 ] ], 3 ],
#          [ 1, [ [ 1, 1 ], [ 2 ] ], 3 ], [ 1, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ 1, [ [  ], [ 2, 1, 1 ] ], 6 ], [ 1, [ [ 1, 1, 1 ], [ 1 ] ], 7 ],
#          [ 1, [ [ 1 ], [ 1, 1, 1 ] ], 7 ],
#          [ 1, [ [  ], [ 1, 1, 1, 1 ] ], 12 ] ] ],
#  [ [ [ 4 ], [  ] ], 0, [ [ 1, [ [ 4 ], [  ] ], 0 ],
#          [ 1, [ [ 1 ], [ 3 ] ], 1 ], [ 1, [ [ 3 ], [ 1 ] ], 1 ],
#          [ 1, [ [ 2 ], [ 2 ] ], 2 ], [ 1, [ [ 3, 1 ], [  ] ], 2 ],
#          [ 1, [ [ 1, 1 ], [ 2 ] ], 3 ], [ 1, [ [ 2 ], [ 1, 1 ] ], 3 ],
#          [ 1, [ [ 1, 1 ], [ 1, 1 ] ], 6 ], [ 1, [ [ 2, 1, 1 ], [  ] ], 6 ],
#          [ 1, [ [ 1 ], [ 1, 1, 1 ] ], 7 ], [ 1, [ [ 1, 1, 1 ], [ 1 ] ], 7 ],
#          [ 1, [ [ 1, 1, 1, 1 ], [  ] ], 12 ] ] ],
#  [ [ [ 1 ], [ 3 ] ], 1, [ [ 1, [ [ 1 ], [ 3 ] ], 1 ],
#          [ 1, [ [ 2 ], [ 2 ] ], 2 ], [ 1, [ [  ], [ 3, 1 ] ], 2 ],
#          [ 1, [ [ 1, 1 ], [ 2 ] ], 3 ], [ 1, [ [ 2 ], [ 1, 1 ] ], 3 ],
#          [ 1, [ [  ], [ 2, 2 ] ], 3 ], [ 1, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ 1, [ [  ], [ 2, 1, 1 ] ], 6 ], [ 1, [ [ 1 ], [ 1, 1, 1 ] ], 7 ] ]
#     ],
#  [ [ [ 3 ], [ 1 ] ], 1, [ [ 1, [ [ 3 ], [ 1 ] ], 1 ], [ 1, [ [ 2 ], [ 2 ] ],
#              2 ], [ 1, [ [ 3, 1 ], [  ] ], 2 ],
#          [ 1, [ [ 2 ], [ 1, 1 ] ], 3 ], [ 1, [ [ 1, 1 ], [ 2 ] ], 3 ],
#          [ 1, [ [ 2, 2 ], [  ] ], 3 ], [ 1, [ [ 1, 1 ], [ 1, 1 ] ], 6 ],
#          [ 1, [ [ 2, 1, 1 ], [  ] ], 6 ], [ 1, [ [ 1, 1, 1 ], [ 1 ] ], 7 ] ]
#     ], [ [ [ 1 ], [ 2, 1 ] ], 3, [ [ 1, [ [ 1 ], [ 2, 1 ] ], 3 ] ] ],
#  [ [ [ 2, 1 ], [ 1 ] ], 3, [ [ 1, [ [ 2, 1 ], [ 1 ] ], 3 ] ] ] ]

verif:=function(liste)local i,j,lo;
  i:=1;
  for i in [1..Length(liste)] do
    lo:=liste[i][3];
    if lo[1][2]<>liste[i][1] then return [liste[i][1],lo[1][2]];
    elif lo[1][1]<>v^0 then return   [liste[i][1],lo[1][2]] ;
    fi;
    j:=2;
    if Length(lo)=1 then i:=i+1; else
      while j<Length(lo)+1 do
	if Value(lo[j][1],0)=0 then j:=j+1;
	else return  [liste[i][1],lo[j][2]] ;fi;
      od;
    fi;
  od;
  return "OK";
end;

Crystal:=function(parametre,e,l,n)
local ensemble,ensemble1,ensemble2,avaleurs,decompo,base,long,
      ensemble4,flotw2,i,j,p,liste,lala,k,la,m,li,flotw,v;
v:=X(Rationals);v.name:="v";ensemble:=PartitionTuples(n,l);
long:=Length(ensemble);
flotw2:=[]; li:=[]; flotw:=[]; ensemble1:=[]; ensemble2:=[]; ensemble4:=[];
avaleurs:=[];
for i in [1..long] do
  Add(ensemble1,[ensemble[i],afonction2(l,ensemble[i],n,parametre,e)]);
od;
SortBy(ensemble1,x->x[2]);
for i in [1..long] do
  Add(ensemble4,ensemble1[i][1]);
  Add(avaleurs,ensemble1[i][2]);
  if Condition1(l,parametre,ensemble1[i][1],e)+Condition2(l,parametre,ensemble1[i][1],e)=2
   then Add(flotw2,ensemble1[i]);Add(ensemble2,[ensemble1[i][1],i]);
                                 Add(flotw,ensemble1[i][1]);
   fi;
od;
decompo:=[];
Base:=Basecanonique(parametre,e,l,v,n,Copy(flotw2));
for j in [1..Length(ensemble2)] do
   m:=ensemble2[j][2]; liste:=[]; lala:=Base[j][3];
  for i in [1..long] do Add(liste,0); od;
  for k in [1..Length(lala)] do
      m:=j; la:=Base[j][3][k][2];
      while ensemble4[m]<>la do m:=m+1; od;
      p:=Base[j][3][k][1];
      if p=0*v^0 then liste[m]:=0;
      elif p=v^0 then liste[m]:=1;
      else liste[m]:=Base[j][3][k][1];
      fi;
  od;
Add(li,liste);
od;
return [ensemble4,avaleurs,flotw,li];
end;

Matricedecompo:=function(parametre,e,l,n)
local ensemble,ensemble1,ensemble2,avaleurs,decompo,base,long,
      ensemble4,flotw2,i,j,p,liste,lala,k,la,m,li,flotw,v;
v:=X(Rationals);v.name:="v"; ensemble:=PartitionTuples(n,l);
long:=Length(ensemble); flotw2:=[]; li:=[]; flotw:=[]; ensemble1:=[];
ensemble2:=[]; ensemble4:=[]; avaleurs:=[];
for i in [1..long] do
    Add(ensemble1,[ensemble[i],afonction2(l,ensemble[i],n,parametre,e)]);
  od;
SortBy(ensemble1,x->x[2]);
for i in [1..long] do
   Add(ensemble4,ensemble1[i][1]);
   Add(avaleurs,ensemble1[i][2]);
   if Condition1(l,parametre,ensemble1[i][1],e)+Condition2(l,parametre,ensemble1[i][1],e)=2
   then Add(flotw2,ensemble1[i]);Add(ensemble2,[ensemble1[i][1],i]);
              Add(flotw,ensemble1[i][1]);
   fi;
od;
decompo:=[];
Base:=decomposition(parametre,e,l,n,Copy(flotw2));
for j in [1..Length(ensemble2)] do
   m:=ensemble2[j][2];
  liste:=[];
  lala:=Base[j][3];
  for i in [1..long] do Add(liste,0); od;
  for k in [1..Length(lala)] do
      m:=j;
      la:=Base[j][3][k][2];
      while ensemble4[m]<>la do m:=m+1; od;
      p:=Base[j][3][k][1];
      if p=0*v^0 then liste[m]:=0;
      elif p=v^0 then liste[m]:=1;
      else liste[m]:=Base[j][3][k][1];
      fi;
  od;
Add(li,liste);
od;
return [ensemble4,avaleurs,flotw,li];
end;

#gap> Matricedecompo([0],4,1,4);
#[ [ [ [ 4 ] ], [ [ 3, 1 ] ], [ [ 2, 2 ] ], [ [ 2, 1, 1 ] ],
#      [ [ 1, 1, 1, 1 ] ] ], [ 0, 1, 2, 3, 6 ],
# [ [ [ 4 ] ], [ [ 3, 1 ] ], [ [ 2, 2 ] ], [ [ 2, 1, 1 ] ] ],
# [ [ 1, 1, 0, 0, 0 ], [ 0, 1, 0, 1, 0 ], [ 0, 0, 1, 0, 0 ],
#     [ 0, 0, 0, 1, 1 ] ] ]

Essai1:=function(parametre,e,l,n)local liste, i,j,k,long,loc;
liste:= Matricedecompo(parametre,e,l,n);
loc:=TransposedMat(liste[4]);
for i in [1..Length(loc)] do
  loc[i]:=Concatenation([ PartitionTupleToString(liste[1][i]),"|"],loc[i]);
  od;
return loc;
end;

Essai2:=function(parametre,e,l,n)local liste, i,j,k,long,loc;
liste:= Crystal(parametre,e,l,n);
loc:=TransposedMat(liste[4]);
for i in [1..Length(loc)] do
  loc[i]:=Concatenation([ PartitionTupleToString(liste[1][i]),"|"],loc[i]);
  od;
return loc;
end;

DisplayFockBasis:=function(n,fbase)local i,j;
  for i in [1..Length(fbase)] do
    for j in [1..2*n-Length(fbase[i][1])] do
      Print(" ");
    od;
    Print(fbase[i][1],"   ");
    if fbase[i][2]<100 then
      Print(" ");
    fi;
    if fbase[i][2]<10 then
      Print(" ");
    fi;
    Print(fbase[i][2],"  ");
    for j in [3..Length(fbase[i])] do
      if fbase[i][j]=0 then
        Print(".");
      else
        Print(fbase[i][j]);
      fi;
    od;
    Print("\n");
 od;
end;

MatrixDecomposition:=function(parametre,e,l,n)
return DisplayFockBasis(n,Essai1(parametre,e,l,n));
end;

CrystalMatrixDecomposition:=function(parametre,e,l,n)
return PrintArray(Essai2(parametre,e,l,n));
end;
