#########################################################################
##
#A  specpie.g                               Meinolf Geck and Gunter Malle
##
#Y  Copyright (C) 1997     Equipe des groupes finis, Universite Paris VII
#Y                         and            IWR der Universitaet Heidelberg
##
##  This file contains functions for computing the polynomials (or, rather,
##  the rational functions) associated with the characters of a finite,
##  real or complex reflection group by the algorithm described in our 
##  paper 'On special pieces in the unipotent variety' (to appear in 
##  Experimental Math.). You need CHEVIE to work with them. 
##
##  Let W be a finite (complex or real) reflection group. Then the function 
##  'GenericLambdaP' returns the pair (Lambda, P), where Lambda and P are 
##  two matrices with rows and columns indexed by the irreducible characters 
##  of W. The algorithm for the computation of these two matrices is 
##  described in the Prop. 2.2 of our paper. The function 'GenericLambdaP' 
##  is written such that it also works for the extension of the basic 
##  algorithm discussed in Section 3 of the paper. Here are some examples:
##
##  (1) The matrices in Example 2.10 of our paper are computed as follows:
##
##      gap> W:=CoxeterGroup("H",3);;         # definition of the group
##      gap> l:=GenericLambdaP(W,W,LowestPowerGenericDegrees(W));;
##      gap> PolSpecialPieces(W,l[1],LowestPowerGenericDegrees(W));
##      # sum of lengths of special pieces: u^30
##      [ [ [ 1 ], [ 5 ], [ 3 ], [ 9 ], [ 4 ], [ 7 ], [ 2 ] ], 
##      [ u^0, u^18 + u^14 + u^10 - u^8 - u^4 - 1, u^20 - u^14 - u^10 + u^4, 
##        u^24 - u^18 - u^14 + u^8, u^26 - u^20 - u^16 + u^10, 
##        u^28 - u^26 - u^22 + u^20 - u^18 + u^16 + u^12 - u^10, 
##        u^30 - u^28 - u^24 + u^22 - u^20 + u^18 + u^14 - u^12 ] ]
##
##      This works similarly for all finite Coxeter groups and certain
##      complex reflection groups. The input data (a-values) for the latter
##      are given at the end of this file.
##
##  (2) The matrices in Example 3.1 of our paper are computed as follows:
##      
##      gap> W:=CoxeterGroup("A",1);
##      gap> W1:=CoxeterGroup("H",4);
##      gap> GenericLambdaP(W,W1,[18,3]);
##
##      The input data for various relative pairs as above are given at
##      the end of this file.
##
#H  written MG, GM, Dec 1997
##
#RequirePackage("chevie");
u:=X(Cyclotomics);u.name:="u";

###########################################################################
##
#F MatFakeDegrees( <W>, <chars>, <fak> )  . . . . . . matrix of fake degrees
##
## 'MatFakeDegree' computes the matrix Omega defined as in Lemma 2.1 of our 
## paper. Here, <W> is a finite reflection group, <chars> is a list which
## determines an ordering of the characters of <W> such that the a-values
## are in decreasing order (this is given by the first component of 
## 'PrepData') and <fak> is the list of fake degrees.
##
MatFakeDegree:=function(W,chars,fak)local R,ct,i,j,val,fakes,det,b;
  ct:=CharTable(W).irreducibles{chars};
  fakes:=fak{chars};
  R:=arr->Sum([1..Length(arr)],i->arr[i]*fakes[i]);
  Print("#I Computing Omega ");
  det:=Position(chars,PositionDet(W));
  det:=Position(ct,List(ct[det],i->GaloisCyc(i,-1)));
  Print(" (compl. conj. det character is ",det,"): ");
  val:=[];
  for i in [1..Length(ct)] do
    Print(".\c");
    val[i]:=[];
    for j in [1..i] do
      val[i][j]:=R(MatScalarProducts(CharTable(W),ct,
               Tensored([ct[det]],Tensored([ct[i]],[(ct[j])])))[1]);
      if j<i then val[j][i]:=val[i][j]; fi;
    od;
  od;
  Print("\n");
  return val;
end;

###########################################################################
##
#F LambdaP( <W>, <W1>, <avals>, <omega>, <q> ) . . . . . . . . the matrices 
## . . . . . . . . . . . . . . . . . .  Lambda and P for a given value of q
##
## 'LambdaP' computes the matrices Lambda and P (as in Prop. 2.2 of our
## article) for a given value of <q>. Here, <W> is a finite reflection
## group, <W1> is another reflection group such that (W,W1) is a relative
## pair as in the variation of the algorithm discussed in Section 3, 
## <avals> is the list of a-values of the characters of <W>, <omega> is
## the result of 'MatFakeDegree' applied to <W>, and <q> is any non-zero
## integer. If we are just interested in the non-relative version of the
## algorithm, we use <W1>=<W>.
##
LambdaP:=function(W,W1,avals,omega,q)
  local J,P,Pi,Lambda,data,x,y,i,j,k;
  Pi:=q^(W.semisimpleRank-W1.semisimpleRank)
     *(q-1)^(W1.semisimpleRank-W.semisimpleRank)
     *Product(ReflectionDegrees(W1),i->(q^i-1)/(q-1))
     /Product(ReflectionDegrees(W),i->(q^i-1)/(q-1))
     *q^(Sum(ReflectionDegrees(W1))-Length(W1.generators))
     *List(omega,x->List(x,y->Value(y,q)));
  data:=Reversed(Collected(avals));
  J:=[]; x:=1; y:=1;
  for k in data do Add(J,[x..x+k[2]-1]);x:=x+k[2]; od;
  P:=List([1..x-1],i->0*[1..x-1]); 
  Lambda:=[]; 
  for i in [1..Length(data)] do 
    x:=Pi{J[i]}{J[i]};
    for k in [1..i-1] do
      x:=x-TransposedMat(P{J[k]}{J[i]})*Lambda[k]*(P{J[k]}{J[i]});
    od;
    Lambda[i]:=q^(-2*data[i][1])*x;
    y:=q^(-data[i][1])*Lambda[i]^(-1);
    P{J[i]}{J[i]}:=q^(data[i][1])*x^0; 
    for j in [i+1..Length(data)] do
      x:=Pi{J[i]}{J[j]};
      for k in [1..i-1] do
        x:=x-TransposedMat(P{J[k]}{J[i]})*Lambda[k]*(P{J[k]}{J[j]});
      od;
      P{J[i]}{J[j]}:=(y)*x;
    od;
  od;
  return [Lambda,TransposedMat(P)];
end;

###########################################################################
##
#F PrepData( <W>, <avals>, <fakes> )  . . . . . . . prepare some basic data
##
## 'PrepData' is a utility function which prepares input data for some of 
## the programs below. It takes as input a finite reflection group <W>, the 
## list of a-values of the irreducible characters of <W> and the corresponding 
## fake degrees. The result contains, for example, the information about an
## ordering of the irreducible characters compatible with the a-values
## and the b-values. 
##
PrepData:=function(W,avals,fakes)local ct,b,a1,i,j,l,x,data;
  data:=[];
  b:=List(fakes,i->i.valuation);
  a1:=Reversed(Set(avals));
  ct:=[1..Length(b)];
  data[1]:=[]; data[2]:=[];
  for i in a1 do
    l:=Filtered([1..Length(avals)],j->avals[j]=i);
    SortBy(l,x->b[x]);
    Append(data[1],ct{l});
    Add(data[2],[i,Length(l),Filtered(l,j->b[j]=b[l[1]])]);
  od;
  return data;
end;

###########################################################################
##
#F GenericLambdaP( <W>, <W1>, <avals> ) . . . . . . . . . . . . the matrices 
## . . . . . . . . . . . . . . . . . . . . . . . Lambda and P for generic q
##
## 'GenericLambdaP' calls 'LambdaP' for a large number of values of <q> and
## then interpolates to find polynomial solutions for the coefficients of
## Lambda and P.  Note that, a priori, it is not known that such polynomial
## solutions exist, but as soon as some polynomials have been found by 
## interpolation, their correctness is verified by simply checking the 
## relations in Prop. 2.2, which amounts to a matrix product calculation. 
## (This verification is done in the last step before returning the result.)
## If we are only interested in the non-relative version of the algorithm, 
## we use <W1>=<W>.
## 
GenericLambdaP:=function(W,W1,avals)
  local q,g,c,lambda,pol,x,i,j,l,max,omega,fakes,P,Pi,lpol,Ppol,Lmat,ver;
  q:=X(Cyclotomics);
  fakes:=FakeDegrees(W,q);
  max:=2*Sum(ReflectionDegrees(W1))+2;
  omega:=MatFakeDegree(W,PrepData(W,avals,fakes)[1],fakes);
  lambda:=[];
  P:=[];
  Print("#I Solving P\*Lambda\*P\^t=Omega (",max," times):  \c");
  for i in [1..max] do
    if IsInt(i/10) then Print(i," \c");fi;
    g:=LambdaP(W,W1,avals,omega,i+1);
    lambda[i]:=g[1];
    P[i]:=g[2];
  od;
  Print("\n#I Now interpolating \n");
  Print("#I Lambda: ");
  lpol:=[];
  for l in [1..Length(lambda[1])] do
    lpol[l]:=[];
    for i in [1..Length(lambda[1][l])] do
    Print(".\c"); 
      lpol[l][i]:=[];
      for j in [1..Length(lambda[1][l])] do
        c:=List([1..max],x->lambda[x][l][i][j]);
        lpol[l][i][j]:=InterpolatedPolynomial(Cyclotomics,[2..max+1],c);
      od;
    od;
  od;
  Print("\n#I      P: ");
  Ppol:=[];
  for l in [1..Length(P[1])] do
    Print(".\c"); 
    Ppol[l]:=[];
    for i in [1..Length(P[1][i])] do
      c:=List([1..max],x->P[x][l][i]);
      Ppol[l][i]:=InterpolatedPolynomial(Cyclotomics,[2..max+1],c);
    od;
  od;
  Print("\n#I Verifying correctness of solutions:  \c");
  l:=[];x:=1;
  for i in [1..Length(lpol)] do 
    l[i]:=[x..x+Length(lpol[i])-1]; x:=x+Length(lpol[i]);
  od;
  Lmat:=IdentityMat(Length(Ppol))*(0*u);
  for i in [1..Length(lpol)] do Lmat{l[i]}{l[i]}:=lpol[i]; od;
  ver:=Ppol*Lmat*TransposedMat(Ppol)=q^(W.semisimpleRank-W1.semisimpleRank)
     *(q-1)^(W1.semisimpleRank-W.semisimpleRank)
     *Product(ReflectionDegrees(W1),i->(q^i-1)/(q-1))
     /Product(ReflectionDegrees(W),i->(q^i-1)/(q-1))
     *q^(Sum(ReflectionDegrees(W1))-Length(W1.generators))*omega;
  if ver=true then 
    Print(ver,"\n"); return [lpol,Ppol];
  else 
    return false;
  fi;
end;

###########################################################################
##
#F PolSpecialPieces( <W>, <Lambda>, <avals> ) . . . . . . . . . . . . . . . 
## . . . . . . . . . . . . . . .  return polynomials for special characters
##
## 'PolSpecialPieces' takes a finite reflection group <W>, the result 
## <Lambda> of 'GenericLambdaP( <W>, <W>, <avals> )[1]' and the list 
## <avals> of a-values, and returns the polynomials on the diagonal of 
## Lambda which correspond to special characters of <W>.
##
PolSpecialPieces:=function(W,genpol,avals)
  local prep,p,i,j,l,dec,families,special,fakes;
  fakes:=FakeDegrees(W,X(Cyclotomics));
  prep:=PrepData(W,avals,fakes);
  dec:=Copy(genpol);
  for l in [1..Length(genpol)] do
    for i in [1..Length(genpol[l])] do
      for j in [1..Length(genpol[l])] do
        if genpol[l][i][j]=0*genpol[l][i][j] then dec[l][i][j]:=0;
        else dec[l][i][j]:=1;
        fi;
      od;
    od;
  od;
  p:=List(prep[2],i->i[3]);
  dec:=List([1..Length(genpol)],i->DecomposedMat(dec[i]));
  #Print("# blocks of matrices: ",p,"\n# ",dec,"\n");
  special:=Flat(List([1..Length(genpol)],
                        i->List([1..Length(p[i])],j->genpol[i][j][j])));
  Print("# sum of lengths of special pieces: ",Sum(special),"\n");
  return [p,special];
end;

#####################################
# Data for complex reflection groups

testspecpie:=function(i)local l,g,g1,ag27,ag29,ag33,ag34,g1ag29;
if i="24" then
  g:=CoxeterGroup("A",1);g1:=ComplexReflectionGroup(24); 
  GenericLambdaP(g,g1,[1,8]);

  g:=ComplexReflectionGroup(24);g1:=ComplexReflectionGroup(24); 
  l:=GenericLambdaP(g,g1,[0,21,8,1,8,1,1,8,6,3,4,4]);
  PolSpecialPieces(g,l[1],[0,21,8,1,8,1,1,8,6,3,4,4]); 
elif i="27" then
  g:=ComplexReflectionGroup(6,1,1);g1:=ComplexReflectionGroup(27); 
# GenericLambdaP(g,g1,[1,16,1,16,12,3]);
# GenericLambdaP(g,g1,[1,16,1,16,6,6]);
  g:=ComplexReflectionGroup(27);
  ag27:=[0,45,1,16,1,16,1,16,1,16,3,3,12,12,16,1,16,1,6,6,6,6,4,9,9,4,9,4,
    12,3,5,8,5,8];
  l:=GenericLambdaP(g,g1,ag27); PolSpecialPieces(g,l[1],ag27);
elif i="29" then
  g:=ComplexReflectionGroup(29);g1:=g;
  ag29:=[0,40,1,21,1,21,1,21,4,12,6,6,6,6,6,18,2,12,4,4,4,12,12,13,3,13,3,
    5,9,9,5,9,5,6,6,6,6];
  l:=GenericLambdaP(g,g1ag29); PolSpecialPieces(g,l[1],ag29);
elif i="33" then
  g:=ComplexReflectionGroup(33);g1:=g;
  ag33:=[0,45,28,1,28,1,18,3,4,13,4,13,12,9,2,23,4,13, 3,18,18,3,4,13,
    4,13,13,4,13,4,10,7,10,7,10,7,8,8,6,11];
  l:=GenericLambdaP(g,g1,ag33); PolSpecialPieces(g,l[1],ag33);
elif i="34" then
  g:=ComplexReflectionGroup(34);g1:=g;
  ag34:=[0,126,1,85,1,85,4,46,4,46,15,15,57,3,2,68,2,68,57,3,28,10,3,57,5,
  41,9,45,5,5,41,41,7,31,31,7,4,46,36,4,46,36,6,4,46,6,13,19,13,19,15,15,4,
  46,4,46,5,41,5,41,18,27,9,19,13,30,12,36,36,6,6,23,11,23,11,28,28,10,10,
  36,6,19,13,10,28,10,28,7,31,7,31,19,13,19,13,8,29,8,29,8,29,8,29,27,9,18,
  7,31,31,7,13,19,13,19,7,31,13,19,15,15,15,15,9,27,18,15,15,11,23,14,20,14,
  20,23,11,23,11,15,15,24,10,10,24,10,24,13,19,13,19,19,13,19,13,23,11,23,
  11,21,12,12,21,13,19,13,19,15,15,15,15];
  l:=GenericLambdaP(g,g1,ag34); PolSpecialPieces(g,l[1],ag34);
##################################################################
# Data for relative pairs (W,W1) where W is a finite Coxeter group.
elif i="B2 in F4" then
  g:=CoxeterGroup("B",2);g1:=CoxeterGroup("F",4);
  l:=GenericLambdaP(g,g1,[4,4,13,1,5]);
elif i="A2 in E6" then
  g:=CoxeterGroup("A",2);g1:=CoxeterGroup("E",6);
  GenericLambdaP(g,g1,[15,7,3]);
elif i="A1 in E7" then
  g:=CoxeterGroup("A",1);g1:=CoxeterGroup("E",7);
  GenericLambdaP(g,g1,[16,7]);
elif i="C3 in E7" then
  g:=CoxeterGroup("C",3);g1:=CoxeterGroup("E",7);
  GenericLambdaP(g,g1,[25,15,16,30,10,8,7,13,3,4]);
elif i="A1 in E8" then
  g:=CoxeterGroup("A",1);g1:=CoxeterGroup("E",8);
  GenericLambdaP(g,g1,[26,11]);
elif i="G2 in E8" then
  g:=CoxeterGroup("G",2);g1:=CoxeterGroup("E",8);
  GenericLambdaP(g,g1,[7,37,8,32,16,16]);
elif i="F4 in E8" then
  g:=CoxeterGroup("F",4);g1:=CoxeterGroup("E",8);
  GenericLambdaP(g,g1,
  [3,42,6,63,12,24,4,52,16,10,25,13,30,16,16,16,7,28,10,37,15,21,8,32,16]);
elif i="A1 in H4" then
  g:=CoxeterGroup("A",1);g1:=CoxeterGroup("H",4);
  GenericLambdaP(g,g1,[18,3]);
elif i="I2(10) in H4" then
  g:=CoxeterGroup("I",2,10);g1:=CoxeterGroup("H",4);
  GenericLambdaP(g,g1,[0,22,2,31,6,6,6,6]);
fi;
end;
# Finally, since it takes quite some time and space, here is the result of 
# W:=CoxeterGroup("E",8);GenericLambdaP(W,W,LowestPowerGenericDegrees(W));
#lambdapE8:=
#[ [ [ [ u^0 ] ], 
#      [ [ u^58 + u^52 + u^48 + u^46 + u^42 + u^40 + u^36 + u^30 - u^28 - u^
#                22 - u^18 - u^16 - u^12 - u^10 - u^6 - 1 ] ], 
#      [ [ u^92 + u^88 + u^86 + u^84 + 2*u^82 + 2*u^80 + u^78 + 3*u^76 + u^
#                74 + 2*u^72 + 2*u^70 + u^66 - 2*u^62 - 3*u^58 - 3*u^56 - 2*u^
#                54 - 5*u^52 - 3*u^50 - 3*u^48 - 5*u^46 - 2*u^44 - 3*u^42 - 
#                3*u^40 - 2*u^36 + u^32 + 2*u^28 + 2*u^26 + u^24 + 3*u^22 + u^
#                20 + 2*u^18 + 2*u^16 + u^14 + u^12 + u^10 + u^6 ] ], 
#      [ [ u^114 + u^112 + u^110 + 2*u^108 + 2*u^106 + 3*u^104 + 3*u^102 + 3*u^
#                100 + 3*u^98 + 3*u^96 + 2*u^94 + u^92 - u^88 - 3*u^86 - 4*u^
#                84 - 5*u^82 - 7*u^80 - 6*u^78 - 8*u^76 - 8*u^74 - 6*u^72 - 
#                7*u^70 - 5*u^68 - 3*u^66 - 3*u^64 + u^62 + 2*u^60 + 2*u^58 + 
#                6*u^56 + 5*u^54 + 6*u^52 + 8*u^50 + 5*u^48 + 7*u^46 + 6*u^
#                44 + 3*u^42 + 4*u^40 + 2*u^38 + u^34 - 2*u^32 - 2*u^30 - u^
#                28 - 3*u^26 - 2*u^24 - 2*u^22 - 2*u^20 - u^18 - u^16 - u^14, 
#              u^113 + u^111 + u^109 + 2*u^107 + 2*u^105 + 2*u^103 + 3*u^101 + 
#                2*u^99 + 2*u^97 + 2*u^95 - u^89 - 3*u^87 - 3*u^85 - 5*u^83 - 
#                6*u^81 - 5*u^79 - 7*u^77 - 6*u^75 - 5*u^73 - 6*u^71 - 3*u^
#                69 - 2*u^67 - 2*u^65 + 2*u^63 + 2*u^61 + 3*u^59 + 6*u^57 + 
#                5*u^55 + 6*u^53 + 7*u^51 + 5*u^49 + 6*u^47 + 5*u^45 + 3*u^
#                43 + 3*u^41 + u^39 - 2*u^33 - 2*u^31 - 2*u^29 - 3*u^27 - 2*u^
#                25 - 2*u^23 - 2*u^21 - u^19 - u^17 - u^15, 
#              u^109 + u^105 + u^103 + 2*u^99 + u^97 + 2*u^93 + 2*u^87 - 2*u^
#                85 - 3*u^79 - 2*u^75 - 3*u^73 - 3*u^69 - 2*u^67 - 3*u^63 - 
#                2*u^57 + 2*u^55 + 2*u^49 + u^45 + 2*u^43 + u^39 + u^37 + u^33 
#             ], 
#          [ u^113 + u^111 + u^109 + 2*u^107 + 2*u^105 + 2*u^103 + 3*u^101 + 
#                2*u^99 + 2*u^97 + 2*u^95 - u^89 - 3*u^87 - 3*u^85 - 5*u^83 - 
#                6*u^81 - 5*u^79 - 7*u^77 - 6*u^75 - 5*u^73 - 6*u^71 - 3*u^
#                69 - 2*u^67 - 2*u^65 + 2*u^63 + 2*u^61 + 3*u^59 + 6*u^57 + 
#                5*u^55 + 6*u^53 + 7*u^51 + 5*u^49 + 6*u^47 + 5*u^45 + 3*u^
#                43 + 3*u^41 + u^39 - 2*u^33 - 2*u^31 - 2*u^29 - 3*u^27 - 2*u^
#                25 - 2*u^23 - 2*u^21 - u^19 - u^17 - u^15, 
#              u^114 + u^112 + u^110 + 2*u^108 + 2*u^106 + 2*u^104 + 3*u^102 + 
#                2*u^100 + 2*u^98 + 2*u^96 - u^90 - 3*u^88 - 3*u^86 - 5*u^84 - 
#                6*u^82 - 5*u^80 - 7*u^78 - 6*u^76 - 5*u^74 - 6*u^72 - 3*u^
#                70 - 2*u^68 - 2*u^66 + 2*u^64 + 2*u^62 + 3*u^60 + 6*u^58 + 
#                5*u^56 + 6*u^54 + 7*u^52 + 5*u^50 + 6*u^48 + 5*u^46 + 3*u^
#                44 + 3*u^42 + u^40 - 2*u^34 - 2*u^32 - 2*u^30 - 3*u^28 - 2*u^
#                26 - 2*u^24 - 2*u^22 - u^20 - u^18 - u^16, 0*u^0 ], 
#          [ u^109 + u^105 + u^103 + 2*u^99 + u^97 + 2*u^93 + 2*u^87 - 2*u^
#                85 - 3*u^79 - 2*u^75 - 3*u^73 - 3*u^69 - 2*u^67 - 3*u^63 - 
#                2*u^57 + 2*u^55 + 2*u^49 + u^45 + 2*u^43 + u^39 + u^37 + u^33,
#              0*u^0, u^114 + u^108 + u^104 + u^102 + u^98 + u^96 + u^92 - u^
#                84 - u^80 - u^78 - u^76 - 2*u^74 - u^72 - u^70 - 2*u^68 - u^
#                66 - u^64 - u^62 - u^58 + u^50 + u^46 + u^44 + u^40 + u^
#                38 + u^34 + u^28 ] ], 
#      [ [ u^136 + u^132 + 2*u^130 + 2*u^128 + 3*u^126 + 4*u^124 + 3*u^122 + 
#                5*u^120 + 4*u^118 + 3*u^116 + 3*u^114 + u^112 - u^110 - 2*u^
#                108 - 5*u^106 - 7*u^104 - 8*u^102 - 11*u^100 - 11*u^98 - 11*u^
#                96 - 12*u^94 - 10*u^92 - 8*u^90 - 8*u^88 - 2*u^86 - u^84 + u^
#                82 + 7*u^80 + 7*u^78 + 9*u^76 + 14*u^74 + 10*u^72 + 13*u^70 + 
#                13*u^68 + 8*u^66 + 10*u^64 + 7*u^62 + 2*u^60 + 4*u^58 - 2*u^
#                56 - 3*u^54 - 2*u^52 - 7*u^50 - 5*u^48 - 5*u^46 - 7*u^44 - 
#                3*u^42 - 4*u^40 - 4*u^38 - 2*u^34 + u^30 + u^26 + u^24 + u^20,
#              u^133 + u^131 + u^129 + 3*u^127 + 2*u^125 + 3*u^123 + 5*u^121 + 
#                2*u^119 + 4*u^117 + 4*u^115 + 3*u^111 - u^109 - 4*u^107 - u^
#                105 - 8*u^103 - 7*u^101 - 6*u^99 - 13*u^97 - 7*u^95 - 9*u^
#                93 - 12*u^91 - 3*u^89 - 8*u^87 - 5*u^85 + 3*u^83 - 3*u^81 + 
#                5*u^79 + 8*u^77 + 3*u^75 + 12*u^73 + 9*u^71 + 7*u^69 + 13*u^
#                67 + 6*u^65 + 7*u^63 + 8*u^61 + u^59 + 4*u^57 + u^55 - 3*u^
#                53 - 4*u^49 - 4*u^47 - 2*u^45 - 5*u^43 - 3*u^41 - 2*u^39 - 
#                3*u^37 - u^35 - u^33 - u^31, 
#              u^132 + u^128 + u^126 + u^124 + u^122 + u^120 + u^116 - u^
#                114 - u^112 - u^110 - 3*u^108 - 2*u^106 - 3*u^104 - 4*u^102 - 
#                2*u^100 - 3*u^98 - 3*u^96 - 2*u^92 + u^90 + 2*u^88 + u^86 + 
#                4*u^84 + 4*u^82 + 3*u^80 + 6*u^78 + 3*u^76 + 4*u^74 + 4*u^
#                72 + u^70 + 2*u^68 + u^66 - 2*u^64 - 3*u^60 - 3*u^58 - 2*u^
#                56 - 4*u^54 - 3*u^52 - 2*u^50 - 3*u^48 - u^46 - u^44 - u^
#                42 + u^40 + u^36 + u^34 + u^32 + u^30 + u^28 + u^24 ], 
#          [ u^133 + u^131 + u^129 + 3*u^127 + 2*u^125 + 3*u^123 + 5*u^121 + 
#                2*u^119 + 4*u^117 + 4*u^115 + 3*u^111 - u^109 - 4*u^107 - u^
#                105 - 8*u^103 - 7*u^101 - 6*u^99 - 13*u^97 - 7*u^95 - 9*u^
#                93 - 12*u^91 - 3*u^89 - 8*u^87 - 5*u^85 + 3*u^83 - 3*u^81 + 
#                5*u^79 + 8*u^77 + 3*u^75 + 12*u^73 + 9*u^71 + 7*u^69 + 13*u^
#                67 + 6*u^65 + 7*u^63 + 8*u^61 + u^59 + 4*u^57 + u^55 - 3*u^
#                53 - 4*u^49 - 4*u^47 - 2*u^45 - 5*u^43 - 3*u^41 - 2*u^39 - 
#                3*u^37 - u^35 - u^33 - u^31, 
#              u^136 + u^132 + 2*u^130 + u^128 + 3*u^126 + 3*u^124 + 2*u^122 + 
#                4*u^120 + 3*u^118 + 2*u^116 + 3*u^114 - u^108 - 4*u^106 - 4*u^
#                104 - 6*u^102 - 8*u^100 - 7*u^98 - 9*u^96 - 9*u^94 - 7*u^92 - 
#                8*u^90 - 6*u^88 - 3*u^86 - 3*u^84 + 3*u^80 + 3*u^78 + 6*u^
#                76 + 8*u^74 + 7*u^72 + 9*u^70 + 9*u^68 + 7*u^66 + 8*u^64 + 
#                6*u^62 + 4*u^60 + 4*u^58 + u^56 - 3*u^50 - 2*u^48 - 3*u^46 - 
#                4*u^44 - 2*u^42 - 3*u^40 - 3*u^38 - u^36 - 2*u^34 - u^32 - u^
#                28, 0*u^0 ], 
#          [ u^132 + u^128 + u^126 + u^124 + u^122 + u^120 + u^116 - u^114 - u^
#                112 - u^110 - 3*u^108 - 2*u^106 - 3*u^104 - 4*u^102 - 2*u^
#                100 - 3*u^98 - 3*u^96 - 2*u^92 + u^90 + 2*u^88 + u^86 + 4*u^
#                84 + 4*u^82 + 3*u^80 + 6*u^78 + 3*u^76 + 4*u^74 + 4*u^72 + u^
#                70 + 2*u^68 + u^66 - 2*u^64 - 3*u^60 - 3*u^58 - 2*u^56 - 4*u^
#                54 - 3*u^52 - 2*u^50 - 3*u^48 - u^46 - u^44 - u^42 + u^40 + u^
#                36 + u^34 + u^32 + u^30 + u^28 + u^24, 0*u^0, 
#              u^136 + u^132 + u^130 + u^128 + u^126 + u^124 + u^120 - u^
#                118 - u^116 - u^114 - 3*u^112 - 2*u^110 - 3*u^108 - 4*u^106 - 
#                2*u^104 - 3*u^102 - 3*u^100 - 2*u^96 + u^94 + 2*u^92 + u^90 + 
#                4*u^88 + 4*u^86 + 3*u^84 + 6*u^82 + 3*u^80 + 4*u^78 + 4*u^
#                76 + u^74 + 2*u^72 + u^70 - 2*u^68 - 3*u^64 - 3*u^62 - 2*u^
#                60 - 4*u^58 - 3*u^56 - 2*u^54 - 3*u^52 - u^50 - u^48 - u^
#                46 + u^44 + u^40 + u^38 + u^36 + u^34 + u^32 + u^28 ] ], 
#      [ [ u^146 + u^144 + 2*u^142 + 3*u^140 + 3*u^138 + 4*u^136 + 4*u^134 + 
#                3*u^132 + 3*u^130 + u^128 - u^126 - 2*u^124 - 6*u^122 - 7*u^
#                120 - 9*u^118 - 12*u^116 - 11*u^114 - 12*u^112 - 12*u^110 - 
#                8*u^108 - 8*u^106 - 4*u^104 + u^102 + 2*u^100 + 8*u^98 + 11*u^
#                96 + 12*u^94 + 17*u^92 + 16*u^90 + 16*u^88 + 17*u^86 + 12*u^
#                84 + 11*u^82 + 8*u^80 + 2*u^78 + u^76 - 4*u^74 - 8*u^72 - 8*u^
#                70 - 12*u^68 - 12*u^66 - 11*u^64 - 12*u^62 - 9*u^60 - 7*u^
#                58 - 6*u^56 - 2*u^54 - u^52 + u^50 + 3*u^48 + 3*u^46 + 4*u^
#                44 + 4*u^42 + 3*u^40 + 3*u^38 + 2*u^36 + u^34 + u^32 ] ], 
#      [ [ u^148 + u^144 + u^142 + u^140 + 2*u^138 + u^136 + u^134 + 2*u^
#                132 + u^128 - 2*u^124 - 3*u^120 - 3*u^118 - 2*u^116 - 5*u^
#                114 - 3*u^112 - 3*u^110 - 5*u^108 - u^106 - 3*u^104 - 2*u^
#                102 + u^100 - u^98 + 2*u^96 + 3*u^94 + u^92 + 5*u^90 + 3*u^
#                88 + 3*u^86 + 5*u^84 + 2*u^82 + 3*u^80 + 3*u^78 + 2*u^74 - u^
#                70 - 2*u^66 - u^64 - u^62 - 2*u^60 - u^58 - u^56 - u^54 - u^
#                50 ] ], 
#      [ [ u^156 + u^154 + 2*u^152 + 2*u^150 + 3*u^148 + 2*u^146 + 3*u^144 + u^
#                142 - u^138 - 3*u^136 - 5*u^134 - 5*u^132 - 8*u^130 - 8*u^
#                128 - 7*u^126 - 9*u^124 - 5*u^122 - 4*u^120 - 3*u^118 + 2*u^
#                116 + 3*u^114 + 5*u^112 + 11*u^110 + 8*u^108 + 12*u^106 + 
#                13*u^104 + 9*u^102 + 12*u^100 + 9*u^98 + 3*u^96 + 6*u^94 - 
#                2*u^92 - 4*u^90 - 3*u^88 - 10*u^86 - 9*u^84 - 8*u^82 - 13*u^
#                80 - 7*u^78 - 8*u^76 - 8*u^74 - 2*u^72 - 3*u^70 + 4*u^66 + 
#                2*u^64 + 5*u^62 + 6*u^60 + 3*u^58 + 6*u^56 + 3*u^54 + 2*u^
#                52 + 2*u^50 - u^46 - 2*u^42 - u^40 - u^38 - u^36, 
#              u^155 + u^153 + u^151 + 2*u^149 + u^147 + u^145 + u^143 - u^
#                141 - u^139 - 2*u^137 - 4*u^135 - 3*u^133 - 5*u^131 - 5*u^
#                129 - 3*u^127 - 5*u^125 - 2*u^123 - u^119 + 4*u^117 + 4*u^
#                115 + 4*u^113 + 9*u^111 + 6*u^109 + 7*u^107 + 9*u^105 + 4*u^
#                103 + 6*u^101 + 4*u^99 - u^97 + u^95 - 4*u^93 - 6*u^91 - 4*u^
#                89 - 9*u^87 - 7*u^85 - 6*u^83 - 9*u^81 - 4*u^79 - 4*u^77 - 
#                4*u^75 + u^73 + 2*u^69 + 5*u^67 + 3*u^65 + 5*u^63 + 5*u^61 + 
#                3*u^59 + 4*u^57 + 2*u^55 + u^53 + u^51 - u^49 - u^47 - u^45 - 
#                2*u^43 - u^41 - u^39 - u^37, 
#              u^154 + 2*u^150 + u^148 + u^146 + 2*u^144 + u^142 + 2*u^138 - 
#                2*u^136 - u^134 - u^132 - 5*u^130 - 2*u^128 - 4*u^126 - 6*u^
#                124 - 2*u^122 - 5*u^120 - 4*u^118 + u^116 - 4*u^114 + 2*u^
#                112 + 3*u^110 + 7*u^106 + 5*u^104 + 3*u^102 + 10*u^100 + 3*u^
#                98 + 5*u^96 + 7*u^94 + 3*u^90 + 2*u^88 - 4*u^86 + u^84 - 4*u^
#                82 - 5*u^80 - 2*u^78 - 6*u^76 - 4*u^74 - 2*u^72 - 5*u^70 - u^
#                68 - u^66 - 2*u^64 + 2*u^62 + u^58 + 2*u^56 + u^54 + u^52 + 
#                2*u^50 + u^46 ], 
#          [ u^155 + u^153 + u^151 + 2*u^149 + u^147 + u^145 + u^143 - u^
#                141 - u^139 - 2*u^137 - 4*u^135 - 3*u^133 - 5*u^131 - 5*u^
#                129 - 3*u^127 - 5*u^125 - 2*u^123 - u^119 + 4*u^117 + 4*u^
#                115 + 4*u^113 + 9*u^111 + 6*u^109 + 7*u^107 + 9*u^105 + 4*u^
#                103 + 6*u^101 + 4*u^99 - u^97 + u^95 - 4*u^93 - 6*u^91 - 4*u^
#                89 - 9*u^87 - 7*u^85 - 6*u^83 - 9*u^81 - 4*u^79 - 4*u^77 - 
#                4*u^75 + u^73 + 2*u^69 + 5*u^67 + 3*u^65 + 5*u^63 + 5*u^61 + 
#                3*u^59 + 4*u^57 + 2*u^55 + u^53 + u^51 - u^49 - u^47 - u^45 - 
#                2*u^43 - u^41 - u^39 - u^37, 
#              u^156 + u^154 + u^152 + 2*u^150 + u^148 + u^146 + u^144 - u^
#                142 - u^140 - 2*u^138 - 4*u^136 - 3*u^134 - 5*u^132 - 5*u^
#                130 - 3*u^128 - 5*u^126 - 2*u^124 - u^120 + 4*u^118 + 4*u^
#                116 + 4*u^114 + 9*u^112 + 6*u^110 + 7*u^108 + 9*u^106 + 4*u^
#                104 + 6*u^102 + 4*u^100 - u^98 + u^96 - 4*u^94 - 6*u^92 - 4*u^
#                90 - 9*u^88 - 7*u^86 - 6*u^84 - 9*u^82 - 4*u^80 - 4*u^78 - 
#                4*u^76 + u^74 + 2*u^70 + 5*u^68 + 3*u^66 + 5*u^64 + 5*u^62 + 
#                3*u^60 + 4*u^58 + 2*u^56 + u^54 + u^52 - u^50 - u^48 - u^46 - 
#                2*u^44 - u^42 - u^40 - u^38, 0*u^0 ], 
#          [ u^154 + 2*u^150 + u^148 + u^146 + 2*u^144 + u^142 + 2*u^138 - 2*u^
#                136 - u^134 - u^132 - 5*u^130 - 2*u^128 - 4*u^126 - 6*u^124 - 
#                2*u^122 - 5*u^120 - 4*u^118 + u^116 - 4*u^114 + 2*u^112 + 3*u^
#                110 + 7*u^106 + 5*u^104 + 3*u^102 + 10*u^100 + 3*u^98 + 5*u^
#                96 + 7*u^94 + 3*u^90 + 2*u^88 - 4*u^86 + u^84 - 4*u^82 - 5*u^
#                80 - 2*u^78 - 6*u^76 - 4*u^74 - 2*u^72 - 5*u^70 - u^68 - u^
#                66 - 2*u^64 + 2*u^62 + u^58 + 2*u^56 + u^54 + u^52 + 2*u^
#                50 + u^46, 0*u^0, u^156 + u^152 + u^150 + u^148 + u^146 + 2*u^
#                144 + u^140 - u^136 - u^134 - 2*u^132 - 3*u^130 - 3*u^128 - 
#                4*u^126 - 4*u^124 - 3*u^122 - 4*u^120 - 2*u^118 - 2*u^116 - u^
#                114 + u^112 + 2*u^110 + 2*u^108 + 5*u^106 + 4*u^104 + 5*u^
#                102 + 6*u^100 + 5*u^98 + 4*u^96 + 5*u^94 + 2*u^92 + 2*u^
#                90 + u^88 - u^86 - 2*u^84 - 2*u^82 - 4*u^80 - 3*u^78 - 4*u^
#                76 - 4*u^74 - 3*u^72 - 3*u^70 - 2*u^68 - u^66 - u^64 + u^60 + 
#                2*u^56 + u^54 + u^52 + u^50 + u^48 + u^44 ] ], 
#      [ [ u^166 + u^164 + 2*u^162 + 4*u^160 + 4*u^158 + 5*u^156 + 6*u^154 + 
#                4*u^152 + 4*u^150 + 2*u^148 - 2*u^146 - 3*u^144 - 7*u^142 - 
#                11*u^140 - 12*u^138 - 16*u^136 - 16*u^134 - 15*u^132 - 16*u^
#                130 - 11*u^128 - 8*u^126 - 5*u^124 + 3*u^122 + 6*u^120 + 11*u^
#                118 + 18*u^116 + 18*u^114 + 22*u^112 + 24*u^110 + 21*u^108 + 
#                21*u^106 + 17*u^104 + 11*u^102 + 8*u^100 + u^98 - 4*u^96 - 
#                8*u^94 - 13*u^92 - 15*u^90 - 17*u^88 - 18*u^86 - 16*u^84 - 
#                15*u^82 - 12*u^80 - 8*u^78 - 5*u^76 - u^74 + 2*u^72 + 4*u^
#                70 + 6*u^68 + 7*u^66 + 7*u^64 + 6*u^62 + 5*u^60 + 4*u^58 + 
#                2*u^56 + u^54 - u^50 - u^48 - u^46 - u^44, 
#              u^165 + 2*u^163 + 3*u^161 + 4*u^159 + 5*u^157 + 5*u^155 + 5*u^
#                153 + 4*u^151 + 2*u^149 - 3*u^145 - 6*u^143 - 9*u^141 - 12*u^
#                139 - 14*u^137 - 15*u^135 - 16*u^133 - 14*u^131 - 12*u^129 - 
#                9*u^127 - 4*u^125 + 5*u^121 + 11*u^119 + 14*u^117 + 18*u^
#                115 + 21*u^113 + 21*u^111 + 22*u^109 + 20*u^107 + 16*u^105 + 
#                13*u^103 + 7*u^101 + 2*u^99 - 2*u^97 - 8*u^95 - 11*u^93 - 
#                14*u^91 - 17*u^89 - 16*u^87 - 16*u^85 - 15*u^83 - 11*u^81 - 
#                9*u^79 - 5*u^77 - u^75 + u^73 + 4*u^71 + 6*u^69 + 6*u^67 + 
#                7*u^65 + 6*u^63 + 5*u^61 + 4*u^59 + 2*u^57 + u^55 - u^51 - u^
#                49 - u^47 - u^45, u^164 + 2*u^162 + 2*u^160 + 4*u^158 + 4*u^
#                156 + 4*u^154 + 5*u^152 + 4*u^150 + 2*u^148 + 2*u^146 - 2*u^
#                144 - 4*u^142 - 6*u^140 - 10*u^138 - 11*u^136 - 13*u^134 - 
#                15*u^132 - 13*u^130 - 13*u^128 - 11*u^126 - 6*u^124 - 5*u^
#                122 + u^120 + 6*u^118 + 8*u^116 + 14*u^114 + 17*u^112 + 17*u^
#                110 + 21*u^108 + 19*u^106 + 17*u^104 + 16*u^102 + 11*u^100 + 
#                7*u^98 + 4*u^96 - 2*u^94 - 5*u^92 - 9*u^90 - 12*u^88 - 13*u^
#                86 - 14*u^84 - 14*u^82 - 12*u^80 - 11*u^78 - 8*u^76 - 5*u^
#                74 - 3*u^72 + 2*u^68 + 3*u^66 + 4*u^64 + 5*u^62 + 4*u^60 + 
#                4*u^58 + 3*u^56 + 2*u^54 + u^52 + u^50, 
#              u^164 + u^162 + u^160 + 2*u^158 + u^156 + u^154 + u^152 - u^
#                150 - u^148 - 2*u^146 - 4*u^144 - 3*u^142 - 5*u^140 - 5*u^
#                138 - 3*u^136 - 5*u^134 - 2*u^132 - u^128 + 4*u^126 + 4*u^
#                124 + 4*u^122 + 9*u^120 + 6*u^118 + 7*u^116 + 9*u^114 + 4*u^
#                112 + 6*u^110 + 4*u^108 - u^106 + u^104 - 4*u^102 - 6*u^100 - 
#                4*u^98 - 9*u^96 - 7*u^94 - 6*u^92 - 9*u^90 - 4*u^88 - 4*u^
#                86 - 4*u^84 + u^82 + 2*u^78 + 5*u^76 + 3*u^74 + 5*u^72 + 5*u^
#                70 + 3*u^68 + 4*u^66 + 2*u^64 + u^62 + u^60 - u^58 - u^56 - u^
#                54 - 2*u^52 - u^50 - u^48 - u^46, 
#              u^154 + u^148 + u^142 - u^140 - 2*u^134 - u^130 - 2*u^128 - 2*u^
#                124 - u^122 + u^120 - 2*u^118 + u^116 + u^114 - u^112 + 3*u^
#                110 + u^108 + 4*u^104 + u^100 + 3*u^98 - u^96 + u^94 + u^92 - 
#                2*u^90 + u^88 - u^86 - 2*u^84 - 2*u^80 - u^78 - 2*u^74 - u^
#                68 + u^66 + u^60 + u^54 ], 
#          [ u^165 + 2*u^163 + 3*u^161 + 4*u^159 + 5*u^157 + 5*u^155 + 5*u^
#                153 + 4*u^151 + 2*u^149 - 3*u^145 - 6*u^143 - 9*u^141 - 12*u^
#                139 - 14*u^137 - 15*u^135 - 16*u^133 - 14*u^131 - 12*u^129 - 
#                9*u^127 - 4*u^125 + 5*u^121 + 11*u^119 + 14*u^117 + 18*u^
#                115 + 21*u^113 + 21*u^111 + 22*u^109 + 20*u^107 + 16*u^105 + 
#                13*u^103 + 7*u^101 + 2*u^99 - 2*u^97 - 8*u^95 - 11*u^93 - 
#                14*u^91 - 17*u^89 - 16*u^87 - 16*u^85 - 15*u^83 - 11*u^81 - 
#                9*u^79 - 5*u^77 - u^75 + u^73 + 4*u^71 + 6*u^69 + 6*u^67 + 
#                7*u^65 + 6*u^63 + 5*u^61 + 4*u^59 + 2*u^57 + u^55 - u^51 - u^
#                49 - u^47 - u^45, u^166 + 2*u^164 + 3*u^162 + 4*u^160 + 5*u^
#                158 + 5*u^156 + 5*u^154 + 4*u^152 + 2*u^150 - 3*u^146 - 6*u^
#                144 - 9*u^142 - 12*u^140 - 14*u^138 - 15*u^136 - 16*u^134 - 
#                14*u^132 - 12*u^130 - 9*u^128 - 4*u^126 + 5*u^122 + 11*u^
#                120 + 14*u^118 + 18*u^116 + 21*u^114 + 21*u^112 + 22*u^110 + 
#                20*u^108 + 16*u^106 + 13*u^104 + 7*u^102 + 2*u^100 - 2*u^98 - 
#                8*u^96 - 11*u^94 - 14*u^92 - 17*u^90 - 16*u^88 - 16*u^86 - 
#                15*u^84 - 11*u^82 - 9*u^80 - 5*u^78 - u^76 + u^74 + 4*u^72 + 
#                6*u^70 + 6*u^68 + 7*u^66 + 6*u^64 + 5*u^62 + 4*u^60 + 2*u^
#                58 + u^56 - u^52 - u^50 - u^48 - u^46, 
#              u^165 + u^163 + 2*u^161 + 3*u^159 + 3*u^157 + 4*u^155 + 4*u^
#                153 + 3*u^151 + 3*u^149 + u^147 - u^145 - 2*u^143 - 6*u^141 - 
#                7*u^139 - 9*u^137 - 12*u^135 - 11*u^133 - 12*u^131 - 12*u^
#                129 - 8*u^127 - 8*u^125 - 4*u^123 + u^121 + 2*u^119 + 8*u^
#                117 + 11*u^115 + 12*u^113 + 17*u^111 + 16*u^109 + 16*u^107 + 
#                17*u^105 + 12*u^103 + 11*u^101 + 8*u^99 + 2*u^97 + u^95 - 4*u^
#                93 - 8*u^91 - 8*u^89 - 12*u^87 - 12*u^85 - 11*u^83 - 12*u^
#                81 - 9*u^79 - 7*u^77 - 6*u^75 - 2*u^73 - u^71 + u^69 + 3*u^
#                67 + 3*u^65 + 4*u^63 + 4*u^61 + 3*u^59 + 3*u^57 + 2*u^55 + u^
#                53 + u^51, u^165 + u^163 + u^161 + 2*u^159 + u^157 + u^
#                155 + u^153 - u^151 - u^149 - 2*u^147 - 4*u^145 - 3*u^143 - 
#                5*u^141 - 5*u^139 - 3*u^137 - 5*u^135 - 2*u^133 - u^129 + 4*u^
#                127 + 4*u^125 + 4*u^123 + 9*u^121 + 6*u^119 + 7*u^117 + 9*u^
#                115 + 4*u^113 + 6*u^111 + 4*u^109 - u^107 + u^105 - 4*u^103 - 
#                6*u^101 - 4*u^99 - 9*u^97 - 7*u^95 - 6*u^93 - 9*u^91 - 4*u^
#                89 - 4*u^87 - 4*u^85 + u^83 + 2*u^79 + 5*u^77 + 3*u^75 + 5*u^
#                73 + 5*u^71 + 3*u^69 + 4*u^67 + 2*u^65 + u^63 + u^61 - u^
#                59 - u^57 - u^55 - 2*u^53 - u^51 - u^49 - u^47, 0*u^0 ], 
#          [ u^164 + 2*u^162 + 2*u^160 + 4*u^158 + 4*u^156 + 4*u^154 + 5*u^
#                152 + 4*u^150 + 2*u^148 + 2*u^146 - 2*u^144 - 4*u^142 - 6*u^
#                140 - 10*u^138 - 11*u^136 - 13*u^134 - 15*u^132 - 13*u^130 - 
#                13*u^128 - 11*u^126 - 6*u^124 - 5*u^122 + u^120 + 6*u^118 + 
#                8*u^116 + 14*u^114 + 17*u^112 + 17*u^110 + 21*u^108 + 19*u^
#                106 + 17*u^104 + 16*u^102 + 11*u^100 + 7*u^98 + 4*u^96 - 2*u^
#                94 - 5*u^92 - 9*u^90 - 12*u^88 - 13*u^86 - 14*u^84 - 14*u^
#                82 - 12*u^80 - 11*u^78 - 8*u^76 - 5*u^74 - 3*u^72 + 2*u^68 + 
#                3*u^66 + 4*u^64 + 5*u^62 + 4*u^60 + 4*u^58 + 3*u^56 + 2*u^
#                54 + u^52 + u^50, u^165 + u^163 + 2*u^161 + 3*u^159 + 3*u^
#                157 + 4*u^155 + 4*u^153 + 3*u^151 + 3*u^149 + u^147 - u^145 - 
#                2*u^143 - 6*u^141 - 7*u^139 - 9*u^137 - 12*u^135 - 11*u^133 - 
#                12*u^131 - 12*u^129 - 8*u^127 - 8*u^125 - 4*u^123 + u^121 + 
#                2*u^119 + 8*u^117 + 11*u^115 + 12*u^113 + 17*u^111 + 16*u^
#                109 + 16*u^107 + 17*u^105 + 12*u^103 + 11*u^101 + 8*u^99 + 
#                2*u^97 + u^95 - 4*u^93 - 8*u^91 - 8*u^89 - 12*u^87 - 12*u^
#                85 - 11*u^83 - 12*u^81 - 9*u^79 - 7*u^77 - 6*u^75 - 2*u^
#                73 - u^71 + u^69 + 3*u^67 + 3*u^65 + 4*u^63 + 4*u^61 + 3*u^
#                59 + 3*u^57 + 2*u^55 + u^53 + u^51, 
#              u^166 + u^164 + 2*u^162 + 3*u^160 + 4*u^158 + 4*u^156 + 6*u^
#                154 + 4*u^152 + 4*u^150 + 3*u^148 - 2*u^144 - 4*u^142 - 9*u^
#                140 - 10*u^138 - 13*u^136 - 16*u^134 - 14*u^132 - 16*u^130 - 
#                14*u^128 - 10*u^126 - 9*u^124 - 3*u^122 + 3*u^120 + 4*u^118 + 
#                13*u^116 + 15*u^114 + 17*u^112 + 23*u^110 + 21*u^108 + 20*u^
#                106 + 22*u^104 + 14*u^102 + 13*u^100 + 9*u^98 + u^96 - u^94 - 
#                6*u^92 - 12*u^90 - 11*u^88 - 16*u^86 - 16*u^84 - 14*u^82 - 
#                15*u^80 - 11*u^78 - 8*u^76 - 7*u^74 - 2*u^72 + u^68 + 5*u^
#                66 + 4*u^64 + 5*u^62 + 5*u^60 + 4*u^58 + 3*u^56 + 3*u^54 + u^
#                52 + u^50, 0*u^0, u^162 + u^158 + u^156 + u^152 + u^150 - u^
#                148 + u^146 - u^144 - 2*u^142 - 3*u^138 - 2*u^136 - u^134 - 
#                4*u^132 - u^130 - u^128 - 3*u^126 + 2*u^124 - u^122 + 4*u^
#                118 + 3*u^114 + 5*u^112 + 5*u^108 + 3*u^106 + 4*u^102 - u^
#                98 + 2*u^96 - 3*u^94 - u^92 - u^90 - 4*u^88 - u^86 - 2*u^84 - 
#                3*u^82 - 2*u^78 - u^76 + u^74 - u^72 + u^70 + u^68 + u^64 + u^
#                62 + u^58 ], 
#          [ u^164 + u^162 + u^160 + 2*u^158 + u^156 + u^154 + u^152 - u^
#                150 - u^148 - 2*u^146 - 4*u^144 - 3*u^142 - 5*u^140 - 5*u^
#                138 - 3*u^136 - 5*u^134 - 2*u^132 - u^128 + 4*u^126 + 4*u^
#                124 + 4*u^122 + 9*u^120 + 6*u^118 + 7*u^116 + 9*u^114 + 4*u^
#                112 + 6*u^110 + 4*u^108 - u^106 + u^104 - 4*u^102 - 6*u^100 - 
#                4*u^98 - 9*u^96 - 7*u^94 - 6*u^92 - 9*u^90 - 4*u^88 - 4*u^
#                86 - 4*u^84 + u^82 + 2*u^78 + 5*u^76 + 3*u^74 + 5*u^72 + 5*u^
#                70 + 3*u^68 + 4*u^66 + 2*u^64 + u^62 + u^60 - u^58 - u^56 - u^
#                54 - 2*u^52 - u^50 - u^48 - u^46, 
#              u^165 + u^163 + u^161 + 2*u^159 + u^157 + u^155 + u^153 - u^
#                151 - u^149 - 2*u^147 - 4*u^145 - 3*u^143 - 5*u^141 - 5*u^
#                139 - 3*u^137 - 5*u^135 - 2*u^133 - u^129 + 4*u^127 + 4*u^
#                125 + 4*u^123 + 9*u^121 + 6*u^119 + 7*u^117 + 9*u^115 + 4*u^
#                113 + 6*u^111 + 4*u^109 - u^107 + u^105 - 4*u^103 - 6*u^101 - 
#                4*u^99 - 9*u^97 - 7*u^95 - 6*u^93 - 9*u^91 - 4*u^89 - 4*u^
#                87 - 4*u^85 + u^83 + 2*u^79 + 5*u^77 + 3*u^75 + 5*u^73 + 5*u^
#                71 + 3*u^69 + 4*u^67 + 2*u^65 + u^63 + u^61 - u^59 - u^57 - u^
#                55 - 2*u^53 - u^51 - u^49 - u^47, 0*u^0, 
#              u^166 + u^164 + u^162 + 2*u^160 + u^158 + u^156 + u^154 - u^
#                152 - u^150 - 2*u^148 - 4*u^146 - 3*u^144 - 5*u^142 - 5*u^
#                140 - 3*u^138 - 5*u^136 - 2*u^134 - u^130 + 4*u^128 + 4*u^
#                126 + 4*u^124 + 9*u^122 + 6*u^120 + 7*u^118 + 9*u^116 + 4*u^
#                114 + 6*u^112 + 4*u^110 - u^108 + u^106 - 4*u^104 - 6*u^102 - 
#                4*u^100 - 9*u^98 - 7*u^96 - 6*u^94 - 9*u^92 - 4*u^90 - 4*u^
#                88 - 4*u^86 + u^84 + 2*u^80 + 5*u^78 + 3*u^76 + 5*u^74 + 5*u^
#                72 + 3*u^70 + 4*u^68 + 2*u^66 + u^64 + u^62 - u^60 - u^58 - u^
#                56 - 2*u^54 - u^52 - u^50 - u^48, 0*u^0 ], 
#          [ u^154 + u^148 + u^142 - u^140 - 2*u^134 - u^130 - 2*u^128 - 2*u^
#                124 - u^122 + u^120 - 2*u^118 + u^116 + u^114 - u^112 + 3*u^
#                110 + u^108 + 4*u^104 + u^100 + 3*u^98 - u^96 + u^94 + u^92 - 
#                2*u^90 + u^88 - u^86 - 2*u^84 - 2*u^80 - u^78 - 2*u^74 - u^
#                68 + u^66 + u^60 + u^54, 0*u^0, 
#              u^162 + u^158 + u^156 + u^152 + u^150 - u^148 + u^146 - u^144 - 
#                2*u^142 - 3*u^138 - 2*u^136 - u^134 - 4*u^132 - u^130 - u^
#                128 - 3*u^126 + 2*u^124 - u^122 + 4*u^118 + 3*u^114 + 5*u^
#                112 + 5*u^108 + 3*u^106 + 4*u^102 - u^98 + 2*u^96 - 3*u^
#                94 - u^92 - u^90 - 4*u^88 - u^86 - 2*u^84 - 3*u^82 - 2*u^
#                78 - u^76 + u^74 - u^72 + u^70 + u^68 + u^64 + u^62 + u^58, 
#              0*u^0, u^166 + u^160 + u^154 - u^152 - 2*u^146 - u^142 - 2*u^
#                140 - 2*u^136 - u^134 + u^132 - 2*u^130 + u^128 + u^126 - u^
#                124 + 3*u^122 + u^120 + 4*u^116 + u^112 + 3*u^110 - u^108 + u^
#                106 + u^104 - 2*u^102 + u^100 - u^98 - 2*u^96 - 2*u^92 - u^
#                90 - 2*u^86 - u^80 + u^78 + u^72 + u^66 ] ], 
#      [ [ u^168 + u^162 + u^156 - u^154 - 2*u^148 - u^144 - 2*u^142 - 2*u^
#                138 - u^136 + u^134 - 2*u^132 + u^130 + u^128 - u^126 + 3*u^
#                124 + u^122 + 4*u^118 + u^114 + 3*u^112 - u^110 + u^108 + u^
#                106 - 2*u^104 + u^102 - u^100 - 2*u^98 - 2*u^94 - u^92 - 2*u^
#                88 - u^82 + u^80 + u^74 + u^68 ] ], 
#      [ [ u^176 + 2*u^172 + 2*u^170 + 3*u^168 + 3*u^166 + 3*u^164 + u^162 + u^
#                160 - 2*u^158 - 4*u^156 - 5*u^154 - 9*u^152 - 8*u^150 - 10*u^
#                148 - 10*u^146 - 8*u^144 - 6*u^142 - 5*u^140 + 2*u^138 + 2*u^
#                136 + 8*u^134 + 12*u^132 + 12*u^130 + 16*u^128 + 17*u^126 + 
#                14*u^124 + 16*u^122 + 11*u^120 + 7*u^118 + 6*u^116 - 2*u^
#                114 - 5*u^112 - 8*u^110 - 15*u^108 - 14*u^106 - 16*u^104 - 
#                18*u^102 - 13*u^100 - 14*u^98 - 11*u^96 - 5*u^94 - 4*u^92 + u^
#                90 + 6*u^88 + 5*u^86 + 11*u^84 + 10*u^82 + 10*u^80 + 10*u^
#                78 + 8*u^76 + 5*u^74 + 5*u^72 - 2*u^66 - 4*u^64 - 3*u^62 - 
#                4*u^60 - 3*u^58 - 2*u^56 - u^54 - u^52 + u^50 + u^46, 
#              u^174 + 2*u^172 + 2*u^170 + 4*u^168 + 3*u^166 + 3*u^164 + 3*u^
#                162 - u^158 - 3*u^156 - 7*u^154 - 7*u^152 - 10*u^150 - 12*u^
#                148 - 9*u^146 - 12*u^144 - 8*u^142 - 4*u^140 - 4*u^138 + 4*u^
#                136 + 7*u^134 + 8*u^132 + 17*u^130 + 15*u^128 + 17*u^126 + 
#                21*u^124 + 14*u^122 + 16*u^120 + 13*u^118 + 4*u^116 + 5*u^
#                114 - 4*u^112 - 10*u^110 - 9*u^108 - 18*u^106 - 17*u^104 - 
#                16*u^102 - 21*u^100 - 14*u^98 - 13*u^96 - 12*u^94 - 3*u^92 - 
#                3*u^90 + 2*u^88 + 8*u^86 + 7*u^84 + 11*u^82 + 12*u^80 + 9*u^
#                78 + 11*u^76 + 7*u^74 + 5*u^72 + 4*u^70 - u^66 - 2*u^64 - 4*u^
#                62 - 3*u^60 - 3*u^58 - 3*u^56 - u^54 - u^52, 
#              u^174 + u^172 + 3*u^170 + 2*u^168 + 3*u^166 + 2*u^164 + u^162 - 
#                2*u^158 - 4*u^156 - 5*u^154 - 7*u^152 - 8*u^150 - 7*u^148 - 
#                9*u^146 - 5*u^144 - 5*u^142 - 2*u^140 + 2*u^138 + 4*u^136 + 
#                7*u^134 + 11*u^132 + 10*u^130 + 14*u^128 + 13*u^126 + 11*u^
#                124 + 12*u^122 + 7*u^120 + 5*u^118 + 3*u^116 - 4*u^114 - 5*u^
#                112 - 9*u^110 - 13*u^108 - 11*u^106 - 14*u^104 - 13*u^102 - 
#                10*u^100 - 11*u^98 - 6*u^96 - 3*u^94 - 2*u^92 + 4*u^90 + 4*u^
#                88 + 7*u^86 + 9*u^84 + 8*u^82 + 8*u^80 + 8*u^78 + 4*u^76 + 
#                5*u^74 + u^72 - u^68 - 3*u^66 - 3*u^64 - 3*u^62 - 3*u^60 - 
#                2*u^58 - u^56 - u^54 + u^52 + u^48, 
#              u^172 + u^168 - u^160 - u^158 - u^156 - 2*u^154 - u^152 - u^
#                150 - 2*u^148 + u^146 - u^144 + u^142 + 2*u^140 + u^138 + 3*u^
#                136 + 3*u^134 + u^132 + 4*u^130 + u^128 + u^126 + 2*u^124 - 
#                2*u^122 - u^118 - 4*u^116 - u^114 - 4*u^112 - 4*u^110 - u^
#                108 - 4*u^106 - u^104 - 2*u^100 + 2*u^98 + u^96 + u^94 + 4*u^
#                92 + u^90 + 3*u^88 + 3*u^86 + u^84 + 2*u^82 + u^80 - u^78 + u^
#                76 - 2*u^74 - u^72 - u^70 - 2*u^68 - u^66 - u^64 - u^62 + u^
#                54 + u^50, u^170 + u^166 + u^164 + u^160 - u^156 - 2*u^152 - 
#                2*u^150 - u^148 - 4*u^146 - u^144 - 2*u^142 - 3*u^140 + u^
#                138 - u^136 + 4*u^132 + 4*u^128 + 5*u^126 + u^124 + 6*u^122 + 
#                3*u^120 + u^118 + 5*u^116 - u^114 + u^110 - 5*u^108 - u^106 - 
#                3*u^104 - 6*u^102 - u^100 - 5*u^98 - 4*u^96 - 4*u^92 + u^
#                88 - u^86 + 3*u^84 + 2*u^82 + u^80 + 4*u^78 + u^76 + 2*u^74 + 
#                2*u^72 + u^68 - u^64 - u^60 - u^58 - u^54 ], 
#          [ u^174 + 2*u^172 + 2*u^170 + 4*u^168 + 3*u^166 + 3*u^164 + 3*u^
#                162 - u^158 - 3*u^156 - 7*u^154 - 7*u^152 - 10*u^150 - 12*u^
#                148 - 9*u^146 - 12*u^144 - 8*u^142 - 4*u^140 - 4*u^138 + 4*u^
#                136 + 7*u^134 + 8*u^132 + 17*u^130 + 15*u^128 + 17*u^126 + 
#                21*u^124 + 14*u^122 + 16*u^120 + 13*u^118 + 4*u^116 + 5*u^
#                114 - 4*u^112 - 10*u^110 - 9*u^108 - 18*u^106 - 17*u^104 - 
#                16*u^102 - 21*u^100 - 14*u^98 - 13*u^96 - 12*u^94 - 3*u^92 - 
#                3*u^90 + 2*u^88 + 8*u^86 + 7*u^84 + 11*u^82 + 12*u^80 + 9*u^
#                78 + 11*u^76 + 7*u^74 + 5*u^72 + 4*u^70 - u^66 - 2*u^64 - 4*u^
#                62 - 3*u^60 - 3*u^58 - 3*u^56 - u^54 - u^52, 
#              u^176 + u^174 + 3*u^172 + 4*u^170 + 4*u^168 + 5*u^166 + 4*u^
#                164 + 2*u^162 + u^160 - 3*u^158 - 6*u^156 - 8*u^154 - 13*u^
#                152 - 13*u^150 - 15*u^148 - 16*u^146 - 12*u^144 - 11*u^142 - 
#                7*u^140 + u^138 + 3*u^136 + 11*u^134 + 17*u^132 + 18*u^130 + 
#                25*u^128 + 25*u^126 + 23*u^124 + 25*u^122 + 18*u^120 + 14*u^
#                118 + 10*u^116 - u^114 - 5*u^112 - 12*u^110 - 20*u^108 - 20*u^
#                106 - 25*u^104 - 26*u^102 - 22*u^100 - 23*u^98 - 17*u^96 - 
#                11*u^94 - 8*u^92 + 5*u^88 + 8*u^86 + 14*u^84 + 14*u^82 + 15*u^
#                80 + 15*u^78 + 12*u^76 + 10*u^74 + 7*u^72 + 3*u^70 + u^68 - 
#                2*u^66 - 4*u^64 - 4*u^62 - 5*u^60 - 4*u^58 - 3*u^56 - 2*u^
#                54 - u^52, u^174 + u^172 + 2*u^170 + 2*u^168 + 2*u^166 + 2*u^
#                164 + u^162 - u^158 - 3*u^156 - 4*u^154 - 5*u^152 - 7*u^150 - 
#                6*u^148 - 7*u^146 - 6*u^144 - 4*u^142 - 3*u^140 + 3*u^136 + 
#                4*u^134 + 8*u^132 + 9*u^130 + 10*u^128 + 12*u^126 + 10*u^
#                124 + 10*u^122 + 9*u^120 + 5*u^118 + 4*u^116 - 4*u^112 - 5*u^
#                110 - 9*u^108 - 10*u^106 - 10*u^104 - 12*u^102 - 10*u^100 - 
#                9*u^98 - 8*u^96 - 4*u^94 - 3*u^92 + 3*u^88 + 4*u^86 + 6*u^
#                84 + 7*u^82 + 6*u^80 + 7*u^78 + 5*u^76 + 4*u^74 + 3*u^72 + u^
#                70 - u^66 - 2*u^64 - 2*u^62 - 2*u^60 - 2*u^58 - u^56 - u^54, 
#              0*u^0, u^174 + u^172 + u^170 + 2*u^168 + u^166 + u^164 + u^
#                162 - u^160 - u^158 - 2*u^156 - 4*u^154 - 3*u^152 - 5*u^150 - 
#                5*u^148 - 3*u^146 - 5*u^144 - 2*u^142 - u^138 + 4*u^136 + 4*u^
#                134 + 4*u^132 + 9*u^130 + 6*u^128 + 7*u^126 + 9*u^124 + 4*u^
#                122 + 6*u^120 + 4*u^118 - u^116 + u^114 - 4*u^112 - 6*u^110 - 
#                4*u^108 - 9*u^106 - 7*u^104 - 6*u^102 - 9*u^100 - 4*u^98 - 
#                4*u^96 - 4*u^94 + u^92 + 2*u^88 + 5*u^86 + 3*u^84 + 5*u^82 + 
#                5*u^80 + 3*u^78 + 4*u^76 + 2*u^74 + u^72 + u^70 - u^68 - u^
#                66 - u^64 - 2*u^62 - u^60 - u^58 - u^56 ], 
#          [ u^174 + u^172 + 3*u^170 + 2*u^168 + 3*u^166 + 2*u^164 + u^162 - 
#                2*u^158 - 4*u^156 - 5*u^154 - 7*u^152 - 8*u^150 - 7*u^148 - 
#                9*u^146 - 5*u^144 - 5*u^142 - 2*u^140 + 2*u^138 + 4*u^136 + 
#                7*u^134 + 11*u^132 + 10*u^130 + 14*u^128 + 13*u^126 + 11*u^
#                124 + 12*u^122 + 7*u^120 + 5*u^118 + 3*u^116 - 4*u^114 - 5*u^
#                112 - 9*u^110 - 13*u^108 - 11*u^106 - 14*u^104 - 13*u^102 - 
#                10*u^100 - 11*u^98 - 6*u^96 - 3*u^94 - 2*u^92 + 4*u^90 + 4*u^
#                88 + 7*u^86 + 9*u^84 + 8*u^82 + 8*u^80 + 8*u^78 + 4*u^76 + 
#                5*u^74 + u^72 - u^68 - 3*u^66 - 3*u^64 - 3*u^62 - 3*u^60 - 
#                2*u^58 - u^56 - u^54 + u^52 + u^48, 
#              u^174 + u^172 + 2*u^170 + 2*u^168 + 2*u^166 + 2*u^164 + u^
#                162 - u^158 - 3*u^156 - 4*u^154 - 5*u^152 - 7*u^150 - 6*u^
#                148 - 7*u^146 - 6*u^144 - 4*u^142 - 3*u^140 + 3*u^136 + 4*u^
#                134 + 8*u^132 + 9*u^130 + 10*u^128 + 12*u^126 + 10*u^124 + 
#                10*u^122 + 9*u^120 + 5*u^118 + 4*u^116 - 4*u^112 - 5*u^110 - 
#                9*u^108 - 10*u^106 - 10*u^104 - 12*u^102 - 10*u^100 - 9*u^
#                98 - 8*u^96 - 4*u^94 - 3*u^92 + 3*u^88 + 4*u^86 + 6*u^84 + 
#                7*u^82 + 6*u^80 + 7*u^78 + 5*u^76 + 4*u^74 + 3*u^72 + u^
#                70 - u^66 - 2*u^64 - 2*u^62 - 2*u^60 - 2*u^58 - u^56 - u^54, 
#              u^176 + u^174 + 3*u^172 + 2*u^170 + 3*u^168 + 2*u^166 + u^164 - 
#                2*u^160 - 4*u^158 - 5*u^156 - 7*u^154 - 8*u^152 - 7*u^150 - 
#                9*u^148 - 5*u^146 - 5*u^144 - 2*u^142 + 2*u^140 + 4*u^138 + 
#                7*u^136 + 11*u^134 + 10*u^132 + 14*u^130 + 13*u^128 + 11*u^
#                126 + 12*u^124 + 7*u^122 + 5*u^120 + 3*u^118 - 4*u^116 - 5*u^
#                114 - 9*u^112 - 13*u^110 - 11*u^108 - 14*u^106 - 13*u^104 - 
#                10*u^102 - 11*u^100 - 6*u^98 - 3*u^96 - 2*u^94 + 4*u^92 + 4*u^
#                90 + 7*u^88 + 9*u^86 + 8*u^84 + 8*u^82 + 8*u^80 + 4*u^78 + 
#                5*u^76 + u^74 - u^70 - 3*u^68 - 3*u^66 - 3*u^64 - 3*u^62 - 
#                2*u^60 - u^58 - u^56 + u^54 + u^50, 
#              u^174 + u^170 - u^162 - u^160 - u^158 - 2*u^156 - u^154 - u^
#                152 - 2*u^150 + u^148 - u^146 + u^144 + 2*u^142 + u^140 + 3*u^
#                138 + 3*u^136 + u^134 + 4*u^132 + u^130 + u^128 + 2*u^126 - 
#                2*u^124 - u^120 - 4*u^118 - u^116 - 4*u^114 - 4*u^112 - u^
#                110 - 4*u^108 - u^106 - 2*u^102 + 2*u^100 + u^98 + u^96 + 4*u^
#                94 + u^92 + 3*u^90 + 3*u^88 + u^86 + 2*u^84 + u^82 - u^80 + u^
#                78 - 2*u^76 - u^74 - u^72 - 2*u^70 - u^68 - u^66 - u^64 + u^
#                56 + u^52, 0*u^0 ], 
#          [ u^172 + u^168 - u^160 - u^158 - u^156 - 2*u^154 - u^152 - u^150 - 
#                2*u^148 + u^146 - u^144 + u^142 + 2*u^140 + u^138 + 3*u^136 + 
#                3*u^134 + u^132 + 4*u^130 + u^128 + u^126 + 2*u^124 - 2*u^
#                122 - u^118 - 4*u^116 - u^114 - 4*u^112 - 4*u^110 - u^108 - 
#                4*u^106 - u^104 - 2*u^100 + 2*u^98 + u^96 + u^94 + 4*u^92 + u^
#                90 + 3*u^88 + 3*u^86 + u^84 + 2*u^82 + u^80 - u^78 + u^76 - 
#                2*u^74 - u^72 - u^70 - 2*u^68 - u^66 - u^64 - u^62 + u^54 + u^
#                50, 0*u^0, u^174 + u^170 - u^162 - u^160 - u^158 - 2*u^
#                156 - u^154 - u^152 - 2*u^150 + u^148 - u^146 + u^144 + 2*u^
#                142 + u^140 + 3*u^138 + 3*u^136 + u^134 + 4*u^132 + u^130 + u^
#                128 + 2*u^126 - 2*u^124 - u^120 - 4*u^118 - u^116 - 4*u^114 - 
#                4*u^112 - u^110 - 4*u^108 - u^106 - 2*u^102 + 2*u^100 + u^
#                98 + u^96 + 4*u^94 + u^92 + 3*u^90 + 3*u^88 + u^86 + 2*u^
#                84 + u^82 - u^80 + u^78 - 2*u^76 - u^74 - u^72 - 2*u^70 - u^
#                68 - u^66 - u^64 + u^56 + u^52, 
#              u^176 + u^172 - u^164 - u^162 - u^160 - 2*u^158 - u^156 - u^
#                154 - 2*u^152 + u^150 - u^148 + u^146 + 2*u^144 + u^142 + 3*u^
#                140 + 3*u^138 + u^136 + 4*u^134 + u^132 + u^130 + 2*u^128 - 
#                2*u^126 - u^122 - 4*u^120 - u^118 - 4*u^116 - 4*u^114 - u^
#                112 - 4*u^110 - u^108 - 2*u^104 + 2*u^102 + u^100 + u^98 + 
#                4*u^96 + u^94 + 3*u^92 + 3*u^90 + u^88 + 2*u^86 + u^84 - u^
#                82 + u^80 - 2*u^78 - u^76 - u^74 - 2*u^72 - u^70 - u^68 - u^
#                66 + u^58 + u^54, 0*u^0 ], 
#          [ u^170 + u^166 + u^164 + u^160 - u^156 - 2*u^152 - 2*u^150 - u^
#                148 - 4*u^146 - u^144 - 2*u^142 - 3*u^140 + u^138 - u^136 + 
#                4*u^132 + 4*u^128 + 5*u^126 + u^124 + 6*u^122 + 3*u^120 + u^
#                118 + 5*u^116 - u^114 + u^110 - 5*u^108 - u^106 - 3*u^104 - 
#                6*u^102 - u^100 - 5*u^98 - 4*u^96 - 4*u^92 + u^88 - u^86 + 
#                3*u^84 + 2*u^82 + u^80 + 4*u^78 + u^76 + 2*u^74 + 2*u^72 + u^
#                68 - u^64 - u^60 - u^58 - u^54, 
#              u^174 + u^172 + u^170 + 2*u^168 + u^166 + u^164 + u^162 - u^
#                160 - u^158 - 2*u^156 - 4*u^154 - 3*u^152 - 5*u^150 - 5*u^
#                148 - 3*u^146 - 5*u^144 - 2*u^142 - u^138 + 4*u^136 + 4*u^
#                134 + 4*u^132 + 9*u^130 + 6*u^128 + 7*u^126 + 9*u^124 + 4*u^
#                122 + 6*u^120 + 4*u^118 - u^116 + u^114 - 4*u^112 - 6*u^110 - 
#                4*u^108 - 9*u^106 - 7*u^104 - 6*u^102 - 9*u^100 - 4*u^98 - 
#                4*u^96 - 4*u^94 + u^92 + 2*u^88 + 5*u^86 + 3*u^84 + 5*u^82 + 
#                5*u^80 + 3*u^78 + 4*u^76 + 2*u^74 + u^72 + u^70 - u^68 - u^
#                66 - u^64 - 2*u^62 - u^60 - u^58 - u^56, 0*u^0, 0*u^0, 
#              u^176 + u^172 + u^170 + u^166 - u^162 - 2*u^158 - 2*u^156 - u^
#                154 - 4*u^152 - u^150 - 2*u^148 - 3*u^146 + u^144 - u^142 + 
#                4*u^138 + 4*u^134 + 5*u^132 + u^130 + 6*u^128 + 3*u^126 + u^
#                124 + 5*u^122 - u^120 + u^116 - 5*u^114 - u^112 - 3*u^110 - 
#                6*u^108 - u^106 - 5*u^104 - 4*u^102 - 4*u^98 + u^94 - u^92 + 
#                3*u^90 + 2*u^88 + u^86 + 4*u^84 + u^82 + 2*u^80 + 2*u^78 + u^
#                74 - u^70 - u^66 - u^64 - u^60 ] ], 
#      [ [ u^178 + u^176 + 2*u^174 + 2*u^172 + 2*u^170 + 2*u^168 + u^166 - u^
#                162 - 3*u^160 - 4*u^158 - 5*u^156 - 7*u^154 - 6*u^152 - 7*u^
#                150 - 6*u^148 - 4*u^146 - 3*u^144 + 3*u^140 + 4*u^138 + 8*u^
#                136 + 9*u^134 + 10*u^132 + 12*u^130 + 10*u^128 + 10*u^126 + 
#                9*u^124 + 5*u^122 + 4*u^120 - 4*u^116 - 5*u^114 - 9*u^112 - 
#                10*u^110 - 10*u^108 - 12*u^106 - 10*u^104 - 9*u^102 - 8*u^
#                100 - 4*u^98 - 3*u^96 + 3*u^92 + 4*u^90 + 6*u^88 + 7*u^86 + 
#                6*u^84 + 7*u^82 + 5*u^80 + 4*u^78 + 3*u^76 + u^74 - u^70 - 
#                2*u^68 - 2*u^66 - 2*u^64 - 2*u^62 - u^60 - u^58 ] ], 
#      [ [ u^180 + u^176 + u^174 + u^172 + 2*u^170 + u^168 + u^166 + u^164 - u^
#                158 - 3*u^156 - 2*u^154 - 4*u^152 - 4*u^150 - 4*u^148 - 5*u^
#                146 - 4*u^144 - 3*u^142 - 3*u^140 - u^138 + u^134 + 3*u^132 + 
#                4*u^130 + 5*u^128 + 6*u^126 + 6*u^124 + 6*u^122 + 6*u^120 + 
#                5*u^118 + 4*u^116 + 3*u^114 + u^112 - u^108 - 3*u^106 - 3*u^
#                104 - 4*u^102 - 5*u^100 - 4*u^98 - 4*u^96 - 4*u^94 - 2*u^92 - 
#                3*u^90 - u^88 + u^82 + u^80 + u^78 + 2*u^76 + u^74 + u^72 + u^
#                70 + u^66, 0*u^0, u^177 + u^175 + u^173 + 2*u^171 + u^169 + 
#                2*u^167 + 2*u^165 + u^161 - u^159 - 2*u^157 - u^155 - 5*u^
#                153 - 4*u^151 - 4*u^149 - 7*u^147 - 3*u^145 - 5*u^143 - 5*u^
#                141 - 3*u^137 + u^135 + 4*u^133 + u^131 + 7*u^129 + 6*u^127 + 
#                5*u^125 + 10*u^123 + 5*u^121 + 6*u^119 + 7*u^117 + u^115 + 
#                4*u^113 + u^111 - 3*u^109 - 5*u^105 - 5*u^103 - 3*u^101 - 7*u^
#                99 - 4*u^97 - 4*u^95 - 5*u^93 - u^91 - 2*u^89 - u^87 + u^85 + 
#                2*u^81 + 2*u^79 + u^77 + 2*u^75 + u^73 + u^71 + u^69 ], 
#          [ 0*u^0, u^180 + u^176 - u^168 - u^166 - u^164 - 2*u^162 - u^
#                160 - u^158 - 2*u^156 + u^154 - u^152 + u^150 + 2*u^148 + u^
#                146 + 3*u^144 + 3*u^142 + u^140 + 4*u^138 + u^136 + u^134 + 
#                2*u^132 - 2*u^130 - u^126 - 4*u^124 - u^122 - 4*u^120 - 4*u^
#                118 - u^116 - 4*u^114 - u^112 - 2*u^108 + 2*u^106 + u^104 + u^
#                102 + 4*u^100 + u^98 + 3*u^96 + 3*u^94 + u^92 + 2*u^90 + u^
#                88 - u^86 + u^84 - 2*u^82 - u^80 - u^78 - 2*u^76 - u^74 - u^
#                72 - u^70 + u^62 + u^58, 0*u^0 ], 
#          [ u^177 + u^175 + u^173 + 2*u^171 + u^169 + 2*u^167 + 2*u^165 + u^
#                161 - u^159 - 2*u^157 - u^155 - 5*u^153 - 4*u^151 - 4*u^149 - 
#                7*u^147 - 3*u^145 - 5*u^143 - 5*u^141 - 3*u^137 + u^135 + 4*u^
#                133 + u^131 + 7*u^129 + 6*u^127 + 5*u^125 + 10*u^123 + 5*u^
#                121 + 6*u^119 + 7*u^117 + u^115 + 4*u^113 + u^111 - 3*u^109 - 
#                5*u^105 - 5*u^103 - 3*u^101 - 7*u^99 - 4*u^97 - 4*u^95 - 5*u^
#                93 - u^91 - 2*u^89 - u^87 + u^85 + 2*u^81 + 2*u^79 + u^77 + 
#                2*u^75 + u^73 + u^71 + u^69, 0*u^0, 
#              u^180 + u^176 + u^174 + u^172 + 2*u^170 + u^168 + u^166 + u^
#                164 - u^158 - 3*u^156 - 2*u^154 - 4*u^152 - 4*u^150 - 4*u^
#                148 - 5*u^146 - 4*u^144 - 3*u^142 - 3*u^140 - u^138 + u^134 + 
#                3*u^132 + 4*u^130 + 5*u^128 + 6*u^126 + 6*u^124 + 6*u^122 + 
#                6*u^120 + 5*u^118 + 4*u^116 + 3*u^114 + u^112 - u^108 - 3*u^
#                106 - 3*u^104 - 4*u^102 - 5*u^100 - 4*u^98 - 4*u^96 - 4*u^
#                94 - 2*u^92 - 3*u^90 - u^88 + u^82 + u^80 + u^78 + 2*u^76 + u^
#                74 + u^72 + u^70 + u^66 ] ], 
#      [ [ u^184 + u^182 + u^180 + 2*u^178 - 3*u^170 - 2*u^168 - 3*u^166 - 5*u^
#                164 - 2*u^162 - 4*u^160 - 3*u^158 + u^156 - 2*u^154 + 3*u^
#                152 + 5*u^150 + 2*u^148 + 9*u^146 + 6*u^144 + 4*u^142 + 10*u^
#                140 + 2*u^138 + 3*u^136 + 5*u^134 - 5*u^132 - 3*u^128 - 10*u^
#                126 - 3*u^124 - 10*u^122 - 10*u^120 - 3*u^118 - 10*u^116 - 
#                3*u^114 - 5*u^110 + 5*u^108 + 3*u^106 + 2*u^104 + 10*u^102 + 
#                4*u^100 + 6*u^98 + 9*u^96 + 2*u^94 + 5*u^92 + 3*u^90 - 2*u^
#                88 + u^86 - 3*u^84 - 4*u^82 - 2*u^80 - 5*u^78 - 3*u^76 - 2*u^
#                74 - 3*u^72 + 2*u^64 + u^62 + u^60 + u^58, 
#              u^183 + u^181 + u^179 + u^177 - u^171 - 2*u^169 - 2*u^167 - 3*u^
#                165 - 3*u^163 - 2*u^161 - 3*u^159 - u^157 + 3*u^151 + 3*u^
#                149 + 4*u^147 + 6*u^145 + 4*u^143 + 5*u^141 + 5*u^139 + 2*u^
#                137 + 3*u^135 - 2*u^131 - u^129 - 5*u^127 - 5*u^125 - 5*u^
#                123 - 8*u^121 - 5*u^119 - 5*u^117 - 5*u^115 - u^113 - 2*u^
#                111 + 3*u^107 + 2*u^105 + 5*u^103 + 5*u^101 + 4*u^99 + 6*u^
#                97 + 4*u^95 + 3*u^93 + 3*u^91 - u^85 - 3*u^83 - 2*u^81 - 3*u^
#                79 - 3*u^77 - 2*u^75 - 2*u^73 - u^71 + u^65 + u^63 + u^61 + u^
#                59, u^181 + u^175 - u^173 - 2*u^167 - u^163 - 2*u^161 + u^
#                159 - 2*u^157 + 2*u^153 - 2*u^151 + 3*u^149 + 2*u^147 - u^
#                145 + 5*u^143 + 5*u^137 - 3*u^135 + u^133 + 2*u^131 - 5*u^
#                129 + 2*u^127 - 2*u^125 - 5*u^123 + 2*u^121 - 5*u^119 - 2*u^
#                117 + 2*u^115 - 5*u^113 + 2*u^111 + u^109 - 3*u^107 + 5*u^
#                105 + 5*u^99 - u^97 + 2*u^95 + 3*u^93 - 2*u^91 + 2*u^89 - 2*u^
#                85 + u^83 - 2*u^81 - u^79 - 2*u^75 - u^69 + u^67 + u^61 ], 
#          [ u^183 + u^181 + u^179 + u^177 - u^171 - 2*u^169 - 2*u^167 - 3*u^
#                165 - 3*u^163 - 2*u^161 - 3*u^159 - u^157 + 3*u^151 + 3*u^
#                149 + 4*u^147 + 6*u^145 + 4*u^143 + 5*u^141 + 5*u^139 + 2*u^
#                137 + 3*u^135 - 2*u^131 - u^129 - 5*u^127 - 5*u^125 - 5*u^
#                123 - 8*u^121 - 5*u^119 - 5*u^117 - 5*u^115 - u^113 - 2*u^
#                111 + 3*u^107 + 2*u^105 + 5*u^103 + 5*u^101 + 4*u^99 + 6*u^
#                97 + 4*u^95 + 3*u^93 + 3*u^91 - u^85 - 3*u^83 - 2*u^81 - 3*u^
#                79 - 3*u^77 - 2*u^75 - 2*u^73 - u^71 + u^65 + u^63 + u^61 + u^
#                59, u^184 + u^182 + u^180 + u^178 - u^172 - 2*u^170 - 2*u^
#                168 - 3*u^166 - 3*u^164 - 2*u^162 - 3*u^160 - u^158 + 3*u^
#                152 + 3*u^150 + 4*u^148 + 6*u^146 + 4*u^144 + 5*u^142 + 5*u^
#                140 + 2*u^138 + 3*u^136 - 2*u^132 - u^130 - 5*u^128 - 5*u^
#                126 - 5*u^124 - 8*u^122 - 5*u^120 - 5*u^118 - 5*u^116 - u^
#                114 - 2*u^112 + 3*u^108 + 2*u^106 + 5*u^104 + 5*u^102 + 4*u^
#                100 + 6*u^98 + 4*u^96 + 3*u^94 + 3*u^92 - u^86 - 3*u^84 - 2*u^
#                82 - 3*u^80 - 3*u^78 - 2*u^76 - 2*u^74 - u^72 + u^66 + u^
#                64 + u^62 + u^60, 0*u^0 ], 
#          [ u^181 + u^175 - u^173 - 2*u^167 - u^163 - 2*u^161 + u^159 - 2*u^
#                157 + 2*u^153 - 2*u^151 + 3*u^149 + 2*u^147 - u^145 + 5*u^
#                143 + 5*u^137 - 3*u^135 + u^133 + 2*u^131 - 5*u^129 + 2*u^
#                127 - 2*u^125 - 5*u^123 + 2*u^121 - 5*u^119 - 2*u^117 + 2*u^
#                115 - 5*u^113 + 2*u^111 + u^109 - 3*u^107 + 5*u^105 + 5*u^
#                99 - u^97 + 2*u^95 + 3*u^93 - 2*u^91 + 2*u^89 - 2*u^85 + u^
#                83 - 2*u^81 - u^79 - 2*u^75 - u^69 + u^67 + u^61, 0*u^0, 
#              u^184 + u^178 - u^176 - 2*u^170 - u^166 - 2*u^164 + u^162 - 2*u^
#                160 + 2*u^156 - 2*u^154 + 3*u^152 + 2*u^150 - u^148 + 5*u^
#                146 + 5*u^140 - 3*u^138 + u^136 + 2*u^134 - 5*u^132 + 2*u^
#                130 - 2*u^128 - 5*u^126 + 2*u^124 - 5*u^122 - 2*u^120 + 2*u^
#                118 - 5*u^116 + 2*u^114 + u^112 - 3*u^110 + 5*u^108 + 5*u^
#                102 - u^100 + 2*u^98 + 3*u^96 - 2*u^94 + 2*u^92 - 2*u^88 + u^
#                86 - 2*u^84 - u^82 - 2*u^78 - u^72 + u^70 + u^64 ] ], 
#      [ [ u^188 + u^186 + 2*u^184 + 3*u^182 + 2*u^180 + 3*u^178 + 2*u^176 - 
#                3*u^170 - 5*u^168 - 5*u^166 - 9*u^164 - 8*u^162 - 8*u^160 - 
#                10*u^158 - 5*u^156 - 5*u^154 - 3*u^152 + 4*u^150 + 3*u^148 + 
#                8*u^146 + 13*u^144 + 10*u^142 + 16*u^140 + 15*u^138 + 11*u^
#                136 + 15*u^134 + 8*u^132 + 5*u^130 + 5*u^128 - 5*u^126 - 5*u^
#                124 - 8*u^122 - 15*u^120 - 11*u^118 - 15*u^116 - 16*u^114 - 
#                10*u^112 - 13*u^110 - 8*u^108 - 3*u^106 - 4*u^104 + 3*u^102 + 
#                5*u^100 + 5*u^98 + 10*u^96 + 8*u^94 + 8*u^92 + 9*u^90 + 5*u^
#                88 + 5*u^86 + 3*u^84 - 2*u^78 - 3*u^76 - 2*u^74 - 3*u^72 - 
#                2*u^70 - u^68 - u^66, u^187 + 2*u^185 + 2*u^183 + 3*u^181 + 
#                3*u^179 + 2*u^177 + 2*u^175 - 2*u^171 - 3*u^169 - 6*u^167 - 
#                7*u^165 - 8*u^163 - 10*u^161 - 8*u^159 - 8*u^157 - 7*u^155 - 
#                2*u^153 - u^151 + 3*u^149 + 8*u^147 + 8*u^145 + 13*u^143 + 
#                15*u^141 + 13*u^139 + 16*u^137 + 13*u^135 + 10*u^133 + 10*u^
#                131 + 3*u^129 - 3*u^125 - 10*u^123 - 10*u^121 - 13*u^119 - 
#                16*u^117 - 13*u^115 - 15*u^113 - 13*u^111 - 8*u^109 - 8*u^
#                107 - 3*u^105 + u^103 + 2*u^101 + 7*u^99 + 8*u^97 + 8*u^95 + 
#                10*u^93 + 8*u^91 + 7*u^89 + 6*u^87 + 3*u^85 + 2*u^83 - 2*u^
#                79 - 2*u^77 - 3*u^75 - 3*u^73 - 2*u^71 - 2*u^69 - u^67 ], 
#          [ u^187 + 2*u^185 + 2*u^183 + 3*u^181 + 3*u^179 + 2*u^177 + 2*u^
#                175 - 2*u^171 - 3*u^169 - 6*u^167 - 7*u^165 - 8*u^163 - 10*u^
#                161 - 8*u^159 - 8*u^157 - 7*u^155 - 2*u^153 - u^151 + 3*u^
#                149 + 8*u^147 + 8*u^145 + 13*u^143 + 15*u^141 + 13*u^139 + 
#                16*u^137 + 13*u^135 + 10*u^133 + 10*u^131 + 3*u^129 - 3*u^
#                125 - 10*u^123 - 10*u^121 - 13*u^119 - 16*u^117 - 13*u^115 - 
#                15*u^113 - 13*u^111 - 8*u^109 - 8*u^107 - 3*u^105 + u^103 + 
#                2*u^101 + 7*u^99 + 8*u^97 + 8*u^95 + 10*u^93 + 8*u^91 + 7*u^
#                89 + 6*u^87 + 3*u^85 + 2*u^83 - 2*u^79 - 2*u^77 - 3*u^75 - 
#                3*u^73 - 2*u^71 - 2*u^69 - u^67, 
#              u^188 + u^186 + 2*u^184 + 3*u^182 + 2*u^180 + 3*u^178 + 2*u^
#                176 - 3*u^170 - 5*u^168 - 5*u^166 - 9*u^164 - 8*u^162 - 8*u^
#                160 - 10*u^158 - 5*u^156 - 5*u^154 - 3*u^152 + 4*u^150 + 3*u^
#                148 + 8*u^146 + 13*u^144 + 10*u^142 + 16*u^140 + 15*u^138 + 
#                11*u^136 + 15*u^134 + 8*u^132 + 5*u^130 + 5*u^128 - 5*u^126 - 
#                5*u^124 - 8*u^122 - 15*u^120 - 11*u^118 - 15*u^116 - 16*u^
#                114 - 10*u^112 - 13*u^110 - 8*u^108 - 3*u^106 - 4*u^104 + 3*u^
#                102 + 5*u^100 + 5*u^98 + 10*u^96 + 8*u^94 + 8*u^92 + 9*u^90 + 
#                5*u^88 + 5*u^86 + 3*u^84 - 2*u^78 - 3*u^76 - 2*u^74 - 3*u^
#                72 - 2*u^70 - u^68 - u^66 ] ], 
#      [ [ u^190 + u^186 + 2*u^184 + 2*u^180 + u^178 - u^176 + u^174 - 2*u^
#                172 - 3*u^170 - u^168 - 6*u^166 - 3*u^164 - 3*u^162 - 7*u^
#                160 - 3*u^156 - 3*u^154 + 5*u^152 - u^150 + 4*u^148 + 9*u^
#                146 + u^144 + 10*u^142 + 8*u^140 + 2*u^138 + 11*u^136 + 2*u^
#                134 + u^132 + 6*u^130 - 6*u^128 - u^126 - 2*u^124 - 11*u^
#                122 - 2*u^120 - 8*u^118 - 10*u^116 - u^114 - 9*u^112 - 4*u^
#                110 + u^108 - 5*u^106 + 3*u^104 + 3*u^102 + 7*u^98 + 3*u^96 + 
#                3*u^94 + 6*u^92 + u^90 + 3*u^88 + 2*u^86 - u^84 + u^82 - u^
#                80 - 2*u^78 - 2*u^74 - u^72 - u^68, 
#              u^187 + u^183 + u^181 + u^177 - u^173 - 2*u^169 - 2*u^167 - u^
#                165 - 4*u^163 - u^161 - 2*u^159 - 3*u^157 + u^155 - u^153 + 
#                4*u^149 + 4*u^145 + 5*u^143 + u^141 + 6*u^139 + 3*u^137 + u^
#                135 + 5*u^133 - u^131 + u^127 - 5*u^125 - u^123 - 3*u^121 - 
#                6*u^119 - u^117 - 5*u^115 - 4*u^113 - 4*u^109 + u^105 - u^
#                103 + 3*u^101 + 2*u^99 + u^97 + 4*u^95 + u^93 + 2*u^91 + 2*u^
#                89 + u^85 - u^81 - u^77 - u^75 - u^71, 
#              u^187 + u^183 + u^181 + u^177 - u^173 - 2*u^169 - 2*u^167 - u^
#                165 - 4*u^163 - u^161 - 2*u^159 - 3*u^157 + u^155 - u^153 + 
#                4*u^149 + 4*u^145 + 5*u^143 + u^141 + 6*u^139 + 3*u^137 + u^
#                135 + 5*u^133 - u^131 + u^127 - 5*u^125 - u^123 - 3*u^121 - 
#                6*u^119 - u^117 - 5*u^115 - 4*u^113 - 4*u^109 + u^105 - u^
#                103 + 3*u^101 + 2*u^99 + u^97 + 4*u^95 + u^93 + 2*u^91 + 2*u^
#                89 + u^85 - u^81 - u^77 - u^75 - u^71 ], 
#          [ u^187 + u^183 + u^181 + u^177 - u^173 - 2*u^169 - 2*u^167 - u^
#                165 - 4*u^163 - u^161 - 2*u^159 - 3*u^157 + u^155 - u^153 + 
#                4*u^149 + 4*u^145 + 5*u^143 + u^141 + 6*u^139 + 3*u^137 + u^
#                135 + 5*u^133 - u^131 + u^127 - 5*u^125 - u^123 - 3*u^121 - 
#                6*u^119 - u^117 - 5*u^115 - 4*u^113 - 4*u^109 + u^105 - u^
#                103 + 3*u^101 + 2*u^99 + u^97 + 4*u^95 + u^93 + 2*u^91 + 2*u^
#                89 + u^85 - u^81 - u^77 - u^75 - u^71, 
#              u^190 + u^186 + u^184 + u^180 - u^176 - 2*u^172 - 2*u^170 - u^
#                168 - 4*u^166 - u^164 - 2*u^162 - 3*u^160 + u^158 - u^156 + 
#                4*u^152 + 4*u^148 + 5*u^146 + u^144 + 6*u^142 + 3*u^140 + u^
#                138 + 5*u^136 - u^134 + u^130 - 5*u^128 - u^126 - 3*u^124 - 
#                6*u^122 - u^120 - 5*u^118 - 4*u^116 - 4*u^112 + u^108 - u^
#                106 + 3*u^104 + 2*u^102 + u^100 + 4*u^98 + u^96 + 2*u^94 + 
#                2*u^92 + u^88 - u^84 - u^80 - u^78 - u^74, 0*u^0 ], 
#          [ u^187 + u^183 + u^181 + u^177 - u^173 - 2*u^169 - 2*u^167 - u^
#                165 - 4*u^163 - u^161 - 2*u^159 - 3*u^157 + u^155 - u^153 + 
#                4*u^149 + 4*u^145 + 5*u^143 + u^141 + 6*u^139 + 3*u^137 + u^
#                135 + 5*u^133 - u^131 + u^127 - 5*u^125 - u^123 - 3*u^121 - 
#                6*u^119 - u^117 - 5*u^115 - 4*u^113 - 4*u^109 + u^105 - u^
#                103 + 3*u^101 + 2*u^99 + u^97 + 4*u^95 + u^93 + 2*u^91 + 2*u^
#                89 + u^85 - u^81 - u^77 - u^75 - u^71, 0*u^0, 
#              u^190 + u^186 + u^184 + u^180 - u^176 - 2*u^172 - 2*u^170 - u^
#                168 - 4*u^166 - u^164 - 2*u^162 - 3*u^160 + u^158 - u^156 + 
#                4*u^152 + 4*u^148 + 5*u^146 + u^144 + 6*u^142 + 3*u^140 + u^
#                138 + 5*u^136 - u^134 + u^130 - 5*u^128 - u^126 - 3*u^124 - 
#                6*u^122 - u^120 - 5*u^118 - 4*u^116 - 4*u^112 + u^108 - u^
#                106 + 3*u^104 + 2*u^102 + u^100 + 4*u^98 + u^96 + 2*u^94 + 
#                2*u^92 + u^88 - u^84 - u^80 - u^78 - u^74 ] ], 
#      [ [ u^192 + u^190 + 2*u^188 + u^186 + u^184 - u^180 - 2*u^178 - 3*u^
#                176 - 4*u^174 - 4*u^172 - 4*u^170 - 4*u^168 - 2*u^166 - 2*u^
#                164 + u^162 + 2*u^160 + 4*u^158 + 6*u^156 + 7*u^154 + 7*u^
#                152 + 8*u^150 + 6*u^148 + 6*u^146 + 4*u^144 + u^142 - 3*u^
#                138 - 5*u^136 - 6*u^134 - 9*u^132 - 9*u^130 - 9*u^128 - 9*u^
#                126 - 6*u^124 - 5*u^122 - 3*u^120 + u^116 + 4*u^114 + 6*u^
#                112 + 6*u^110 + 8*u^108 + 7*u^106 + 7*u^104 + 6*u^102 + 4*u^
#                100 + 2*u^98 + u^96 - 2*u^94 - 2*u^92 - 4*u^90 - 4*u^88 - 4*u^
#                86 - 4*u^84 - 3*u^82 - 2*u^80 - u^78 + u^74 + u^72 + 2*u^
#                70 + u^68 + u^66, u^191 + u^189 + u^187 + u^185 - u^179 - 2*u^
#                177 - 2*u^175 - 3*u^173 - 3*u^171 - 2*u^169 - 3*u^167 - u^
#                165 + 3*u^159 + 3*u^157 + 4*u^155 + 6*u^153 + 4*u^151 + 5*u^
#                149 + 5*u^147 + 2*u^145 + 3*u^143 - 2*u^139 - u^137 - 5*u^
#                135 - 5*u^133 - 5*u^131 - 8*u^129 - 5*u^127 - 5*u^125 - 5*u^
#                123 - u^121 - 2*u^119 + 3*u^115 + 2*u^113 + 5*u^111 + 5*u^
#                109 + 4*u^107 + 6*u^105 + 4*u^103 + 3*u^101 + 3*u^99 - u^93 - 
#                3*u^91 - 2*u^89 - 3*u^87 - 3*u^85 - 2*u^83 - 2*u^81 - u^
#                79 + u^73 + u^71 + u^69 + u^67, 
#              u^190 + u^186 - u^178 - u^176 - u^174 - 2*u^172 - u^170 - u^
#                168 - 2*u^166 + u^164 - u^162 + u^160 + 2*u^158 + u^156 + 3*u^
#                154 + 3*u^152 + u^150 + 4*u^148 + u^146 + u^144 + 2*u^142 - 
#                2*u^140 - u^136 - 4*u^134 - u^132 - 4*u^130 - 4*u^128 - u^
#                126 - 4*u^124 - u^122 - 2*u^118 + 2*u^116 + u^114 + u^112 + 
#                4*u^110 + u^108 + 3*u^106 + 3*u^104 + u^102 + 2*u^100 + u^
#                98 - u^96 + u^94 - 2*u^92 - u^90 - u^88 - 2*u^86 - u^84 - u^
#                82 - u^80 + u^72 + u^68 ], 
#          [ u^191 + u^189 + u^187 + u^185 - u^179 - 2*u^177 - 2*u^175 - 3*u^
#                173 - 3*u^171 - 2*u^169 - 3*u^167 - u^165 + 3*u^159 + 3*u^
#                157 + 4*u^155 + 6*u^153 + 4*u^151 + 5*u^149 + 5*u^147 + 2*u^
#                145 + 3*u^143 - 2*u^139 - u^137 - 5*u^135 - 5*u^133 - 5*u^
#                131 - 8*u^129 - 5*u^127 - 5*u^125 - 5*u^123 - u^121 - 2*u^
#                119 + 3*u^115 + 2*u^113 + 5*u^111 + 5*u^109 + 4*u^107 + 6*u^
#                105 + 4*u^103 + 3*u^101 + 3*u^99 - u^93 - 3*u^91 - 2*u^89 - 
#                3*u^87 - 3*u^85 - 2*u^83 - 2*u^81 - u^79 + u^73 + u^71 + u^
#                69 + u^67, u^192 + u^190 + u^188 + u^186 - u^180 - 2*u^178 - 
#                2*u^176 - 3*u^174 - 3*u^172 - 2*u^170 - 3*u^168 - u^166 + 3*u^
#                160 + 3*u^158 + 4*u^156 + 6*u^154 + 4*u^152 + 5*u^150 + 5*u^
#                148 + 2*u^146 + 3*u^144 - 2*u^140 - u^138 - 5*u^136 - 5*u^
#                134 - 5*u^132 - 8*u^130 - 5*u^128 - 5*u^126 - 5*u^124 - u^
#                122 - 2*u^120 + 3*u^116 + 2*u^114 + 5*u^112 + 5*u^110 + 4*u^
#                108 + 6*u^106 + 4*u^104 + 3*u^102 + 3*u^100 - u^94 - 3*u^92 - 
#                2*u^90 - 3*u^88 - 3*u^86 - 2*u^84 - 2*u^82 - u^80 + u^74 + u^
#                72 + u^70 + u^68, 0*u^0 ], 
#          [ u^190 + u^186 - u^178 - u^176 - u^174 - 2*u^172 - u^170 - u^168 - 
#                2*u^166 + u^164 - u^162 + u^160 + 2*u^158 + u^156 + 3*u^154 + 
#                3*u^152 + u^150 + 4*u^148 + u^146 + u^144 + 2*u^142 - 2*u^
#                140 - u^136 - 4*u^134 - u^132 - 4*u^130 - 4*u^128 - u^126 - 
#                4*u^124 - u^122 - 2*u^118 + 2*u^116 + u^114 + u^112 + 4*u^
#                110 + u^108 + 3*u^106 + 3*u^104 + u^102 + 2*u^100 + u^98 - u^
#                96 + u^94 - 2*u^92 - u^90 - u^88 - 2*u^86 - u^84 - u^82 - u^
#                80 + u^72 + u^68, 0*u^0, u^192 + u^188 - u^180 - u^178 - u^
#                176 - 2*u^174 - u^172 - u^170 - 2*u^168 + u^166 - u^164 + u^
#                162 + 2*u^160 + u^158 + 3*u^156 + 3*u^154 + u^152 + 4*u^
#                150 + u^148 + u^146 + 2*u^144 - 2*u^142 - u^138 - 4*u^136 - u^
#                134 - 4*u^132 - 4*u^130 - u^128 - 4*u^126 - u^124 - 2*u^120 + 
#                2*u^118 + u^116 + u^114 + 4*u^112 + u^110 + 3*u^108 + 3*u^
#                106 + u^104 + 2*u^102 + u^100 - u^98 + u^96 - 2*u^94 - u^
#                92 - u^90 - 2*u^88 - u^86 - u^84 - u^82 + u^74 + u^70 ] ], 
#      [ [ u^194 + u^192 + u^190 + u^188 - u^182 - 2*u^180 - 2*u^178 - 3*u^
#                176 - 3*u^174 - 2*u^172 - 3*u^170 - u^168 + 3*u^162 + 3*u^
#                160 + 4*u^158 + 6*u^156 + 4*u^154 + 5*u^152 + 5*u^150 + 2*u^
#                148 + 3*u^146 - 2*u^142 - u^140 - 5*u^138 - 5*u^136 - 5*u^
#                134 - 8*u^132 - 5*u^130 - 5*u^128 - 5*u^126 - u^124 - 2*u^
#                122 + 3*u^118 + 2*u^116 + 5*u^114 + 5*u^112 + 4*u^110 + 6*u^
#                108 + 4*u^106 + 3*u^104 + 3*u^102 - u^96 - 3*u^94 - 2*u^92 - 
#                3*u^90 - 3*u^88 - 2*u^86 - 2*u^84 - u^82 + u^76 + u^74 + u^
#                72 + u^70 ] ], 
#      [ [ u^196 + u^194 + u^192 + u^190 - u^184 - 2*u^182 - 2*u^180 - 3*u^
#                178 - 3*u^176 - 2*u^174 - 3*u^172 - u^170 + 3*u^164 + 3*u^
#                162 + 4*u^160 + 6*u^158 + 4*u^156 + 5*u^154 + 5*u^152 + 2*u^
#                150 + 3*u^148 - 2*u^144 - u^142 - 5*u^140 - 5*u^138 - 5*u^
#                136 - 8*u^134 - 5*u^132 - 5*u^130 - 5*u^128 - u^126 - 2*u^
#                124 + 3*u^120 + 2*u^118 + 5*u^116 + 5*u^114 + 4*u^112 + 6*u^
#                110 + 4*u^108 + 3*u^106 + 3*u^104 - u^98 - 3*u^96 - 2*u^94 - 
#                3*u^92 - 3*u^90 - 2*u^88 - 2*u^86 - u^84 + u^78 + u^76 + u^
#                74 + u^72, 0*u^0 ], 
#          [ 0*u^0, u^196 - u^188 - u^184 - u^182 - u^178 + u^174 - u^172 + 
#                2*u^170 + u^168 + 3*u^164 + u^160 + 2*u^158 - 2*u^156 + u^
#                154 - 3*u^150 + u^148 - 3*u^146 - 2*u^144 + u^142 - 4*u^140 - 
#                3*u^134 + 3*u^132 + 4*u^126 - u^124 + 2*u^122 + 3*u^120 - u^
#                118 + 3*u^116 - u^112 + 2*u^110 - 2*u^108 - u^106 - 3*u^
#                102 - u^98 - 2*u^96 + u^94 - u^92 + u^88 + u^84 + u^82 + u^
#                78 - u^70 ] ], 
#      [ [ u^198 + u^196 + u^194 + 2*u^192 + u^190 + u^188 + u^186 - u^184 - u^
#                182 - 2*u^180 - 4*u^178 - 3*u^176 - 5*u^174 - 5*u^172 - 3*u^
#                170 - 5*u^168 - 2*u^166 - u^162 + 4*u^160 + 4*u^158 + 4*u^
#                156 + 9*u^154 + 6*u^152 + 7*u^150 + 9*u^148 + 4*u^146 + 6*u^
#                144 + 4*u^142 - u^140 + u^138 - 4*u^136 - 6*u^134 - 4*u^132 - 
#                9*u^130 - 7*u^128 - 6*u^126 - 9*u^124 - 4*u^122 - 4*u^120 - 
#                4*u^118 + u^116 + 2*u^112 + 5*u^110 + 3*u^108 + 5*u^106 + 5*u^
#                104 + 3*u^102 + 4*u^100 + 2*u^98 + u^96 + u^94 - u^92 - u^
#                90 - u^88 - 2*u^86 - u^84 - u^82 - u^80, 0*u^0, 
#              u^197 + u^195 + u^193 + 2*u^191 + u^189 + u^187 + u^185 - u^
#                183 - u^181 - 2*u^179 - 4*u^177 - 3*u^175 - 5*u^173 - 5*u^
#                171 - 3*u^169 - 5*u^167 - 2*u^165 - u^161 + 4*u^159 + 4*u^
#                157 + 4*u^155 + 9*u^153 + 6*u^151 + 7*u^149 + 9*u^147 + 4*u^
#                145 + 6*u^143 + 4*u^141 - u^139 + u^137 - 4*u^135 - 6*u^133 - 
#                4*u^131 - 9*u^129 - 7*u^127 - 6*u^125 - 9*u^123 - 4*u^121 - 
#                4*u^119 - 4*u^117 + u^115 + 2*u^111 + 5*u^109 + 3*u^107 + 5*u^
#                105 + 5*u^103 + 3*u^101 + 4*u^99 + 2*u^97 + u^95 + u^93 - u^
#                91 - u^89 - u^87 - 2*u^85 - u^83 - u^81 - u^79, 0*u^0 ], 
#          [ 0*u^0, u^198 + u^192 - u^190 - 2*u^184 - u^180 - 2*u^178 + u^
#                176 - 2*u^174 + 2*u^170 - 2*u^168 + 3*u^166 + 2*u^164 - u^
#                162 + 5*u^160 + 5*u^154 - 3*u^152 + u^150 + 2*u^148 - 5*u^
#                146 + 2*u^144 - 2*u^142 - 5*u^140 + 2*u^138 - 5*u^136 - 2*u^
#                134 + 2*u^132 - 5*u^130 + 2*u^128 + u^126 - 3*u^124 + 5*u^
#                122 + 5*u^116 - u^114 + 2*u^112 + 3*u^110 - 2*u^108 + 2*u^
#                106 - 2*u^102 + u^100 - 2*u^98 - u^96 - 2*u^92 - u^86 + u^
#                84 + u^78, 0*u^0, 0*u^0 ], 
#          [ u^197 + u^195 + u^193 + 2*u^191 + u^189 + u^187 + u^185 - u^
#                183 - u^181 - 2*u^179 - 4*u^177 - 3*u^175 - 5*u^173 - 5*u^
#                171 - 3*u^169 - 5*u^167 - 2*u^165 - u^161 + 4*u^159 + 4*u^
#                157 + 4*u^155 + 9*u^153 + 6*u^151 + 7*u^149 + 9*u^147 + 4*u^
#                145 + 6*u^143 + 4*u^141 - u^139 + u^137 - 4*u^135 - 6*u^133 - 
#                4*u^131 - 9*u^129 - 7*u^127 - 6*u^125 - 9*u^123 - 4*u^121 - 
#                4*u^119 - 4*u^117 + u^115 + 2*u^111 + 5*u^109 + 3*u^107 + 5*u^
#                105 + 5*u^103 + 3*u^101 + 4*u^99 + 2*u^97 + u^95 + u^93 - u^
#                91 - u^89 - u^87 - 2*u^85 - u^83 - u^81 - u^79, 0*u^0, 
#              u^198 + u^196 + u^194 + 2*u^192 + u^190 + u^188 + u^186 - u^
#                184 - u^182 - 2*u^180 - 4*u^178 - 3*u^176 - 5*u^174 - 5*u^
#                172 - 3*u^170 - 5*u^168 - 2*u^166 - u^162 + 4*u^160 + 4*u^
#                158 + 4*u^156 + 9*u^154 + 6*u^152 + 7*u^150 + 9*u^148 + 4*u^
#                146 + 6*u^144 + 4*u^142 - u^140 + u^138 - 4*u^136 - 6*u^134 - 
#                4*u^132 - 9*u^130 - 7*u^128 - 6*u^126 - 9*u^124 - 4*u^122 - 
#                4*u^120 - 4*u^118 + u^116 + 2*u^112 + 5*u^110 + 3*u^108 + 5*u^
#                106 + 5*u^104 + 3*u^102 + 4*u^100 + 2*u^98 + u^96 + u^94 - u^
#                92 - u^90 - u^88 - 2*u^86 - u^84 - u^82 - u^80, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, u^198 + u^192 - u^190 - 2*u^184 - u^180 - 
#                2*u^178 + u^176 - 2*u^174 + 2*u^170 - 2*u^168 + 3*u^166 + 2*u^
#                164 - u^162 + 5*u^160 + 5*u^154 - 3*u^152 + u^150 + 2*u^148 - 
#                5*u^146 + 2*u^144 - 2*u^142 - 5*u^140 + 2*u^138 - 5*u^136 - 
#                2*u^134 + 2*u^132 - 5*u^130 + 2*u^128 + u^126 - 3*u^124 + 5*u^
#                122 + 5*u^116 - u^114 + 2*u^112 + 3*u^110 - 2*u^108 + 2*u^
#                106 - 2*u^102 + u^100 - 2*u^98 - u^96 - 2*u^92 - u^86 + u^
#                84 + u^78 ] ], 
#      [ [ u^200 + u^196 + u^194 + u^190 - u^186 - 2*u^182 - 2*u^180 - u^178 - 
#                4*u^176 - u^174 - 2*u^172 - 3*u^170 + u^168 - u^166 + 4*u^
#                162 + 4*u^158 + 5*u^156 + u^154 + 6*u^152 + 3*u^150 + u^148 + 
#                5*u^146 - u^144 + u^140 - 5*u^138 - u^136 - 3*u^134 - 6*u^
#                132 - u^130 - 5*u^128 - 4*u^126 - 4*u^122 + u^118 - u^116 + 
#                3*u^114 + 2*u^112 + u^110 + 4*u^108 + u^106 + 2*u^104 + 2*u^
#                102 + u^98 - u^94 - u^90 - u^88 - u^84 ] ], 
#      [ [ u^208 + 2*u^204 + 2*u^202 + 2*u^200 + u^198 - u^196 - 2*u^194 - 3*u^
#                192 - 5*u^190 - 5*u^188 - 5*u^186 - 6*u^184 - 2*u^182 - 3*u^
#                180 + 3*u^176 + 4*u^174 + 7*u^172 + 9*u^170 + 8*u^168 + 11*u^
#                166 + 8*u^164 + 6*u^162 + 6*u^160 + u^158 - 3*u^154 - 8*u^
#                152 - 7*u^150 - 10*u^148 - 12*u^146 - 10*u^144 - 12*u^142 - 
#                9*u^140 - 5*u^138 - 4*u^136 + u^134 + 3*u^132 + 4*u^130 + 8*u^
#                128 + 8*u^126 + 10*u^124 + 10*u^122 + 8*u^120 + 7*u^118 + 6*u^
#                116 + 2*u^114 + u^112 - 3*u^110 - 4*u^108 - 5*u^106 - 6*u^
#                104 - 6*u^102 - 5*u^100 - 4*u^98 - 3*u^96 - u^94 + 2*u^90 + 
#                2*u^88 + 3*u^86 + 2*u^84 + 2*u^82 + u^80 - u^76 - u^74, 
#              u^207 + 2*u^205 + 4*u^203 + 4*u^201 + 2*u^199 + u^197 - 2*u^
#                195 - 4*u^193 - 6*u^191 - 9*u^189 - 9*u^187 - 9*u^185 - 9*u^
#                183 - 5*u^181 - 4*u^179 + 6*u^175 + 7*u^173 + 13*u^171 + 16*u^
#                169 + 15*u^167 + 18*u^165 + 14*u^163 + 11*u^161 + 11*u^159 + 
#                2*u^157 - u^155 - 5*u^153 - 13*u^151 - 12*u^149 - 18*u^147 - 
#                22*u^145 - 18*u^143 - 21*u^141 - 15*u^139 - 8*u^137 - 8*u^
#                135 + u^133 + 4*u^131 + 6*u^129 + 15*u^127 + 14*u^125 + 16*u^
#                123 + 18*u^121 + 13*u^119 + 14*u^117 + 10*u^115 + 3*u^113 + 
#                2*u^111 - 4*u^109 - 7*u^107 - 7*u^105 - 11*u^103 - 9*u^101 - 
#                8*u^99 - 8*u^97 - 4*u^95 - 2*u^93 + 3*u^89 + 3*u^87 + 4*u^
#                85 + 4*u^83 + 2*u^81 + u^79 - u^77 - u^75, 
#              u^206 + 2*u^204 + 2*u^202 + u^200 - u^196 - 2*u^194 - 3*u^192 - 
#                4*u^190 - 4*u^188 - 4*u^186 - 3*u^184 - 2*u^182 - u^180 + u^
#                178 + 3*u^176 + 4*u^174 + 6*u^172 + 7*u^170 + 7*u^168 + 7*u^
#                166 + 5*u^164 + 4*u^162 + 3*u^160 - 2*u^156 - 4*u^154 - 6*u^
#                152 - 6*u^150 - 8*u^148 - 9*u^146 - 8*u^144 - 8*u^142 - 5*u^
#                140 - 2*u^138 - u^136 + 2*u^134 + 3*u^132 + 4*u^130 + 7*u^
#                128 + 7*u^126 + 7*u^124 + 7*u^122 + 5*u^120 + 5*u^118 + 3*u^
#                116 - u^112 - 3*u^110 - 4*u^108 - 4*u^106 - 5*u^104 - 4*u^
#                102 - 3*u^100 - 3*u^98 - u^96 + u^92 + 2*u^90 + 2*u^88 + 2*u^
#                86 + 2*u^84 + u^82 - u^78 - u^76, 
#              u^206 + 2*u^204 + 3*u^202 + u^200 + u^198 - u^196 - 2*u^194 - 
#                3*u^192 - 5*u^190 - 5*u^188 - 5*u^186 - 5*u^184 - 3*u^182 - 
#                2*u^180 - u^178 + 4*u^176 + 3*u^174 + 7*u^172 + 9*u^170 + 8*u^
#                168 + 10*u^166 + 8*u^164 + 5*u^162 + 7*u^160 + u^158 - u^
#                156 - 2*u^154 - 8*u^152 - 6*u^150 - 9*u^148 - 13*u^146 - 9*u^
#                144 - 12*u^142 - 9*u^140 - 3*u^138 - 5*u^136 + u^134 + 3*u^
#                132 + 2*u^130 + 9*u^128 + 8*u^126 + 8*u^124 + 11*u^122 + 6*u^
#                120 + 8*u^118 + 6*u^116 + u^114 + u^112 - 2*u^110 - 5*u^108 - 
#                3*u^106 - 7*u^104 - 5*u^102 - 4*u^100 - 5*u^98 - 2*u^96 - u^
#                94 + 2*u^90 + 2*u^88 + 2*u^86 + 3*u^84 + u^82 + u^80 - u^
#                78 - u^76, u^206 + 2*u^204 + 3*u^202 + 2*u^200 + u^198 - 2*u^
#                194 - 3*u^192 - 5*u^190 - 6*u^188 - 6*u^186 - 6*u^184 - 5*u^
#                182 - 3*u^180 - 2*u^178 + 2*u^176 + 4*u^174 + 6*u^172 + 10*u^
#                170 + 10*u^168 + 11*u^166 + 11*u^164 + 8*u^162 + 8*u^160 + 
#                5*u^158 - u^154 - 6*u^152 - 8*u^150 - 9*u^148 - 14*u^146 - 
#                13*u^144 - 13*u^142 - 13*u^140 - 7*u^138 - 6*u^136 - 3*u^
#                134 + 2*u^132 + 2*u^130 + 7*u^128 + 10*u^126 + 9*u^124 + 12*u^
#                122 + 10*u^120 + 9*u^118 + 9*u^116 + 4*u^114 + 2*u^112 - 4*u^
#                108 - 4*u^106 - 6*u^104 - 7*u^102 - 5*u^100 - 6*u^98 - 4*u^
#                96 - 2*u^94 - u^92 + u^90 + 2*u^88 + 2*u^86 + 3*u^84 + 2*u^
#                82 + u^80 - u^76, u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^
#                198 - u^194 - 3*u^192 - 4*u^190 - 5*u^188 - 6*u^186 - 5*u^
#                184 - 5*u^182 - 4*u^180 - u^178 + 3*u^174 + 6*u^172 + 7*u^
#                170 + 10*u^168 + 10*u^166 + 9*u^164 + 10*u^162 + 7*u^160 + 
#                5*u^158 + 3*u^156 - 2*u^154 - 3*u^152 - 6*u^150 - 10*u^148 - 
#                10*u^146 - 13*u^144 - 13*u^142 - 10*u^140 - 10*u^138 - 6*u^
#                136 - 3*u^134 - 2*u^132 + 3*u^130 + 5*u^128 + 7*u^126 + 10*u^
#                124 + 9*u^122 + 10*u^120 + 10*u^118 + 7*u^116 + 6*u^114 + 3*u^
#                112 - u^108 - 4*u^106 - 5*u^104 - 5*u^102 - 6*u^100 - 5*u^
#                98 - 4*u^96 - 3*u^94 - u^92 + u^88 + 2*u^86 + 2*u^84 + 2*u^
#                82 + u^80, 2*u^205 + 3*u^203 + 2*u^201 + 2*u^199 - u^195 - 
#                2*u^193 - 5*u^191 - 5*u^189 - 6*u^187 - 7*u^185 - 4*u^183 - 
#                5*u^181 - 3*u^179 + 2*u^177 + u^175 + 6*u^173 + 9*u^171 + 8*u^
#                169 + 13*u^167 + 10*u^165 + 8*u^163 + 11*u^161 + 4*u^159 + 
#                3*u^157 + u^155 - 7*u^153 - 4*u^151 - 9*u^149 - 14*u^147 - 
#                10*u^145 - 16*u^143 - 13*u^141 - 7*u^139 - 10*u^137 - 2*u^
#                135 - u^131 + 8*u^129 + 7*u^127 + 9*u^125 + 13*u^123 + 8*u^
#                121 + 11*u^119 + 10*u^117 + 4*u^115 + 5*u^113 - 3*u^109 - 2*u^
#                107 - 7*u^105 - 6*u^103 - 5*u^101 - 7*u^99 - 4*u^97 - 3*u^
#                95 - 2*u^93 + u^91 + u^89 + 2*u^87 + 3*u^85 + 2*u^83 + 2*u^
#                81 - u^77, u^205 + 2*u^203 + u^201 + u^199 - u^195 - u^193 - 
#                3*u^191 - 3*u^189 - 3*u^187 - 4*u^185 - 2*u^183 - 2*u^181 - 
#                2*u^179 + 2*u^177 + u^175 + 3*u^173 + 6*u^171 + 4*u^169 + 7*u^
#                167 + 6*u^165 + 3*u^163 + 6*u^161 + 2*u^159 + u^155 - 5*u^
#                153 - 3*u^151 - 4*u^149 - 9*u^147 - 5*u^145 - 8*u^143 - 8*u^
#                141 - 2*u^139 - 5*u^137 - u^135 + 2*u^133 - u^131 + 5*u^129 + 
#                5*u^127 + 4*u^125 + 8*u^123 + 4*u^121 + 5*u^119 + 6*u^117 + u^
#                115 + 2*u^113 - 3*u^109 - u^107 - 4*u^105 - 4*u^103 - 2*u^
#                101 - 4*u^99 - 2*u^97 - u^95 - u^93 + u^91 + u^89 + u^87 + 
#                2*u^85 + u^83 + u^81 - u^77, 
#              u^205 + u^203 - u^197 - u^195 - u^193 - 2*u^191 - u^189 - u^
#                187 - u^185 + u^183 + u^179 + 3*u^177 + u^175 + 3*u^173 + 3*u^
#                171 + u^169 + 3*u^167 - u^163 + u^161 - 3*u^159 - 2*u^157 - 
#                2*u^155 - 5*u^153 - u^151 - 3*u^149 - 4*u^147 - 3*u^143 + 3*u^
#                139 + 4*u^135 + 3*u^133 + u^131 + 5*u^129 + 2*u^127 + 2*u^
#                125 + 3*u^123 - u^121 + u^119 - 3*u^115 - u^113 - 3*u^111 - 
#                3*u^109 - u^107 - 3*u^105 - u^103 - u^99 + u^97 + u^95 + u^
#                93 + 2*u^91 + u^89 + u^87 + u^85 - u^79 - u^77, 
#              u^204 - u^196 - u^192 - u^190 - u^186 + u^182 - u^180 + 2*u^
#                178 + u^176 + 3*u^172 + u^168 + 2*u^166 - 2*u^164 + u^162 - 
#                3*u^158 + u^156 - 3*u^154 - 2*u^152 + u^150 - 4*u^148 - 3*u^
#                142 + 3*u^140 + 4*u^134 - u^132 + 2*u^130 + 3*u^128 - u^126 + 
#                3*u^124 - u^120 + 2*u^118 - 2*u^116 - u^114 - 3*u^110 - u^
#                106 - 2*u^104 + u^102 - u^100 + u^96 + u^92 + u^90 + u^86 - u^
#                78, u^204 - u^196 - u^192 - u^190 - u^186 + u^182 - u^180 + 
#                2*u^178 + u^176 + 3*u^172 + u^168 + 2*u^166 - 2*u^164 + u^
#                162 - 3*u^158 + u^156 - 3*u^154 - 2*u^152 + u^150 - 4*u^148 - 
#                3*u^142 + 3*u^140 + 4*u^134 - u^132 + 2*u^130 + 3*u^128 - u^
#                126 + 3*u^124 - u^120 + 2*u^118 - 2*u^116 - u^114 - 3*u^
#                110 - u^106 - 2*u^104 + u^102 - u^100 + u^96 + u^92 + u^
#                90 + u^86 - u^78, 2*u^204 + u^202 + 2*u^200 + u^198 - 2*u^
#                192 - 3*u^190 - 3*u^188 - 5*u^186 - 4*u^184 - 3*u^182 - 5*u^
#                180 - u^176 + u^174 + 5*u^172 + 4*u^170 + 7*u^168 + 9*u^166 + 
#                5*u^164 + 9*u^162 + 6*u^160 + 3*u^158 + 5*u^156 - 2*u^154 - 
#                2*u^152 - 2*u^150 - 9*u^148 - 6*u^146 - 9*u^144 - 12*u^142 - 
#                6*u^140 - 9*u^138 - 6*u^136 - u^134 - 4*u^132 + 2*u^130 + 4*u^
#                128 + 3*u^126 + 9*u^124 + 6*u^122 + 7*u^120 + 9*u^118 + 5*u^
#                116 + 5*u^114 + 4*u^112 - u^110 + u^108 - 3*u^106 - 4*u^104 - 
#                3*u^102 - 5*u^100 - 4*u^98 - 3*u^96 - 3*u^94 - u^92 + 2*u^
#                86 + u^84 + 2*u^82 + u^80, u^204 + u^200 - u^192 - u^190 - u^
#                188 - 2*u^186 - u^184 - u^182 - 2*u^180 + u^178 - u^176 + u^
#                174 + 2*u^172 + u^170 + 3*u^168 + 3*u^166 + u^164 + 4*u^
#                162 + u^160 + u^158 + 2*u^156 - 2*u^154 - u^150 - 4*u^148 - u^
#                146 - 4*u^144 - 4*u^142 - u^140 - 4*u^138 - u^136 - 2*u^132 + 
#                2*u^130 + u^128 + u^126 + 4*u^124 + u^122 + 3*u^120 + 3*u^
#                118 + u^116 + 2*u^114 + u^112 - u^110 + u^108 - 2*u^106 - u^
#                104 - u^102 - 2*u^100 - u^98 - u^96 - u^94 + u^86 + u^82, 
#              u^202 + u^198 - u^190 - u^188 - u^186 - 2*u^184 - u^182 - u^
#                180 - 2*u^178 + u^176 - u^174 + u^172 + 2*u^170 + u^168 + 3*u^
#                166 + 3*u^164 + u^162 + 4*u^160 + u^158 + u^156 + 2*u^154 - 
#                2*u^152 - u^148 - 4*u^146 - u^144 - 4*u^142 - 4*u^140 - u^
#                138 - 4*u^136 - u^134 - 2*u^130 + 2*u^128 + u^126 + u^124 + 
#                4*u^122 + u^120 + 3*u^118 + 3*u^116 + u^114 + 2*u^112 + u^
#                110 - u^108 + u^106 - 2*u^104 - u^102 - u^100 - 2*u^98 - u^
#                96 - u^94 - u^92 + u^84 + u^80, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ u^207 + 2*u^205 + 4*u^203 + 4*u^201 + 2*u^199 + u^197 - 2*u^195 - 
#                4*u^193 - 6*u^191 - 9*u^189 - 9*u^187 - 9*u^185 - 9*u^183 - 
#                5*u^181 - 4*u^179 + 6*u^175 + 7*u^173 + 13*u^171 + 16*u^169 + 
#                15*u^167 + 18*u^165 + 14*u^163 + 11*u^161 + 11*u^159 + 2*u^
#                157 - u^155 - 5*u^153 - 13*u^151 - 12*u^149 - 18*u^147 - 22*u^
#                145 - 18*u^143 - 21*u^141 - 15*u^139 - 8*u^137 - 8*u^135 + u^
#                133 + 4*u^131 + 6*u^129 + 15*u^127 + 14*u^125 + 16*u^123 + 
#                18*u^121 + 13*u^119 + 14*u^117 + 10*u^115 + 3*u^113 + 2*u^
#                111 - 4*u^109 - 7*u^107 - 7*u^105 - 11*u^103 - 9*u^101 - 8*u^
#                99 - 8*u^97 - 4*u^95 - 2*u^93 + 3*u^89 + 3*u^87 + 4*u^85 + 
#                4*u^83 + 2*u^81 + u^79 - u^77 - u^75, 
#              u^208 + 3*u^206 + 7*u^204 + 7*u^202 + 5*u^200 + 3*u^198 - 2*u^
#                196 - 5*u^194 - 10*u^192 - 15*u^190 - 16*u^188 - 18*u^186 - 
#                17*u^184 - 12*u^182 - 11*u^180 - 2*u^178 + 6*u^176 + 10*u^
#                174 + 22*u^172 + 26*u^170 + 29*u^168 + 34*u^166 + 27*u^164 + 
#                26*u^162 + 23*u^160 + 9*u^158 + 5*u^156 - 7*u^154 - 18*u^
#                152 - 19*u^150 - 33*u^148 - 37*u^146 - 36*u^144 - 42*u^142 - 
#                30*u^140 - 23*u^138 - 19*u^136 - 3*u^134 + 9*u^130 + 23*u^
#                128 + 23*u^126 + 31*u^124 + 32*u^122 + 27*u^120 + 30*u^118 + 
#                21*u^116 + 12*u^114 + 8*u^112 - 4*u^110 - 8*u^108 - 12*u^
#                106 - 19*u^104 - 16*u^102 - 17*u^100 - 16*u^98 - 10*u^96 - 
#                7*u^94 - 2*u^92 + 3*u^90 + 4*u^88 + 7*u^86 + 7*u^84 + 5*u^
#                82 + 3*u^80 - u^78 - u^76, u^207 + 3*u^205 + 3*u^203 + 2*u^
#                201 + u^199 - u^197 - 2*u^195 - 4*u^193 - 6*u^191 - 6*u^189 - 
#                7*u^187 - 6*u^185 - 4*u^183 - 4*u^181 + 3*u^177 + 4*u^175 + 
#                9*u^173 + 10*u^171 + 11*u^169 + 13*u^167 + 9*u^165 + 9*u^
#                163 + 8*u^161 + 2*u^159 + u^157 - 4*u^155 - 8*u^153 - 7*u^
#                151 - 13*u^149 - 14*u^147 - 13*u^145 - 16*u^143 - 10*u^141 - 
#                7*u^139 - 6*u^137 + u^135 + u^133 + 4*u^131 + 10*u^129 + 9*u^
#                127 + 12*u^125 + 12*u^123 + 9*u^121 + 11*u^119 + 7*u^117 + 
#                3*u^115 + 2*u^113 - 3*u^111 - 4*u^109 - 5*u^107 - 8*u^105 - 
#                6*u^103 - 6*u^101 - 6*u^99 - 3*u^97 - 2*u^95 + 2*u^91 + 2*u^
#                89 + 3*u^87 + 3*u^85 + 2*u^83 + u^81 - u^79 - u^77, 
#              u^207 + 4*u^205 + 4*u^203 + 3*u^201 + 2*u^199 - u^197 - 2*u^
#                195 - 5*u^193 - 8*u^191 - 8*u^189 - 10*u^187 - 9*u^185 - 6*u^
#                183 - 7*u^181 - u^179 + 3*u^177 + 4*u^175 + 12*u^173 + 13*u^
#                171 + 15*u^169 + 19*u^167 + 13*u^165 + 14*u^163 + 13*u^161 + 
#                4*u^159 + 4*u^157 - 4*u^155 - 10*u^153 - 8*u^151 - 18*u^149 - 
#                19*u^147 - 18*u^145 - 24*u^143 - 15*u^141 - 12*u^139 - 11*u^
#                137 - u^133 + 4*u^131 + 13*u^129 + 11*u^127 + 17*u^125 + 17*u^
#                123 + 13*u^121 + 17*u^119 + 11*u^117 + 6*u^115 + 5*u^113 - 
#                3*u^111 - 4*u^109 - 6*u^107 - 11*u^105 - 8*u^103 - 9*u^101 - 
#                9*u^99 - 5*u^97 - 4*u^95 - u^93 + 2*u^91 + 2*u^89 + 4*u^87 + 
#                4*u^85 + 3*u^83 + 2*u^81 - u^79 - u^77, 
#              u^207 + 3*u^205 + 5*u^203 + 4*u^201 + 3*u^199 + u^197 - 2*u^
#                195 - 4*u^193 - 8*u^191 - 10*u^189 - 11*u^187 - 12*u^185 - 
#                10*u^183 - 8*u^181 - 6*u^179 + u^177 + 4*u^175 + 9*u^173 + 
#                16*u^171 + 17*u^169 + 21*u^167 + 21*u^165 + 17*u^163 + 18*u^
#                161 + 12*u^159 + 5*u^157 + 2*u^155 - 8*u^153 - 11*u^151 - 
#                15*u^149 - 24*u^147 - 23*u^145 - 26*u^143 - 26*u^141 - 17*u^
#                139 - 16*u^137 - 9*u^135 - u^133 + 10*u^129 + 15*u^127 + 16*u^
#                125 + 22*u^123 + 19*u^121 + 19*u^119 + 19*u^117 + 11*u^115 + 
#                8*u^113 + 3*u^111 - 4*u^109 - 5*u^107 - 10*u^105 - 12*u^103 - 
#                10*u^101 - 12*u^99 - 9*u^97 - 6*u^95 - 4*u^93 + 2*u^89 + 3*u^
#                87 + 5*u^85 + 4*u^83 + 3*u^81 + u^79 - u^77, 
#              u^207 + 3*u^205 + 4*u^203 + 4*u^201 + 3*u^199 + u^197 - u^195 - 
#                4*u^193 - 7*u^191 - 9*u^189 - 11*u^187 - 11*u^185 - 10*u^
#                183 - 9*u^181 - 5*u^179 - u^177 + 3*u^175 + 9*u^173 + 13*u^
#                171 + 17*u^169 + 20*u^167 + 19*u^165 + 19*u^163 + 17*u^161 + 
#                12*u^159 + 8*u^157 + u^155 - 5*u^153 - 9*u^151 - 16*u^149 - 
#                20*u^147 - 23*u^145 - 26*u^143 - 23*u^141 - 20*u^139 - 16*u^
#                137 - 9*u^135 - 5*u^133 + u^131 + 8*u^129 + 12*u^127 + 17*u^
#                125 + 19*u^123 + 19*u^121 + 20*u^119 + 17*u^117 + 13*u^115 + 
#                9*u^113 + 3*u^111 - u^109 - 5*u^107 - 9*u^105 - 10*u^103 - 
#                11*u^101 - 11*u^99 - 9*u^97 - 7*u^95 - 4*u^93 - u^91 + u^89 + 
#                3*u^87 + 4*u^85 + 4*u^83 + 3*u^81 + u^79, 
#              3*u^206 + 5*u^204 + 4*u^202 + 4*u^200 + u^198 - u^196 - 3*u^
#                194 - 8*u^192 - 9*u^190 - 11*u^188 - 13*u^186 - 9*u^184 - 
#                10*u^182 - 7*u^180 + u^178 + u^176 + 9*u^174 + 15*u^172 + 
#                15*u^170 + 23*u^168 + 20*u^166 + 17*u^164 + 21*u^162 + 11*u^
#                160 + 8*u^158 + 4*u^156 - 9*u^154 - 7*u^152 - 15*u^150 - 24*u^
#                148 - 20*u^146 - 29*u^144 - 26*u^142 - 17*u^140 - 20*u^138 - 
#                8*u^136 - 3*u^134 - 3*u^132 + 11*u^130 + 12*u^128 + 16*u^
#                126 + 23*u^124 + 17*u^122 + 21*u^120 + 20*u^118 + 11*u^116 + 
#                11*u^114 + 3*u^112 - 3*u^110 - 3*u^108 - 11*u^106 - 11*u^
#                104 - 10*u^102 - 13*u^100 - 9*u^98 - 7*u^96 - 5*u^94 + u^90 + 
#                3*u^88 + 5*u^86 + 4*u^84 + 4*u^82 + u^80 - u^78, 
#              2*u^206 + 3*u^204 + 2*u^202 + 2*u^200 - u^196 - 2*u^194 - 5*u^
#                192 - 5*u^190 - 6*u^188 - 7*u^186 - 4*u^184 - 5*u^182 - 3*u^
#                180 + 2*u^178 + u^176 + 6*u^174 + 9*u^172 + 8*u^170 + 13*u^
#                168 + 10*u^166 + 8*u^164 + 11*u^162 + 4*u^160 + 3*u^158 + u^
#                156 - 7*u^154 - 4*u^152 - 9*u^150 - 14*u^148 - 10*u^146 - 
#                16*u^144 - 13*u^142 - 7*u^140 - 10*u^138 - 2*u^136 - u^132 + 
#                8*u^130 + 7*u^128 + 9*u^126 + 13*u^124 + 8*u^122 + 11*u^120 + 
#                10*u^118 + 4*u^116 + 5*u^114 - 3*u^110 - 2*u^108 - 7*u^106 - 
#                6*u^104 - 5*u^102 - 7*u^100 - 4*u^98 - 3*u^96 - 2*u^94 + u^
#                92 + u^90 + 2*u^88 + 3*u^86 + 2*u^84 + 2*u^82 - u^78, 
#              u^206 + u^204 - u^198 - u^196 - u^194 - 2*u^192 - u^190 - u^
#                188 - u^186 + u^184 + u^180 + 3*u^178 + u^176 + 3*u^174 + 3*u^
#                172 + u^170 + 3*u^168 - u^164 + u^162 - 3*u^160 - 2*u^158 - 
#                2*u^156 - 5*u^154 - u^152 - 3*u^150 - 4*u^148 - 3*u^144 + 3*u^
#                140 + 4*u^136 + 3*u^134 + u^132 + 5*u^130 + 2*u^128 + 2*u^
#                126 + 3*u^124 - u^122 + u^120 - 3*u^116 - u^114 - 3*u^112 - 
#                3*u^110 - u^108 - 3*u^106 - u^104 - u^100 + u^98 + u^96 + u^
#                94 + 2*u^92 + u^90 + u^88 + u^86 - u^80 - u^78, 
#              u^205 - u^197 - u^193 - u^191 - u^187 + u^183 - u^181 + 2*u^
#                179 + u^177 + 3*u^173 + u^169 + 2*u^167 - 2*u^165 + u^163 - 
#                3*u^159 + u^157 - 3*u^155 - 2*u^153 + u^151 - 4*u^149 - 3*u^
#                143 + 3*u^141 + 4*u^135 - u^133 + 2*u^131 + 3*u^129 - u^127 + 
#                3*u^125 - u^121 + 2*u^119 - 2*u^117 - u^115 - 3*u^111 - u^
#                107 - 2*u^105 + u^103 - u^101 + u^97 + u^93 + u^91 + u^87 - u^
#                79, u^205 - u^197 - u^193 - u^191 - u^187 + u^183 - u^181 + 
#                2*u^179 + u^177 + 3*u^173 + u^169 + 2*u^167 - 2*u^165 + u^
#                163 - 3*u^159 + u^157 - 3*u^155 - 2*u^153 + u^151 - 4*u^149 - 
#                3*u^143 + 3*u^141 + 4*u^135 - u^133 + 2*u^131 + 3*u^129 - u^
#                127 + 3*u^125 - u^121 + 2*u^119 - 2*u^117 - u^115 - 3*u^
#                111 - u^107 - 2*u^105 + u^103 - u^101 + u^97 + u^93 + u^
#                91 + u^87 - u^79, u^207 + 3*u^205 + 3*u^203 + 3*u^201 + 2*u^
#                199 - u^195 - 4*u^193 - 6*u^191 - 7*u^189 - 9*u^187 - 8*u^
#                185 - 7*u^183 - 7*u^181 - 2*u^179 + 3*u^175 + 9*u^173 + 10*u^
#                171 + 14*u^169 + 16*u^167 + 13*u^165 + 15*u^163 + 12*u^161 + 
#                7*u^159 + 6*u^157 - 2*u^155 - 5*u^153 - 7*u^151 - 15*u^149 - 
#                15*u^147 - 18*u^145 - 21*u^143 - 15*u^141 - 15*u^139 - 11*u^
#                137 - 4*u^135 - 4*u^133 + 3*u^131 + 8*u^129 + 9*u^127 + 15*u^
#                125 + 14*u^123 + 14*u^121 + 16*u^119 + 11*u^117 + 9*u^115 + 
#                6*u^113 - u^109 - 5*u^107 - 8*u^105 - 7*u^103 - 9*u^101 - 8*u^
#                99 - 6*u^97 - 5*u^95 - 2*u^93 + u^89 + 3*u^87 + 3*u^85 + 3*u^
#                83 + 2*u^81, u^205 + u^203 + u^201 + u^199 - u^193 - 2*u^
#                191 - 2*u^189 - 3*u^187 - 3*u^185 - 2*u^183 - 3*u^181 - u^
#                179 + 3*u^173 + 3*u^171 + 4*u^169 + 6*u^167 + 4*u^165 + 5*u^
#                163 + 5*u^161 + 2*u^159 + 3*u^157 - 2*u^153 - u^151 - 5*u^
#                149 - 5*u^147 - 5*u^145 - 8*u^143 - 5*u^141 - 5*u^139 - 5*u^
#                137 - u^135 - 2*u^133 + 3*u^129 + 2*u^127 + 5*u^125 + 5*u^
#                123 + 4*u^121 + 6*u^119 + 4*u^117 + 3*u^115 + 3*u^113 - u^
#                107 - 3*u^105 - 2*u^103 - 3*u^101 - 3*u^99 - 2*u^97 - 2*u^
#                95 - u^93 + u^87 + u^85 + u^83 + u^81, 
#              u^205 + u^203 + u^201 + u^199 - u^193 - 2*u^191 - 2*u^189 - 3*u^
#                187 - 3*u^185 - 2*u^183 - 3*u^181 - u^179 + 3*u^173 + 3*u^
#                171 + 4*u^169 + 6*u^167 + 4*u^165 + 5*u^163 + 5*u^161 + 2*u^
#                159 + 3*u^157 - 2*u^153 - u^151 - 5*u^149 - 5*u^147 - 5*u^
#                145 - 8*u^143 - 5*u^141 - 5*u^139 - 5*u^137 - u^135 - 2*u^
#                133 + 3*u^129 + 2*u^127 + 5*u^125 + 5*u^123 + 4*u^121 + 6*u^
#                119 + 4*u^117 + 3*u^115 + 3*u^113 - u^107 - 3*u^105 - 2*u^
#                103 - 3*u^101 - 3*u^99 - 2*u^97 - 2*u^95 - u^93 + u^87 + u^
#                85 + u^83 + u^81, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ u^206 + 2*u^204 + 2*u^202 + u^200 - u^196 - 2*u^194 - 3*u^192 - 
#                4*u^190 - 4*u^188 - 4*u^186 - 3*u^184 - 2*u^182 - u^180 + u^
#                178 + 3*u^176 + 4*u^174 + 6*u^172 + 7*u^170 + 7*u^168 + 7*u^
#                166 + 5*u^164 + 4*u^162 + 3*u^160 - 2*u^156 - 4*u^154 - 6*u^
#                152 - 6*u^150 - 8*u^148 - 9*u^146 - 8*u^144 - 8*u^142 - 5*u^
#                140 - 2*u^138 - u^136 + 2*u^134 + 3*u^132 + 4*u^130 + 7*u^
#                128 + 7*u^126 + 7*u^124 + 7*u^122 + 5*u^120 + 5*u^118 + 3*u^
#                116 - u^112 - 3*u^110 - 4*u^108 - 4*u^106 - 5*u^104 - 4*u^
#                102 - 3*u^100 - 3*u^98 - u^96 + u^92 + 2*u^90 + 2*u^88 + 2*u^
#                86 + 2*u^84 + u^82 - u^78 - u^76, 
#              u^207 + 3*u^205 + 3*u^203 + 2*u^201 + u^199 - u^197 - 2*u^195 - 
#                4*u^193 - 6*u^191 - 6*u^189 - 7*u^187 - 6*u^185 - 4*u^183 - 
#                4*u^181 + 3*u^177 + 4*u^175 + 9*u^173 + 10*u^171 + 11*u^169 + 
#                13*u^167 + 9*u^165 + 9*u^163 + 8*u^161 + 2*u^159 + u^157 - 
#                4*u^155 - 8*u^153 - 7*u^151 - 13*u^149 - 14*u^147 - 13*u^
#                145 - 16*u^143 - 10*u^141 - 7*u^139 - 6*u^137 + u^135 + u^
#                133 + 4*u^131 + 10*u^129 + 9*u^127 + 12*u^125 + 12*u^123 + 
#                9*u^121 + 11*u^119 + 7*u^117 + 3*u^115 + 2*u^113 - 3*u^111 - 
#                4*u^109 - 5*u^107 - 8*u^105 - 6*u^103 - 6*u^101 - 6*u^99 - 
#                3*u^97 - 2*u^95 + 2*u^91 + 2*u^89 + 3*u^87 + 3*u^85 + 2*u^
#                83 + u^81 - u^79 - u^77, 
#              u^208 + 2*u^206 + 2*u^204 + u^202 - u^198 - 2*u^196 - 3*u^194 - 
#                4*u^192 - 4*u^190 - 4*u^188 - 3*u^186 - 2*u^184 - u^182 + u^
#                180 + 3*u^178 + 4*u^176 + 6*u^174 + 7*u^172 + 7*u^170 + 7*u^
#                168 + 5*u^166 + 4*u^164 + 3*u^162 - 2*u^158 - 4*u^156 - 6*u^
#                154 - 6*u^152 - 8*u^150 - 9*u^148 - 8*u^146 - 8*u^144 - 5*u^
#                142 - 2*u^140 - u^138 + 2*u^136 + 3*u^134 + 4*u^132 + 7*u^
#                130 + 7*u^128 + 7*u^126 + 7*u^124 + 5*u^122 + 5*u^120 + 3*u^
#                118 - u^114 - 3*u^112 - 4*u^110 - 4*u^108 - 5*u^106 - 4*u^
#                104 - 3*u^102 - 3*u^100 - u^98 + u^94 + 2*u^92 + 2*u^90 + 2*u^
#                88 + 2*u^86 + u^84 - u^80 - u^78, 
#              2*u^206 + 2*u^204 + u^202 + u^200 - u^198 - u^196 - 2*u^194 - 
#                4*u^192 - 3*u^190 - 4*u^188 - 4*u^186 - u^184 - 3*u^182 + 3*u^
#                178 + u^176 + 6*u^174 + 6*u^172 + 5*u^170 + 9*u^168 + 4*u^
#                166 + 4*u^164 + 6*u^162 - u^160 + u^158 - 2*u^156 - 7*u^154 - 
#                2*u^152 - 8*u^150 - 9*u^148 - 5*u^146 - 11*u^144 - 5*u^142 - 
#                2*u^140 - 5*u^138 + 3*u^136 + u^134 + u^132 + 8*u^130 + 4*u^
#                128 + 7*u^126 + 8*u^124 + 3*u^122 + 7*u^120 + 4*u^118 + 2*u^
#                114 - 3*u^112 - 3*u^110 - 2*u^108 - 6*u^106 - 3*u^104 - 3*u^
#                102 - 4*u^100 - u^98 - u^96 + 2*u^92 + u^90 + 2*u^88 + 2*u^
#                86 + u^84 + u^82 - u^80 - u^78, 
#              u^206 + 2*u^204 + u^202 + u^200 - u^196 - u^194 - 3*u^192 - 3*u^
#                190 - 3*u^188 - 4*u^186 - 2*u^184 - 2*u^182 - 2*u^180 + 2*u^
#                178 + u^176 + 3*u^174 + 6*u^172 + 4*u^170 + 7*u^168 + 6*u^
#                166 + 3*u^164 + 6*u^162 + 2*u^160 + u^156 - 5*u^154 - 3*u^
#                152 - 4*u^150 - 9*u^148 - 5*u^146 - 8*u^144 - 8*u^142 - 2*u^
#                140 - 5*u^138 - u^136 + 2*u^134 - u^132 + 5*u^130 + 5*u^128 + 
#                4*u^126 + 8*u^124 + 4*u^122 + 5*u^120 + 6*u^118 + u^116 + 2*u^
#                114 - 3*u^110 - u^108 - 4*u^106 - 4*u^104 - 2*u^102 - 4*u^
#                100 - 2*u^98 - u^96 - u^94 + u^92 + u^90 + u^88 + 2*u^86 + u^
#                84 + u^82 - u^78, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^
#                192 - 2*u^190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^
#                180 + 3*u^174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^
#                164 + 5*u^162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^
#                150 - 5*u^148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^
#                138 - u^136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^
#                124 + 4*u^122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^
#                108 - 3*u^106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^
#                96 - u^94 + u^88 + u^86 + u^84 + u^82, 
#              u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^195 - 3*u^193 - 3*u^
#                191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^183 - 2*u^181 + 2*u^
#                179 + u^177 + 3*u^175 + 6*u^173 + 4*u^171 + 7*u^169 + 6*u^
#                167 + 3*u^165 + 6*u^163 + 2*u^161 + u^157 - 5*u^155 - 3*u^
#                153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^145 - 8*u^143 - 2*u^
#                141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 5*u^131 + 5*u^129 + 
#                4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 6*u^119 + u^117 + 2*u^
#                115 - 3*u^111 - u^109 - 4*u^107 - 4*u^105 - 2*u^103 - 4*u^
#                101 - 2*u^99 - u^97 - u^95 + u^93 + u^91 + u^89 + 2*u^87 + u^
#                85 + u^83 - u^79, u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^
#                195 - 3*u^193 - 3*u^191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^
#                183 - 2*u^181 + 2*u^179 + u^177 + 3*u^175 + 6*u^173 + 4*u^
#                171 + 7*u^169 + 6*u^167 + 3*u^165 + 6*u^163 + 2*u^161 + u^
#                157 - 5*u^155 - 3*u^153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^
#                145 - 8*u^143 - 2*u^141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 
#                5*u^131 + 5*u^129 + 4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 
#                6*u^119 + u^117 + 2*u^115 - 3*u^111 - u^109 - 4*u^107 - 4*u^
#                105 - 2*u^103 - 4*u^101 - 2*u^99 - u^97 - u^95 + u^93 + u^
#                91 + u^89 + 2*u^87 + u^85 + u^83 - u^79, 
#              u^207 + u^205 - u^199 - u^197 - u^195 - 2*u^193 - u^191 - u^
#                189 - u^187 + u^185 + u^181 + 3*u^179 + u^177 + 3*u^175 + 3*u^
#                173 + u^171 + 3*u^169 - u^165 + u^163 - 3*u^161 - 2*u^159 - 
#                2*u^157 - 5*u^155 - u^153 - 3*u^151 - 4*u^149 - 3*u^145 + 3*u^
#                141 + 4*u^137 + 3*u^135 + u^133 + 5*u^131 + 2*u^129 + 2*u^
#                127 + 3*u^125 - u^123 + u^121 - 3*u^117 - u^115 - 3*u^113 - 
#                3*u^111 - u^109 - 3*u^107 - u^105 - u^101 + u^99 + u^97 + u^
#                95 + 2*u^93 + u^91 + u^89 + u^87 - u^81 - u^79, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 
#                2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^
#                164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 
#                3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^
#                128 + 3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^
#                112 - u^108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^
#                92 + u^88 - u^80, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^
#                192 - 2*u^190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^
#                180 + 3*u^174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^
#                164 + 5*u^162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^
#                150 - 5*u^148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^
#                138 - u^136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^
#                124 + 4*u^122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^
#                108 - 3*u^106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^
#                96 - u^94 + u^88 + u^86 + u^84 + u^82, 0*u^0, 0*u^0, 0*u^0, 
#              0*u^0, 0*u^0 ], 
#          [ u^206 + 2*u^204 + 3*u^202 + u^200 + u^198 - u^196 - 2*u^194 - 3*u^
#                192 - 5*u^190 - 5*u^188 - 5*u^186 - 5*u^184 - 3*u^182 - 2*u^
#                180 - u^178 + 4*u^176 + 3*u^174 + 7*u^172 + 9*u^170 + 8*u^
#                168 + 10*u^166 + 8*u^164 + 5*u^162 + 7*u^160 + u^158 - u^
#                156 - 2*u^154 - 8*u^152 - 6*u^150 - 9*u^148 - 13*u^146 - 9*u^
#                144 - 12*u^142 - 9*u^140 - 3*u^138 - 5*u^136 + u^134 + 3*u^
#                132 + 2*u^130 + 9*u^128 + 8*u^126 + 8*u^124 + 11*u^122 + 6*u^
#                120 + 8*u^118 + 6*u^116 + u^114 + u^112 - 2*u^110 - 5*u^108 - 
#                3*u^106 - 7*u^104 - 5*u^102 - 4*u^100 - 5*u^98 - 2*u^96 - u^
#                94 + 2*u^90 + 2*u^88 + 2*u^86 + 3*u^84 + u^82 + u^80 - u^
#                78 - u^76, u^207 + 4*u^205 + 4*u^203 + 3*u^201 + 2*u^199 - u^
#                197 - 2*u^195 - 5*u^193 - 8*u^191 - 8*u^189 - 10*u^187 - 9*u^
#                185 - 6*u^183 - 7*u^181 - u^179 + 3*u^177 + 4*u^175 + 12*u^
#                173 + 13*u^171 + 15*u^169 + 19*u^167 + 13*u^165 + 14*u^163 + 
#                13*u^161 + 4*u^159 + 4*u^157 - 4*u^155 - 10*u^153 - 8*u^151 - 
#                18*u^149 - 19*u^147 - 18*u^145 - 24*u^143 - 15*u^141 - 12*u^
#                139 - 11*u^137 - u^133 + 4*u^131 + 13*u^129 + 11*u^127 + 17*u^
#                125 + 17*u^123 + 13*u^121 + 17*u^119 + 11*u^117 + 6*u^115 + 
#                5*u^113 - 3*u^111 - 4*u^109 - 6*u^107 - 11*u^105 - 8*u^103 - 
#                9*u^101 - 9*u^99 - 5*u^97 - 4*u^95 - u^93 + 2*u^91 + 2*u^89 + 
#                4*u^87 + 4*u^85 + 3*u^83 + 2*u^81 - u^79 - u^77, 
#              2*u^206 + 2*u^204 + u^202 + u^200 - u^198 - u^196 - 2*u^194 - 
#                4*u^192 - 3*u^190 - 4*u^188 - 4*u^186 - u^184 - 3*u^182 + 3*u^
#                178 + u^176 + 6*u^174 + 6*u^172 + 5*u^170 + 9*u^168 + 4*u^
#                166 + 4*u^164 + 6*u^162 - u^160 + u^158 - 2*u^156 - 7*u^154 - 
#                2*u^152 - 8*u^150 - 9*u^148 - 5*u^146 - 11*u^144 - 5*u^142 - 
#                2*u^140 - 5*u^138 + 3*u^136 + u^134 + u^132 + 8*u^130 + 4*u^
#                128 + 7*u^126 + 8*u^124 + 3*u^122 + 7*u^120 + 4*u^118 + 2*u^
#                114 - 3*u^112 - 3*u^110 - 2*u^108 - 6*u^106 - 3*u^104 - 3*u^
#                102 - 4*u^100 - u^98 - u^96 + 2*u^92 + u^90 + 2*u^88 + 2*u^
#                86 + u^84 + u^82 - u^80 - u^78, 
#              u^208 + 2*u^206 + 4*u^204 + u^202 + u^200 - u^198 - 3*u^196 - 
#                3*u^194 - 6*u^192 - 6*u^190 - 5*u^188 - 6*u^186 - 3*u^184 - u^
#                182 - 2*u^180 + 6*u^178 + 4*u^176 + 7*u^174 + 12*u^172 + 8*u^
#                170 + 11*u^168 + 10*u^166 + 3*u^164 + 8*u^162 + u^160 - 4*u^
#                158 - u^156 - 11*u^154 - 8*u^152 - 8*u^150 - 17*u^148 - 9*u^
#                146 - 12*u^144 - 12*u^142 - 5*u^138 + u^136 + 7*u^134 + u^
#                132 + 11*u^130 + 11*u^128 + 7*u^126 + 14*u^124 + 6*u^122 + 
#                7*u^120 + 8*u^118 - u^116 - 2*u^112 - 8*u^110 - 3*u^108 - 8*u^
#                106 - 7*u^104 - 3*u^102 - 6*u^100 - 2*u^98 + 3*u^92 + 3*u^
#                90 + 2*u^88 + 4*u^86 + u^84 + u^82 - u^80 - 2*u^78, 
#              u^206 + 3*u^204 + 2*u^202 + 2*u^200 + u^198 - u^196 - u^194 - 
#                4*u^192 - 5*u^190 - 5*u^188 - 7*u^186 - 5*u^184 - 4*u^182 - 
#                5*u^180 + u^178 + u^176 + 3*u^174 + 9*u^172 + 7*u^170 + 11*u^
#                168 + 12*u^166 + 7*u^164 + 11*u^162 + 7*u^160 + 2*u^158 + 4*u^
#                156 - 5*u^154 - 5*u^152 - 5*u^150 - 14*u^148 - 10*u^146 - 
#                13*u^144 - 16*u^142 - 7*u^140 - 10*u^138 - 6*u^136 + u^134 - 
#                3*u^132 + 5*u^130 + 8*u^128 + 6*u^126 + 13*u^124 + 9*u^122 + 
#                9*u^120 + 12*u^118 + 5*u^116 + 5*u^114 + 3*u^112 - 3*u^
#                110 - u^108 - 5*u^106 - 7*u^104 - 4*u^102 - 7*u^100 - 5*u^
#                98 - 3*u^96 - 3*u^94 + u^90 + u^88 + 3*u^86 + 2*u^84 + 2*u^
#                82 + u^80 - u^78, u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^
#                198 - u^194 - 3*u^192 - 4*u^190 - 5*u^188 - 6*u^186 - 5*u^
#                184 - 5*u^182 - 4*u^180 - u^178 + 3*u^174 + 6*u^172 + 7*u^
#                170 + 10*u^168 + 10*u^166 + 9*u^164 + 10*u^162 + 7*u^160 + 
#                5*u^158 + 3*u^156 - 2*u^154 - 3*u^152 - 6*u^150 - 10*u^148 - 
#                10*u^146 - 13*u^144 - 13*u^142 - 10*u^140 - 10*u^138 - 6*u^
#                136 - 3*u^134 - 2*u^132 + 3*u^130 + 5*u^128 + 7*u^126 + 10*u^
#                124 + 9*u^122 + 10*u^120 + 10*u^118 + 7*u^116 + 6*u^114 + 3*u^
#                112 - u^108 - 4*u^106 - 5*u^104 - 5*u^102 - 6*u^100 - 5*u^
#                98 - 4*u^96 - 3*u^94 - u^92 + u^88 + 2*u^86 + 2*u^84 + 2*u^
#                82 + u^80, u^207 + 3*u^205 + 2*u^203 + 2*u^201 + u^199 - u^
#                197 - u^195 - 4*u^193 - 5*u^191 - 5*u^189 - 7*u^187 - 5*u^
#                185 - 4*u^183 - 5*u^181 + u^179 + u^177 + 3*u^175 + 9*u^173 + 
#                7*u^171 + 11*u^169 + 12*u^167 + 7*u^165 + 11*u^163 + 7*u^
#                161 + 2*u^159 + 4*u^157 - 5*u^155 - 5*u^153 - 5*u^151 - 14*u^
#                149 - 10*u^147 - 13*u^145 - 16*u^143 - 7*u^141 - 10*u^139 - 
#                6*u^137 + u^135 - 3*u^133 + 5*u^131 + 8*u^129 + 6*u^127 + 
#                13*u^125 + 9*u^123 + 9*u^121 + 12*u^119 + 5*u^117 + 5*u^115 + 
#                3*u^113 - 3*u^111 - u^109 - 5*u^107 - 7*u^105 - 4*u^103 - 7*u^
#                101 - 5*u^99 - 3*u^97 - 3*u^95 + u^91 + u^89 + 3*u^87 + 2*u^
#                85 + 2*u^83 + u^81 - u^79, 
#              u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^195 - 3*u^193 - 3*u^
#                191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^183 - 2*u^181 + 2*u^
#                179 + u^177 + 3*u^175 + 6*u^173 + 4*u^171 + 7*u^169 + 6*u^
#                167 + 3*u^165 + 6*u^163 + 2*u^161 + u^157 - 5*u^155 - 3*u^
#                153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^145 - 8*u^143 - 2*u^
#                141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 5*u^131 + 5*u^129 + 
#                4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 6*u^119 + u^117 + 2*u^
#                115 - 3*u^111 - u^109 - 4*u^107 - 4*u^105 - 2*u^103 - 4*u^
#                101 - 2*u^99 - u^97 - u^95 + u^93 + u^91 + u^89 + 2*u^87 + u^
#                85 + u^83 - u^79, u^207 + 2*u^205 - u^199 - 2*u^197 - u^195 - 
#                3*u^193 - 2*u^191 - u^189 - 2*u^187 + u^185 + u^183 + 5*u^
#                179 + 2*u^177 + 3*u^175 + 6*u^173 + u^171 + 4*u^169 + 2*u^
#                167 - 3*u^165 + 2*u^163 - 3*u^161 - 5*u^159 - u^157 - 8*u^
#                155 - 3*u^153 - 2*u^151 - 8*u^149 - 3*u^145 - 3*u^143 + 6*u^
#                141 + 4*u^137 + 7*u^135 + 7*u^131 + 5*u^129 + u^127 + 6*u^
#                125 - u^123 + 2*u^119 - 5*u^117 - 2*u^115 - 3*u^113 - 6*u^
#                111 - u^109 - 4*u^107 - 3*u^105 + u^103 - 2*u^101 + u^99 + 
#                2*u^97 + u^95 + 3*u^93 + 2*u^91 + u^89 + 2*u^87 - u^81 - 2*u^
#                79, u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 
#                2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^
#                164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 
#                3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^
#                128 + 3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^
#                112 - u^108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^
#                92 + u^88 - u^80, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 
#              2*u^206 + u^204 + 2*u^202 + u^200 - 2*u^194 - 3*u^192 - 3*u^
#                190 - 5*u^188 - 4*u^186 - 3*u^184 - 5*u^182 - u^178 + u^176 + 
#                5*u^174 + 4*u^172 + 7*u^170 + 9*u^168 + 5*u^166 + 9*u^164 + 
#                6*u^162 + 3*u^160 + 5*u^158 - 2*u^156 - 2*u^154 - 2*u^152 - 
#                9*u^150 - 6*u^148 - 9*u^146 - 12*u^144 - 6*u^142 - 9*u^140 - 
#                6*u^138 - u^136 - 4*u^134 + 2*u^132 + 4*u^130 + 3*u^128 + 9*u^
#                126 + 6*u^124 + 7*u^122 + 9*u^120 + 5*u^118 + 5*u^116 + 4*u^
#                114 - u^112 + u^110 - 3*u^108 - 4*u^106 - 3*u^104 - 5*u^102 - 
#                4*u^100 - 3*u^98 - 3*u^96 - u^94 + 2*u^88 + u^86 + 2*u^84 + u^
#                82, u^206 + u^202 - u^194 - u^192 - u^190 - 2*u^188 - u^
#                186 - u^184 - 2*u^182 + u^180 - u^178 + u^176 + 2*u^174 + u^
#                172 + 3*u^170 + 3*u^168 + u^166 + 4*u^164 + u^162 + u^160 + 
#                2*u^158 - 2*u^156 - u^152 - 4*u^150 - u^148 - 4*u^146 - 4*u^
#                144 - u^142 - 4*u^140 - u^138 - 2*u^134 + 2*u^132 + u^130 + u^
#                128 + 4*u^126 + u^124 + 3*u^122 + 3*u^120 + u^118 + 2*u^
#                116 + u^114 - u^112 + u^110 - 2*u^108 - u^106 - u^104 - 2*u^
#                102 - u^100 - u^98 - u^96 + u^88 + u^84, 
#              u^204 + u^200 - u^192 - u^190 - u^188 - 2*u^186 - u^184 - u^
#                182 - 2*u^180 + u^178 - u^176 + u^174 + 2*u^172 + u^170 + 3*u^
#                168 + 3*u^166 + u^164 + 4*u^162 + u^160 + u^158 + 2*u^156 - 
#                2*u^154 - u^150 - 4*u^148 - u^146 - 4*u^144 - 4*u^142 - u^
#                140 - 4*u^138 - u^136 - 2*u^132 + 2*u^130 + u^128 + u^126 + 
#                4*u^124 + u^122 + 3*u^120 + 3*u^118 + u^116 + 2*u^114 + u^
#                112 - u^110 + u^108 - 2*u^106 - u^104 - u^102 - 2*u^100 - u^
#                98 - u^96 - u^94 + u^86 + u^82, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, 0*u^0, 0*u^0 ], 
#          [ u^206 + 2*u^204 + 3*u^202 + 2*u^200 + u^198 - 2*u^194 - 3*u^192 - 
#                5*u^190 - 6*u^188 - 6*u^186 - 6*u^184 - 5*u^182 - 3*u^180 - 
#                2*u^178 + 2*u^176 + 4*u^174 + 6*u^172 + 10*u^170 + 10*u^168 + 
#                11*u^166 + 11*u^164 + 8*u^162 + 8*u^160 + 5*u^158 - u^154 - 
#                6*u^152 - 8*u^150 - 9*u^148 - 14*u^146 - 13*u^144 - 13*u^
#                142 - 13*u^140 - 7*u^138 - 6*u^136 - 3*u^134 + 2*u^132 + 2*u^
#                130 + 7*u^128 + 10*u^126 + 9*u^124 + 12*u^122 + 10*u^120 + 
#                9*u^118 + 9*u^116 + 4*u^114 + 2*u^112 - 4*u^108 - 4*u^106 - 
#                6*u^104 - 7*u^102 - 5*u^100 - 6*u^98 - 4*u^96 - 2*u^94 - u^
#                92 + u^90 + 2*u^88 + 2*u^86 + 3*u^84 + 2*u^82 + u^80 - u^76, 
#              u^207 + 3*u^205 + 5*u^203 + 4*u^201 + 3*u^199 + u^197 - 2*u^
#                195 - 4*u^193 - 8*u^191 - 10*u^189 - 11*u^187 - 12*u^185 - 
#                10*u^183 - 8*u^181 - 6*u^179 + u^177 + 4*u^175 + 9*u^173 + 
#                16*u^171 + 17*u^169 + 21*u^167 + 21*u^165 + 17*u^163 + 18*u^
#                161 + 12*u^159 + 5*u^157 + 2*u^155 - 8*u^153 - 11*u^151 - 
#                15*u^149 - 24*u^147 - 23*u^145 - 26*u^143 - 26*u^141 - 17*u^
#                139 - 16*u^137 - 9*u^135 - u^133 + 10*u^129 + 15*u^127 + 16*u^
#                125 + 22*u^123 + 19*u^121 + 19*u^119 + 19*u^117 + 11*u^115 + 
#                8*u^113 + 3*u^111 - 4*u^109 - 5*u^107 - 10*u^105 - 12*u^103 - 
#                10*u^101 - 12*u^99 - 9*u^97 - 6*u^95 - 4*u^93 + 2*u^89 + 3*u^
#                87 + 5*u^85 + 4*u^83 + 3*u^81 + u^79 - u^77, 
#              u^206 + 2*u^204 + u^202 + u^200 - u^196 - u^194 - 3*u^192 - 3*u^
#                190 - 3*u^188 - 4*u^186 - 2*u^184 - 2*u^182 - 2*u^180 + 2*u^
#                178 + u^176 + 3*u^174 + 6*u^172 + 4*u^170 + 7*u^168 + 6*u^
#                166 + 3*u^164 + 6*u^162 + 2*u^160 + u^156 - 5*u^154 - 3*u^
#                152 - 4*u^150 - 9*u^148 - 5*u^146 - 8*u^144 - 8*u^142 - 2*u^
#                140 - 5*u^138 - u^136 + 2*u^134 - u^132 + 5*u^130 + 5*u^128 + 
#                4*u^126 + 8*u^124 + 4*u^122 + 5*u^120 + 6*u^118 + u^116 + 2*u^
#                114 - 3*u^110 - u^108 - 4*u^106 - 4*u^104 - 2*u^102 - 4*u^
#                100 - 2*u^98 - u^96 - u^94 + u^92 + u^90 + u^88 + 2*u^86 + u^
#                84 + u^82 - u^78, u^206 + 3*u^204 + 2*u^202 + 2*u^200 + u^
#                198 - u^196 - u^194 - 4*u^192 - 5*u^190 - 5*u^188 - 7*u^186 - 
#                5*u^184 - 4*u^182 - 5*u^180 + u^178 + u^176 + 3*u^174 + 9*u^
#                172 + 7*u^170 + 11*u^168 + 12*u^166 + 7*u^164 + 11*u^162 + 
#                7*u^160 + 2*u^158 + 4*u^156 - 5*u^154 - 5*u^152 - 5*u^150 - 
#                14*u^148 - 10*u^146 - 13*u^144 - 16*u^142 - 7*u^140 - 10*u^
#                138 - 6*u^136 + u^134 - 3*u^132 + 5*u^130 + 8*u^128 + 6*u^
#                126 + 13*u^124 + 9*u^122 + 9*u^120 + 12*u^118 + 5*u^116 + 5*u^
#                114 + 3*u^112 - 3*u^110 - u^108 - 5*u^106 - 7*u^104 - 4*u^
#                102 - 7*u^100 - 5*u^98 - 3*u^96 - 3*u^94 + u^90 + u^88 + 3*u^
#                86 + 2*u^84 + 2*u^82 + u^80 - u^78, 
#              u^208 + u^206 + 4*u^204 + 4*u^202 + 2*u^200 + 2*u^198 - 2*u^
#                196 - 3*u^194 - 5*u^192 - 9*u^190 - 8*u^188 - 9*u^186 - 10*u^
#                184 - 4*u^182 - 6*u^180 - u^178 + 6*u^176 + 4*u^174 + 13*u^
#                172 + 15*u^170 + 13*u^168 + 20*u^166 + 13*u^164 + 11*u^162 + 
#                14*u^160 + u^158 + 2*u^156 - 3*u^154 - 14*u^152 - 8*u^150 - 
#                18*u^148 - 22*u^146 - 15*u^144 - 24*u^142 - 15*u^140 - 8*u^
#                138 - 12*u^136 + 2*u^134 + 2*u^132 + 3*u^130 + 16*u^128 + 
#                11*u^126 + 16*u^124 + 19*u^122 + 11*u^120 + 16*u^118 + 11*u^
#                116 + 3*u^114 + 5*u^112 - 4*u^110 - 6*u^108 - 5*u^106 - 12*u^
#                104 - 8*u^102 - 8*u^100 - 9*u^98 - 4*u^96 - 3*u^94 - u^92 + 
#                3*u^90 + 2*u^88 + 4*u^86 + 4*u^84 + 2*u^82 + 2*u^80 - u^
#                78 - u^76, 2*u^206 + 3*u^204 + 3*u^202 + 3*u^200 + u^198 - 
#                2*u^194 - 5*u^192 - 6*u^190 - 8*u^188 - 9*u^186 - 7*u^184 - 
#                8*u^182 - 5*u^180 - u^178 + 6*u^174 + 9*u^172 + 11*u^170 + 
#                16*u^168 + 14*u^166 + 14*u^164 + 15*u^162 + 9*u^160 + 8*u^
#                158 + 3*u^156 - 4*u^154 - 4*u^152 - 11*u^150 - 15*u^148 - 
#                15*u^146 - 21*u^144 - 18*u^142 - 15*u^140 - 15*u^138 - 7*u^
#                136 - 5*u^134 - 2*u^132 + 6*u^130 + 7*u^128 + 12*u^126 + 15*u^
#                124 + 13*u^122 + 16*u^120 + 14*u^118 + 10*u^116 + 9*u^114 + 
#                3*u^112 - 2*u^108 - 7*u^106 - 7*u^104 - 8*u^102 - 9*u^100 - 
#                7*u^98 - 6*u^96 - 4*u^94 - u^92 + 2*u^88 + 3*u^86 + 3*u^84 + 
#                3*u^82 + u^80, u^207 + 4*u^205 + 4*u^203 + 3*u^201 + 2*u^
#                199 - u^197 - 2*u^195 - 5*u^193 - 8*u^191 - 8*u^189 - 10*u^
#                187 - 9*u^185 - 6*u^183 - 7*u^181 - u^179 + 3*u^177 + 4*u^
#                175 + 12*u^173 + 13*u^171 + 15*u^169 + 19*u^167 + 13*u^165 + 
#                14*u^163 + 13*u^161 + 4*u^159 + 4*u^157 - 4*u^155 - 10*u^
#                153 - 8*u^151 - 18*u^149 - 19*u^147 - 18*u^145 - 24*u^143 - 
#                15*u^141 - 12*u^139 - 11*u^137 - u^133 + 4*u^131 + 13*u^129 + 
#                11*u^127 + 17*u^125 + 17*u^123 + 13*u^121 + 17*u^119 + 11*u^
#                117 + 6*u^115 + 5*u^113 - 3*u^111 - 4*u^109 - 6*u^107 - 11*u^
#                105 - 8*u^103 - 9*u^101 - 9*u^99 - 5*u^97 - 4*u^95 - u^93 + 
#                2*u^91 + 2*u^89 + 4*u^87 + 4*u^85 + 3*u^83 + 2*u^81 - u^
#                79 - u^77, u^205 + u^203 + u^201 + u^199 - u^193 - 2*u^191 - 
#                2*u^189 - 3*u^187 - 3*u^185 - 2*u^183 - 3*u^181 - u^179 + 3*u^
#                173 + 3*u^171 + 4*u^169 + 6*u^167 + 4*u^165 + 5*u^163 + 5*u^
#                161 + 2*u^159 + 3*u^157 - 2*u^153 - u^151 - 5*u^149 - 5*u^
#                147 - 5*u^145 - 8*u^143 - 5*u^141 - 5*u^139 - 5*u^137 - u^
#                135 - 2*u^133 + 3*u^129 + 2*u^127 + 5*u^125 + 5*u^123 + 4*u^
#                121 + 6*u^119 + 4*u^117 + 3*u^115 + 3*u^113 - u^107 - 3*u^
#                105 - 2*u^103 - 3*u^101 - 3*u^99 - 2*u^97 - 2*u^95 - u^93 + u^
#                87 + u^85 + u^83 + u^81, u^205 - u^197 - u^193 - u^191 - u^
#                187 + u^183 - u^181 + 2*u^179 + u^177 + 3*u^173 + u^169 + 2*u^
#                167 - 2*u^165 + u^163 - 3*u^159 + u^157 - 3*u^155 - 2*u^
#                153 + u^151 - 4*u^149 - 3*u^143 + 3*u^141 + 4*u^135 - u^133 + 
#                2*u^131 + 3*u^129 - u^127 + 3*u^125 - u^121 + 2*u^119 - 2*u^
#                117 - u^115 - 3*u^111 - u^107 - 2*u^105 + u^103 - u^101 + u^
#                97 + u^93 + u^91 + u^87 - u^79, 
#              u^206 + u^204 - u^198 - u^196 - u^194 - 2*u^192 - u^190 - u^
#                188 - u^186 + u^184 + u^180 + 3*u^178 + u^176 + 3*u^174 + 3*u^
#                172 + u^170 + 3*u^168 - u^164 + u^162 - 3*u^160 - 2*u^158 - 
#                2*u^156 - 5*u^154 - u^152 - 3*u^150 - 4*u^148 - 3*u^144 + 3*u^
#                140 + 4*u^136 + 3*u^134 + u^132 + 5*u^130 + 2*u^128 + 2*u^
#                126 + 3*u^124 - u^122 + u^120 - 3*u^116 - u^114 - 3*u^112 - 
#                3*u^110 - u^108 - 3*u^106 - u^104 - u^100 + u^98 + u^96 + u^
#                94 + 2*u^92 + u^90 + u^88 + u^86 - u^80 - u^78, 0*u^0, 
#              u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^198 - u^194 - 3*u^192 - 
#                4*u^190 - 5*u^188 - 6*u^186 - 5*u^184 - 5*u^182 - 4*u^180 - u^
#                178 + 3*u^174 + 6*u^172 + 7*u^170 + 10*u^168 + 10*u^166 + 9*u^
#                164 + 10*u^162 + 7*u^160 + 5*u^158 + 3*u^156 - 2*u^154 - 3*u^
#                152 - 6*u^150 - 10*u^148 - 10*u^146 - 13*u^144 - 13*u^142 - 
#                10*u^140 - 10*u^138 - 6*u^136 - 3*u^134 - 2*u^132 + 3*u^130 + 
#                5*u^128 + 7*u^126 + 10*u^124 + 9*u^122 + 10*u^120 + 10*u^
#                118 + 7*u^116 + 6*u^114 + 3*u^112 - u^108 - 4*u^106 - 5*u^
#                104 - 5*u^102 - 6*u^100 - 5*u^98 - 4*u^96 - 3*u^94 - u^92 + u^
#                88 + 2*u^86 + 2*u^84 + 2*u^82 + u^80, 
#              u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^192 - 2*u^190 - 3*u^
#                188 - 3*u^186 - 2*u^184 - 3*u^182 - u^180 + 3*u^174 + 3*u^
#                172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^164 + 5*u^162 + 2*u^
#                160 + 3*u^158 - 2*u^154 - u^152 - 5*u^150 - 5*u^148 - 5*u^
#                146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^138 - u^136 - 2*u^
#                134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^124 + 4*u^122 + 6*u^
#                120 + 4*u^118 + 3*u^116 + 3*u^114 - u^108 - 3*u^106 - 2*u^
#                104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^96 - u^94 + u^88 + u^
#                86 + u^84 + u^82, u^206 + 2*u^204 + u^202 + u^200 - u^196 - u^
#                194 - 3*u^192 - 3*u^190 - 3*u^188 - 4*u^186 - 2*u^184 - 2*u^
#                182 - 2*u^180 + 2*u^178 + u^176 + 3*u^174 + 6*u^172 + 4*u^
#                170 + 7*u^168 + 6*u^166 + 3*u^164 + 6*u^162 + 2*u^160 + u^
#                156 - 5*u^154 - 3*u^152 - 4*u^150 - 9*u^148 - 5*u^146 - 8*u^
#                144 - 8*u^142 - 2*u^140 - 5*u^138 - u^136 + 2*u^134 - u^132 + 
#                5*u^130 + 5*u^128 + 4*u^126 + 8*u^124 + 4*u^122 + 5*u^120 + 
#                6*u^118 + u^116 + 2*u^114 - 3*u^110 - u^108 - 4*u^106 - 4*u^
#                104 - 2*u^102 - 4*u^100 - 2*u^98 - u^96 - u^94 + u^92 + u^
#                90 + u^88 + 2*u^86 + u^84 + u^82 - u^78, 0*u^0, 
#              u^205 - u^197 - u^193 - u^191 - u^187 + u^183 - u^181 + 2*u^
#                179 + u^177 + 3*u^173 + u^169 + 2*u^167 - 2*u^165 + u^163 - 
#                3*u^159 + u^157 - 3*u^155 - 2*u^153 + u^151 - 4*u^149 - 3*u^
#                143 + 3*u^141 + 4*u^135 - u^133 + 2*u^131 + 3*u^129 - u^127 + 
#                3*u^125 - u^121 + 2*u^119 - 2*u^117 - u^115 - 3*u^111 - u^
#                107 - 2*u^105 + u^103 - u^101 + u^97 + u^93 + u^91 + u^87 - u^
#                79, 0*u^0 ], 
#          [ u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^198 - u^194 - 3*u^192 - 
#                4*u^190 - 5*u^188 - 6*u^186 - 5*u^184 - 5*u^182 - 4*u^180 - u^
#                178 + 3*u^174 + 6*u^172 + 7*u^170 + 10*u^168 + 10*u^166 + 9*u^
#                164 + 10*u^162 + 7*u^160 + 5*u^158 + 3*u^156 - 2*u^154 - 3*u^
#                152 - 6*u^150 - 10*u^148 - 10*u^146 - 13*u^144 - 13*u^142 - 
#                10*u^140 - 10*u^138 - 6*u^136 - 3*u^134 - 2*u^132 + 3*u^130 + 
#                5*u^128 + 7*u^126 + 10*u^124 + 9*u^122 + 10*u^120 + 10*u^
#                118 + 7*u^116 + 6*u^114 + 3*u^112 - u^108 - 4*u^106 - 5*u^
#                104 - 5*u^102 - 6*u^100 - 5*u^98 - 4*u^96 - 3*u^94 - u^92 + u^
#                88 + 2*u^86 + 2*u^84 + 2*u^82 + u^80, 
#              u^207 + 3*u^205 + 4*u^203 + 4*u^201 + 3*u^199 + u^197 - u^195 - 
#                4*u^193 - 7*u^191 - 9*u^189 - 11*u^187 - 11*u^185 - 10*u^
#                183 - 9*u^181 - 5*u^179 - u^177 + 3*u^175 + 9*u^173 + 13*u^
#                171 + 17*u^169 + 20*u^167 + 19*u^165 + 19*u^163 + 17*u^161 + 
#                12*u^159 + 8*u^157 + u^155 - 5*u^153 - 9*u^151 - 16*u^149 - 
#                20*u^147 - 23*u^145 - 26*u^143 - 23*u^141 - 20*u^139 - 16*u^
#                137 - 9*u^135 - 5*u^133 + u^131 + 8*u^129 + 12*u^127 + 17*u^
#                125 + 19*u^123 + 19*u^121 + 20*u^119 + 17*u^117 + 13*u^115 + 
#                9*u^113 + 3*u^111 - u^109 - 5*u^107 - 9*u^105 - 10*u^103 - 
#                11*u^101 - 11*u^99 - 9*u^97 - 7*u^95 - 4*u^93 - u^91 + u^89 + 
#                3*u^87 + 4*u^85 + 4*u^83 + 3*u^81 + u^79, 
#              u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^192 - 2*u^190 - 3*u^
#                188 - 3*u^186 - 2*u^184 - 3*u^182 - u^180 + 3*u^174 + 3*u^
#                172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^164 + 5*u^162 + 2*u^
#                160 + 3*u^158 - 2*u^154 - u^152 - 5*u^150 - 5*u^148 - 5*u^
#                146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^138 - u^136 - 2*u^
#                134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^124 + 4*u^122 + 6*u^
#                120 + 4*u^118 + 3*u^116 + 3*u^114 - u^108 - 3*u^106 - 2*u^
#                104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^96 - u^94 + u^88 + u^
#                86 + u^84 + u^82, u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^
#                198 - u^194 - 3*u^192 - 4*u^190 - 5*u^188 - 6*u^186 - 5*u^
#                184 - 5*u^182 - 4*u^180 - u^178 + 3*u^174 + 6*u^172 + 7*u^
#                170 + 10*u^168 + 10*u^166 + 9*u^164 + 10*u^162 + 7*u^160 + 
#                5*u^158 + 3*u^156 - 2*u^154 - 3*u^152 - 6*u^150 - 10*u^148 - 
#                10*u^146 - 13*u^144 - 13*u^142 - 10*u^140 - 10*u^138 - 6*u^
#                136 - 3*u^134 - 2*u^132 + 3*u^130 + 5*u^128 + 7*u^126 + 10*u^
#                124 + 9*u^122 + 10*u^120 + 10*u^118 + 7*u^116 + 6*u^114 + 3*u^
#                112 - u^108 - 4*u^106 - 5*u^104 - 5*u^102 - 6*u^100 - 5*u^
#                98 - 4*u^96 - 3*u^94 - u^92 + u^88 + 2*u^86 + 2*u^84 + 2*u^
#                82 + u^80, 2*u^206 + 3*u^204 + 3*u^202 + 3*u^200 + u^198 - 
#                2*u^194 - 5*u^192 - 6*u^190 - 8*u^188 - 9*u^186 - 7*u^184 - 
#                8*u^182 - 5*u^180 - u^178 + 6*u^174 + 9*u^172 + 11*u^170 + 
#                16*u^168 + 14*u^166 + 14*u^164 + 15*u^162 + 9*u^160 + 8*u^
#                158 + 3*u^156 - 4*u^154 - 4*u^152 - 11*u^150 - 15*u^148 - 
#                15*u^146 - 21*u^144 - 18*u^142 - 15*u^140 - 15*u^138 - 7*u^
#                136 - 5*u^134 - 2*u^132 + 6*u^130 + 7*u^128 + 12*u^126 + 15*u^
#                124 + 13*u^122 + 16*u^120 + 14*u^118 + 10*u^116 + 9*u^114 + 
#                3*u^112 - 2*u^108 - 7*u^106 - 7*u^104 - 8*u^102 - 9*u^100 - 
#                7*u^98 - 6*u^96 - 4*u^94 - u^92 + 2*u^88 + 3*u^86 + 3*u^84 + 
#                3*u^82 + u^80, u^208 + u^206 + 3*u^204 + 3*u^202 + 2*u^200 + 
#                2*u^198 - u^196 - 2*u^194 - 4*u^192 - 7*u^190 - 7*u^188 - 8*u^
#                186 - 9*u^184 - 5*u^182 - 6*u^180 - 2*u^178 + 3*u^176 + 3*u^
#                174 + 10*u^172 + 12*u^170 + 12*u^168 + 17*u^166 + 13*u^164 + 
#                12*u^162 + 13*u^160 + 4*u^158 + 4*u^156 - u^154 - 9*u^152 - 
#                7*u^150 - 15*u^148 - 18*u^146 - 15*u^144 - 21*u^142 - 15*u^
#                140 - 11*u^138 - 12*u^136 - 2*u^134 - u^132 + 2*u^130 + 11*u^
#                128 + 9*u^126 + 14*u^124 + 16*u^122 + 12*u^120 + 15*u^118 + 
#                11*u^116 + 6*u^114 + 6*u^112 - u^110 - 3*u^108 - 4*u^106 - 
#                9*u^104 - 7*u^102 - 8*u^100 - 8*u^98 - 5*u^96 - 4*u^94 - 2*u^
#                92 + u^90 + u^88 + 3*u^86 + 3*u^84 + 2*u^82 + 2*u^80, 
#              u^207 + 3*u^205 + 3*u^203 + 3*u^201 + 2*u^199 - u^195 - 4*u^
#                193 - 6*u^191 - 7*u^189 - 9*u^187 - 8*u^185 - 7*u^183 - 7*u^
#                181 - 2*u^179 + 3*u^175 + 9*u^173 + 10*u^171 + 14*u^169 + 
#                16*u^167 + 13*u^165 + 15*u^163 + 12*u^161 + 7*u^159 + 6*u^
#                157 - 2*u^155 - 5*u^153 - 7*u^151 - 15*u^149 - 15*u^147 - 
#                18*u^145 - 21*u^143 - 15*u^141 - 15*u^139 - 11*u^137 - 4*u^
#                135 - 4*u^133 + 3*u^131 + 8*u^129 + 9*u^127 + 15*u^125 + 14*u^
#                123 + 14*u^121 + 16*u^119 + 11*u^117 + 9*u^115 + 6*u^113 - u^
#                109 - 5*u^107 - 8*u^105 - 7*u^103 - 9*u^101 - 8*u^99 - 6*u^
#                97 - 5*u^95 - 2*u^93 + u^89 + 3*u^87 + 3*u^85 + 3*u^83 + 2*u^
#                81, u^205 + u^203 + u^201 + u^199 - u^193 - 2*u^191 - 2*u^
#                189 - 3*u^187 - 3*u^185 - 2*u^183 - 3*u^181 - u^179 + 3*u^
#                173 + 3*u^171 + 4*u^169 + 6*u^167 + 4*u^165 + 5*u^163 + 5*u^
#                161 + 2*u^159 + 3*u^157 - 2*u^153 - u^151 - 5*u^149 - 5*u^
#                147 - 5*u^145 - 8*u^143 - 5*u^141 - 5*u^139 - 5*u^137 - u^
#                135 - 2*u^133 + 3*u^129 + 2*u^127 + 5*u^125 + 5*u^123 + 4*u^
#                121 + 6*u^119 + 4*u^117 + 3*u^115 + 3*u^113 - u^107 - 3*u^
#                105 - 2*u^103 - 3*u^101 - 3*u^99 - 2*u^97 - 2*u^95 - u^93 + u^
#                87 + u^85 + u^83 + u^81, 0*u^0, 0*u^0, 0*u^0, 
#              u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^198 - u^194 - 3*u^192 - 
#                4*u^190 - 5*u^188 - 6*u^186 - 5*u^184 - 5*u^182 - 4*u^180 - u^
#                178 + 3*u^174 + 6*u^172 + 7*u^170 + 10*u^168 + 10*u^166 + 9*u^
#                164 + 10*u^162 + 7*u^160 + 5*u^158 + 3*u^156 - 2*u^154 - 3*u^
#                152 - 6*u^150 - 10*u^148 - 10*u^146 - 13*u^144 - 13*u^142 - 
#                10*u^140 - 10*u^138 - 6*u^136 - 3*u^134 - 2*u^132 + 3*u^130 + 
#                5*u^128 + 7*u^126 + 10*u^124 + 9*u^122 + 10*u^120 + 10*u^
#                118 + 7*u^116 + 6*u^114 + 3*u^112 - u^108 - 4*u^106 - 5*u^
#                104 - 5*u^102 - 6*u^100 - 5*u^98 - 4*u^96 - 3*u^94 - u^92 + u^
#                88 + 2*u^86 + 2*u^84 + 2*u^82 + u^80, 
#              u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^192 - 2*u^190 - 3*u^
#                188 - 3*u^186 - 2*u^184 - 3*u^182 - u^180 + 3*u^174 + 3*u^
#                172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^164 + 5*u^162 + 2*u^
#                160 + 3*u^158 - 2*u^154 - u^152 - 5*u^150 - 5*u^148 - 5*u^
#                146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^138 - u^136 - 2*u^
#                134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^124 + 4*u^122 + 6*u^
#                120 + 4*u^118 + 3*u^116 + 3*u^114 - u^108 - 3*u^106 - 2*u^
#                104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^96 - u^94 + u^88 + u^
#                86 + u^84 + u^82, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^
#                192 - 2*u^190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^
#                180 + 3*u^174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^
#                164 + 5*u^162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^
#                150 - 5*u^148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^
#                138 - u^136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^
#                124 + 4*u^122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^
#                108 - 3*u^106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^
#                96 - u^94 + u^88 + u^86 + u^84 + u^82, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ 2*u^205 + 3*u^203 + 2*u^201 + 2*u^199 - u^195 - 2*u^193 - 5*u^
#                191 - 5*u^189 - 6*u^187 - 7*u^185 - 4*u^183 - 5*u^181 - 3*u^
#                179 + 2*u^177 + u^175 + 6*u^173 + 9*u^171 + 8*u^169 + 13*u^
#                167 + 10*u^165 + 8*u^163 + 11*u^161 + 4*u^159 + 3*u^157 + u^
#                155 - 7*u^153 - 4*u^151 - 9*u^149 - 14*u^147 - 10*u^145 - 
#                16*u^143 - 13*u^141 - 7*u^139 - 10*u^137 - 2*u^135 - u^131 + 
#                8*u^129 + 7*u^127 + 9*u^125 + 13*u^123 + 8*u^121 + 11*u^119 + 
#                10*u^117 + 4*u^115 + 5*u^113 - 3*u^109 - 2*u^107 - 7*u^105 - 
#                6*u^103 - 5*u^101 - 7*u^99 - 4*u^97 - 3*u^95 - 2*u^93 + u^
#                91 + u^89 + 2*u^87 + 3*u^85 + 2*u^83 + 2*u^81 - u^77, 
#              3*u^206 + 5*u^204 + 4*u^202 + 4*u^200 + u^198 - u^196 - 3*u^
#                194 - 8*u^192 - 9*u^190 - 11*u^188 - 13*u^186 - 9*u^184 - 
#                10*u^182 - 7*u^180 + u^178 + u^176 + 9*u^174 + 15*u^172 + 
#                15*u^170 + 23*u^168 + 20*u^166 + 17*u^164 + 21*u^162 + 11*u^
#                160 + 8*u^158 + 4*u^156 - 9*u^154 - 7*u^152 - 15*u^150 - 24*u^
#                148 - 20*u^146 - 29*u^144 - 26*u^142 - 17*u^140 - 20*u^138 - 
#                8*u^136 - 3*u^134 - 3*u^132 + 11*u^130 + 12*u^128 + 16*u^
#                126 + 23*u^124 + 17*u^122 + 21*u^120 + 20*u^118 + 11*u^116 + 
#                11*u^114 + 3*u^112 - 3*u^110 - 3*u^108 - 11*u^106 - 11*u^
#                104 - 10*u^102 - 13*u^100 - 9*u^98 - 7*u^96 - 5*u^94 + u^90 + 
#                3*u^88 + 5*u^86 + 4*u^84 + 4*u^82 + u^80 - u^78, 
#              u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^195 - 3*u^193 - 3*u^
#                191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^183 - 2*u^181 + 2*u^
#                179 + u^177 + 3*u^175 + 6*u^173 + 4*u^171 + 7*u^169 + 6*u^
#                167 + 3*u^165 + 6*u^163 + 2*u^161 + u^157 - 5*u^155 - 3*u^
#                153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^145 - 8*u^143 - 2*u^
#                141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 5*u^131 + 5*u^129 + 
#                4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 6*u^119 + u^117 + 2*u^
#                115 - 3*u^111 - u^109 - 4*u^107 - 4*u^105 - 2*u^103 - 4*u^
#                101 - 2*u^99 - u^97 - u^95 + u^93 + u^91 + u^89 + 2*u^87 + u^
#                85 + u^83 - u^79, u^207 + 3*u^205 + 2*u^203 + 2*u^201 + u^
#                199 - u^197 - u^195 - 4*u^193 - 5*u^191 - 5*u^189 - 7*u^187 - 
#                5*u^185 - 4*u^183 - 5*u^181 + u^179 + u^177 + 3*u^175 + 9*u^
#                173 + 7*u^171 + 11*u^169 + 12*u^167 + 7*u^165 + 11*u^163 + 
#                7*u^161 + 2*u^159 + 4*u^157 - 5*u^155 - 5*u^153 - 5*u^151 - 
#                14*u^149 - 10*u^147 - 13*u^145 - 16*u^143 - 7*u^141 - 10*u^
#                139 - 6*u^137 + u^135 - 3*u^133 + 5*u^131 + 8*u^129 + 6*u^
#                127 + 13*u^125 + 9*u^123 + 9*u^121 + 12*u^119 + 5*u^117 + 5*u^
#                115 + 3*u^113 - 3*u^111 - u^109 - 5*u^107 - 7*u^105 - 4*u^
#                103 - 7*u^101 - 5*u^99 - 3*u^97 - 3*u^95 + u^91 + u^89 + 3*u^
#                87 + 2*u^85 + 2*u^83 + u^81 - u^79, 
#              u^207 + 4*u^205 + 4*u^203 + 3*u^201 + 2*u^199 - u^197 - 2*u^
#                195 - 5*u^193 - 8*u^191 - 8*u^189 - 10*u^187 - 9*u^185 - 6*u^
#                183 - 7*u^181 - u^179 + 3*u^177 + 4*u^175 + 12*u^173 + 13*u^
#                171 + 15*u^169 + 19*u^167 + 13*u^165 + 14*u^163 + 13*u^161 + 
#                4*u^159 + 4*u^157 - 4*u^155 - 10*u^153 - 8*u^151 - 18*u^149 - 
#                19*u^147 - 18*u^145 - 24*u^143 - 15*u^141 - 12*u^139 - 11*u^
#                137 - u^133 + 4*u^131 + 13*u^129 + 11*u^127 + 17*u^125 + 17*u^
#                123 + 13*u^121 + 17*u^119 + 11*u^117 + 6*u^115 + 5*u^113 - 
#                3*u^111 - 4*u^109 - 6*u^107 - 11*u^105 - 8*u^103 - 9*u^101 - 
#                9*u^99 - 5*u^97 - 4*u^95 - u^93 + 2*u^91 + 2*u^89 + 4*u^87 + 
#                4*u^85 + 3*u^83 + 2*u^81 - u^79 - u^77, 
#              u^207 + 3*u^205 + 3*u^203 + 3*u^201 + 2*u^199 - u^195 - 4*u^
#                193 - 6*u^191 - 7*u^189 - 9*u^187 - 8*u^185 - 7*u^183 - 7*u^
#                181 - 2*u^179 + 3*u^175 + 9*u^173 + 10*u^171 + 14*u^169 + 
#                16*u^167 + 13*u^165 + 15*u^163 + 12*u^161 + 7*u^159 + 6*u^
#                157 - 2*u^155 - 5*u^153 - 7*u^151 - 15*u^149 - 15*u^147 - 
#                18*u^145 - 21*u^143 - 15*u^141 - 15*u^139 - 11*u^137 - 4*u^
#                135 - 4*u^133 + 3*u^131 + 8*u^129 + 9*u^127 + 15*u^125 + 14*u^
#                123 + 14*u^121 + 16*u^119 + 11*u^117 + 9*u^115 + 6*u^113 - u^
#                109 - 5*u^107 - 8*u^105 - 7*u^103 - 9*u^101 - 8*u^99 - 6*u^
#                97 - 5*u^95 - 2*u^93 + u^89 + 3*u^87 + 3*u^85 + 3*u^83 + 2*u^
#                81, u^208 + 4*u^206 + 4*u^204 + 3*u^202 + 2*u^200 - u^198 - 
#                2*u^196 - 5*u^194 - 8*u^192 - 8*u^190 - 10*u^188 - 9*u^186 - 
#                6*u^184 - 7*u^182 - u^180 + 3*u^178 + 4*u^176 + 12*u^174 + 
#                13*u^172 + 15*u^170 + 19*u^168 + 13*u^166 + 14*u^164 + 13*u^
#                162 + 4*u^160 + 4*u^158 - 4*u^156 - 10*u^154 - 8*u^152 - 18*u^
#                150 - 19*u^148 - 18*u^146 - 24*u^144 - 15*u^142 - 12*u^140 - 
#                11*u^138 - u^134 + 4*u^132 + 13*u^130 + 11*u^128 + 17*u^126 + 
#                17*u^124 + 13*u^122 + 17*u^120 + 11*u^118 + 6*u^116 + 5*u^
#                114 - 3*u^112 - 4*u^110 - 6*u^108 - 11*u^106 - 8*u^104 - 9*u^
#                102 - 9*u^100 - 5*u^98 - 4*u^96 - u^94 + 2*u^92 + 2*u^90 + 
#                4*u^88 + 4*u^86 + 3*u^84 + 2*u^82 - u^80 - u^78, 
#              u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^192 - 2*u^190 - 3*u^
#                188 - 3*u^186 - 2*u^184 - 3*u^182 - u^180 + 3*u^174 + 3*u^
#                172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^164 + 5*u^162 + 2*u^
#                160 + 3*u^158 - 2*u^154 - u^152 - 5*u^150 - 5*u^148 - 5*u^
#                146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^138 - u^136 - 2*u^
#                134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^124 + 4*u^122 + 6*u^
#                120 + 4*u^118 + 3*u^116 + 3*u^114 - u^108 - 3*u^106 - 2*u^
#                104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^96 - u^94 + u^88 + u^
#                86 + u^84 + u^82, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 
#              u^207 + u^205 - u^199 - u^197 - u^195 - 2*u^193 - u^191 - u^
#                189 - u^187 + u^185 + u^181 + 3*u^179 + u^177 + 3*u^175 + 3*u^
#                173 + u^171 + 3*u^169 - u^165 + u^163 - 3*u^161 - 2*u^159 - 
#                2*u^157 - 5*u^155 - u^153 - 3*u^151 - 4*u^149 - 3*u^145 + 3*u^
#                141 + 4*u^137 + 3*u^135 + u^133 + 5*u^131 + 2*u^129 + 2*u^
#                127 + 3*u^125 - u^123 + u^121 - 3*u^117 - u^115 - 3*u^113 - 
#                3*u^111 - u^109 - 3*u^107 - u^105 - u^101 + u^99 + u^97 + u^
#                95 + 2*u^93 + u^91 + u^89 + u^87 - u^81 - u^79, 0*u^0, 
#              u^207 + 2*u^205 + 2*u^203 + 2*u^201 + u^199 - u^195 - 3*u^193 - 
#                4*u^191 - 5*u^189 - 6*u^187 - 5*u^185 - 5*u^183 - 4*u^181 - u^
#                179 + 3*u^175 + 6*u^173 + 7*u^171 + 10*u^169 + 10*u^167 + 9*u^
#                165 + 10*u^163 + 7*u^161 + 5*u^159 + 3*u^157 - 2*u^155 - 3*u^
#                153 - 6*u^151 - 10*u^149 - 10*u^147 - 13*u^145 - 13*u^143 - 
#                10*u^141 - 10*u^139 - 6*u^137 - 3*u^135 - 2*u^133 + 3*u^131 + 
#                5*u^129 + 7*u^127 + 10*u^125 + 9*u^123 + 10*u^121 + 10*u^
#                119 + 7*u^117 + 6*u^115 + 3*u^113 - u^109 - 4*u^107 - 5*u^
#                105 - 5*u^103 - 6*u^101 - 5*u^99 - 4*u^97 - 3*u^95 - u^93 + u^
#                89 + 2*u^87 + 2*u^85 + 2*u^83 + u^81, 
#              u^207 + u^205 + u^203 + u^201 - u^195 - 2*u^193 - 2*u^191 - 3*u^
#                189 - 3*u^187 - 2*u^185 - 3*u^183 - u^181 + 3*u^175 + 3*u^
#                173 + 4*u^171 + 6*u^169 + 4*u^167 + 5*u^165 + 5*u^163 + 2*u^
#                161 + 3*u^159 - 2*u^155 - u^153 - 5*u^151 - 5*u^149 - 5*u^
#                147 - 8*u^145 - 5*u^143 - 5*u^141 - 5*u^139 - u^137 - 2*u^
#                135 + 3*u^131 + 2*u^129 + 5*u^127 + 5*u^125 + 4*u^123 + 6*u^
#                121 + 4*u^119 + 3*u^117 + 3*u^115 - u^109 - 3*u^107 - 2*u^
#                105 - 3*u^103 - 3*u^101 - 2*u^99 - 2*u^97 - u^95 + u^89 + u^
#                87 + u^85 + u^83, u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^
#                195 - 3*u^193 - 3*u^191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^
#                183 - 2*u^181 + 2*u^179 + u^177 + 3*u^175 + 6*u^173 + 4*u^
#                171 + 7*u^169 + 6*u^167 + 3*u^165 + 6*u^163 + 2*u^161 + u^
#                157 - 5*u^155 - 3*u^153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^
#                145 - 8*u^143 - 2*u^141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 
#                5*u^131 + 5*u^129 + 4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 
#                6*u^119 + u^117 + 2*u^115 - 3*u^111 - u^109 - 4*u^107 - 4*u^
#                105 - 2*u^103 - 4*u^101 - 2*u^99 - u^97 - u^95 + u^93 + u^
#                91 + u^89 + 2*u^87 + u^85 + u^83 - u^79, 0*u^0, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, 0*u^0 ], 
#          [ u^205 + 2*u^203 + u^201 + u^199 - u^195 - u^193 - 3*u^191 - 3*u^
#                189 - 3*u^187 - 4*u^185 - 2*u^183 - 2*u^181 - 2*u^179 + 2*u^
#                177 + u^175 + 3*u^173 + 6*u^171 + 4*u^169 + 7*u^167 + 6*u^
#                165 + 3*u^163 + 6*u^161 + 2*u^159 + u^155 - 5*u^153 - 3*u^
#                151 - 4*u^149 - 9*u^147 - 5*u^145 - 8*u^143 - 8*u^141 - 2*u^
#                139 - 5*u^137 - u^135 + 2*u^133 - u^131 + 5*u^129 + 5*u^127 + 
#                4*u^125 + 8*u^123 + 4*u^121 + 5*u^119 + 6*u^117 + u^115 + 2*u^
#                113 - 3*u^109 - u^107 - 4*u^105 - 4*u^103 - 2*u^101 - 4*u^
#                99 - 2*u^97 - u^95 - u^93 + u^91 + u^89 + u^87 + 2*u^85 + u^
#                83 + u^81 - u^77, 2*u^206 + 3*u^204 + 2*u^202 + 2*u^200 - u^
#                196 - 2*u^194 - 5*u^192 - 5*u^190 - 6*u^188 - 7*u^186 - 4*u^
#                184 - 5*u^182 - 3*u^180 + 2*u^178 + u^176 + 6*u^174 + 9*u^
#                172 + 8*u^170 + 13*u^168 + 10*u^166 + 8*u^164 + 11*u^162 + 
#                4*u^160 + 3*u^158 + u^156 - 7*u^154 - 4*u^152 - 9*u^150 - 
#                14*u^148 - 10*u^146 - 16*u^144 - 13*u^142 - 7*u^140 - 10*u^
#                138 - 2*u^136 - u^132 + 8*u^130 + 7*u^128 + 9*u^126 + 13*u^
#                124 + 8*u^122 + 11*u^120 + 10*u^118 + 4*u^116 + 5*u^114 - 3*u^
#                110 - 2*u^108 - 7*u^106 - 6*u^104 - 5*u^102 - 7*u^100 - 4*u^
#                98 - 3*u^96 - 2*u^94 + u^92 + u^90 + 2*u^88 + 3*u^86 + 2*u^
#                84 + 2*u^82 - u^78, u^207 + 2*u^205 + u^203 + u^201 - u^
#                197 - u^195 - 3*u^193 - 3*u^191 - 3*u^189 - 4*u^187 - 2*u^
#                185 - 2*u^183 - 2*u^181 + 2*u^179 + u^177 + 3*u^175 + 6*u^
#                173 + 4*u^171 + 7*u^169 + 6*u^167 + 3*u^165 + 6*u^163 + 2*u^
#                161 + u^157 - 5*u^155 - 3*u^153 - 4*u^151 - 9*u^149 - 5*u^
#                147 - 8*u^145 - 8*u^143 - 2*u^141 - 5*u^139 - u^137 + 2*u^
#                135 - u^133 + 5*u^131 + 5*u^129 + 4*u^127 + 8*u^125 + 4*u^
#                123 + 5*u^121 + 6*u^119 + u^117 + 2*u^115 - 3*u^111 - u^109 - 
#                4*u^107 - 4*u^105 - 2*u^103 - 4*u^101 - 2*u^99 - u^97 - u^
#                95 + u^93 + u^91 + u^89 + 2*u^87 + u^85 + u^83 - u^79, 
#              u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^195 - 3*u^193 - 3*u^
#                191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^183 - 2*u^181 + 2*u^
#                179 + u^177 + 3*u^175 + 6*u^173 + 4*u^171 + 7*u^169 + 6*u^
#                167 + 3*u^165 + 6*u^163 + 2*u^161 + u^157 - 5*u^155 - 3*u^
#                153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^145 - 8*u^143 - 2*u^
#                141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 5*u^131 + 5*u^129 + 
#                4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 6*u^119 + u^117 + 2*u^
#                115 - 3*u^111 - u^109 - 4*u^107 - 4*u^105 - 2*u^103 - 4*u^
#                101 - 2*u^99 - u^97 - u^95 + u^93 + u^91 + u^89 + 2*u^87 + u^
#                85 + u^83 - u^79, u^205 + u^203 + u^201 + u^199 - u^193 - 2*u^
#                191 - 2*u^189 - 3*u^187 - 3*u^185 - 2*u^183 - 3*u^181 - u^
#                179 + 3*u^173 + 3*u^171 + 4*u^169 + 6*u^167 + 4*u^165 + 5*u^
#                163 + 5*u^161 + 2*u^159 + 3*u^157 - 2*u^153 - u^151 - 5*u^
#                149 - 5*u^147 - 5*u^145 - 8*u^143 - 5*u^141 - 5*u^139 - 5*u^
#                137 - u^135 - 2*u^133 + 3*u^129 + 2*u^127 + 5*u^125 + 5*u^
#                123 + 4*u^121 + 6*u^119 + 4*u^117 + 3*u^115 + 3*u^113 - u^
#                107 - 3*u^105 - 2*u^103 - 3*u^101 - 3*u^99 - 2*u^97 - 2*u^
#                95 - u^93 + u^87 + u^85 + u^83 + u^81, 
#              u^205 + u^203 + u^201 + u^199 - u^193 - 2*u^191 - 2*u^189 - 3*u^
#                187 - 3*u^185 - 2*u^183 - 3*u^181 - u^179 + 3*u^173 + 3*u^
#                171 + 4*u^169 + 6*u^167 + 4*u^165 + 5*u^163 + 5*u^161 + 2*u^
#                159 + 3*u^157 - 2*u^153 - u^151 - 5*u^149 - 5*u^147 - 5*u^
#                145 - 8*u^143 - 5*u^141 - 5*u^139 - 5*u^137 - u^135 - 2*u^
#                133 + 3*u^129 + 2*u^127 + 5*u^125 + 5*u^123 + 4*u^121 + 6*u^
#                119 + 4*u^117 + 3*u^115 + 3*u^113 - u^107 - 3*u^105 - 2*u^
#                103 - 3*u^101 - 3*u^99 - 2*u^97 - 2*u^95 - u^93 + u^87 + u^
#                85 + u^83 + u^81, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^
#                192 - 2*u^190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^
#                180 + 3*u^174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^
#                164 + 5*u^162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^
#                150 - 5*u^148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^
#                138 - u^136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^
#                124 + 4*u^122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^
#                108 - 3*u^106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^
#                96 - u^94 + u^88 + u^86 + u^84 + u^82, 
#              u^208 + 2*u^206 + u^204 + u^202 - u^198 - u^196 - 3*u^194 - 3*u^
#                192 - 3*u^190 - 4*u^188 - 2*u^186 - 2*u^184 - 2*u^182 + 2*u^
#                180 + u^178 + 3*u^176 + 6*u^174 + 4*u^172 + 7*u^170 + 6*u^
#                168 + 3*u^166 + 6*u^164 + 2*u^162 + u^158 - 5*u^156 - 3*u^
#                154 - 4*u^152 - 9*u^150 - 5*u^148 - 8*u^146 - 8*u^144 - 2*u^
#                142 - 5*u^140 - u^138 + 2*u^136 - u^134 + 5*u^132 + 5*u^130 + 
#                4*u^128 + 8*u^126 + 4*u^124 + 5*u^122 + 6*u^120 + u^118 + 2*u^
#                116 - 3*u^112 - u^110 - 4*u^108 - 4*u^106 - 2*u^104 - 4*u^
#                102 - 2*u^100 - u^98 - u^96 + u^94 + u^92 + u^90 + 2*u^88 + u^
#                86 + u^84 - u^80, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 0*u^0, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81, u^207 + u^205 + u^203 + u^201 - u^195 - 2*u^193 - 2*u^
#                191 - 3*u^189 - 3*u^187 - 2*u^185 - 3*u^183 - u^181 + 3*u^
#                175 + 3*u^173 + 4*u^171 + 6*u^169 + 4*u^167 + 5*u^165 + 5*u^
#                163 + 2*u^161 + 3*u^159 - 2*u^155 - u^153 - 5*u^151 - 5*u^
#                149 - 5*u^147 - 8*u^145 - 5*u^143 - 5*u^141 - 5*u^139 - u^
#                137 - 2*u^135 + 3*u^131 + 2*u^129 + 5*u^127 + 5*u^125 + 4*u^
#                123 + 6*u^121 + 4*u^119 + 3*u^117 + 3*u^115 - u^109 - 3*u^
#                107 - 2*u^105 - 3*u^103 - 3*u^101 - 2*u^99 - 2*u^97 - u^
#                95 + u^89 + u^87 + u^85 + u^83, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#              0*u^0 ], 
#          [ u^205 + u^203 - u^197 - u^195 - u^193 - 2*u^191 - u^189 - u^
#                187 - u^185 + u^183 + u^179 + 3*u^177 + u^175 + 3*u^173 + 3*u^
#                171 + u^169 + 3*u^167 - u^163 + u^161 - 3*u^159 - 2*u^157 - 
#                2*u^155 - 5*u^153 - u^151 - 3*u^149 - 4*u^147 - 3*u^143 + 3*u^
#                139 + 4*u^135 + 3*u^133 + u^131 + 5*u^129 + 2*u^127 + 2*u^
#                125 + 3*u^123 - u^121 + u^119 - 3*u^115 - u^113 - 3*u^111 - 
#                3*u^109 - u^107 - 3*u^105 - u^103 - u^99 + u^97 + u^95 + u^
#                93 + 2*u^91 + u^89 + u^87 + u^85 - u^79 - u^77, 
#              u^206 + u^204 - u^198 - u^196 - u^194 - 2*u^192 - u^190 - u^
#                188 - u^186 + u^184 + u^180 + 3*u^178 + u^176 + 3*u^174 + 3*u^
#                172 + u^170 + 3*u^168 - u^164 + u^162 - 3*u^160 - 2*u^158 - 
#                2*u^156 - 5*u^154 - u^152 - 3*u^150 - 4*u^148 - 3*u^144 + 3*u^
#                140 + 4*u^136 + 3*u^134 + u^132 + 5*u^130 + 2*u^128 + 2*u^
#                126 + 3*u^124 - u^122 + u^120 - 3*u^116 - u^114 - 3*u^112 - 
#                3*u^110 - u^108 - 3*u^106 - u^104 - u^100 + u^98 + u^96 + u^
#                94 + 2*u^92 + u^90 + u^88 + u^86 - u^80 - u^78, 
#              u^207 + u^205 - u^199 - u^197 - u^195 - 2*u^193 - u^191 - u^
#                189 - u^187 + u^185 + u^181 + 3*u^179 + u^177 + 3*u^175 + 3*u^
#                173 + u^171 + 3*u^169 - u^165 + u^163 - 3*u^161 - 2*u^159 - 
#                2*u^157 - 5*u^155 - u^153 - 3*u^151 - 4*u^149 - 3*u^145 + 3*u^
#                141 + 4*u^137 + 3*u^135 + u^133 + 5*u^131 + 2*u^129 + 2*u^
#                127 + 3*u^125 - u^123 + u^121 - 3*u^117 - u^115 - 3*u^113 - 
#                3*u^111 - u^109 - 3*u^107 - u^105 - u^101 + u^99 + u^97 + u^
#                95 + 2*u^93 + u^91 + u^89 + u^87 - u^81 - u^79, 
#              u^207 + 2*u^205 - u^199 - 2*u^197 - u^195 - 3*u^193 - 2*u^
#                191 - u^189 - 2*u^187 + u^185 + u^183 + 5*u^179 + 2*u^177 + 
#                3*u^175 + 6*u^173 + u^171 + 4*u^169 + 2*u^167 - 3*u^165 + 2*u^
#                163 - 3*u^161 - 5*u^159 - u^157 - 8*u^155 - 3*u^153 - 2*u^
#                151 - 8*u^149 - 3*u^145 - 3*u^143 + 6*u^141 + 4*u^137 + 7*u^
#                135 + 7*u^131 + 5*u^129 + u^127 + 6*u^125 - u^123 + 2*u^119 - 
#                5*u^117 - 2*u^115 - 3*u^113 - 6*u^111 - u^109 - 4*u^107 - 3*u^
#                105 + u^103 - 2*u^101 + u^99 + 2*u^97 + u^95 + 3*u^93 + 2*u^
#                91 + u^89 + 2*u^87 - u^81 - 2*u^79, 
#              u^205 - u^197 - u^193 - u^191 - u^187 + u^183 - u^181 + 2*u^
#                179 + u^177 + 3*u^173 + u^169 + 2*u^167 - 2*u^165 + u^163 - 
#                3*u^159 + u^157 - 3*u^155 - 2*u^153 + u^151 - 4*u^149 - 3*u^
#                143 + 3*u^141 + 4*u^135 - u^133 + 2*u^131 + 3*u^129 - u^127 + 
#                3*u^125 - u^121 + 2*u^119 - 2*u^117 - u^115 - 3*u^111 - u^
#                107 - 2*u^105 + u^103 - u^101 + u^97 + u^93 + u^91 + u^87 - u^
#                79, 0*u^0, u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^
#                182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^
#                166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 
#                4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^
#                130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 
#                3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^
#                92 + u^88 - u^80, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 
#              u^208 + 2*u^206 - u^200 - 2*u^198 - u^196 - 3*u^194 - 2*u^
#                192 - u^190 - 2*u^188 + u^186 + u^184 + 5*u^180 + 2*u^178 + 
#                3*u^176 + 6*u^174 + u^172 + 4*u^170 + 2*u^168 - 3*u^166 + 2*u^
#                164 - 3*u^162 - 5*u^160 - u^158 - 8*u^156 - 3*u^154 - 2*u^
#                152 - 8*u^150 - 3*u^146 - 3*u^144 + 6*u^142 + 4*u^138 + 7*u^
#                136 + 7*u^132 + 5*u^130 + u^128 + 6*u^126 - u^124 + 2*u^120 - 
#                5*u^118 - 2*u^116 - 3*u^114 - 6*u^112 - u^110 - 4*u^108 - 3*u^
#                106 + u^104 - 2*u^102 + u^100 + 2*u^98 + u^96 + 3*u^94 + 2*u^
#                92 + u^90 + 2*u^88 - u^82 - 2*u^80, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81, u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 
#                2*u^181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^
#                165 - 3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 
#                3*u^145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^
#                129 + 3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^
#                113 - u^109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^
#                93 + u^89 - u^81, 0*u^0, 0*u^0, 0*u^0, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81, 0*u^0, 0*u^0 ], 
#          [ u^204 - u^196 - u^192 - u^190 - u^186 + u^182 - u^180 + 2*u^
#                178 + u^176 + 3*u^172 + u^168 + 2*u^166 - 2*u^164 + u^162 - 
#                3*u^158 + u^156 - 3*u^154 - 2*u^152 + u^150 - 4*u^148 - 3*u^
#                142 + 3*u^140 + 4*u^134 - u^132 + 2*u^130 + 3*u^128 - u^126 + 
#                3*u^124 - u^120 + 2*u^118 - 2*u^116 - u^114 - 3*u^110 - u^
#                106 - 2*u^104 + u^102 - u^100 + u^96 + u^92 + u^90 + u^86 - u^
#                78, u^205 - u^197 - u^193 - u^191 - u^187 + u^183 - u^181 + 
#                2*u^179 + u^177 + 3*u^173 + u^169 + 2*u^167 - 2*u^165 + u^
#                163 - 3*u^159 + u^157 - 3*u^155 - 2*u^153 + u^151 - 4*u^149 - 
#                3*u^143 + 3*u^141 + 4*u^135 - u^133 + 2*u^131 + 3*u^129 - u^
#                127 + 3*u^125 - u^121 + 2*u^119 - 2*u^117 - u^115 - 3*u^
#                111 - u^107 - 2*u^105 + u^103 - u^101 + u^97 + u^93 + u^
#                91 + u^87 - u^79, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, u^206 + u^204 - u^198 - u^196 - u^194 - 2*u^192 - u^
#                190 - u^188 - u^186 + u^184 + u^180 + 3*u^178 + u^176 + 3*u^
#                174 + 3*u^172 + u^170 + 3*u^168 - u^164 + u^162 - 3*u^160 - 
#                2*u^158 - 2*u^156 - 5*u^154 - u^152 - 3*u^150 - 4*u^148 - 3*u^
#                144 + 3*u^140 + 4*u^136 + 3*u^134 + u^132 + 5*u^130 + 2*u^
#                128 + 2*u^126 + 3*u^124 - u^122 + u^120 - 3*u^116 - u^114 - 
#                3*u^112 - 3*u^110 - u^108 - 3*u^106 - u^104 - u^100 + u^
#                98 + u^96 + u^94 + 2*u^92 + u^90 + u^88 + u^86 - u^80 - u^78, 
#              0*u^0, u^207 + u^205 - u^199 - u^197 - u^195 - 2*u^193 - u^
#                191 - u^189 - u^187 + u^185 + u^181 + 3*u^179 + u^177 + 3*u^
#                175 + 3*u^173 + u^171 + 3*u^169 - u^165 + u^163 - 3*u^161 - 
#                2*u^159 - 2*u^157 - 5*u^155 - u^153 - 3*u^151 - 4*u^149 - 3*u^
#                145 + 3*u^141 + 4*u^137 + 3*u^135 + u^133 + 5*u^131 + 2*u^
#                129 + 2*u^127 + 3*u^125 - u^123 + u^121 - 3*u^117 - u^115 - 
#                3*u^113 - 3*u^111 - u^109 - 3*u^107 - u^105 - u^101 + u^
#                99 + u^97 + u^95 + 2*u^93 + u^91 + u^89 + u^87 - u^81 - u^79, 
#              0*u^0, u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 
#                2*u^181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^
#                165 - 3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 
#                3*u^145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^
#                129 + 3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^
#                113 - u^109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^
#                93 + u^89 - u^81, u^208 + u^206 - u^200 - u^198 - u^196 - 2*u^
#                194 - u^192 - u^190 - u^188 + u^186 + u^182 + 3*u^180 + u^
#                178 + 3*u^176 + 3*u^174 + u^172 + 3*u^170 - u^166 + u^164 - 
#                3*u^162 - 2*u^160 - 2*u^158 - 5*u^156 - u^154 - 3*u^152 - 4*u^
#                150 - 3*u^146 + 3*u^142 + 4*u^138 + 3*u^136 + u^134 + 5*u^
#                132 + 2*u^130 + 2*u^128 + 3*u^126 - u^124 + u^122 - 3*u^
#                118 - u^116 - 3*u^114 - 3*u^112 - u^110 - 3*u^108 - u^106 - u^
#                102 + u^100 + u^98 + u^96 + 2*u^94 + u^92 + u^90 + u^88 - u^
#                82 - u^80, 0*u^0, 0*u^0, 0*u^0, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, 0*u^0, u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^
#                183 + 2*u^181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^
#                167 + u^165 - 3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 
#                4*u^151 - 3*u^145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^
#                131 - u^129 + 3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 
#                3*u^113 - u^109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^
#                93 + u^89 - u^81, 0*u^0 ], 
#          [ u^204 - u^196 - u^192 - u^190 - u^186 + u^182 - u^180 + 2*u^
#                178 + u^176 + 3*u^172 + u^168 + 2*u^166 - 2*u^164 + u^162 - 
#                3*u^158 + u^156 - 3*u^154 - 2*u^152 + u^150 - 4*u^148 - 3*u^
#                142 + 3*u^140 + 4*u^134 - u^132 + 2*u^130 + 3*u^128 - u^126 + 
#                3*u^124 - u^120 + 2*u^118 - 2*u^116 - u^114 - 3*u^110 - u^
#                106 - 2*u^104 + u^102 - u^100 + u^96 + u^92 + u^90 + u^86 - u^
#                78, u^205 - u^197 - u^193 - u^191 - u^187 + u^183 - u^181 + 
#                2*u^179 + u^177 + 3*u^173 + u^169 + 2*u^167 - 2*u^165 + u^
#                163 - 3*u^159 + u^157 - 3*u^155 - 2*u^153 + u^151 - 4*u^149 - 
#                3*u^143 + 3*u^141 + 4*u^135 - u^133 + 2*u^131 + 3*u^129 - u^
#                127 + 3*u^125 - u^121 + 2*u^119 - 2*u^117 - u^115 - 3*u^
#                111 - u^107 - 2*u^105 + u^103 - u^101 + u^97 + u^93 + u^
#                91 + u^87 - u^79, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, 0*u^0, 0*u^0, 0*u^0, u^207 - u^199 - u^195 - u^193 - u^
#                189 + u^185 - u^183 + 2*u^181 + u^179 + 3*u^175 + u^171 + 2*u^
#                169 - 2*u^167 + u^165 - 3*u^161 + u^159 - 3*u^157 - 2*u^
#                155 + u^153 - 4*u^151 - 3*u^145 + 3*u^143 + 4*u^137 - u^135 + 
#                2*u^133 + 3*u^131 - u^129 + 3*u^127 - u^123 + 2*u^121 - 2*u^
#                119 - u^117 - 3*u^113 - u^109 - 2*u^107 + u^105 - u^103 + u^
#                99 + u^95 + u^93 + u^89 - u^81, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81, 0*u^0, u^208 - u^200 - u^196 - u^194 - u^190 + u^186 - u^
#                184 + 2*u^182 + u^180 + 3*u^176 + u^172 + 2*u^170 - 2*u^
#                168 + u^166 - 3*u^162 + u^160 - 3*u^158 - 2*u^156 + u^154 - 
#                4*u^152 - 3*u^146 + 3*u^144 + 4*u^138 - u^136 + 2*u^134 + 3*u^
#                132 - u^130 + 3*u^128 - u^124 + 2*u^122 - 2*u^120 - u^118 - 
#                3*u^114 - u^110 - 2*u^108 + u^106 - u^104 + u^100 + u^96 + u^
#                94 + u^90 - u^82, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ 2*u^204 + u^202 + 2*u^200 + u^198 - 2*u^192 - 3*u^190 - 3*u^188 - 
#                5*u^186 - 4*u^184 - 3*u^182 - 5*u^180 - u^176 + u^174 + 5*u^
#                172 + 4*u^170 + 7*u^168 + 9*u^166 + 5*u^164 + 9*u^162 + 6*u^
#                160 + 3*u^158 + 5*u^156 - 2*u^154 - 2*u^152 - 2*u^150 - 9*u^
#                148 - 6*u^146 - 9*u^144 - 12*u^142 - 6*u^140 - 9*u^138 - 6*u^
#                136 - u^134 - 4*u^132 + 2*u^130 + 4*u^128 + 3*u^126 + 9*u^
#                124 + 6*u^122 + 7*u^120 + 9*u^118 + 5*u^116 + 5*u^114 + 4*u^
#                112 - u^110 + u^108 - 3*u^106 - 4*u^104 - 3*u^102 - 5*u^100 - 
#                4*u^98 - 3*u^96 - 3*u^94 - u^92 + 2*u^86 + u^84 + 2*u^82 + u^
#                80, u^207 + 3*u^205 + 3*u^203 + 3*u^201 + 2*u^199 - u^195 - 
#                4*u^193 - 6*u^191 - 7*u^189 - 9*u^187 - 8*u^185 - 7*u^183 - 
#                7*u^181 - 2*u^179 + 3*u^175 + 9*u^173 + 10*u^171 + 14*u^169 + 
#                16*u^167 + 13*u^165 + 15*u^163 + 12*u^161 + 7*u^159 + 6*u^
#                157 - 2*u^155 - 5*u^153 - 7*u^151 - 15*u^149 - 15*u^147 - 
#                18*u^145 - 21*u^143 - 15*u^141 - 15*u^139 - 11*u^137 - 4*u^
#                135 - 4*u^133 + 3*u^131 + 8*u^129 + 9*u^127 + 15*u^125 + 14*u^
#                123 + 14*u^121 + 16*u^119 + 11*u^117 + 9*u^115 + 6*u^113 - u^
#                109 - 5*u^107 - 8*u^105 - 7*u^103 - 9*u^101 - 8*u^99 - 6*u^
#                97 - 5*u^95 - 2*u^93 + u^89 + 3*u^87 + 3*u^85 + 3*u^83 + 2*u^
#                81, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^192 - 2*u^
#                190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^180 + 3*u^
#                174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^164 + 5*u^
#                162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^150 - 5*u^
#                148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^138 - u^
#                136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^124 + 4*u^
#                122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^108 - 3*u^
#                106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^96 - u^
#                94 + u^88 + u^86 + u^84 + u^82, 
#              2*u^206 + u^204 + 2*u^202 + u^200 - 2*u^194 - 3*u^192 - 3*u^
#                190 - 5*u^188 - 4*u^186 - 3*u^184 - 5*u^182 - u^178 + u^176 + 
#                5*u^174 + 4*u^172 + 7*u^170 + 9*u^168 + 5*u^166 + 9*u^164 + 
#                6*u^162 + 3*u^160 + 5*u^158 - 2*u^156 - 2*u^154 - 2*u^152 - 
#                9*u^150 - 6*u^148 - 9*u^146 - 12*u^144 - 6*u^142 - 9*u^140 - 
#                6*u^138 - u^136 - 4*u^134 + 2*u^132 + 4*u^130 + 3*u^128 + 9*u^
#                126 + 6*u^124 + 7*u^122 + 9*u^120 + 5*u^118 + 5*u^116 + 4*u^
#                114 - u^112 + u^110 - 3*u^108 - 4*u^106 - 3*u^104 - 5*u^102 - 
#                4*u^100 - 3*u^98 - 3*u^96 - u^94 + 2*u^88 + u^86 + 2*u^84 + u^
#                82, u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^198 - u^194 - 3*u^
#                192 - 4*u^190 - 5*u^188 - 6*u^186 - 5*u^184 - 5*u^182 - 4*u^
#                180 - u^178 + 3*u^174 + 6*u^172 + 7*u^170 + 10*u^168 + 10*u^
#                166 + 9*u^164 + 10*u^162 + 7*u^160 + 5*u^158 + 3*u^156 - 2*u^
#                154 - 3*u^152 - 6*u^150 - 10*u^148 - 10*u^146 - 13*u^144 - 
#                13*u^142 - 10*u^140 - 10*u^138 - 6*u^136 - 3*u^134 - 2*u^
#                132 + 3*u^130 + 5*u^128 + 7*u^126 + 10*u^124 + 9*u^122 + 10*u^
#                120 + 10*u^118 + 7*u^116 + 6*u^114 + 3*u^112 - u^108 - 4*u^
#                106 - 5*u^104 - 5*u^102 - 6*u^100 - 5*u^98 - 4*u^96 - 3*u^
#                94 - u^92 + u^88 + 2*u^86 + 2*u^84 + 2*u^82 + u^80, 
#              u^206 + 2*u^204 + 2*u^202 + 2*u^200 + u^198 - u^194 - 3*u^192 - 
#                4*u^190 - 5*u^188 - 6*u^186 - 5*u^184 - 5*u^182 - 4*u^180 - u^
#                178 + 3*u^174 + 6*u^172 + 7*u^170 + 10*u^168 + 10*u^166 + 9*u^
#                164 + 10*u^162 + 7*u^160 + 5*u^158 + 3*u^156 - 2*u^154 - 3*u^
#                152 - 6*u^150 - 10*u^148 - 10*u^146 - 13*u^144 - 13*u^142 - 
#                10*u^140 - 10*u^138 - 6*u^136 - 3*u^134 - 2*u^132 + 3*u^130 + 
#                5*u^128 + 7*u^126 + 10*u^124 + 9*u^122 + 10*u^120 + 10*u^
#                118 + 7*u^116 + 6*u^114 + 3*u^112 - u^108 - 4*u^106 - 5*u^
#                104 - 5*u^102 - 6*u^100 - 5*u^98 - 4*u^96 - 3*u^94 - u^92 + u^
#                88 + 2*u^86 + 2*u^84 + 2*u^82 + u^80, 
#              u^207 + 2*u^205 + 2*u^203 + 2*u^201 + u^199 - u^195 - 3*u^193 - 
#                4*u^191 - 5*u^189 - 6*u^187 - 5*u^185 - 5*u^183 - 4*u^181 - u^
#                179 + 3*u^175 + 6*u^173 + 7*u^171 + 10*u^169 + 10*u^167 + 9*u^
#                165 + 10*u^163 + 7*u^161 + 5*u^159 + 3*u^157 - 2*u^155 - 3*u^
#                153 - 6*u^151 - 10*u^149 - 10*u^147 - 13*u^145 - 13*u^143 - 
#                10*u^141 - 10*u^139 - 6*u^137 - 3*u^135 - 2*u^133 + 3*u^131 + 
#                5*u^129 + 7*u^127 + 10*u^125 + 9*u^123 + 10*u^121 + 10*u^
#                119 + 7*u^117 + 6*u^115 + 3*u^113 - u^109 - 4*u^107 - 5*u^
#                105 - 5*u^103 - 6*u^101 - 5*u^99 - 4*u^97 - 3*u^95 - u^93 + u^
#                89 + 2*u^87 + 2*u^85 + 2*u^83 + u^81, 
#              u^207 + u^205 + u^203 + u^201 - u^195 - 2*u^193 - 2*u^191 - 3*u^
#                189 - 3*u^187 - 2*u^185 - 3*u^183 - u^181 + 3*u^175 + 3*u^
#                173 + 4*u^171 + 6*u^169 + 4*u^167 + 5*u^165 + 5*u^163 + 2*u^
#                161 + 3*u^159 - 2*u^155 - u^153 - 5*u^151 - 5*u^149 - 5*u^
#                147 - 8*u^145 - 5*u^143 - 5*u^141 - 5*u^139 - u^137 - 2*u^
#                135 + 3*u^131 + 2*u^129 + 5*u^127 + 5*u^125 + 4*u^123 + 6*u^
#                121 + 4*u^119 + 3*u^117 + 3*u^115 - u^109 - 3*u^107 - 2*u^
#                105 - 3*u^103 - 3*u^101 - 2*u^99 - 2*u^97 - u^95 + u^89 + u^
#                87 + u^85 + u^83, 0*u^0, 0*u^0, 0*u^0, 
#              u^208 + u^206 + 2*u^204 + u^202 + u^200 - u^196 - 2*u^194 - 3*u^
#                192 - 4*u^190 - 4*u^188 - 4*u^186 - 4*u^184 - 2*u^182 - 2*u^
#                180 + u^178 + 2*u^176 + 4*u^174 + 6*u^172 + 7*u^170 + 7*u^
#                168 + 8*u^166 + 6*u^164 + 6*u^162 + 4*u^160 + u^158 - 3*u^
#                154 - 5*u^152 - 6*u^150 - 9*u^148 - 9*u^146 - 9*u^144 - 9*u^
#                142 - 6*u^140 - 5*u^138 - 3*u^136 + u^132 + 4*u^130 + 6*u^
#                128 + 6*u^126 + 8*u^124 + 7*u^122 + 7*u^120 + 6*u^118 + 4*u^
#                116 + 2*u^114 + u^112 - 2*u^110 - 2*u^108 - 4*u^106 - 4*u^
#                104 - 4*u^102 - 4*u^100 - 3*u^98 - 2*u^96 - u^94 + u^90 + u^
#                88 + 2*u^86 + u^84 + u^82, u^204 + u^200 - u^192 - u^190 - u^
#                188 - 2*u^186 - u^184 - u^182 - 2*u^180 + u^178 - u^176 + u^
#                174 + 2*u^172 + u^170 + 3*u^168 + 3*u^166 + u^164 + 4*u^
#                162 + u^160 + u^158 + 2*u^156 - 2*u^154 - u^150 - 4*u^148 - u^
#                146 - 4*u^144 - 4*u^142 - u^140 - 4*u^138 - u^136 - 2*u^132 + 
#                2*u^130 + u^128 + u^126 + 4*u^124 + u^122 + 3*u^120 + 3*u^
#                118 + u^116 + 2*u^114 + u^112 - u^110 + u^108 - 2*u^106 - u^
#                104 - u^102 - 2*u^100 - u^98 - u^96 - u^94 + u^86 + u^82, 
#              u^206 + u^202 - u^194 - u^192 - u^190 - 2*u^188 - u^186 - u^
#                184 - 2*u^182 + u^180 - u^178 + u^176 + 2*u^174 + u^172 + 3*u^
#                170 + 3*u^168 + u^166 + 4*u^164 + u^162 + u^160 + 2*u^158 - 
#                2*u^156 - u^152 - 4*u^150 - u^148 - 4*u^146 - 4*u^144 - u^
#                142 - 4*u^140 - u^138 - 2*u^134 + 2*u^132 + u^130 + u^128 + 
#                4*u^126 + u^124 + 3*u^122 + 3*u^120 + u^118 + 2*u^116 + u^
#                114 - u^112 + u^110 - 2*u^108 - u^106 - u^104 - 2*u^102 - u^
#                100 - u^98 - u^96 + u^88 + u^84, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ u^204 + u^200 - u^192 - u^190 - u^188 - 2*u^186 - u^184 - u^182 - 
#                2*u^180 + u^178 - u^176 + u^174 + 2*u^172 + u^170 + 3*u^168 + 
#                3*u^166 + u^164 + 4*u^162 + u^160 + u^158 + 2*u^156 - 2*u^
#                154 - u^150 - 4*u^148 - u^146 - 4*u^144 - 4*u^142 - u^140 - 
#                4*u^138 - u^136 - 2*u^132 + 2*u^130 + u^128 + u^126 + 4*u^
#                124 + u^122 + 3*u^120 + 3*u^118 + u^116 + 2*u^114 + u^112 - u^
#                110 + u^108 - 2*u^106 - u^104 - u^102 - 2*u^100 - u^98 - u^
#                96 - u^94 + u^86 + u^82, u^205 + u^203 + u^201 + u^199 - u^
#                193 - 2*u^191 - 2*u^189 - 3*u^187 - 3*u^185 - 2*u^183 - 3*u^
#                181 - u^179 + 3*u^173 + 3*u^171 + 4*u^169 + 6*u^167 + 4*u^
#                165 + 5*u^163 + 5*u^161 + 2*u^159 + 3*u^157 - 2*u^153 - u^
#                151 - 5*u^149 - 5*u^147 - 5*u^145 - 8*u^143 - 5*u^141 - 5*u^
#                139 - 5*u^137 - u^135 - 2*u^133 + 3*u^129 + 2*u^127 + 5*u^
#                125 + 5*u^123 + 4*u^121 + 6*u^119 + 4*u^117 + 3*u^115 + 3*u^
#                113 - u^107 - 3*u^105 - 2*u^103 - 3*u^101 - 3*u^99 - 2*u^97 - 
#                2*u^95 - u^93 + u^87 + u^85 + u^83 + u^81, 0*u^0, 
#              u^206 + u^202 - u^194 - u^192 - u^190 - 2*u^188 - u^186 - u^
#                184 - 2*u^182 + u^180 - u^178 + u^176 + 2*u^174 + u^172 + 3*u^
#                170 + 3*u^168 + u^166 + 4*u^164 + u^162 + u^160 + 2*u^158 - 
#                2*u^156 - u^152 - 4*u^150 - u^148 - 4*u^146 - 4*u^144 - u^
#                142 - 4*u^140 - u^138 - 2*u^134 + 2*u^132 + u^130 + u^128 + 
#                4*u^126 + u^124 + 3*u^122 + 3*u^120 + u^118 + 2*u^116 + u^
#                114 - u^112 + u^110 - 2*u^108 - u^106 - u^104 - 2*u^102 - u^
#                100 - u^98 - u^96 + u^88 + u^84, 
#              u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^192 - 2*u^190 - 3*u^
#                188 - 3*u^186 - 2*u^184 - 3*u^182 - u^180 + 3*u^174 + 3*u^
#                172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^164 + 5*u^162 + 2*u^
#                160 + 3*u^158 - 2*u^154 - u^152 - 5*u^150 - 5*u^148 - 5*u^
#                146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^138 - u^136 - 2*u^
#                134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^124 + 4*u^122 + 6*u^
#                120 + 4*u^118 + 3*u^116 + 3*u^114 - u^108 - 3*u^106 - 2*u^
#                104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^96 - u^94 + u^88 + u^
#                86 + u^84 + u^82, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^
#                192 - 2*u^190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^
#                180 + 3*u^174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^
#                164 + 5*u^162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^
#                150 - 5*u^148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^
#                138 - u^136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^
#                124 + 4*u^122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^
#                108 - 3*u^106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^
#                96 - u^94 + u^88 + u^86 + u^84 + u^82, 
#              u^207 + u^205 + u^203 + u^201 - u^195 - 2*u^193 - 2*u^191 - 3*u^
#                189 - 3*u^187 - 2*u^185 - 3*u^183 - u^181 + 3*u^175 + 3*u^
#                173 + 4*u^171 + 6*u^169 + 4*u^167 + 5*u^165 + 5*u^163 + 2*u^
#                161 + 3*u^159 - 2*u^155 - u^153 - 5*u^151 - 5*u^149 - 5*u^
#                147 - 8*u^145 - 5*u^143 - 5*u^141 - 5*u^139 - u^137 - 2*u^
#                135 + 3*u^131 + 2*u^129 + 5*u^127 + 5*u^125 + 4*u^123 + 6*u^
#                121 + 4*u^119 + 3*u^117 + 3*u^115 - u^109 - 3*u^107 - 2*u^
#                105 - 3*u^103 - 3*u^101 - 2*u^99 - 2*u^97 - u^95 + u^89 + u^
#                87 + u^85 + u^83, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#              u^204 + u^200 - u^192 - u^190 - u^188 - 2*u^186 - u^184 - u^
#                182 - 2*u^180 + u^178 - u^176 + u^174 + 2*u^172 + u^170 + 3*u^
#                168 + 3*u^166 + u^164 + 4*u^162 + u^160 + u^158 + 2*u^156 - 
#                2*u^154 - u^150 - 4*u^148 - u^146 - 4*u^144 - 4*u^142 - u^
#                140 - 4*u^138 - u^136 - 2*u^132 + 2*u^130 + u^128 + u^126 + 
#                4*u^124 + u^122 + 3*u^120 + 3*u^118 + u^116 + 2*u^114 + u^
#                112 - u^110 + u^108 - 2*u^106 - u^104 - u^102 - 2*u^100 - u^
#                98 - u^96 - u^94 + u^86 + u^82, 
#              u^208 + u^204 - u^196 - u^194 - u^192 - 2*u^190 - u^188 - u^
#                186 - 2*u^184 + u^182 - u^180 + u^178 + 2*u^176 + u^174 + 3*u^
#                172 + 3*u^170 + u^168 + 4*u^166 + u^164 + u^162 + 2*u^160 - 
#                2*u^158 - u^154 - 4*u^152 - u^150 - 4*u^148 - 4*u^146 - u^
#                144 - 4*u^142 - u^140 - 2*u^136 + 2*u^134 + u^132 + u^130 + 
#                4*u^128 + u^126 + 3*u^124 + 3*u^122 + u^120 + 2*u^118 + u^
#                116 - u^114 + u^112 - 2*u^110 - u^108 - u^106 - 2*u^104 - u^
#                102 - u^100 - u^98 + u^90 + u^86, 
#              u^206 + u^202 - u^194 - u^192 - u^190 - 2*u^188 - u^186 - u^
#                184 - 2*u^182 + u^180 - u^178 + u^176 + 2*u^174 + u^172 + 3*u^
#                170 + 3*u^168 + u^166 + 4*u^164 + u^162 + u^160 + 2*u^158 - 
#                2*u^156 - u^152 - 4*u^150 - u^148 - 4*u^146 - 4*u^144 - u^
#                142 - 4*u^140 - u^138 - 2*u^134 + 2*u^132 + u^130 + u^128 + 
#                4*u^126 + u^124 + 3*u^122 + 3*u^120 + u^118 + 2*u^116 + u^
#                114 - u^112 + u^110 - 2*u^108 - u^106 - u^104 - 2*u^102 - u^
#                100 - u^98 - u^96 + u^88 + u^84, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ u^202 + u^198 - u^190 - u^188 - u^186 - 2*u^184 - u^182 - u^180 - 
#                2*u^178 + u^176 - u^174 + u^172 + 2*u^170 + u^168 + 3*u^166 + 
#                3*u^164 + u^162 + 4*u^160 + u^158 + u^156 + 2*u^154 - 2*u^
#                152 - u^148 - 4*u^146 - u^144 - 4*u^142 - 4*u^140 - u^138 - 
#                4*u^136 - u^134 - 2*u^130 + 2*u^128 + u^126 + u^124 + 4*u^
#                122 + u^120 + 3*u^118 + 3*u^116 + u^114 + 2*u^112 + u^110 - u^
#                108 + u^106 - 2*u^104 - u^102 - u^100 - 2*u^98 - u^96 - u^
#                94 - u^92 + u^84 + u^80, u^205 + u^203 + u^201 + u^199 - u^
#                193 - 2*u^191 - 2*u^189 - 3*u^187 - 3*u^185 - 2*u^183 - 3*u^
#                181 - u^179 + 3*u^173 + 3*u^171 + 4*u^169 + 6*u^167 + 4*u^
#                165 + 5*u^163 + 5*u^161 + 2*u^159 + 3*u^157 - 2*u^153 - u^
#                151 - 5*u^149 - 5*u^147 - 5*u^145 - 8*u^143 - 5*u^141 - 5*u^
#                139 - 5*u^137 - u^135 - 2*u^133 + 3*u^129 + 2*u^127 + 5*u^
#                125 + 5*u^123 + 4*u^121 + 6*u^119 + 4*u^117 + 3*u^115 + 3*u^
#                113 - u^107 - 3*u^105 - 2*u^103 - 3*u^101 - 3*u^99 - 2*u^97 - 
#                2*u^95 - u^93 + u^87 + u^85 + u^83 + u^81, 0*u^0, 
#              u^204 + u^200 - u^192 - u^190 - u^188 - 2*u^186 - u^184 - u^
#                182 - 2*u^180 + u^178 - u^176 + u^174 + 2*u^172 + u^170 + 3*u^
#                168 + 3*u^166 + u^164 + 4*u^162 + u^160 + u^158 + 2*u^156 - 
#                2*u^154 - u^150 - 4*u^148 - u^146 - 4*u^144 - 4*u^142 - u^
#                140 - 4*u^138 - u^136 - 2*u^132 + 2*u^130 + u^128 + u^126 + 
#                4*u^124 + u^122 + 3*u^120 + 3*u^118 + u^116 + 2*u^114 + u^
#                112 - u^110 + u^108 - 2*u^106 - u^104 - u^102 - 2*u^100 - u^
#                98 - u^96 - u^94 + u^86 + u^82, 
#              u^206 + 2*u^204 + u^202 + u^200 - u^196 - u^194 - 3*u^192 - 3*u^
#                190 - 3*u^188 - 4*u^186 - 2*u^184 - 2*u^182 - 2*u^180 + 2*u^
#                178 + u^176 + 3*u^174 + 6*u^172 + 4*u^170 + 7*u^168 + 6*u^
#                166 + 3*u^164 + 6*u^162 + 2*u^160 + u^156 - 5*u^154 - 3*u^
#                152 - 4*u^150 - 9*u^148 - 5*u^146 - 8*u^144 - 8*u^142 - 2*u^
#                140 - 5*u^138 - u^136 + 2*u^134 - u^132 + 5*u^130 + 5*u^128 + 
#                4*u^126 + 8*u^124 + 4*u^122 + 5*u^120 + 6*u^118 + u^116 + 2*u^
#                114 - 3*u^110 - u^108 - 4*u^106 - 4*u^104 - 2*u^102 - 4*u^
#                100 - 2*u^98 - u^96 - u^94 + u^92 + u^90 + u^88 + 2*u^86 + u^
#                84 + u^82 - u^78, u^206 + u^204 + u^202 + u^200 - u^194 - 2*u^
#                192 - 2*u^190 - 3*u^188 - 3*u^186 - 2*u^184 - 3*u^182 - u^
#                180 + 3*u^174 + 3*u^172 + 4*u^170 + 6*u^168 + 4*u^166 + 5*u^
#                164 + 5*u^162 + 2*u^160 + 3*u^158 - 2*u^154 - u^152 - 5*u^
#                150 - 5*u^148 - 5*u^146 - 8*u^144 - 5*u^142 - 5*u^140 - 5*u^
#                138 - u^136 - 2*u^134 + 3*u^130 + 2*u^128 + 5*u^126 + 5*u^
#                124 + 4*u^122 + 6*u^120 + 4*u^118 + 3*u^116 + 3*u^114 - u^
#                108 - 3*u^106 - 2*u^104 - 3*u^102 - 3*u^100 - 2*u^98 - 2*u^
#                96 - u^94 + u^88 + u^86 + u^84 + u^82, 
#              u^207 + 2*u^205 + u^203 + u^201 - u^197 - u^195 - 3*u^193 - 3*u^
#                191 - 3*u^189 - 4*u^187 - 2*u^185 - 2*u^183 - 2*u^181 + 2*u^
#                179 + u^177 + 3*u^175 + 6*u^173 + 4*u^171 + 7*u^169 + 6*u^
#                167 + 3*u^165 + 6*u^163 + 2*u^161 + u^157 - 5*u^155 - 3*u^
#                153 - 4*u^151 - 9*u^149 - 5*u^147 - 8*u^145 - 8*u^143 - 2*u^
#                141 - 5*u^139 - u^137 + 2*u^135 - u^133 + 5*u^131 + 5*u^129 + 
#                4*u^127 + 8*u^125 + 4*u^123 + 5*u^121 + 6*u^119 + u^117 + 2*u^
#                115 - 3*u^111 - u^109 - 4*u^107 - 4*u^105 - 2*u^103 - 4*u^
#                101 - 2*u^99 - u^97 - u^95 + u^93 + u^91 + u^89 + 2*u^87 + u^
#                85 + u^83 - u^79, 0*u^0, 0*u^0, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, 0*u^0, u^206 + u^202 - u^194 - u^192 - u^190 - 2*u^
#                188 - u^186 - u^184 - 2*u^182 + u^180 - u^178 + u^176 + 2*u^
#                174 + u^172 + 3*u^170 + 3*u^168 + u^166 + 4*u^164 + u^162 + u^
#                160 + 2*u^158 - 2*u^156 - u^152 - 4*u^150 - u^148 - 4*u^146 - 
#                4*u^144 - u^142 - 4*u^140 - u^138 - 2*u^134 + 2*u^132 + u^
#                130 + u^128 + 4*u^126 + u^124 + 3*u^122 + 3*u^120 + u^118 + 
#                2*u^116 + u^114 - u^112 + u^110 - 2*u^108 - u^106 - u^104 - 
#                2*u^102 - u^100 - u^98 - u^96 + u^88 + u^84, 
#              u^206 + u^202 - u^194 - u^192 - u^190 - 2*u^188 - u^186 - u^
#                184 - 2*u^182 + u^180 - u^178 + u^176 + 2*u^174 + u^172 + 3*u^
#                170 + 3*u^168 + u^166 + 4*u^164 + u^162 + u^160 + 2*u^158 - 
#                2*u^156 - u^152 - 4*u^150 - u^148 - 4*u^146 - 4*u^144 - u^
#                142 - 4*u^140 - u^138 - 2*u^134 + 2*u^132 + u^130 + u^128 + 
#                4*u^126 + u^124 + 3*u^122 + 3*u^120 + u^118 + 2*u^116 + u^
#                114 - u^112 + u^110 - 2*u^108 - u^106 - u^104 - 2*u^102 - u^
#                100 - u^98 - u^96 + u^88 + u^84, 
#              u^208 + u^206 + 2*u^204 - u^198 - 2*u^196 - 2*u^194 - 3*u^192 - 
#                3*u^190 - 2*u^188 - 2*u^186 - u^184 + u^182 + 4*u^178 + 3*u^
#                176 + 4*u^174 + 6*u^172 + 4*u^170 + 4*u^168 + 4*u^166 + 2*u^
#                162 - u^160 - 4*u^158 - 2*u^156 - 6*u^154 - 5*u^152 - 4*u^
#                150 - 8*u^148 - 4*u^146 - 4*u^144 - 4*u^142 + 2*u^140 + 2*u^
#                136 + 5*u^134 + 2*u^132 + 6*u^130 + 6*u^128 + 3*u^126 + 6*u^
#                124 + 2*u^122 + 2*u^120 + 2*u^118 - 2*u^116 - 2*u^114 - 2*u^
#                112 - 5*u^110 - 2*u^108 - 4*u^106 - 3*u^104 - u^102 - 2*u^
#                100 + u^96 + u^94 + 2*u^92 + 2*u^90 + u^88 + 2*u^86 - u^
#                80 - u^78, 0*u^0, u^207 + u^205 - u^199 - u^197 - u^195 - 2*u^
#                193 - u^191 - u^189 - u^187 + u^185 + u^181 + 3*u^179 + u^
#                177 + 3*u^175 + 3*u^173 + u^171 + 3*u^169 - u^165 + u^163 - 
#                3*u^161 - 2*u^159 - 2*u^157 - 5*u^155 - u^153 - 3*u^151 - 4*u^
#                149 - 3*u^145 + 3*u^141 + 4*u^137 + 3*u^135 + u^133 + 5*u^
#                131 + 2*u^129 + 2*u^127 + 3*u^125 - u^123 + u^121 - 3*u^
#                117 - u^115 - 3*u^113 - 3*u^111 - u^109 - 3*u^107 - u^105 - u^
#                101 + u^99 + u^97 + u^95 + 2*u^93 + u^91 + u^89 + u^87 - u^
#                81 - u^79, u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^
#                182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^
#                166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 
#                4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^
#                130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 
#                3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^
#                92 + u^88 - u^80 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, u^206 - u^198 - u^194 - u^192 - u^188 + u^
#                184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 
#                2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^
#                152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 2*u^
#                132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#              u^208 - u^200 - u^196 - u^194 - u^190 + u^186 - u^184 + 2*u^
#                182 + u^180 + 3*u^176 + u^172 + 2*u^170 - 2*u^168 + u^166 - 
#                3*u^162 + u^160 - 3*u^158 - 2*u^156 + u^154 - 4*u^152 - 3*u^
#                146 + 3*u^144 + 4*u^138 - u^136 + 2*u^134 + 3*u^132 - u^130 + 
#                3*u^128 - u^124 + 2*u^122 - 2*u^120 - u^118 - 3*u^114 - u^
#                110 - 2*u^108 + u^106 - u^104 + u^100 + u^96 + u^94 + u^
#                90 - u^82, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^205 - u^197 - u^193 - u^191 - u^
#                187 + u^183 - u^181 + 2*u^179 + u^177 + 3*u^173 + u^169 + 2*u^
#                167 - 2*u^165 + u^163 - 3*u^159 + u^157 - 3*u^155 - 2*u^
#                153 + u^151 - 4*u^149 - 3*u^143 + 3*u^141 + 4*u^135 - u^133 + 
#                2*u^131 + 3*u^129 - u^127 + 3*u^125 - u^121 + 2*u^119 - 2*u^
#                117 - u^115 - 3*u^111 - u^107 - 2*u^105 + u^103 - u^101 + u^
#                97 + u^93 + u^91 + u^87 - u^79, 0*u^0, 
#              u^206 - u^198 - u^194 - u^192 - u^188 + u^184 - u^182 + 2*u^
#                180 + u^178 + 3*u^174 + u^170 + 2*u^168 - 2*u^166 + u^164 - 
#                3*u^160 + u^158 - 3*u^156 - 2*u^154 + u^152 - 4*u^150 - 3*u^
#                144 + 3*u^142 + 4*u^136 - u^134 + 2*u^132 + 3*u^130 - u^128 + 
#                3*u^126 - u^122 + 2*u^120 - 2*u^118 - u^116 - 3*u^112 - u^
#                108 - 2*u^106 + u^104 - u^102 + u^98 + u^94 + u^92 + u^88 - u^
#                80, 0*u^0, 0*u^0, u^207 - u^199 - u^195 - u^193 - u^189 + u^
#                185 - u^183 + 2*u^181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 
#                2*u^167 + u^165 - 3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^
#                153 - 4*u^151 - 3*u^145 + 3*u^143 + 4*u^137 - u^135 + 2*u^
#                133 + 3*u^131 - u^129 + 3*u^127 - u^123 + 2*u^121 - 2*u^
#                119 - u^117 - 3*u^113 - u^109 - 2*u^107 + u^105 - u^103 + u^
#                99 + u^95 + u^93 + u^89 - u^81, 0*u^0, 0*u^0, 0*u^0, 
#              u^207 + u^205 - u^199 - u^197 - u^195 - 2*u^193 - u^191 - u^
#                189 - u^187 + u^185 + u^181 + 3*u^179 + u^177 + 3*u^175 + 3*u^
#                173 + u^171 + 3*u^169 - u^165 + u^163 - 3*u^161 - 2*u^159 - 
#                2*u^157 - 5*u^155 - u^153 - 3*u^151 - 4*u^149 - 3*u^145 + 3*u^
#                141 + 4*u^137 + 3*u^135 + u^133 + 5*u^131 + 2*u^129 + 2*u^
#                127 + 3*u^125 - u^123 + u^121 - 3*u^117 - u^115 - 3*u^113 - 
#                3*u^111 - u^109 - 3*u^107 - u^105 - u^101 + u^99 + u^97 + u^
#                95 + 2*u^93 + u^91 + u^89 + u^87 - u^81 - u^79, 0*u^0, 
#              u^208 + u^206 - u^200 - u^198 - u^196 - 2*u^194 - u^192 - u^
#                190 - u^188 + u^186 + u^182 + 3*u^180 + u^178 + 3*u^176 + 3*u^
#                174 + u^172 + 3*u^170 - u^166 + u^164 - 3*u^162 - 2*u^160 - 
#                2*u^158 - 5*u^156 - u^154 - 3*u^152 - 4*u^150 - 3*u^146 + 3*u^
#                142 + 4*u^138 + 3*u^136 + u^134 + 5*u^132 + 2*u^130 + 2*u^
#                128 + 3*u^126 - u^124 + u^122 - 3*u^118 - u^116 - 3*u^114 - 
#                3*u^112 - u^110 - 3*u^108 - u^106 - u^102 + u^100 + u^98 + u^
#                96 + 2*u^94 + u^92 + u^90 + u^88 - u^82 - u^80, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#              0*u^0, 0*u^0, 0*u^0, 0*u^0, u^206 - u^198 - u^194 - u^192 - u^
#                188 + u^184 - u^182 + 2*u^180 + u^178 + 3*u^174 + u^170 + 2*u^
#                168 - 2*u^166 + u^164 - 3*u^160 + u^158 - 3*u^156 - 2*u^
#                154 + u^152 - 4*u^150 - 3*u^144 + 3*u^142 + 4*u^136 - u^134 + 
#                2*u^132 + 3*u^130 - u^128 + 3*u^126 - u^122 + 2*u^120 - 2*u^
#                118 - u^116 - 3*u^112 - u^108 - 2*u^106 + u^104 - u^102 + u^
#                98 + u^94 + u^92 + u^88 - u^80, 0*u^0, 
#              u^207 - u^199 - u^195 - u^193 - u^189 + u^185 - u^183 + 2*u^
#                181 + u^179 + 3*u^175 + u^171 + 2*u^169 - 2*u^167 + u^165 - 
#                3*u^161 + u^159 - 3*u^157 - 2*u^155 + u^153 - 4*u^151 - 3*u^
#                145 + 3*u^143 + 4*u^137 - u^135 + 2*u^133 + 3*u^131 - u^129 + 
#                3*u^127 - u^123 + 2*u^121 - 2*u^119 - u^117 - 3*u^113 - u^
#                109 - 2*u^107 + u^105 - u^103 + u^99 + u^95 + u^93 + u^89 - u^
#                81, u^208 - u^200 - u^196 - u^194 - u^190 + u^186 - u^184 + 
#                2*u^182 + u^180 + 3*u^176 + u^172 + 2*u^170 - 2*u^168 + u^
#                166 - 3*u^162 + u^160 - 3*u^158 - 2*u^156 + u^154 - 4*u^152 - 
#                3*u^146 + 3*u^144 + 4*u^138 - u^136 + 2*u^134 + 3*u^132 - u^
#                130 + 3*u^128 - u^124 + 2*u^122 - 2*u^120 - u^118 - 3*u^
#                114 - u^110 - 2*u^108 + u^106 - u^104 + u^100 + u^96 + u^
#                94 + u^90 - u^82 ] ], 
#      [ [ u^210 + u^208 + 2*u^206 + u^204 + u^202 - u^198 - 2*u^196 - 3*u^
#                194 - 4*u^192 - 4*u^190 - 4*u^188 - 4*u^186 - 2*u^184 - 2*u^
#                182 + u^180 + 2*u^178 + 4*u^176 + 6*u^174 + 7*u^172 + 7*u^
#                170 + 8*u^168 + 6*u^166 + 6*u^164 + 4*u^162 + u^160 - 3*u^
#                156 - 5*u^154 - 6*u^152 - 9*u^150 - 9*u^148 - 9*u^146 - 9*u^
#                144 - 6*u^142 - 5*u^140 - 3*u^138 + u^134 + 4*u^132 + 6*u^
#                130 + 6*u^128 + 8*u^126 + 7*u^124 + 7*u^122 + 6*u^120 + 4*u^
#                118 + 2*u^116 + u^114 - 2*u^112 - 2*u^110 - 4*u^108 - 4*u^
#                106 - 4*u^104 - 4*u^102 - 3*u^100 - 2*u^98 - u^96 + u^92 + u^
#                90 + 2*u^88 + u^86 + u^84, 0*u^0, 
#              u^209 + u^207 + u^205 + u^203 - u^197 - 2*u^195 - 2*u^193 - 3*u^
#                191 - 3*u^189 - 2*u^187 - 3*u^185 - u^183 + 3*u^177 + 3*u^
#                175 + 4*u^173 + 6*u^171 + 4*u^169 + 5*u^167 + 5*u^165 + 2*u^
#                163 + 3*u^161 - 2*u^157 - u^155 - 5*u^153 - 5*u^151 - 5*u^
#                149 - 8*u^147 - 5*u^145 - 5*u^143 - 5*u^141 - u^139 - 2*u^
#                137 + 3*u^133 + 2*u^131 + 5*u^129 + 5*u^127 + 4*u^125 + 6*u^
#                123 + 4*u^121 + 3*u^119 + 3*u^117 - u^111 - 3*u^109 - 2*u^
#                107 - 3*u^105 - 3*u^103 - 2*u^101 - 2*u^99 - u^97 + u^91 + u^
#                89 + u^87 + u^85, u^208 + u^204 - u^196 - u^194 - u^192 - 2*u^
#                190 - u^188 - u^186 - 2*u^184 + u^182 - u^180 + u^178 + 2*u^
#                176 + u^174 + 3*u^172 + 3*u^170 + u^168 + 4*u^166 + u^164 + u^
#                162 + 2*u^160 - 2*u^158 - u^154 - 4*u^152 - u^150 - 4*u^148 - 
#                4*u^146 - u^144 - 4*u^142 - u^140 - 2*u^136 + 2*u^134 + u^
#                132 + u^130 + 4*u^128 + u^126 + 3*u^124 + 3*u^122 + u^120 + 
#                2*u^118 + u^116 - u^114 + u^112 - 2*u^110 - u^108 - u^106 - 
#                2*u^104 - u^102 - u^100 - u^98 + u^90 + u^86 ], 
#          [ 0*u^0, u^210 + u^208 + u^206 + u^204 - u^198 - 2*u^196 - 2*u^
#                194 - 3*u^192 - 3*u^190 - 2*u^188 - 3*u^186 - u^184 + 3*u^
#                178 + 3*u^176 + 4*u^174 + 6*u^172 + 4*u^170 + 5*u^168 + 5*u^
#                166 + 2*u^164 + 3*u^162 - 2*u^158 - u^156 - 5*u^154 - 5*u^
#                152 - 5*u^150 - 8*u^148 - 5*u^146 - 5*u^144 - 5*u^142 - u^
#                140 - 2*u^138 + 3*u^134 + 2*u^132 + 5*u^130 + 5*u^128 + 4*u^
#                126 + 6*u^124 + 4*u^122 + 3*u^120 + 3*u^118 - u^112 - 3*u^
#                110 - 2*u^108 - 3*u^106 - 3*u^104 - 2*u^102 - 2*u^100 - u^
#                98 + u^92 + u^90 + u^88 + u^86, 0*u^0, 0*u^0 ], 
#          [ u^209 + u^207 + u^205 + u^203 - u^197 - 2*u^195 - 2*u^193 - 3*u^
#                191 - 3*u^189 - 2*u^187 - 3*u^185 - u^183 + 3*u^177 + 3*u^
#                175 + 4*u^173 + 6*u^171 + 4*u^169 + 5*u^167 + 5*u^165 + 2*u^
#                163 + 3*u^161 - 2*u^157 - u^155 - 5*u^153 - 5*u^151 - 5*u^
#                149 - 8*u^147 - 5*u^145 - 5*u^143 - 5*u^141 - u^139 - 2*u^
#                137 + 3*u^133 + 2*u^131 + 5*u^129 + 5*u^127 + 4*u^125 + 6*u^
#                123 + 4*u^121 + 3*u^119 + 3*u^117 - u^111 - 3*u^109 - 2*u^
#                107 - 3*u^105 - 3*u^103 - 2*u^101 - 2*u^99 - u^97 + u^91 + u^
#                89 + u^87 + u^85, 0*u^0, u^210 + u^208 + u^206 + u^204 - u^
#                198 - 2*u^196 - 2*u^194 - 3*u^192 - 3*u^190 - 2*u^188 - 3*u^
#                186 - u^184 + 3*u^178 + 3*u^176 + 4*u^174 + 6*u^172 + 4*u^
#                170 + 5*u^168 + 5*u^166 + 2*u^164 + 3*u^162 - 2*u^158 - u^
#                156 - 5*u^154 - 5*u^152 - 5*u^150 - 8*u^148 - 5*u^146 - 5*u^
#                144 - 5*u^142 - u^140 - 2*u^138 + 3*u^134 + 2*u^132 + 5*u^
#                130 + 5*u^128 + 4*u^126 + 6*u^124 + 4*u^122 + 3*u^120 + 3*u^
#                118 - u^112 - 3*u^110 - 2*u^108 - 3*u^106 - 3*u^104 - 2*u^
#                102 - 2*u^100 - u^98 + u^92 + u^90 + u^88 + u^86, 0*u^0 ], 
#          [ u^208 + u^204 - u^196 - u^194 - u^192 - 2*u^190 - u^188 - u^186 - 
#                2*u^184 + u^182 - u^180 + u^178 + 2*u^176 + u^174 + 3*u^172 + 
#                3*u^170 + u^168 + 4*u^166 + u^164 + u^162 + 2*u^160 - 2*u^
#                158 - u^154 - 4*u^152 - u^150 - 4*u^148 - 4*u^146 - u^144 - 
#                4*u^142 - u^140 - 2*u^136 + 2*u^134 + u^132 + u^130 + 4*u^
#                128 + u^126 + 3*u^124 + 3*u^122 + u^120 + 2*u^118 + u^116 - u^
#                114 + u^112 - 2*u^110 - u^108 - u^106 - 2*u^104 - u^102 - u^
#                100 - u^98 + u^90 + u^86, 0*u^0, 0*u^0, 
#              u^210 + u^206 - u^198 - u^196 - u^194 - 2*u^192 - u^190 - u^
#                188 - 2*u^186 + u^184 - u^182 + u^180 + 2*u^178 + u^176 + 3*u^
#                174 + 3*u^172 + u^170 + 4*u^168 + u^166 + u^164 + 2*u^162 - 
#                2*u^160 - u^156 - 4*u^154 - u^152 - 4*u^150 - 4*u^148 - u^
#                146 - 4*u^144 - u^142 - 2*u^138 + 2*u^136 + u^134 + u^132 + 
#                4*u^130 + u^128 + 3*u^126 + 3*u^124 + u^122 + 2*u^120 + u^
#                118 - u^116 + u^114 - 2*u^112 - u^110 - u^108 - 2*u^106 - u^
#                104 - u^102 - u^100 + u^92 + u^88 ] ], 
#      [ [ u^212 - u^204 - u^200 - u^198 - u^194 + u^190 - u^188 + 2*u^186 + u^
#                184 + 3*u^180 + u^176 + 2*u^174 - 2*u^172 + u^170 - 3*u^
#                166 + u^164 - 3*u^162 - 2*u^160 + u^158 - 4*u^156 - 3*u^150 + 
#                3*u^148 + 4*u^142 - u^140 + 2*u^138 + 3*u^136 - u^134 + 3*u^
#                132 - u^128 + 2*u^126 - 2*u^124 - u^122 - 3*u^118 - u^114 - 
#                2*u^112 + u^110 - u^108 + u^104 + u^100 + u^98 + u^94 - u^86, 
#              0*u^0 ], 
#          [ 0*u^0, u^212 - u^204 - u^200 - u^198 - u^194 + u^190 - u^188 + 
#                2*u^186 + u^184 + 3*u^180 + u^176 + 2*u^174 - 2*u^172 + u^
#                170 - 3*u^166 + u^164 - 3*u^162 - 2*u^160 + u^158 - 4*u^156 - 
#                3*u^150 + 3*u^148 + 4*u^142 - u^140 + 2*u^138 + 3*u^136 - u^
#                134 + 3*u^132 - u^128 + 2*u^126 - 2*u^124 - u^122 - 3*u^
#                118 - u^114 - 2*u^112 + u^110 - u^108 + u^104 + u^100 + u^
#                98 + u^94 - u^86 ] ], 
#      [ [ u^214 - u^206 - u^202 - u^200 - u^196 + u^192 - u^190 + 2*u^188 + u^
#                186 + 3*u^182 + u^178 + 2*u^176 - 2*u^174 + u^172 - 3*u^
#                168 + u^166 - 3*u^164 - 2*u^162 + u^160 - 4*u^158 - 3*u^152 + 
#                3*u^150 + 4*u^144 - u^142 + 2*u^140 + 3*u^138 - u^136 + 3*u^
#                134 - u^130 + 2*u^128 - 2*u^126 - u^124 - 3*u^120 - u^116 - 
#                2*u^114 + u^112 - u^110 + u^106 + u^102 + u^100 + u^96 - u^88,
#              0*u^0, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, u^214 + u^208 - u^206 - 2*u^200 - u^196 - 2*u^194 + u^
#                192 - 2*u^190 + 2*u^186 - 2*u^184 + 3*u^182 + 2*u^180 - u^
#                178 + 5*u^176 + 5*u^170 - 3*u^168 + u^166 + 2*u^164 - 5*u^
#                162 + 2*u^160 - 2*u^158 - 5*u^156 + 2*u^154 - 5*u^152 - 2*u^
#                150 + 2*u^148 - 5*u^146 + 2*u^144 + u^142 - 3*u^140 + 5*u^
#                138 + 5*u^132 - u^130 + 2*u^128 + 3*u^126 - 2*u^124 + 2*u^
#                122 - 2*u^118 + u^116 - 2*u^114 - u^112 - 2*u^108 - u^102 + u^
#                100 + u^94, u^211 + u^205 - u^203 - 2*u^197 - u^193 - 2*u^
#                191 + u^189 - 2*u^187 + 2*u^183 - 2*u^181 + 3*u^179 + 2*u^
#                177 - u^175 + 5*u^173 + 5*u^167 - 3*u^165 + u^163 + 2*u^161 - 
#                5*u^159 + 2*u^157 - 2*u^155 - 5*u^153 + 2*u^151 - 5*u^149 - 
#                2*u^147 + 2*u^145 - 5*u^143 + 2*u^141 + u^139 - 3*u^137 + 5*u^
#                135 + 5*u^129 - u^127 + 2*u^125 + 3*u^123 - 2*u^121 + 2*u^
#                119 - 2*u^115 + u^113 - 2*u^111 - u^109 - 2*u^105 - u^99 + u^
#                97 + u^91, 0*u^0 ], 
#          [ 0*u^0, u^211 + u^205 - u^203 - 2*u^197 - u^193 - 2*u^191 + u^
#                189 - 2*u^187 + 2*u^183 - 2*u^181 + 3*u^179 + 2*u^177 - u^
#                175 + 5*u^173 + 5*u^167 - 3*u^165 + u^163 + 2*u^161 - 5*u^
#                159 + 2*u^157 - 2*u^155 - 5*u^153 + 2*u^151 - 5*u^149 - 2*u^
#                147 + 2*u^145 - 5*u^143 + 2*u^141 + u^139 - 3*u^137 + 5*u^
#                135 + 5*u^129 - u^127 + 2*u^125 + 3*u^123 - 2*u^121 + 2*u^
#                119 - 2*u^115 + u^113 - 2*u^111 - u^109 - 2*u^105 - u^99 + u^
#                97 + u^91, u^214 + u^208 - u^206 - 2*u^200 - u^196 - 2*u^
#                194 + u^192 - 2*u^190 + 2*u^186 - 2*u^184 + 3*u^182 + 2*u^
#                180 - u^178 + 5*u^176 + 5*u^170 - 3*u^168 + u^166 + 2*u^164 - 
#                5*u^162 + 2*u^160 - 2*u^158 - 5*u^156 + 2*u^154 - 5*u^152 - 
#                2*u^150 + 2*u^148 - 5*u^146 + 2*u^144 + u^142 - 3*u^140 + 5*u^
#                138 + 5*u^132 - u^130 + 2*u^128 + 3*u^126 - 2*u^124 + 2*u^
#                122 - 2*u^118 + u^116 - 2*u^114 - u^112 - 2*u^108 - u^102 + u^
#                100 + u^94, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, u^214 - u^206 - u^202 - u^200 - u^196 + u^
#                192 - u^190 + 2*u^188 + u^186 + 3*u^182 + u^178 + 2*u^176 - 
#                2*u^174 + u^172 - 3*u^168 + u^166 - 3*u^164 - 2*u^162 + u^
#                160 - 4*u^158 - 3*u^152 + 3*u^150 + 4*u^144 - u^142 + 2*u^
#                140 + 3*u^138 - u^136 + 3*u^134 - u^130 + 2*u^128 - 2*u^
#                126 - u^124 - 3*u^120 - u^116 - 2*u^114 + u^112 - u^110 + u^
#                106 + u^102 + u^100 + u^96 - u^88 ] ], 
#      [ [ u^216 - u^208 - u^204 - u^202 - u^198 + u^194 - u^192 + 2*u^190 + u^
#                188 + 3*u^184 + u^180 + 2*u^178 - 2*u^176 + u^174 - 3*u^
#                170 + u^168 - 3*u^166 - 2*u^164 + u^162 - 4*u^160 - 3*u^154 + 
#                3*u^152 + 4*u^146 - u^144 + 2*u^142 + 3*u^140 - u^138 + 3*u^
#                136 - u^132 + 2*u^130 - 2*u^128 - u^126 - 3*u^122 - u^118 - 
#                2*u^116 + u^114 - u^112 + u^108 + u^104 + u^102 + u^98 - u^90,
#              0*u^0, u^215 - u^207 - u^203 - u^201 - u^197 + u^193 - u^191 + 
#                2*u^189 + u^187 + 3*u^183 + u^179 + 2*u^177 - 2*u^175 + u^
#                173 - 3*u^169 + u^167 - 3*u^165 - 2*u^163 + u^161 - 4*u^159 - 
#                3*u^153 + 3*u^151 + 4*u^145 - u^143 + 2*u^141 + 3*u^139 - u^
#                137 + 3*u^135 - u^131 + 2*u^129 - 2*u^127 - u^125 - 3*u^
#                121 - u^117 - 2*u^115 + u^113 - u^111 + u^107 + u^103 + u^
#                101 + u^97 - u^89, 0*u^0 ], 
#          [ 0*u^0, u^216 + u^210 - u^208 - 2*u^202 - u^198 - 2*u^196 + u^
#                194 - 2*u^192 + 2*u^188 - 2*u^186 + 3*u^184 + 2*u^182 - u^
#                180 + 5*u^178 + 5*u^172 - 3*u^170 + u^168 + 2*u^166 - 5*u^
#                164 + 2*u^162 - 2*u^160 - 5*u^158 + 2*u^156 - 5*u^154 - 2*u^
#                152 + 2*u^150 - 5*u^148 + 2*u^146 + u^144 - 3*u^142 + 5*u^
#                140 + 5*u^134 - u^132 + 2*u^130 + 3*u^128 - 2*u^126 + 2*u^
#                124 - 2*u^120 + u^118 - 2*u^116 - u^114 - 2*u^110 - u^104 + u^
#                102 + u^96, 0*u^0, 0*u^0 ], 
#          [ u^215 - u^207 - u^203 - u^201 - u^197 + u^193 - u^191 + 2*u^
#                189 + u^187 + 3*u^183 + u^179 + 2*u^177 - 2*u^175 + u^173 - 
#                3*u^169 + u^167 - 3*u^165 - 2*u^163 + u^161 - 4*u^159 - 3*u^
#                153 + 3*u^151 + 4*u^145 - u^143 + 2*u^141 + 3*u^139 - u^137 + 
#                3*u^135 - u^131 + 2*u^129 - 2*u^127 - u^125 - 3*u^121 - u^
#                117 - 2*u^115 + u^113 - u^111 + u^107 + u^103 + u^101 + u^
#                97 - u^89, 0*u^0, u^216 - u^208 - u^204 - u^202 - u^198 + u^
#                194 - u^192 + 2*u^190 + u^188 + 3*u^184 + u^180 + 2*u^178 - 
#                2*u^176 + u^174 - 3*u^170 + u^168 - 3*u^166 - 2*u^164 + u^
#                162 - 4*u^160 - 3*u^154 + 3*u^152 + 4*u^146 - u^144 + 2*u^
#                142 + 3*u^140 - u^138 + 3*u^136 - u^132 + 2*u^130 - 2*u^
#                128 - u^126 - 3*u^122 - u^118 - 2*u^116 + u^114 - u^112 + u^
#                108 + u^104 + u^102 + u^98 - u^90, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, u^216 - u^214 - u^208 + u^206 - u^204 + u^
#                200 - u^198 + u^196 + u^194 - 2*u^192 + 3*u^190 - u^188 - u^
#                186 + 3*u^184 - 3*u^182 + u^180 + u^178 - 4*u^176 + 3*u^
#                174 - u^172 - 3*u^170 + 4*u^168 - 4*u^166 + u^164 + 3*u^162 - 
#                5*u^160 + 4*u^158 - 3*u^154 + 6*u^152 - 3*u^150 + 4*u^146 - 
#                5*u^144 + 3*u^142 + u^140 - 4*u^138 + 4*u^136 - 3*u^134 - u^
#                132 + 3*u^130 - 4*u^128 + u^126 + u^124 - 3*u^122 + 3*u^
#                120 - u^118 - u^116 + 3*u^114 - 2*u^112 + u^110 + u^108 - u^
#                106 + u^104 - u^100 + u^98 - u^96 - u^90 + u^88 ] ], 
#      [ [ u^218 - u^210 - u^206 - u^204 - u^200 + u^196 - u^194 + 2*u^192 + u^
#                190 + 3*u^186 + u^182 + 2*u^180 - 2*u^178 + u^176 - 3*u^
#                172 + u^170 - 3*u^168 - 2*u^166 + u^164 - 4*u^162 - 3*u^156 + 
#                3*u^154 + 4*u^148 - u^146 + 2*u^144 + 3*u^142 - u^140 + 3*u^
#                138 - u^134 + 2*u^132 - 2*u^130 - u^128 - 3*u^124 - u^120 - 
#                2*u^118 + u^116 - u^114 + u^110 + u^106 + u^104 + u^100 - u^92
#                , u^217 - u^209 - u^205 - u^203 - u^199 + u^195 - u^193 + 2*u^
#                191 + u^189 + 3*u^185 + u^181 + 2*u^179 - 2*u^177 + u^175 - 
#                3*u^171 + u^169 - 3*u^167 - 2*u^165 + u^163 - 4*u^161 - 3*u^
#                155 + 3*u^153 + 4*u^147 - u^145 + 2*u^143 + 3*u^141 - u^139 + 
#                3*u^137 - u^133 + 2*u^131 - 2*u^129 - u^127 - 3*u^123 - u^
#                119 - 2*u^117 + u^115 - u^113 + u^109 + u^105 + u^103 + u^
#                99 - u^91 ], 
#          [ u^217 - u^209 - u^205 - u^203 - u^199 + u^195 - u^193 + 2*u^
#                191 + u^189 + 3*u^185 + u^181 + 2*u^179 - 2*u^177 + u^175 - 
#                3*u^171 + u^169 - 3*u^167 - 2*u^165 + u^163 - 4*u^161 - 3*u^
#                155 + 3*u^153 + 4*u^147 - u^145 + 2*u^143 + 3*u^141 - u^139 + 
#                3*u^137 - u^133 + 2*u^131 - 2*u^129 - u^127 - 3*u^123 - u^
#                119 - 2*u^117 + u^115 - u^113 + u^109 + u^105 + u^103 + u^
#                99 - u^91, u^218 - u^210 - u^206 - u^204 - u^200 + u^196 - u^
#                194 + 2*u^192 + u^190 + 3*u^186 + u^182 + 2*u^180 - 2*u^
#                178 + u^176 - 3*u^172 + u^170 - 3*u^168 - 2*u^166 + u^164 - 
#                4*u^162 - 3*u^156 + 3*u^154 + 4*u^148 - u^146 + 2*u^144 + 3*u^
#                142 - u^140 + 3*u^138 - u^134 + 2*u^132 - 2*u^130 - u^128 - 
#                3*u^124 - u^120 - 2*u^118 + u^116 - u^114 + u^110 + u^106 + u^
#                104 + u^100 - u^92 ] ], 
#      [ [ u^220 + u^216 - u^208 - u^206 - u^204 - 2*u^202 - u^200 - u^198 - 
#                2*u^196 + u^194 - u^192 + u^190 + 2*u^188 + u^186 + 3*u^184 + 
#                3*u^182 + u^180 + 4*u^178 + u^176 + u^174 + 2*u^172 - 2*u^
#                170 - u^166 - 4*u^164 - u^162 - 4*u^160 - 4*u^158 - u^156 - 
#                4*u^154 - u^152 - 2*u^148 + 2*u^146 + u^144 + u^142 + 4*u^
#                140 + u^138 + 3*u^136 + 3*u^134 + u^132 + 2*u^130 + u^128 - u^
#                126 + u^124 - 2*u^122 - u^120 - u^118 - 2*u^116 - u^114 - u^
#                112 - u^110 + u^102 + u^98, 0*u^0, 0*u^0, 
#              u^218 + u^214 - u^206 - u^204 - u^202 - 2*u^200 - u^198 - u^
#                196 - 2*u^194 + u^192 - u^190 + u^188 + 2*u^186 + u^184 + 3*u^
#                182 + 3*u^180 + u^178 + 4*u^176 + u^174 + u^172 + 2*u^170 - 
#                2*u^168 - u^164 - 4*u^162 - u^160 - 4*u^158 - 4*u^156 - u^
#                154 - 4*u^152 - u^150 - 2*u^146 + 2*u^144 + u^142 + u^140 + 
#                4*u^138 + u^136 + 3*u^134 + 3*u^132 + u^130 + 2*u^128 + u^
#                126 - u^124 + u^122 - 2*u^120 - u^118 - u^116 - 2*u^114 - u^
#                112 - u^110 - u^108 + u^100 + u^96, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, u^220 - u^212 - u^208 - u^206 - u^202 + u^198 - u^196 + 
#                2*u^194 + u^192 + 3*u^188 + u^184 + 2*u^182 - 2*u^180 + u^
#                178 - 3*u^174 + u^172 - 3*u^170 - 2*u^168 + u^166 - 4*u^164 - 
#                3*u^158 + 3*u^156 + 4*u^150 - u^148 + 2*u^146 + 3*u^144 - u^
#                142 + 3*u^140 - u^136 + 2*u^134 - 2*u^132 - u^130 - 3*u^
#                126 - u^122 - 2*u^120 + u^118 - u^116 + u^112 + u^108 + u^
#                106 + u^102 - u^94, u^219 - u^211 - u^207 - u^205 - u^201 + u^
#                197 - u^195 + 2*u^193 + u^191 + 3*u^187 + u^183 + 2*u^181 - 
#                2*u^179 + u^177 - 3*u^173 + u^171 - 3*u^169 - 2*u^167 + u^
#                165 - 4*u^163 - 3*u^157 + 3*u^155 + 4*u^149 - u^147 + 2*u^
#                145 + 3*u^143 - u^141 + 3*u^139 - u^135 + 2*u^133 - 2*u^
#                131 - u^129 - 3*u^125 - u^121 - 2*u^119 + u^117 - u^115 + u^
#                111 + u^107 + u^105 + u^101 - u^93, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, u^219 - u^211 - u^207 - u^205 - u^201 + u^197 - u^195 + 
#                2*u^193 + u^191 + 3*u^187 + u^183 + 2*u^181 - 2*u^179 + u^
#                177 - 3*u^173 + u^171 - 3*u^169 - 2*u^167 + u^165 - 4*u^163 - 
#                3*u^157 + 3*u^155 + 4*u^149 - u^147 + 2*u^145 + 3*u^143 - u^
#                141 + 3*u^139 - u^135 + 2*u^133 - 2*u^131 - u^129 - 3*u^
#                125 - u^121 - 2*u^119 + u^117 - u^115 + u^111 + u^107 + u^
#                105 + u^101 - u^93, u^220 - u^212 - u^208 - u^206 - u^202 + u^
#                198 - u^196 + 2*u^194 + u^192 + 3*u^188 + u^184 + 2*u^182 - 
#                2*u^180 + u^178 - 3*u^174 + u^172 - 3*u^170 - 2*u^168 + u^
#                166 - 4*u^164 - 3*u^158 + 3*u^156 + 4*u^150 - u^148 + 2*u^
#                146 + 3*u^144 - u^142 + 3*u^140 - u^136 + 2*u^134 - 2*u^
#                132 - u^130 - 3*u^126 - u^122 - 2*u^120 + u^118 - u^116 + u^
#                112 + u^108 + u^106 + u^102 - u^94, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ u^218 + u^214 - u^206 - u^204 - u^202 - 2*u^200 - u^198 - u^196 - 
#                2*u^194 + u^192 - u^190 + u^188 + 2*u^186 + u^184 + 3*u^182 + 
#                3*u^180 + u^178 + 4*u^176 + u^174 + u^172 + 2*u^170 - 2*u^
#                168 - u^164 - 4*u^162 - u^160 - 4*u^158 - 4*u^156 - u^154 - 
#                4*u^152 - u^150 - 2*u^146 + 2*u^144 + u^142 + u^140 + 4*u^
#                138 + u^136 + 3*u^134 + 3*u^132 + u^130 + 2*u^128 + u^126 - u^
#                124 + u^122 - 2*u^120 - u^118 - u^116 - 2*u^114 - u^112 - u^
#                110 - u^108 + u^100 + u^96, 0*u^0, 0*u^0, 
#              u^220 + u^216 - u^208 - u^206 - u^204 - 2*u^202 - u^200 - u^
#                198 - 2*u^196 + u^194 - u^192 + u^190 + 2*u^188 + u^186 + 3*u^
#                184 + 3*u^182 + u^180 + 4*u^178 + u^176 + u^174 + 2*u^172 - 
#                2*u^170 - u^166 - 4*u^164 - u^162 - 4*u^160 - 4*u^158 - u^
#                156 - 4*u^154 - u^152 - 2*u^148 + 2*u^146 + u^144 + u^142 + 
#                4*u^140 + u^138 + 3*u^136 + 3*u^134 + u^132 + 2*u^130 + u^
#                128 - u^126 + u^124 - 2*u^122 - u^120 - u^118 - 2*u^116 - u^
#                114 - u^112 - u^110 + u^102 + u^98, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^220 - u^212 - u^208 - u^206 - u^
#                202 + u^198 - u^196 + 2*u^194 + u^192 + 3*u^188 + u^184 + 2*u^
#                182 - 2*u^180 + u^178 - 3*u^174 + u^172 - 3*u^170 - 2*u^
#                168 + u^166 - 4*u^164 - 3*u^158 + 3*u^156 + 4*u^150 - u^148 + 
#                2*u^146 + 3*u^144 - u^142 + 3*u^140 - u^136 + 2*u^134 - 2*u^
#                132 - u^130 - 3*u^126 - u^122 - 2*u^120 + u^118 - u^116 + u^
#                112 + u^108 + u^106 + u^102 - u^94, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#              u^220 - u^218 - u^212 + u^210 - u^208 + u^204 - u^202 + u^
#                200 + u^198 - 2*u^196 + 3*u^194 - u^192 - u^190 + 3*u^188 - 
#                3*u^186 + u^184 + u^182 - 4*u^180 + 3*u^178 - u^176 - 3*u^
#                174 + 4*u^172 - 4*u^170 + u^168 + 3*u^166 - 5*u^164 + 4*u^
#                162 - 3*u^158 + 6*u^156 - 3*u^154 + 4*u^150 - 5*u^148 + 3*u^
#                146 + u^144 - 4*u^142 + 4*u^140 - 3*u^138 - u^136 + 3*u^134 - 
#                4*u^132 + u^130 + u^128 - 3*u^126 + 3*u^124 - u^122 - u^120 + 
#                3*u^118 - 2*u^116 + u^114 + u^112 - u^110 + u^108 - u^104 + u^
#                102 - u^100 - u^94 + u^92 ] ], 
#      [ [ u^222 - u^214 - u^210 - u^208 - u^204 + u^200 - u^198 + 2*u^196 + u^
#                194 + 3*u^190 + u^186 + 2*u^184 - 2*u^182 + u^180 - 3*u^
#                176 + u^174 - 3*u^172 - 2*u^170 + u^168 - 4*u^166 - 3*u^160 + 
#                3*u^158 + 4*u^152 - u^150 + 2*u^148 + 3*u^146 - u^144 + 3*u^
#                142 - u^138 + 2*u^136 - 2*u^134 - u^132 - 3*u^128 - u^124 - 
#                2*u^122 + u^120 - u^118 + u^114 + u^110 + u^108 + u^104 - u^
#                96 ] ], 
#      [ [ u^224 - u^222 - u^216 + u^214 - u^212 + u^208 - u^206 + u^204 + u^
#                202 - 2*u^200 + 3*u^198 - u^196 - u^194 + 3*u^192 - 3*u^
#                190 + u^188 + u^186 - 4*u^184 + 3*u^182 - u^180 - 3*u^178 + 
#                4*u^176 - 4*u^174 + u^172 + 3*u^170 - 5*u^168 + 4*u^166 - 3*u^
#                162 + 6*u^160 - 3*u^158 + 4*u^154 - 5*u^152 + 3*u^150 + u^
#                148 - 4*u^146 + 4*u^144 - 3*u^142 - u^140 + 3*u^138 - 4*u^
#                136 + u^134 + u^132 - 3*u^130 + 3*u^128 - u^126 - u^124 + 3*u^
#                122 - 2*u^120 + u^118 + u^116 - u^114 + u^112 - u^108 + u^
#                106 - u^104 - u^98 + u^96, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, u^224 - u^222 - u^216 + u^214 - u^212 + u^208 - u^206 + u^
#                204 + u^202 - 2*u^200 + 3*u^198 - u^196 - u^194 + 3*u^192 - 
#                3*u^190 + u^188 + u^186 - 4*u^184 + 3*u^182 - u^180 - 3*u^
#                178 + 4*u^176 - 4*u^174 + u^172 + 3*u^170 - 5*u^168 + 4*u^
#                166 - 3*u^162 + 6*u^160 - 3*u^158 + 4*u^154 - 5*u^152 + 3*u^
#                150 + u^148 - 4*u^146 + 4*u^144 - 3*u^142 - u^140 + 3*u^138 - 
#                4*u^136 + u^134 + u^132 - 3*u^130 + 3*u^128 - u^126 - u^124 + 
#                3*u^122 - 2*u^120 + u^118 + u^116 - u^114 + u^112 - u^108 + u^
#                106 - u^104 - u^98 + u^96, 0*u^0, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, u^224 - u^222 - u^216 + u^214 - u^212 + u^208 - u^
#                206 + u^204 + u^202 - 2*u^200 + 3*u^198 - u^196 - u^194 + 3*u^
#                192 - 3*u^190 + u^188 + u^186 - 4*u^184 + 3*u^182 - u^180 - 
#                3*u^178 + 4*u^176 - 4*u^174 + u^172 + 3*u^170 - 5*u^168 + 4*u^
#                166 - 3*u^162 + 6*u^160 - 3*u^158 + 4*u^154 - 5*u^152 + 3*u^
#                150 + u^148 - 4*u^146 + 4*u^144 - 3*u^142 - u^140 + 3*u^138 - 
#                4*u^136 + u^134 + u^132 - 3*u^130 + 3*u^128 - u^126 - u^124 + 
#                3*u^122 - 2*u^120 + u^118 + u^116 - u^114 + u^112 - u^108 + u^
#                106 - u^104 - u^98 + u^96, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, u^224 - u^222 - u^216 + u^214 - u^212 + u^
#                208 - u^206 + u^204 + u^202 - 2*u^200 + 3*u^198 - u^196 - u^
#                194 + 3*u^192 - 3*u^190 + u^188 + u^186 - 4*u^184 + 3*u^
#                182 - u^180 - 3*u^178 + 4*u^176 - 4*u^174 + u^172 + 3*u^170 - 
#                5*u^168 + 4*u^166 - 3*u^162 + 6*u^160 - 3*u^158 + 4*u^154 - 
#                5*u^152 + 3*u^150 + u^148 - 4*u^146 + 4*u^144 - 3*u^142 - u^
#                140 + 3*u^138 - 4*u^136 + u^134 + u^132 - 3*u^130 + 3*u^
#                128 - u^126 - u^124 + 3*u^122 - 2*u^120 + u^118 + u^116 - u^
#                114 + u^112 - u^108 + u^106 - u^104 - u^98 + u^96, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^224 - u^222 - u^216 + u^214 - u^
#                212 + u^208 - u^206 + u^204 + u^202 - 2*u^200 + 3*u^198 - u^
#                196 - u^194 + 3*u^192 - 3*u^190 + u^188 + u^186 - 4*u^184 + 
#                3*u^182 - u^180 - 3*u^178 + 4*u^176 - 4*u^174 + u^172 + 3*u^
#                170 - 5*u^168 + 4*u^166 - 3*u^162 + 6*u^160 - 3*u^158 + 4*u^
#                154 - 5*u^152 + 3*u^150 + u^148 - 4*u^146 + 4*u^144 - 3*u^
#                142 - u^140 + 3*u^138 - 4*u^136 + u^134 + u^132 - 3*u^130 + 
#                3*u^128 - u^126 - u^124 + 3*u^122 - 2*u^120 + u^118 + u^
#                116 - u^114 + u^112 - u^108 + u^106 - u^104 - u^98 + u^96 ] ],
#      [ [ u^226 + u^222 - u^218 - 2*u^214 - u^212 - u^210 - 2*u^208 - u^202 + 
#                3*u^200 + 2*u^196 + 4*u^194 + 4*u^190 + 2*u^188 - u^186 + 3*u^
#                184 - 2*u^182 - 2*u^180 + u^178 - 6*u^176 - u^174 - 2*u^172 - 
#                6*u^170 + u^168 - 4*u^166 - 3*u^164 + 3*u^162 - 3*u^160 + 3*u^
#                158 + 4*u^156 - u^154 + 6*u^152 + 2*u^150 + u^148 + 6*u^
#                146 - u^144 + 2*u^142 + 2*u^140 - 3*u^138 + u^136 - 2*u^134 - 
#                4*u^132 - 4*u^128 - 2*u^126 - 3*u^122 + u^120 + 2*u^114 + u^
#                112 + u^110 + 2*u^108 + u^104 - u^100 - u^96, 
#              u^225 + u^223 - u^217 - u^215 - u^213 - 2*u^211 - u^209 - u^
#                207 - u^205 + u^203 + u^199 + 3*u^197 + u^195 + 3*u^193 + 3*u^
#                191 + u^189 + 3*u^187 - u^183 + u^181 - 3*u^179 - 2*u^177 - 
#                2*u^175 - 5*u^173 - u^171 - 3*u^169 - 4*u^167 - 3*u^163 + 3*u^
#                159 + 4*u^155 + 3*u^153 + u^151 + 5*u^149 + 2*u^147 + 2*u^
#                145 + 3*u^143 - u^141 + u^139 - 3*u^135 - u^133 - 3*u^131 - 
#                3*u^129 - u^127 - 3*u^125 - u^123 - u^119 + u^117 + u^115 + u^
#                113 + 2*u^111 + u^109 + u^107 + u^105 - u^99 - u^97, 
#              u^224 - u^216 - u^212 - u^210 - u^206 + u^202 - u^200 + 2*u^
#                198 + u^196 + 3*u^192 + u^188 + 2*u^186 - 2*u^184 + u^182 - 
#                3*u^178 + u^176 - 3*u^174 - 2*u^172 + u^170 - 4*u^168 - 3*u^
#                162 + 3*u^160 + 4*u^154 - u^152 + 2*u^150 + 3*u^148 - u^146 + 
#                3*u^144 - u^140 + 2*u^138 - 2*u^136 - u^134 - 3*u^130 - u^
#                126 - 2*u^124 + u^122 - u^120 + u^116 + u^112 + u^110 + u^
#                106 - u^98, u^224 - u^216 - u^212 - u^210 - u^206 + u^202 - u^
#                200 + 2*u^198 + u^196 + 3*u^192 + u^188 + 2*u^186 - 2*u^
#                184 + u^182 - 3*u^178 + u^176 - 3*u^174 - 2*u^172 + u^170 - 
#                4*u^168 - 3*u^162 + 3*u^160 + 4*u^154 - u^152 + 2*u^150 + 3*u^
#                148 - u^146 + 3*u^144 - u^140 + 2*u^138 - 2*u^136 - u^134 - 
#                3*u^130 - u^126 - 2*u^124 + u^122 - u^120 + u^116 + u^112 + u^
#                110 + u^106 - u^98, 0*u^0 ], 
#          [ u^225 + u^223 - u^217 - u^215 - u^213 - 2*u^211 - u^209 - u^
#                207 - u^205 + u^203 + u^199 + 3*u^197 + u^195 + 3*u^193 + 3*u^
#                191 + u^189 + 3*u^187 - u^183 + u^181 - 3*u^179 - 2*u^177 - 
#                2*u^175 - 5*u^173 - u^171 - 3*u^169 - 4*u^167 - 3*u^163 + 3*u^
#                159 + 4*u^155 + 3*u^153 + u^151 + 5*u^149 + 2*u^147 + 2*u^
#                145 + 3*u^143 - u^141 + u^139 - 3*u^135 - u^133 - 3*u^131 - 
#                3*u^129 - u^127 - 3*u^125 - u^123 - u^119 + u^117 + u^115 + u^
#                113 + 2*u^111 + u^109 + u^107 + u^105 - u^99 - u^97, 
#              u^226 + u^224 - u^218 - u^216 - u^214 - 2*u^212 - u^210 - u^
#                208 - u^206 + u^204 + u^200 + 3*u^198 + u^196 + 3*u^194 + 3*u^
#                192 + u^190 + 3*u^188 - u^184 + u^182 - 3*u^180 - 2*u^178 - 
#                2*u^176 - 5*u^174 - u^172 - 3*u^170 - 4*u^168 - 3*u^164 + 3*u^
#                160 + 4*u^156 + 3*u^154 + u^152 + 5*u^150 + 2*u^148 + 2*u^
#                146 + 3*u^144 - u^142 + u^140 - 3*u^136 - u^134 - 3*u^132 - 
#                3*u^130 - u^128 - 3*u^126 - u^124 - u^120 + u^118 + u^116 + u^
#                114 + 2*u^112 + u^110 + u^108 + u^106 - u^100 - u^98, 
#              u^225 - u^217 - u^213 - u^211 - u^207 + u^203 - u^201 + 2*u^
#                199 + u^197 + 3*u^193 + u^189 + 2*u^187 - 2*u^185 + u^183 - 
#                3*u^179 + u^177 - 3*u^175 - 2*u^173 + u^171 - 4*u^169 - 3*u^
#                163 + 3*u^161 + 4*u^155 - u^153 + 2*u^151 + 3*u^149 - u^147 + 
#                3*u^145 - u^141 + 2*u^139 - 2*u^137 - u^135 - 3*u^131 - u^
#                127 - 2*u^125 + u^123 - u^121 + u^117 + u^113 + u^111 + u^
#                107 - u^99, u^225 - u^217 - u^213 - u^211 - u^207 + u^203 - u^
#                201 + 2*u^199 + u^197 + 3*u^193 + u^189 + 2*u^187 - 2*u^
#                185 + u^183 - 3*u^179 + u^177 - 3*u^175 - 2*u^173 + u^171 - 
#                4*u^169 - 3*u^163 + 3*u^161 + 4*u^155 - u^153 + 2*u^151 + 3*u^
#                149 - u^147 + 3*u^145 - u^141 + 2*u^139 - 2*u^137 - u^135 - 
#                3*u^131 - u^127 - 2*u^125 + u^123 - u^121 + u^117 + u^113 + u^
#                111 + u^107 - u^99, 0*u^0 ], 
#          [ u^224 - u^216 - u^212 - u^210 - u^206 + u^202 - u^200 + 2*u^
#                198 + u^196 + 3*u^192 + u^188 + 2*u^186 - 2*u^184 + u^182 - 
#                3*u^178 + u^176 - 3*u^174 - 2*u^172 + u^170 - 4*u^168 - 3*u^
#                162 + 3*u^160 + 4*u^154 - u^152 + 2*u^150 + 3*u^148 - u^146 + 
#                3*u^144 - u^140 + 2*u^138 - 2*u^136 - u^134 - 3*u^130 - u^
#                126 - 2*u^124 + u^122 - u^120 + u^116 + u^112 + u^110 + u^
#                106 - u^98, u^225 - u^217 - u^213 - u^211 - u^207 + u^203 - u^
#                201 + 2*u^199 + u^197 + 3*u^193 + u^189 + 2*u^187 - 2*u^
#                185 + u^183 - 3*u^179 + u^177 - 3*u^175 - 2*u^173 + u^171 - 
#                4*u^169 - 3*u^163 + 3*u^161 + 4*u^155 - u^153 + 2*u^151 + 3*u^
#                149 - u^147 + 3*u^145 - u^141 + 2*u^139 - 2*u^137 - u^135 - 
#                3*u^131 - u^127 - 2*u^125 + u^123 - u^121 + u^117 + u^113 + u^
#                111 + u^107 - u^99, u^226 - u^218 - u^214 - u^212 - u^208 + u^
#                204 - u^202 + 2*u^200 + u^198 + 3*u^194 + u^190 + 2*u^188 - 
#                2*u^186 + u^184 - 3*u^180 + u^178 - 3*u^176 - 2*u^174 + u^
#                172 - 4*u^170 - 3*u^164 + 3*u^162 + 4*u^156 - u^154 + 2*u^
#                152 + 3*u^150 - u^148 + 3*u^146 - u^142 + 2*u^140 - 2*u^
#                138 - u^136 - 3*u^132 - u^128 - 2*u^126 + u^124 - u^122 + u^
#                118 + u^114 + u^112 + u^108 - u^100, 0*u^0, 0*u^0 ], 
#          [ u^224 - u^216 - u^212 - u^210 - u^206 + u^202 - u^200 + 2*u^
#                198 + u^196 + 3*u^192 + u^188 + 2*u^186 - 2*u^184 + u^182 - 
#                3*u^178 + u^176 - 3*u^174 - 2*u^172 + u^170 - 4*u^168 - 3*u^
#                162 + 3*u^160 + 4*u^154 - u^152 + 2*u^150 + 3*u^148 - u^146 + 
#                3*u^144 - u^140 + 2*u^138 - 2*u^136 - u^134 - 3*u^130 - u^
#                126 - 2*u^124 + u^122 - u^120 + u^116 + u^112 + u^110 + u^
#                106 - u^98, u^225 - u^217 - u^213 - u^211 - u^207 + u^203 - u^
#                201 + 2*u^199 + u^197 + 3*u^193 + u^189 + 2*u^187 - 2*u^
#                185 + u^183 - 3*u^179 + u^177 - 3*u^175 - 2*u^173 + u^171 - 
#                4*u^169 - 3*u^163 + 3*u^161 + 4*u^155 - u^153 + 2*u^151 + 3*u^
#                149 - u^147 + 3*u^145 - u^141 + 2*u^139 - 2*u^137 - u^135 - 
#                3*u^131 - u^127 - 2*u^125 + u^123 - u^121 + u^117 + u^113 + u^
#                111 + u^107 - u^99, 0*u^0, u^226 - u^218 - u^214 - u^212 - u^
#                208 + u^204 - u^202 + 2*u^200 + u^198 + 3*u^194 + u^190 + 2*u^
#                188 - 2*u^186 + u^184 - 3*u^180 + u^178 - 3*u^176 - 2*u^
#                174 + u^172 - 4*u^170 - 3*u^164 + 3*u^162 + 4*u^156 - u^154 + 
#                2*u^152 + 3*u^150 - u^148 + 3*u^146 - u^142 + 2*u^140 - 2*u^
#                138 - u^136 - 3*u^132 - u^128 - 2*u^126 + u^124 - u^122 + u^
#                118 + u^114 + u^112 + u^108 - u^100, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^226 - u^224 - u^218 + u^216 - u^
#                214 + u^210 - u^208 + u^206 + u^204 - 2*u^202 + 3*u^200 - u^
#                198 - u^196 + 3*u^194 - 3*u^192 + u^190 + u^188 - 4*u^186 + 
#                3*u^184 - u^182 - 3*u^180 + 4*u^178 - 4*u^176 + u^174 + 3*u^
#                172 - 5*u^170 + 4*u^168 - 3*u^164 + 6*u^162 - 3*u^160 + 4*u^
#                156 - 5*u^154 + 3*u^152 + u^150 - 4*u^148 + 4*u^146 - 3*u^
#                144 - u^142 + 3*u^140 - 4*u^138 + u^136 + u^134 - 3*u^132 + 
#                3*u^130 - u^128 - u^126 + 3*u^124 - 2*u^122 + u^120 + u^
#                118 - u^116 + u^114 - u^110 + u^108 - u^106 - u^100 + u^98 ] ]
#        , 
#      [ [ u^228 - u^220 - u^216 - u^214 - u^210 + u^206 - u^204 + 2*u^202 + u^
#                200 + 3*u^196 + u^192 + 2*u^190 - 2*u^188 + u^186 - 3*u^
#                182 + u^180 - 3*u^178 - 2*u^176 + u^174 - 4*u^172 - 3*u^166 + 
#                3*u^164 + 4*u^158 - u^156 + 2*u^154 + 3*u^152 - u^150 + 3*u^
#                148 - u^144 + 2*u^142 - 2*u^140 - u^138 - 3*u^134 - u^130 - 
#                2*u^128 + u^126 - u^124 + u^120 + u^116 + u^114 + u^110 - u^
#                102, 0*u^0, u^227 - u^219 - u^215 - u^213 - u^209 + u^205 - u^
#                203 + 2*u^201 + u^199 + 3*u^195 + u^191 + 2*u^189 - 2*u^
#                187 + u^185 - 3*u^181 + u^179 - 3*u^177 - 2*u^175 + u^173 - 
#                4*u^171 - 3*u^165 + 3*u^163 + 4*u^157 - u^155 + 2*u^153 + 3*u^
#                151 - u^149 + 3*u^147 - u^143 + 2*u^141 - 2*u^139 - u^137 - 
#                3*u^133 - u^129 - 2*u^127 + u^125 - u^123 + u^119 + u^115 + u^
#                113 + u^109 - u^101, 0*u^0 ], 
#          [ 0*u^0, u^228 - u^220 - u^216 - u^214 - u^210 + u^206 - u^204 + 
#                2*u^202 + u^200 + 3*u^196 + u^192 + 2*u^190 - 2*u^188 + u^
#                186 - 3*u^182 + u^180 - 3*u^178 - 2*u^176 + u^174 - 4*u^172 - 
#                3*u^166 + 3*u^164 + 4*u^158 - u^156 + 2*u^154 + 3*u^152 - u^
#                150 + 3*u^148 - u^144 + 2*u^142 - 2*u^140 - u^138 - 3*u^
#                134 - u^130 - 2*u^128 + u^126 - u^124 + u^120 + u^116 + u^
#                114 + u^110 - u^102, 0*u^0, 0*u^0 ], 
#          [ u^227 - u^219 - u^215 - u^213 - u^209 + u^205 - u^203 + 2*u^
#                201 + u^199 + 3*u^195 + u^191 + 2*u^189 - 2*u^187 + u^185 - 
#                3*u^181 + u^179 - 3*u^177 - 2*u^175 + u^173 - 4*u^171 - 3*u^
#                165 + 3*u^163 + 4*u^157 - u^155 + 2*u^153 + 3*u^151 - u^149 + 
#                3*u^147 - u^143 + 2*u^141 - 2*u^139 - u^137 - 3*u^133 - u^
#                129 - 2*u^127 + u^125 - u^123 + u^119 + u^115 + u^113 + u^
#                109 - u^101, 0*u^0, u^228 - u^220 - u^216 - u^214 - u^210 + u^
#                206 - u^204 + 2*u^202 + u^200 + 3*u^196 + u^192 + 2*u^190 - 
#                2*u^188 + u^186 - 3*u^182 + u^180 - 3*u^178 - 2*u^176 + u^
#                174 - 4*u^172 - 3*u^166 + 3*u^164 + 4*u^158 - u^156 + 2*u^
#                154 + 3*u^152 - u^150 + 3*u^148 - u^144 + 2*u^142 - 2*u^
#                140 - u^138 - 3*u^134 - u^130 - 2*u^128 + u^126 - u^124 + u^
#                120 + u^116 + u^114 + u^110 - u^102, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, 0*u^0, u^228 - u^226 - u^220 + u^218 - u^216 + u^
#                212 - u^210 + u^208 + u^206 - 2*u^204 + 3*u^202 - u^200 - u^
#                198 + 3*u^196 - 3*u^194 + u^192 + u^190 - 4*u^188 + 3*u^
#                186 - u^184 - 3*u^182 + 4*u^180 - 4*u^178 + u^176 + 3*u^174 - 
#                5*u^172 + 4*u^170 - 3*u^166 + 6*u^164 - 3*u^162 + 4*u^158 - 
#                5*u^156 + 3*u^154 + u^152 - 4*u^150 + 4*u^148 - 3*u^146 - u^
#                144 + 3*u^142 - 4*u^140 + u^138 + u^136 - 3*u^134 + 3*u^
#                132 - u^130 - u^128 + 3*u^126 - 2*u^124 + u^122 + u^120 - u^
#                118 + u^116 - u^112 + u^110 - u^108 - u^102 + u^100 ] ], 
#      [ [ u^230 - u^228 - u^222 + u^220 - u^218 + u^214 - u^212 + u^210 + u^
#                208 - 2*u^206 + 3*u^204 - u^202 - u^200 + 3*u^198 - 3*u^
#                196 + u^194 + u^192 - 4*u^190 + 3*u^188 - u^186 - 3*u^184 + 
#                4*u^182 - 4*u^180 + u^178 + 3*u^176 - 5*u^174 + 4*u^172 - 3*u^
#                168 + 6*u^166 - 3*u^164 + 4*u^160 - 5*u^158 + 3*u^156 + u^
#                154 - 4*u^152 + 4*u^150 - 3*u^148 - u^146 + 3*u^144 - 4*u^
#                142 + u^140 + u^138 - 3*u^136 + 3*u^134 - u^132 - u^130 + 3*u^
#                128 - 2*u^126 + u^124 + u^122 - u^120 + u^118 - u^114 + u^
#                112 - u^110 - u^104 + u^102 ] ], 
#      [ [ u^232 - u^230 - u^224 + u^222 - u^220 + u^216 - u^214 + u^212 + u^
#                210 - 2*u^208 + 3*u^206 - u^204 - u^202 + 3*u^200 - 3*u^
#                198 + u^196 + u^194 - 4*u^192 + 3*u^190 - u^188 - 3*u^186 + 
#                4*u^184 - 4*u^182 + u^180 + 3*u^178 - 5*u^176 + 4*u^174 - 3*u^
#                170 + 6*u^168 - 3*u^166 + 4*u^162 - 5*u^160 + 3*u^158 + u^
#                156 - 4*u^154 + 4*u^152 - 3*u^150 - u^148 + 3*u^146 - 4*u^
#                144 + u^142 + u^140 - 3*u^138 + 3*u^136 - u^134 - u^132 + 3*u^
#                130 - 2*u^128 + u^126 + u^124 - u^122 + u^120 - u^116 + u^
#                114 - u^112 - u^106 + u^104, 0*u^0, 0*u^0 ], 
#          [ 0*u^0, u^232 - u^230 - u^224 + u^222 - u^220 + u^216 - u^214 + u^
#                212 + u^210 - 2*u^208 + 3*u^206 - u^204 - u^202 + 3*u^200 - 
#                3*u^198 + u^196 + u^194 - 4*u^192 + 3*u^190 - u^188 - 3*u^
#                186 + 4*u^184 - 4*u^182 + u^180 + 3*u^178 - 5*u^176 + 4*u^
#                174 - 3*u^170 + 6*u^168 - 3*u^166 + 4*u^162 - 5*u^160 + 3*u^
#                158 + u^156 - 4*u^154 + 4*u^152 - 3*u^150 - u^148 + 3*u^146 - 
#                4*u^144 + u^142 + u^140 - 3*u^138 + 3*u^136 - u^134 - u^132 + 
#                3*u^130 - 2*u^128 + u^126 + u^124 - u^122 + u^120 - u^116 + u^
#                114 - u^112 - u^106 + u^104, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, u^232 - u^230 - u^224 + u^222 - u^220 + u^216 - u^
#                214 + u^212 + u^210 - 2*u^208 + 3*u^206 - u^204 - u^202 + 3*u^
#                200 - 3*u^198 + u^196 + u^194 - 4*u^192 + 3*u^190 - u^188 - 
#                3*u^186 + 4*u^184 - 4*u^182 + u^180 + 3*u^178 - 5*u^176 + 4*u^
#                174 - 3*u^170 + 6*u^168 - 3*u^166 + 4*u^162 - 5*u^160 + 3*u^
#                158 + u^156 - 4*u^154 + 4*u^152 - 3*u^150 - u^148 + 3*u^146 - 
#                4*u^144 + u^142 + u^140 - 3*u^138 + 3*u^136 - u^134 - u^132 + 
#                3*u^130 - 2*u^128 + u^126 + u^124 - u^122 + u^120 - u^116 + u^
#                114 - u^112 - u^106 + u^104 ] ], 
#      [ [ u^234 - u^226 - u^222 - u^220 - u^216 + u^212 - u^210 + 2*u^208 + u^
#                206 + 3*u^202 + u^198 + 2*u^196 - 2*u^194 + u^192 - 3*u^
#                188 + u^186 - 3*u^184 - 2*u^182 + u^180 - 4*u^178 - 3*u^172 + 
#                3*u^170 + 4*u^164 - u^162 + 2*u^160 + 3*u^158 - u^156 + 3*u^
#                154 - u^150 + 2*u^148 - 2*u^146 - u^144 - 3*u^140 - u^136 - 
#                2*u^134 + u^132 - u^130 + u^126 + u^122 + u^120 + u^116 - u^
#                108, u^233 - u^225 - u^221 - u^219 - u^215 + u^211 - u^209 + 
#                2*u^207 + u^205 + 3*u^201 + u^197 + 2*u^195 - 2*u^193 + u^
#                191 - 3*u^187 + u^185 - 3*u^183 - 2*u^181 + u^179 - 4*u^177 - 
#                3*u^171 + 3*u^169 + 4*u^163 - u^161 + 2*u^159 + 3*u^157 - u^
#                155 + 3*u^153 - u^149 + 2*u^147 - 2*u^145 - u^143 - 3*u^
#                139 - u^135 - 2*u^133 + u^131 - u^129 + u^125 + u^121 + u^
#                119 + u^115 - u^107, 0*u^0 ], 
#          [ u^233 - u^225 - u^221 - u^219 - u^215 + u^211 - u^209 + 2*u^
#                207 + u^205 + 3*u^201 + u^197 + 2*u^195 - 2*u^193 + u^191 - 
#                3*u^187 + u^185 - 3*u^183 - 2*u^181 + u^179 - 4*u^177 - 3*u^
#                171 + 3*u^169 + 4*u^163 - u^161 + 2*u^159 + 3*u^157 - u^155 + 
#                3*u^153 - u^149 + 2*u^147 - 2*u^145 - u^143 - 3*u^139 - u^
#                135 - 2*u^133 + u^131 - u^129 + u^125 + u^121 + u^119 + u^
#                115 - u^107, u^234 - u^226 - u^222 - u^220 - u^216 + u^
#                212 - u^210 + 2*u^208 + u^206 + 3*u^202 + u^198 + 2*u^196 - 
#                2*u^194 + u^192 - 3*u^188 + u^186 - 3*u^184 - 2*u^182 + u^
#                180 - 4*u^178 - 3*u^172 + 3*u^170 + 4*u^164 - u^162 + 2*u^
#                160 + 3*u^158 - u^156 + 3*u^154 - u^150 + 2*u^148 - 2*u^
#                146 - u^144 - 3*u^140 - u^136 - 2*u^134 + u^132 - u^130 + u^
#                126 + u^122 + u^120 + u^116 - u^108, 0*u^0 ], 
#          [ 0*u^0, 0*u^0, u^234 - u^232 - u^226 + u^224 - u^222 + u^218 - u^
#                216 + u^214 + u^212 - 2*u^210 + 3*u^208 - u^206 - u^204 + 3*u^
#                202 - 3*u^200 + u^198 + u^196 - 4*u^194 + 3*u^192 - u^190 - 
#                3*u^188 + 4*u^186 - 4*u^184 + u^182 + 3*u^180 - 5*u^178 + 4*u^
#                176 - 3*u^172 + 6*u^170 - 3*u^168 + 4*u^164 - 5*u^162 + 3*u^
#                160 + u^158 - 4*u^156 + 4*u^154 - 3*u^152 - u^150 + 3*u^148 - 
#                4*u^146 + u^144 + u^142 - 3*u^140 + 3*u^138 - u^136 - u^134 + 
#                3*u^132 - 2*u^130 + u^128 + u^126 - u^124 + u^122 - u^118 + u^
#                116 - u^114 - u^108 + u^106 ] ], 
#      [ [ u^236 - u^234 - u^228 + u^226 - u^224 + u^220 - u^218 + u^216 + u^
#                214 - 2*u^212 + 3*u^210 - u^208 - u^206 + 3*u^204 - 3*u^
#                202 + u^200 + u^198 - 4*u^196 + 3*u^194 - u^192 - 3*u^190 + 
#                4*u^188 - 4*u^186 + u^184 + 3*u^182 - 5*u^180 + 4*u^178 - 3*u^
#                174 + 6*u^172 - 3*u^170 + 4*u^166 - 5*u^164 + 3*u^162 + u^
#                160 - 4*u^158 + 4*u^156 - 3*u^154 - u^152 + 3*u^150 - 4*u^
#                148 + u^146 + u^144 - 3*u^142 + 3*u^140 - u^138 - u^136 + 3*u^
#                134 - 2*u^132 + u^130 + u^128 - u^126 + u^124 - u^120 + u^
#                118 - u^116 - u^110 + u^108 ] ], 
#      [ [ u^238 - u^236 - u^230 + u^228 - u^226 + u^222 - u^220 + u^218 + u^
#                216 - 2*u^214 + 3*u^212 - u^210 - u^208 + 3*u^206 - 3*u^
#                204 + u^202 + u^200 - 4*u^198 + 3*u^196 - u^194 - 3*u^192 + 
#                4*u^190 - 4*u^188 + u^186 + 3*u^184 - 5*u^182 + 4*u^180 - 3*u^
#                176 + 6*u^174 - 3*u^172 + 4*u^168 - 5*u^166 + 3*u^164 + u^
#                162 - 4*u^160 + 4*u^158 - 3*u^156 - u^154 + 3*u^152 - 4*u^
#                150 + u^148 + u^146 - 3*u^144 + 3*u^142 - u^140 - u^138 + 3*u^
#                136 - 2*u^134 + u^132 + u^130 - u^128 + u^126 - u^122 + u^
#                120 - u^118 - u^112 + u^110 ] ], 
#      [ [ u^240 - u^238 - u^232 + u^230 - u^228 + u^224 - u^222 + u^220 + u^
#                218 - 2*u^216 + 3*u^214 - u^212 - u^210 + 3*u^208 - 3*u^
#                206 + u^204 + u^202 - 4*u^200 + 3*u^198 - u^196 - 3*u^194 + 
#                4*u^192 - 4*u^190 + u^188 + 3*u^186 - 5*u^184 + 4*u^182 - 3*u^
#                178 + 6*u^176 - 3*u^174 + 4*u^170 - 5*u^168 + 3*u^166 + u^
#                164 - 4*u^162 + 4*u^160 - 3*u^158 - u^156 + 3*u^154 - 4*u^
#                152 + u^150 + u^148 - 3*u^146 + 3*u^144 - u^142 - u^140 + 3*u^
#                138 - 2*u^136 + u^134 + u^132 - u^130 + u^128 - u^124 + u^
#                122 - u^120 - u^114 + u^112 ] ] ], 
#  [ [ u^120, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^119 + u^113 + u^109 + u^107 + u^103 + u^101 + u^97 + u^91, u^91, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^118 + u^114 + u^112 + u^110 + 2*u^108 + 2*u^106 + u^104 + 3*u^102 + 
#            2*u^100 + 2*u^98 + 3*u^96 + 2*u^94 + 2*u^92 + 3*u^90 + u^88 + 2*u^
#            86 + 2*u^84 + u^82 + u^80 + u^78 + u^74, 
#          u^90 + u^86 + u^84 + u^82 + u^80 + u^78 + u^74, u^74, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^117 + u^115 + u^113 + 2*u^111 + 2*u^109 + 3*u^107 + 4*u^105 + 4*u^
#            103 + 5*u^101 + 6*u^99 + 6*u^97 + 7*u^95 + 7*u^93 + 7*u^91 + 7*u^
#            89 + 7*u^87 + 7*u^85 + 6*u^83 + 6*u^81 + 5*u^79 + 4*u^77 + 4*u^
#            75 + 3*u^73 + 2*u^71 + 2*u^69 + u^67 + u^65 + u^63, 
#          u^89 + u^87 + 2*u^85 + 2*u^83 + 3*u^81 + 3*u^79 + 3*u^77 + 3*u^75 + 
#            3*u^73 + 2*u^71 + 2*u^69 + u^67 + u^65 + u^63, 
#          u^73 + u^71 + u^69 + u^67 + u^65 + u^63, u^63, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^116 + u^114 + u^112 + 2*u^110 + 2*u^108 + 2*u^106 + 4*u^104 + 3*u^
#            102 + 4*u^100 + 5*u^98 + 4*u^96 + 5*u^94 + 6*u^92 + 4*u^90 + 6*u^
#            88 + 5*u^86 + 4*u^84 + 5*u^82 + 4*u^80 + 3*u^78 + 4*u^76 + 2*u^
#            74 + 2*u^72 + 2*u^70 + u^68 + u^66 + u^64, 
#          u^88 + u^86 + u^84 + 2*u^82 + 2*u^80 + 2*u^78 + 3*u^76 + 2*u^74 + 
#            2*u^72 + 2*u^70 + u^68 + u^66 + u^64, 
#          u^72 + u^70 + u^68 + u^66 + u^64, 0*u^0, u^63, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^112 + u^108 + u^106 + 2*u^102 + 2*u^100 + 3*u^96 + u^94 + u^92 + 
#            4*u^90 + u^88 + u^86 + 3*u^84 + 2*u^80 + 2*u^78 + u^74 + u^72 + u^
#            68, u^90 + u^84 + u^80 + u^78 + u^74 + u^72 + u^68, u^68, 0*u^0, 
#          0*u^0, u^63, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^116 + u^112 + 2*u^110 + 2*u^108 + 3*u^106 + 5*u^104 + 4*u^102 + 7*u^
#            100 + 8*u^98 + 8*u^96 + 10*u^94 + 12*u^92 + 10*u^90 + 13*u^88 + 
#            13*u^86 + 12*u^84 + 13*u^82 + 13*u^80 + 10*u^78 + 12*u^76 + 10*u^
#            74 + 8*u^72 + 8*u^70 + 7*u^68 + 4*u^66 + 5*u^64 + 3*u^62 + 2*u^
#            60 + 2*u^58 + u^56 + u^52, u^88 + u^86 + 2*u^84 + 3*u^82 + 4*u^
#            80 + 4*u^78 + 6*u^76 + 6*u^74 + 6*u^72 + 6*u^70 + 6*u^68 + 4*u^
#            66 + 5*u^64 + 3*u^62 + 2*u^60 + 2*u^58 + u^56 + u^52, 
#          u^72 + u^70 + 2*u^68 + 2*u^66 + 3*u^64 + 2*u^62 + 2*u^60 + 2*u^
#            58 + u^56 + u^52, u^62 + u^58 + u^56 + u^52, 0*u^0, u^59 + u^55, 
#          u^52, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^113 + u^111 + u^109 + 3*u^107 + 2*u^105 + 3*u^103 + 6*u^101 + 4*u^
#            99 + 6*u^97 + 9*u^95 + 6*u^93 + 9*u^91 + 11*u^89 + 7*u^87 + 11*u^
#            85 + 11*u^83 + 7*u^81 + 11*u^79 + 9*u^77 + 6*u^75 + 9*u^73 + 6*u^
#            71 + 4*u^69 + 6*u^67 + 3*u^65 + 2*u^63 + 3*u^61 + u^59 + u^57 + u^
#            55, u^89 + 2*u^85 + 2*u^83 + 2*u^81 + 4*u^79 + 4*u^77 + 3*u^75 + 
#            6*u^73 + 4*u^71 + 4*u^69 + 5*u^67 + 3*u^65 + 2*u^63 + 3*u^61 + u^
#            59 + u^57 + u^55, u^73 + u^69 + 2*u^67 + u^65 + 2*u^63 + 2*u^
#            61 + u^59 + u^57 + u^55, u^59 + u^55, 0*u^0, 
#          u^62 + u^58 + u^56 + u^52, 0*u^0, u^52, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^112 + u^108 + u^106 + u^104 + u^102 + 2*u^100 + u^98 + 3*u^96 + 2*u^
#            94 + 2*u^92 + 3*u^90 + 3*u^88 + 2*u^86 + 4*u^84 + 2*u^82 + 3*u^
#            80 + 3*u^78 + 2*u^76 + 2*u^74 + 3*u^72 + u^70 + 2*u^68 + u^66 + u^
#            64 + u^62 + u^60 + u^56, u^84 + u^80 + u^78 + u^76 + u^74 + 2*u^
#            72 + u^70 + 2*u^68 + u^66 + u^64 + u^62 + u^60 + u^56, 
#          u^68 + u^64 + u^62 + u^60 + u^56, 0*u^0, u^59 + u^55, 0*u^0, 0*u^0, 
#          0*u^0, u^52, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^115 + u^113 + 2*u^111 + 4*u^109 + 4*u^107 + 7*u^105 + 10*u^103 + 
#            10*u^101 + 15*u^99 + 18*u^97 + 18*u^95 + 25*u^93 + 26*u^91 + 26*u^
#            89 + 33*u^87 + 31*u^85 + 31*u^83 + 36*u^81 + 31*u^79 + 31*u^77 + 
#            33*u^75 + 26*u^73 + 26*u^71 + 25*u^69 + 18*u^67 + 18*u^65 + 15*u^
#            63 + 10*u^61 + 10*u^59 + 7*u^57 + 4*u^55 + 4*u^53 + 2*u^51 + u^
#            49 + u^47, 2*u^87 + 2*u^85 + 4*u^83 + 7*u^81 + 8*u^79 + 11*u^77 + 
#            15*u^75 + 14*u^73 + 17*u^71 + 18*u^69 + 15*u^67 + 16*u^65 + 14*u^
#            63 + 10*u^61 + 10*u^59 + 7*u^57 + 4*u^55 + 4*u^53 + 2*u^51 + u^
#            49 + u^47, 2*u^71 + 3*u^69 + 4*u^67 + 7*u^65 + 7*u^63 + 7*u^61 + 
#            8*u^59 + 6*u^57 + 4*u^55 + 4*u^53 + 2*u^51 + u^49 + u^47, 
#          u^61 + u^59 + 2*u^57 + 2*u^55 + 2*u^53 + 2*u^51 + u^49 + u^47, 
#          u^62 + u^58 + u^56 + u^52, u^60 + u^58 + u^56 + 2*u^54 + u^52 + u^
#            50 + u^48, u^51 + u^49 + u^47, u^50 + u^48, u^49, u^47, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^114 + u^112 + 2*u^110 + 3*u^108 + 5*u^106 + 7*u^104 + 9*u^102 + 
#            11*u^100 + 15*u^98 + 17*u^96 + 21*u^94 + 23*u^92 + 26*u^90 + 29*u^
#            88 + 31*u^86 + 32*u^84 + 34*u^82 + 33*u^80 + 34*u^78 + 32*u^76 + 
#            31*u^74 + 29*u^72 + 26*u^70 + 23*u^68 + 21*u^66 + 17*u^64 + 15*u^
#            62 + 11*u^60 + 9*u^58 + 7*u^56 + 5*u^54 + 3*u^52 + 2*u^50 + u^
#            48 + u^46, u^88 + 2*u^86 + 3*u^84 + 5*u^82 + 7*u^80 + 10*u^78 + 
#            12*u^76 + 14*u^74 + 16*u^72 + 17*u^70 + 17*u^68 + 17*u^66 + 15*u^
#            64 + 14*u^62 + 11*u^60 + 9*u^58 + 7*u^56 + 5*u^54 + 3*u^52 + 2*u^
#            50 + u^48 + u^46, u^72 + 2*u^70 + 3*u^68 + 5*u^66 + 6*u^64 + 8*u^
#            62 + 7*u^60 + 7*u^58 + 6*u^56 + 5*u^54 + 3*u^52 + 2*u^50 + u^
#            48 + u^46, u^62 + u^60 + 2*u^58 + 2*u^56 + 3*u^54 + 2*u^52 + 2*u^
#            50 + u^48 + u^46, 0*u^0, u^61 + u^59 + 2*u^57 + 2*u^55 + 2*u^53 + 
#            2*u^51 + u^49 + u^47, u^50 + u^48 + u^46, u^51 + u^49 + u^47, 
#          0*u^0, u^46, u^46, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^114 + u^112 + 2*u^110 + 3*u^108 + 5*u^106 + 6*u^104 + 10*u^102 + 
#            11*u^100 + 15*u^98 + 18*u^96 + 22*u^94 + 24*u^92 + 30*u^90 + 31*u^
#            88 + 35*u^86 + 37*u^84 + 39*u^82 + 39*u^80 + 42*u^78 + 39*u^76 + 
#            39*u^74 + 37*u^72 + 35*u^70 + 31*u^68 + 30*u^66 + 24*u^64 + 22*u^
#            62 + 18*u^60 + 15*u^58 + 11*u^56 + 10*u^54 + 6*u^52 + 5*u^50 + 
#            3*u^48 + 2*u^46 + u^44 + u^42, u^86 + 2*u^84 + 4*u^82 + 6*u^80 + 
#            10*u^78 + 12*u^76 + 16*u^74 + 18*u^72 + 21*u^70 + 21*u^68 + 23*u^
#            66 + 20*u^64 + 20*u^62 + 17*u^60 + 15*u^58 + 11*u^56 + 10*u^54 + 
#            6*u^52 + 5*u^50 + 3*u^48 + 2*u^46 + u^44 + u^42, 
#          2*u^70 + 3*u^68 + 6*u^66 + 7*u^64 + 10*u^62 + 10*u^60 + 11*u^58 + 
#            9*u^56 + 9*u^54 + 6*u^52 + 5*u^50 + 3*u^48 + 2*u^46 + u^44 + u^42,
#          u^60 + u^58 + 2*u^56 + 3*u^54 + 3*u^52 + 3*u^50 + 3*u^48 + 2*u^
#            46 + u^44 + u^42, u^61 + u^59 + 2*u^57 + u^55 + 2*u^53 + u^51 + u^
#            49, u^57 + u^55 + 2*u^53 + 2*u^51 + 2*u^49 + u^47 + u^45, 
#          u^50 + 2*u^48 + 2*u^46 + u^44 + u^42, u^49 + u^47 + u^45, 
#          u^50 + u^48 + u^46, u^46 + u^44 + u^42, 0*u^0, u^42, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^113 + u^111 + u^109 + 3*u^107 + 3*u^105 + 4*u^103 + 7*u^101 + 7*u^
#            99 + 9*u^97 + 13*u^95 + 12*u^93 + 15*u^91 + 19*u^89 + 17*u^87 + 
#            21*u^85 + 23*u^83 + 20*u^81 + 24*u^79 + 24*u^77 + 20*u^75 + 23*u^
#            73 + 21*u^71 + 17*u^69 + 19*u^67 + 15*u^65 + 12*u^63 + 13*u^61 + 
#            9*u^59 + 7*u^57 + 7*u^55 + 4*u^53 + 3*u^51 + 3*u^49 + u^47 + u^
#            45 + u^43, u^85 + 2*u^83 + 2*u^81 + 5*u^79 + 6*u^77 + 7*u^75 + 
#            10*u^73 + 11*u^71 + 11*u^69 + 14*u^67 + 12*u^65 + 11*u^63 + 12*u^
#            61 + 9*u^59 + 7*u^57 + 7*u^55 + 4*u^53 + 3*u^51 + 3*u^49 + u^
#            47 + u^45 + u^43, u^69 + 3*u^67 + 3*u^65 + 5*u^63 + 7*u^61 + 6*u^
#            59 + 6*u^57 + 6*u^55 + 4*u^53 + 3*u^51 + 3*u^49 + u^47 + u^45 + u^
#            43, u^59 + u^57 + u^55 + 2*u^53 + u^51 + 2*u^49 + u^47 + u^45 + u^
#            43, u^60 + u^58 + u^56 + 2*u^54 + u^52 + u^50 + u^48, 
#          u^56 + u^52 + u^50 + u^46, u^49 + u^47 + u^45 + u^43, u^46, 
#          u^51 + u^49 + u^47 + u^45, u^45 + u^43, 0*u^0, 0*u^0, u^42, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^112 + 2*u^108 + 2*u^106 + 2*u^104 + 4*u^102 + 6*u^100 + 4*u^98 + 
#            10*u^96 + 9*u^94 + 9*u^92 + 14*u^90 + 14*u^88 + 12*u^86 + 20*u^
#            84 + 15*u^82 + 16*u^80 + 20*u^78 + 16*u^76 + 15*u^74 + 20*u^72 + 
#            12*u^70 + 14*u^68 + 14*u^66 + 9*u^64 + 9*u^62 + 10*u^60 + 4*u^
#            58 + 6*u^56 + 4*u^54 + 2*u^52 + 2*u^50 + 2*u^48 + u^44, 
#          2*u^84 + u^82 + 3*u^80 + 5*u^78 + 5*u^76 + 6*u^74 + 10*u^72 + 7*u^
#            70 + 10*u^68 + 10*u^66 + 8*u^64 + 8*u^62 + 9*u^60 + 4*u^58 + 6*u^
#            56 + 4*u^54 + 2*u^52 + 2*u^50 + 2*u^48 + u^44, 
#          u^72 + 2*u^68 + 2*u^66 + 3*u^64 + 4*u^62 + 5*u^60 + 3*u^58 + 5*u^
#            56 + 3*u^54 + 2*u^52 + 2*u^50 + 2*u^48 + u^44, 
#          u^58 + 2*u^54 + u^52 + u^50 + 2*u^48 + u^44, u^59 + u^55, 
#          u^57 + u^55 + 3*u^51 + u^47 + u^45, u^48 + u^44, u^51 + u^47 + u^45,
#          0*u^0, u^44, 0*u^0, 0*u^0, 0*u^0, u^42, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^113 + u^111 + 2*u^109 + 5*u^107 + 6*u^105 + 9*u^103 + 15*u^101 + 
#            17*u^99 + 23*u^97 + 32*u^95 + 34*u^93 + 43*u^91 + 53*u^89 + 54*u^
#            87 + 64*u^85 + 72*u^83 + 70*u^81 + 79*u^79 + 82*u^77 + 76*u^75 + 
#            82*u^73 + 79*u^71 + 70*u^69 + 72*u^67 + 64*u^65 + 54*u^63 + 53*u^
#            61 + 43*u^59 + 34*u^57 + 32*u^55 + 23*u^53 + 17*u^51 + 15*u^49 + 
#            9*u^47 + 6*u^45 + 5*u^43 + 2*u^41 + u^39 + u^37, 
#          2*u^85 + 4*u^83 + 6*u^81 + 12*u^79 + 17*u^77 + 21*u^75 + 30*u^73 + 
#            35*u^71 + 38*u^69 + 45*u^67 + 45*u^65 + 43*u^63 + 45*u^61 + 39*u^
#            59 + 33*u^57 + 31*u^55 + 23*u^53 + 17*u^51 + 15*u^49 + 9*u^47 + 
#            6*u^45 + 5*u^43 + 2*u^41 + u^39 + u^37, 
#          u^71 + 3*u^69 + 6*u^67 + 9*u^65 + 14*u^63 + 19*u^61 + 20*u^59 + 
#            22*u^57 + 23*u^55 + 19*u^53 + 16*u^51 + 14*u^49 + 9*u^47 + 6*u^
#            45 + 5*u^43 + 2*u^41 + u^39 + u^37, 
#          u^61 + 2*u^59 + 3*u^57 + 5*u^55 + 8*u^53 + 7*u^51 + 9*u^49 + 8*u^
#            47 + 5*u^45 + 5*u^43 + 2*u^41 + u^39 + u^37, 
#          u^60 + u^58 + u^56 + 2*u^54 + u^52 + u^50 + u^48, 
#          u^58 + 3*u^56 + 3*u^54 + 5*u^52 + 7*u^50 + 4*u^48 + 5*u^46 + 3*u^
#            44 + u^42 + u^40, 2*u^49 + 4*u^47 + 4*u^45 + 4*u^43 + 2*u^41 + u^
#            39 + u^37, u^50 + 2*u^48 + 4*u^46 + 3*u^44 + u^42 + u^40, 
#          u^47 + u^45, 2*u^45 + 3*u^43 + 2*u^41 + u^39 + u^37, 
#          u^45 + u^43 + u^41 + u^39 + u^37, u^41 + u^37, 0*u^0, u^41 + u^37, 
#          u^37, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^112 + 2*u^110 + 3*u^108 + 5*u^106 + 8*u^104 + 11*u^102 + 15*u^100 + 
#            21*u^98 + 25*u^96 + 31*u^94 + 40*u^92 + 44*u^90 + 51*u^88 + 60*u^
#            86 + 61*u^84 + 69*u^82 + 75*u^80 + 72*u^78 + 78*u^76 + 78*u^74 + 
#            72*u^72 + 75*u^70 + 69*u^68 + 61*u^66 + 60*u^64 + 51*u^62 + 44*u^
#            60 + 40*u^58 + 31*u^56 + 25*u^54 + 21*u^52 + 15*u^50 + 11*u^48 + 
#            8*u^46 + 5*u^44 + 3*u^42 + 2*u^40 + u^38, 
#          u^86 + 2*u^84 + 5*u^82 + 9*u^80 + 12*u^78 + 19*u^76 + 25*u^74 + 
#            29*u^72 + 37*u^70 + 40*u^68 + 41*u^66 + 45*u^64 + 42*u^62 + 39*u^
#            60 + 37*u^58 + 30*u^56 + 25*u^54 + 21*u^52 + 15*u^50 + 11*u^48 + 
#            8*u^46 + 5*u^44 + 3*u^42 + 2*u^40 + u^38, 
#          2*u^70 + 4*u^68 + 7*u^66 + 12*u^64 + 15*u^62 + 19*u^60 + 22*u^58 + 
#            21*u^56 + 20*u^54 + 18*u^52 + 14*u^50 + 11*u^48 + 8*u^46 + 5*u^
#            44 + 3*u^42 + 2*u^40 + u^38, u^60 + 2*u^58 + 4*u^56 + 5*u^54 + 
#            7*u^52 + 8*u^50 + 7*u^48 + 7*u^46 + 5*u^44 + 3*u^42 + 2*u^40 + u^
#            38, u^61 + u^59 + 2*u^57 + 2*u^55 + 2*u^53 + 2*u^51 + u^49 + u^47,
#          u^59 + u^57 + 3*u^55 + 5*u^53 + 4*u^51 + 6*u^49 + 5*u^47 + 3*u^45 + 
#            3*u^43 + u^41, u^50 + 2*u^48 + 4*u^46 + 4*u^44 + 3*u^42 + 2*u^
#            40 + u^38, 2*u^49 + 3*u^47 + 3*u^45 + 3*u^43 + u^41, 
#          u^48 + u^46 + u^44, u^46 + 2*u^44 + 3*u^42 + 2*u^40 + u^38, 
#          u^44 + u^42 + u^40 + u^38, u^40 + u^38, 0*u^0, u^40 + u^38, 0*u^0, 
#          u^37, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^111 + 2*u^109 + 2*u^107 + 5*u^105 + 7*u^103 + 8*u^101 + 14*u^99 + 
#            17*u^97 + 19*u^95 + 28*u^93 + 31*u^91 + 34*u^89 + 44*u^87 + 45*u^
#            85 + 48*u^83 + 57*u^81 + 55*u^79 + 56*u^77 + 62*u^75 + 56*u^73 + 
#            55*u^71 + 57*u^69 + 48*u^67 + 45*u^65 + 44*u^63 + 34*u^61 + 31*u^
#            59 + 28*u^57 + 19*u^55 + 17*u^53 + 14*u^51 + 8*u^49 + 7*u^47 + 
#            5*u^45 + 2*u^43 + 2*u^41 + u^39, u^87 + u^85 + 3*u^83 + 6*u^81 + 
#            8*u^79 + 12*u^77 + 18*u^75 + 20*u^73 + 25*u^71 + 30*u^69 + 30*u^
#            67 + 32*u^65 + 34*u^63 + 29*u^61 + 28*u^59 + 26*u^57 + 19*u^55 + 
#            17*u^53 + 14*u^51 + 8*u^49 + 7*u^47 + 5*u^45 + 2*u^43 + 2*u^
#            41 + u^39, u^71 + 2*u^69 + 4*u^67 + 7*u^65 + 10*u^63 + 12*u^61 + 
#            15*u^59 + 16*u^57 + 14*u^55 + 14*u^53 + 12*u^51 + 8*u^49 + 7*u^
#            47 + 5*u^45 + 2*u^43 + 2*u^41 + u^39, 
#          u^61 + u^59 + 4*u^57 + 4*u^55 + 5*u^53 + 7*u^51 + 5*u^49 + 5*u^47 + 
#            5*u^45 + 2*u^43 + 2*u^41 + u^39, 0*u^0, 
#          u^60 + 2*u^58 + 2*u^56 + 5*u^54 + 4*u^52 + 4*u^50 + 5*u^48 + 3*u^
#            46 + 2*u^44 + 2*u^42, u^51 + u^49 + 2*u^47 + 3*u^45 + 2*u^43 + 
#            2*u^41 + u^39, u^50 + 2*u^48 + 2*u^46 + 2*u^44 + 2*u^42, 0*u^0, 
#          u^45 + u^43 + 2*u^41 + u^39, u^45 + u^43 + 2*u^41 + u^39 + u^37, 
#          u^39, 0*u^0, u^39, 0*u^0, 0*u^0, u^37, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^111 + u^109 + u^107 + 3*u^105 + 3*u^103 + 4*u^101 + 7*u^99 + 7*u^
#            97 + 9*u^95 + 13*u^93 + 13*u^91 + 16*u^89 + 19*u^87 + 19*u^85 + 
#            22*u^83 + 24*u^81 + 24*u^79 + 25*u^77 + 26*u^75 + 25*u^73 + 24*u^
#            71 + 24*u^69 + 22*u^67 + 19*u^65 + 19*u^63 + 16*u^61 + 13*u^59 + 
#            13*u^57 + 9*u^55 + 7*u^53 + 7*u^51 + 4*u^49 + 3*u^47 + 3*u^45 + u^
#            43 + u^41 + u^39, u^83 + 2*u^81 + 3*u^79 + 5*u^77 + 7*u^75 + 9*u^
#            73 + 11*u^71 + 13*u^69 + 14*u^67 + 14*u^65 + 15*u^63 + 14*u^61 + 
#            12*u^59 + 12*u^57 + 9*u^55 + 7*u^53 + 7*u^51 + 4*u^49 + 3*u^47 + 
#            3*u^45 + u^43 + u^41 + u^39, u^69 + 2*u^67 + 3*u^65 + 5*u^63 + 
#            6*u^61 + 7*u^59 + 8*u^57 + 7*u^55 + 6*u^53 + 6*u^51 + 4*u^49 + 
#            3*u^47 + 3*u^45 + u^43 + u^41 + u^39, 
#          u^57 + u^55 + u^53 + 3*u^51 + 2*u^49 + 2*u^47 + 3*u^45 + u^43 + u^
#            41 + u^39, u^60 + u^58 + u^56 + 2*u^54 + u^52 + u^50 + u^48, 
#          u^54 + u^52 + u^50 + 2*u^48 + u^46 + u^44 + u^42, 
#          u^47 + 2*u^45 + u^43 + u^41 + u^39, u^48 + u^46 + u^44 + u^42, 
#          u^47 + u^45, u^45 + u^43 + u^41 + u^39, 0*u^0, u^39, 0*u^0, u^39, 
#          0*u^0, 0*u^0, 0*u^0, u^37, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^101 + u^99 + 2*u^95 + u^93 + u^91 + 4*u^89 + u^87 + 2*u^85 + 5*u^
#            83 + u^81 + 4*u^79 + 5*u^77 + 5*u^73 + 4*u^71 + u^69 + 5*u^67 + 
#            2*u^65 + u^63 + 4*u^61 + u^59 + u^57 + 2*u^55 + u^51 + u^49, 
#          u^83 + u^79 + u^77 + 2*u^73 + 2*u^71 + 3*u^67 + u^65 + u^63 + 3*u^
#            61 + u^59 + u^57 + 2*u^55 + u^51 + u^49, 
#          u^67 + u^61 + u^57 + u^55 + u^51 + u^49, u^49, 0*u^0, 
#          u^62 + u^56 + u^52 + u^50 + u^46, 0*u^0, u^46, 0*u^0, 0*u^0, u^41, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^37, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^108 + u^106 + u^104 + 3*u^102 + 4*u^100 + 4*u^98 + 8*u^96 + 9*u^
#            94 + 10*u^92 + 15*u^90 + 16*u^88 + 17*u^86 + 24*u^84 + 23*u^82 + 
#            24*u^80 + 30*u^78 + 28*u^76 + 28*u^74 + 33*u^72 + 28*u^70 + 28*u^
#            68 + 30*u^66 + 24*u^64 + 23*u^62 + 24*u^60 + 17*u^58 + 16*u^56 + 
#            15*u^54 + 10*u^52 + 9*u^50 + 8*u^48 + 4*u^46 + 4*u^44 + 3*u^
#            42 + u^40 + u^38 + u^36, u^84 + u^82 + 2*u^80 + 4*u^78 + 5*u^76 + 
#            7*u^74 + 11*u^72 + 11*u^70 + 14*u^68 + 17*u^66 + 16*u^64 + 17*u^
#            62 + 19*u^60 + 15*u^58 + 15*u^56 + 14*u^54 + 10*u^52 + 9*u^50 + 
#            8*u^48 + 4*u^46 + 4*u^44 + 3*u^42 + u^40 + u^38 + u^36, 
#          u^68 + 2*u^66 + 3*u^64 + 5*u^62 + 7*u^60 + 7*u^58 + 9*u^56 + 9*u^
#            54 + 8*u^52 + 8*u^50 + 7*u^48 + 4*u^46 + 4*u^44 + 3*u^42 + u^
#            40 + u^38 + u^36, u^60 + u^58 + u^56 + 3*u^54 + 3*u^52 + 3*u^50 + 
#            5*u^48 + 3*u^46 + 3*u^44 + 3*u^42 + u^40 + u^38 + u^36, 0*u^0, 
#          u^57 + u^55 + u^53 + 3*u^51 + 2*u^49 + 2*u^47 + 3*u^45 + u^43 + u^
#            41 + u^39, u^48 + u^46 + 2*u^44 + 2*u^42 + u^40 + u^38 + u^36, 
#          u^47 + 2*u^45 + u^43 + u^41 + u^39, 0*u^0, 
#          u^44 + u^42 + u^40 + u^38 + u^36, u^44 + u^42 + u^40 + u^38 + u^36, 
#          u^36, 0*u^0, u^36, u^36, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^36, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^112 + 2*u^108 + 3*u^106 + 4*u^104 + 7*u^102 + 12*u^100 + 11*u^98 + 
#            21*u^96 + 24*u^94 + 28*u^92 + 38*u^90 + 45*u^88 + 45*u^86 + 62*u^
#            84 + 61*u^82 + 66*u^80 + 77*u^78 + 76*u^76 + 74*u^74 + 86*u^72 + 
#            74*u^70 + 76*u^68 + 77*u^66 + 66*u^64 + 61*u^62 + 62*u^60 + 45*u^
#            58 + 45*u^56 + 38*u^54 + 28*u^52 + 24*u^50 + 21*u^48 + 11*u^46 + 
#            12*u^44 + 7*u^42 + 4*u^40 + 3*u^38 + 2*u^36 + u^32, 
#          2*u^84 + 2*u^82 + 5*u^80 + 10*u^78 + 13*u^76 + 18*u^74 + 28*u^72 + 
#            29*u^70 + 38*u^68 + 44*u^66 + 44*u^64 + 46*u^62 + 50*u^60 + 40*u^
#            58 + 42*u^56 + 36*u^54 + 28*u^52 + 24*u^50 + 21*u^48 + 11*u^46 + 
#            12*u^44 + 7*u^42 + 4*u^40 + 3*u^38 + 2*u^36 + u^32, 
#          3*u^68 + 4*u^66 + 8*u^64 + 13*u^62 + 19*u^60 + 19*u^58 + 27*u^56 + 
#            24*u^54 + 23*u^52 + 21*u^50 + 19*u^48 + 11*u^46 + 12*u^44 + 7*u^
#            42 + 4*u^40 + 3*u^38 + 2*u^36 + u^32, 
#          u^60 + 2*u^58 + 3*u^56 + 6*u^54 + 8*u^52 + 8*u^50 + 12*u^48 + 8*u^
#            46 + 9*u^44 + 7*u^42 + 4*u^40 + 3*u^38 + 2*u^36 + u^32, 
#          u^59 + 2*u^55 + 2*u^53 + 2*u^51 + 2*u^49 + 2*u^47 + u^43, 
#          u^57 + 3*u^55 + 2*u^53 + 7*u^51 + 5*u^49 + 5*u^47 + 6*u^45 + 3*u^
#            43 + 2*u^41 + 2*u^39, 3*u^48 + 4*u^46 + 6*u^44 + 5*u^42 + 4*u^
#            40 + 3*u^38 + 2*u^36 + u^32, 2*u^47 + 4*u^45 + 2*u^43 + 2*u^41 + 
#            2*u^39, u^48 + 2*u^46 + 2*u^44 + u^40, 
#          3*u^44 + 3*u^42 + 3*u^40 + 3*u^38 + 2*u^36 + u^32, 
#          u^44 + u^42 + 2*u^40 + 2*u^38 + 2*u^36 + u^34 + u^32, 
#          u^40 + 2*u^36 + u^32, u^41 + u^37, u^36, u^36 + u^32, 0*u^0, u^34, 
#          0*u^0, 0*u^0, 0*u^0, u^32, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^110 + 2*u^108 + 2*u^106 + 6*u^104 + 7*u^102 + 10*u^100 + 17*u^98 + 
#            20*u^96 + 25*u^94 + 37*u^92 + 38*u^90 + 48*u^88 + 60*u^86 + 61*u^
#            84 + 71*u^82 + 82*u^80 + 77*u^78 + 89*u^76 + 92*u^74 + 85*u^72 + 
#            92*u^70 + 89*u^68 + 77*u^66 + 82*u^64 + 71*u^62 + 61*u^60 + 60*u^
#            58 + 48*u^56 + 38*u^54 + 37*u^52 + 25*u^50 + 20*u^48 + 17*u^46 + 
#            10*u^44 + 7*u^42 + 6*u^40 + 2*u^38 + 2*u^36 + u^34, 
#          u^86 + u^84 + 4*u^82 + 7*u^80 + 9*u^78 + 17*u^76 + 23*u^74 + 27*u^
#            72 + 38*u^70 + 43*u^68 + 45*u^66 + 54*u^64 + 52*u^62 + 50*u^60 + 
#            52*u^58 + 44*u^56 + 37*u^54 + 36*u^52 + 25*u^50 + 20*u^48 + 17*u^
#            46 + 10*u^44 + 7*u^42 + 6*u^40 + 2*u^38 + 2*u^36 + u^34, 
#          u^70 + 2*u^68 + 5*u^66 + 11*u^64 + 13*u^62 + 20*u^60 + 25*u^58 + 
#            25*u^56 + 26*u^54 + 28*u^52 + 21*u^50 + 19*u^48 + 16*u^46 + 10*u^
#            44 + 7*u^42 + 6*u^40 + 2*u^38 + 2*u^36 + u^34, 
#          u^60 + 2*u^58 + 5*u^56 + 6*u^54 + 9*u^52 + 11*u^50 + 10*u^48 + 11*u^
#            46 + 9*u^44 + 6*u^42 + 6*u^40 + 2*u^38 + 2*u^36 + u^34, 
#          u^57 + u^55 + 2*u^51 + u^47 + u^45, 
#          2*u^59 + 2*u^57 + 4*u^55 + 7*u^53 + 6*u^51 + 8*u^49 + 9*u^47 + 5*u^
#            45 + 6*u^43 + 3*u^41 + u^39 + u^37, 
#          u^50 + 2*u^48 + 4*u^46 + 5*u^44 + 5*u^42 + 5*u^40 + 2*u^38 + 2*u^
#            36 + u^34, 2*u^49 + 3*u^47 + 3*u^45 + 5*u^43 + 3*u^41 + u^39 + u^
#            37, u^48 + u^44 + u^42, u^46 + u^44 + 3*u^42 + 4*u^40 + 2*u^38 + 
#            2*u^36 + u^34, u^44 + 2*u^42 + 3*u^40 + 3*u^38 + 3*u^36 + 2*u^
#            34 + u^32, u^38 + u^36 + u^34, u^39, u^38 + u^34, u^34, 0*u^0, 
#          u^36 + u^34 + u^32, 0*u^0, u^34, 0*u^0, 0*u^0, u^32, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^110 + u^108 + 3*u^106 + 3*u^104 + 6*u^102 + 8*u^100 + 11*u^98 + 
#            14*u^96 + 20*u^94 + 21*u^92 + 29*u^90 + 32*u^88 + 37*u^86 + 43*u^
#            84 + 48*u^82 + 48*u^80 + 57*u^78 + 55*u^76 + 58*u^74 + 60*u^72 + 
#            58*u^70 + 55*u^68 + 57*u^66 + 48*u^64 + 48*u^62 + 43*u^60 + 37*u^
#            58 + 32*u^56 + 29*u^54 + 21*u^52 + 20*u^50 + 14*u^48 + 11*u^46 + 
#            8*u^44 + 6*u^42 + 3*u^40 + 3*u^38 + u^36 + u^34, 
#          u^84 + 2*u^82 + 3*u^80 + 7*u^78 + 9*u^76 + 14*u^74 + 19*u^72 + 23*u^
#            70 + 27*u^68 + 33*u^66 + 32*u^64 + 36*u^62 + 35*u^60 + 33*u^58 + 
#            30*u^56 + 28*u^54 + 21*u^52 + 20*u^50 + 14*u^48 + 11*u^46 + 8*u^
#            44 + 6*u^42 + 3*u^40 + 3*u^38 + u^36 + u^34, 
#          u^68 + 4*u^66 + 5*u^64 + 11*u^62 + 13*u^60 + 17*u^58 + 18*u^56 + 
#            20*u^54 + 17*u^52 + 18*u^50 + 13*u^48 + 11*u^46 + 8*u^44 + 6*u^
#            42 + 3*u^40 + 3*u^38 + u^36 + u^34, 
#          u^58 + 2*u^56 + 3*u^54 + 5*u^52 + 6*u^50 + 7*u^48 + 7*u^46 + 6*u^
#            44 + 6*u^42 + 3*u^40 + 3*u^38 + u^36 + u^34, 
#          u^59 + 2*u^57 + 2*u^55 + 4*u^53 + 2*u^51 + 4*u^49 + 2*u^47 + 2*u^
#            45 + u^43, u^57 + u^55 + 2*u^53 + 3*u^51 + 4*u^49 + 3*u^47 + 5*u^
#            45 + 2*u^43 + 2*u^41 + u^39, u^48 + 3*u^46 + 4*u^44 + 5*u^42 + 
#            3*u^40 + 3*u^38 + u^36 + u^34, 
#          u^47 + 3*u^45 + 2*u^43 + 2*u^41 + u^39, 
#          u^50 + u^48 + 3*u^46 + 2*u^44 + 2*u^42 + u^40, 
#          2*u^44 + 3*u^42 + 3*u^40 + 3*u^38 + u^36 + u^34, 
#          u^42 + u^40 + u^38 + u^36 + u^34, u^38 + u^36 + u^34, 
#          u^41 + u^39 + u^37, u^36, 0*u^0, u^35 + u^33, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^32, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^108 + u^104 + u^102 + u^100 + 2*u^98 + 3*u^96 + 2*u^94 + 5*u^92 + 
#            4*u^90 + 5*u^88 + 7*u^86 + 7*u^84 + 7*u^82 + 10*u^80 + 7*u^78 + 
#            10*u^76 + 10*u^74 + 9*u^72 + 10*u^70 + 10*u^68 + 7*u^66 + 10*u^
#            64 + 7*u^62 + 7*u^60 + 7*u^58 + 5*u^56 + 4*u^54 + 5*u^52 + 2*u^
#            50 + 3*u^48 + 2*u^46 + u^44 + u^42 + u^40 + u^36, 
#          u^80 + 2*u^76 + 2*u^74 + 3*u^72 + 4*u^70 + 5*u^68 + 4*u^66 + 7*u^
#            64 + 5*u^62 + 6*u^60 + 6*u^58 + 5*u^56 + 4*u^54 + 5*u^52 + 2*u^
#            50 + 3*u^48 + 2*u^46 + u^44 + u^42 + u^40 + u^36, 
#          2*u^64 + u^62 + 3*u^60 + 3*u^58 + 3*u^56 + 3*u^54 + 4*u^52 + 2*u^
#            50 + 3*u^48 + 2*u^46 + u^44 + u^42 + u^40 + u^36, 
#          u^54 + u^50 + u^48 + u^46 + u^44 + u^42 + u^40 + u^36, 
#          u^57 + u^55 + 2*u^51 + u^47 + u^45, u^47 + u^43, 
#          u^44 + u^42 + u^40 + u^36, u^43, u^48 + u^44 + u^42, 
#          u^42 + u^40 + u^36, 0*u^0, u^36, u^39, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^34, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^32, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^106 + 2*u^102 + 3*u^100 + 2*u^98 + 5*u^96 + 7*u^94 + 5*u^92 + 11*u^
#            90 + 11*u^88 + 10*u^86 + 17*u^84 + 16*u^82 + 14*u^80 + 23*u^78 + 
#            18*u^76 + 18*u^74 + 24*u^72 + 18*u^70 + 18*u^68 + 23*u^66 + 14*u^
#            64 + 16*u^62 + 17*u^60 + 10*u^58 + 11*u^56 + 11*u^54 + 5*u^52 + 
#            7*u^50 + 5*u^48 + 2*u^46 + 3*u^44 + 2*u^42 + u^38, 
#          u^84 + u^82 + u^80 + 4*u^78 + 3*u^76 + 5*u^74 + 8*u^72 + 7*u^70 + 
#            9*u^68 + 13*u^66 + 9*u^64 + 12*u^62 + 13*u^60 + 9*u^58 + 10*u^
#            56 + 10*u^54 + 5*u^52 + 7*u^50 + 5*u^48 + 2*u^46 + 3*u^44 + 2*u^
#            42 + u^38, u^68 + 2*u^66 + u^64 + 4*u^62 + 4*u^60 + 4*u^58 + 6*u^
#            56 + 6*u^54 + 4*u^52 + 6*u^50 + 4*u^48 + 2*u^46 + 3*u^44 + 2*u^
#            42 + u^38, u^58 + 2*u^54 + u^52 + 2*u^50 + 3*u^48 + u^46 + 2*u^
#            44 + 2*u^42 + u^38, 0*u^0, u^61 + 2*u^57 + 2*u^55 + u^53 + 4*u^
#            51 + 2*u^49 + 2*u^47 + 3*u^45 + u^43 + u^41 + u^39, 
#          u^44 + u^42 + u^38, u^51 + u^47 + 2*u^45 + u^43 + u^41 + u^39, 
#          0*u^0, u^42 + u^38, u^42 + u^40 + u^38 + u^36 + u^34, 0*u^0, 0*u^0, 
#          u^36, 0*u^0, 0*u^0, u^34, 0*u^0, u^36 + u^32, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^32, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^111 + 2*u^109 + 4*u^107 + 7*u^105 + 11*u^103 + 17*u^101 + 24*u^99 + 
#            33*u^97 + 44*u^95 + 56*u^93 + 70*u^91 + 85*u^89 + 100*u^87 + 
#            116*u^85 + 131*u^83 + 145*u^81 + 158*u^79 + 168*u^77 + 176*u^75 + 
#            181*u^73 + 182*u^71 + 181*u^69 + 176*u^67 + 168*u^65 + 158*u^63 + 
#            145*u^61 + 131*u^59 + 116*u^57 + 100*u^55 + 85*u^53 + 70*u^51 + 
#            56*u^49 + 44*u^47 + 33*u^45 + 24*u^43 + 17*u^41 + 11*u^39 + 7*u^
#            37 + 4*u^35 + 2*u^33 + u^31, u^85 + 3*u^83 + 7*u^81 + 13*u^79 + 
#            22*u^77 + 33*u^75 + 47*u^73 + 61*u^71 + 76*u^69 + 89*u^67 + 100*u^
#            65 + 107*u^63 + 110*u^61 + 108*u^59 + 103*u^57 + 93*u^55 + 82*u^
#            53 + 69*u^51 + 56*u^49 + 44*u^47 + 33*u^45 + 24*u^43 + 17*u^41 + 
#            11*u^39 + 7*u^37 + 4*u^35 + 2*u^33 + u^31, 
#          2*u^69 + 6*u^67 + 13*u^65 + 22*u^63 + 33*u^61 + 44*u^59 + 53*u^57 + 
#            58*u^55 + 59*u^53 + 56*u^51 + 49*u^49 + 41*u^47 + 32*u^45 + 24*u^
#            43 + 17*u^41 + 11*u^39 + 7*u^37 + 4*u^35 + 2*u^33 + u^31, 
#          2*u^59 + 4*u^57 + 9*u^55 + 13*u^53 + 19*u^51 + 22*u^49 + 24*u^47 + 
#            23*u^45 + 20*u^43 + 16*u^41 + 11*u^39 + 7*u^37 + 4*u^35 + 2*u^
#            33 + u^31, u^60 + 2*u^58 + 4*u^56 + 5*u^54 + 6*u^52 + 6*u^50 + 
#            5*u^48 + 4*u^46 + 2*u^44 + u^42, 
#          u^58 + 3*u^56 + 6*u^54 + 10*u^52 + 12*u^50 + 15*u^48 + 14*u^46 + 
#            13*u^44 + 9*u^42 + 6*u^40 + 3*u^38 + u^36, 
#          2*u^49 + 6*u^47 + 11*u^45 + 13*u^43 + 13*u^41 + 10*u^39 + 7*u^37 + 
#            4*u^35 + 2*u^33 + u^31, 3*u^48 + 6*u^46 + 9*u^44 + 8*u^42 + 6*u^
#            40 + 3*u^38 + u^36, u^49 + 3*u^47 + 4*u^45 + 4*u^43 + 2*u^41 + u^
#            39, 2*u^45 + 6*u^43 + 9*u^41 + 9*u^39 + 7*u^37 + 4*u^35 + 2*u^
#            33 + u^31, 2*u^43 + 3*u^41 + 5*u^39 + 5*u^37 + 5*u^35 + 3*u^33 + 
#            2*u^31, u^41 + 2*u^39 + 3*u^37 + 3*u^35 + 2*u^33 + u^31, 
#          u^40 + u^38 + u^36, u^39 + 2*u^37 + 2*u^35 + u^33, 
#          u^35 + u^33 + u^31, u^36 + u^34 + u^32, u^35 + u^33 + u^31, 0*u^0, 
#          0*u^0, 0*u^0, u^31, u^31, u^31, 0*u^0, 0*u^0, u^31, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^110 + u^108 + 3*u^106 + 4*u^104 + 8*u^102 + 11*u^100 + 17*u^98 + 
#            22*u^96 + 31*u^94 + 38*u^92 + 50*u^90 + 58*u^88 + 71*u^86 + 80*u^
#            84 + 93*u^82 + 101*u^80 + 112*u^78 + 117*u^76 + 125*u^74 + 126*u^
#            72 + 130*u^70 + 126*u^68 + 125*u^66 + 117*u^64 + 112*u^62 + 101*u^
#            60 + 93*u^58 + 80*u^56 + 71*u^54 + 58*u^52 + 50*u^50 + 38*u^48 + 
#            31*u^46 + 22*u^44 + 17*u^42 + 11*u^40 + 8*u^38 + 4*u^36 + 3*u^
#            34 + u^32 + u^30, u^84 + 3*u^82 + 6*u^80 + 11*u^78 + 17*u^76 + 
#            26*u^74 + 35*u^72 + 46*u^70 + 55*u^68 + 65*u^66 + 71*u^64 + 77*u^
#            62 + 77*u^60 + 77*u^58 + 71*u^56 + 66*u^54 + 56*u^52 + 49*u^50 + 
#            38*u^48 + 31*u^46 + 22*u^44 + 17*u^42 + 11*u^40 + 8*u^38 + 4*u^
#            36 + 3*u^34 + u^32 + u^30, u^70 + 2*u^68 + 6*u^66 + 10*u^64 + 
#            18*u^62 + 24*u^60 + 33*u^58 + 37*u^56 + 42*u^54 + 40*u^52 + 40*u^
#            50 + 33*u^48 + 29*u^46 + 21*u^44 + 17*u^42 + 11*u^40 + 8*u^38 + 
#            4*u^36 + 3*u^34 + u^32 + u^30, u^60 + 2*u^58 + 5*u^56 + 8*u^54 + 
#            12*u^52 + 16*u^50 + 17*u^48 + 19*u^46 + 16*u^44 + 15*u^42 + 10*u^
#            40 + 8*u^38 + 4*u^36 + 3*u^34 + u^32 + u^30, 
#          u^57 + u^55 + 2*u^53 + u^51 + 2*u^49 + u^47 + u^45, 
#          2*u^57 + 3*u^55 + 7*u^53 + 8*u^51 + 12*u^49 + 11*u^47 + 12*u^45 + 
#            9*u^43 + 7*u^41 + 4*u^39 + 2*u^37 + u^35, 
#          u^50 + 2*u^48 + 6*u^46 + 8*u^44 + 10*u^42 + 8*u^40 + 7*u^38 + 4*u^
#            36 + 3*u^34 + u^32 + u^30, u^49 + 3*u^47 + 6*u^45 + 7*u^43 + 6*u^
#            41 + 4*u^39 + 2*u^37 + u^35, u^46 + u^44 + u^42, 
#          2*u^44 + 5*u^42 + 6*u^40 + 6*u^38 + 4*u^36 + 3*u^34 + u^32 + u^30, 
#          u^44 + 3*u^42 + 4*u^40 + 5*u^38 + 5*u^36 + 4*u^34 + 2*u^32 + 2*u^30,
#          u^40 + 2*u^38 + 2*u^36 + 2*u^34 + u^32 + u^30, 0*u^0, 
#          u^40 + u^38 + 2*u^36 + u^34 + u^32, u^36 + u^34 + u^32 + u^30, 
#          0*u^0, u^36 + u^34 + u^32 + u^30, 0*u^0, 0*u^0, 0*u^0, u^30, u^30, 
#          0*u^0, 0*u^0, 0*u^0, u^30, u^30, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^108 + u^106 + 3*u^104 + 3*u^102 + 6*u^100 + 7*u^98 + 11*u^96 + 13*u^
#            94 + 18*u^92 + 20*u^90 + 27*u^88 + 29*u^86 + 36*u^84 + 38*u^82 + 
#            44*u^80 + 46*u^78 + 51*u^76 + 51*u^74 + 55*u^72 + 52*u^70 + 55*u^
#            68 + 51*u^66 + 51*u^64 + 46*u^62 + 44*u^60 + 38*u^58 + 36*u^56 + 
#            29*u^54 + 27*u^52 + 20*u^50 + 18*u^48 + 13*u^46 + 11*u^44 + 7*u^
#            42 + 6*u^40 + 3*u^38 + 3*u^36 + u^34 + u^32, 
#          u^82 + 2*u^80 + 4*u^78 + 7*u^76 + 10*u^74 + 15*u^72 + 18*u^70 + 
#            24*u^68 + 27*u^66 + 31*u^64 + 32*u^62 + 34*u^60 + 32*u^58 + 32*u^
#            56 + 27*u^54 + 26*u^52 + 20*u^50 + 18*u^48 + 13*u^46 + 11*u^44 + 
#            7*u^42 + 6*u^40 + 3*u^38 + 3*u^36 + u^34 + u^32, 
#          u^68 + 2*u^66 + 5*u^64 + 7*u^62 + 12*u^60 + 13*u^58 + 18*u^56 + 
#            17*u^54 + 20*u^52 + 16*u^50 + 16*u^48 + 12*u^46 + 11*u^44 + 7*u^
#            42 + 6*u^40 + 3*u^38 + 3*u^36 + u^34 + u^32, 
#          u^56 + 2*u^54 + 4*u^52 + 5*u^50 + 7*u^48 + 7*u^46 + 8*u^44 + 6*u^
#            42 + 6*u^40 + 3*u^38 + 3*u^36 + u^34 + u^32, 
#          u^59 + u^57 + 3*u^55 + 2*u^53 + 4*u^51 + 2*u^49 + 3*u^47 + u^45 + u^
#            43, u^55 + u^53 + 3*u^51 + 3*u^49 + 5*u^47 + 4*u^45 + 4*u^43 + 
#            3*u^41 + 2*u^39 + u^37, u^48 + 2*u^46 + 4*u^44 + 4*u^42 + 5*u^
#            40 + 3*u^38 + 3*u^36 + u^34 + u^32, 
#          u^47 + 2*u^45 + 3*u^43 + 3*u^41 + 2*u^39 + u^37, 
#          u^48 + u^46 + 2*u^44 + u^42 + u^40, 
#          u^44 + 2*u^42 + 4*u^40 + 3*u^38 + 3*u^36 + u^34 + u^32, 
#          u^40 + u^38 + 2*u^36 + u^34 + u^32, u^40 + u^38 + 2*u^36 + u^34 + u^
#            32, 0*u^0, u^38 + u^36 + u^34, 0*u^0, u^35 + u^33 + u^31, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^30, 0*u^0, 0*u^0, u^30, 0*u^0, 
#          u^30, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^107 + 2*u^105 + 3*u^103 + 5*u^101 + 7*u^99 + 11*u^97 + 15*u^95 + 
#            19*u^93 + 25*u^91 + 30*u^89 + 37*u^87 + 44*u^85 + 49*u^83 + 56*u^
#            81 + 61*u^79 + 66*u^77 + 71*u^75 + 72*u^73 + 74*u^71 + 74*u^69 + 
#            72*u^67 + 71*u^65 + 66*u^63 + 61*u^61 + 56*u^59 + 49*u^57 + 44*u^
#            55 + 37*u^53 + 30*u^51 + 25*u^49 + 19*u^47 + 15*u^45 + 11*u^43 + 
#            7*u^41 + 5*u^39 + 3*u^37 + 2*u^35 + u^33, 
#          u^85 + u^83 + 3*u^81 + 5*u^79 + 9*u^77 + 13*u^75 + 18*u^73 + 23*u^
#            71 + 30*u^69 + 34*u^67 + 40*u^65 + 42*u^63 + 44*u^61 + 44*u^59 + 
#            42*u^57 + 39*u^55 + 35*u^53 + 29*u^51 + 25*u^49 + 19*u^47 + 15*u^
#            45 + 11*u^43 + 7*u^41 + 5*u^39 + 3*u^37 + 2*u^35 + u^33, 
#          u^69 + 2*u^67 + 5*u^65 + 8*u^63 + 12*u^61 + 16*u^59 + 19*u^57 + 
#            22*u^55 + 23*u^53 + 22*u^51 + 20*u^49 + 17*u^47 + 14*u^45 + 11*u^
#            43 + 7*u^41 + 5*u^39 + 3*u^37 + 2*u^35 + u^33, 
#          u^59 + 2*u^57 + 4*u^55 + 6*u^53 + 8*u^51 + 9*u^49 + 10*u^47 + 9*u^
#            45 + 9*u^43 + 6*u^41 + 5*u^39 + 3*u^37 + 2*u^35 + u^33, 0*u^0, 
#          u^60 + 2*u^58 + 3*u^56 + 5*u^54 + 6*u^52 + 8*u^50 + 8*u^48 + 8*u^
#            46 + 7*u^44 + 5*u^42 + 4*u^40 + 2*u^38 + u^36, 
#          u^49 + 2*u^47 + 3*u^45 + 4*u^43 + 4*u^41 + 4*u^39 + 3*u^37 + 2*u^
#            35 + u^33, u^50 + 2*u^48 + 3*u^46 + 4*u^44 + 4*u^42 + 4*u^40 + 
#            2*u^38 + u^36, 0*u^0, u^45 + u^43 + 2*u^41 + 3*u^39 + 3*u^37 + 
#            2*u^35 + u^33, u^45 + u^43 + 2*u^41 + 4*u^39 + 4*u^37 + 4*u^35 + 
#            3*u^33 + u^31, u^37 + u^35 + u^33, 0*u^0, u^37 + u^35 + u^33, 
#          u^33, 0*u^0, u^35 + 2*u^33 + u^31, 0*u^0, u^35 + u^33 + u^31, 
#          0*u^0, 0*u^0, u^31, 0*u^0, 0*u^0, u^31, 0*u^0, 0*u^0, 0*u^0, u^30, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^110 + u^108 + 2*u^106 + 5*u^104 + 6*u^102 + 10*u^100 + 16*u^98 + 
#            19*u^96 + 27*u^94 + 37*u^92 + 42*u^90 + 54*u^88 + 66*u^86 + 72*u^
#            84 + 86*u^82 + 97*u^80 + 101*u^78 + 113*u^76 + 120*u^74 + 119*u^
#            72 + 126*u^70 + 126*u^68 + 119*u^66 + 120*u^64 + 113*u^62 + 101*u^
#            60 + 97*u^58 + 86*u^56 + 72*u^54 + 66*u^52 + 54*u^50 + 42*u^48 + 
#            37*u^46 + 27*u^44 + 19*u^42 + 16*u^40 + 10*u^38 + 6*u^36 + 5*u^
#            34 + 2*u^32 + u^30 + u^28, 2*u^82 + 4*u^80 + 7*u^78 + 14*u^76 + 
#            21*u^74 + 29*u^72 + 41*u^70 + 51*u^68 + 59*u^66 + 70*u^64 + 75*u^
#            62 + 76*u^60 + 79*u^58 + 75*u^56 + 67*u^54 + 63*u^52 + 53*u^50 + 
#            42*u^48 + 37*u^46 + 27*u^44 + 19*u^42 + 16*u^40 + 10*u^38 + 6*u^
#            36 + 5*u^34 + 2*u^32 + u^30 + u^28, 
#          u^68 + 4*u^66 + 9*u^64 + 15*u^62 + 23*u^60 + 32*u^58 + 38*u^56 + 
#            42*u^54 + 45*u^52 + 42*u^50 + 37*u^48 + 34*u^46 + 26*u^44 + 19*u^
#            42 + 16*u^40 + 10*u^38 + 6*u^36 + 5*u^34 + 2*u^32 + u^30 + u^28, 
#          u^58 + 3*u^56 + 5*u^54 + 9*u^52 + 14*u^50 + 15*u^48 + 19*u^46 + 
#            19*u^44 + 15*u^42 + 15*u^40 + 10*u^38 + 6*u^36 + 5*u^34 + 2*u^
#            32 + u^30 + u^28, u^59 + 2*u^57 + 4*u^55 + 5*u^53 + 6*u^51 + 6*u^
#            49 + 5*u^47 + 4*u^45 + 2*u^43 + u^41, 
#          u^55 + 4*u^53 + 5*u^51 + 8*u^49 + 11*u^47 + 9*u^45 + 10*u^43 + 8*u^
#            41 + 4*u^39 + 3*u^37 + u^35, u^48 + 5*u^46 + 9*u^44 + 10*u^42 + 
#            12*u^40 + 9*u^38 + 6*u^36 + 5*u^34 + 2*u^32 + u^30 + u^28, 
#          2*u^47 + 4*u^45 + 7*u^43 + 7*u^41 + 4*u^39 + 3*u^37 + u^35, 
#          u^48 + 3*u^46 + 4*u^44 + 4*u^42 + 2*u^40 + u^38, 
#          2*u^44 + 5*u^42 + 9*u^40 + 8*u^38 + 6*u^36 + 5*u^34 + 2*u^32 + u^
#            30 + u^28, u^42 + 2*u^40 + 3*u^38 + 4*u^36 + 4*u^34 + 3*u^32 + 
#            2*u^30 + u^28, u^40 + 3*u^38 + 2*u^36 + 4*u^34 + 2*u^32 + u^
#            30 + u^28, u^39 + u^37 + u^35, 2*u^38 + u^36 + 2*u^34 + u^32, 
#          u^34 + u^32 + u^30 + u^28, u^35 + u^33 + u^31, u^32 + u^30, 
#          u^36 + u^32, 0*u^0, 0*u^0, u^30 + u^28, u^30, u^30, u^30, 0*u^0, 
#          u^30 + u^28, 0*u^0, 0*u^0, 0*u^0, u^28, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^109 + u^107 + 2*u^105 + 4*u^103 + 5*u^101 + 8*u^99 + 12*u^97 + 14*u^
#            95 + 20*u^93 + 26*u^91 + 29*u^89 + 38*u^87 + 44*u^85 + 48*u^83 + 
#            58*u^81 + 62*u^79 + 65*u^77 + 74*u^75 + 74*u^73 + 75*u^71 + 80*u^
#            69 + 75*u^67 + 74*u^65 + 74*u^63 + 65*u^61 + 62*u^59 + 58*u^57 + 
#            48*u^55 + 44*u^53 + 38*u^51 + 29*u^49 + 26*u^47 + 20*u^45 + 14*u^
#            43 + 12*u^41 + 8*u^39 + 5*u^37 + 4*u^35 + 2*u^33 + u^31 + u^29, 
#          2*u^81 + 3*u^79 + 6*u^77 + 11*u^75 + 15*u^73 + 21*u^71 + 29*u^69 + 
#            33*u^67 + 40*u^65 + 46*u^63 + 46*u^61 + 49*u^59 + 49*u^57 + 44*u^
#            55 + 42*u^53 + 37*u^51 + 29*u^49 + 26*u^47 + 20*u^45 + 14*u^43 + 
#            12*u^41 + 8*u^39 + 5*u^37 + 4*u^35 + 2*u^33 + u^31 + u^29, 
#          u^67 + 4*u^65 + 7*u^63 + 11*u^61 + 18*u^59 + 22*u^57 + 25*u^55 + 
#            29*u^53 + 28*u^51 + 25*u^49 + 24*u^47 + 19*u^45 + 14*u^43 + 12*u^
#            41 + 8*u^39 + 5*u^37 + 4*u^35 + 2*u^33 + u^31 + u^29, 
#          u^57 + 2*u^55 + 4*u^53 + 6*u^51 + 9*u^49 + 10*u^47 + 12*u^45 + 11*u^
#            43 + 10*u^41 + 8*u^39 + 5*u^37 + 4*u^35 + 2*u^33 + u^31 + u^29, 
#          u^58 + 3*u^56 + 3*u^54 + 5*u^52 + 6*u^50 + 4*u^48 + 5*u^46 + 3*u^
#            44 + u^42 + u^40, u^54 + 2*u^52 + 3*u^50 + 5*u^48 + 6*u^46 + 5*u^
#            44 + 6*u^42 + 3*u^40 + 2*u^38 + u^36, 
#          u^47 + 4*u^45 + 7*u^43 + 8*u^41 + 7*u^39 + 5*u^37 + 4*u^35 + 2*u^
#            33 + u^31 + u^29, u^46 + 3*u^44 + 5*u^42 + 3*u^40 + 2*u^38 + u^36,
#          u^49 + 2*u^47 + 3*u^45 + 5*u^43 + 3*u^41 + u^39 + u^37, 
#          2*u^43 + 6*u^41 + 6*u^39 + 5*u^37 + 4*u^35 + 2*u^33 + u^31 + u^29, 
#          u^41 + u^39 + 2*u^37 + 2*u^35 + 2*u^33 + u^31 + u^29, 
#          u^39 + 2*u^37 + 2*u^35 + 2*u^33 + u^31 + u^29, 
#          u^40 + 2*u^38 + u^36 + u^34, u^37 + u^35 + u^33, u^33 + u^29, 
#          u^34 + u^32 + u^30, 0*u^0, u^35 + u^33, 0*u^0, 0*u^0, u^29, 0*u^0, 
#          u^31 + u^29, u^31, 0*u^0, u^29, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^28, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^107 + u^105 + u^103 + 4*u^101 + 4*u^99 + 5*u^97 + 11*u^95 + 10*u^
#            93 + 13*u^91 + 22*u^89 + 19*u^87 + 25*u^85 + 35*u^83 + 29*u^81 + 
#            38*u^79 + 46*u^77 + 37*u^75 + 48*u^73 + 51*u^71 + 40*u^69 + 51*u^
#            67 + 48*u^65 + 37*u^63 + 46*u^61 + 38*u^59 + 29*u^57 + 35*u^55 + 
#            25*u^53 + 19*u^51 + 22*u^49 + 13*u^47 + 10*u^45 + 11*u^43 + 5*u^
#            41 + 4*u^39 + 4*u^37 + u^35 + u^33 + u^31, 
#          u^83 + 3*u^79 + 5*u^77 + 5*u^75 + 11*u^73 + 15*u^71 + 14*u^69 + 
#            24*u^67 + 25*u^65 + 24*u^63 + 32*u^61 + 29*u^59 + 25*u^57 + 31*u^
#            55 + 23*u^53 + 19*u^51 + 21*u^49 + 13*u^47 + 10*u^45 + 11*u^43 + 
#            5*u^41 + 4*u^39 + 4*u^37 + u^35 + u^33 + u^31, 
#          2*u^67 + 2*u^65 + 4*u^63 + 9*u^61 + 9*u^59 + 12*u^57 + 17*u^55 + 
#            14*u^53 + 15*u^51 + 17*u^49 + 11*u^47 + 10*u^45 + 10*u^43 + 5*u^
#            41 + 4*u^39 + 4*u^37 + u^35 + u^33 + u^31, 
#          u^57 + u^55 + 4*u^53 + 3*u^51 + 6*u^49 + 7*u^47 + 5*u^45 + 8*u^43 + 
#            5*u^41 + 3*u^39 + 4*u^37 + u^35 + u^33 + u^31, 
#          u^58 + 2*u^54 + u^52 + u^50 + 2*u^48 + u^44, 
#          2*u^56 + u^54 + 2*u^52 + 6*u^50 + 2*u^48 + 6*u^46 + 6*u^44 + 2*u^
#            42 + 4*u^40 + 2*u^38 + u^34, u^47 + u^45 + 4*u^43 + 3*u^41 + 3*u^
#            39 + 3*u^37 + u^35 + u^33 + u^31, 
#          u^50 + 3*u^46 + 3*u^44 + 2*u^42 + 4*u^40 + 2*u^38 + u^34, 
#          u^45 + u^41, 2*u^43 + u^41 + 3*u^39 + 3*u^37 + u^35 + u^33 + u^31, 
#          u^41 + u^39 + 2*u^37 + 2*u^35 + 2*u^33 + u^31 + u^29, 
#          u^37 + u^35 + u^31, 0*u^0, u^41 + u^37 + 2*u^35 + u^31, u^31, u^34, 
#          u^33 + u^29, 0*u^0, u^31, 0*u^0, 0*u^0, u^29, 0*u^0, 0*u^0, u^31, 
#          u^29, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^28, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^108 + 3*u^106 + 4*u^104 + 8*u^102 + 14*u^100 + 18*u^98 + 29*u^96 + 
#            41*u^94 + 49*u^92 + 69*u^90 + 86*u^88 + 98*u^86 + 126*u^84 + 
#            142*u^82 + 155*u^80 + 185*u^78 + 194*u^76 + 204*u^74 + 227*u^72 + 
#            223*u^70 + 226*u^68 + 236*u^66 + 219*u^64 + 213*u^62 + 208*u^60 + 
#            182*u^58 + 170*u^56 + 156*u^54 + 127*u^52 + 113*u^50 + 96*u^48 + 
#            72*u^46 + 61*u^44 + 47*u^42 + 31*u^40 + 25*u^38 + 17*u^36 + 9*u^
#            34 + 7*u^32 + 3*u^30 + u^28 + u^26, 
#          u^84 + 2*u^82 + 5*u^80 + 13*u^78 + 20*u^76 + 32*u^74 + 51*u^72 + 
#            64*u^70 + 84*u^68 + 107*u^66 + 117*u^64 + 133*u^62 + 145*u^60 + 
#            141*u^58 + 143*u^56 + 138*u^54 + 119*u^52 + 109*u^50 + 94*u^48 + 
#            72*u^46 + 61*u^44 + 47*u^42 + 31*u^40 + 25*u^38 + 17*u^36 + 9*u^
#            34 + 7*u^32 + 3*u^30 + u^28 + u^26, 
#          2*u^68 + 6*u^66 + 12*u^64 + 24*u^62 + 37*u^60 + 49*u^58 + 66*u^56 + 
#            76*u^54 + 78*u^52 + 82*u^50 + 76*u^48 + 64*u^46 + 57*u^44 + 45*u^
#            42 + 31*u^40 + 25*u^38 + 17*u^36 + 9*u^34 + 7*u^32 + 3*u^30 + u^
#            28 + u^26, 2*u^58 + 5*u^56 + 11*u^54 + 18*u^52 + 24*u^50 + 33*u^
#            48 + 35*u^46 + 36*u^44 + 36*u^42 + 27*u^40 + 23*u^38 + 17*u^36 + 
#            9*u^34 + 7*u^32 + 3*u^30 + u^28 + u^26, 
#          u^59 + u^57 + 3*u^55 + 5*u^53 + 4*u^51 + 6*u^49 + 5*u^47 + 3*u^45 + 
#            3*u^43 + u^41, 2*u^57 + 5*u^55 + 7*u^53 + 15*u^51 + 19*u^49 + 
#            20*u^47 + 26*u^45 + 21*u^43 + 17*u^41 + 15*u^39 + 7*u^37 + 4*u^
#            35 + 2*u^33, 3*u^48 + 7*u^46 + 13*u^44 + 19*u^42 + 19*u^40 + 19*u^
#            38 + 15*u^36 + 9*u^34 + 7*u^32 + 3*u^30 + u^28 + u^26, 
#          u^49 + 4*u^47 + 10*u^45 + 12*u^43 + 14*u^41 + 14*u^39 + 7*u^37 + 
#            4*u^35 + 2*u^33, 2*u^46 + 3*u^44 + 3*u^42 + 3*u^40 + u^38, 
#          3*u^44 + 7*u^42 + 11*u^40 + 16*u^38 + 14*u^36 + 9*u^34 + 7*u^32 + 
#            3*u^30 + u^28 + u^26, u^44 + 3*u^42 + 6*u^40 + 10*u^38 + 12*u^
#            36 + 12*u^34 + 10*u^32 + 6*u^30 + 3*u^28 + u^26, 
#          u^40 + 2*u^38 + 6*u^36 + 5*u^34 + 5*u^32 + 3*u^30 + u^28 + u^26, 
#          u^37 + u^35, u^40 + u^38 + 5*u^36 + 4*u^34 + 3*u^32 + 2*u^30, 
#          u^36 + u^34 + 3*u^32 + 3*u^30 + u^28 + u^26, u^35 + 2*u^33 + u^31, 
#          2*u^34 + 3*u^32 + 3*u^30 + 2*u^28, 0*u^0, u^32 + u^30, 0*u^0, 
#          u^30 + u^28 + u^26, 2*u^30 + 2*u^28, u^30, 0*u^0, u^30, 
#          u^30 + 2*u^28 + u^26, u^28 + u^26, 0*u^0, u^29 + u^27, u^26, 0*u^0, 
#          u^27, u^26, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^109 + u^107 + 3*u^105 + 7*u^103 + 9*u^101 + 17*u^99 + 25*u^97 + 
#            31*u^95 + 47*u^93 + 61*u^91 + 72*u^89 + 96*u^87 + 113*u^85 + 
#            127*u^83 + 156*u^81 + 170*u^79 + 182*u^77 + 208*u^75 + 213*u^73 + 
#            219*u^71 + 236*u^69 + 226*u^67 + 223*u^65 + 227*u^63 + 204*u^61 + 
#            194*u^59 + 185*u^57 + 155*u^55 + 142*u^53 + 126*u^51 + 98*u^49 + 
#            86*u^47 + 69*u^45 + 49*u^43 + 41*u^41 + 29*u^39 + 18*u^37 + 14*u^
#            35 + 8*u^33 + 4*u^31 + 3*u^29 + u^27, 
#          u^83 + 5*u^81 + 8*u^79 + 15*u^77 + 28*u^75 + 39*u^73 + 56*u^71 + 
#            78*u^69 + 92*u^67 + 111*u^65 + 131*u^63 + 135*u^61 + 144*u^59 + 
#            148*u^57 + 135*u^55 + 130*u^53 + 119*u^51 + 96*u^49 + 85*u^47 + 
#            69*u^45 + 49*u^43 + 41*u^41 + 29*u^39 + 18*u^37 + 14*u^35 + 8*u^
#            33 + 4*u^31 + 3*u^29 + u^27, u^69 + 3*u^67 + 9*u^65 + 18*u^63 + 
#            28*u^61 + 45*u^59 + 59*u^57 + 68*u^55 + 80*u^53 + 82*u^51 + 76*u^
#            49 + 73*u^47 + 62*u^45 + 47*u^43 + 40*u^41 + 29*u^39 + 18*u^37 + 
#            14*u^35 + 8*u^33 + 4*u^31 + 3*u^29 + u^27, 
#          u^59 + 3*u^57 + 8*u^55 + 13*u^53 + 22*u^51 + 29*u^49 + 33*u^47 + 
#            39*u^45 + 34*u^43 + 32*u^41 + 27*u^39 + 17*u^37 + 14*u^35 + 8*u^
#            33 + 4*u^31 + 3*u^29 + u^27, u^58 + 3*u^56 + 3*u^54 + 5*u^52 + 
#            6*u^50 + 4*u^48 + 5*u^46 + 3*u^44 + u^42 + u^40, 
#          u^58 + 2*u^56 + 7*u^54 + 12*u^52 + 14*u^50 + 23*u^48 + 23*u^46 + 
#            21*u^44 + 23*u^42 + 14*u^40 + 10*u^38 + 7*u^36 + 2*u^34 + u^32, 
#          u^49 + 4*u^47 + 11*u^45 + 16*u^43 + 20*u^41 + 20*u^39 + 15*u^37 + 
#            13*u^35 + 8*u^33 + 4*u^31 + 3*u^29 + u^27, 
#          3*u^48 + 6*u^46 + 11*u^44 + 16*u^42 + 12*u^40 + 10*u^38 + 7*u^36 + 
#            2*u^34 + u^32, u^47 + 2*u^45 + 4*u^43 + 3*u^41 + u^39 + u^37, 
#          u^45 + 4*u^43 + 11*u^41 + 14*u^39 + 13*u^37 + 13*u^35 + 8*u^33 + 
#            4*u^31 + 3*u^29 + u^27, 2*u^43 + 5*u^41 + 8*u^39 + 11*u^37 + 12*u^
#            35 + 11*u^33 + 8*u^31 + 5*u^29 + 2*u^27, 
#          2*u^39 + 4*u^37 + 5*u^35 + 6*u^33 + 3*u^31 + 3*u^29 + u^27, 
#          u^38 + u^34, 2*u^39 + 3*u^37 + 3*u^35 + 5*u^33 + 2*u^31 + u^29, 
#          u^35 + 3*u^33 + 2*u^31 + 3*u^29 + u^27, u^36 + u^34 + u^32 + u^30, 
#          u^35 + 2*u^33 + 4*u^31 + 2*u^29 + u^27, 0*u^0, u^33 + u^29, 0*u^0, 
#          2*u^29 + u^27, u^31 + 2*u^29 + u^27, u^29, 0*u^0, u^29, 
#          2*u^29 + 2*u^27, u^29 + u^27, 0*u^0, u^28 + u^26, u^27, 0*u^0, 
#          u^26, 0*u^0, u^26, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^107 + u^105 + 2*u^103 + 6*u^101 + 7*u^99 + 11*u^97 + 20*u^95 + 23*u^
#            93 + 32*u^91 + 47*u^89 + 51*u^87 + 66*u^85 + 85*u^83 + 88*u^81 + 
#            107*u^79 + 125*u^77 + 124*u^75 + 143*u^73 + 154*u^71 + 146*u^69 + 
#            161*u^67 + 161*u^65 + 146*u^63 + 154*u^61 + 143*u^59 + 124*u^57 + 
#            125*u^55 + 107*u^53 + 88*u^51 + 85*u^49 + 66*u^47 + 51*u^45 + 
#            47*u^43 + 32*u^41 + 23*u^39 + 20*u^37 + 11*u^35 + 7*u^33 + 6*u^
#            31 + 2*u^29 + u^27 + u^25, u^83 + u^81 + 5*u^79 + 9*u^77 + 13*u^
#            75 + 24*u^73 + 35*u^71 + 42*u^69 + 61*u^67 + 72*u^65 + 79*u^63 + 
#            95*u^61 + 99*u^59 + 96*u^57 + 104*u^55 + 94*u^53 + 83*u^51 + 81*u^
#            49 + 65*u^47 + 51*u^45 + 47*u^43 + 32*u^41 + 23*u^39 + 20*u^37 + 
#            11*u^35 + 7*u^33 + 6*u^31 + 2*u^29 + u^27 + u^25, 
#          2*u^67 + 4*u^65 + 9*u^63 + 18*u^61 + 25*u^59 + 35*u^57 + 48*u^55 + 
#            51*u^53 + 55*u^51 + 60*u^49 + 52*u^47 + 46*u^45 + 43*u^43 + 31*u^
#            41 + 23*u^39 + 20*u^37 + 11*u^35 + 7*u^33 + 6*u^31 + 2*u^29 + u^
#            27 + u^25, u^59 + 2*u^57 + 5*u^55 + 11*u^53 + 13*u^51 + 21*u^49 + 
#            26*u^47 + 25*u^45 + 29*u^43 + 26*u^41 + 19*u^39 + 19*u^37 + 11*u^
#            35 + 7*u^33 + 6*u^31 + 2*u^29 + u^27 + u^25, 
#          u^54 + u^52 + u^50 + 2*u^48 + u^46 + u^44 + u^42, 
#          3*u^56 + 3*u^54 + 7*u^52 + 13*u^50 + 12*u^48 + 17*u^46 + 19*u^44 + 
#            13*u^42 + 14*u^40 + 10*u^38 + 4*u^36 + 4*u^34 + u^32, 
#          u^49 + 3*u^47 + 6*u^45 + 12*u^43 + 14*u^41 + 14*u^39 + 15*u^37 + 
#            10*u^35 + 7*u^33 + 6*u^31 + 2*u^29 + u^27 + u^25, 
#          4*u^46 + 7*u^44 + 8*u^42 + 11*u^40 + 9*u^38 + 4*u^36 + 4*u^34 + u^32
#            , u^45 + u^43 + u^41 + u^39, 3*u^43 + 5*u^41 + 9*u^39 + 12*u^37 + 
#            9*u^35 + 7*u^33 + 6*u^31 + 2*u^29 + u^27 + u^25, 
#          2*u^43 + 4*u^41 + 7*u^39 + 10*u^37 + 11*u^35 + 10*u^33 + 8*u^31 + 
#            5*u^29 + 2*u^27 + u^25, u^39 + 2*u^37 + 5*u^35 + 3*u^33 + 5*u^
#            31 + 2*u^29 + u^27 + u^25, u^36, 
#          u^37 + 4*u^35 + u^33 + 3*u^31 + u^29, 
#          2*u^35 + u^33 + 4*u^31 + 2*u^29 + u^27 + u^25, 0*u^0, 
#          u^35 + 3*u^33 + 2*u^31 + 3*u^29 + u^27, 0*u^0, u^31, 
#          u^35 + u^31 + u^29 + u^25, u^31 + u^29 + u^27 + u^25, 
#          2*u^29 + u^27, 0*u^0, 0*u^0, 0*u^0, u^29 + u^27 + u^25, 
#          u^29 + u^27 + u^25, 0*u^0, u^28 + u^26, u^25, 0*u^0, 0*u^0, u^25, 
#          0*u^0, u^25, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ 2*u^104 + 2*u^102 + 4*u^100 + 8*u^98 + 10*u^96 + 14*u^94 + 24*u^92 + 
#            26*u^90 + 36*u^88 + 48*u^86 + 52*u^84 + 64*u^82 + 80*u^80 + 80*u^
#            78 + 96*u^76 + 106*u^74 + 104*u^72 + 116*u^70 + 122*u^68 + 112*u^
#            66 + 122*u^64 + 116*u^62 + 104*u^60 + 106*u^58 + 96*u^56 + 80*u^
#            54 + 80*u^52 + 64*u^50 + 52*u^48 + 48*u^46 + 36*u^44 + 26*u^42 + 
#            24*u^40 + 14*u^38 + 10*u^36 + 8*u^34 + 4*u^32 + 2*u^30 + 2*u^28, 
#          u^82 + 3*u^80 + 4*u^78 + 10*u^76 + 15*u^74 + 20*u^72 + 31*u^70 + 
#            41*u^68 + 46*u^66 + 61*u^64 + 66*u^62 + 69*u^60 + 77*u^58 + 76*u^
#            56 + 69*u^54 + 72*u^52 + 60*u^50 + 51*u^48 + 47*u^46 + 36*u^44 + 
#            26*u^42 + 24*u^40 + 14*u^38 + 10*u^36 + 8*u^34 + 4*u^32 + 2*u^
#            30 + 2*u^28, u^68 + 2*u^66 + 6*u^64 + 9*u^62 + 16*u^60 + 24*u^
#            58 + 30*u^56 + 35*u^54 + 43*u^52 + 40*u^50 + 40*u^48 + 39*u^46 + 
#            32*u^44 + 25*u^42 + 23*u^40 + 14*u^38 + 10*u^36 + 8*u^34 + 4*u^
#            32 + 2*u^30 + 2*u^28, u^58 + 3*u^56 + 5*u^54 + 9*u^52 + 14*u^50 + 
#            15*u^48 + 20*u^46 + 21*u^44 + 17*u^42 + 19*u^40 + 13*u^38 + 9*u^
#            36 + 8*u^34 + 4*u^32 + 2*u^30 + 2*u^28, u^51, 
#          u^59 + u^57 + 3*u^55 + 7*u^53 + 8*u^51 + 11*u^49 + 16*u^47 + 12*u^
#            45 + 15*u^43 + 13*u^41 + 8*u^39 + 7*u^37 + 4*u^35 + u^33 + u^31, 
#          u^48 + 3*u^46 + 6*u^44 + 7*u^42 + 11*u^40 + 9*u^38 + 8*u^36 + 7*u^
#            34 + 4*u^32 + 2*u^30 + 2*u^28, u^49 + 3*u^47 + 4*u^45 + 8*u^43 + 
#            9*u^41 + 7*u^39 + 7*u^37 + 4*u^35 + u^33 + u^31, 0*u^0, 
#          u^44 + 2*u^42 + 6*u^40 + 6*u^38 + 7*u^36 + 7*u^34 + 4*u^32 + 2*u^
#            30 + 2*u^28, u^44 + 2*u^42 + 5*u^40 + 7*u^38 + 9*u^36 + 9*u^34 + 
#            8*u^32 + 5*u^30 + 3*u^28 + u^26, 
#          u^38 + u^36 + 3*u^34 + 3*u^32 + u^30 + 2*u^28, 0*u^0, 
#          2*u^38 + u^36 + 3*u^34 + 3*u^32 + u^30 + u^28, 
#          u^34 + 2*u^32 + u^30 + 2*u^28, 0*u^0, 
#          u^36 + u^34 + 3*u^32 + 3*u^30 + u^28 + u^26, 0*u^0, 
#          u^34 + u^32 + u^30 + u^28, u^32 + u^28, u^28, u^30 + u^28 + u^26, 
#          0*u^0, 0*u^0, u^30 + u^28, u^28 + u^26, u^28 + u^26, 0*u^0, 
#          u^29 + u^27 + u^25, 0*u^0, 0*u^0, u^25, 0*u^0, u^25, 0*u^0, u^25, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^104 + u^102 + 2*u^100 + 3*u^98 + 4*u^96 + 5*u^94 + 9*u^92 + 9*u^
#            90 + 13*u^88 + 16*u^86 + 18*u^84 + 21*u^82 + 27*u^80 + 26*u^78 + 
#            32*u^76 + 34*u^74 + 34*u^72 + 37*u^70 + 40*u^68 + 36*u^66 + 40*u^
#            64 + 37*u^62 + 34*u^60 + 34*u^58 + 32*u^56 + 26*u^54 + 27*u^52 + 
#            21*u^50 + 18*u^48 + 16*u^46 + 13*u^44 + 9*u^42 + 9*u^40 + 5*u^
#            38 + 4*u^36 + 3*u^34 + 2*u^32 + u^30 + u^28, 
#          u^80 + u^78 + 3*u^76 + 4*u^74 + 6*u^72 + 9*u^70 + 13*u^68 + 14*u^
#            66 + 20*u^64 + 21*u^62 + 23*u^60 + 25*u^58 + 26*u^56 + 23*u^54 + 
#            25*u^52 + 20*u^50 + 18*u^48 + 16*u^46 + 13*u^44 + 9*u^42 + 9*u^
#            40 + 5*u^38 + 4*u^36 + 3*u^34 + 2*u^32 + u^30 + u^28, 
#          2*u^64 + 2*u^62 + 5*u^60 + 8*u^58 + 11*u^56 + 12*u^54 + 16*u^52 + 
#            14*u^50 + 15*u^48 + 14*u^46 + 12*u^44 + 9*u^42 + 9*u^40 + 5*u^
#            38 + 4*u^36 + 3*u^34 + 2*u^32 + u^30 + u^28, 
#          u^56 + u^54 + 2*u^52 + 4*u^50 + 4*u^48 + 6*u^46 + 7*u^44 + 6*u^42 + 
#            7*u^40 + 5*u^38 + 4*u^36 + 3*u^34 + 2*u^32 + u^30 + u^28, 
#          u^55 + u^53 + 3*u^51 + 2*u^49 + 3*u^47 + 2*u^45 + 2*u^43 + u^41 + u^
#            39, u^53 + u^51 + 2*u^49 + 3*u^47 + 2*u^45 + 4*u^43 + 3*u^41 + 
#            2*u^39 + 2*u^37 + u^35, u^46 + 2*u^44 + 3*u^42 + 5*u^40 + 4*u^
#            38 + 4*u^36 + 3*u^34 + 2*u^32 + u^30 + u^28, 
#          2*u^43 + 2*u^41 + 2*u^39 + 2*u^37 + u^35, 
#          u^48 + u^46 + 2*u^44 + 2*u^42 + 2*u^40 + u^38 + u^36, 
#          u^42 + 3*u^40 + 3*u^38 + 4*u^36 + 3*u^34 + 2*u^32 + u^30 + u^28, 
#          u^40 + u^38 + 2*u^36 + 2*u^34 + 2*u^32 + u^30 + u^28, 
#          u^36 + u^34 + 2*u^32 + u^30 + u^28, u^39 + u^37 + u^35 + u^33, 
#          u^34 + u^32, u^32 + u^28, u^33 + u^31 + u^29, 0*u^0, 0*u^0, 0*u^0, 
#          u^32 + u^28, u^28, 0*u^0, u^30 + u^28, 0*u^0, 0*u^0, u^28, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^27, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^25, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^108 + u^106 + 4*u^104 + 6*u^102 + 10*u^100 + 16*u^98 + 24*u^96 + 
#            31*u^94 + 47*u^92 + 57*u^90 + 74*u^88 + 93*u^86 + 110*u^84 + 
#            128*u^82 + 154*u^80 + 165*u^78 + 188*u^76 + 204*u^74 + 213*u^72 + 
#            225*u^70 + 235*u^68 + 228*u^66 + 235*u^64 + 225*u^62 + 213*u^60 + 
#            204*u^58 + 188*u^56 + 165*u^54 + 154*u^52 + 128*u^50 + 110*u^48 + 
#            93*u^46 + 74*u^44 + 57*u^42 + 47*u^40 + 31*u^38 + 24*u^36 + 16*u^
#            34 + 10*u^32 + 6*u^30 + 4*u^28 + u^26 + u^24, 
#          u^82 + 4*u^80 + 7*u^78 + 16*u^76 + 26*u^74 + 39*u^72 + 57*u^70 + 
#            77*u^68 + 93*u^66 + 117*u^64 + 130*u^62 + 142*u^60 + 151*u^58 + 
#            152*u^56 + 144*u^54 + 141*u^52 + 122*u^50 + 108*u^48 + 92*u^46 + 
#            74*u^44 + 57*u^42 + 47*u^40 + 31*u^38 + 24*u^36 + 16*u^34 + 10*u^
#            32 + 6*u^30 + 4*u^28 + u^26 + u^24, 
#          u^68 + 3*u^66 + 10*u^64 + 17*u^62 + 32*u^60 + 47*u^58 + 63*u^56 + 
#            75*u^54 + 88*u^52 + 86*u^50 + 87*u^48 + 79*u^46 + 68*u^44 + 55*u^
#            42 + 46*u^40 + 31*u^38 + 24*u^36 + 16*u^34 + 10*u^32 + 6*u^30 + 
#            4*u^28 + u^26 + u^24, u^58 + 3*u^56 + 8*u^54 + 14*u^52 + 23*u^
#            50 + 31*u^48 + 37*u^46 + 42*u^44 + 39*u^42 + 38*u^40 + 29*u^38 + 
#            23*u^36 + 16*u^34 + 10*u^32 + 6*u^30 + 4*u^28 + u^26 + u^24, 
#          2*u^57 + 4*u^55 + 5*u^53 + 10*u^51 + 8*u^49 + 9*u^47 + 8*u^45 + 5*u^
#            43 + 3*u^41 + 2*u^39, 2*u^55 + 5*u^53 + 10*u^51 + 14*u^49 + 22*u^
#            47 + 21*u^45 + 24*u^43 + 21*u^41 + 15*u^39 + 11*u^37 + 6*u^35 + 
#            2*u^33 + u^31, u^48 + 5*u^46 + 14*u^44 + 20*u^42 + 25*u^40 + 23*u^
#            38 + 21*u^36 + 15*u^34 + 10*u^32 + 6*u^30 + 4*u^28 + u^26 + u^24, 
#          2*u^47 + 6*u^45 + 13*u^43 + 15*u^41 + 13*u^39 + 11*u^37 + 6*u^35 + 
#            2*u^33 + u^31, u^48 + 2*u^46 + 6*u^44 + 7*u^42 + 5*u^40 + 3*u^
#            38 + 2*u^36, u^44 + 7*u^42 + 15*u^40 + 17*u^38 + 19*u^36 + 15*u^
#            34 + 10*u^32 + 6*u^30 + 4*u^28 + u^26 + u^24, 
#          2*u^42 + 5*u^40 + 8*u^38 + 12*u^36 + 13*u^34 + 12*u^32 + 9*u^30 + 
#            6*u^28 + 2*u^26 + u^24, u^40 + 3*u^38 + 6*u^36 + 7*u^34 + 8*u^
#            32 + 5*u^30 + 4*u^28 + u^26 + u^24, 
#          u^39 + 2*u^37 + 2*u^35 + 2*u^33, 
#          2*u^38 + 3*u^36 + 5*u^34 + 4*u^32 + 2*u^30 + u^28, 
#          2*u^34 + 3*u^32 + 3*u^30 + 4*u^28 + u^26 + u^24, 
#          u^35 + 2*u^33 + 3*u^31 + 2*u^29, 
#          u^34 + 2*u^32 + 3*u^30 + 2*u^28 + u^26, u^34, 0*u^0, 0*u^0, 
#          u^30 + 3*u^28 + u^26 + u^24, u^30 + 2*u^28 + u^26, 2*u^30 + 2*u^28, 
#          0*u^0, 0*u^0, u^30 + 3*u^28 + 2*u^26 + u^24, u^28 + u^26 + u^24, 
#          0*u^0, u^27 + u^25, u^26 + u^24, u^27, u^25, u^24, u^25, 0*u^0, 
#          0*u^0, 0*u^0, u^24, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^107 + 2*u^105 + 3*u^103 + 7*u^101 + 10*u^99 + 14*u^97 + 24*u^95 + 
#            30*u^93 + 39*u^91 + 56*u^89 + 64*u^87 + 79*u^85 + 101*u^83 + 
#            108*u^81 + 127*u^79 + 148*u^77 + 150*u^75 + 169*u^73 + 182*u^71 + 
#            176*u^69 + 190*u^67 + 190*u^65 + 176*u^63 + 182*u^61 + 169*u^59 + 
#            150*u^57 + 148*u^55 + 127*u^53 + 108*u^51 + 101*u^49 + 79*u^47 + 
#            64*u^45 + 56*u^43 + 39*u^41 + 30*u^39 + 24*u^37 + 14*u^35 + 10*u^
#            33 + 7*u^31 + 3*u^29 + 2*u^27 + u^25, 
#          u^83 + u^81 + 5*u^79 + 10*u^77 + 15*u^75 + 27*u^73 + 40*u^71 + 50*u^
#            69 + 71*u^67 + 85*u^65 + 95*u^63 + 113*u^61 + 118*u^59 + 117*u^
#            57 + 124*u^55 + 113*u^53 + 102*u^51 + 97*u^49 + 78*u^47 + 64*u^
#            45 + 56*u^43 + 39*u^41 + 30*u^39 + 24*u^37 + 14*u^35 + 10*u^33 + 
#            7*u^31 + 3*u^29 + 2*u^27 + u^25, 2*u^67 + 4*u^65 + 10*u^63 + 21*u^
#            61 + 30*u^59 + 43*u^57 + 58*u^55 + 63*u^53 + 69*u^51 + 73*u^49 + 
#            64*u^47 + 58*u^45 + 52*u^43 + 38*u^41 + 30*u^39 + 24*u^37 + 14*u^
#            35 + 10*u^33 + 7*u^31 + 3*u^29 + 2*u^27 + u^25, 
#          2*u^57 + 4*u^55 + 10*u^53 + 15*u^51 + 22*u^49 + 29*u^47 + 30*u^45 + 
#            34*u^43 + 31*u^41 + 25*u^39 + 23*u^37 + 14*u^35 + 10*u^33 + 7*u^
#            31 + 3*u^29 + 2*u^27 + u^25, u^58 + u^56 + 4*u^54 + 5*u^52 + 5*u^
#            50 + 8*u^48 + 5*u^46 + 5*u^44 + 4*u^42 + u^40 + u^38, 
#          2*u^56 + 3*u^54 + 6*u^52 + 13*u^50 + 13*u^48 + 18*u^46 + 21*u^44 + 
#            15*u^42 + 16*u^40 + 11*u^38 + 5*u^36 + 4*u^34 + u^32, 
#          3*u^47 + 7*u^45 + 14*u^43 + 18*u^41 + 19*u^39 + 19*u^37 + 13*u^35 + 
#            10*u^33 + 7*u^31 + 3*u^29 + 2*u^27 + u^25, 
#          4*u^46 + 8*u^44 + 9*u^42 + 13*u^40 + 10*u^38 + 5*u^36 + 4*u^34 + u^
#            32, u^47 + 4*u^45 + 4*u^43 + 5*u^41 + 4*u^39 + u^37 + u^35, 
#          4*u^43 + 7*u^41 + 13*u^39 + 16*u^37 + 12*u^35 + 10*u^33 + 7*u^31 + 
#            3*u^29 + 2*u^27 + u^25, u^43 + 3*u^41 + 6*u^39 + 9*u^37 + 11*u^
#            35 + 11*u^33 + 9*u^31 + 6*u^29 + 3*u^27 + u^25, 
#          u^39 + 3*u^37 + 6*u^35 + 5*u^33 + 6*u^31 + 3*u^29 + 2*u^27 + u^25, 
#          u^40 + u^38 + 3*u^36 + u^34 + u^32, 
#          u^37 + 4*u^35 + 2*u^33 + 3*u^31 + u^29, 
#          u^35 + u^33 + 4*u^31 + 2*u^29 + 2*u^27 + u^25, 
#          2*u^34 + 2*u^32 + u^30 + u^28, 2*u^33 + 2*u^31 + 3*u^29 + u^27, 
#          0*u^0, u^31, 0*u^0, u^31 + u^29 + 2*u^27 + u^25, 2*u^29 + u^27, 
#          u^31 + u^29 + u^27, 0*u^0, 0*u^0, u^29 + 2*u^27 + u^25, 
#          u^27 + u^25, 0*u^0, u^28 + u^26 + u^24, u^25, u^26, 0*u^0, u^25, 
#          u^24, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^24, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^106 + 2*u^102 + 2*u^100 + 3*u^98 + 5*u^96 + 8*u^94 + 7*u^92 + 14*u^
#            90 + 14*u^88 + 18*u^86 + 23*u^84 + 27*u^82 + 27*u^80 + 37*u^78 + 
#            35*u^76 + 40*u^74 + 44*u^72 + 45*u^70 + 43*u^68 + 50*u^66 + 43*u^
#            64 + 45*u^62 + 44*u^60 + 40*u^58 + 35*u^56 + 37*u^54 + 27*u^52 + 
#            27*u^50 + 23*u^48 + 18*u^46 + 14*u^44 + 14*u^42 + 7*u^40 + 8*u^
#            38 + 5*u^36 + 3*u^34 + 2*u^32 + 2*u^30 + u^26, 
#          2*u^78 + 2*u^76 + 5*u^74 + 8*u^72 + 11*u^70 + 14*u^68 + 21*u^66 + 
#            21*u^64 + 27*u^62 + 29*u^60 + 30*u^58 + 29*u^56 + 32*u^54 + 25*u^
#            52 + 26*u^50 + 22*u^48 + 18*u^46 + 14*u^44 + 14*u^42 + 7*u^40 + 
#            8*u^38 + 5*u^36 + 3*u^34 + 2*u^32 + 2*u^30 + u^26, 
#          u^66 + u^64 + 5*u^62 + 6*u^60 + 10*u^58 + 12*u^56 + 17*u^54 + 15*u^
#            52 + 20*u^50 + 17*u^48 + 16*u^46 + 13*u^44 + 13*u^42 + 7*u^40 + 
#            8*u^38 + 5*u^36 + 3*u^34 + 2*u^32 + 2*u^30 + u^26, 
#          u^54 + 2*u^52 + 3*u^50 + 6*u^48 + 7*u^46 + 7*u^44 + 10*u^42 + 6*u^
#            40 + 7*u^38 + 5*u^36 + 3*u^34 + 2*u^32 + 2*u^30 + u^26, 
#          u^57 + u^55 + 3*u^53 + 2*u^51 + 4*u^49 + 2*u^47 + 3*u^45 + u^43 + u^
#            41, u^51 + 3*u^49 + 2*u^47 + 6*u^45 + 4*u^43 + 4*u^41 + 4*u^39 + 
#            2*u^37 + u^35 + u^33, u^46 + 2*u^44 + 5*u^42 + 4*u^40 + 6*u^38 + 
#            4*u^36 + 3*u^34 + 2*u^32 + 2*u^30 + u^26, 
#          2*u^45 + 2*u^43 + 3*u^41 + 4*u^39 + 2*u^37 + u^35 + u^33, 
#          u^46 + u^44 + 2*u^42 + u^40 + u^38, 
#          2*u^42 + 2*u^40 + 5*u^38 + 4*u^36 + 3*u^34 + 2*u^32 + 2*u^30 + u^26,
#          u^38 + u^36 + 2*u^34 + 2*u^32 + 2*u^30 + u^28 + u^26, 
#          u^38 + 2*u^36 + 2*u^34 + u^32 + 2*u^30 + u^26, 0*u^0, 
#          u^40 + 2*u^36 + u^34 + u^32 + u^30, u^30 + u^26, u^33 + u^31 + u^29,
#          u^28, u^34, 0*u^0, 0*u^0, u^26, u^28, u^28, 0*u^0, 0*u^0, 
#          2*u^28 + u^26, 0*u^0, u^28, 0*u^0, u^26, 0*u^0, u^27, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^24, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^107 + 2*u^105 + 4*u^103 + 7*u^101 + 11*u^99 + 18*u^97 + 26*u^95 + 
#            36*u^93 + 49*u^91 + 63*u^89 + 81*u^87 + 100*u^85 + 119*u^83 + 
#            141*u^81 + 161*u^79 + 182*u^77 + 202*u^75 + 217*u^73 + 232*u^71 + 
#            242*u^69 + 248*u^67 + 252*u^65 + 248*u^63 + 242*u^61 + 232*u^59 + 
#            217*u^57 + 202*u^55 + 182*u^53 + 161*u^51 + 141*u^49 + 119*u^47 + 
#            100*u^45 + 81*u^43 + 63*u^41 + 49*u^39 + 36*u^37 + 26*u^35 + 18*u^
#            33 + 11*u^31 + 7*u^29 + 4*u^27 + 2*u^25 + u^23, 
#          u^81 + 4*u^79 + 9*u^77 + 18*u^75 + 29*u^73 + 45*u^71 + 64*u^69 + 
#            85*u^67 + 107*u^65 + 127*u^63 + 144*u^61 + 158*u^59 + 164*u^57 + 
#            166*u^55 + 160*u^53 + 149*u^51 + 135*u^49 + 117*u^47 + 99*u^45 + 
#            81*u^43 + 63*u^41 + 49*u^39 + 36*u^37 + 26*u^35 + 18*u^33 + 11*u^
#            31 + 7*u^29 + 4*u^27 + 2*u^25 + u^23, 
#          u^67 + 5*u^65 + 12*u^63 + 23*u^61 + 38*u^59 + 55*u^57 + 73*u^55 + 
#            87*u^53 + 96*u^51 + 99*u^49 + 95*u^47 + 87*u^45 + 75*u^43 + 61*u^
#            41 + 48*u^39 + 36*u^37 + 26*u^35 + 18*u^33 + 11*u^31 + 7*u^29 + 
#            4*u^27 + 2*u^25 + u^23, u^57 + 4*u^55 + 9*u^53 + 17*u^51 + 26*u^
#            49 + 36*u^47 + 43*u^45 + 47*u^43 + 46*u^41 + 41*u^39 + 34*u^37 + 
#            25*u^35 + 18*u^33 + 11*u^31 + 7*u^29 + 4*u^27 + 2*u^25 + u^23, 
#          u^58 + 3*u^56 + 6*u^54 + 9*u^52 + 11*u^50 + 12*u^48 + 11*u^46 + 9*u^
#            44 + 6*u^42 + 3*u^40 + u^38, 2*u^54 + 5*u^52 + 11*u^50 + 17*u^
#            48 + 22*u^46 + 26*u^44 + 25*u^42 + 23*u^40 + 17*u^38 + 11*u^36 + 
#            6*u^34 + 2*u^32 + u^30, 2*u^47 + 8*u^45 + 17*u^43 + 25*u^41 + 
#            29*u^39 + 28*u^37 + 23*u^35 + 17*u^33 + 11*u^31 + 7*u^29 + 4*u^
#            27 + 2*u^25 + u^23, 3*u^46 + 9*u^44 + 15*u^42 + 18*u^40 + 16*u^
#            38 + 11*u^36 + 6*u^34 + 2*u^32 + u^30, 
#          u^47 + 4*u^45 + 7*u^43 + 8*u^41 + 6*u^39 + 3*u^37 + u^35, 
#          3*u^43 + 10*u^41 + 19*u^39 + 23*u^37 + 22*u^35 + 17*u^33 + 11*u^
#            31 + 7*u^29 + 4*u^27 + 2*u^25 + u^23, 
#          2*u^41 + 6*u^39 + 10*u^37 + 14*u^35 + 15*u^33 + 13*u^31 + 10*u^29 + 
#            6*u^27 + 3*u^25 + u^23, 2*u^39 + 5*u^37 + 9*u^35 + 10*u^33 + 9*u^
#            31 + 6*u^29 + 4*u^27 + 2*u^25 + u^23, 
#          u^38 + 2*u^36 + 2*u^34 + u^32, u^39 + 3*u^37 + 6*u^35 + 6*u^33 + 
#            5*u^31 + 2*u^29 + u^27, u^35 + 2*u^33 + 4*u^31 + 4*u^29 + 4*u^
#            27 + 2*u^25 + u^23, 2*u^34 + 4*u^32 + 4*u^30 + 2*u^28, 
#          u^33 + 2*u^31 + 3*u^29 + 2*u^27 + u^25, u^35 + u^33 + u^31, 0*u^0, 
#          0*u^0, 2*u^29 + 3*u^27 + 2*u^25 + u^23, 2*u^29 + 2*u^27 + u^25, 
#          3*u^29 + 2*u^27, u^29, 0*u^0, 3*u^29 + 5*u^27 + 3*u^25 + u^23, 
#          u^27 + 2*u^25 + u^23, u^29 + u^27, u^26 + u^24, 
#          u^27 + 2*u^25 + u^23, u^26, u^26 + u^24, u^25 + u^23, u^24, 0*u^0, 
#          0*u^0, 0*u^0, u^23, 0*u^0, u^23, u^23, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^106 + 2*u^104 + 5*u^102 + 8*u^100 + 15*u^98 + 22*u^96 + 34*u^94 + 
#            46*u^92 + 65*u^90 + 83*u^88 + 108*u^86 + 131*u^84 + 161*u^82 + 
#            187*u^80 + 218*u^78 + 243*u^76 + 272*u^74 + 292*u^72 + 314*u^70 + 
#            325*u^68 + 337*u^66 + 337*u^64 + 337*u^62 + 325*u^60 + 314*u^58 + 
#            292*u^56 + 272*u^54 + 243*u^52 + 218*u^50 + 187*u^48 + 161*u^46 + 
#            131*u^44 + 108*u^42 + 83*u^40 + 65*u^38 + 46*u^36 + 34*u^34 + 
#            22*u^32 + 15*u^30 + 8*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          u^82 + 3*u^80 + 8*u^78 + 15*u^76 + 28*u^74 + 44*u^72 + 67*u^70 + 
#            91*u^68 + 121*u^66 + 148*u^64 + 177*u^62 + 197*u^60 + 216*u^58 + 
#            222*u^56 + 225*u^54 + 214*u^52 + 202*u^50 + 179*u^48 + 158*u^46 + 
#            130*u^44 + 108*u^42 + 83*u^40 + 65*u^38 + 46*u^36 + 34*u^34 + 
#            22*u^32 + 15*u^30 + 8*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          3*u^66 + 8*u^64 + 20*u^62 + 34*u^60 + 57*u^58 + 78*u^56 + 103*u^
#            54 + 118*u^52 + 132*u^50 + 132*u^48 + 129*u^46 + 114*u^44 + 100*u^
#            42 + 80*u^40 + 64*u^38 + 46*u^36 + 34*u^34 + 22*u^32 + 15*u^30 + 
#            8*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          u^58 + 4*u^56 + 10*u^54 + 19*u^52 + 31*u^50 + 44*u^48 + 56*u^46 + 
#            63*u^44 + 66*u^42 + 62*u^40 + 55*u^38 + 43*u^36 + 33*u^34 + 22*u^
#            32 + 15*u^30 + 8*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          u^57 + 2*u^55 + 5*u^53 + 6*u^51 + 9*u^49 + 8*u^47 + 9*u^45 + 6*u^
#            43 + 5*u^41 + 2*u^39 + u^37, u^57 + 3*u^55 + 8*u^53 + 14*u^51 + 
#            23*u^49 + 30*u^47 + 37*u^45 + 38*u^43 + 37*u^41 + 30*u^39 + 23*u^
#            37 + 14*u^35 + 8*u^33 + 3*u^31 + u^29, 
#          u^48 + 6*u^46 + 15*u^44 + 27*u^42 + 35*u^40 + 39*u^38 + 35*u^36 + 
#            30*u^34 + 21*u^32 + 15*u^30 + 8*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          2*u^47 + 8*u^45 + 16*u^43 + 23*u^41 + 24*u^39 + 21*u^37 + 14*u^35 + 
#            8*u^33 + 3*u^31 + u^29, 2*u^46 + 4*u^44 + 7*u^42 + 6*u^40 + 5*u^
#            38 + 2*u^36 + u^34, u^44 + 7*u^42 + 16*u^40 + 26*u^38 + 29*u^36 + 
#            28*u^34 + 21*u^32 + 15*u^30 + 8*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          3*u^42 + 8*u^40 + 16*u^38 + 22*u^36 + 27*u^34 + 26*u^32 + 22*u^30 + 
#            14*u^28 + 8*u^26 + 3*u^24 + u^22, 
#          3*u^38 + 6*u^36 + 12*u^34 + 12*u^32 + 12*u^30 + 7*u^28 + 5*u^26 + 
#            2*u^24 + u^22, u^39 + 2*u^37 + 3*u^35 + 2*u^33 + u^31, 
#          u^38 + 4*u^36 + 7*u^34 + 8*u^32 + 6*u^30 + 3*u^28 + u^26, 
#          3*u^34 + 5*u^32 + 8*u^30 + 6*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          u^35 + 2*u^33 + 3*u^31 + 2*u^29 + u^27, 
#          2*u^34 + 5*u^32 + 7*u^30 + 6*u^28 + 3*u^26 + u^24, 0*u^0, 
#          u^32 + u^30 + u^28, u^34 + u^32 + 2*u^30 + u^28 + 2*u^26 + u^24 + u^
#            22, 2*u^30 + 3*u^28 + 4*u^26 + 2*u^24 + u^22, 
#          2*u^30 + 4*u^28 + 3*u^26 + u^24, u^30 + 2*u^28 + u^26, 0*u^0, u^28, 
#          u^30 + 3*u^28 + 5*u^26 + 3*u^24 + u^22, 
#          u^28 + 3*u^26 + 3*u^24 + u^22, 0*u^0, 2*u^27 + 3*u^25 + 2*u^23, 
#          u^26 + 2*u^24 + u^22, u^25, u^25 + u^23, 2*u^24 + u^22, 
#          u^25 + 2*u^23, u^24 + u^22, u^23, u^23, u^22, u^23, 0*u^0, u^22, 
#          u^22, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^106 + u^104 + 3*u^102 + 4*u^100 + 8*u^98 + 11*u^96 + 17*u^94 + 22*u^
#            92 + 32*u^90 + 39*u^88 + 52*u^86 + 61*u^84 + 76*u^82 + 87*u^80 + 
#            102*u^78 + 112*u^76 + 127*u^74 + 134*u^72 + 146*u^70 + 149*u^68 + 
#            156*u^66 + 155*u^64 + 156*u^62 + 149*u^60 + 146*u^58 + 134*u^56 + 
#            127*u^54 + 112*u^52 + 102*u^50 + 87*u^48 + 76*u^46 + 61*u^44 + 
#            52*u^42 + 39*u^40 + 32*u^38 + 22*u^36 + 17*u^34 + 11*u^32 + 8*u^
#            30 + 4*u^28 + 3*u^26 + u^24 + u^22, 
#          u^80 + 3*u^78 + 6*u^76 + 12*u^74 + 19*u^72 + 30*u^70 + 41*u^68 + 
#            55*u^66 + 68*u^64 + 82*u^62 + 91*u^60 + 101*u^58 + 103*u^56 + 
#            106*u^54 + 100*u^52 + 95*u^50 + 84*u^48 + 75*u^46 + 61*u^44 + 
#            52*u^42 + 39*u^40 + 32*u^38 + 22*u^36 + 17*u^34 + 11*u^32 + 8*u^
#            30 + 4*u^28 + 3*u^26 + u^24 + u^22, 
#          u^66 + 3*u^64 + 9*u^62 + 15*u^60 + 27*u^58 + 36*u^56 + 50*u^54 + 
#            56*u^52 + 64*u^50 + 63*u^48 + 63*u^46 + 54*u^44 + 49*u^42 + 38*u^
#            40 + 32*u^38 + 22*u^36 + 17*u^34 + 11*u^32 + 8*u^30 + 4*u^28 + 
#            3*u^26 + u^24 + u^22, u^56 + 3*u^54 + 7*u^52 + 12*u^50 + 18*u^
#            48 + 25*u^46 + 28*u^44 + 32*u^42 + 29*u^40 + 28*u^38 + 21*u^36 + 
#            17*u^34 + 11*u^32 + 8*u^30 + 4*u^28 + 3*u^26 + u^24 + u^22, 
#          u^57 + 2*u^55 + 5*u^53 + 6*u^51 + 9*u^49 + 8*u^47 + 9*u^45 + 6*u^
#            43 + 5*u^41 + 2*u^39 + u^37, 2*u^53 + 3*u^51 + 8*u^49 + 10*u^47 + 
#            15*u^45 + 15*u^43 + 16*u^41 + 13*u^39 + 10*u^37 + 6*u^35 + 3*u^
#            33 + u^31, 2*u^46 + 6*u^44 + 14*u^42 + 18*u^40 + 21*u^38 + 18*u^
#            36 + 16*u^34 + 11*u^32 + 8*u^30 + 4*u^28 + 3*u^26 + u^24 + u^22, 
#          2*u^45 + 6*u^43 + 10*u^41 + 11*u^39 + 9*u^37 + 6*u^35 + 3*u^33 + u^
#            31, 2*u^46 + 4*u^44 + 7*u^42 + 6*u^40 + 5*u^38 + 2*u^36 + u^34, 
#          3*u^42 + 9*u^40 + 15*u^38 + 16*u^36 + 15*u^34 + 11*u^32 + 8*u^30 + 
#            4*u^28 + 3*u^26 + u^24 + u^22, 2*u^40 + 5*u^38 + 7*u^36 + 10*u^
#            34 + 10*u^32 + 9*u^30 + 6*u^28 + 4*u^26 + u^24 + u^22, 
#          2*u^38 + 4*u^36 + 7*u^34 + 7*u^32 + 7*u^30 + 4*u^28 + 3*u^26 + u^
#            24 + u^22, u^39 + 2*u^37 + 3*u^35 + 2*u^33 + u^31, 
#          2*u^36 + 3*u^34 + 4*u^32 + 2*u^30 + u^28, 
#          u^34 + 2*u^32 + 3*u^30 + 3*u^28 + 3*u^26 + u^24 + u^22, 
#          2*u^33 + 3*u^31 + 3*u^29 + u^27, u^32 + u^30 + 2*u^28 + u^26, 
#          u^34 + u^32 + u^30, 0*u^0, 0*u^0, u^30 + 2*u^28 + 3*u^26 + u^24 + u^
#            22, u^28 + u^26, u^30 + 3*u^28 + u^26, u^30 + u^28, 0*u^0, 
#          3*u^28 + 3*u^26 + u^24 + u^22, u^26 + u^24 + u^22, u^28, 
#          u^25 + u^23, u^26 + u^24 + u^22, u^27 + u^25, 0*u^0, u^24 + u^22, 
#          u^23, 0*u^0, 0*u^0, 0*u^0, u^22, u^23, 0*u^0, u^22, 0*u^0, u^22, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^105 + 2*u^103 + 4*u^101 + 8*u^99 + 13*u^97 + 20*u^95 + 31*u^93 + 
#            43*u^91 + 58*u^89 + 78*u^87 + 98*u^85 + 121*u^83 + 148*u^81 + 
#            173*u^79 + 199*u^77 + 227*u^75 + 249*u^73 + 270*u^71 + 290*u^69 + 
#            301*u^67 + 309*u^65 + 314*u^63 + 309*u^61 + 301*u^59 + 290*u^57 + 
#            270*u^55 + 249*u^53 + 227*u^51 + 199*u^49 + 173*u^47 + 148*u^45 + 
#            121*u^43 + 98*u^41 + 78*u^39 + 58*u^37 + 43*u^35 + 31*u^33 + 20*u^
#            31 + 13*u^29 + 8*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^81 + 3*u^79 + 8*u^77 + 16*u^75 + 28*u^73 + 44*u^71 + 66*u^69 + 
#            89*u^67 + 116*u^65 + 143*u^63 + 167*u^61 + 187*u^59 + 203*u^57 + 
#            208*u^55 + 208*u^53 + 201*u^51 + 185*u^49 + 166*u^47 + 145*u^45 + 
#            120*u^43 + 98*u^41 + 78*u^39 + 58*u^37 + 43*u^35 + 31*u^33 + 20*u^
#            31 + 13*u^29 + 8*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^67 + 4*u^65 + 10*u^63 + 21*u^61 + 37*u^59 + 57*u^57 + 78*u^55 + 
#            99*u^53 + 115*u^51 + 123*u^49 + 125*u^47 + 119*u^45 + 106*u^43 + 
#            91*u^41 + 75*u^39 + 57*u^37 + 43*u^35 + 31*u^33 + 20*u^31 + 13*u^
#            29 + 8*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          2*u^57 + 5*u^55 + 11*u^53 + 22*u^51 + 33*u^49 + 45*u^47 + 57*u^45 + 
#            62*u^43 + 63*u^41 + 60*u^39 + 50*u^37 + 40*u^35 + 30*u^33 + 20*u^
#            31 + 13*u^29 + 8*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^56 + 2*u^54 + 4*u^52 + 5*u^50 + 6*u^48 + 6*u^46 + 5*u^44 + 4*u^
#            42 + 2*u^40 + u^38, u^56 + 4*u^54 + 9*u^52 + 16*u^50 + 25*u^48 + 
#            32*u^46 + 37*u^44 + 39*u^42 + 35*u^40 + 29*u^38 + 21*u^36 + 13*u^
#            34 + 7*u^32 + 3*u^30 + u^28, 2*u^47 + 8*u^45 + 17*u^43 + 28*u^
#            41 + 35*u^39 + 36*u^37 + 33*u^35 + 27*u^33 + 19*u^31 + 13*u^29 + 
#            8*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^48 + 4*u^46 + 11*u^44 + 19*u^42 + 24*u^40 + 24*u^38 + 20*u^36 + 
#            13*u^34 + 7*u^32 + 3*u^30 + u^28, 
#          u^45 + 3*u^43 + 4*u^41 + 4*u^39 + 2*u^37 + u^35, 
#          2*u^43 + 9*u^41 + 18*u^39 + 26*u^37 + 28*u^35 + 26*u^33 + 19*u^31 + 
#            13*u^29 + 8*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^43 + 5*u^41 + 11*u^39 + 18*u^37 + 25*u^35 + 28*u^33 + 26*u^31 + 
#            21*u^29 + 14*u^27 + 7*u^25 + 3*u^23 + u^21, 
#          u^39 + 4*u^37 + 8*u^35 + 12*u^33 + 12*u^31 + 10*u^29 + 7*u^27 + 4*u^
#            25 + 2*u^23 + u^21, u^36 + u^34 + u^32, 
#          u^39 + 3*u^37 + 6*u^35 + 9*u^33 + 8*u^31 + 6*u^29 + 3*u^27 + u^25, 
#          u^35 + 4*u^33 + 6*u^31 + 8*u^29 + 7*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^34 + 2*u^32 + 2*u^30 + u^28, u^35 + 3*u^33 + 6*u^31 + 8*u^29 + 
#            6*u^27 + 3*u^25 + u^23, 0*u^0, u^31 + u^29 + u^27, 
#          u^33 + u^31 + u^29 + 2*u^27 + u^25 + u^23 + u^21, 
#          2*u^29 + 4*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          3*u^29 + 5*u^27 + 3*u^25 + u^23, u^29 + u^27, 0*u^0, u^29 + u^27, 
#          2*u^29 + 5*u^27 + 5*u^25 + 3*u^23 + u^21, 
#          u^29 + 3*u^27 + 4*u^25 + 3*u^23 + u^21, 0*u^0, 
#          u^28 + 3*u^26 + 3*u^24 + u^22, 2*u^25 + 2*u^23 + u^21, 0*u^0, 
#          u^26 + 2*u^24 + u^22, u^25 + 2*u^23 + u^21, 2*u^24 + u^22, 
#          u^23 + u^21, u^24 + u^22, 0*u^0, u^23 + u^21, 0*u^0, 0*u^0, u^21, 
#          u^21, 0*u^0, u^21, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^105 + 2*u^103 + 3*u^101 + 7*u^99 + 11*u^97 + 15*u^95 + 25*u^93 + 
#            34*u^91 + 43*u^89 + 61*u^87 + 75*u^85 + 89*u^83 + 114*u^81 + 
#            130*u^79 + 146*u^77 + 173*u^75 + 185*u^73 + 198*u^71 + 220*u^69 + 
#            222*u^67 + 227*u^65 + 238*u^63 + 227*u^61 + 222*u^59 + 220*u^57 + 
#            198*u^55 + 185*u^53 + 173*u^51 + 146*u^49 + 130*u^47 + 114*u^45 + 
#            89*u^43 + 75*u^41 + 61*u^39 + 43*u^37 + 34*u^35 + 25*u^33 + 15*u^
#            31 + 11*u^29 + 7*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          u^81 + 2*u^79 + 5*u^77 + 12*u^75 + 19*u^73 + 31*u^71 + 49*u^69 + 
#            64*u^67 + 84*u^65 + 108*u^63 + 122*u^61 + 139*u^59 + 154*u^57 + 
#            154*u^55 + 156*u^53 + 154*u^51 + 137*u^49 + 126*u^47 + 112*u^45 + 
#            89*u^43 + 75*u^41 + 61*u^39 + 43*u^37 + 34*u^35 + 25*u^33 + 15*u^
#            31 + 11*u^29 + 7*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          2*u^65 + 7*u^63 + 14*u^61 + 27*u^59 + 43*u^57 + 58*u^55 + 76*u^53 + 
#            89*u^51 + 93*u^49 + 97*u^47 + 93*u^45 + 80*u^43 + 71*u^41 + 59*u^
#            39 + 43*u^37 + 34*u^35 + 25*u^33 + 15*u^31 + 11*u^29 + 7*u^27 + 
#            3*u^25 + 2*u^23 + u^21, u^57 + 3*u^55 + 7*u^53 + 15*u^51 + 23*u^
#            49 + 31*u^47 + 43*u^45 + 45*u^43 + 47*u^41 + 48*u^39 + 38*u^37 + 
#            32*u^35 + 25*u^33 + 15*u^31 + 11*u^29 + 7*u^27 + 3*u^25 + 2*u^
#            23 + u^21, u^56 + 3*u^54 + 5*u^52 + 8*u^50 + 9*u^48 + 10*u^46 + 
#            9*u^44 + 7*u^42 + 5*u^40 + 2*u^38 + u^36, 
#          2*u^54 + 5*u^52 + 8*u^50 + 16*u^48 + 20*u^46 + 22*u^44 + 28*u^42 + 
#            23*u^40 + 19*u^38 + 16*u^36 + 8*u^34 + 4*u^32 + 2*u^30, 
#          u^47 + 6*u^45 + 13*u^43 + 22*u^41 + 30*u^39 + 29*u^37 + 28*u^35 + 
#            23*u^33 + 15*u^31 + 11*u^29 + 7*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          u^46 + 5*u^44 + 13*u^42 + 15*u^40 + 16*u^38 + 15*u^36 + 8*u^34 + 
#            4*u^32 + 2*u^30, u^47 + 3*u^45 + 6*u^43 + 8*u^41 + 7*u^39 + 5*u^
#            37 + 2*u^35 + u^33, u^43 + 7*u^41 + 16*u^39 + 21*u^37 + 25*u^35 + 
#            22*u^33 + 15*u^31 + 11*u^29 + 7*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          2*u^41 + 6*u^39 + 11*u^37 + 16*u^35 + 19*u^33 + 18*u^31 + 15*u^29 + 
#            10*u^27 + 5*u^25 + 2*u^23 + u^21, 
#          u^39 + 3*u^37 + 6*u^35 + 12*u^33 + 10*u^31 + 9*u^29 + 7*u^27 + 3*u^
#            25 + 2*u^23 + u^21, 2*u^38 + 3*u^36 + 4*u^34 + 2*u^32 + u^30, 
#          u^37 + 2*u^35 + 7*u^33 + 5*u^31 + 3*u^29 + 2*u^27, 
#          3*u^33 + 4*u^31 + 6*u^29 + 6*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          u^34 + 3*u^32 + 4*u^30 + 2*u^28 + u^26, 
#          u^33 + 4*u^31 + 4*u^29 + 3*u^27 + 2*u^25, u^31, 0*u^0, 
#          u^33 + u^31 + u^29 + 2*u^27 + u^25 + u^23 + u^21, 
#          3*u^29 + 4*u^27 + 3*u^25 + 2*u^23 + u^21, 2*u^29 + 2*u^27 + 2*u^25, 
#          3*u^29 + 2*u^27 + u^25, u^29, 0*u^0, 
#          2*u^29 + 4*u^27 + 4*u^25 + 2*u^23 + u^21, 
#          u^27 + 2*u^25 + 2*u^23 + u^21, 0*u^0, u^26 + 2*u^24 + u^22, 
#          u^27 + u^25 + 2*u^23 + u^21, u^26 + u^24, u^24, 2*u^23 + u^21, 
#          u^24 + u^22, u^23 + u^21, 0*u^0, u^24 + u^22, u^23 + u^21, u^22, 
#          0*u^0, u^21, u^21, u^21, 0*u^0, u^21, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^104 + 2*u^102 + 3*u^100 + 7*u^98 + 10*u^96 + 14*u^94 + 23*u^92 + 
#            29*u^90 + 38*u^88 + 53*u^86 + 62*u^84 + 76*u^82 + 95*u^80 + 104*u^
#            78 + 121*u^76 + 139*u^74 + 144*u^72 + 160*u^70 + 171*u^68 + 169*u^
#            66 + 179*u^64 + 179*u^62 + 169*u^60 + 171*u^58 + 160*u^56 + 144*u^
#            54 + 139*u^52 + 121*u^50 + 104*u^48 + 95*u^46 + 76*u^44 + 62*u^
#            42 + 53*u^40 + 38*u^38 + 29*u^36 + 23*u^34 + 14*u^32 + 10*u^30 + 
#            7*u^28 + 3*u^26 + 2*u^24 + u^22, 
#          u^80 + 2*u^78 + 6*u^76 + 12*u^74 + 18*u^72 + 31*u^70 + 44*u^68 + 
#            56*u^66 + 75*u^64 + 89*u^62 + 99*u^60 + 114*u^58 + 118*u^56 + 
#            117*u^54 + 120*u^52 + 110*u^50 + 99*u^48 + 92*u^46 + 75*u^44 + 
#            62*u^42 + 53*u^40 + 38*u^38 + 29*u^36 + 23*u^34 + 14*u^32 + 10*u^
#            30 + 7*u^28 + 3*u^26 + 2*u^24 + u^22, 
#          u^66 + 4*u^64 + 8*u^62 + 16*u^60 + 28*u^58 + 38*u^56 + 51*u^54 + 
#            64*u^52 + 68*u^50 + 72*u^48 + 73*u^46 + 64*u^44 + 57*u^42 + 50*u^
#            40 + 37*u^38 + 29*u^36 + 23*u^34 + 14*u^32 + 10*u^30 + 7*u^28 + 
#            3*u^26 + 2*u^24 + u^22, u^56 + 3*u^54 + 7*u^52 + 14*u^50 + 20*u^
#            48 + 28*u^46 + 34*u^44 + 35*u^42 + 37*u^40 + 32*u^38 + 26*u^36 + 
#            22*u^34 + 14*u^32 + 10*u^30 + 7*u^28 + 3*u^26 + 2*u^24 + u^22, 
#          u^57 + 2*u^55 + 3*u^53 + 6*u^51 + 6*u^49 + 6*u^47 + 7*u^45 + 4*u^
#            43 + 3*u^41 + 2*u^39, 3*u^53 + 5*u^51 + 9*u^49 + 16*u^47 + 17*u^
#            45 + 21*u^43 + 22*u^41 + 17*u^39 + 15*u^37 + 10*u^35 + 5*u^33 + 
#            3*u^31 + u^29, 2*u^46 + 7*u^44 + 13*u^42 + 19*u^40 + 21*u^38 + 
#            21*u^36 + 19*u^34 + 13*u^32 + 10*u^30 + 7*u^28 + 3*u^26 + 2*u^
#            24 + u^22, u^47 + 3*u^45 + 9*u^43 + 13*u^41 + 13*u^39 + 14*u^37 + 
#            10*u^35 + 5*u^33 + 3*u^31 + u^29, 
#          2*u^44 + 4*u^42 + 3*u^40 + 3*u^38 + 2*u^36, 
#          3*u^42 + 9*u^40 + 13*u^38 + 17*u^36 + 18*u^34 + 13*u^32 + 10*u^30 + 
#            7*u^28 + 3*u^26 + 2*u^24 + u^22, 
#          u^42 + 3*u^40 + 6*u^38 + 10*u^36 + 13*u^34 + 14*u^32 + 13*u^30 + 
#            10*u^28 + 6*u^26 + 3*u^24 + u^22, 
#          2*u^38 + 4*u^36 + 7*u^34 + 8*u^32 + 7*u^30 + 6*u^28 + 3*u^26 + 2*u^
#            24 + u^22, u^33, 2*u^38 + 3*u^36 + 5*u^34 + 6*u^32 + 4*u^30 + 3*u^
#            28 + u^26, u^34 + 2*u^32 + 3*u^30 + 5*u^28 + 3*u^26 + 2*u^24 + u^
#            22, u^35 + 2*u^33 + 3*u^31 + 3*u^29 + u^27, 
#          u^32 + 3*u^30 + 3*u^28 + 3*u^26 + u^24, 0*u^0, u^28, 0*u^0, 
#          2*u^28 + 2*u^26 + 2*u^24 + u^22, 2*u^28 + 3*u^26 + u^24, 
#          2*u^28 + u^26, 0*u^0, u^28, 3*u^28 + 5*u^26 + 3*u^24 + u^22, 
#          u^28 + 2*u^26 + 2*u^24 + u^22, u^28 + u^26, u^27 + 2*u^25 + u^23, 
#          u^26 + 2*u^24 + u^22, 0*u^0, 2*u^25 + u^23, u^24 + u^22, 
#          u^25 + u^23, 0*u^0, 0*u^0, 0*u^0, u^22, 0*u^0, u^22, u^22, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^21, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^103 + u^101 + 3*u^99 + 6*u^97 + 7*u^95 + 13*u^93 + 19*u^91 + 22*u^
#            89 + 34*u^87 + 43*u^85 + 48*u^83 + 66*u^81 + 75*u^79 + 81*u^77 + 
#            102*u^75 + 107*u^73 + 112*u^71 + 131*u^69 + 128*u^67 + 130*u^65 + 
#            142*u^63 + 130*u^61 + 128*u^59 + 131*u^57 + 112*u^55 + 107*u^53 + 
#            102*u^51 + 81*u^49 + 75*u^47 + 66*u^45 + 48*u^43 + 43*u^41 + 34*u^
#            39 + 22*u^37 + 19*u^35 + 13*u^33 + 7*u^31 + 6*u^29 + 3*u^27 + u^
#            25 + u^23, u^81 + 2*u^79 + 3*u^77 + 9*u^75 + 12*u^73 + 19*u^71 + 
#            31*u^69 + 38*u^67 + 49*u^65 + 65*u^63 + 69*u^61 + 80*u^59 + 90*u^
#            57 + 86*u^55 + 89*u^53 + 89*u^51 + 75*u^49 + 72*u^47 + 64*u^45 + 
#            48*u^43 + 43*u^41 + 34*u^39 + 22*u^37 + 19*u^35 + 13*u^33 + 7*u^
#            31 + 6*u^29 + 3*u^27 + u^25 + u^23, 
#          2*u^65 + 5*u^63 + 8*u^61 + 17*u^59 + 25*u^57 + 31*u^55 + 43*u^53 + 
#            49*u^51 + 49*u^49 + 54*u^47 + 51*u^45 + 42*u^43 + 40*u^41 + 32*u^
#            39 + 22*u^37 + 19*u^35 + 13*u^33 + 7*u^31 + 6*u^29 + 3*u^27 + u^
#            25 + u^23, u^57 + 3*u^55 + 5*u^53 + 11*u^51 + 15*u^49 + 19*u^47 + 
#            26*u^45 + 25*u^43 + 26*u^41 + 26*u^39 + 19*u^37 + 17*u^35 + 13*u^
#            33 + 7*u^31 + 6*u^29 + 3*u^27 + u^25 + u^23, 
#          u^50 + u^46 + u^44 + u^40, u^58 + u^56 + 4*u^54 + 7*u^52 + 8*u^50 + 
#            15*u^48 + 16*u^46 + 16*u^44 + 20*u^42 + 15*u^40 + 12*u^38 + 11*u^
#            36 + 5*u^34 + 3*u^32 + 2*u^30, u^47 + 4*u^45 + 7*u^43 + 11*u^41 + 
#            14*u^39 + 13*u^37 + 14*u^35 + 11*u^33 + 7*u^31 + 6*u^29 + 3*u^
#            27 + u^25 + u^23, u^48 + 2*u^46 + 5*u^44 + 10*u^42 + 9*u^40 + 
#            10*u^38 + 10*u^36 + 5*u^34 + 3*u^32 + 2*u^30, u^43 + u^41 + u^37, 
#          u^43 + 4*u^41 + 7*u^39 + 8*u^37 + 12*u^35 + 10*u^33 + 7*u^31 + 6*u^
#            29 + 3*u^27 + u^25 + u^23, u^43 + 3*u^41 + 6*u^39 + 10*u^37 + 
#            13*u^35 + 14*u^33 + 13*u^31 + 10*u^29 + 6*u^27 + 3*u^25 + u^23, 
#          u^37 + 2*u^35 + 5*u^33 + 4*u^31 + 4*u^29 + 3*u^27 + u^25 + u^23, 
#          u^38 + u^34, u^37 + u^35 + 4*u^33 + 3*u^31 + 2*u^29 + 2*u^27, 
#          2*u^33 + 2*u^31 + 3*u^29 + 3*u^27 + u^25 + u^23, 0*u^0, 
#          u^35 + 2*u^33 + 5*u^31 + 4*u^29 + 3*u^27 + 2*u^25, 0*u^0, 
#          2*u^33 + u^31 + 2*u^29 + u^27, u^33 + u^31 + u^29 + 2*u^27 + u^
#            25 + u^23 + u^21, u^29 + u^27 + u^25 + u^23, 
#          u^31 + u^29 + 2*u^27 + 2*u^25, 0*u^0, 0*u^0, u^29 + u^27, 
#          u^27 + 2*u^25 + u^23, u^27 + 2*u^25 + u^23, 0*u^0, 
#          u^28 + 2*u^26 + 2*u^24 + u^22, u^23, 0*u^0, u^24, u^23, 
#          u^24 + u^22, u^23 + u^21, u^24 + u^22, 0*u^0, 0*u^0, u^22, 0*u^0, 
#          0*u^0, u^21, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^21, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^100 + u^98 + 3*u^96 + 5*u^94 + 7*u^92 + 11*u^90 + 17*u^88 + 20*u^
#            86 + 30*u^84 + 37*u^82 + 44*u^80 + 56*u^78 + 66*u^76 + 72*u^74 + 
#            88*u^72 + 93*u^70 + 100*u^68 + 111*u^66 + 113*u^64 + 114*u^62 + 
#            122*u^60 + 114*u^58 + 113*u^56 + 111*u^54 + 100*u^52 + 93*u^50 + 
#            88*u^48 + 72*u^46 + 66*u^44 + 56*u^42 + 44*u^40 + 37*u^38 + 30*u^
#            36 + 20*u^34 + 17*u^32 + 11*u^30 + 7*u^28 + 5*u^26 + 3*u^24 + u^
#            22 + u^20, u^78 + 2*u^76 + 4*u^74 + 9*u^72 + 13*u^70 + 20*u^68 + 
#            30*u^66 + 38*u^64 + 48*u^62 + 61*u^60 + 66*u^58 + 75*u^56 + 81*u^
#            54 + 80*u^52 + 80*u^50 + 79*u^48 + 68*u^46 + 64*u^44 + 55*u^42 + 
#            44*u^40 + 37*u^38 + 30*u^36 + 20*u^34 + 17*u^32 + 11*u^30 + 7*u^
#            28 + 5*u^26 + 3*u^24 + u^22 + u^20, 
#          u^64 + 3*u^62 + 7*u^60 + 11*u^58 + 20*u^56 + 27*u^54 + 35*u^52 + 
#            43*u^50 + 49*u^48 + 48*u^46 + 51*u^44 + 46*u^42 + 40*u^40 + 35*u^
#            38 + 29*u^36 + 20*u^34 + 17*u^32 + 11*u^30 + 7*u^28 + 5*u^26 + 
#            3*u^24 + u^22 + u^20, u^56 + 2*u^54 + 5*u^52 + 8*u^50 + 14*u^48 + 
#            18*u^46 + 23*u^44 + 27*u^42 + 27*u^40 + 26*u^38 + 25*u^36 + 18*u^
#            34 + 16*u^32 + 11*u^30 + 7*u^28 + 5*u^26 + 3*u^24 + u^22 + u^20, 
#          u^47 + u^43, u^55 + u^53 + 4*u^51 + 7*u^49 + 9*u^47 + 14*u^45 + 
#            16*u^43 + 15*u^41 + 17*u^39 + 13*u^37 + 10*u^35 + 8*u^33 + 4*u^
#            31 + 2*u^29 + u^27, u^46 + 3*u^44 + 7*u^42 + 11*u^40 + 14*u^38 + 
#            16*u^36 + 14*u^34 + 14*u^32 + 10*u^30 + 7*u^28 + 5*u^26 + 3*u^
#            24 + u^22 + u^20, 2*u^45 + 4*u^43 + 7*u^41 + 11*u^39 + 10*u^37 + 
#            9*u^35 + 8*u^33 + 4*u^31 + 2*u^29 + u^27, u^40, 
#          u^42 + 3*u^40 + 7*u^38 + 11*u^36 + 11*u^34 + 13*u^32 + 10*u^30 + 
#            7*u^28 + 5*u^26 + 3*u^24 + u^22 + u^20, 
#          u^42 + 3*u^40 + 6*u^38 + 10*u^36 + 13*u^34 + 15*u^32 + 14*u^30 + 
#            12*u^28 + 8*u^26 + 5*u^24 + 2*u^22 + u^20, 
#          2*u^36 + 3*u^34 + 5*u^32 + 6*u^30 + 5*u^28 + 4*u^26 + 3*u^24 + u^
#            22 + u^20, 0*u^0, u^36 + 3*u^34 + 3*u^32 + 5*u^30 + 3*u^28 + 2*u^
#            26 + u^24, u^34 + 2*u^32 + 4*u^30 + 4*u^28 + 4*u^26 + 3*u^24 + u^
#            22 + u^20, 0*u^0, u^34 + 2*u^32 + 3*u^30 + 5*u^28 + 3*u^26 + 2*u^
#            24 + u^22, 0*u^0, u^30 + u^26, u^34 + u^32 + 2*u^30 + 2*u^28 + 
#            2*u^26 + 2*u^24 + u^22 + u^20, u^28 + 2*u^26 + 2*u^24 + u^22 + u^
#            20, 2*u^28 + 2*u^26 + 2*u^24 + u^22, 0*u^0, 0*u^0, u^26, 
#          u^28 + 2*u^26 + 3*u^24 + 2*u^22 + u^20, 
#          u^28 + 2*u^26 + 3*u^24 + 2*u^22 + u^20, 0*u^0, 
#          u^27 + 2*u^25 + 2*u^23 + u^21, u^24 + u^22 + u^20, 0*u^0, 
#          u^23 + u^21, u^24 + u^22 + u^20, u^23 + u^21, u^24 + u^22 + u^20, 
#          u^23 + u^21, 0*u^0, u^20, 0*u^0, 0*u^0, u^20, u^20, 0*u^0, u^20, 
#          0*u^0, 0*u^0, 0*u^0, u^20, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^104 + 2*u^100 + 4*u^98 + 7*u^96 + 11*u^94 + 19*u^92 + 23*u^90 + 
#            37*u^88 + 48*u^86 + 63*u^84 + 78*u^82 + 101*u^80 + 114*u^78 + 
#            141*u^76 + 160*u^74 + 179*u^72 + 197*u^70 + 220*u^68 + 223*u^66 + 
#            243*u^64 + 246*u^62 + 246*u^60 + 246*u^58 + 243*u^56 + 223*u^54 + 
#            220*u^52 + 197*u^50 + 179*u^48 + 160*u^46 + 141*u^44 + 114*u^42 + 
#            101*u^40 + 78*u^38 + 63*u^36 + 48*u^34 + 37*u^32 + 23*u^30 + 19*u^
#            28 + 11*u^26 + 7*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          u^78 + 4*u^76 + 8*u^74 + 15*u^72 + 26*u^70 + 42*u^68 + 57*u^66 + 
#            81*u^64 + 102*u^62 + 123*u^60 + 144*u^58 + 162*u^56 + 166*u^54 + 
#            178*u^52 + 171*u^50 + 164*u^48 + 152*u^46 + 137*u^44 + 113*u^42 + 
#            101*u^40 + 78*u^38 + 63*u^36 + 48*u^34 + 37*u^32 + 23*u^30 + 19*u^
#            28 + 11*u^26 + 7*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          2*u^64 + 5*u^62 + 13*u^60 + 23*u^58 + 41*u^56 + 56*u^54 + 81*u^52 + 
#            92*u^50 + 107*u^48 + 110*u^46 + 111*u^44 + 98*u^42 + 93*u^40 + 
#            74*u^38 + 62*u^36 + 48*u^34 + 37*u^32 + 23*u^30 + 19*u^28 + 11*u^
#            26 + 7*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          u^56 + 2*u^54 + 7*u^52 + 15*u^50 + 23*u^48 + 36*u^46 + 49*u^44 + 
#            52*u^42 + 62*u^40 + 57*u^38 + 53*u^36 + 44*u^34 + 36*u^32 + 23*u^
#            30 + 19*u^28 + 11*u^26 + 7*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          u^55 + 2*u^53 + 5*u^51 + 6*u^49 + 10*u^47 + 8*u^45 + 10*u^43 + 6*u^
#            41 + 5*u^39 + 2*u^37 + u^35, u^53 + 4*u^51 + 8*u^49 + 18*u^47 + 
#            20*u^45 + 30*u^43 + 31*u^41 + 30*u^39 + 27*u^37 + 21*u^35 + 13*u^
#            33 + 8*u^31 + 4*u^29 + u^27, u^46 + 6*u^44 + 13*u^42 + 27*u^40 + 
#            33*u^38 + 38*u^36 + 36*u^34 + 32*u^32 + 22*u^30 + 19*u^28 + 11*u^
#            26 + 7*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          u^45 + 7*u^43 + 14*u^41 + 19*u^39 + 22*u^37 + 19*u^35 + 13*u^33 + 
#            8*u^31 + 4*u^29 + u^27, 2*u^44 + 4*u^42 + 8*u^40 + 6*u^38 + 5*u^
#            36 + 2*u^34 + u^32, u^42 + 8*u^40 + 17*u^38 + 27*u^36 + 31*u^34 + 
#            30*u^32 + 22*u^30 + 19*u^28 + 11*u^26 + 7*u^24 + 4*u^22 + 2*u^
#            20 + u^16, 3*u^40 + 7*u^38 + 15*u^36 + 21*u^34 + 27*u^32 + 26*u^
#            30 + 25*u^28 + 17*u^26 + 11*u^24 + 5*u^22 + 2*u^20 + u^16, 
#          u^38 + 4*u^36 + 9*u^34 + 15*u^32 + 13*u^30 + 15*u^28 + 10*u^26 + 
#            7*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          u^37 + 2*u^35 + 3*u^33 + 2*u^31 + u^29, 
#          u^38 + u^36 + 6*u^34 + 8*u^32 + 8*u^30 + 6*u^28 + 4*u^26 + u^24, 
#          u^34 + 4*u^32 + 6*u^30 + 10*u^28 + 9*u^26 + 7*u^24 + 4*u^22 + 2*u^
#            20 + u^16, u^33 + 3*u^31 + 4*u^29 + 3*u^27 + u^25, 
#          2*u^32 + 5*u^30 + 6*u^28 + 6*u^26 + 4*u^24 + u^22, u^32 + u^30, 
#          0*u^0, u^32 + 2*u^28 + 2*u^26 + 2*u^24 + 2*u^22 + 2*u^20 + u^16, 
#          4*u^28 + 5*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          3*u^28 + 4*u^26 + 4*u^24 + u^22, 2*u^28 + 3*u^26 + u^24, u^28, 
#          0*u^0, 3*u^28 + 7*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^16, 
#          u^28 + 2*u^26 + 5*u^24 + 5*u^22 + 2*u^20 + u^16, u^26, 
#          2*u^25 + 4*u^23 + 2*u^21, u^26 + 3*u^24 + 4*u^22 + 2*u^20 + u^16, 
#          u^25 + u^23, u^25 + 2*u^23 + u^21, u^24 + 4*u^22 + 2*u^20 + u^16, 
#          2*u^23 + 2*u^21, 2*u^22 + 2*u^20 + u^16, u^21, u^21, 
#          2*u^22 + 2*u^20 + u^16, u^21, u^22, u^22 + 2*u^20 + u^16, 
#          2*u^20 + u^16, u^20 + u^16, u^20 + u^16, u^20 + u^16, 0*u^0, 0*u^0, 
#          0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^103 + 2*u^101 + 4*u^99 + 9*u^97 + 14*u^95 + 22*u^93 + 35*u^91 + 
#            48*u^89 + 66*u^87 + 90*u^85 + 112*u^83 + 142*u^81 + 175*u^79 + 
#            203*u^77 + 240*u^75 + 273*u^73 + 300*u^71 + 334*u^69 + 356*u^67 + 
#            372*u^65 + 392*u^63 + 394*u^61 + 394*u^59 + 392*u^57 + 372*u^55 + 
#            356*u^53 + 334*u^51 + 300*u^49 + 273*u^47 + 240*u^45 + 203*u^43 + 
#            175*u^41 + 142*u^39 + 112*u^37 + 90*u^35 + 66*u^33 + 48*u^31 + 
#            35*u^29 + 22*u^27 + 14*u^25 + 9*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          u^79 + 3*u^77 + 9*u^75 + 18*u^73 + 32*u^71 + 54*u^69 + 79*u^67 + 
#            110*u^65 + 147*u^63 + 180*u^61 + 214*u^59 + 245*u^57 + 262*u^55 + 
#            277*u^53 + 280*u^51 + 268*u^49 + 255*u^47 + 231*u^45 + 200*u^43 + 
#            174*u^41 + 142*u^39 + 112*u^37 + 90*u^35 + 66*u^33 + 48*u^31 + 
#            35*u^29 + 22*u^27 + 14*u^25 + 9*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          u^65 + 5*u^63 + 13*u^61 + 29*u^59 + 51*u^57 + 77*u^55 + 110*u^53 + 
#            138*u^51 + 159*u^49 + 176*u^47 + 177*u^45 + 168*u^43 + 156*u^41 + 
#            133*u^39 + 109*u^37 + 89*u^35 + 66*u^33 + 48*u^31 + 35*u^29 + 
#            22*u^27 + 14*u^25 + 9*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          2*u^55 + 7*u^53 + 15*u^51 + 30*u^49 + 47*u^47 + 65*u^45 + 82*u^43 + 
#            92*u^41 + 95*u^39 + 89*u^37 + 79*u^35 + 63*u^33 + 47*u^31 + 35*u^
#            29 + 22*u^27 + 14*u^25 + 9*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          u^56 + 2*u^54 + 6*u^52 + 11*u^50 + 12*u^48 + 16*u^46 + 16*u^44 + 
#            12*u^42 + 11*u^40 + 6*u^38 + 2*u^36 + u^34, 
#          u^54 + 5*u^52 + 10*u^50 + 20*u^48 + 32*u^46 + 40*u^44 + 50*u^42 + 
#            51*u^40 + 46*u^38 + 40*u^36 + 27*u^34 + 17*u^32 + 9*u^30 + 3*u^
#            28 + u^26, 4*u^45 + 15*u^43 + 31*u^41 + 47*u^39 + 58*u^37 + 61*u^
#            35 + 54*u^33 + 44*u^31 + 34*u^29 + 22*u^27 + 14*u^25 + 9*u^23 + 
#            4*u^21 + 2*u^19 + u^17, u^46 + 6*u^44 + 17*u^42 + 27*u^40 + 34*u^
#            38 + 35*u^36 + 26*u^34 + 17*u^32 + 9*u^30 + 3*u^28 + u^26, 
#          u^45 + 6*u^43 + 10*u^41 + 11*u^39 + 11*u^37 + 6*u^35 + 2*u^33 + u^31
#            , 7*u^41 + 19*u^39 + 35*u^37 + 49*u^35 + 49*u^33 + 43*u^31 + 34*u^
#            29 + 22*u^27 + 14*u^25 + 9*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          2*u^41 + 7*u^39 + 17*u^37 + 29*u^35 + 39*u^33 + 45*u^31 + 42*u^29 + 
#            33*u^27 + 22*u^25 + 12*u^23 + 5*u^21 + 2*u^19 + u^17, 
#          4*u^37 + 10*u^35 + 18*u^33 + 24*u^31 + 24*u^29 + 19*u^27 + 13*u^
#            25 + 9*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          u^38 + 2*u^36 + 5*u^34 + 5*u^32 + 2*u^30 + u^28, 
#          2*u^37 + 5*u^35 + 11*u^33 + 14*u^31 + 12*u^29 + 8*u^27 + 3*u^25 + u^
#            23, 3*u^33 + 7*u^31 + 13*u^29 + 15*u^27 + 12*u^25 + 9*u^23 + 4*u^
#            21 + 2*u^19 + u^17, u^34 + 4*u^32 + 8*u^30 + 7*u^28 + 3*u^26 + u^
#            24, u^33 + 5*u^31 + 9*u^29 + 11*u^27 + 8*u^25 + 3*u^23 + u^21, 
#          u^33 + u^29, u^29 + u^27, u^31 + 2*u^29 + 2*u^27 + 3*u^25 + 3*u^
#            23 + 2*u^21 + 2*u^19 + u^17, 2*u^29 + 7*u^27 + 9*u^25 + 8*u^23 + 
#            4*u^21 + 2*u^19 + u^17, u^29 + 7*u^27 + 7*u^25 + 3*u^23 + u^21, 
#          2*u^29 + 6*u^27 + 3*u^25 + u^23, u^27, u^27, 
#          u^29 + 8*u^27 + 13*u^25 + 10*u^23 + 5*u^21 + 2*u^19 + u^17, 
#          2*u^27 + 6*u^25 + 8*u^23 + 5*u^21 + 2*u^19 + u^17, u^27 + u^25, 
#          3*u^26 + 6*u^24 + 5*u^22 + 2*u^20, 
#          3*u^25 + 6*u^23 + 4*u^21 + 2*u^19 + u^17, u^26 + 2*u^24 + u^22, 
#          3*u^24 + 2*u^22 + u^20, 5*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          2*u^24 + 4*u^22 + 2*u^20, u^23 + 2*u^21 + 2*u^19 + u^17, 
#          u^22 + u^20, u^22 + u^20, u^23 + 3*u^21 + 2*u^19 + u^17, 
#          2*u^22 + u^20, u^21, 3*u^21 + 2*u^19 + u^17, u^21 + 2*u^19 + u^17, 
#          u^21 + u^19 + u^17, u^19 + u^17, u^19 + u^17, u^20 + u^16, 0*u^0, 
#          0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^102 + 2*u^100 + 2*u^98 + 6*u^96 + 8*u^94 + 12*u^92 + 19*u^90 + 26*u^
#            88 + 32*u^86 + 47*u^84 + 55*u^82 + 68*u^80 + 85*u^78 + 98*u^76 + 
#            109*u^74 + 130*u^72 + 137*u^70 + 150*u^68 + 163*u^66 + 167*u^64 + 
#            169*u^62 + 178*u^60 + 169*u^58 + 167*u^56 + 163*u^54 + 150*u^52 + 
#            137*u^50 + 130*u^48 + 109*u^46 + 98*u^44 + 85*u^42 + 68*u^40 + 
#            55*u^38 + 47*u^36 + 32*u^34 + 26*u^32 + 19*u^30 + 12*u^28 + 8*u^
#            26 + 6*u^24 + 2*u^22 + 2*u^20 + u^18, 
#          u^78 + 2*u^76 + 5*u^74 + 11*u^72 + 17*u^70 + 28*u^68 + 42*u^66 + 
#            55*u^64 + 70*u^62 + 89*u^60 + 99*u^58 + 112*u^56 + 121*u^54 + 
#            122*u^52 + 120*u^50 + 119*u^48 + 104*u^46 + 96*u^44 + 84*u^42 + 
#            68*u^40 + 55*u^38 + 47*u^36 + 32*u^34 + 26*u^32 + 19*u^30 + 12*u^
#            28 + 8*u^26 + 6*u^24 + 2*u^22 + 2*u^20 + u^18, 
#          u^64 + 3*u^62 + 9*u^60 + 16*u^58 + 29*u^56 + 41*u^54 + 55*u^52 + 
#            66*u^50 + 77*u^48 + 76*u^46 + 79*u^44 + 73*u^42 + 63*u^40 + 53*u^
#            38 + 46*u^36 + 32*u^34 + 26*u^32 + 19*u^30 + 12*u^28 + 8*u^26 + 
#            6*u^24 + 2*u^22 + 2*u^20 + u^18, 2*u^54 + 4*u^52 + 8*u^50 + 17*u^
#            48 + 23*u^46 + 31*u^44 + 40*u^42 + 41*u^40 + 40*u^38 + 40*u^36 + 
#            30*u^34 + 25*u^32 + 19*u^30 + 12*u^28 + 8*u^26 + 6*u^24 + 2*u^
#            22 + 2*u^20 + u^18, u^55 + 2*u^53 + 5*u^51 + 6*u^49 + 9*u^47 + 
#            8*u^45 + 9*u^43 + 6*u^41 + 5*u^39 + 2*u^37 + u^35, 
#          3*u^51 + 5*u^49 + 9*u^47 + 15*u^45 + 19*u^43 + 20*u^41 + 23*u^39 + 
#            18*u^37 + 14*u^35 + 10*u^33 + 5*u^31 + 2*u^29 + u^27, 
#          3*u^44 + 10*u^42 + 17*u^40 + 24*u^38 + 29*u^36 + 25*u^34 + 23*u^
#            32 + 18*u^30 + 12*u^28 + 8*u^26 + 6*u^24 + 2*u^22 + 2*u^20 + u^18,
#          u^45 + 4*u^43 + 9*u^41 + 15*u^39 + 15*u^37 + 13*u^35 + 10*u^33 + 
#            5*u^31 + 2*u^29 + u^27, 2*u^44 + 4*u^42 + 7*u^40 + 6*u^38 + 5*u^
#            36 + 2*u^34 + u^32, u^42 + 5*u^40 + 13*u^38 + 21*u^36 + 22*u^34 + 
#            22*u^32 + 18*u^30 + 12*u^28 + 8*u^26 + 6*u^24 + 2*u^22 + 2*u^
#            20 + u^18, u^40 + 4*u^38 + 9*u^36 + 13*u^34 + 18*u^32 + 18*u^30 + 
#            16*u^28 + 12*u^26 + 8*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          4*u^36 + 6*u^34 + 10*u^32 + 12*u^30 + 10*u^28 + 7*u^26 + 6*u^24 + 
#            2*u^22 + 2*u^20 + u^18, u^37 + 2*u^35 + 3*u^33 + 2*u^31 + u^29, 
#          2*u^36 + 4*u^34 + 5*u^32 + 7*u^30 + 4*u^28 + 2*u^26 + u^24, 
#          2*u^32 + 5*u^30 + 6*u^28 + 6*u^26 + 6*u^24 + 2*u^22 + 2*u^20 + u^18,
#          u^33 + 3*u^31 + 4*u^29 + 3*u^27 + u^25, 
#          u^32 + 2*u^30 + 4*u^28 + 4*u^26 + 2*u^24 + u^22, 
#          u^34 + u^32 + u^30 + u^28, 0*u^0, u^30 + u^28 + 2*u^24 + u^20 + u^18
#            , 2*u^28 + 4*u^26 + 5*u^24 + 2*u^22 + 2*u^20 + u^18, 
#          2*u^28 + 3*u^26 + 2*u^24 + u^22, 2*u^28 + 3*u^26 + u^24, 
#          u^28 + u^26, 0*u^0, 2*u^28 + 6*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^
#            18, u^26 + 3*u^24 + 3*u^22 + 2*u^20 + u^18, u^26, 
#          u^25 + 2*u^23 + 2*u^21, u^26 + 3*u^24 + 2*u^22 + 2*u^20 + u^18, 
#          u^25 + u^23, u^25 + u^23 + u^21, u^24 + 2*u^22 + 2*u^20 + u^18, 
#          u^23 + 2*u^21, u^20 + u^18, u^21, u^21, u^22 + 2*u^20 + u^18, u^21, 
#          u^22, u^22 + 2*u^20 + u^18, u^20 + u^18, u^20 + u^18, u^18, u^18, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^102 + 2*u^100 + 4*u^98 + 6*u^96 + 12*u^94 + 16*u^92 + 25*u^90 + 
#            34*u^88 + 46*u^86 + 59*u^84 + 77*u^82 + 90*u^80 + 113*u^78 + 
#            130*u^76 + 149*u^74 + 168*u^72 + 187*u^70 + 198*u^68 + 217*u^66 + 
#            221*u^64 + 229*u^62 + 232*u^60 + 229*u^58 + 221*u^56 + 217*u^54 + 
#            198*u^52 + 187*u^50 + 168*u^48 + 149*u^46 + 130*u^44 + 113*u^42 + 
#            90*u^40 + 77*u^38 + 59*u^36 + 46*u^34 + 34*u^32 + 25*u^30 + 16*u^
#            28 + 12*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^78 + 3*u^76 + 7*u^74 + 14*u^72 + 24*u^70 + 37*u^68 + 56*u^66 + 
#            73*u^64 + 95*u^62 + 116*u^60 + 134*u^58 + 148*u^56 + 161*u^54 + 
#            161*u^52 + 163*u^50 + 154*u^48 + 142*u^46 + 127*u^44 + 112*u^42 + 
#            90*u^40 + 77*u^38 + 59*u^36 + 46*u^34 + 34*u^32 + 25*u^30 + 16*u^
#            28 + 12*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^64 + 5*u^62 + 11*u^60 + 23*u^58 + 37*u^56 + 56*u^54 + 71*u^52 + 
#            91*u^50 + 98*u^48 + 105*u^46 + 103*u^44 + 98*u^42 + 83*u^40 + 
#            74*u^38 + 58*u^36 + 46*u^34 + 34*u^32 + 25*u^30 + 16*u^28 + 12*u^
#            26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          2*u^54 + 5*u^52 + 12*u^50 + 20*u^48 + 32*u^46 + 41*u^44 + 52*u^42 + 
#            54*u^40 + 56*u^38 + 50*u^36 + 43*u^34 + 33*u^32 + 25*u^30 + 16*u^
#            28 + 12*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^55 + 4*u^53 + 5*u^51 + 10*u^49 + 11*u^47 + 12*u^45 + 11*u^43 + 
#            10*u^41 + 5*u^39 + 4*u^37 + u^35, 
#          u^53 + 3*u^51 + 8*u^49 + 12*u^47 + 20*u^45 + 25*u^43 + 28*u^41 + 
#            29*u^39 + 25*u^37 + 19*u^35 + 13*u^33 + 7*u^31 + 3*u^29 + u^27, 
#          u^46 + 4*u^44 + 13*u^42 + 23*u^40 + 33*u^38 + 36*u^36 + 36*u^34 + 
#            30*u^32 + 24*u^30 + 16*u^28 + 12*u^26 + 6*u^24 + 4*u^22 + 2*u^
#            20 + u^18, u^45 + 5*u^43 + 12*u^41 + 19*u^39 + 20*u^37 + 18*u^
#            35 + 13*u^33 + 7*u^31 + 3*u^29 + u^27, 
#          u^46 + 2*u^44 + 6*u^42 + 9*u^40 + 9*u^38 + 5*u^36 + 4*u^34 + u^32, 
#          u^42 + 6*u^40 + 18*u^38 + 26*u^36 + 31*u^34 + 29*u^32 + 24*u^30 + 
#            16*u^28 + 12*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^40 + 6*u^38 + 11*u^36 + 19*u^34 + 23*u^32 + 25*u^30 + 21*u^28 + 
#            17*u^26 + 9*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^38 + 4*u^36 + 10*u^34 + 13*u^32 + 16*u^30 + 13*u^28 + 11*u^26 + 
#            6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^37 + 3*u^35 + 3*u^33 + 3*u^31 + u^29, 
#          2*u^36 + 5*u^34 + 8*u^32 + 8*u^30 + 6*u^28 + 3*u^26 + u^24, 
#          2*u^32 + 6*u^30 + 7*u^28 + 9*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18,
#          2*u^33 + 5*u^31 + 6*u^29 + 5*u^27 + 2*u^25, 
#          u^32 + 3*u^30 + 5*u^28 + 5*u^26 + 3*u^24 + u^22, u^32 + u^30, 
#          0*u^0, u^30 + u^28 + 2*u^26 + u^24 + 2*u^22 + u^20 + u^18, 
#          2*u^28 + 6*u^26 + 5*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          2*u^28 + 4*u^26 + 3*u^24 + u^22, 4*u^28 + 4*u^26 + 2*u^24, u^28, 
#          0*u^0, 3*u^28 + 8*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          2*u^26 + 4*u^24 + 4*u^22 + 2*u^20 + u^18, u^28 + u^26 + u^24, 
#          2*u^25 + 3*u^23 + 2*u^21, u^26 + 3*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          2*u^25 + u^23, u^25 + 2*u^23 + u^21, u^24 + 3*u^22 + 2*u^20 + u^18, 
#          2*u^23 + 2*u^21, u^22 + u^20 + u^18, 0*u^0, u^23 + u^21 + u^19, 
#          2*u^22 + 2*u^20 + u^18, u^21, u^22 + u^20, u^22 + 2*u^20 + u^18, 
#          u^20 + u^18, u^20 + u^18, 0*u^0, u^20 + u^18 + u^16, u^19 + u^17, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^102 + 2*u^100 + 5*u^98 + 8*u^96 + 12*u^94 + 23*u^92 + 30*u^90 + 
#            42*u^88 + 63*u^86 + 76*u^84 + 98*u^82 + 129*u^80 + 144*u^78 + 
#            175*u^76 + 208*u^74 + 222*u^72 + 254*u^70 + 280*u^68 + 283*u^66 + 
#            308*u^64 + 318*u^62 + 308*u^60 + 318*u^58 + 308*u^56 + 283*u^54 + 
#            280*u^52 + 254*u^50 + 222*u^48 + 208*u^46 + 175*u^44 + 144*u^42 + 
#            129*u^40 + 98*u^38 + 76*u^36 + 63*u^34 + 42*u^32 + 30*u^30 + 23*u^
#            28 + 12*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^80 + u^78 + 5*u^76 + 12*u^74 + 19*u^72 + 36*u^70 + 55*u^68 + 73*u^
#            66 + 105*u^64 + 132*u^62 + 154*u^60 + 186*u^58 + 203*u^56 + 210*u^
#            54 + 225*u^52 + 218*u^50 + 203*u^48 + 196*u^46 + 170*u^44 + 143*u^
#            42 + 128*u^40 + 98*u^38 + 76*u^36 + 63*u^34 + 42*u^32 + 30*u^30 + 
#            23*u^28 + 12*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          3*u^64 + 7*u^62 + 16*u^60 + 33*u^58 + 50*u^56 + 72*u^54 + 100*u^
#            52 + 116*u^50 + 130*u^48 + 141*u^46 + 134*u^44 + 124*u^42 + 116*u^
#            40 + 93*u^38 + 75*u^36 + 62*u^34 + 42*u^32 + 30*u^30 + 23*u^28 + 
#            12*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^56 + 4*u^54 + 10*u^52 + 21*u^50 + 32*u^48 + 47*u^46 + 62*u^44 + 
#            67*u^42 + 76*u^40 + 73*u^38 + 62*u^36 + 57*u^34 + 41*u^32 + 29*u^
#            30 + 23*u^28 + 12*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^55 + u^53 + 4*u^51 + 5*u^49 + 6*u^47 + 8*u^45 + 6*u^43 + 5*u^41 + 
#            4*u^39 + u^37 + u^35, u^55 + 5*u^53 + 8*u^51 + 15*u^49 + 27*u^
#            47 + 31*u^45 + 40*u^43 + 46*u^41 + 38*u^39 + 37*u^37 + 29*u^35 + 
#            16*u^33 + 12*u^31 + 5*u^29 + u^27 + u^25, 
#          2*u^46 + 8*u^44 + 17*u^42 + 31*u^40 + 40*u^38 + 43*u^36 + 45*u^34 + 
#            36*u^32 + 28*u^30 + 22*u^28 + 12*u^26 + 8*u^24 + 5*u^22 + 2*u^
#            20 + u^18, u^47 + 3*u^45 + 11*u^43 + 20*u^41 + 24*u^39 + 30*u^
#            37 + 26*u^35 + 16*u^33 + 12*u^31 + 5*u^29 + u^27 + u^25, 
#          u^44 + 4*u^42 + 5*u^40 + 5*u^38 + 4*u^36 + u^34 + u^32, 
#          2*u^42 + 10*u^40 + 19*u^38 + 30*u^36 + 38*u^34 + 33*u^32 + 28*u^
#            30 + 22*u^28 + 12*u^26 + 8*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^42 + 5*u^40 + 12*u^38 + 22*u^36 + 32*u^34 + 38*u^32 + 38*u^30 + 
#            32*u^28 + 22*u^26 + 13*u^24 + 6*u^22 + 3*u^20 + u^18, 
#          u^38 + 4*u^36 + 10*u^34 + 16*u^32 + 15*u^30 + 17*u^28 + 11*u^26 + 
#            7*u^24 + 5*u^22 + 2*u^20 + u^18, u^37 + u^35 + 3*u^33 + u^31 + u^
#            29, u^38 + 2*u^36 + 6*u^34 + 12*u^32 + 9*u^30 + 9*u^28 + 5*u^
#            26 + u^24 + u^22, u^34 + 5*u^32 + 7*u^30 + 13*u^28 + 10*u^26 + 
#            7*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^33 + 3*u^31 + 3*u^29 + u^27 + u^25, 
#          u^34 + 4*u^32 + 10*u^30 + 10*u^28 + 10*u^26 + 5*u^24 + u^22 + u^20, 
#          0*u^0, u^32 + u^30 + 3*u^28 + u^26, u^32 + u^30 + 3*u^28 + 4*u^26 + 
#            3*u^24 + 4*u^22 + 3*u^20 + u^18 + u^16, 
#          4*u^28 + 5*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^30 + 4*u^28 + 7*u^26 + 5*u^24 + u^22 + u^20, 2*u^28 + u^26 + u^24,
#          0*u^0, 2*u^28 + u^26, 3*u^28 + 7*u^26 + 8*u^24 + 5*u^22 + 3*u^
#            20 + u^18, u^28 + 4*u^26 + 6*u^24 + 5*u^22 + 3*u^20 + u^18, 
#          0*u^0, 2*u^27 + 5*u^25 + 5*u^23 + 3*u^21 + u^19, 
#          u^26 + 2*u^24 + 4*u^22 + 2*u^20 + u^18, u^23, 
#          2*u^25 + 2*u^23 + u^21 + u^19, u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^25 + 2*u^23 + 3*u^21 + u^19, 3*u^22 + 3*u^20 + u^18 + u^16, 
#          u^23 + 2*u^21 + u^19, u^21, u^22 + 2*u^20 + u^18, 2*u^21, 0*u^0, 
#          2*u^20 + u^18, 3*u^20 + u^18 + u^16, u^20, u^20 + u^18 + u^16, 
#          u^20 + u^16, 0*u^0, u^20 + u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^102 + 2*u^100 + 3*u^98 + 6*u^96 + 12*u^94 + 15*u^92 + 26*u^90 + 
#            36*u^88 + 45*u^86 + 64*u^84 + 82*u^82 + 93*u^80 + 124*u^78 + 
#            141*u^76 + 156*u^74 + 188*u^72 + 202*u^70 + 212*u^68 + 243*u^66 + 
#            240*u^64 + 246*u^62 + 262*u^60 + 246*u^58 + 240*u^56 + 243*u^54 + 
#            212*u^52 + 202*u^50 + 188*u^48 + 156*u^46 + 141*u^44 + 124*u^42 + 
#            93*u^40 + 82*u^38 + 64*u^36 + 45*u^34 + 36*u^32 + 26*u^30 + 15*u^
#            28 + 12*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          2*u^78 + 4*u^76 + 8*u^74 + 18*u^72 + 27*u^70 + 41*u^68 + 65*u^66 + 
#            80*u^64 + 103*u^62 + 131*u^60 + 143*u^58 + 160*u^56 + 178*u^54 + 
#            171*u^52 + 175*u^50 + 170*u^48 + 148*u^46 + 137*u^44 + 122*u^42 + 
#            93*u^40 + 82*u^38 + 64*u^36 + 45*u^34 + 36*u^32 + 26*u^30 + 15*u^
#            28 + 12*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          u^66 + u^64 + 6*u^62 + 14*u^60 + 24*u^58 + 41*u^56 + 61*u^54 + 74*u^
#            52 + 96*u^50 + 106*u^48 + 107*u^46 + 110*u^44 + 104*u^42 + 85*u^
#            40 + 78*u^38 + 62*u^36 + 45*u^34 + 36*u^32 + 26*u^30 + 15*u^28 + 
#            12*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          u^56 + 3*u^54 + 9*u^52 + 15*u^50 + 28*u^48 + 38*u^46 + 47*u^44 + 
#            60*u^42 + 57*u^40 + 58*u^38 + 54*u^36 + 41*u^34 + 34*u^32 + 26*u^
#            30 + 15*u^28 + 12*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          2*u^53 + 2*u^51 + 4*u^49 + 5*u^47 + 4*u^45 + 5*u^43 + 4*u^41 + 2*u^
#            39 + 2*u^37, u^55 + 2*u^53 + 7*u^51 + 14*u^49 + 16*u^47 + 30*u^
#            45 + 31*u^43 + 32*u^41 + 36*u^39 + 27*u^37 + 20*u^35 + 17*u^33 + 
#            7*u^31 + 4*u^29 + 2*u^27, 2*u^46 + 6*u^44 + 16*u^42 + 23*u^40 + 
#            33*u^38 + 37*u^36 + 33*u^34 + 30*u^32 + 24*u^30 + 15*u^28 + 12*u^
#            26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          4*u^45 + 7*u^43 + 14*u^41 + 23*u^39 + 21*u^37 + 19*u^35 + 16*u^33 + 
#            7*u^31 + 4*u^29 + 2*u^27, u^44 + 2*u^42 + 4*u^40 + 4*u^38 + 2*u^
#            36 + 2*u^34, 2*u^42 + 6*u^40 + 17*u^38 + 25*u^36 + 27*u^34 + 29*u^
#            32 + 23*u^30 + 15*u^28 + 12*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18
#            , u^42 + 4*u^40 + 10*u^38 + 18*u^36 + 26*u^34 + 30*u^32 + 30*u^
#            30 + 25*u^28 + 18*u^26 + 10*u^24 + 5*u^22 + 2*u^20 + u^18, 
#          u^38 + 4*u^36 + 8*u^34 + 11*u^32 + 16*u^30 + 11*u^28 + 10*u^26 + 
#            6*u^24 + 3*u^22 + 2*u^20 + u^18, 2*u^35 + u^33 + 2*u^31, 
#          3*u^36 + 5*u^34 + 6*u^32 + 10*u^30 + 6*u^28 + 3*u^26 + 2*u^24, 
#          u^34 + 3*u^32 + 9*u^30 + 8*u^28 + 9*u^26 + 6*u^24 + 3*u^22 + 2*u^
#            20 + u^18, u^33 + u^31 + 2*u^29 + 2*u^27, 
#          u^34 + 4*u^32 + 6*u^30 + 10*u^28 + 6*u^26 + 4*u^24 + 2*u^22, 0*u^0, 
#          2*u^30 + u^26, 2*u^30 + 2*u^28 + 2*u^26 + 3*u^24 + 2*u^22 + 2*u^
#            20 + 2*u^18, u^30 + 2*u^28 + 6*u^26 + 4*u^24 + 3*u^22 + 2*u^
#            20 + u^18, 4*u^28 + 5*u^26 + 3*u^24 + 2*u^22, u^28 + 2*u^26, 
#          0*u^0, u^26, 2*u^28 + 6*u^26 + 6*u^24 + 4*u^22 + 2*u^20 + u^18, 
#          u^28 + 4*u^26 + 5*u^24 + 4*u^22 + 2*u^20 + u^18, 0*u^0, 
#          2*u^27 + 4*u^25 + 4*u^23 + 2*u^21 + u^19, 
#          3*u^24 + 2*u^22 + 2*u^20 + u^18, u^25, 2*u^23 + u^21, 
#          2*u^24 + 2*u^22 + 2*u^20 + u^18, 3*u^23 + u^21 + u^19, 
#          u^24 + u^22 + 2*u^20 + 2*u^18, u^23 + u^21 + u^19, 0*u^0, 
#          u^22 + u^20 + u^18, u^23 + u^19, 0*u^0, u^20 + u^18, u^20 + 2*u^18, 
#          u^18, u^20 + u^18 + u^16, u^18, 0*u^0, u^18, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ 2*u^101 + 3*u^99 + 5*u^97 + 12*u^95 + 16*u^93 + 24*u^91 + 40*u^89 + 
#            49*u^87 + 66*u^85 + 93*u^83 + 106*u^81 + 133*u^79 + 167*u^77 + 
#            180*u^75 + 213*u^73 + 245*u^71 + 251*u^69 + 282*u^67 + 302*u^65 + 
#            295*u^63 + 316*u^61 + 316*u^59 + 295*u^57 + 302*u^55 + 282*u^53 + 
#            251*u^51 + 245*u^49 + 213*u^47 + 180*u^45 + 167*u^43 + 133*u^41 + 
#            106*u^39 + 93*u^37 + 66*u^35 + 49*u^33 + 40*u^31 + 24*u^29 + 16*u^
#            27 + 12*u^25 + 5*u^23 + 3*u^21 + 2*u^19, 
#          u^79 + 4*u^77 + 6*u^75 + 16*u^73 + 28*u^71 + 40*u^69 + 65*u^67 + 
#            90*u^65 + 110*u^63 + 146*u^61 + 170*u^59 + 185*u^57 + 212*u^55 + 
#            217*u^53 + 211*u^51 + 217*u^49 + 197*u^47 + 174*u^45 + 163*u^43 + 
#            132*u^41 + 106*u^39 + 93*u^37 + 66*u^35 + 49*u^33 + 40*u^31 + 
#            24*u^29 + 16*u^27 + 12*u^25 + 5*u^23 + 3*u^21 + 2*u^19, 
#          u^65 + 4*u^63 + 12*u^61 + 22*u^59 + 39*u^57 + 64*u^55 + 83*u^53 + 
#            105*u^51 + 128*u^49 + 132*u^47 + 134*u^45 + 135*u^43 + 116*u^41 + 
#            100*u^39 + 89*u^37 + 65*u^35 + 49*u^33 + 40*u^31 + 24*u^29 + 16*u^
#            27 + 12*u^25 + 5*u^23 + 3*u^21 + 2*u^19, 
#          2*u^55 + 7*u^53 + 13*u^51 + 25*u^49 + 41*u^47 + 50*u^45 + 66*u^43 + 
#            73*u^41 + 69*u^39 + 71*u^37 + 59*u^35 + 45*u^33 + 39*u^31 + 24*u^
#            29 + 16*u^27 + 12*u^25 + 5*u^23 + 3*u^21 + 2*u^19, 
#          2*u^54 + 3*u^52 + 4*u^50 + 8*u^48 + 7*u^46 + 7*u^44 + 8*u^42 + 4*u^
#            40 + 3*u^38 + 2*u^36, u^56 + 2*u^54 + 5*u^52 + 14*u^50 + 18*u^
#            48 + 28*u^46 + 40*u^44 + 38*u^42 + 43*u^40 + 41*u^38 + 28*u^36 + 
#            24*u^34 + 15*u^32 + 6*u^30 + 4*u^28 + u^26, 
#          u^47 + 3*u^45 + 12*u^43 + 24*u^41 + 34*u^39 + 44*u^37 + 43*u^35 + 
#            39*u^33 + 35*u^31 + 23*u^29 + 16*u^27 + 12*u^25 + 5*u^23 + 3*u^
#            21 + 2*u^19, 2*u^46 + 7*u^44 + 13*u^42 + 25*u^40 + 29*u^38 + 24*u^
#            36 + 23*u^34 + 15*u^32 + 6*u^30 + 4*u^28 + u^26, 
#          u^45 + 2*u^43 + 5*u^41 + 7*u^39 + 4*u^37 + 3*u^35 + 2*u^33, 
#          u^43 + 4*u^41 + 14*u^39 + 28*u^37 + 32*u^35 + 35*u^33 + 34*u^31 + 
#            23*u^29 + 16*u^27 + 12*u^25 + 5*u^23 + 3*u^21 + 2*u^19, 
#          2*u^41 + 8*u^39 + 17*u^37 + 26*u^35 + 35*u^33 + 38*u^31 + 34*u^29 + 
#            27*u^27 + 18*u^25 + 9*u^23 + 4*u^21 + 2*u^19, 
#          2*u^37 + 8*u^35 + 11*u^33 + 17*u^31 + 17*u^29 + 12*u^27 + 11*u^25 + 
#            5*u^23 + 3*u^21 + 2*u^19, 2*u^36 + 2*u^34 + 2*u^32 + 2*u^30, 
#          u^37 + 6*u^35 + 7*u^33 + 11*u^31 + 11*u^29 + 5*u^27 + 4*u^25 + u^23,
#          u^33 + 7*u^31 + 10*u^29 + 9*u^27 + 11*u^25 + 5*u^23 + 3*u^21 + 2*u^
#            19, u^34 + 3*u^32 + 3*u^30 + 3*u^28 + 2*u^26, 
#          3*u^33 + 6*u^31 + 10*u^29 + 11*u^27 + 6*u^25 + 4*u^23 + u^21, 
#          0*u^0, 2*u^31 + 2*u^29 + u^27 + u^25, 
#          u^31 + 3*u^29 + 2*u^27 + 4*u^25 + 4*u^23 + 2*u^21 + 3*u^19 + u^17, 
#          u^29 + 4*u^27 + 7*u^25 + 4*u^23 + 3*u^21 + 2*u^19, 
#          3*u^29 + 6*u^27 + 5*u^25 + 4*u^23 + u^21, u^29 + 2*u^27 + 2*u^25, 
#          0*u^0, u^29 + u^27 + u^25, u^29 + 5*u^27 + 8*u^25 + 7*u^23 + 4*u^
#            21 + 2*u^19, 2*u^27 + 5*u^25 + 6*u^23 + 4*u^21 + 2*u^19, 0*u^0, 
#          u^28 + 3*u^26 + 5*u^24 + 5*u^22 + 2*u^20, 
#          2*u^25 + 3*u^23 + 3*u^21 + 2*u^19, u^24, 
#          u^26 + u^24 + 3*u^22 + u^20, u^25 + 2*u^23 + 3*u^21 + 2*u^19, 
#          u^24 + 4*u^22 + 2*u^20, u^23 + 2*u^21 + 3*u^19 + u^17, 
#          2*u^22 + 2*u^20, u^22 + u^20, 2*u^21 + 2*u^19, u^22 + u^20, 0*u^0, 
#          u^21 + 2*u^19, u^21 + 3*u^19 + u^17, u^19, u^19 + u^17, 
#          u^19 + u^17, u^18, u^19 + u^17, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^101 + 2*u^99 + 2*u^97 + 5*u^95 + 7*u^93 + 9*u^91 + 16*u^89 + 19*u^
#            87 + 24*u^85 + 35*u^83 + 39*u^81 + 48*u^79 + 61*u^77 + 64*u^75 + 
#            76*u^73 + 88*u^71 + 89*u^69 + 100*u^67 + 107*u^65 + 104*u^63 + 
#            112*u^61 + 112*u^59 + 104*u^57 + 107*u^55 + 100*u^53 + 89*u^51 + 
#            88*u^49 + 76*u^47 + 64*u^45 + 61*u^43 + 48*u^41 + 39*u^39 + 35*u^
#            37 + 24*u^35 + 19*u^33 + 16*u^31 + 9*u^29 + 7*u^27 + 5*u^25 + 2*u^
#            23 + 2*u^21 + u^19, u^77 + u^75 + 5*u^73 + 9*u^71 + 13*u^69 + 
#            22*u^67 + 31*u^65 + 38*u^63 + 52*u^61 + 60*u^59 + 66*u^57 + 76*u^
#            55 + 78*u^53 + 76*u^51 + 79*u^49 + 71*u^47 + 63*u^45 + 60*u^43 + 
#            48*u^41 + 39*u^39 + 35*u^37 + 24*u^35 + 19*u^33 + 16*u^31 + 9*u^
#            29 + 7*u^27 + 5*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          u^63 + 4*u^61 + 7*u^59 + 14*u^57 + 23*u^55 + 30*u^53 + 39*u^51 + 
#            48*u^49 + 49*u^47 + 50*u^45 + 51*u^43 + 43*u^41 + 38*u^39 + 34*u^
#            37 + 24*u^35 + 19*u^33 + 16*u^31 + 9*u^29 + 7*u^27 + 5*u^25 + 2*u^
#            23 + 2*u^21 + u^19, u^53 + 3*u^51 + 6*u^49 + 12*u^47 + 16*u^45 + 
#            22*u^43 + 26*u^41 + 26*u^39 + 27*u^37 + 23*u^35 + 18*u^33 + 16*u^
#            31 + 9*u^29 + 7*u^27 + 5*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          2*u^54 + 3*u^52 + 4*u^50 + 8*u^48 + 7*u^46 + 7*u^44 + 8*u^42 + 4*u^
#            40 + 3*u^38 + 2*u^36, 2*u^50 + 3*u^48 + 6*u^46 + 11*u^44 + 11*u^
#            42 + 13*u^40 + 14*u^38 + 9*u^36 + 8*u^34 + 5*u^32 + u^30 + u^28, 
#          3*u^43 + 9*u^41 + 14*u^39 + 18*u^37 + 18*u^35 + 17*u^33 + 15*u^31 + 
#            9*u^29 + 7*u^27 + 5*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          u^44 + 3*u^42 + 8*u^40 + 10*u^38 + 8*u^36 + 8*u^34 + 5*u^32 + u^
#            30 + u^28, u^45 + 2*u^43 + 5*u^41 + 7*u^39 + 4*u^37 + 3*u^35 + 
#            2*u^33, u^41 + 6*u^39 + 13*u^37 + 14*u^35 + 16*u^33 + 15*u^31 + 
#            9*u^29 + 7*u^27 + 5*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          u^39 + 3*u^37 + 5*u^35 + 9*u^33 + 11*u^31 + 10*u^29 + 9*u^27 + 6*u^
#            25 + 3*u^23 + 2*u^21 + u^19, u^37 + 4*u^35 + 6*u^33 + 8*u^31 + 
#            8*u^29 + 6*u^27 + 5*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          2*u^36 + 2*u^34 + 2*u^32 + 2*u^30, 
#          2*u^35 + 3*u^33 + 4*u^31 + 4*u^29 + u^27 + u^25, 
#          2*u^31 + 3*u^29 + 3*u^27 + 5*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          2*u^32 + 3*u^30 + 4*u^28 + 3*u^26, u^29 + 2*u^27 + u^25 + u^23, 
#          u^33 + 2*u^31 + u^29, 0*u^0, 0*u^0, 
#          2*u^27 + 4*u^25 + 2*u^23 + 2*u^21 + u^19, u^27 + u^25 + u^23, 
#          u^29 + 3*u^27 + 3*u^25, u^29 + u^27, 0*u^0, 
#          3*u^27 + 5*u^25 + 3*u^23 + 2*u^21 + u^19, 
#          u^25 + 2*u^23 + 2*u^21 + u^19, u^27 + u^25, u^24 + 2*u^22, 
#          u^25 + 2*u^23 + 2*u^21 + u^19, u^26 + 2*u^24, u^22, 
#          u^23 + 2*u^21 + u^19, 2*u^22, 0*u^0, 0*u^0, 0*u^0, 2*u^21 + u^19, 
#          u^22, u^21, 2*u^21 + u^19, 0*u^0, u^21 + u^19, 0*u^0, 0*u^0, u^18, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^101 + u^99 + u^97 + 4*u^95 + 4*u^93 + 6*u^91 + 12*u^89 + 11*u^87 + 
#            16*u^85 + 25*u^83 + 24*u^81 + 33*u^79 + 42*u^77 + 40*u^75 + 52*u^
#            73 + 60*u^71 + 57*u^69 + 68*u^67 + 72*u^65 + 67*u^63 + 76*u^61 + 
#            76*u^59 + 67*u^57 + 72*u^55 + 68*u^53 + 57*u^51 + 60*u^49 + 52*u^
#            47 + 40*u^45 + 42*u^43 + 33*u^41 + 24*u^39 + 25*u^37 + 16*u^35 + 
#            11*u^33 + 12*u^31 + 6*u^29 + 4*u^27 + 4*u^25 + u^23 + u^21 + u^19,
#          u^77 + 4*u^73 + 6*u^71 + 8*u^69 + 15*u^67 + 21*u^65 + 24*u^63 + 
#            36*u^61 + 40*u^59 + 43*u^57 + 51*u^55 + 53*u^53 + 49*u^51 + 54*u^
#            49 + 48*u^47 + 40*u^45 + 41*u^43 + 33*u^41 + 24*u^39 + 25*u^37 + 
#            16*u^35 + 11*u^33 + 12*u^31 + 6*u^29 + 4*u^27 + 4*u^25 + u^23 + u^
#            21 + u^19, 3*u^61 + 5*u^59 + 9*u^57 + 16*u^55 + 20*u^53 + 25*u^
#            51 + 33*u^49 + 33*u^47 + 32*u^45 + 35*u^43 + 29*u^41 + 24*u^39 + 
#            24*u^37 + 16*u^35 + 11*u^33 + 12*u^31 + 6*u^29 + 4*u^27 + 4*u^
#            25 + u^23 + u^21 + u^19, u^53 + 2*u^51 + 4*u^49 + 9*u^47 + 9*u^
#            45 + 15*u^43 + 18*u^41 + 16*u^39 + 19*u^37 + 16*u^35 + 10*u^33 + 
#            12*u^31 + 6*u^29 + 4*u^27 + 4*u^25 + u^23 + u^21 + u^19, 
#          u^54 + 2*u^52 + 3*u^50 + 5*u^48 + 5*u^46 + 5*u^44 + 5*u^42 + 3*u^
#            40 + 2*u^38 + u^36, 2*u^50 + u^48 + 4*u^46 + 8*u^44 + 6*u^42 + 
#            9*u^40 + 10*u^38 + 5*u^36 + 6*u^34 + 4*u^32 + u^28, 
#          3*u^43 + 6*u^41 + 8*u^39 + 13*u^37 + 12*u^35 + 10*u^33 + 11*u^31 + 
#            6*u^29 + 4*u^27 + 4*u^25 + u^23 + u^21 + u^19, 
#          u^44 + u^42 + 5*u^40 + 7*u^38 + 5*u^36 + 6*u^34 + 4*u^32 + u^28, 
#          u^45 + 2*u^43 + 3*u^41 + 4*u^39 + 3*u^37 + 2*u^35 + u^33, 
#          u^41 + 3*u^39 + 9*u^37 + 9*u^35 + 10*u^33 + 11*u^31 + 6*u^29 + 4*u^
#            27 + 4*u^25 + u^23 + u^21 + u^19, 
#          2*u^37 + 3*u^35 + 6*u^33 + 8*u^31 + 7*u^29 + 6*u^27 + 4*u^25 + 2*u^
#            23 + u^21 + u^19, u^37 + 3*u^35 + 3*u^33 + 6*u^31 + 6*u^29 + 3*u^
#            27 + 4*u^25 + u^23 + u^21 + u^19, u^36 + u^34 + u^32 + u^30, 
#          2*u^35 + u^33 + 3*u^31 + 4*u^29 + u^25, 
#          u^31 + 3*u^29 + u^27 + 4*u^25 + u^23 + u^21 + u^19, 
#          u^32 + 2*u^30 + 3*u^28 + 2*u^26, u^29 + 2*u^27 + u^23, u^33 + u^29, 
#          0*u^0, u^29 + u^25 + u^23 + u^19, u^27 + 3*u^25 + u^23 + u^21 + u^19
#            , 2*u^27 + u^23, 2*u^27 + 2*u^25, u^27, 0*u^0, 
#          3*u^27 + 4*u^25 + 2*u^23 + u^21 + u^19, u^23 + u^21 + u^19, 
#          u^27 + u^25, u^22, 2*u^25 + u^23 + u^21 + u^19, u^24, u^26 + u^22, 
#          u^21 + u^19, u^22, u^19, 0*u^0, u^22 + u^20, u^21 + u^19, 0*u^0, 
#          u^23 + u^21, u^21 + u^19, u^19, u^19, 0*u^0, u^19 + u^17, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^100 + 2*u^96 + 3*u^94 + 3*u^92 + 7*u^90 + 11*u^88 + 8*u^86 + 18*u^
#            84 + 21*u^82 + 21*u^80 + 33*u^78 + 37*u^76 + 33*u^74 + 53*u^72 + 
#            49*u^70 + 49*u^68 + 66*u^66 + 59*u^64 + 56*u^62 + 74*u^60 + 56*u^
#            58 + 59*u^56 + 66*u^54 + 49*u^52 + 49*u^50 + 53*u^48 + 33*u^46 + 
#            37*u^44 + 33*u^42 + 21*u^40 + 21*u^38 + 18*u^36 + 8*u^34 + 11*u^
#            32 + 7*u^30 + 3*u^28 + 3*u^26 + 2*u^24 + u^20, 
#          u^78 + u^76 + u^74 + 6*u^72 + 6*u^70 + 9*u^68 + 18*u^66 + 19*u^64 + 
#            24*u^62 + 37*u^60 + 32*u^58 + 40*u^56 + 48*u^54 + 40*u^52 + 43*u^
#            50 + 47*u^48 + 32*u^46 + 36*u^44 + 32*u^42 + 21*u^40 + 21*u^38 + 
#            18*u^36 + 8*u^34 + 11*u^32 + 7*u^30 + 3*u^28 + 3*u^26 + 2*u^
#            24 + u^20, u^62 + 5*u^60 + 4*u^58 + 12*u^56 + 16*u^54 + 17*u^52 + 
#            24*u^50 + 29*u^48 + 23*u^46 + 30*u^44 + 26*u^42 + 20*u^40 + 20*u^
#            38 + 17*u^36 + 8*u^34 + 11*u^32 + 7*u^30 + 3*u^28 + 3*u^26 + 2*u^
#            24 + u^20, u^54 + 2*u^52 + 3*u^50 + 7*u^48 + 9*u^46 + 10*u^44 + 
#            16*u^42 + 13*u^40 + 13*u^38 + 16*u^36 + 7*u^34 + 10*u^32 + 7*u^
#            30 + 3*u^28 + 3*u^26 + 2*u^24 + u^20, 
#          u^53 + u^51 + u^49 + 3*u^47 + 3*u^43 + u^41 + u^39 + u^37, 
#          u^55 + 3*u^51 + 4*u^49 + 3*u^47 + 9*u^45 + 9*u^43 + 5*u^41 + 13*u^
#            39 + 5*u^37 + 5*u^35 + 6*u^33 + u^31 + u^29 + u^27, 
#          u^44 + 4*u^42 + 5*u^40 + 7*u^38 + 10*u^36 + 6*u^34 + 9*u^32 + 6*u^
#            30 + 3*u^28 + 3*u^26 + 2*u^24 + u^20, 
#          2*u^45 + 2*u^43 + 3*u^41 + 8*u^39 + 4*u^37 + 5*u^35 + 6*u^33 + u^
#            31 + u^29 + u^27, u^44 + 2*u^40 + u^38 + u^36 + u^34, 
#          u^42 + u^40 + 5*u^38 + 6*u^36 + 5*u^34 + 9*u^32 + 6*u^30 + 3*u^28 + 
#            3*u^26 + 2*u^24 + u^20, u^40 + 2*u^38 + 4*u^36 + 6*u^34 + 8*u^
#            32 + 8*u^30 + 7*u^28 + 4*u^26 + 3*u^24 + u^22 + u^20, 
#          u^36 + u^34 + 2*u^32 + 5*u^30 + 2*u^28 + 2*u^26 + 2*u^24 + u^20, 
#          u^35 + u^31, u^36 + 2*u^34 + 5*u^30 + u^28 + u^26 + u^24, 
#          3*u^30 + u^28 + 2*u^26 + 2*u^24 + u^20, u^33 + u^29 + u^27, 
#          u^32 + u^30 + 4*u^28 + u^26 + u^24 + u^22, 0*u^0, 2*u^30 + u^26, 
#          u^30 + u^28 + 2*u^24 + u^20 + u^18, u^26 + u^24 + u^20, 
#          2*u^28 + u^26 + u^24 + u^22, u^26, 0*u^0, u^30 + u^26, 
#          u^28 + u^26 + 2*u^24 + u^22 + u^20, u^24 + u^22 + u^20, 0*u^0, 
#          u^25 + u^23 + u^21, u^24 + u^20, 0*u^0, u^27 + u^23 + u^21, u^20, 
#          u^23 + u^21, u^20 + u^18, u^21, u^21, u^20, 0*u^0, 0*u^0, u^20, 
#          u^20 + u^18, 0*u^0, 0*u^0, u^18, 0*u^0, u^18, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^16, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^100 + u^96 + u^94 + 2*u^92 + 2*u^90 + 5*u^88 + 3*u^86 + 7*u^84 + 
#            8*u^82 + 9*u^80 + 11*u^78 + 15*u^76 + 11*u^74 + 20*u^72 + 18*u^
#            70 + 18*u^68 + 23*u^66 + 22*u^64 + 19*u^62 + 28*u^60 + 19*u^58 + 
#            22*u^56 + 23*u^54 + 18*u^52 + 18*u^50 + 20*u^48 + 11*u^46 + 15*u^
#            44 + 11*u^42 + 9*u^40 + 8*u^38 + 7*u^36 + 3*u^34 + 5*u^32 + 2*u^
#            30 + 2*u^28 + u^26 + u^24 + u^20, 
#          2*u^72 + 2*u^70 + 3*u^68 + 6*u^66 + 7*u^64 + 8*u^62 + 14*u^60 + 
#            11*u^58 + 15*u^56 + 17*u^54 + 15*u^52 + 16*u^50 + 18*u^48 + 11*u^
#            46 + 15*u^44 + 11*u^42 + 9*u^40 + 8*u^38 + 7*u^36 + 3*u^34 + 5*u^
#            32 + 2*u^30 + 2*u^28 + u^26 + u^24 + u^20, 
#          2*u^60 + u^58 + 5*u^56 + 5*u^54 + 7*u^52 + 9*u^50 + 12*u^48 + 8*u^
#            46 + 13*u^44 + 9*u^42 + 9*u^40 + 8*u^38 + 7*u^36 + 3*u^34 + 5*u^
#            32 + 2*u^30 + 2*u^28 + u^26 + u^24 + u^20, 
#          u^50 + u^48 + 3*u^46 + 3*u^44 + 5*u^42 + 6*u^40 + 5*u^38 + 7*u^36 + 
#            3*u^34 + 5*u^32 + 2*u^30 + 2*u^28 + u^26 + u^24 + u^20, 
#          u^53 + u^51 + u^49 + 4*u^47 + 4*u^43 + u^41 + u^39 + u^37, 
#          u^45 + 3*u^43 + 5*u^39 + u^37 + 2*u^35 + 2*u^33, 
#          u^42 + 3*u^40 + 3*u^38 + 5*u^36 + 3*u^34 + 5*u^32 + 2*u^30 + 2*u^
#            28 + u^26 + u^24 + u^20, 3*u^39 + u^37 + 2*u^35 + 2*u^33, 
#          u^44 + 3*u^40 + u^38 + u^36 + u^34, 
#          3*u^38 + 3*u^36 + 3*u^34 + 5*u^32 + 2*u^30 + 2*u^28 + u^26 + u^
#            24 + u^20, u^34 + 2*u^32 + 2*u^30 + 2*u^28 + u^26 + u^24 + u^20, 
#          u^36 + u^34 + 2*u^32 + 2*u^30 + 2*u^28 + u^26 + u^24 + u^20, 
#          u^35 + u^31, u^34 + 2*u^30, u^30 + u^26 + u^24 + u^20, 
#          u^29 + 2*u^27, 0*u^0, u^32 + u^30, 0*u^0, 0*u^0, u^26 + u^24 + u^20,
#          0*u^0, 2*u^26, u^28, 0*u^0, 2*u^26 + u^24 + u^20, u^20, u^26, 
#          0*u^0, u^24 + u^20, u^25, 0*u^0, u^20, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^20, 0*u^0, u^22, u^20, 0*u^0, u^20, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ 2*u^100 + u^98 + 5*u^96 + 7*u^94 + 10*u^92 + 15*u^90 + 24*u^88 + 25*u^
#            86 + 41*u^84 + 48*u^82 + 57*u^80 + 72*u^78 + 86*u^76 + 89*u^74 + 
#            115*u^72 + 117*u^70 + 126*u^68 + 141*u^66 + 143*u^64 + 141*u^62 + 
#            158*u^60 + 141*u^58 + 143*u^56 + 141*u^54 + 126*u^52 + 117*u^50 + 
#            115*u^48 + 89*u^46 + 86*u^44 + 72*u^42 + 57*u^40 + 48*u^38 + 41*u^
#            36 + 25*u^34 + 24*u^32 + 15*u^30 + 10*u^28 + 7*u^26 + 5*u^24 + u^
#            22 + 2*u^20, u^78 + 2*u^76 + 4*u^74 + 11*u^72 + 15*u^70 + 24*u^
#            68 + 37*u^66 + 47*u^64 + 59*u^62 + 79*u^60 + 82*u^58 + 96*u^56 + 
#            104*u^54 + 102*u^52 + 102*u^50 + 104*u^48 + 85*u^46 + 84*u^44 + 
#            71*u^42 + 57*u^40 + 48*u^38 + 41*u^36 + 25*u^34 + 24*u^32 + 15*u^
#            30 + 10*u^28 + 7*u^26 + 5*u^24 + u^22 + 2*u^20, 
#          u^64 + 3*u^62 + 9*u^60 + 12*u^58 + 26*u^56 + 34*u^54 + 46*u^52 + 
#            56*u^50 + 67*u^48 + 61*u^46 + 69*u^44 + 60*u^42 + 53*u^40 + 46*u^
#            38 + 40*u^36 + 25*u^34 + 24*u^32 + 15*u^30 + 10*u^28 + 7*u^26 + 
#            5*u^24 + u^22 + 2*u^20, u^54 + 4*u^52 + 7*u^50 + 15*u^48 + 20*u^
#            46 + 27*u^44 + 33*u^42 + 35*u^40 + 33*u^38 + 35*u^36 + 23*u^34 + 
#            23*u^32 + 15*u^30 + 10*u^28 + 7*u^26 + 5*u^24 + u^22 + 2*u^20, 
#          u^55 + 2*u^53 + 3*u^51 + 4*u^49 + 8*u^47 + 4*u^45 + 8*u^43 + 4*u^
#            41 + 3*u^39 + 2*u^37 + u^35, u^55 + 4*u^51 + 6*u^49 + 9*u^47 + 
#            15*u^45 + 19*u^43 + 16*u^41 + 22*u^39 + 15*u^37 + 12*u^35 + 10*u^
#            33 + 4*u^31 + 2*u^29 + u^27, 3*u^44 + 8*u^42 + 15*u^40 + 19*u^
#            38 + 24*u^36 + 19*u^34 + 21*u^32 + 14*u^30 + 10*u^28 + 7*u^26 + 
#            5*u^24 + u^22 + 2*u^20, 2*u^45 + 4*u^43 + 8*u^41 + 14*u^39 + 12*u^
#            37 + 11*u^35 + 10*u^33 + 4*u^31 + 2*u^29 + u^27, 
#          2*u^44 + 2*u^42 + 6*u^40 + 4*u^38 + 3*u^36 + 2*u^34 + u^32, 
#          u^42 + 4*u^40 + 11*u^38 + 17*u^36 + 16*u^34 + 20*u^32 + 14*u^30 + 
#            10*u^28 + 7*u^26 + 5*u^24 + u^22 + 2*u^20, 
#          2*u^40 + 4*u^38 + 9*u^36 + 12*u^34 + 17*u^32 + 16*u^30 + 15*u^28 + 
#            10*u^26 + 7*u^24 + 2*u^22 + 2*u^20, 
#          3*u^36 + 5*u^34 + 8*u^32 + 9*u^30 + 8*u^28 + 6*u^26 + 5*u^24 + u^
#            22 + 2*u^20, u^37 + 2*u^35 + u^33 + 2*u^31 + u^29, 
#          u^36 + 4*u^34 + 3*u^32 + 7*u^30 + 3*u^28 + 2*u^26 + u^24, 
#          u^32 + 4*u^30 + 4*u^28 + 5*u^26 + 5*u^24 + u^22 + 2*u^20, 
#          2*u^33 + 3*u^31 + 3*u^29 + 3*u^27 + u^25, 
#          u^32 + 2*u^30 + 5*u^28 + 3*u^26 + 2*u^24 + u^22, 0*u^0, 
#          u^30 + u^26, u^28 + 2*u^24 + u^20, u^28 + 3*u^26 + 4*u^24 + u^22 + 
#            2*u^20, 2*u^28 + 2*u^26 + 2*u^24 + u^22, 
#          u^30 + u^28 + 3*u^26 + u^24, 0*u^0, u^26, 
#          u^28 + 4*u^26 + 5*u^24 + 2*u^22 + 2*u^20, 
#          u^26 + 3*u^24 + 2*u^22 + 2*u^20, u^26, 
#          u^27 + 2*u^25 + 3*u^23 + 2*u^21, 2*u^24 + u^22 + 2*u^20, 
#          u^25 + u^23, u^23 + u^21, u^24 + u^22 + 2*u^20, 2*u^23 + 2*u^21, 
#          u^20, u^21, u^21, 2*u^20, u^23 + u^21, 0*u^0, 2*u^20, u^20, u^20, 
#          0*u^0, 0*u^0, u^19 + u^17, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^16, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^100 + 3*u^96 + 3*u^94 + 4*u^92 + 8*u^90 + 12*u^88 + 11*u^86 + 22*u^
#            84 + 24*u^82 + 28*u^80 + 39*u^78 + 44*u^76 + 45*u^74 + 63*u^72 + 
#            60*u^70 + 65*u^68 + 77*u^66 + 75*u^64 + 73*u^62 + 86*u^60 + 73*u^
#            58 + 75*u^56 + 77*u^54 + 65*u^52 + 60*u^50 + 63*u^48 + 45*u^46 + 
#            44*u^44 + 39*u^42 + 28*u^40 + 24*u^38 + 22*u^36 + 11*u^34 + 12*u^
#            32 + 8*u^30 + 4*u^28 + 3*u^26 + 3*u^24 + u^20, 
#          u^78 + u^76 + 2*u^74 + 7*u^72 + 8*u^70 + 13*u^68 + 21*u^66 + 25*u^
#            64 + 31*u^62 + 43*u^60 + 42*u^58 + 50*u^56 + 56*u^54 + 52*u^52 + 
#            52*u^50 + 56*u^48 + 43*u^46 + 43*u^44 + 38*u^42 + 28*u^40 + 24*u^
#            38 + 22*u^36 + 11*u^34 + 12*u^32 + 8*u^30 + 4*u^28 + 3*u^26 + 3*u^
#            24 + u^20, u^64 + 2*u^62 + 5*u^60 + 6*u^58 + 14*u^56 + 18*u^54 + 
#            23*u^52 + 28*u^50 + 35*u^48 + 30*u^46 + 35*u^44 + 31*u^42 + 26*u^
#            40 + 23*u^38 + 21*u^36 + 11*u^34 + 12*u^32 + 8*u^30 + 4*u^28 + 
#            3*u^26 + 3*u^24 + u^20, 2*u^54 + 3*u^52 + 4*u^50 + 11*u^48 + 11*u^
#            46 + 15*u^44 + 19*u^42 + 18*u^40 + 16*u^38 + 19*u^36 + 10*u^34 + 
#            11*u^32 + 8*u^30 + 4*u^28 + 3*u^26 + 3*u^24 + u^20, u^47 + u^43, 
#          u^55 + 4*u^51 + 4*u^49 + 6*u^47 + 10*u^45 + 11*u^43 + 9*u^41 + 13*u^
#            39 + 8*u^37 + 6*u^35 + 6*u^33 + 2*u^31 + u^29 + u^27, 
#          u^48 + 2*u^44 + 5*u^42 + 8*u^40 + 9*u^38 + 12*u^36 + 8*u^34 + 10*u^
#            32 + 7*u^30 + 4*u^28 + 3*u^26 + 3*u^24 + u^20, 
#          u^45 + 2*u^43 + 5*u^41 + 8*u^39 + 6*u^37 + 5*u^35 + 6*u^33 + 2*u^
#            31 + u^29 + u^27, u^40, 2*u^40 + 5*u^38 + 8*u^36 + 6*u^34 + 9*u^
#            32 + 7*u^30 + 4*u^28 + 3*u^26 + 3*u^24 + u^20, 
#          2*u^40 + 4*u^38 + 7*u^36 + 8*u^34 + 11*u^32 + 9*u^30 + 8*u^28 + 5*u^
#            26 + 4*u^24 + u^22 + u^20, 2*u^36 + 2*u^34 + 3*u^32 + 5*u^30 + 
#            3*u^28 + 2*u^26 + 3*u^24 + u^20, 0*u^0, 
#          u^36 + 2*u^34 + u^32 + 4*u^30 + u^28 + u^26 + u^24, 
#          u^32 + 3*u^30 + 2*u^28 + 2*u^26 + 3*u^24 + u^20, 0*u^0, 
#          u^34 + 2*u^32 + 2*u^30 + 4*u^28 + 2*u^26 + u^24 + u^22, 0*u^0, 
#          u^30, u^30 + u^28 + 3*u^24 + u^20 + u^18, 
#          u^28 + u^26 + 2*u^24 + u^20, 2*u^28 + u^26 + u^24 + u^22, 0*u^0, 
#          0*u^0, 0*u^0, u^28 + u^26 + 2*u^24 + u^22 + u^20, 
#          u^28 + u^26 + 2*u^24 + u^22 + u^20, 0*u^0, 
#          u^27 + u^25 + u^23 + u^21, u^24 + u^20, 0*u^0, u^21, u^24 + u^20, 
#          u^21, u^24 + u^20 + u^18, u^21, 0*u^0, u^20, 0*u^0, 0*u^0, 0*u^0, 
#          u^18, 0*u^0, u^18, u^18, 0*u^0, u^18, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^98 + 2*u^96 + 4*u^94 + 5*u^92 + 9*u^90 + 12*u^88 + 17*u^86 + 23*u^
#            84 + 29*u^82 + 35*u^80 + 45*u^78 + 51*u^76 + 61*u^74 + 68*u^72 + 
#            75*u^70 + 82*u^68 + 88*u^66 + 91*u^64 + 95*u^62 + 94*u^60 + 95*u^
#            58 + 91*u^56 + 88*u^54 + 82*u^52 + 75*u^50 + 68*u^48 + 61*u^46 + 
#            51*u^44 + 45*u^42 + 35*u^40 + 29*u^38 + 23*u^36 + 17*u^34 + 12*u^
#            32 + 9*u^30 + 5*u^28 + 4*u^26 + 2*u^24 + u^22, 
#          u^78 + 2*u^76 + 4*u^74 + 7*u^72 + 11*u^70 + 17*u^68 + 24*u^66 + 
#            31*u^64 + 40*u^62 + 47*u^60 + 55*u^58 + 60*u^56 + 64*u^54 + 65*u^
#            52 + 64*u^50 + 61*u^48 + 57*u^46 + 49*u^44 + 44*u^42 + 35*u^40 + 
#            29*u^38 + 23*u^36 + 17*u^34 + 12*u^32 + 9*u^30 + 5*u^28 + 4*u^
#            26 + 2*u^24 + u^22, u^64 + 3*u^62 + 5*u^60 + 10*u^58 + 15*u^56 + 
#            22*u^54 + 28*u^52 + 34*u^50 + 37*u^48 + 40*u^46 + 38*u^44 + 37*u^
#            42 + 31*u^40 + 27*u^38 + 22*u^36 + 17*u^34 + 12*u^32 + 9*u^30 + 
#            5*u^28 + 4*u^26 + 2*u^24 + u^22, 2*u^54 + 3*u^52 + 7*u^50 + 10*u^
#            48 + 15*u^46 + 17*u^44 + 21*u^42 + 20*u^40 + 20*u^38 + 18*u^36 + 
#            15*u^34 + 11*u^32 + 9*u^30 + 5*u^28 + 4*u^26 + 2*u^24 + u^22, 
#          0*u^0, u^57 + u^55 + 3*u^53 + 5*u^51 + 8*u^49 + 11*u^47 + 13*u^45 + 
#            15*u^43 + 15*u^41 + 14*u^39 + 12*u^37 + 9*u^35 + 6*u^33 + 4*u^
#            31 + 2*u^29 + u^27, u^46 + 2*u^44 + 5*u^42 + 7*u^40 + 10*u^38 + 
#            11*u^36 + 11*u^34 + 9*u^32 + 8*u^30 + 5*u^28 + 4*u^26 + 2*u^
#            24 + u^22, u^47 + 2*u^45 + 5*u^43 + 7*u^41 + 9*u^39 + 9*u^37 + 
#            8*u^35 + 6*u^33 + 4*u^31 + 2*u^29 + u^27, 0*u^0, 
#          u^42 + 2*u^40 + 5*u^38 + 7*u^36 + 9*u^34 + 8*u^32 + 8*u^30 + 5*u^
#            28 + 4*u^26 + 2*u^24 + u^22, u^42 + 2*u^40 + 6*u^38 + 9*u^36 + 
#            12*u^34 + 13*u^32 + 13*u^30 + 10*u^28 + 8*u^26 + 4*u^24 + 2*u^22, 
#          u^36 + 2*u^34 + 3*u^32 + 4*u^30 + 3*u^28 + 3*u^26 + 2*u^24 + u^22, 
#          0*u^0, u^36 + 2*u^34 + 3*u^32 + 3*u^30 + 3*u^28 + 2*u^26 + u^24, 
#          u^32 + 2*u^30 + 2*u^28 + 3*u^26 + 2*u^24 + u^22, 0*u^0, 
#          u^34 + 3*u^32 + 4*u^30 + 4*u^28 + 4*u^26 + 2*u^24 + u^22, 0*u^0, 
#          u^34 + 2*u^32 + 2*u^30 + 2*u^28 + u^26 + u^24, 
#          u^30 + u^28 + 2*u^26 + 2*u^24 + 2*u^22 + u^20 + u^18, 
#          u^26 + u^24 + u^22, u^30 + u^28 + 2*u^26 + 2*u^24 + u^22, 0*u^0, 
#          0*u^0, u^30 + u^28 + u^26 + u^24, u^26 + 2*u^24 + 2*u^22, 
#          u^26 + 2*u^24 + 2*u^22, 0*u^0, u^29 + u^27 + 2*u^25 + 3*u^23 + 2*u^
#            21, u^22, 0*u^0, u^23 + u^21, u^22, u^23 + 2*u^21, 
#          u^22 + u^20 + u^18, u^23 + 2*u^21 + u^19, 0*u^0, 0*u^0, u^21, 
#          0*u^0, 0*u^0, u^20 + u^18, 0*u^0, u^18, 0*u^0, 0*u^0, 
#          u^20 + u^18 + u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^16, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^96 + u^92 + u^90 + 2*u^88 + u^86 + 3*u^84 + 3*u^82 + 4*u^80 + 4*u^
#            78 + 6*u^76 + 4*u^74 + 9*u^72 + 7*u^70 + 7*u^68 + 9*u^66 + 9*u^
#            64 + 7*u^62 + 12*u^60 + 7*u^58 + 9*u^56 + 9*u^54 + 7*u^52 + 7*u^
#            50 + 9*u^48 + 4*u^46 + 6*u^44 + 4*u^42 + 4*u^40 + 3*u^38 + 3*u^
#            36 + u^34 + 2*u^32 + u^30 + u^28 + u^24, 
#          u^72 + u^70 + u^68 + 2*u^66 + 3*u^64 + 3*u^62 + 6*u^60 + 4*u^58 + 
#            6*u^56 + 7*u^54 + 6*u^52 + 6*u^50 + 8*u^48 + 4*u^46 + 6*u^44 + 
#            4*u^42 + 4*u^40 + 3*u^38 + 3*u^36 + u^34 + 2*u^32 + u^30 + u^
#            28 + u^24, u^60 + 2*u^56 + 2*u^54 + 3*u^52 + 3*u^50 + 6*u^48 + 
#            3*u^46 + 5*u^44 + 3*u^42 + 4*u^40 + 3*u^38 + 3*u^36 + u^34 + 2*u^
#            32 + u^30 + u^28 + u^24, u^48 + u^46 + u^44 + 2*u^42 + 2*u^40 + 
#            2*u^38 + 3*u^36 + u^34 + 2*u^32 + u^30 + u^28 + u^24, 
#          u^51 + 2*u^47 + 2*u^43 + u^39, u^45 + u^43 + 2*u^39 + u^35 + u^33, 
#          u^42 + u^40 + u^38 + 2*u^36 + u^34 + 2*u^32 + u^30 + u^28 + u^24, 
#          u^39 + u^35 + u^33, u^44 + u^40 + u^36, 
#          u^38 + u^36 + u^34 + 2*u^32 + u^30 + u^28 + u^24, 
#          u^32 + u^30 + u^28 + u^24, u^36 + u^32 + u^30 + u^28 + u^24, 0*u^0, 
#          u^30, u^24, u^29 + u^27, 0*u^0, 0*u^0, 0*u^0, u^24, u^24, 0*u^0, 
#          u^26, 0*u^0, 0*u^0, u^26 + u^24, 0*u^0, u^26, 0*u^0, u^24, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^21, 0*u^0, 0*u^0, u^22, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^18, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^95 + u^93 + u^91 + 4*u^89 + 3*u^87 + 4*u^85 + 9*u^83 + 7*u^81 + 
#            10*u^79 + 16*u^77 + 12*u^75 + 17*u^73 + 23*u^71 + 17*u^69 + 24*u^
#            67 + 27*u^65 + 20*u^63 + 28*u^61 + 28*u^59 + 20*u^57 + 27*u^55 + 
#            24*u^53 + 17*u^51 + 23*u^49 + 17*u^47 + 12*u^45 + 16*u^43 + 10*u^
#            41 + 7*u^39 + 9*u^37 + 4*u^35 + 3*u^33 + 4*u^31 + u^29 + u^27 + u^
#            25, u^77 + 2*u^73 + 3*u^71 + 3*u^69 + 6*u^67 + 9*u^65 + 7*u^63 + 
#            14*u^61 + 14*u^59 + 13*u^57 + 18*u^55 + 18*u^53 + 14*u^51 + 20*u^
#            49 + 15*u^47 + 12*u^45 + 15*u^43 + 10*u^41 + 7*u^39 + 9*u^37 + 
#            4*u^35 + 3*u^33 + 4*u^31 + u^29 + u^27 + u^25, 
#          2*u^61 + 2*u^59 + 3*u^57 + 6*u^55 + 6*u^53 + 7*u^51 + 11*u^49 + 9*u^
#            47 + 9*u^45 + 12*u^43 + 8*u^41 + 7*u^39 + 8*u^37 + 4*u^35 + 3*u^
#            33 + 4*u^31 + u^29 + u^27 + u^25, u^53 + u^51 + 2*u^49 + 4*u^47 + 
#            2*u^45 + 6*u^43 + 5*u^41 + 4*u^39 + 6*u^37 + 4*u^35 + 2*u^33 + 
#            4*u^31 + u^29 + u^27 + u^25, 0*u^0, 
#          u^56 + u^54 + u^52 + 4*u^50 + 2*u^48 + 4*u^46 + 6*u^44 + 3*u^42 + 
#            5*u^40 + 5*u^38 + 2*u^36 + 3*u^34 + 2*u^32 + u^28, 
#          u^43 + u^41 + u^39 + 3*u^37 + 2*u^35 + 2*u^33 + 3*u^31 + u^29 + u^
#            27 + u^25, u^46 + 2*u^44 + u^42 + 3*u^40 + 3*u^38 + 2*u^36 + 3*u^
#            34 + 2*u^32 + u^28, 0*u^0, u^41 + 2*u^37 + u^35 + 2*u^33 + 3*u^
#            31 + u^29 + u^27 + u^25, u^41 + u^39 + 2*u^37 + 3*u^35 + 4*u^33 + 
#            4*u^31 + 4*u^29 + 3*u^27 + u^25 + u^23, u^31 + u^29 + u^25, 
#          0*u^0, u^35 + u^31 + 2*u^29 + u^25, u^29 + u^25, 0*u^0, 
#          u^33 + 2*u^29 + 2*u^27 + u^23, 0*u^0, 
#          u^35 + 2*u^31 + u^29 + u^27 + u^25, u^29 + u^25 + u^23 + u^19, 
#          0*u^0, u^27 + u^23, 0*u^0, 0*u^0, u^31 + u^27 + u^25, u^23, u^23, 
#          0*u^0, u^26 + u^24 + u^22, 0*u^0, 0*u^0, u^22, 0*u^0, u^22, u^19, 
#          u^22 + u^20, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^19, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^19 + u^17, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^88 + u^84 + 2*u^82 + u^80 + 2*u^78 + 3*u^76 + 5*u^72 + 4*u^70 + u^
#            68 + 6*u^66 + 3*u^64 + 2*u^62 + 8*u^60 + 2*u^58 + 3*u^56 + 6*u^
#            54 + u^52 + 4*u^50 + 5*u^48 + 3*u^44 + 2*u^42 + u^40 + 2*u^38 + u^
#            36 + u^32, u^72 + u^70 + 2*u^66 + u^64 + u^62 + 4*u^60 + u^58 + 
#            2*u^56 + 4*u^54 + u^52 + 3*u^50 + 4*u^48 + 3*u^44 + 2*u^42 + u^
#            40 + 2*u^38 + u^36 + u^32, u^60 + u^56 + u^54 + 2*u^50 + 2*u^48 + 
#            2*u^44 + u^42 + u^40 + 2*u^38 + u^36 + u^32, 
#          u^48 + u^42 + u^38 + u^36 + u^32, 0*u^0, 
#          u^55 + u^51 + u^49 + 2*u^45 + u^43 + 2*u^39 + u^35 + u^33, u^32, 
#          u^45 + u^39 + u^35 + u^33, 0*u^0, u^32, 
#          u^40 + u^34 + u^32 + u^30 + u^28, 0*u^0, 0*u^0, u^30, 0*u^0, 0*u^0, 
#          u^28, 0*u^0, u^36 + u^30 + u^26, u^24, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^26, 0*u^0, 0*u^0, 0*u^0, u^25, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^21, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^18, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^16, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^99 + 2*u^97 + 4*u^95 + 8*u^93 + 13*u^91 + 20*u^89 + 31*u^87 + 43*u^
#            85 + 58*u^83 + 78*u^81 + 98*u^79 + 121*u^77 + 148*u^75 + 173*u^
#            73 + 199*u^71 + 227*u^69 + 249*u^67 + 270*u^65 + 290*u^63 + 301*u^
#            61 + 309*u^59 + 314*u^57 + 309*u^55 + 301*u^53 + 290*u^51 + 270*u^
#            49 + 249*u^47 + 227*u^45 + 199*u^43 + 173*u^41 + 148*u^39 + 121*u^
#            37 + 98*u^35 + 78*u^33 + 58*u^31 + 43*u^29 + 31*u^27 + 20*u^25 + 
#            13*u^23 + 8*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^77 + 3*u^75 + 7*u^73 + 14*u^71 + 26*u^69 + 41*u^67 + 62*u^65 + 
#            87*u^63 + 114*u^61 + 142*u^59 + 171*u^57 + 193*u^55 + 212*u^53 + 
#            224*u^51 + 226*u^49 + 221*u^47 + 211*u^45 + 191*u^43 + 170*u^41 + 
#            147*u^39 + 121*u^37 + 98*u^35 + 78*u^33 + 58*u^31 + 43*u^29 + 
#            31*u^27 + 20*u^25 + 13*u^23 + 8*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^63 + 5*u^61 + 12*u^59 + 25*u^57 + 43*u^55 + 66*u^53 + 91*u^51 + 
#            114*u^49 + 133*u^47 + 145*u^45 + 147*u^43 + 142*u^41 + 131*u^39 + 
#            113*u^37 + 95*u^35 + 77*u^33 + 58*u^31 + 43*u^29 + 31*u^27 + 20*u^
#            25 + 13*u^23 + 8*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^55 + 3*u^53 + 10*u^51 + 18*u^49 + 32*u^47 + 48*u^45 + 63*u^43 + 
#            75*u^41 + 84*u^39 + 83*u^37 + 78*u^35 + 69*u^33 + 55*u^31 + 42*u^
#            29 + 31*u^27 + 20*u^25 + 13*u^23 + 8*u^21 + 4*u^19 + 2*u^17 + u^15
#            , u^52 + 2*u^50 + 4*u^48 + 5*u^46 + 6*u^44 + 6*u^42 + 5*u^40 + 
#            4*u^38 + 2*u^36 + u^34, u^54 + 3*u^52 + 7*u^50 + 15*u^48 + 23*u^
#            46 + 33*u^44 + 41*u^42 + 46*u^40 + 45*u^38 + 42*u^36 + 33*u^34 + 
#            24*u^32 + 15*u^30 + 8*u^28 + 3*u^26 + u^24, 
#          2*u^45 + 7*u^43 + 18*u^41 + 32*u^39 + 44*u^37 + 51*u^35 + 53*u^33 + 
#            47*u^31 + 39*u^29 + 30*u^27 + 20*u^25 + 13*u^23 + 8*u^21 + 4*u^
#            19 + 2*u^17 + u^15, 3*u^44 + 9*u^42 + 19*u^40 + 27*u^38 + 32*u^
#            36 + 29*u^34 + 23*u^32 + 15*u^30 + 8*u^28 + 3*u^26 + u^24, 
#          u^43 + 3*u^41 + 5*u^39 + 5*u^37 + 4*u^35 + 2*u^33 + u^31, 
#          2*u^41 + 9*u^39 + 21*u^37 + 34*u^35 + 43*u^33 + 43*u^31 + 38*u^29 + 
#            30*u^27 + 20*u^25 + 13*u^23 + 8*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^41 + 6*u^39 + 14*u^37 + 27*u^35 + 39*u^33 + 48*u^31 + 49*u^29 + 
#            44*u^27 + 32*u^25 + 21*u^23 + 11*u^21 + 5*u^19 + 2*u^17 + u^15, 
#          u^37 + 4*u^35 + 12*u^33 + 18*u^31 + 22*u^29 + 22*u^27 + 17*u^25 + 
#            12*u^23 + 8*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^36 + 2*u^34 + 3*u^32 + 2*u^30 + u^28, 
#          2*u^35 + 7*u^33 + 11*u^31 + 13*u^29 + 11*u^27 + 7*u^25 + 3*u^23 + u^
#            21, 2*u^33 + 6*u^31 + 12*u^29 + 16*u^27 + 15*u^25 + 12*u^23 + 8*u^
#            21 + 4*u^19 + 2*u^17 + u^15, u^32 + 2*u^30 + 3*u^28 + 2*u^26 + u^
#            24, u^33 + 6*u^31 + 11*u^29 + 14*u^27 + 12*u^25 + 8*u^23 + 3*u^
#            21 + u^19, 0*u^0, u^31 + 2*u^29 + 2*u^27 + u^25, 
#          u^33 + 2*u^31 + 4*u^29 + 6*u^27 + 7*u^25 + 7*u^23 + 7*u^21 + 5*u^
#            19 + 3*u^17 + 2*u^15, u^29 + 5*u^27 + 8*u^25 + 9*u^23 + 7*u^21 + 
#            4*u^19 + 2*u^17 + u^15, u^29 + 6*u^27 + 8*u^25 + 7*u^23 + 3*u^
#            21 + u^19, 2*u^27 + 2*u^25 + u^23, 0*u^0, u^27 + u^25, 
#          4*u^27 + 9*u^25 + 12*u^23 + 9*u^21 + 5*u^19 + 2*u^17 + u^15, 
#          2*u^27 + 6*u^25 + 10*u^23 + 9*u^21 + 5*u^19 + 2*u^17 + u^15, 0*u^0, 
#          3*u^26 + 7*u^24 + 8*u^22 + 5*u^20 + 2*u^18, 
#          u^25 + 4*u^23 + 6*u^21 + 4*u^19 + 2*u^17 + u^15, u^24 + u^22, 
#          u^24 + 3*u^22 + 2*u^20 + u^18, 
#          3*u^23 + 6*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^24 + 4*u^22 + 4*u^20 + 2*u^18, 
#          2*u^23 + 5*u^21 + 5*u^19 + 3*u^17 + 2*u^15, 
#          2*u^22 + 3*u^20 + 2*u^18, u^22 + u^20, 
#          2*u^21 + 3*u^19 + 2*u^17 + u^15, u^22 + 2*u^20 + u^18, 0*u^0, 
#          u^21 + 3*u^19 + 2*u^17 + u^15, u^21 + 4*u^19 + 3*u^17 + 2*u^15, 
#          u^19 + u^17 + u^15, 2*u^19 + 2*u^17 + 3*u^15, u^19 + u^17 + 2*u^15, 
#          0*u^0, u^19 + u^17 + u^15, u^19 + u^17 + u^15, u^15, 0*u^0, 0*u^0, 
#          0*u^0, u^15, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^99 + 2*u^97 + 3*u^95 + 7*u^93 + 11*u^91 + 15*u^89 + 25*u^87 + 34*u^
#            85 + 43*u^83 + 61*u^81 + 75*u^79 + 89*u^77 + 114*u^75 + 130*u^
#            73 + 146*u^71 + 173*u^69 + 185*u^67 + 198*u^65 + 220*u^63 + 222*u^
#            61 + 227*u^59 + 238*u^57 + 227*u^55 + 222*u^53 + 220*u^51 + 198*u^
#            49 + 185*u^47 + 173*u^45 + 146*u^43 + 130*u^41 + 114*u^39 + 89*u^
#            37 + 75*u^35 + 61*u^33 + 43*u^31 + 34*u^29 + 25*u^27 + 15*u^25 + 
#            11*u^23 + 7*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          2*u^75 + 4*u^73 + 9*u^71 + 19*u^69 + 29*u^67 + 44*u^65 + 66*u^63 + 
#            83*u^61 + 105*u^59 + 130*u^57 + 143*u^55 + 158*u^53 + 171*u^51 + 
#            167*u^49 + 166*u^47 + 161*u^45 + 141*u^43 + 128*u^41 + 113*u^39 + 
#            89*u^37 + 75*u^35 + 61*u^33 + 43*u^31 + 34*u^29 + 25*u^27 + 15*u^
#            25 + 11*u^23 + 7*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          u^63 + 3*u^61 + 9*u^59 + 19*u^57 + 31*u^55 + 50*u^53 + 70*u^51 + 
#            85*u^49 + 102*u^47 + 112*u^45 + 110*u^43 + 109*u^41 + 101*u^39 + 
#            84*u^37 + 73*u^35 + 60*u^33 + 43*u^31 + 34*u^29 + 25*u^27 + 15*u^
#            25 + 11*u^23 + 7*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          u^53 + 5*u^51 + 11*u^49 + 19*u^47 + 34*u^45 + 44*u^43 + 54*u^41 + 
#            65*u^39 + 62*u^37 + 60*u^35 + 55*u^33 + 41*u^31 + 33*u^29 + 25*u^
#            27 + 15*u^25 + 11*u^23 + 7*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          u^54 + 2*u^52 + 5*u^50 + 7*u^48 + 9*u^46 + 10*u^44 + 9*u^42 + 8*u^
#            40 + 5*u^38 + 3*u^36 + u^34, u^52 + 2*u^50 + 8*u^48 + 14*u^46 + 
#            19*u^44 + 30*u^42 + 32*u^40 + 32*u^38 + 33*u^36 + 24*u^34 + 17*u^
#            32 + 12*u^30 + 5*u^28 + 2*u^26 + u^24, 
#          u^45 + 4*u^43 + 12*u^41 + 25*u^39 + 33*u^37 + 41*u^35 + 43*u^33 + 
#            36*u^31 + 31*u^29 + 24*u^27 + 15*u^25 + 11*u^23 + 7*u^21 + 3*u^
#            19 + 2*u^17 + u^15, u^44 + 7*u^42 + 13*u^40 + 20*u^38 + 26*u^36 + 
#            22*u^34 + 17*u^32 + 12*u^30 + 5*u^28 + 2*u^26 + u^24, 
#          u^43 + 4*u^41 + 6*u^39 + 7*u^37 + 5*u^35 + 3*u^33 + u^31, 
#          u^41 + 7*u^39 + 16*u^37 + 29*u^35 + 36*u^33 + 34*u^31 + 31*u^29 + 
#            24*u^27 + 15*u^25 + 11*u^23 + 7*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          2*u^39 + 7*u^37 + 15*u^35 + 24*u^33 + 31*u^31 + 34*u^29 + 31*u^27 + 
#            24*u^25 + 16*u^23 + 9*u^21 + 4*u^19 + 2*u^17 + u^15, 
#          u^37 + 4*u^35 + 12*u^33 + 15*u^31 + 18*u^29 + 19*u^27 + 13*u^25 + 
#            10*u^23 + 7*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          u^34 + 2*u^32 + 2*u^30 + u^28, u^37 + 2*u^35 + 8*u^33 + 10*u^31 + 
#            10*u^29 + 10*u^27 + 5*u^25 + 2*u^23 + u^21, 
#          u^33 + 3*u^31 + 7*u^29 + 12*u^27 + 11*u^25 + 10*u^23 + 7*u^21 + 3*u^
#            19 + 2*u^17 + u^15, 2*u^32 + 5*u^30 + 6*u^28 + 5*u^26 + 2*u^24, 
#          2*u^31 + 5*u^29 + 7*u^27 + 9*u^25 + 5*u^23 + 2*u^21 + u^19, u^31, 
#          u^27, 2*u^27 + 2*u^25 + 2*u^23 + 4*u^21 + 2*u^19 + 2*u^17 + 2*u^15, 
#          3*u^27 + 6*u^25 + 8*u^23 + 6*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          3*u^27 + 7*u^25 + 5*u^23 + 2*u^21 + u^19, 3*u^27 + 4*u^25 + 2*u^23, 
#          0*u^0, u^27, 4*u^27 + 11*u^25 + 12*u^23 + 8*u^21 + 4*u^19 + 2*u^
#            17 + u^15, u^27 + 4*u^25 + 7*u^23 + 7*u^21 + 4*u^19 + 2*u^17 + u^
#            15, u^27 + 2*u^25 + u^23, u^26 + 4*u^24 + 5*u^22 + 3*u^20 + u^18, 
#          u^25 + 5*u^23 + 6*u^21 + 3*u^19 + 2*u^17 + u^15, u^24 + u^22, 
#          3*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          2*u^23 + 5*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          u^24 + 3*u^22 + 3*u^20 + u^18, 2*u^21 + 2*u^19 + 2*u^17 + 2*u^15, 
#          u^20 + u^18, u^20 + u^18, 3*u^21 + 3*u^19 + 2*u^17 + u^15, u^20, 
#          2*u^21 + u^19, 2*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          2*u^19 + 2*u^17 + 2*u^15, u^19 + u^17 + u^15, u^19 + u^17 + 2*u^15, 
#          u^19 + u^17 + 3*u^15, u^20 + u^18 + u^16, u^15, 0*u^0, u^15, 0*u^0, 
#          0*u^0, u^15, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^98 + 2*u^96 + 3*u^94 + 7*u^92 + 10*u^90 + 14*u^88 + 23*u^86 + 29*u^
#            84 + 38*u^82 + 53*u^80 + 62*u^78 + 76*u^76 + 95*u^74 + 104*u^72 + 
#            121*u^70 + 139*u^68 + 144*u^66 + 160*u^64 + 171*u^62 + 169*u^60 + 
#            179*u^58 + 179*u^56 + 169*u^54 + 171*u^52 + 160*u^50 + 144*u^48 + 
#            139*u^46 + 121*u^44 + 104*u^42 + 95*u^40 + 76*u^38 + 62*u^36 + 
#            53*u^34 + 38*u^32 + 29*u^30 + 23*u^28 + 14*u^26 + 10*u^24 + 7*u^
#            22 + 3*u^20 + 2*u^18 + u^16, u^76 + 3*u^74 + 5*u^72 + 11*u^70 + 
#            19*u^68 + 27*u^66 + 42*u^64 + 57*u^62 + 70*u^60 + 90*u^58 + 104*u^
#            56 + 113*u^54 + 127*u^52 + 129*u^50 + 126*u^48 + 127*u^46 + 115*u^
#            44 + 102*u^42 + 94*u^40 + 76*u^38 + 62*u^36 + 53*u^34 + 38*u^32 + 
#            29*u^30 + 23*u^28 + 14*u^26 + 10*u^24 + 7*u^22 + 3*u^20 + 2*u^
#            18 + u^16, u^62 + 4*u^60 + 10*u^58 + 18*u^56 + 30*u^54 + 46*u^
#            52 + 58*u^50 + 71*u^48 + 83*u^46 + 84*u^44 + 84*u^42 + 82*u^40 + 
#            70*u^38 + 60*u^36 + 52*u^34 + 38*u^32 + 29*u^30 + 23*u^28 + 14*u^
#            26 + 10*u^24 + 7*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          u^54 + 3*u^52 + 7*u^50 + 13*u^48 + 21*u^46 + 31*u^44 + 38*u^42 + 
#            46*u^40 + 49*u^38 + 46*u^36 + 45*u^34 + 36*u^32 + 28*u^30 + 23*u^
#            28 + 14*u^26 + 10*u^24 + 7*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          2*u^51 + 3*u^49 + 4*u^47 + 7*u^45 + 6*u^43 + 6*u^41 + 6*u^39 + 3*u^
#            37 + 2*u^35 + u^33, u^53 + 2*u^51 + 4*u^49 + 10*u^47 + 13*u^45 + 
#            18*u^43 + 25*u^41 + 23*u^39 + 24*u^37 + 22*u^35 + 14*u^33 + 11*u^
#            31 + 6*u^29 + 2*u^27 + u^25, 2*u^44 + 6*u^42 + 14*u^40 + 23*u^
#            38 + 29*u^36 + 33*u^34 + 30*u^32 + 26*u^30 + 22*u^28 + 14*u^26 + 
#            10*u^24 + 7*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          2*u^43 + 7*u^41 + 11*u^39 + 17*u^37 + 18*u^35 + 13*u^33 + 11*u^31 + 
#            6*u^29 + 2*u^27 + u^25, u^44 + 3*u^42 + 4*u^40 + 6*u^38 + 6*u^
#            36 + 3*u^34 + 2*u^32 + u^30, 3*u^40 + 8*u^38 + 17*u^36 + 26*u^
#            34 + 26*u^32 + 25*u^30 + 22*u^28 + 14*u^26 + 10*u^24 + 7*u^22 + 
#            3*u^20 + 2*u^18 + u^16, u^40 + 4*u^38 + 10*u^36 + 17*u^34 + 23*u^
#            32 + 27*u^30 + 26*u^28 + 21*u^26 + 15*u^24 + 9*u^22 + 4*u^20 + 
#            2*u^18 + u^16, u^36 + 4*u^34 + 10*u^32 + 12*u^30 + 15*u^28 + 12*u^
#            26 + 9*u^24 + 7*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          u^37 + 2*u^35 + 4*u^33 + 3*u^31 + 2*u^29 + u^27, 
#          u^34 + 6*u^32 + 6*u^30 + 7*u^28 + 5*u^26 + 2*u^24 + u^22, 
#          2*u^32 + 4*u^30 + 9*u^28 + 9*u^26 + 8*u^24 + 7*u^22 + 3*u^20 + 2*u^
#            18 + u^16, 2*u^31 + 4*u^29 + 3*u^27 + 2*u^25 + u^23, 
#          u^32 + 4*u^30 + 5*u^28 + 7*u^26 + 5*u^24 + 2*u^22 + u^20, 0*u^0, 
#          u^28, u^32 + 2*u^30 + 3*u^28 + 4*u^26 + 4*u^24 + 4*u^22 + 3*u^20 + 
#            2*u^18 + u^16, 2*u^28 + 4*u^26 + 6*u^24 + 6*u^22 + 3*u^20 + 2*u^
#            18 + u^16, u^28 + 4*u^26 + 4*u^24 + 2*u^22 + u^20, 
#          2*u^28 + 3*u^26 + 2*u^24 + u^22, 0*u^0, 0*u^0, 
#          u^28 + 4*u^26 + 7*u^24 + 7*u^22 + 4*u^20 + 2*u^18 + u^16, 
#          u^26 + 4*u^24 + 6*u^22 + 4*u^20 + 2*u^18 + u^16, 0*u^0, 
#          2*u^25 + 4*u^23 + 4*u^21 + 2*u^19, 
#          u^24 + 4*u^22 + 3*u^20 + 2*u^18 + u^16, u^25 + 2*u^23 + u^21, 
#          u^23 + u^21 + u^19, 4*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          u^23 + 3*u^21 + 2*u^19, 2*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          u^21 + u^19, u^23 + 2*u^21 + u^19, u^22 + 2*u^20 + 2*u^18 + u^16, 
#          2*u^21 + u^19, 0*u^0, 2*u^20 + 2*u^18 + u^16, 
#          2*u^20 + 2*u^18 + u^16, u^20 + u^18 + u^16, u^18 + u^16, 
#          u^20 + u^18 + u^16, u^15, 0*u^0, u^18 + u^16, 0*u^0, u^15, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^97 + u^95 + 3*u^93 + 6*u^91 + 7*u^89 + 13*u^87 + 19*u^85 + 22*u^
#            83 + 34*u^81 + 43*u^79 + 48*u^77 + 66*u^75 + 75*u^73 + 81*u^71 + 
#            102*u^69 + 107*u^67 + 112*u^65 + 131*u^63 + 128*u^61 + 130*u^59 + 
#            142*u^57 + 130*u^55 + 128*u^53 + 131*u^51 + 112*u^49 + 107*u^47 + 
#            102*u^45 + 81*u^43 + 75*u^41 + 66*u^39 + 48*u^37 + 43*u^35 + 34*u^
#            33 + 22*u^31 + 19*u^29 + 13*u^27 + 7*u^25 + 6*u^23 + 3*u^21 + u^
#            19 + u^17, 2*u^75 + 3*u^73 + 6*u^71 + 13*u^69 + 18*u^67 + 26*u^
#            65 + 41*u^63 + 48*u^61 + 61*u^59 + 77*u^57 + 81*u^55 + 90*u^53 + 
#            100*u^51 + 93*u^49 + 95*u^47 + 93*u^45 + 78*u^43 + 73*u^41 + 65*u^
#            39 + 48*u^37 + 43*u^35 + 34*u^33 + 22*u^31 + 19*u^29 + 13*u^27 + 
#            7*u^25 + 6*u^23 + 3*u^21 + u^19 + u^17, 
#          u^63 + 2*u^61 + 6*u^59 + 12*u^57 + 18*u^55 + 29*u^53 + 40*u^51 + 
#            46*u^49 + 57*u^47 + 62*u^45 + 59*u^43 + 61*u^41 + 56*u^39 + 45*u^
#            37 + 41*u^35 + 33*u^33 + 22*u^31 + 19*u^29 + 13*u^27 + 7*u^25 + 
#            6*u^23 + 3*u^21 + u^19 + u^17, 2*u^53 + 4*u^51 + 9*u^49 + 13*u^
#            47 + 22*u^45 + 26*u^43 + 31*u^41 + 37*u^39 + 33*u^37 + 32*u^35 + 
#            30*u^33 + 20*u^31 + 18*u^29 + 13*u^27 + 7*u^25 + 6*u^23 + 3*u^
#            21 + u^19 + u^17, u^50 + u^46 + u^44 + u^40, 
#          u^54 + 3*u^52 + 4*u^50 + 9*u^48 + 14*u^46 + 14*u^44 + 23*u^42 + 
#            22*u^40 + 20*u^38 + 22*u^36 + 15*u^34 + 10*u^32 + 9*u^30 + 3*u^
#            28 + 2*u^26 + u^24, u^45 + 3*u^43 + 6*u^41 + 13*u^39 + 16*u^37 + 
#            20*u^35 + 21*u^33 + 17*u^31 + 16*u^29 + 12*u^27 + 7*u^25 + 6*u^
#            23 + 3*u^21 + u^19 + u^17, u^46 + 2*u^44 + 7*u^42 + 9*u^40 + 13*u^
#            38 + 17*u^36 + 13*u^34 + 10*u^32 + 9*u^30 + 3*u^28 + 2*u^26 + u^24
#            , u^37, u^41 + 4*u^39 + 7*u^37 + 14*u^35 + 17*u^33 + 15*u^31 + 
#            16*u^29 + 12*u^27 + 7*u^25 + 6*u^23 + 3*u^21 + u^19 + u^17, 
#          u^41 + 3*u^39 + 8*u^37 + 13*u^35 + 19*u^33 + 22*u^31 + 23*u^29 + 
#            19*u^27 + 15*u^25 + 9*u^23 + 5*u^21 + 2*u^19 + u^17, 
#          u^35 + 5*u^33 + 6*u^31 + 7*u^29 + 9*u^27 + 5*u^25 + 5*u^23 + 3*u^
#            21 + u^19 + u^17, 0*u^0, u^37 + u^35 + 5*u^33 + 6*u^31 + 5*u^29 + 
#            7*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          u^33 + 2*u^31 + 4*u^29 + 7*u^27 + 5*u^25 + 5*u^23 + 3*u^21 + u^
#            19 + u^17, u^30, u^33 + 4*u^31 + 6*u^29 + 6*u^27 + 8*u^25 + 3*u^
#            23 + 2*u^21 + u^19, 0*u^0, u^33 + u^31 + 2*u^29 + 3*u^27 + u^
#            25 + u^23, u^31 + u^29 + 3*u^27 + 3*u^25 + 3*u^23 + 4*u^21 + 2*u^
#            19 + 2*u^17 + u^15, u^27 + 2*u^25 + 3*u^23 + 2*u^21 + u^19 + u^17,
#          u^29 + 2*u^27 + 5*u^25 + 3*u^23 + 2*u^21 + u^19, 0*u^0, 0*u^0, 
#          u^29 + 2*u^27 + u^25 + u^23, u^27 + 4*u^25 + 4*u^23 + 4*u^21 + 2*u^
#            19 + u^17, u^27 + 3*u^25 + 4*u^23 + 4*u^21 + 2*u^19 + u^17, 
#          0*u^0, u^28 + 2*u^26 + 4*u^24 + 4*u^22 + 3*u^20 + u^18, 
#          u^23 + 2*u^21 + u^19 + u^17, 0*u^0, 2*u^24 + u^22 + 2*u^20 + u^18, 
#          u^23 + 2*u^21 + u^19 + u^17, u^24 + u^22 + 3*u^20 + u^18, 
#          2*u^21 + 2*u^19 + 2*u^17 + u^15, u^24 + u^22 + 3*u^20 + 2*u^18, 
#          0*u^0, u^19 + u^17, u^20, 0*u^0, u^19 + u^17, 
#          2*u^19 + 2*u^17 + u^15, 0*u^0, u^19 + 2*u^17 + u^15, u^15, 0*u^0, 
#          u^19 + u^17 + 2*u^15, u^17, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^98 + 2*u^96 + 5*u^94 + 8*u^92 + 15*u^90 + 22*u^88 + 34*u^86 + 46*u^
#            84 + 65*u^82 + 83*u^80 + 108*u^78 + 131*u^76 + 161*u^74 + 187*u^
#            72 + 218*u^70 + 243*u^68 + 272*u^66 + 292*u^64 + 314*u^62 + 325*u^
#            60 + 337*u^58 + 337*u^56 + 337*u^54 + 325*u^52 + 314*u^50 + 292*u^
#            48 + 272*u^46 + 243*u^44 + 218*u^42 + 187*u^40 + 161*u^38 + 131*u^
#            36 + 108*u^34 + 83*u^32 + 65*u^30 + 46*u^28 + 34*u^26 + 22*u^24 + 
#            15*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^76 + 3*u^74 + 8*u^72 + 16*u^70 + 29*u^68 + 47*u^66 + 70*u^64 + 
#            98*u^62 + 128*u^60 + 160*u^58 + 189*u^56 + 216*u^54 + 234*u^52 + 
#            247*u^50 + 248*u^48 + 244*u^46 + 228*u^44 + 210*u^42 + 184*u^40 + 
#            160*u^38 + 131*u^36 + 108*u^34 + 83*u^32 + 65*u^30 + 46*u^28 + 
#            34*u^26 + 22*u^24 + 15*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          2*u^62 + 6*u^60 + 16*u^58 + 30*u^56 + 53*u^54 + 77*u^52 + 107*u^
#            50 + 130*u^48 + 153*u^46 + 161*u^44 + 166*u^42 + 156*u^40 + 145*u^
#            38 + 123*u^36 + 105*u^34 + 82*u^32 + 65*u^30 + 46*u^28 + 34*u^
#            26 + 22*u^24 + 15*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^54 + 4*u^52 + 11*u^50 + 22*u^48 + 38*u^46 + 55*u^44 + 74*u^42 + 
#            86*u^40 + 96*u^38 + 93*u^36 + 89*u^34 + 74*u^32 + 62*u^30 + 45*u^
#            28 + 34*u^26 + 22*u^24 + 15*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14
#            , u^53 + 2*u^51 + 5*u^49 + 6*u^47 + 9*u^45 + 8*u^43 + 9*u^41 + 
#            6*u^39 + 5*u^37 + 2*u^35 + u^33, 
#          u^53 + 3*u^51 + 9*u^49 + 16*u^47 + 28*u^45 + 37*u^43 + 48*u^41 + 
#            51*u^39 + 52*u^37 + 45*u^35 + 36*u^33 + 25*u^31 + 15*u^29 + 8*u^
#            27 + 3*u^25 + u^23, 2*u^44 + 10*u^42 + 23*u^40 + 40*u^38 + 52*u^
#            36 + 61*u^34 + 59*u^32 + 54*u^30 + 42*u^28 + 33*u^26 + 22*u^24 + 
#            15*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^45 + 5*u^43 + 14*u^41 + 25*u^39 + 34*u^37 + 37*u^35 + 33*u^33 + 
#            25*u^31 + 15*u^29 + 8*u^27 + 3*u^25 + u^23, 
#          2*u^42 + 4*u^40 + 7*u^38 + 6*u^36 + 5*u^34 + 2*u^32 + u^30, 
#          4*u^40 + 14*u^38 + 29*u^36 + 44*u^34 + 51*u^32 + 51*u^30 + 42*u^
#            28 + 33*u^26 + 22*u^24 + 15*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14
#            , 2*u^40 + 8*u^38 + 18*u^36 + 33*u^34 + 46*u^32 + 56*u^30 + 55*u^
#            28 + 49*u^26 + 35*u^24 + 23*u^22 + 11*u^20 + 6*u^18 + 2*u^16 + u^
#            14, 2*u^36 + 8*u^34 + 16*u^32 + 24*u^30 + 26*u^28 + 25*u^26 + 
#            19*u^24 + 14*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^35 + 2*u^33 + 3*u^31 + 2*u^29 + u^27, 
#          2*u^36 + 5*u^34 + 12*u^32 + 15*u^30 + 17*u^28 + 12*u^26 + 8*u^24 + 
#            3*u^22 + u^20, 3*u^32 + 8*u^30 + 15*u^28 + 19*u^26 + 18*u^24 + 
#            14*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^33 + 3*u^31 + 5*u^29 + 5*u^27 + 3*u^25 + u^23, 
#          2*u^32 + 7*u^30 + 13*u^28 + 16*u^26 + 13*u^24 + 8*u^22 + 3*u^20 + u^
#            18, 0*u^0, u^30 + 2*u^28 + 2*u^26 + u^24, 
#          u^32 + 2*u^30 + 4*u^28 + 6*u^26 + 7*u^24 + 8*u^22 + 6*u^20 + 6*u^
#            18 + 3*u^16 + 2*u^14, u^28 + 7*u^26 + 10*u^24 + 11*u^22 + 7*u^
#            20 + 5*u^18 + 2*u^16 + u^14, 2*u^28 + 9*u^26 + 10*u^24 + 8*u^22 + 
#            3*u^20 + u^18, u^28 + 3*u^26 + 3*u^24 + u^22, 0*u^0, 
#          u^28 + 2*u^26 + u^24, u^28 + 8*u^26 + 15*u^24 + 16*u^22 + 10*u^20 + 
#            6*u^18 + 2*u^16 + u^14, 4*u^26 + 9*u^24 + 13*u^22 + 10*u^20 + 6*u^
#            18 + 2*u^16 + u^14, u^24, u^27 + 5*u^25 + 9*u^23 + 9*u^21 + 5*u^
#            19 + u^17, 3*u^24 + 7*u^22 + 7*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^23 + u^21, u^25 + 4*u^23 + 5*u^21 + 3*u^19 + u^17, 
#          u^24 + 5*u^22 + 7*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          3*u^23 + 6*u^21 + 5*u^19 + u^17, 
#          2*u^22 + 5*u^20 + 6*u^18 + 3*u^16 + 2*u^14, 
#          u^23 + 3*u^21 + 4*u^19 + u^17, u^21 + 2*u^19, 
#          u^22 + 4*u^20 + 5*u^18 + 2*u^16 + u^14, u^21 + 2*u^19, u^20, 
#          3*u^20 + 5*u^18 + 2*u^16 + u^14, 2*u^20 + 6*u^18 + 3*u^16 + 2*u^14, 
#          2*u^18 + u^16 + u^14, u^20 + 3*u^18 + 3*u^16 + 2*u^14, 
#          2*u^18 + 2*u^16 + 2*u^14, u^19 + u^17 + u^15, u^18 + u^16 + u^14, 
#          u^18 + u^16 + u^14, u^14, u^15, 0*u^0, 0*u^0, u^14, 0*u^0, u^15, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^14, u^14, 0*u^0, u^14, u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^98 + u^96 + 3*u^94 + 4*u^92 + 8*u^90 + 11*u^88 + 17*u^86 + 22*u^
#            84 + 32*u^82 + 39*u^80 + 52*u^78 + 61*u^76 + 76*u^74 + 87*u^72 + 
#            102*u^70 + 112*u^68 + 127*u^66 + 134*u^64 + 146*u^62 + 149*u^60 + 
#            156*u^58 + 155*u^56 + 156*u^54 + 149*u^52 + 146*u^50 + 134*u^48 + 
#            127*u^46 + 112*u^44 + 102*u^42 + 87*u^40 + 76*u^38 + 61*u^36 + 
#            52*u^34 + 39*u^32 + 32*u^30 + 22*u^28 + 17*u^26 + 11*u^24 + 8*u^
#            22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^74 + 3*u^72 + 7*u^70 + 12*u^68 + 21*u^66 + 31*u^64 + 45*u^62 + 
#            58*u^60 + 74*u^58 + 87*u^56 + 101*u^54 + 108*u^52 + 116*u^50 + 
#            115*u^48 + 115*u^46 + 106*u^44 + 99*u^42 + 86*u^40 + 76*u^38 + 
#            61*u^36 + 52*u^34 + 39*u^32 + 32*u^30 + 22*u^28 + 17*u^26 + 11*u^
#            24 + 8*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^62 + 2*u^60 + 7*u^58 + 13*u^56 + 25*u^54 + 35*u^52 + 51*u^50 + 
#            61*u^48 + 74*u^46 + 76*u^44 + 80*u^42 + 74*u^40 + 70*u^38 + 58*u^
#            36 + 51*u^34 + 39*u^32 + 32*u^30 + 22*u^28 + 17*u^26 + 11*u^24 + 
#            8*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^52 + 3*u^50 + 8*u^48 + 15*u^46 + 23*u^44 + 33*u^42 + 39*u^40 + 
#            46*u^38 + 44*u^36 + 44*u^34 + 36*u^32 + 31*u^30 + 22*u^28 + 17*u^
#            26 + 11*u^24 + 8*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^53 + 2*u^51 + 5*u^49 + 6*u^47 + 9*u^45 + 8*u^43 + 9*u^41 + 6*u^
#            39 + 5*u^37 + 2*u^35 + u^33, 2*u^49 + 4*u^47 + 10*u^45 + 13*u^
#            43 + 20*u^41 + 21*u^39 + 23*u^37 + 20*u^35 + 16*u^33 + 11*u^31 + 
#            6*u^29 + 3*u^27 + u^25, 4*u^42 + 10*u^40 + 20*u^38 + 26*u^36 + 
#            32*u^34 + 30*u^32 + 28*u^30 + 21*u^28 + 17*u^26 + 11*u^24 + 8*u^
#            22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^43 + 5*u^41 + 10*u^39 + 15*u^37 + 17*u^35 + 15*u^33 + 11*u^31 + 
#            6*u^29 + 3*u^27 + u^25, 2*u^42 + 4*u^40 + 7*u^38 + 6*u^36 + 5*u^
#            34 + 2*u^32 + u^30, u^40 + 7*u^38 + 15*u^36 + 24*u^34 + 27*u^32 + 
#            27*u^30 + 21*u^28 + 17*u^26 + 11*u^24 + 8*u^22 + 4*u^20 + 3*u^
#            18 + u^16 + u^14, 2*u^38 + 5*u^36 + 11*u^34 + 17*u^32 + 23*u^30 + 
#            23*u^28 + 22*u^26 + 16*u^24 + 11*u^22 + 5*u^20 + 3*u^18 + u^
#            16 + u^14, u^36 + 5*u^34 + 9*u^32 + 14*u^30 + 14*u^28 + 14*u^26 + 
#            10*u^24 + 8*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^35 + 2*u^33 + 3*u^31 + 2*u^29 + u^27, 
#          u^36 + 2*u^34 + 6*u^32 + 7*u^30 + 8*u^28 + 5*u^26 + 3*u^24 + u^22, 
#          u^32 + 3*u^30 + 6*u^28 + 9*u^26 + 9*u^24 + 8*u^22 + 4*u^20 + 3*u^
#            18 + u^16 + u^14, 2*u^31 + 5*u^29 + 6*u^27 + 4*u^25 + u^23, 
#          u^30 + 3*u^28 + 5*u^26 + 5*u^24 + 3*u^22 + u^20, u^32 + u^30 + u^28,
#          0*u^0, u^26 + u^24 + 2*u^22 + u^20 + 2*u^18 + u^16 + u^14, 
#          4*u^26 + 6*u^24 + 7*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          3*u^26 + 4*u^24 + 3*u^22 + u^20, u^28 + 4*u^26 + 4*u^24 + u^22, 
#          u^26, 0*u^0, 5*u^26 + 10*u^24 + 9*u^22 + 5*u^20 + 3*u^18 + u^16 + u^
#            14, u^26 + 3*u^24 + 6*u^22 + 5*u^20 + 3*u^18 + u^16 + u^14, 
#          u^26 + 2*u^24, u^25 + 3*u^23 + 4*u^21 + 2*u^19, 
#          2*u^24 + 5*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^25 + 2*u^23 + u^21, 2*u^23 + 2*u^21 + u^19, 
#          3*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14, u^23 + 3*u^21 + 2*u^19, 
#          u^20 + 2*u^18 + u^16 + u^14, u^19, u^19, 
#          u^22 + 3*u^20 + 3*u^18 + u^16 + u^14, u^21 + u^19, 2*u^20, 
#          3*u^20 + 3*u^18 + u^16 + u^14, 2*u^18 + u^16 + u^14, 
#          u^20 + 2*u^18 + u^16 + u^14, u^18 + u^16 + u^14, u^18 + u^16 + u^14,
#          u^19 + u^17 + u^15, 0*u^0, 0*u^0, u^14, u^15, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^97 + 2*u^95 + 4*u^93 + 7*u^91 + 11*u^89 + 18*u^87 + 26*u^85 + 36*u^
#            83 + 49*u^81 + 63*u^79 + 81*u^77 + 100*u^75 + 119*u^73 + 141*u^
#            71 + 161*u^69 + 182*u^67 + 202*u^65 + 217*u^63 + 232*u^61 + 242*u^
#            59 + 248*u^57 + 252*u^55 + 248*u^53 + 242*u^51 + 232*u^49 + 217*u^
#            47 + 202*u^45 + 182*u^43 + 161*u^41 + 141*u^39 + 119*u^37 + 100*u^
#            35 + 81*u^33 + 63*u^31 + 49*u^29 + 36*u^27 + 26*u^25 + 18*u^23 + 
#            11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          u^75 + 2*u^73 + 6*u^71 + 12*u^69 + 22*u^67 + 36*u^65 + 53*u^63 + 
#            74*u^61 + 98*u^59 + 121*u^57 + 145*u^55 + 163*u^53 + 178*u^51 + 
#            187*u^49 + 188*u^47 + 184*u^45 + 173*u^43 + 157*u^41 + 140*u^39 + 
#            119*u^37 + 100*u^35 + 81*u^33 + 63*u^31 + 49*u^29 + 36*u^27 + 
#            26*u^25 + 18*u^23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          u^61 + 5*u^59 + 12*u^57 + 25*u^55 + 42*u^53 + 63*u^51 + 85*u^49 + 
#            104*u^47 + 120*u^45 + 128*u^43 + 128*u^41 + 122*u^39 + 110*u^37 + 
#            96*u^35 + 80*u^33 + 63*u^31 + 49*u^29 + 36*u^27 + 26*u^25 + 18*u^
#            23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          u^53 + 3*u^51 + 9*u^49 + 17*u^47 + 30*u^45 + 44*u^43 + 58*u^41 + 
#            69*u^39 + 75*u^37 + 75*u^35 + 70*u^33 + 59*u^31 + 48*u^29 + 36*u^
#            27 + 26*u^25 + 18*u^23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13
#            , u^52 + 3*u^50 + 6*u^48 + 9*u^46 + 11*u^44 + 12*u^42 + 11*u^40 + 
#            9*u^38 + 6*u^36 + 3*u^34 + u^32, 2*u^50 + 5*u^48 + 11*u^46 + 18*u^
#            44 + 26*u^42 + 33*u^40 + 37*u^38 + 36*u^36 + 32*u^34 + 24*u^32 + 
#            17*u^30 + 9*u^28 + 4*u^26 + u^24, 
#          2*u^43 + 9*u^41 + 22*u^39 + 36*u^37 + 47*u^35 + 52*u^33 + 50*u^31 + 
#            44*u^29 + 35*u^27 + 26*u^25 + 18*u^23 + 11*u^21 + 7*u^19 + 4*u^
#            17 + 2*u^15 + u^13, 3*u^42 + 10*u^40 + 19*u^38 + 25*u^36 + 27*u^
#            34 + 23*u^32 + 17*u^30 + 9*u^28 + 4*u^26 + u^24, 
#          u^43 + 4*u^41 + 8*u^39 + 10*u^37 + 9*u^35 + 6*u^33 + 3*u^31 + u^29, 
#          4*u^39 + 15*u^37 + 29*u^35 + 41*u^33 + 45*u^31 + 43*u^29 + 35*u^
#            27 + 26*u^25 + 18*u^23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13
#            , u^39 + 6*u^37 + 14*u^35 + 25*u^33 + 35*u^31 + 43*u^29 + 42*u^
#            27 + 37*u^25 + 26*u^23 + 15*u^21 + 8*u^19 + 4*u^17 + 2*u^15 + u^13
#            , 3*u^35 + 9*u^33 + 17*u^31 + 23*u^29 + 25*u^27 + 22*u^25 + 17*u^
#            23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          u^36 + 3*u^34 + 5*u^32 + 5*u^30 + 3*u^28 + u^26, 
#          u^35 + 4*u^33 + 9*u^31 + 12*u^29 + 12*u^27 + 8*u^25 + 4*u^23 + u^21,
#          3*u^31 + 8*u^29 + 14*u^27 + 17*u^25 + 16*u^23 + 11*u^21 + 7*u^19 + 
#            4*u^17 + 2*u^15 + u^13, u^32 + 4*u^30 + 7*u^28 + 7*u^26 + 4*u^
#            24 + u^22, u^31 + 5*u^29 + 8*u^27 + 11*u^25 + 8*u^23 + 4*u^21 + u^
#            19, u^31 + u^29 + u^27, 0*u^0, u^31 + 2*u^29 + 4*u^27 + 5*u^25 + 
#            6*u^23 + 6*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          3*u^27 + 9*u^25 + 12*u^23 + 10*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^
#            13, 2*u^27 + 7*u^25 + 7*u^23 + 4*u^21 + u^19, 
#          3*u^27 + 6*u^25 + 4*u^23 + u^21, u^27 + u^25, 0*u^0, 
#          2*u^27 + 11*u^25 + 16*u^23 + 13*u^21 + 8*u^19 + 4*u^17 + 2*u^15 + u^
#            13, 3*u^25 + 8*u^23 + 11*u^21 + 8*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          u^25 + u^23, 3*u^24 + 7*u^22 + 7*u^20 + 2*u^18, 
#          u^25 + 5*u^23 + 8*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          2*u^24 + 3*u^22 + u^20, u^24 + 3*u^22 + 3*u^20 + u^18, 
#          u^23 + 7*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          3*u^22 + 6*u^20 + 2*u^18, 3*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          2*u^20 + u^18, u^22 + 3*u^20 + 2*u^18, 
#          3*u^21 + 6*u^19 + 4*u^17 + 2*u^15 + u^13, 3*u^20 + u^18, 
#          u^21 + u^19, u^21 + 6*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          4*u^19 + 4*u^17 + 2*u^15 + u^13, 3*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^19 + 2*u^17 + 2*u^15 + u^13, 2*u^19 + 3*u^17 + 3*u^15 + u^13, 
#          u^18 + u^16 + u^14, 0*u^0, u^17 + u^15 + u^13, u^15 + u^13, u^14, 
#          u^15, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, u^13, u^14, 0*u^0, u^13, 
#          u^13, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^95 + u^93 + 2*u^91 + 6*u^89 + 7*u^87 + 11*u^85 + 20*u^83 + 23*u^
#            81 + 32*u^79 + 47*u^77 + 51*u^75 + 66*u^73 + 85*u^71 + 88*u^69 + 
#            107*u^67 + 125*u^65 + 124*u^63 + 143*u^61 + 154*u^59 + 146*u^57 + 
#            161*u^55 + 161*u^53 + 146*u^51 + 154*u^49 + 143*u^47 + 124*u^45 + 
#            125*u^43 + 107*u^41 + 88*u^39 + 85*u^37 + 66*u^35 + 51*u^33 + 
#            47*u^31 + 32*u^29 + 23*u^27 + 20*u^25 + 11*u^23 + 7*u^21 + 6*u^
#            19 + 2*u^17 + u^15 + u^13, u^73 + 4*u^71 + 5*u^69 + 13*u^67 + 
#            21*u^65 + 28*u^63 + 44*u^61 + 59*u^59 + 67*u^57 + 89*u^55 + 100*u^
#            53 + 104*u^51 + 119*u^49 + 119*u^47 + 111*u^45 + 116*u^43 + 102*u^
#            41 + 87*u^39 + 84*u^37 + 66*u^35 + 51*u^33 + 47*u^31 + 32*u^29 + 
#            23*u^27 + 20*u^25 + 11*u^23 + 7*u^21 + 6*u^19 + 2*u^17 + u^15 + u^
#            13, u^61 + 3*u^59 + 6*u^57 + 15*u^55 + 23*u^53 + 34*u^51 + 51*u^
#            49 + 60*u^47 + 69*u^45 + 81*u^43 + 78*u^41 + 74*u^39 + 75*u^37 + 
#            61*u^35 + 50*u^33 + 46*u^31 + 32*u^29 + 23*u^27 + 20*u^25 + 11*u^
#            23 + 7*u^21 + 6*u^19 + 2*u^17 + u^15 + u^13, 
#          u^53 + 2*u^51 + 6*u^49 + 13*u^47 + 18*u^45 + 29*u^43 + 39*u^41 + 
#            40*u^39 + 50*u^37 + 48*u^35 + 41*u^33 + 41*u^31 + 31*u^29 + 22*u^
#            27 + 20*u^25 + 11*u^23 + 7*u^21 + 6*u^19 + 2*u^17 + u^15 + u^13, 
#          u^48 + u^46 + u^44 + 2*u^42 + u^40 + u^38 + u^36, 
#          3*u^50 + 4*u^48 + 8*u^46 + 17*u^44 + 17*u^42 + 24*u^40 + 29*u^38 + 
#            23*u^36 + 24*u^34 + 20*u^32 + 11*u^30 + 9*u^28 + 5*u^26 + u^
#            24 + u^22, 2*u^43 + 6*u^41 + 11*u^39 + 21*u^37 + 26*u^35 + 28*u^
#            33 + 32*u^31 + 26*u^29 + 21*u^27 + 19*u^25 + 11*u^23 + 7*u^21 + 
#            6*u^19 + 2*u^17 + u^15 + u^13, u^44 + 2*u^42 + 8*u^40 + 14*u^38 + 
#            15*u^36 + 20*u^34 + 18*u^32 + 11*u^30 + 9*u^28 + 5*u^26 + u^
#            24 + u^22, u^39 + u^37 + u^35 + u^33, 
#          2*u^39 + 8*u^37 + 14*u^35 + 20*u^33 + 28*u^31 + 24*u^29 + 21*u^27 + 
#            19*u^25 + 11*u^23 + 7*u^21 + 6*u^19 + 2*u^17 + u^15 + u^13, 
#          2*u^39 + 6*u^37 + 12*u^35 + 20*u^33 + 28*u^31 + 31*u^29 + 31*u^27 + 
#            27*u^25 + 19*u^23 + 12*u^21 + 7*u^19 + 3*u^17 + u^15 + u^13, 
#          2*u^35 + 3*u^33 + 9*u^31 + 13*u^29 + 12*u^27 + 14*u^25 + 10*u^23 + 
#            6*u^21 + 6*u^19 + 2*u^17 + u^15 + u^13, u^30, 
#          2*u^35 + 2*u^33 + 6*u^31 + 10*u^29 + 7*u^27 + 7*u^25 + 5*u^23 + u^
#            21 + u^19, 3*u^31 + 6*u^29 + 8*u^27 + 12*u^25 + 10*u^23 + 6*u^
#            21 + 6*u^19 + 2*u^17 + u^15 + u^13, u^28 + u^26, 
#          u^33 + 2*u^31 + 5*u^29 + 10*u^27 + 8*u^25 + 8*u^23 + 5*u^21 + u^
#            19 + u^17, 0*u^0, u^29 + 2*u^25, 
#          u^31 + 3*u^29 + 3*u^27 + 5*u^25 + 7*u^23 + 4*u^21 + 6*u^19 + 4*u^
#            17 + u^15 + 2*u^13, u^27 + 5*u^25 + 5*u^23 + 5*u^21 + 5*u^19 + 
#            2*u^17 + u^15 + u^13, 2*u^27 + 5*u^25 + 6*u^23 + 5*u^21 + u^
#            19 + u^17, u^25, 0*u^0, 2*u^25, u^27 + 5*u^25 + 8*u^23 + 8*u^21 + 
#            6*u^19 + 3*u^17 + u^15 + u^13, u^27 + 4*u^25 + 7*u^23 + 8*u^21 + 
#            6*u^19 + 3*u^17 + u^15 + u^13, 0*u^0, 
#          u^26 + 4*u^24 + 6*u^22 + 5*u^20 + 2*u^18 + u^16, 
#          2*u^23 + 3*u^21 + 5*u^19 + 2*u^17 + u^15 + u^13, 0*u^0, 
#          3*u^22 + 3*u^20 + u^18 + u^16, u^23 + 3*u^21 + 5*u^19 + 2*u^17 + u^
#            15 + u^13, 3*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          u^23 + 2*u^21 + 5*u^19 + 4*u^17 + u^15 + 2*u^13, 
#          2*u^22 + 3*u^20 + 2*u^18 + u^16, 0*u^0, 
#          u^21 + 3*u^19 + 2*u^17 + u^15 + u^13, u^18, 0*u^0, 
#          2*u^19 + 2*u^17 + u^15 + u^13, 2*u^19 + 4*u^17 + u^15 + 2*u^13, 
#          u^17 + u^13, 2*u^19 + 3*u^17 + 2*u^15 + 2*u^13, 2*u^17 + 2*u^13, 
#          0*u^0, 2*u^17 + u^13, u^19 + u^17 + u^15 + u^13, u^13, 0*u^0, 
#          0*u^0, 0*u^0, u^13, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, u^13, 0*u^0, u^13, u^13, 
#          0*u^0, 0*u^0, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ 2*u^92 + 2*u^90 + 4*u^88 + 8*u^86 + 10*u^84 + 14*u^82 + 24*u^80 + 
#            26*u^78 + 36*u^76 + 48*u^74 + 52*u^72 + 64*u^70 + 80*u^68 + 80*u^
#            66 + 96*u^64 + 106*u^62 + 104*u^60 + 116*u^58 + 122*u^56 + 112*u^
#            54 + 122*u^52 + 116*u^50 + 104*u^48 + 106*u^46 + 96*u^44 + 80*u^
#            42 + 80*u^40 + 64*u^38 + 52*u^36 + 48*u^34 + 36*u^32 + 26*u^30 + 
#            24*u^28 + 14*u^26 + 10*u^24 + 8*u^22 + 4*u^20 + 2*u^18 + 2*u^16, 
#          u^74 + u^72 + 4*u^70 + 8*u^68 + 11*u^66 + 20*u^64 + 29*u^62 + 35*u^
#            60 + 50*u^58 + 61*u^56 + 66*u^54 + 81*u^52 + 85*u^50 + 84*u^48 + 
#            91*u^46 + 86*u^44 + 76*u^42 + 77*u^40 + 63*u^38 + 52*u^36 + 48*u^
#            34 + 36*u^32 + 26*u^30 + 24*u^28 + 14*u^26 + 10*u^24 + 8*u^22 + 
#            4*u^20 + 2*u^18 + 2*u^16, u^60 + 4*u^58 + 8*u^56 + 13*u^54 + 24*u^
#            52 + 30*u^50 + 40*u^48 + 51*u^46 + 55*u^44 + 56*u^42 + 62*u^40 + 
#            53*u^38 + 48*u^36 + 45*u^34 + 35*u^32 + 26*u^30 + 24*u^28 + 14*u^
#            26 + 10*u^24 + 8*u^22 + 4*u^20 + 2*u^18 + 2*u^16, 
#          u^52 + 4*u^50 + 6*u^48 + 12*u^46 + 20*u^44 + 22*u^42 + 32*u^40 + 
#            34*u^38 + 33*u^36 + 35*u^34 + 31*u^32 + 23*u^30 + 23*u^28 + 14*u^
#            26 + 10*u^24 + 8*u^22 + 4*u^20 + 2*u^18 + 2*u^16, u^39, 
#          u^53 + 2*u^51 + 3*u^49 + 9*u^47 + 9*u^45 + 15*u^43 + 20*u^41 + 19*u^
#            39 + 21*u^37 + 21*u^35 + 14*u^33 + 13*u^31 + 9*u^29 + 4*u^27 + 
#            3*u^25 + u^23, u^44 + 2*u^42 + 7*u^40 + 11*u^38 + 16*u^36 + 21*u^
#            34 + 21*u^32 + 19*u^30 + 20*u^28 + 13*u^26 + 10*u^24 + 8*u^22 + 
#            4*u^20 + 2*u^18 + 2*u^16, 2*u^43 + 5*u^41 + 7*u^39 + 13*u^37 + 
#            15*u^35 + 12*u^33 + 12*u^31 + 9*u^29 + 4*u^27 + 3*u^25 + u^23, 
#          u^36, u^40 + 3*u^38 + 7*u^36 + 14*u^34 + 16*u^32 + 17*u^30 + 19*u^
#            28 + 13*u^26 + 10*u^24 + 8*u^22 + 4*u^20 + 2*u^18 + 2*u^16, 
#          u^40 + 3*u^38 + 8*u^36 + 14*u^34 + 20*u^32 + 24*u^30 + 26*u^28 + 
#            22*u^26 + 18*u^24 + 12*u^22 + 7*u^20 + 3*u^18 + 2*u^16, 
#          u^34 + 5*u^32 + 6*u^30 + 10*u^28 + 9*u^26 + 7*u^24 + 7*u^22 + 4*u^
#            20 + 2*u^18 + 2*u^16, u^33, u^34 + 4*u^32 + 4*u^30 + 6*u^28 + 7*u^
#            26 + 3*u^24 + 3*u^22 + u^20, u^32 + 2*u^30 + 6*u^28 + 7*u^26 + 
#            6*u^24 + 7*u^22 + 4*u^20 + 2*u^18 + 2*u^16, 0*u^0, 
#          u^32 + 5*u^30 + 6*u^28 + 8*u^26 + 8*u^24 + 4*u^22 + 3*u^20 + u^18, 
#          0*u^0, u^32 + u^30 + 3*u^28 + 2*u^26 + u^24 + u^22, 
#          u^32 + u^30 + 3*u^28 + 5*u^26 + 4*u^24 + 6*u^22 + 6*u^20 + 3*u^18 + 
#            4*u^16 + u^14, u^28 + u^26 + 3*u^24 + 4*u^22 + 3*u^20 + 2*u^18 + 
#            2*u^16, u^28 + 3*u^26 + 5*u^24 + 3*u^22 + 3*u^20 + u^18, 0*u^0, 
#          0*u^0, u^28 + u^26 + u^24 + u^22, u^26 + 4*u^24 + 5*u^22 + 5*u^20 + 
#            3*u^18 + 2*u^16, u^26 + 4*u^24 + 5*u^22 + 5*u^20 + 3*u^18 + 2*u^16
#            , 0*u^0, u^27 + 3*u^25 + 5*u^23 + 5*u^21 + 4*u^19 + 2*u^17, 
#          2*u^22 + 2*u^20 + 2*u^18 + 2*u^16, 0*u^0, 
#          u^23 + u^21 + 2*u^19 + u^17, 2*u^22 + 2*u^20 + 2*u^18 + 2*u^16, 
#          u^23 + 2*u^21 + 3*u^19 + 2*u^17, 
#          2*u^22 + 3*u^20 + 3*u^18 + 4*u^16 + u^14, 
#          u^23 + 2*u^21 + 3*u^19 + 3*u^17, 0*u^0, u^18 + 2*u^16, 
#          u^21 + u^19 + u^17, 0*u^0, u^18 + 2*u^16, 
#          u^20 + 2*u^18 + 4*u^16 + u^14, u^16, u^18 + 3*u^16 + 2*u^14, 
#          u^16 + u^14, 0*u^0, u^20 + u^18 + 2*u^16 + 2*u^14, 
#          u^18 + 2*u^16 + u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, u^14, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 
#          0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^92 + u^90 + 2*u^88 + 3*u^86 + 4*u^84 + 5*u^82 + 9*u^80 + 9*u^78 + 
#            13*u^76 + 16*u^74 + 18*u^72 + 21*u^70 + 27*u^68 + 26*u^66 + 32*u^
#            64 + 34*u^62 + 34*u^60 + 37*u^58 + 40*u^56 + 36*u^54 + 40*u^52 + 
#            37*u^50 + 34*u^48 + 34*u^46 + 32*u^44 + 26*u^42 + 27*u^40 + 21*u^
#            38 + 18*u^36 + 16*u^34 + 13*u^32 + 9*u^30 + 9*u^28 + 5*u^26 + 4*u^
#            24 + 3*u^22 + 2*u^20 + u^18 + u^16, 
#          u^70 + 2*u^68 + 3*u^66 + 6*u^64 + 9*u^62 + 11*u^60 + 16*u^58 + 20*u^
#            56 + 22*u^54 + 27*u^52 + 28*u^50 + 28*u^48 + 30*u^46 + 29*u^44 + 
#            25*u^42 + 26*u^40 + 21*u^38 + 18*u^36 + 16*u^34 + 13*u^32 + 9*u^
#            30 + 9*u^28 + 5*u^26 + 4*u^24 + 3*u^22 + 2*u^20 + u^18 + u^16, 
#          u^58 + 3*u^56 + 4*u^54 + 8*u^52 + 10*u^50 + 14*u^48 + 17*u^46 + 
#            20*u^44 + 19*u^42 + 22*u^40 + 18*u^38 + 17*u^36 + 15*u^34 + 13*u^
#            32 + 9*u^30 + 9*u^28 + 5*u^26 + 4*u^24 + 3*u^22 + 2*u^20 + u^
#            18 + u^16, u^48 + 2*u^46 + 5*u^44 + 6*u^42 + 10*u^40 + 11*u^38 + 
#            12*u^36 + 12*u^34 + 12*u^32 + 8*u^30 + 9*u^28 + 5*u^26 + 4*u^24 + 
#            3*u^22 + 2*u^20 + u^18 + u^16, u^51 + u^49 + 2*u^47 + 2*u^45 + 
#            3*u^43 + 2*u^41 + 3*u^39 + u^37 + u^35, 
#          u^47 + u^45 + 3*u^43 + 5*u^41 + 5*u^39 + 6*u^37 + 7*u^35 + 4*u^33 + 
#            4*u^31 + 3*u^29 + u^27 + u^25, 
#          2*u^40 + 4*u^38 + 6*u^36 + 8*u^34 + 9*u^32 + 7*u^30 + 8*u^28 + 5*u^
#            26 + 4*u^24 + 3*u^22 + 2*u^20 + u^18 + u^16, 
#          u^41 + 2*u^39 + 4*u^37 + 5*u^35 + 4*u^33 + 4*u^31 + 3*u^29 + u^
#            27 + u^25, u^40 + u^38 + 2*u^36 + u^34 + u^32, 
#          u^38 + 3*u^36 + 6*u^34 + 7*u^32 + 7*u^30 + 8*u^28 + 5*u^26 + 4*u^
#            24 + 3*u^22 + 2*u^20 + u^18 + u^16, 
#          u^36 + 2*u^34 + 4*u^32 + 5*u^30 + 7*u^28 + 6*u^26 + 6*u^24 + 4*u^
#            22 + 3*u^20 + u^18 + u^16, u^34 + 3*u^32 + 3*u^30 + 5*u^28 + 4*u^
#            26 + 3*u^24 + 3*u^22 + 2*u^20 + u^18 + u^16, 0*u^0, 
#          u^34 + 2*u^32 + 2*u^30 + 2*u^28 + 3*u^26 + u^24 + u^22, 
#          u^28 + 2*u^26 + 2*u^24 + 3*u^22 + 2*u^20 + u^18 + u^16, 
#          u^31 + 2*u^29 + 2*u^27 + 2*u^25 + u^23, u^26 + 2*u^24 + u^22 + u^20,
#          0*u^0, 0*u^0, u^20 + u^16, u^24 + 2*u^22 + 2*u^20 + u^18 + u^16, 
#          2*u^24 + u^22 + u^20, u^26 + u^24 + u^22, 0*u^0, 0*u^0, 
#          u^26 + 3*u^24 + 4*u^22 + 3*u^20 + u^18 + u^16, 
#          u^24 + 2*u^22 + 2*u^20 + u^18 + u^16, u^26 + u^24 + u^22, 
#          u^23 + u^21 + u^19, 2*u^22 + 2*u^20 + u^18 + u^16, 0*u^0, 
#          u^23 + u^21 + u^19, u^22 + u^20 + u^18 + u^16, u^21 + u^19, u^16, 
#          0*u^0, u^17, u^20 + u^18 + u^16, 0*u^0, u^20 + u^18, 
#          u^20 + u^18 + u^16, u^16, u^16, 0*u^0, u^16 + u^14, 
#          u^19 + u^17 + u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^96 + u^94 + 4*u^92 + 6*u^90 + 10*u^88 + 16*u^86 + 24*u^84 + 31*u^
#            82 + 47*u^80 + 57*u^78 + 74*u^76 + 93*u^74 + 110*u^72 + 128*u^
#            70 + 154*u^68 + 165*u^66 + 188*u^64 + 204*u^62 + 213*u^60 + 225*u^
#            58 + 235*u^56 + 228*u^54 + 235*u^52 + 225*u^50 + 213*u^48 + 204*u^
#            46 + 188*u^44 + 165*u^42 + 154*u^40 + 128*u^38 + 110*u^36 + 93*u^
#            34 + 74*u^32 + 57*u^30 + 47*u^28 + 31*u^26 + 24*u^24 + 16*u^22 + 
#            10*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          u^74 + 2*u^72 + 6*u^70 + 13*u^68 + 21*u^66 + 36*u^64 + 53*u^62 + 
#            71*u^60 + 95*u^58 + 118*u^56 + 135*u^54 + 158*u^52 + 168*u^50 + 
#            174*u^48 + 178*u^46 + 172*u^44 + 158*u^42 + 150*u^40 + 127*u^38 + 
#            110*u^36 + 93*u^34 + 74*u^32 + 57*u^30 + 47*u^28 + 31*u^26 + 24*u^
#            24 + 16*u^22 + 10*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          2*u^60 + 6*u^58 + 14*u^56 + 25*u^54 + 45*u^52 + 61*u^50 + 84*u^48 + 
#            102*u^46 + 115*u^44 + 119*u^42 + 124*u^40 + 111*u^38 + 103*u^36 + 
#            89*u^34 + 73*u^32 + 57*u^30 + 47*u^28 + 31*u^26 + 24*u^24 + 16*u^
#            22 + 10*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          u^52 + 4*u^50 + 10*u^48 + 19*u^46 + 33*u^44 + 45*u^42 + 60*u^40 + 
#            68*u^38 + 73*u^36 + 71*u^34 + 66*u^32 + 53*u^30 + 46*u^28 + 31*u^
#            26 + 24*u^24 + 16*u^22 + 10*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          2*u^51 + 3*u^49 + 5*u^47 + 8*u^45 + 9*u^43 + 8*u^41 + 10*u^39 + 5*u^
#            37 + 4*u^35 + 2*u^33, u^51 + 2*u^49 + 8*u^47 + 12*u^45 + 21*u^
#            43 + 29*u^41 + 34*u^39 + 36*u^37 + 37*u^35 + 28*u^33 + 23*u^31 + 
#            15*u^29 + 7*u^27 + 4*u^25 + u^23, 
#          u^44 + 3*u^42 + 12*u^40 + 24*u^38 + 38*u^36 + 46*u^34 + 50*u^32 + 
#            46*u^30 + 42*u^28 + 30*u^26 + 24*u^24 + 16*u^22 + 10*u^20 + 6*u^
#            18 + 4*u^16 + u^14 + u^12, u^43 + 5*u^41 + 12*u^39 + 22*u^37 + 
#            27*u^35 + 25*u^33 + 22*u^31 + 15*u^29 + 7*u^27 + 4*u^25 + u^23, 
#          u^42 + 3*u^40 + 6*u^38 + 9*u^36 + 5*u^34 + 4*u^32 + 2*u^30, 
#          u^40 + 5*u^38 + 18*u^36 + 32*u^34 + 40*u^32 + 43*u^30 + 41*u^28 + 
#            30*u^26 + 24*u^24 + 16*u^22 + 10*u^20 + 6*u^18 + 4*u^16 + u^
#            14 + u^12, 2*u^38 + 9*u^36 + 18*u^34 + 30*u^32 + 40*u^30 + 45*u^
#            28 + 41*u^26 + 35*u^24 + 23*u^22 + 14*u^20 + 7*u^18 + 4*u^16 + u^
#            14 + u^12, u^36 + 4*u^34 + 12*u^32 + 18*u^30 + 24*u^28 + 23*u^
#            26 + 20*u^24 + 15*u^22 + 10*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          3*u^33 + 3*u^31 + 3*u^29 + 2*u^27, 2*u^34 + 6*u^32 + 11*u^30 + 12*u^
#            28 + 12*u^26 + 6*u^24 + 4*u^22 + u^20, 
#          u^32 + 4*u^30 + 11*u^28 + 15*u^26 + 16*u^24 + 15*u^22 + 10*u^20 + 
#            6*u^18 + 4*u^16 + u^14 + u^12, 
#          2*u^31 + 5*u^29 + 6*u^27 + 5*u^25 + 3*u^23, 
#          3*u^30 + 7*u^28 + 11*u^26 + 11*u^24 + 7*u^22 + 4*u^20 + u^18, u^28, 
#          u^28 + u^26, u^30 + 2*u^28 + 4*u^26 + 5*u^24 + 6*u^22 + 7*u^20 + 
#            5*u^18 + 5*u^16 + 2*u^14 + u^12, 
#          u^28 + 4*u^26 + 10*u^24 + 11*u^22 + 9*u^20 + 6*u^18 + 4*u^16 + u^
#            14 + u^12, 5*u^26 + 8*u^24 + 6*u^22 + 4*u^20 + u^18, 
#          4*u^26 + 4*u^24 + 3*u^22, u^26, u^26, 
#          5*u^26 + 13*u^24 + 15*u^22 + 12*u^20 + 7*u^18 + 4*u^16 + u^14 + u^12
#            , u^26 + 6*u^24 + 10*u^22 + 11*u^20 + 7*u^18 + 4*u^16 + u^14 + u^
#            12, u^26 + u^24 + u^22, 2*u^25 + 6*u^23 + 8*u^21 + 6*u^19 + 2*u^17
#            , 2*u^24 + 6*u^22 + 8*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          2*u^23 + 2*u^21, 2*u^23 + 3*u^21 + 3*u^19 + u^17, 
#          4*u^22 + 7*u^20 + 6*u^18 + 4*u^16 + u^14 + u^12, 
#          u^23 + 5*u^21 + 5*u^19 + 2*u^17, u^22 + 4*u^20 + 5*u^18 + 5*u^16 + 
#            2*u^14 + u^12, u^21 + 2*u^19 + 2*u^17, u^21 + 2*u^19 + u^17, 
#          4*u^20 + 5*u^18 + 4*u^16 + u^14 + u^12, 2*u^21 + 2*u^19 + u^17, 
#          u^20 + u^18, 3*u^20 + 5*u^18 + 4*u^16 + u^14 + u^12, 
#          u^20 + 4*u^18 + 5*u^16 + 2*u^14 + u^12, 
#          u^20 + 2*u^18 + 3*u^16 + u^14 + u^12, 
#          2*u^18 + 3*u^16 + 3*u^14 + u^12, 2*u^18 + 3*u^16 + 3*u^14 + u^12, 
#          u^19 + u^17 + 2*u^15, u^16 + u^14, u^16 + u^14 + u^12, u^14 + u^12, 
#          u^15, 0*u^0, u^14, u^14, u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14 + u^12, u^14 + u^12, 
#          0*u^0, 0*u^0, u^12, u^12, u^12, 0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^84 + u^82 + u^80 + 3*u^78 + 4*u^76 + 4*u^74 + 8*u^72 + 9*u^70 + 
#            10*u^68 + 15*u^66 + 16*u^64 + 17*u^62 + 24*u^60 + 23*u^58 + 24*u^
#            56 + 30*u^54 + 28*u^52 + 28*u^50 + 33*u^48 + 28*u^46 + 28*u^44 + 
#            30*u^42 + 24*u^40 + 23*u^38 + 24*u^36 + 17*u^34 + 16*u^32 + 15*u^
#            30 + 10*u^28 + 9*u^26 + 8*u^24 + 4*u^22 + 4*u^20 + 3*u^18 + u^
#            16 + u^14 + u^12, u^66 + u^64 + 2*u^62 + 5*u^60 + 6*u^58 + 8*u^
#            56 + 13*u^54 + 14*u^52 + 17*u^50 + 22*u^48 + 21*u^46 + 23*u^44 + 
#            26*u^42 + 22*u^40 + 22*u^38 + 23*u^36 + 17*u^34 + 16*u^32 + 15*u^
#            30 + 10*u^28 + 9*u^26 + 8*u^24 + 4*u^22 + 4*u^20 + 3*u^18 + u^
#            16 + u^14 + u^12, u^54 + 2*u^52 + 4*u^50 + 7*u^48 + 8*u^46 + 12*u^
#            44 + 15*u^42 + 15*u^40 + 17*u^38 + 19*u^36 + 15*u^34 + 15*u^32 + 
#            14*u^30 + 10*u^28 + 9*u^26 + 8*u^24 + 4*u^22 + 4*u^20 + 3*u^
#            18 + u^16 + u^14 + u^12, u^48 + u^46 + 2*u^44 + 5*u^42 + 6*u^40 + 
#            7*u^38 + 12*u^36 + 10*u^34 + 11*u^32 + 12*u^30 + 9*u^28 + 8*u^
#            26 + 8*u^24 + 4*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14 + u^12, 
#          0*u^0, u^45 + 2*u^43 + 2*u^41 + 5*u^39 + 5*u^37 + 5*u^35 + 7*u^33 + 
#            5*u^31 + 4*u^29 + 4*u^27 + 2*u^25 + u^23 + u^21, 
#          u^38 + 4*u^36 + 4*u^34 + 6*u^32 + 8*u^30 + 7*u^28 + 7*u^26 + 7*u^
#            24 + 4*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14 + u^12, 
#          u^39 + 2*u^37 + 3*u^35 + 5*u^33 + 4*u^31 + 4*u^29 + 4*u^27 + 2*u^
#            25 + u^23 + u^21, 0*u^0, u^36 + 2*u^34 + 4*u^32 + 6*u^30 + 6*u^
#            28 + 7*u^26 + 7*u^24 + 4*u^22 + 4*u^20 + 3*u^18 + u^16 + u^14 + u^
#            12, u^36 + 2*u^34 + 4*u^32 + 6*u^30 + 8*u^28 + 9*u^26 + 9*u^24 + 
#            7*u^22 + 6*u^20 + 4*u^18 + 2*u^16 + u^14 + u^12, 
#          2*u^30 + 2*u^28 + 3*u^26 + 5*u^24 + 3*u^22 + 3*u^20 + 3*u^18 + u^
#            16 + u^14 + u^12, 0*u^0, 2*u^30 + 2*u^28 + 2*u^26 + 3*u^24 + 2*u^
#            22 + u^20 + u^18, u^30 + u^28 + 2*u^26 + 4*u^24 + 3*u^22 + 3*u^
#            20 + 3*u^18 + u^16 + u^14 + u^12, 0*u^0, 
#          u^28 + 2*u^26 + 2*u^24 + 3*u^22 + 2*u^20 + u^18 + u^16, 0*u^0, 
#          u^24, u^30 + u^28 + u^26 + 3*u^24 + 2*u^22 + 2*u^20 + 3*u^18 + u^
#            16 + u^14 + 2*u^12, u^24 + u^22 + 2*u^20 + 2*u^18 + u^16 + u^
#            14 + u^12, u^24 + 2*u^22 + 2*u^20 + u^18 + u^16, 0*u^0, 0*u^0, 
#          u^24, u^24 + 2*u^22 + 3*u^20 + 3*u^18 + 2*u^16 + u^14 + u^12, 
#          u^24 + 2*u^22 + 3*u^20 + 3*u^18 + 2*u^16 + u^14 + u^12, 0*u^0, 
#          u^23 + 2*u^21 + 2*u^19 + u^17 + u^15, 
#          u^20 + 2*u^18 + u^16 + u^14 + u^12, 0*u^0, 
#          u^21 + u^19 + u^17 + u^15, u^20 + 2*u^18 + u^16 + u^14 + u^12, 
#          u^21 + u^19 + u^17 + u^15, u^20 + 2*u^18 + u^16 + u^14 + 2*u^12, 
#          u^21 + u^19 + u^17 + u^15, 0*u^0, u^18 + u^16 + u^14 + u^12, 0*u^0, 
#          0*u^0, u^18 + u^16 + u^14 + u^12, u^18 + u^16 + u^14 + 2*u^12, 
#          u^12, u^18 + u^16 + u^14 + 2*u^12, 2*u^12, 0*u^0, u^12, 
#          u^18 + u^16 + u^14 + u^12, u^12, 0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^12, u^12, 0*u^0, u^12, u^12, 0*u^0, 0*u^0, u^12, 
#          0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^95 + 2*u^93 + 3*u^91 + 7*u^89 + 10*u^87 + 14*u^85 + 24*u^83 + 30*u^
#            81 + 39*u^79 + 56*u^77 + 64*u^75 + 79*u^73 + 101*u^71 + 108*u^
#            69 + 127*u^67 + 148*u^65 + 150*u^63 + 169*u^61 + 182*u^59 + 176*u^
#            57 + 190*u^55 + 190*u^53 + 176*u^51 + 182*u^49 + 169*u^47 + 150*u^
#            45 + 148*u^43 + 127*u^41 + 108*u^39 + 101*u^37 + 79*u^35 + 64*u^
#            33 + 56*u^31 + 39*u^29 + 30*u^27 + 24*u^25 + 14*u^23 + 10*u^21 + 
#            7*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^73 + 4*u^71 + 6*u^69 + 14*u^67 + 24*u^65 + 33*u^63 + 51*u^61 + 
#            69*u^59 + 81*u^57 + 105*u^55 + 119*u^53 + 126*u^51 + 142*u^49 + 
#            142*u^47 + 135*u^45 + 138*u^43 + 122*u^41 + 107*u^39 + 100*u^37 + 
#            79*u^35 + 64*u^33 + 56*u^31 + 39*u^29 + 30*u^27 + 24*u^25 + 14*u^
#            23 + 10*u^21 + 7*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^61 + 3*u^59 + 7*u^57 + 17*u^55 + 27*u^53 + 41*u^51 + 61*u^49 + 
#            73*u^47 + 85*u^45 + 98*u^43 + 95*u^41 + 92*u^39 + 90*u^37 + 74*u^
#            35 + 63*u^33 + 55*u^31 + 39*u^29 + 30*u^27 + 24*u^25 + 14*u^23 + 
#            10*u^21 + 7*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          2*u^51 + 5*u^49 + 12*u^47 + 20*u^45 + 31*u^43 + 44*u^41 + 49*u^39 + 
#            59*u^37 + 58*u^35 + 52*u^33 + 50*u^31 + 38*u^29 + 29*u^27 + 24*u^
#            25 + 14*u^23 + 10*u^21 + 7*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^52 + u^50 + 4*u^48 + 5*u^46 + 5*u^44 + 8*u^42 + 5*u^40 + 5*u^38 + 
#            4*u^36 + u^34 + u^32, 2*u^50 + 4*u^48 + 7*u^46 + 17*u^44 + 19*u^
#            42 + 26*u^40 + 33*u^38 + 27*u^36 + 28*u^34 + 23*u^32 + 13*u^30 + 
#            10*u^28 + 5*u^26 + u^24 + u^22, 
#          u^43 + 6*u^41 + 13*u^39 + 25*u^37 + 33*u^35 + 37*u^33 + 40*u^31 + 
#            33*u^29 + 28*u^27 + 23*u^25 + 14*u^23 + 10*u^21 + 7*u^19 + 3*u^
#            17 + 2*u^15 + u^13, u^44 + 2*u^42 + 8*u^40 + 16*u^38 + 18*u^36 + 
#            24*u^34 + 21*u^32 + 13*u^30 + 10*u^28 + 5*u^26 + u^24 + u^22, 
#          u^41 + 4*u^39 + 4*u^37 + 5*u^35 + 4*u^33 + u^31 + u^29, 
#          2*u^39 + 10*u^37 + 18*u^35 + 28*u^33 + 36*u^31 + 31*u^29 + 28*u^
#            27 + 23*u^25 + 14*u^23 + 10*u^21 + 7*u^19 + 3*u^17 + 2*u^15 + u^13
#            , u^39 + 4*u^37 + 11*u^35 + 20*u^33 + 29*u^31 + 35*u^29 + 36*u^
#            27 + 31*u^25 + 23*u^23 + 15*u^21 + 8*u^19 + 4*u^17 + 2*u^15 + u^13
#            , 2*u^35 + 5*u^33 + 12*u^31 + 17*u^29 + 17*u^27 + 18*u^25 + 13*u^
#            23 + 9*u^21 + 7*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^34 + u^32 + 3*u^30 + u^28 + u^26, 2*u^35 + 3*u^33 + 7*u^31 + 12*u^
#            29 + 9*u^27 + 8*u^25 + 5*u^23 + u^21 + u^19, 
#          2*u^31 + 6*u^29 + 9*u^27 + 14*u^25 + 12*u^23 + 9*u^21 + 7*u^19 + 
#            3*u^17 + 2*u^15 + u^13, u^32 + 2*u^30 + 5*u^28 + 5*u^26 + 2*u^
#            24 + u^22, u^31 + 4*u^29 + 9*u^27 + 8*u^25 + 9*u^23 + 5*u^21 + u^
#            19 + u^17, 0*u^0, u^29 + 2*u^25, 
#          u^29 + 2*u^27 + 3*u^25 + 5*u^23 + 4*u^21 + 5*u^19 + 4*u^17 + 2*u^
#            15 + 2*u^13, u^27 + 6*u^25 + 7*u^23 + 8*u^21 + 6*u^19 + 3*u^17 + 
#            2*u^15 + u^13, 2*u^27 + 5*u^25 + 7*u^23 + 5*u^21 + u^19 + u^17, 
#          u^27 + 4*u^25 + 2*u^23 + u^21, 0*u^0, 2*u^25, 
#          u^27 + 7*u^25 + 12*u^23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13,
#          3*u^25 + 7*u^23 + 9*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          u^25 + u^23, u^26 + 4*u^24 + 7*u^22 + 6*u^20 + 3*u^18 + u^16, 
#          3*u^23 + 5*u^21 + 6*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^24 + u^22 + u^20, 4*u^22 + 3*u^20 + u^18 + u^16, 
#          u^23 + 4*u^21 + 6*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          4*u^22 + 4*u^20 + 3*u^18 + u^16, u^21 + 4*u^19 + 4*u^17 + 2*u^15 + 
#            2*u^13, u^22 + 2*u^20 + 2*u^18 + u^16, u^20 + u^18, 
#          u^21 + 4*u^19 + 3*u^17 + 2*u^15 + u^13, u^22 + u^20 + 2*u^18, u^19, 
#          4*u^19 + 3*u^17 + 2*u^15 + u^13, 2*u^19 + 4*u^17 + 2*u^15 + 2*u^13, 
#          u^19 + 2*u^17 + u^15 + u^13, u^19 + 2*u^17 + 3*u^15 + 2*u^13, 
#          2*u^17 + u^15 + 2*u^13, 2*u^18 + u^16 + u^14, u^17 + u^13, 
#          u^15 + u^13, u^13, u^14, 0*u^0, 0*u^0, u^13, u^15, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, 
#          2*u^13, 0*u^0, u^13, u^13, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^94 + 2*u^90 + 2*u^88 + 3*u^86 + 5*u^84 + 8*u^82 + 7*u^80 + 14*u^
#            78 + 14*u^76 + 18*u^74 + 23*u^72 + 27*u^70 + 27*u^68 + 37*u^66 + 
#            35*u^64 + 40*u^62 + 44*u^60 + 45*u^58 + 43*u^56 + 50*u^54 + 43*u^
#            52 + 45*u^50 + 44*u^48 + 40*u^46 + 35*u^44 + 37*u^42 + 27*u^40 + 
#            27*u^38 + 23*u^36 + 18*u^34 + 14*u^32 + 14*u^30 + 7*u^28 + 8*u^
#            26 + 5*u^24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          u^72 + u^70 + 2*u^68 + 5*u^66 + 6*u^64 + 10*u^62 + 15*u^60 + 18*u^
#            58 + 22*u^56 + 29*u^54 + 29*u^52 + 34*u^50 + 36*u^48 + 35*u^46 + 
#            33*u^44 + 35*u^42 + 27*u^40 + 27*u^38 + 23*u^36 + 18*u^34 + 14*u^
#            32 + 14*u^30 + 7*u^28 + 8*u^26 + 5*u^24 + 3*u^22 + 2*u^20 + 2*u^
#            18 + u^14, u^58 + 2*u^56 + 6*u^54 + 7*u^52 + 14*u^50 + 16*u^48 + 
#            21*u^46 + 22*u^44 + 27*u^42 + 22*u^40 + 25*u^38 + 21*u^36 + 18*u^
#            34 + 14*u^32 + 14*u^30 + 7*u^28 + 8*u^26 + 5*u^24 + 3*u^22 + 2*u^
#            20 + 2*u^18 + u^14, u^50 + 2*u^48 + 4*u^46 + 5*u^44 + 10*u^42 + 
#            10*u^40 + 14*u^38 + 15*u^36 + 15*u^34 + 12*u^32 + 14*u^30 + 7*u^
#            28 + 8*u^26 + 5*u^24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          u^49 + u^47 + 3*u^45 + 2*u^43 + 4*u^41 + 2*u^39 + 3*u^37 + u^35 + u^
#            33, u^49 + 3*u^45 + 3*u^43 + 4*u^41 + 6*u^39 + 7*u^37 + 5*u^35 + 
#            7*u^33 + 4*u^31 + 2*u^29 + 2*u^27, 
#          u^42 + 2*u^40 + 5*u^38 + 8*u^36 + 10*u^34 + 10*u^32 + 12*u^30 + 7*u^
#            28 + 8*u^26 + 5*u^24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          2*u^39 + 3*u^37 + 4*u^35 + 6*u^33 + 4*u^31 + 2*u^29 + 2*u^27, 
#          u^42 + u^40 + 3*u^38 + 2*u^36 + 3*u^34 + u^32 + u^30, 
#          u^38 + 4*u^36 + 6*u^34 + 9*u^32 + 11*u^30 + 7*u^28 + 8*u^26 + 5*u^
#            24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          u^36 + 3*u^34 + 5*u^32 + 8*u^30 + 8*u^28 + 9*u^26 + 6*u^24 + 5*u^
#            22 + 2*u^20 + 2*u^18 + u^14, u^34 + u^32 + 6*u^30 + 4*u^28 + 6*u^
#            26 + 5*u^24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          u^35 + u^33 + 2*u^31 + u^29 + u^27, 2*u^30 + 3*u^28 + u^26 + 2*u^24,
#          u^30 + 2*u^28 + 3*u^26 + 4*u^24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          u^29 + 2*u^27 + 2*u^25 + u^23, 2*u^28 + u^26 + u^24 + 2*u^22, u^28, 
#          0*u^0, u^28 + u^26 + 2*u^24 + 2*u^22 + u^20 + 2*u^18 + u^14, 
#          2*u^26 + 2*u^24 + 3*u^22 + 2*u^20 + 2*u^18 + u^14, u^26 + 2*u^22, 
#          u^26 + 2*u^24 + u^22, u^26, 0*u^0, 
#          u^26 + 3*u^24 + 4*u^22 + 2*u^20 + 2*u^18 + u^14, 
#          2*u^22 + 2*u^20 + 2*u^18 + u^14, u^24, 2*u^21 + u^19, 
#          u^24 + u^22 + 2*u^20 + 2*u^18 + u^14, u^23 + u^21, u^21, 
#          2*u^20 + 2*u^18 + u^14, u^21 + u^19, u^20 + 2*u^18 + u^14, 0*u^0, 
#          u^21 + 2*u^19, u^20 + 2*u^18 + u^14, u^19, u^20, 
#          u^20 + 2*u^18 + u^14, 2*u^18 + u^14, 2*u^18 + u^14, u^14, 
#          2*u^18 + u^16 + u^14, 0*u^0, 0*u^0, u^14, u^14, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, 0*u^0, 
#          0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^93 + 3*u^91 + 4*u^89 + 8*u^87 + 14*u^85 + 18*u^83 + 29*u^81 + 41*u^
#            79 + 49*u^77 + 69*u^75 + 86*u^73 + 98*u^71 + 126*u^69 + 142*u^
#            67 + 155*u^65 + 185*u^63 + 194*u^61 + 204*u^59 + 227*u^57 + 223*u^
#            55 + 226*u^53 + 236*u^51 + 219*u^49 + 213*u^47 + 208*u^45 + 182*u^
#            43 + 170*u^41 + 156*u^39 + 127*u^37 + 113*u^35 + 96*u^33 + 72*u^
#            31 + 61*u^29 + 47*u^27 + 31*u^25 + 25*u^23 + 17*u^21 + 9*u^19 + 
#            7*u^17 + 3*u^15 + u^13 + u^11, u^73 + 2*u^71 + 7*u^69 + 12*u^67 + 
#            20*u^65 + 37*u^63 + 50*u^61 + 69*u^59 + 96*u^57 + 112*u^55 + 
#            134*u^53 + 158*u^51 + 163*u^49 + 174*u^47 + 180*u^45 + 167*u^43 + 
#            162*u^41 + 151*u^39 + 126*u^37 + 113*u^35 + 96*u^33 + 72*u^31 + 
#            61*u^29 + 47*u^27 + 31*u^25 + 25*u^23 + 17*u^21 + 9*u^19 + 7*u^
#            17 + 3*u^15 + u^13 + u^11, 2*u^59 + 7*u^57 + 13*u^55 + 27*u^53 + 
#            45*u^51 + 61*u^49 + 85*u^47 + 103*u^45 + 111*u^43 + 123*u^41 + 
#            123*u^39 + 111*u^37 + 105*u^35 + 91*u^33 + 71*u^31 + 61*u^29 + 
#            47*u^27 + 31*u^25 + 25*u^23 + 17*u^21 + 9*u^19 + 7*u^17 + 3*u^
#            15 + u^13 + u^11, 2*u^51 + 6*u^49 + 12*u^47 + 24*u^45 + 36*u^43 + 
#            49*u^41 + 65*u^39 + 70*u^37 + 74*u^35 + 75*u^33 + 63*u^31 + 56*u^
#            29 + 46*u^27 + 31*u^25 + 25*u^23 + 17*u^21 + 9*u^19 + 7*u^17 + 
#            3*u^15 + u^13 + u^11, u^50 + u^48 + 3*u^46 + 5*u^44 + 4*u^42 + 
#            6*u^40 + 5*u^38 + 3*u^36 + 3*u^34 + u^32, 
#          u^52 + u^50 + 5*u^48 + 11*u^46 + 15*u^44 + 27*u^42 + 33*u^40 + 35*u^
#            38 + 42*u^36 + 36*u^34 + 29*u^32 + 25*u^30 + 14*u^28 + 8*u^26 + 
#            5*u^24 + u^22, u^43 + 5*u^41 + 15*u^39 + 26*u^37 + 39*u^35 + 48*u^
#            33 + 48*u^31 + 48*u^29 + 41*u^27 + 30*u^25 + 25*u^23 + 17*u^21 + 
#            9*u^19 + 7*u^17 + 3*u^15 + u^13 + u^11, 
#          3*u^42 + 7*u^40 + 15*u^38 + 26*u^36 + 27*u^34 + 26*u^32 + 24*u^30 + 
#            14*u^28 + 8*u^26 + 5*u^24 + u^22, 
#          u^41 + 2*u^39 + 5*u^37 + 5*u^35 + 3*u^33 + 3*u^31 + u^29, 
#          2*u^39 + 7*u^37 + 21*u^35 + 33*u^33 + 39*u^31 + 45*u^29 + 40*u^27 + 
#            30*u^25 + 25*u^23 + 17*u^21 + 9*u^19 + 7*u^17 + 3*u^15 + u^13 + u^
#            11, u^39 + 5*u^37 + 14*u^35 + 26*u^33 + 39*u^31 + 49*u^29 + 51*u^
#            27 + 46*u^25 + 37*u^23 + 25*u^21 + 14*u^19 + 8*u^17 + 3*u^15 + u^
#            13 + u^11, u^35 + 5*u^33 + 12*u^31 + 18*u^29 + 25*u^27 + 22*u^
#            25 + 20*u^23 + 16*u^21 + 9*u^19 + 7*u^17 + 3*u^15 + u^13 + u^11, 
#          u^34 + 3*u^32 + 2*u^30 + 3*u^28 + u^26, 
#          3*u^33 + 8*u^31 + 10*u^29 + 15*u^27 + 11*u^25 + 7*u^23 + 5*u^21 + u^
#            19, 2*u^31 + 6*u^29 + 15*u^27 + 16*u^25 + 18*u^23 + 16*u^21 + 9*u^
#            19 + 7*u^17 + 3*u^15 + u^13 + u^11, 
#          2*u^30 + 3*u^28 + 3*u^26 + 3*u^24 + u^22, 
#          2*u^31 + 7*u^29 + 11*u^27 + 16*u^25 + 12*u^23 + 8*u^21 + 5*u^19 + u^
#            17, 0*u^0, u^29 + 3*u^27 + u^25 + u^23, 
#          u^31 + 2*u^29 + 6*u^27 + 8*u^25 + 9*u^23 + 12*u^21 + 9*u^19 + 8*u^
#            17 + 6*u^15 + 2*u^13 + u^11, 2*u^27 + 5*u^25 + 11*u^23 + 11*u^
#            21 + 8*u^19 + 7*u^17 + 3*u^15 + u^13 + u^11, 
#          2*u^27 + 8*u^25 + 9*u^23 + 7*u^21 + 5*u^19 + u^17, 
#          u^27 + 2*u^25 + 3*u^23 + u^21, 0*u^0, u^27 + u^25 + u^23, 
#          u^27 + 6*u^25 + 13*u^23 + 15*u^21 + 12*u^19 + 8*u^17 + 3*u^15 + u^
#            13 + u^11, 3*u^25 + 9*u^23 + 13*u^21 + 12*u^19 + 8*u^17 + 3*u^
#            15 + u^13 + u^11, 0*u^0, u^26 + 5*u^24 + 9*u^22 + 10*u^20 + 7*u^
#            18 + 2*u^16, 2*u^23 + 7*u^21 + 7*u^19 + 7*u^17 + 3*u^15 + u^
#            13 + u^11, 2*u^22 + u^20, u^24 + 2*u^22 + 4*u^20 + 4*u^18 + u^16, 
#          u^23 + 6*u^21 + 7*u^19 + 7*u^17 + 3*u^15 + u^13 + u^11, 
#          2*u^22 + 7*u^20 + 6*u^18 + 2*u^16, 
#          4*u^21 + 7*u^19 + 8*u^17 + 6*u^15 + 2*u^13 + u^11, 
#          u^22 + 4*u^20 + 5*u^18 + 2*u^16, 2*u^20 + 2*u^18, 
#          u^21 + 4*u^19 + 6*u^17 + 3*u^15 + u^13 + u^11, 
#          3*u^20 + 2*u^18 + u^16, 0*u^0, 3*u^19 + 6*u^17 + 3*u^15 + u^13 + u^
#            11, 4*u^19 + 7*u^17 + 6*u^15 + 2*u^13 + u^11, 
#          u^19 + 2*u^17 + 2*u^15 + u^13 + u^11, 
#          u^19 + 4*u^17 + 5*u^15 + 3*u^13 + u^11, 
#          u^19 + 2*u^17 + 4*u^15 + 2*u^13 + u^11, u^16 + u^14, 
#          u^19 + u^17 + 3*u^15 + u^13, 2*u^17 + 3*u^15 + 2*u^13 + u^11, 
#          u^13 + u^11, u^14, 0*u^0, 0*u^0, u^15 + u^13, u^13, u^14, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          2*u^13 + u^11, u^13 + u^11, u^14, u^13, u^13 + u^11, u^11, u^11, 
#          u^11, u^12, 0*u^0, u^11, 0*u^0, 0*u^0, 0*u^0, u^11, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^94 + u^92 + 3*u^90 + 7*u^88 + 9*u^86 + 17*u^84 + 25*u^82 + 31*u^
#            80 + 47*u^78 + 61*u^76 + 72*u^74 + 96*u^72 + 113*u^70 + 127*u^
#            68 + 156*u^66 + 170*u^64 + 182*u^62 + 208*u^60 + 213*u^58 + 219*u^
#            56 + 236*u^54 + 226*u^52 + 223*u^50 + 227*u^48 + 204*u^46 + 194*u^
#            44 + 185*u^42 + 155*u^40 + 142*u^38 + 126*u^36 + 98*u^34 + 86*u^
#            32 + 69*u^30 + 49*u^28 + 41*u^26 + 29*u^24 + 18*u^22 + 14*u^20 + 
#            8*u^18 + 4*u^16 + 3*u^14 + u^12, 2*u^72 + 4*u^70 + 8*u^68 + 18*u^
#            66 + 27*u^64 + 41*u^62 + 63*u^60 + 80*u^58 + 102*u^56 + 129*u^
#            54 + 142*u^52 + 159*u^50 + 176*u^48 + 172*u^46 + 174*u^44 + 172*u^
#            42 + 150*u^40 + 140*u^38 + 125*u^36 + 98*u^34 + 86*u^32 + 69*u^
#            30 + 49*u^28 + 41*u^26 + 29*u^24 + 18*u^22 + 14*u^20 + 8*u^18 + 
#            4*u^16 + 3*u^14 + u^12, u^60 + 3*u^58 + 10*u^56 + 21*u^54 + 33*u^
#            52 + 54*u^50 + 75*u^48 + 90*u^46 + 110*u^44 + 121*u^42 + 118*u^
#            40 + 120*u^38 + 112*u^36 + 93*u^34 + 84*u^32 + 68*u^30 + 49*u^
#            28 + 41*u^26 + 29*u^24 + 18*u^22 + 14*u^20 + 8*u^18 + 4*u^16 + 
#            3*u^14 + u^12, u^52 + 3*u^50 + 9*u^48 + 18*u^46 + 28*u^44 + 45*u^
#            42 + 56*u^40 + 66*u^38 + 77*u^36 + 72*u^34 + 70*u^32 + 63*u^30 + 
#            47*u^28 + 40*u^26 + 29*u^24 + 18*u^22 + 14*u^20 + 8*u^18 + 4*u^
#            16 + 3*u^14 + u^12, u^49 + 3*u^47 + 3*u^45 + 5*u^43 + 6*u^41 + 
#            4*u^39 + 5*u^37 + 3*u^35 + u^33 + u^31, 
#          u^51 + 4*u^49 + 6*u^47 + 14*u^45 + 22*u^43 + 26*u^41 + 38*u^39 + 
#            39*u^37 + 36*u^35 + 37*u^33 + 26*u^31 + 18*u^29 + 13*u^27 + 5*u^
#            25 + 2*u^23 + u^21, 3*u^42 + 9*u^40 + 20*u^38 + 34*u^36 + 42*u^
#            34 + 50*u^32 + 50*u^30 + 42*u^28 + 38*u^26 + 28*u^24 + 18*u^22 + 
#            14*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^43 + 4*u^41 + 13*u^39 + 19*u^37 + 25*u^35 + 31*u^33 + 24*u^31 + 
#            18*u^29 + 13*u^27 + 5*u^25 + 2*u^23 + u^21, 
#          2*u^40 + 4*u^38 + 4*u^36 + 5*u^34 + 3*u^32 + u^30 + u^28, 
#          5*u^38 + 14*u^36 + 24*u^34 + 39*u^32 + 44*u^30 + 40*u^28 + 38*u^
#            26 + 28*u^24 + 18*u^22 + 14*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12
#            , 3*u^38 + 9*u^36 + 19*u^34 + 33*u^32 + 44*u^30 + 51*u^28 + 51*u^
#            26 + 42*u^24 + 30*u^22 + 19*u^20 + 10*u^18 + 5*u^16 + 3*u^14 + u^
#            12, 3*u^34 + 7*u^32 + 17*u^30 + 21*u^28 + 24*u^26 + 23*u^24 + 
#            16*u^22 + 13*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^35 + u^33 + 3*u^31 + 3*u^29 + u^27 + u^25, 
#          2*u^34 + 3*u^32 + 11*u^30 + 13*u^28 + 12*u^26 + 11*u^24 + 5*u^22 + 
#            2*u^20 + u^18, 5*u^30 + 9*u^28 + 16*u^26 + 20*u^24 + 15*u^22 + 
#            13*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^31 + 2*u^29 + 4*u^27 + 3*u^25 + u^23 + u^21, 
#          u^32 + 3*u^30 + 10*u^28 + 14*u^26 + 14*u^24 + 12*u^22 + 5*u^20 + 
#            2*u^18 + u^16, 0*u^0, u^30 + u^28 + 2*u^26 + 2*u^24, 
#          2*u^30 + 4*u^28 + 6*u^26 + 10*u^24 + 10*u^22 + 10*u^20 + 10*u^18 + 
#            6*u^16 + 4*u^14 + 2*u^12, 4*u^26 + 9*u^24 + 10*u^22 + 11*u^20 + 
#            7*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^28 + 4*u^26 + 9*u^24 + 10*u^22 + 5*u^20 + 2*u^18 + u^16, 
#          2*u^26 + 3*u^24 + u^22 + u^20, 0*u^0, u^26 + 2*u^24, 
#          3*u^26 + 10*u^24 + 15*u^22 + 14*u^20 + 9*u^18 + 5*u^16 + 3*u^14 + u^
#            12, u^26 + 6*u^24 + 12*u^22 + 13*u^20 + 9*u^18 + 5*u^16 + 3*u^
#            14 + u^12, 0*u^0, 2*u^25 + 8*u^23 + 11*u^21 + 8*u^19 + 4*u^17 + u^
#            15, u^24 + 4*u^22 + 8*u^20 + 7*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^23 + u^21 + u^19, u^23 + 5*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          2*u^22 + 8*u^20 + 7*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^23 + 6*u^21 + 6*u^19 + 4*u^17 + u^15, 
#          u^22 + 6*u^20 + 9*u^18 + 6*u^16 + 4*u^14 + 2*u^12, 
#          3*u^21 + 5*u^19 + 3*u^17 + u^15, u^21 + 2*u^19 + u^17, 
#          3*u^20 + 5*u^18 + 4*u^16 + 3*u^14 + u^12, u^21 + 3*u^19 + 2*u^17, 
#          0*u^0, u^20 + 5*u^18 + 4*u^16 + 3*u^14 + u^12, 
#          u^20 + 7*u^18 + 6*u^16 + 4*u^14 + 2*u^12, 
#          2*u^18 + 2*u^16 + 2*u^14 + u^12, 3*u^18 + 4*u^16 + 5*u^14 + 2*u^12, 
#          2*u^18 + 3*u^16 + 3*u^14 + 2*u^12, u^17 + u^13, 
#          2*u^18 + 2*u^16 + u^14 + u^12, u^18 + 2*u^16 + 3*u^14 + 2*u^12, 
#          u^14 + u^12, u^13, 0*u^0, 0*u^0, u^14 + u^12, u^14, u^15, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^14 + 2*u^12, 2*u^12, u^13, u^12, 2*u^12, u^12, u^12, u^12, u^11, 
#          0*u^0, 0*u^0, 0*u^0, u^11, 0*u^0, 0*u^0, u^11, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^90 + u^88 + 3*u^86 + 4*u^84 + 8*u^82 + 11*u^80 + 17*u^78 + 22*u^
#            76 + 31*u^74 + 38*u^72 + 50*u^70 + 58*u^68 + 71*u^66 + 80*u^64 + 
#            93*u^62 + 101*u^60 + 112*u^58 + 117*u^56 + 125*u^54 + 126*u^52 + 
#            130*u^50 + 126*u^48 + 125*u^46 + 117*u^44 + 112*u^42 + 101*u^40 + 
#            93*u^38 + 80*u^36 + 71*u^34 + 58*u^32 + 50*u^30 + 38*u^28 + 31*u^
#            26 + 22*u^24 + 17*u^22 + 11*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^
#            12 + u^10, u^70 + 2*u^68 + 5*u^66 + 9*u^64 + 16*u^62 + 24*u^60 + 
#            35*u^58 + 46*u^56 + 60*u^54 + 71*u^52 + 84*u^50 + 91*u^48 + 99*u^
#            46 + 100*u^44 + 101*u^42 + 95*u^40 + 90*u^38 + 79*u^36 + 71*u^
#            34 + 58*u^32 + 50*u^30 + 38*u^28 + 31*u^26 + 22*u^24 + 17*u^22 + 
#            11*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^58 + 2*u^56 + 7*u^54 + 12*u^52 + 23*u^50 + 32*u^48 + 46*u^46 + 
#            55*u^44 + 66*u^42 + 69*u^40 + 73*u^38 + 68*u^36 + 65*u^34 + 55*u^
#            32 + 49*u^30 + 38*u^28 + 31*u^26 + 22*u^24 + 17*u^22 + 11*u^20 + 
#            8*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^50 + 3*u^48 + 7*u^46 + 13*u^44 + 21*u^42 + 29*u^40 + 38*u^38 + 
#            42*u^36 + 47*u^34 + 44*u^32 + 43*u^30 + 35*u^28 + 30*u^26 + 22*u^
#            24 + 17*u^22 + 11*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^45 + u^43 + 2*u^41 + u^39 + 2*u^37 + u^35 + u^33, 
#          u^49 + 2*u^47 + 6*u^45 + 9*u^43 + 15*u^41 + 18*u^39 + 23*u^37 + 
#            23*u^35 + 23*u^33 + 19*u^31 + 15*u^29 + 10*u^27 + 6*u^25 + 3*u^
#            23 + u^21, u^42 + 3*u^40 + 9*u^38 + 15*u^36 + 24*u^34 + 28*u^32 + 
#            32*u^30 + 29*u^28 + 27*u^26 + 21*u^24 + 17*u^22 + 11*u^20 + 8*u^
#            18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^41 + 4*u^39 + 9*u^37 + 14*u^35 + 17*u^33 + 17*u^31 + 14*u^29 + 
#            10*u^27 + 6*u^25 + 3*u^23 + u^21, u^38 + u^36 + 2*u^34 + u^32 + u^
#            30, u^38 + 4*u^36 + 12*u^34 + 19*u^32 + 26*u^30 + 27*u^28 + 26*u^
#            26 + 21*u^24 + 17*u^22 + 11*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^
#            12 + u^10, u^38 + 3*u^36 + 10*u^34 + 17*u^32 + 27*u^30 + 32*u^
#            28 + 35*u^26 + 31*u^24 + 26*u^22 + 17*u^20 + 11*u^18 + 5*u^16 + 
#            3*u^14 + u^12 + u^10, u^34 + 3*u^32 + 8*u^30 + 12*u^28 + 16*u^
#            26 + 15*u^24 + 14*u^22 + 10*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^
#            12 + u^10, u^31 + u^29 + u^27, 
#          2*u^32 + 4*u^30 + 8*u^28 + 8*u^26 + 8*u^24 + 5*u^22 + 3*u^20 + u^18,
#          2*u^30 + 5*u^28 + 10*u^26 + 12*u^24 + 13*u^22 + 10*u^20 + 8*u^18 + 
#            4*u^16 + 3*u^14 + u^12 + u^10, u^27 + u^25 + u^23, 
#          2*u^30 + 5*u^28 + 9*u^26 + 10*u^24 + 9*u^22 + 6*u^20 + 3*u^18 + u^16
#            , 0*u^0, u^28 + u^26 + u^24, u^30 + 2*u^28 + 5*u^26 + 6*u^24 + 
#            9*u^22 + 8*u^20 + 9*u^18 + 6*u^16 + 5*u^14 + 2*u^12 + u^10, 
#          2*u^26 + 4*u^24 + 8*u^22 + 7*u^20 + 7*u^18 + 4*u^16 + 3*u^14 + u^
#            12 + u^10, 2*u^26 + 5*u^24 + 7*u^22 + 5*u^20 + 3*u^18 + u^16, 
#          u^24 + u^22, 0*u^0, u^24, u^26 + 4*u^24 + 10*u^22 + 10*u^20 + 9*u^
#            18 + 5*u^16 + 3*u^14 + u^12 + u^10, 
#          u^26 + 3*u^24 + 8*u^22 + 10*u^20 + 9*u^18 + 5*u^16 + 3*u^14 + u^
#            12 + u^10, 0*u^0, u^25 + 4*u^23 + 7*u^21 + 7*u^19 + 4*u^17 + 2*u^
#            15, 2*u^22 + 5*u^20 + 6*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^21, 2*u^21 + 3*u^19 + 2*u^17 + u^15, 
#          u^22 + 5*u^20 + 6*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          2*u^21 + 5*u^19 + 3*u^17 + 2*u^15, 
#          u^22 + 4*u^20 + 7*u^18 + 6*u^16 + 5*u^14 + 2*u^12 + u^10, 
#          u^21 + 4*u^19 + 3*u^17 + 2*u^15, u^19, 
#          u^20 + 4*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 2*u^19 + u^17 + u^15,
#          0*u^0, 3*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          4*u^18 + 5*u^16 + 5*u^14 + 2*u^12 + u^10, 
#          u^18 + u^16 + 2*u^14 + u^12 + u^10, 
#          2*u^18 + 4*u^16 + 5*u^14 + 3*u^12 + u^10, 
#          u^18 + 2*u^16 + 3*u^14 + 2*u^12 + u^10, 0*u^0, 
#          u^18 + 2*u^16 + 2*u^14 + u^12, u^18 + 2*u^16 + 3*u^14 + 2*u^12 + u^
#            10, u^14 + u^12 + u^10, 0*u^0, 0*u^0, 0*u^0, u^14 + u^12, 
#          u^14 + u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^14 + 2*u^12 + u^10, u^12 + u^10, 0*u^0, 
#          u^12, u^12 + u^10, u^10, u^10, u^12 + u^10, u^11, 0*u^0, u^10, 
#          0*u^0, 0*u^0, 0*u^0, u^10, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^92 + u^90 + 2*u^88 + 5*u^86 + 6*u^84 + 10*u^82 + 16*u^80 + 19*u^
#            78 + 27*u^76 + 37*u^74 + 42*u^72 + 54*u^70 + 66*u^68 + 72*u^66 + 
#            86*u^64 + 97*u^62 + 101*u^60 + 113*u^58 + 120*u^56 + 119*u^54 + 
#            126*u^52 + 126*u^50 + 119*u^48 + 120*u^46 + 113*u^44 + 101*u^42 + 
#            97*u^40 + 86*u^38 + 72*u^36 + 66*u^34 + 54*u^32 + 42*u^30 + 37*u^
#            28 + 27*u^26 + 19*u^24 + 16*u^22 + 10*u^20 + 6*u^18 + 5*u^16 + 
#            2*u^14 + u^12 + u^10, u^70 + 3*u^68 + 5*u^66 + 11*u^64 + 18*u^
#            62 + 25*u^60 + 38*u^58 + 50*u^56 + 60*u^54 + 75*u^52 + 85*u^50 + 
#            90*u^48 + 99*u^46 + 99*u^44 + 94*u^42 + 93*u^40 + 84*u^38 + 72*u^
#            36 + 66*u^34 + 54*u^32 + 42*u^30 + 37*u^28 + 27*u^26 + 19*u^24 + 
#            16*u^22 + 10*u^20 + 6*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          u^58 + 3*u^56 + 7*u^54 + 15*u^52 + 24*u^50 + 35*u^48 + 49*u^46 + 
#            58*u^44 + 65*u^42 + 72*u^40 + 70*u^38 + 65*u^36 + 62*u^34 + 52*u^
#            32 + 42*u^30 + 37*u^28 + 27*u^26 + 19*u^24 + 16*u^22 + 10*u^20 + 
#            6*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          u^50 + 2*u^48 + 6*u^46 + 13*u^44 + 18*u^42 + 29*u^40 + 37*u^38 + 
#            40*u^36 + 46*u^34 + 44*u^32 + 38*u^30 + 35*u^28 + 27*u^26 + 19*u^
#            24 + 16*u^22 + 10*u^20 + 6*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          u^49 + 2*u^47 + 4*u^45 + 5*u^43 + 6*u^41 + 6*u^39 + 5*u^37 + 4*u^
#            35 + 2*u^33 + u^31, 2*u^47 + 3*u^45 + 7*u^43 + 13*u^41 + 15*u^
#            39 + 20*u^37 + 22*u^35 + 19*u^33 + 17*u^31 + 13*u^29 + 7*u^27 + 
#            4*u^25 + 2*u^23, 3*u^40 + 9*u^38 + 16*u^36 + 26*u^34 + 30*u^32 + 
#            31*u^30 + 31*u^28 + 25*u^26 + 19*u^24 + 16*u^22 + 10*u^20 + 6*u^
#            18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          u^41 + 3*u^39 + 9*u^37 + 14*u^35 + 15*u^33 + 16*u^31 + 13*u^29 + 
#            7*u^27 + 4*u^25 + 2*u^23, u^40 + 3*u^38 + 5*u^36 + 5*u^34 + 4*u^
#            32 + 2*u^30 + u^28, u^38 + 5*u^36 + 15*u^34 + 22*u^32 + 27*u^30 + 
#            30*u^28 + 25*u^26 + 19*u^24 + 16*u^22 + 10*u^20 + 6*u^18 + 5*u^
#            16 + 2*u^14 + u^12 + u^10, 2*u^36 + 7*u^34 + 13*u^32 + 21*u^30 + 
#            28*u^28 + 29*u^26 + 27*u^24 + 22*u^22 + 14*u^20 + 8*u^18 + 5*u^
#            16 + 2*u^14 + u^12 + u^10, u^34 + 5*u^32 + 8*u^30 + 15*u^28 + 
#            17*u^26 + 15*u^24 + 14*u^22 + 10*u^20 + 6*u^18 + 5*u^16 + 2*u^
#            14 + u^12 + u^10, u^33 + 2*u^31 + 3*u^29 + 2*u^27 + u^25, 
#          3*u^32 + 4*u^30 + 8*u^28 + 9*u^26 + 6*u^24 + 4*u^22 + 2*u^20, 
#          u^30 + 5*u^28 + 9*u^26 + 11*u^24 + 13*u^22 + 10*u^20 + 6*u^18 + 5*u^
#            16 + 2*u^14 + u^12 + u^10, 2*u^29 + 4*u^27 + 5*u^25 + 3*u^23 + u^
#            21, u^30 + 2*u^28 + 5*u^26 + 8*u^24 + 6*u^22 + 4*u^20 + 2*u^18, 
#          u^30 + u^26, 0*u^0, u^28 + 3*u^26 + 3*u^24 + 5*u^22 + 6*u^20 + 4*u^
#            18 + 5*u^16 + 3*u^14 + u^12 + u^10, 
#          u^26 + 5*u^24 + 9*u^22 + 8*u^20 + 6*u^18 + 5*u^16 + 2*u^14 + u^
#            12 + u^10, u^26 + 5*u^24 + 5*u^22 + 4*u^20 + 2*u^18, 
#          u^26 + 4*u^24 + 3*u^22 + u^20, u^24, 0*u^0, 
#          u^26 + 7*u^24 + 12*u^22 + 11*u^20 + 8*u^18 + 5*u^16 + 2*u^14 + u^
#            12 + u^10, 2*u^24 + 6*u^22 + 9*u^20 + 8*u^18 + 5*u^16 + 2*u^
#            14 + u^12 + u^10, u^24 + u^22, 2*u^23 + 5*u^21 + 6*u^19 + 3*u^17, 
#          4*u^22 + 6*u^20 + 6*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          u^23 + 2*u^21 + u^19, u^23 + 2*u^21 + 3*u^19 + 2*u^17, 
#          u^22 + 5*u^20 + 6*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          2*u^21 + 5*u^19 + 3*u^17, 2*u^20 + 4*u^18 + 5*u^16 + 3*u^14 + u^
#            12 + u^10, 2*u^19 + 2*u^17, 2*u^19 + 2*u^17, 
#          2*u^20 + 5*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 2*u^19 + u^17, 
#          u^20 + u^18, u^20 + 5*u^18 + 5*u^16 + 2*u^14 + u^12 + u^10, 
#          3*u^18 + 5*u^16 + 3*u^14 + u^12 + u^10, 
#          2*u^18 + 3*u^16 + 2*u^14 + u^12 + u^10, 
#          u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^18 + 3*u^16 + 4*u^14 + u^12 + u^10, u^17 + u^15 + u^13, u^14, 
#          u^16 + u^14 + u^12 + u^10, u^14 + u^12 + u^10, u^13, u^14, u^14, 
#          u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12 + u^10, u^14 + u^12 + u^10, u^13, 
#          0*u^0, u^12 + u^10, u^12 + u^10, u^12 + u^10, u^10, 0*u^0, 0*u^0, 
#          u^10, 0*u^0, 0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^91 + u^89 + 2*u^87 + 4*u^85 + 5*u^83 + 8*u^81 + 12*u^79 + 14*u^77 + 
#            20*u^75 + 26*u^73 + 29*u^71 + 38*u^69 + 44*u^67 + 48*u^65 + 58*u^
#            63 + 62*u^61 + 65*u^59 + 74*u^57 + 74*u^55 + 75*u^53 + 80*u^51 + 
#            75*u^49 + 74*u^47 + 74*u^45 + 65*u^43 + 62*u^41 + 58*u^39 + 48*u^
#            37 + 44*u^35 + 38*u^33 + 29*u^31 + 26*u^29 + 20*u^27 + 14*u^25 + 
#            12*u^23 + 8*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          u^69 + 2*u^67 + 4*u^65 + 9*u^63 + 13*u^61 + 19*u^59 + 28*u^57 + 
#            34*u^55 + 42*u^53 + 51*u^51 + 54*u^49 + 59*u^47 + 63*u^45 + 59*u^
#            43 + 59*u^41 + 56*u^39 + 48*u^37 + 44*u^35 + 38*u^33 + 29*u^31 + 
#            26*u^29 + 20*u^27 + 14*u^25 + 12*u^23 + 8*u^21 + 5*u^19 + 4*u^
#            17 + 2*u^15 + u^13 + u^11, u^57 + 3*u^55 + 7*u^53 + 12*u^51 + 
#            18*u^49 + 27*u^47 + 34*u^45 + 38*u^43 + 44*u^41 + 45*u^39 + 42*u^
#            37 + 41*u^35 + 36*u^33 + 29*u^31 + 26*u^29 + 20*u^27 + 14*u^25 + 
#            12*u^23 + 8*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          2*u^47 + 4*u^45 + 9*u^43 + 13*u^41 + 20*u^39 + 24*u^37 + 27*u^35 + 
#            29*u^33 + 26*u^31 + 24*u^29 + 20*u^27 + 14*u^25 + 12*u^23 + 8*u^
#            21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          u^50 + u^48 + 3*u^46 + 5*u^44 + 4*u^42 + 6*u^40 + 5*u^38 + 3*u^36 + 
#            3*u^34 + u^32, u^46 + 2*u^44 + 5*u^42 + 9*u^40 + 9*u^38 + 14*u^
#            36 + 13*u^34 + 11*u^32 + 10*u^30 + 6*u^28 + 3*u^26 + 2*u^24, 
#          3*u^39 + 8*u^37 + 13*u^35 + 18*u^33 + 20*u^31 + 21*u^29 + 18*u^27 + 
#            14*u^25 + 12*u^23 + 8*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^
#            11, u^40 + 3*u^38 + 8*u^36 + 9*u^34 + 10*u^32 + 10*u^30 + 6*u^
#            28 + 3*u^26 + 2*u^24, u^39 + 4*u^37 + 4*u^35 + 3*u^33 + 3*u^
#            31 + u^29, u^37 + 7*u^35 + 12*u^33 + 16*u^31 + 20*u^29 + 18*u^
#            27 + 14*u^25 + 12*u^23 + 8*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^
#            13 + u^11, 2*u^35 + 4*u^33 + 9*u^31 + 14*u^29 + 17*u^27 + 17*u^
#            25 + 16*u^23 + 11*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          2*u^33 + 5*u^31 + 8*u^29 + 11*u^27 + 11*u^25 + 10*u^23 + 8*u^21 + 
#            5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          u^32 + u^30 + 2*u^28 + u^26, u^33 + 3*u^31 + 4*u^29 + 6*u^27 + 5*u^
#            25 + 3*u^23 + 2*u^21, u^29 + 4*u^27 + 5*u^25 + 8*u^23 + 8*u^21 + 
#            5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          u^30 + 3*u^28 + 4*u^26 + 5*u^24 + 2*u^22, 
#          u^27 + 3*u^25 + 4*u^23 + 3*u^21 + 2*u^19, u^29 + u^27, 0*u^0, 
#          u^25 + u^23 + 2*u^21 + 2*u^19 + 2*u^17 + 2*u^15 + u^13 + u^11, 
#          u^25 + 5*u^23 + 6*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          u^25 + 3*u^23 + 3*u^21 + 2*u^19, 2*u^25 + 4*u^23 + 2*u^21, u^25, 
#          0*u^0, 2*u^25 + 8*u^23 + 9*u^21 + 7*u^19 + 4*u^17 + 2*u^15 + u^
#            13 + u^11, 2*u^23 + 5*u^21 + 6*u^19 + 4*u^17 + 2*u^15 + u^13 + u^
#            11, u^25 + 2*u^23 + u^21, 2*u^22 + 4*u^20 + 3*u^18, 
#          u^23 + 4*u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          2*u^22 + u^20, u^22 + 2*u^20 + 2*u^18, 
#          2*u^21 + 4*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 3*u^20 + 3*u^18, 
#          u^19 + 2*u^17 + 2*u^15 + u^13 + u^11, u^18, u^18 + u^16, 
#          3*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, u^20 + u^18, 2*u^19 + u^17, 
#          3*u^19 + 4*u^17 + 2*u^15 + u^13 + u^11, 
#          2*u^17 + 2*u^15 + u^13 + u^11, u^19 + 2*u^17 + 2*u^15 + u^13 + u^11,
#          u^17 + u^15 + u^13 + u^11, u^17 + 2*u^15 + 2*u^13 + u^11, 
#          u^18 + 2*u^16 + 2*u^14, 0*u^0, u^11, u^13 + u^11, u^14, 0*u^0, 
#          u^13, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^11, 2*u^13 + u^11, 0*u^0, 0*u^0, 
#          u^11, u^13 + u^11, u^11, 0*u^0, 0*u^0, u^12, u^11, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^88 + u^86 + 3*u^84 + 3*u^82 + 6*u^80 + 7*u^78 + 11*u^76 + 13*u^74 + 
#            18*u^72 + 20*u^70 + 27*u^68 + 29*u^66 + 36*u^64 + 38*u^62 + 44*u^
#            60 + 46*u^58 + 51*u^56 + 51*u^54 + 55*u^52 + 52*u^50 + 55*u^48 + 
#            51*u^46 + 51*u^44 + 46*u^42 + 44*u^40 + 38*u^38 + 36*u^36 + 29*u^
#            34 + 27*u^32 + 20*u^30 + 18*u^28 + 13*u^26 + 11*u^24 + 7*u^22 + 
#            6*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, 
#          u^68 + 2*u^66 + 4*u^64 + 6*u^62 + 10*u^60 + 14*u^58 + 20*u^56 + 
#            24*u^54 + 31*u^52 + 34*u^50 + 40*u^48 + 41*u^46 + 44*u^44 + 42*u^
#            42 + 42*u^40 + 37*u^38 + 36*u^36 + 29*u^34 + 27*u^32 + 20*u^30 + 
#            18*u^28 + 13*u^26 + 11*u^24 + 7*u^22 + 6*u^20 + 3*u^18 + 3*u^
#            16 + u^14 + u^12, u^56 + 2*u^54 + 6*u^52 + 8*u^50 + 15*u^48 + 
#            18*u^46 + 26*u^44 + 27*u^42 + 32*u^40 + 30*u^38 + 32*u^36 + 27*u^
#            34 + 26*u^32 + 20*u^30 + 18*u^28 + 13*u^26 + 11*u^24 + 7*u^22 + 
#            6*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, 
#          u^48 + 2*u^46 + 5*u^44 + 7*u^42 + 12*u^40 + 14*u^38 + 19*u^36 + 
#            19*u^34 + 21*u^32 + 18*u^30 + 17*u^28 + 13*u^26 + 11*u^24 + 7*u^
#            22 + 6*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, 
#          u^47 + u^45 + 3*u^43 + 2*u^41 + 4*u^39 + 2*u^37 + 3*u^35 + u^33 + u^
#            31, u^47 + u^45 + 3*u^43 + 4*u^41 + 7*u^39 + 8*u^37 + 9*u^35 + 
#            9*u^33 + 8*u^31 + 6*u^29 + 4*u^27 + 2*u^25 + u^23, 
#          u^40 + 3*u^38 + 7*u^36 + 10*u^34 + 14*u^32 + 14*u^30 + 15*u^28 + 
#            12*u^26 + 11*u^24 + 7*u^22 + 6*u^20 + 3*u^18 + 3*u^16 + u^14 + u^
#            12, u^39 + 3*u^37 + 5*u^35 + 7*u^33 + 7*u^31 + 6*u^29 + 4*u^27 + 
#            2*u^25 + u^23, u^40 + u^38 + 3*u^36 + 2*u^34 + 3*u^32 + u^30 + u^
#            28, 2*u^36 + 5*u^34 + 10*u^32 + 12*u^30 + 14*u^28 + 12*u^26 + 
#            11*u^24 + 7*u^22 + 6*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, 
#          u^36 + 2*u^34 + 6*u^32 + 8*u^30 + 13*u^28 + 13*u^26 + 14*u^24 + 
#            10*u^22 + 8*u^20 + 4*u^18 + 3*u^16 + u^14 + u^12, 
#          2*u^32 + 3*u^30 + 7*u^28 + 7*u^26 + 9*u^24 + 6*u^22 + 6*u^20 + 3*u^
#            18 + 3*u^16 + u^14 + u^12, u^33 + u^31 + 2*u^29 + u^27 + u^25, 
#          2*u^30 + 3*u^28 + 4*u^26 + 3*u^24 + 2*u^22 + u^20, 
#          2*u^28 + 3*u^26 + 6*u^24 + 5*u^22 + 6*u^20 + 3*u^18 + 3*u^16 + u^
#            14 + u^12, u^29 + 2*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          u^28 + 2*u^26 + 3*u^24 + 3*u^22 + 2*u^20 + u^18, 0*u^0, 0*u^0, 
#          u^28 + u^26 + 3*u^24 + 2*u^22 + 4*u^20 + 2*u^18 + 3*u^16 + u^14 + u^
#            12, 3*u^24 + 3*u^22 + 5*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, 
#          2*u^24 + 2*u^22 + 2*u^20 + u^18, u^26 + 2*u^24 + 2*u^22 + u^20, 
#          0*u^0, 0*u^0, 3*u^24 + 5*u^22 + 6*u^20 + 4*u^18 + 3*u^16 + u^14 + u^
#            12, u^24 + 2*u^22 + 4*u^20 + 4*u^18 + 3*u^16 + u^14 + u^12, u^22, 
#          u^23 + 2*u^21 + 3*u^19 + 2*u^17, 
#          u^22 + 3*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, u^23 + u^21 + u^19, 
#          u^21 + u^19 + u^17, 2*u^20 + 3*u^18 + 3*u^16 + u^14 + u^12, 
#          u^21 + 2*u^19 + 2*u^17, u^20 + 2*u^18 + 3*u^16 + u^14 + u^12, u^17, 
#          u^21 + u^19 + 2*u^17, u^20 + 2*u^18 + 3*u^16 + u^14 + u^12, 
#          u^19 + u^17, u^18, 2*u^18 + 3*u^16 + u^14 + u^12, 
#          u^18 + 3*u^16 + u^14 + u^12, u^18 + 2*u^16 + u^14 + u^12, 
#          u^16 + u^14 + u^12, u^18 + 2*u^16 + 2*u^14 + u^12, 
#          u^17 + u^15 + u^13, 0*u^0, u^16 + u^14 + u^12, u^12, u^13, 0*u^0, 
#          u^14, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12, u^12, u^13, 0*u^0, u^12, 
#          u^12, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^87 + 2*u^85 + 3*u^83 + 5*u^81 + 7*u^79 + 11*u^77 + 15*u^75 + 19*u^
#            73 + 25*u^71 + 30*u^69 + 37*u^67 + 44*u^65 + 49*u^63 + 56*u^61 + 
#            61*u^59 + 66*u^57 + 71*u^55 + 72*u^53 + 74*u^51 + 74*u^49 + 72*u^
#            47 + 71*u^45 + 66*u^43 + 61*u^41 + 56*u^39 + 49*u^37 + 44*u^35 + 
#            37*u^33 + 30*u^31 + 25*u^29 + 19*u^27 + 15*u^25 + 11*u^23 + 7*u^
#            21 + 5*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^69 + 2*u^67 + 5*u^65 + 7*u^63 + 12*u^61 + 17*u^59 + 24*u^57 + 
#            31*u^55 + 38*u^53 + 44*u^51 + 51*u^49 + 54*u^47 + 58*u^45 + 57*u^
#            43 + 56*u^41 + 53*u^39 + 48*u^37 + 43*u^35 + 37*u^33 + 30*u^31 + 
#            25*u^29 + 19*u^27 + 15*u^25 + 11*u^23 + 7*u^21 + 5*u^19 + 3*u^
#            17 + 2*u^15 + u^13, u^57 + 3*u^55 + 6*u^53 + 10*u^51 + 16*u^49 + 
#            22*u^47 + 29*u^45 + 34*u^43 + 38*u^41 + 40*u^39 + 39*u^37 + 38*u^
#            35 + 34*u^33 + 29*u^31 + 24*u^29 + 19*u^27 + 15*u^25 + 11*u^23 + 
#            7*u^21 + 5*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          u^49 + 3*u^47 + 6*u^45 + 10*u^43 + 14*u^41 + 19*u^39 + 22*u^37 + 
#            25*u^35 + 25*u^33 + 24*u^31 + 21*u^29 + 18*u^27 + 14*u^25 + 11*u^
#            23 + 7*u^21 + 5*u^19 + 3*u^17 + 2*u^15 + u^13, 0*u^0, 
#          u^50 + 2*u^48 + 4*u^46 + 6*u^44 + 9*u^42 + 12*u^40 + 14*u^38 + 15*u^
#            36 + 15*u^34 + 13*u^32 + 11*u^30 + 8*u^28 + 5*u^26 + 3*u^24 + u^
#            22 + u^20, u^41 + 3*u^39 + 6*u^37 + 10*u^35 + 13*u^33 + 15*u^31 + 
#            16*u^29 + 15*u^27 + 13*u^25 + 10*u^23 + 7*u^21 + 5*u^19 + 3*u^
#            17 + 2*u^15 + u^13, u^42 + 3*u^40 + 5*u^38 + 8*u^36 + 10*u^34 + 
#            11*u^32 + 10*u^30 + 8*u^28 + 5*u^26 + 3*u^24 + u^22 + u^20, 
#          0*u^0, 2*u^37 + 4*u^35 + 8*u^33 + 11*u^31 + 14*u^29 + 14*u^27 + 
#            13*u^25 + 10*u^23 + 7*u^21 + 5*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          2*u^37 + 5*u^35 + 9*u^33 + 14*u^31 + 19*u^29 + 21*u^27 + 21*u^25 + 
#            17*u^23 + 12*u^21 + 8*u^19 + 4*u^17 + 3*u^15 + u^13, 
#          u^33 + 2*u^31 + 5*u^29 + 6*u^27 + 8*u^25 + 7*u^23 + 6*u^21 + 4*u^
#            19 + 3*u^17 + 2*u^15 + u^13, 0*u^0, 
#          u^33 + 2*u^31 + 4*u^29 + 5*u^27 + 6*u^25 + 4*u^23 + 3*u^21 + u^
#            19 + u^17, 2*u^29 + 3*u^27 + 6*u^25 + 6*u^23 + 6*u^21 + 4*u^19 + 
#            3*u^17 + 2*u^15 + u^13, 0*u^0, u^31 + 3*u^29 + 5*u^27 + 7*u^25 + 
#            7*u^23 + 5*u^21 + 3*u^19 + u^17 + u^15, 0*u^0, 
#          u^31 + u^29 + 2*u^27 + 2*u^25 + 2*u^23 + u^21, 
#          u^29 + 2*u^27 + 4*u^25 + 5*u^23 + 6*u^21 + 6*u^19 + 5*u^17 + 4*u^
#            15 + 2*u^13 + u^11, u^25 + 2*u^23 + 3*u^21 + 3*u^19 + 2*u^17 + 
#            2*u^15 + u^13, u^27 + 2*u^25 + 4*u^23 + 4*u^21 + 3*u^19 + u^
#            17 + u^15, 0*u^0, 0*u^0, u^27 + u^25 + 2*u^23 + u^21, 
#          u^25 + 3*u^23 + 5*u^21 + 5*u^19 + 3*u^17 + 3*u^15 + u^13, 
#          u^25 + 3*u^23 + 5*u^21 + 5*u^19 + 3*u^17 + 3*u^15 + u^13, 0*u^0, 
#          u^26 + 2*u^24 + 5*u^22 + 6*u^20 + 4*u^18 + 2*u^16 + u^14, 
#          u^21 + 2*u^19 + 2*u^17 + 2*u^15 + u^13, 0*u^0, 
#          u^22 + 2*u^20 + 2*u^18 + u^16 + u^14, 
#          u^21 + 2*u^19 + 2*u^17 + 2*u^15 + u^13, 
#          u^22 + 3*u^20 + 3*u^18 + 2*u^16 + u^14, 
#          u^21 + 3*u^19 + 4*u^17 + 4*u^15 + 2*u^13 + u^11, 
#          u^22 + 3*u^20 + 4*u^18 + 3*u^16 + u^14, 0*u^0, u^17 + 2*u^15 + u^13,
#          u^20 + u^18 + u^16, 0*u^0, u^17 + 2*u^15 + u^13, 
#          u^19 + 3*u^17 + 4*u^15 + 2*u^13 + u^11, u^15, 
#          2*u^17 + 3*u^15 + 3*u^13 + u^11, u^15 + u^13 + u^11, 0*u^0, 
#          u^19 + 2*u^17 + 3*u^15 + 2*u^13 + u^11, 
#          u^17 + 2*u^15 + 2*u^13 + u^11, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^13 + u^11, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^15 + u^13, 0*u^0, 0*u^0, 0*u^0, u^13 + u^11, u^11, 0*u^0, 
#          u^13 + u^11, u^11, 0*u^0, 0*u^0, u^11, u^12 + u^10, 0*u^0, 0*u^0, 
#          0*u^0, u^10, 0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^10, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^89 + u^87 + u^85 + 4*u^83 + 4*u^81 + 5*u^79 + 11*u^77 + 10*u^75 + 
#            13*u^73 + 22*u^71 + 19*u^69 + 25*u^67 + 35*u^65 + 29*u^63 + 38*u^
#            61 + 46*u^59 + 37*u^57 + 48*u^55 + 51*u^53 + 40*u^51 + 51*u^49 + 
#            48*u^47 + 37*u^45 + 46*u^43 + 38*u^41 + 29*u^39 + 35*u^37 + 25*u^
#            35 + 19*u^33 + 22*u^31 + 13*u^29 + 10*u^27 + 11*u^25 + 5*u^23 + 
#            4*u^21 + 4*u^19 + u^17 + u^15 + u^13, 
#          u^71 + 2*u^67 + 4*u^65 + 4*u^63 + 9*u^61 + 14*u^59 + 13*u^57 + 23*u^
#            55 + 27*u^53 + 26*u^51 + 36*u^49 + 37*u^47 + 32*u^45 + 41*u^43 + 
#            35*u^41 + 29*u^39 + 34*u^37 + 25*u^35 + 19*u^33 + 22*u^31 + 13*u^
#            29 + 10*u^27 + 11*u^25 + 5*u^23 + 4*u^21 + 4*u^19 + u^17 + u^
#            15 + u^13, 3*u^55 + 4*u^53 + 6*u^51 + 13*u^49 + 15*u^47 + 18*u^
#            45 + 26*u^43 + 24*u^41 + 24*u^39 + 29*u^37 + 22*u^35 + 19*u^33 + 
#            21*u^31 + 13*u^29 + 10*u^27 + 11*u^25 + 5*u^23 + 4*u^21 + 4*u^
#            19 + u^17 + u^15 + u^13, u^49 + 2*u^47 + 4*u^45 + 6*u^43 + 11*u^
#            41 + 10*u^39 + 16*u^37 + 17*u^35 + 13*u^33 + 18*u^31 + 13*u^29 + 
#            9*u^27 + 11*u^25 + 5*u^23 + 4*u^21 + 4*u^19 + u^17 + u^15 + u^13, 
#          u^46 + 2*u^42 + u^40 + u^38 + 2*u^36 + u^32, 
#          u^50 + u^48 + u^46 + 5*u^44 + 4*u^42 + 5*u^40 + 11*u^38 + 5*u^36 + 
#            9*u^34 + 9*u^32 + 3*u^30 + 5*u^28 + 3*u^26 + u^22, 
#          u^41 + 2*u^39 + 5*u^37 + 8*u^35 + 8*u^33 + 13*u^31 + 10*u^29 + 9*u^
#            27 + 10*u^25 + 5*u^23 + 4*u^21 + 4*u^19 + u^17 + u^15 + u^13, 
#          u^40 + 4*u^38 + 2*u^36 + 7*u^34 + 7*u^32 + 3*u^30 + 5*u^28 + 3*u^
#            26 + u^22, u^39 + u^37 + u^35 + 2*u^33 + u^29, 
#          2*u^37 + 3*u^35 + 5*u^33 + 11*u^31 + 8*u^29 + 9*u^27 + 10*u^25 + 
#            5*u^23 + 4*u^21 + 4*u^19 + u^17 + u^15 + u^13, 
#          u^37 + 2*u^35 + 5*u^33 + 8*u^31 + 11*u^29 + 12*u^27 + 12*u^25 + 9*u^
#            23 + 7*u^21 + 4*u^19 + 2*u^17 + u^15 + u^13, 
#          2*u^31 + 5*u^29 + 3*u^27 + 7*u^25 + 5*u^23 + 3*u^21 + 4*u^19 + u^
#            17 + u^15 + u^13, u^34 + 2*u^30 + u^26, 
#          4*u^29 + u^27 + 3*u^25 + 3*u^23 + u^19, 
#          2*u^29 + u^27 + 5*u^25 + 4*u^23 + 3*u^21 + 4*u^19 + u^17 + u^15 + u^
#            13, u^28 + u^26 + u^22, u^29 + 4*u^27 + 2*u^25 + 4*u^23 + 3*u^
#            21 + u^17, 0*u^0, u^29 + 2*u^25, u^29 + u^27 + 2*u^25 + 4*u^23 + 
#            2*u^21 + 4*u^19 + 3*u^17 + u^15 + 2*u^13, 
#          2*u^25 + u^23 + 3*u^21 + 3*u^19 + u^17 + u^15 + u^13, 
#          u^27 + u^25 + 2*u^23 + 3*u^21 + u^17, u^25 + u^21, 0*u^0, u^25, 
#          u^25 + 2*u^23 + 4*u^21 + 3*u^19 + 2*u^17 + u^15 + u^13, 
#          u^23 + 3*u^21 + 3*u^19 + 2*u^17 + u^15 + u^13, 0*u^0, 
#          u^24 + 2*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          u^23 + 3*u^19 + u^17 + u^15 + u^13, u^20, u^22 + u^20 + u^16, 
#          3*u^19 + u^17 + u^15 + u^13, u^22 + u^20 + 2*u^18 + u^16, 
#          3*u^19 + 3*u^17 + u^15 + 2*u^13, u^20 + u^18 + u^16, u^20 + u^18, 
#          u^19 + u^17 + u^15 + u^13, 2*u^18, 0*u^0, u^19 + u^17 + u^15 + u^13,
#          u^19 + 3*u^17 + u^15 + 2*u^13, u^17 + u^13, u^15 + 2*u^13, 
#          2*u^17 + 2*u^13, 0*u^0, 2*u^17 + u^13, u^15 + 2*u^13, u^13, 0*u^0, 
#          0*u^0, 0*u^0, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^15, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, 0*u^0, u^12, 0*u^0, 0*u^0, 
#          0*u^0, u^11, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^11, 0*u^0, 
#          u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^89 + 2*u^87 + 4*u^85 + 7*u^83 + 11*u^81 + 17*u^79 + 24*u^77 + 33*u^
#            75 + 44*u^73 + 56*u^71 + 70*u^69 + 85*u^67 + 100*u^65 + 116*u^
#            63 + 131*u^61 + 145*u^59 + 158*u^57 + 168*u^55 + 176*u^53 + 181*u^
#            51 + 182*u^49 + 181*u^47 + 176*u^45 + 168*u^43 + 158*u^41 + 145*u^
#            39 + 131*u^37 + 116*u^35 + 100*u^33 + 85*u^31 + 70*u^29 + 56*u^
#            27 + 44*u^25 + 33*u^23 + 24*u^21 + 17*u^19 + 11*u^17 + 7*u^15 + 
#            4*u^13 + 2*u^11 + u^9, u^69 + 3*u^67 + 7*u^65 + 13*u^63 + 23*u^
#            61 + 35*u^59 + 51*u^57 + 68*u^55 + 87*u^53 + 105*u^51 + 121*u^
#            49 + 134*u^47 + 143*u^45 + 146*u^43 + 145*u^41 + 138*u^39 + 128*u^
#            37 + 115*u^35 + 100*u^33 + 85*u^31 + 70*u^29 + 56*u^27 + 44*u^
#            25 + 33*u^23 + 24*u^21 + 17*u^19 + 11*u^17 + 7*u^15 + 4*u^13 + 
#            2*u^11 + u^9, u^57 + 4*u^55 + 10*u^53 + 20*u^51 + 34*u^49 + 51*u^
#            47 + 69*u^45 + 85*u^43 + 98*u^41 + 105*u^39 + 106*u^37 + 102*u^
#            35 + 93*u^33 + 82*u^31 + 69*u^29 + 56*u^27 + 44*u^25 + 33*u^23 + 
#            24*u^21 + 17*u^19 + 11*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          u^49 + 4*u^47 + 10*u^45 + 19*u^43 + 31*u^41 + 44*u^39 + 56*u^37 + 
#            65*u^35 + 69*u^33 + 68*u^31 + 62*u^29 + 53*u^27 + 43*u^25 + 33*u^
#            23 + 24*u^21 + 17*u^19 + 11*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9,
#          u^48 + 2*u^46 + 4*u^44 + 5*u^42 + 6*u^40 + 6*u^38 + 5*u^36 + 4*u^
#            34 + 2*u^32 + u^30, u^48 + 3*u^46 + 7*u^44 + 13*u^42 + 20*u^40 + 
#            27*u^38 + 32*u^36 + 34*u^34 + 32*u^32 + 27*u^30 + 20*u^28 + 13*u^
#            26 + 7*u^24 + 3*u^22 + u^20, u^41 + 5*u^39 + 14*u^37 + 26*u^35 + 
#            38*u^33 + 46*u^31 + 49*u^29 + 46*u^27 + 40*u^25 + 32*u^23 + 24*u^
#            21 + 17*u^19 + 11*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          2*u^40 + 7*u^38 + 15*u^36 + 22*u^34 + 26*u^32 + 25*u^30 + 20*u^28 + 
#            13*u^26 + 7*u^24 + 3*u^22 + u^20, 
#          u^39 + 3*u^37 + 5*u^35 + 5*u^33 + 4*u^31 + 2*u^29 + u^27, 
#          2*u^37 + 9*u^35 + 22*u^33 + 34*u^31 + 43*u^29 + 44*u^27 + 40*u^25 + 
#            32*u^23 + 24*u^21 + 17*u^19 + 11*u^17 + 7*u^15 + 4*u^13 + 2*u^
#            11 + u^9, u^37 + 5*u^35 + 14*u^33 + 26*u^31 + 39*u^29 + 49*u^27 + 
#            51*u^25 + 46*u^23 + 36*u^21 + 24*u^19 + 14*u^17 + 8*u^15 + 4*u^
#            13 + 2*u^11 + u^9, 2*u^33 + 7*u^31 + 15*u^29 + 22*u^27 + 26*u^
#            25 + 25*u^23 + 21*u^21 + 16*u^19 + 11*u^17 + 7*u^15 + 4*u^13 + 
#            2*u^11 + u^9, u^32 + 2*u^30 + 3*u^28 + 2*u^26 + u^24, 
#          u^33 + 4*u^31 + 9*u^29 + 13*u^27 + 14*u^25 + 11*u^23 + 7*u^21 + 3*u^
#            19 + u^17, 3*u^29 + 9*u^27 + 16*u^25 + 20*u^23 + 20*u^21 + 16*u^
#            19 + 11*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          u^30 + 3*u^28 + 5*u^26 + 5*u^24 + 3*u^22 + u^20, 
#          2*u^29 + 7*u^27 + 12*u^25 + 14*u^23 + 12*u^21 + 7*u^19 + 3*u^17 + u^
#            15, 0*u^0, u^27 + u^25 + u^23, u^29 + 3*u^27 + 6*u^25 + 9*u^23 + 
#            11*u^21 + 12*u^19 + 11*u^17 + 9*u^15 + 6*u^13 + 3*u^11 + u^9, 
#          3*u^25 + 9*u^23 + 13*u^21 + 13*u^19 + 10*u^17 + 7*u^15 + 4*u^13 + 
#            2*u^11 + u^9, 4*u^25 + 9*u^23 + 10*u^21 + 7*u^19 + 3*u^17 + u^15, 
#          2*u^25 + 4*u^23 + 3*u^21 + u^19, 0*u^0, u^25 + u^23, 
#          3*u^25 + 12*u^23 + 19*u^21 + 18*u^19 + 13*u^17 + 8*u^15 + 4*u^13 + 
#            2*u^11 + u^9, u^25 + 6*u^23 + 13*u^21 + 16*u^19 + 13*u^17 + 8*u^
#            15 + 4*u^13 + 2*u^11 + u^9, u^23 + u^21, 
#          2*u^24 + 7*u^22 + 11*u^20 + 10*u^18 + 5*u^16 + u^14, 
#          u^23 + 6*u^21 + 10*u^19 + 10*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          u^22 + 2*u^20 + u^18, 2*u^22 + 5*u^20 + 5*u^18 + 3*u^16 + u^14, 
#          3*u^21 + 9*u^19 + 10*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          u^22 + 6*u^20 + 8*u^18 + 5*u^16 + u^14, 
#          u^21 + 6*u^19 + 10*u^17 + 9*u^15 + 6*u^13 + 3*u^11 + u^9, 
#          2*u^20 + 5*u^18 + 4*u^16 + u^14, u^20 + 3*u^18 + 2*u^16, 
#          4*u^19 + 8*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          u^20 + 3*u^18 + 2*u^16, u^19 + u^17, 
#          2*u^19 + 8*u^17 + 7*u^15 + 4*u^13 + 2*u^11 + u^9, 
#          u^19 + 8*u^17 + 9*u^15 + 6*u^13 + 3*u^11 + u^9, 
#          3*u^17 + 4*u^15 + 3*u^13 + 2*u^11 + u^9, 
#          4*u^17 + 6*u^15 + 7*u^13 + 3*u^11 + u^9, 
#          3*u^17 + 5*u^15 + 6*u^13 + 3*u^11 + u^9, 
#          u^18 + 2*u^16 + 2*u^14 + u^12, u^17 + 2*u^15 + 2*u^13 + u^11, 
#          u^17 + 3*u^15 + 4*u^13 + 3*u^11 + u^9, 2*u^13 + 2*u^11 + u^9, 
#          u^14 + u^12, 0*u^0, u^13, 2*u^13 + u^11, u^13, u^14, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          2*u^13 + 3*u^11 + u^9, 2*u^13 + 3*u^11 + u^9, u^12, u^13 + u^11, 
#          u^13 + 3*u^11 + u^9, 2*u^11 + u^9, 2*u^11 + u^9, u^11 + u^9, u^10, 
#          0*u^0, u^11 + u^9, 0*u^0, u^10, 0*u^0, u^9, u^10, u^9, u^9, 0*u^0, 
#          u^9, 0*u^0, 0*u^0, u^9, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^88 + 2*u^84 + 3*u^82 + 4*u^80 + 7*u^78 + 12*u^76 + 11*u^74 + 21*u^
#            72 + 24*u^70 + 28*u^68 + 38*u^66 + 45*u^64 + 45*u^62 + 62*u^60 + 
#            61*u^58 + 66*u^56 + 77*u^54 + 76*u^52 + 74*u^50 + 86*u^48 + 74*u^
#            46 + 76*u^44 + 77*u^42 + 66*u^40 + 61*u^38 + 62*u^36 + 45*u^34 + 
#            45*u^32 + 38*u^30 + 28*u^28 + 24*u^26 + 21*u^24 + 11*u^22 + 12*u^
#            20 + 7*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          2*u^66 + 3*u^64 + 5*u^62 + 12*u^60 + 15*u^58 + 22*u^56 + 33*u^54 + 
#            38*u^52 + 45*u^50 + 58*u^48 + 56*u^46 + 63*u^44 + 67*u^42 + 61*u^
#            40 + 59*u^38 + 60*u^36 + 45*u^34 + 45*u^32 + 38*u^30 + 28*u^28 + 
#            24*u^26 + 21*u^24 + 11*u^22 + 12*u^20 + 7*u^18 + 4*u^16 + 3*u^
#            14 + 2*u^12 + u^8, u^56 + 2*u^54 + 5*u^52 + 9*u^50 + 18*u^48 + 
#            21*u^46 + 34*u^44 + 39*u^42 + 43*u^40 + 46*u^38 + 50*u^36 + 40*u^
#            34 + 43*u^32 + 36*u^30 + 28*u^28 + 24*u^26 + 21*u^24 + 11*u^22 + 
#            12*u^20 + 7*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          u^48 + 2*u^46 + 5*u^44 + 10*u^42 + 16*u^40 + 19*u^38 + 30*u^36 + 
#            27*u^34 + 32*u^32 + 31*u^30 + 26*u^28 + 22*u^26 + 21*u^24 + 11*u^
#            22 + 12*u^20 + 7*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          u^47 + 2*u^43 + 2*u^41 + 2*u^39 + 2*u^37 + 2*u^35 + u^31, 
#          u^45 + 5*u^43 + 4*u^41 + 12*u^39 + 13*u^37 + 13*u^35 + 17*u^33 + 
#            14*u^31 + 10*u^29 + 10*u^27 + 5*u^25 + 2*u^23 + 2*u^21, 
#          u^40 + 2*u^38 + 9*u^36 + 12*u^34 + 19*u^32 + 21*u^30 + 21*u^28 + 
#            20*u^26 + 19*u^24 + 11*u^22 + 12*u^20 + 7*u^18 + 4*u^16 + 3*u^
#            14 + 2*u^12 + u^8, 2*u^39 + 4*u^37 + 7*u^35 + 13*u^33 + 11*u^31 + 
#            10*u^29 + 10*u^27 + 5*u^25 + 2*u^23 + 2*u^21, 
#          u^36 + 2*u^34 + 2*u^32 + u^28, 2*u^36 + 4*u^34 + 13*u^32 + 17*u^
#            30 + 18*u^28 + 20*u^26 + 19*u^24 + 11*u^22 + 12*u^20 + 7*u^18 + 
#            4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          u^36 + 3*u^34 + 9*u^32 + 13*u^30 + 20*u^28 + 23*u^26 + 24*u^24 + 
#            19*u^22 + 17*u^20 + 9*u^18 + 6*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          u^32 + 5*u^30 + 8*u^28 + 9*u^26 + 14*u^24 + 9*u^22 + 10*u^20 + 7*u^
#            18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, u^29 + u^25, 
#          u^34 + 4*u^30 + 5*u^28 + 6*u^26 + 7*u^24 + 5*u^22 + 2*u^20 + 2*u^18,
#          u^30 + 2*u^28 + 5*u^26 + 10*u^24 + 8*u^22 + 10*u^20 + 7*u^18 + 4*u^
#            16 + 3*u^14 + 2*u^12 + u^8, 2*u^27 + 2*u^25 + u^23 + u^21, 
#          2*u^28 + 4*u^26 + 5*u^24 + 8*u^22 + 5*u^20 + 2*u^18 + 2*u^16, 
#          0*u^0, u^24, u^28 + u^26 + 4*u^24 + 4*u^22 + 5*u^20 + 6*u^18 + 5*u^
#            16 + 3*u^14 + 4*u^12 + u^8, 3*u^24 + 3*u^22 + 8*u^20 + 5*u^18 + 
#            4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          3*u^24 + 5*u^22 + 5*u^20 + 2*u^18 + 2*u^16, u^24 + u^22 + u^20, 
#          0*u^0, u^24, 3*u^24 + 7*u^22 + 10*u^20 + 7*u^18 + 6*u^16 + 3*u^14 + 
#            2*u^12 + u^8, 2*u^24 + 4*u^22 + 8*u^20 + 7*u^18 + 6*u^16 + 3*u^
#            14 + 2*u^12 + u^8, u^22, 2*u^23 + 4*u^21 + 5*u^19 + 4*u^17 + 2*u^
#            15, u^22 + 4*u^20 + 5*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          u^19, 3*u^21 + 2*u^19 + 2*u^17 + 2*u^15, 
#          3*u^20 + 5*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          2*u^21 + 2*u^19 + 4*u^17 + 2*u^15, 
#          u^20 + 4*u^18 + 5*u^16 + 3*u^14 + 4*u^12 + u^8, 
#          u^21 + u^19 + 3*u^17 + 2*u^15, u^17, 
#          u^20 + 2*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 2*u^17, u^18, 
#          2*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + u^8, 
#          u^18 + 5*u^16 + 3*u^14 + 4*u^12 + u^8, 2*u^16 + u^14 + 2*u^12 + u^8,
#          u^18 + 3*u^16 + 4*u^14 + 4*u^12 + u^8, 2*u^16 + u^14 + 4*u^12 + u^8,
#          u^17 + u^13, u^16 + 2*u^12, u^16 + 2*u^14 + 2*u^12 + u^8, 
#          2*u^12 + u^8, u^13, 0*u^0, 0*u^0, 2*u^12, u^14, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          2*u^12 + u^8, 3*u^12 + u^8, 0*u^0, 2*u^12, 2*u^12 + u^8, 
#          u^12 + u^8, u^8, u^12 + u^8, 0*u^0, 0*u^0, u^8, 0*u^0, u^11, 0*u^0, 
#          u^8, 0*u^0, u^8, u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^8, u^8, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^86 + 2*u^84 + 2*u^82 + 6*u^80 + 7*u^78 + 10*u^76 + 17*u^74 + 20*u^
#            72 + 25*u^70 + 37*u^68 + 38*u^66 + 48*u^64 + 60*u^62 + 61*u^60 + 
#            71*u^58 + 82*u^56 + 77*u^54 + 89*u^52 + 92*u^50 + 85*u^48 + 92*u^
#            46 + 89*u^44 + 77*u^42 + 82*u^40 + 71*u^38 + 61*u^36 + 60*u^34 + 
#            48*u^32 + 38*u^30 + 37*u^28 + 25*u^26 + 20*u^24 + 17*u^22 + 10*u^
#            20 + 7*u^18 + 6*u^16 + 2*u^14 + 2*u^12 + u^10, 
#          u^68 + u^66 + 4*u^64 + 8*u^62 + 11*u^60 + 19*u^58 + 28*u^56 + 32*u^
#            54 + 46*u^52 + 54*u^50 + 58*u^48 + 69*u^46 + 72*u^44 + 68*u^42 + 
#            75*u^40 + 67*u^38 + 60*u^36 + 59*u^34 + 48*u^32 + 38*u^30 + 37*u^
#            28 + 25*u^26 + 20*u^24 + 17*u^22 + 10*u^20 + 7*u^18 + 6*u^16 + 
#            2*u^14 + 2*u^12 + u^10, u^56 + 2*u^54 + 7*u^52 + 11*u^50 + 18*u^
#            48 + 28*u^46 + 35*u^44 + 41*u^42 + 52*u^40 + 50*u^38 + 51*u^36 + 
#            52*u^34 + 44*u^32 + 37*u^30 + 36*u^28 + 25*u^26 + 20*u^24 + 17*u^
#            22 + 10*u^20 + 7*u^18 + 6*u^16 + 2*u^14 + 2*u^12 + u^10, 
#          u^48 + 3*u^46 + 7*u^44 + 11*u^42 + 18*u^40 + 25*u^38 + 28*u^36 + 
#            34*u^34 + 35*u^32 + 30*u^30 + 32*u^28 + 24*u^26 + 19*u^24 + 17*u^
#            22 + 10*u^20 + 7*u^18 + 6*u^16 + 2*u^14 + 2*u^12 + u^10, 
#          u^45 + u^43 + 2*u^39 + u^35 + u^33, 2*u^47 + 3*u^45 + 5*u^43 + 11*u^
#            41 + 12*u^39 + 15*u^37 + 21*u^35 + 16*u^33 + 17*u^31 + 15*u^29 + 
#            8*u^27 + 7*u^25 + 4*u^23 + u^21 + u^19, 
#          u^40 + 4*u^38 + 8*u^36 + 14*u^34 + 19*u^32 + 21*u^30 + 25*u^28 + 
#            20*u^26 + 18*u^24 + 16*u^22 + 10*u^20 + 7*u^18 + 6*u^16 + 2*u^
#            14 + 2*u^12 + u^10, u^41 + 2*u^39 + 6*u^37 + 11*u^35 + 11*u^33 + 
#            15*u^31 + 14*u^29 + 8*u^27 + 7*u^25 + 4*u^23 + u^21 + u^19, 
#          u^36 + u^32 + u^30, 2*u^36 + 7*u^34 + 11*u^32 + 16*u^30 + 23*u^28 + 
#            19*u^26 + 18*u^24 + 16*u^22 + 10*u^20 + 7*u^18 + 6*u^16 + 2*u^
#            14 + 2*u^12 + u^10, 2*u^36 + 5*u^34 + 11*u^32 + 18*u^30 + 25*u^
#            28 + 27*u^26 + 28*u^24 + 23*u^22 + 17*u^20 + 11*u^18 + 7*u^16 + 
#            3*u^14 + 2*u^12 + u^10, 2*u^32 + 3*u^30 + 8*u^28 + 11*u^26 + 11*u^
#            24 + 12*u^22 + 9*u^20 + 6*u^18 + 6*u^16 + 2*u^14 + 2*u^12 + u^10, 
#          u^27, 2*u^32 + 3*u^30 + 5*u^28 + 9*u^26 + 6*u^24 + 6*u^22 + 4*u^
#            20 + u^18 + u^16, 3*u^28 + 6*u^26 + 8*u^24 + 11*u^22 + 9*u^20 + 
#            6*u^18 + 6*u^16 + 2*u^14 + 2*u^12 + u^10, u^29 + u^25 + u^23, 
#          u^30 + 3*u^28 + 6*u^26 + 10*u^24 + 7*u^22 + 7*u^20 + 4*u^18 + u^
#            16 + u^14, 0*u^0, u^28 + 2*u^26 + u^24 + 2*u^22, 
#          u^28 + 3*u^26 + 5*u^24 + 6*u^22 + 9*u^20 + 6*u^18 + 8*u^16 + 5*u^
#            14 + 3*u^12 + 2*u^10, 2*u^24 + 5*u^22 + 5*u^20 + 5*u^18 + 5*u^
#            16 + 2*u^14 + 2*u^12 + u^10, u^26 + 4*u^24 + 5*u^22 + 6*u^20 + 
#            4*u^18 + u^16 + u^14, u^22, 0*u^0, u^26 + u^24 + 2*u^22, 
#          3*u^24 + 6*u^22 + 9*u^20 + 8*u^18 + 6*u^16 + 3*u^14 + 2*u^12 + u^10,
#          2*u^24 + 5*u^22 + 8*u^20 + 8*u^18 + 6*u^16 + 3*u^14 + 2*u^12 + u^10,
#          0*u^0, u^25 + 3*u^23 + 6*u^21 + 7*u^19 + 5*u^17 + 2*u^15 + u^13, 
#          u^22 + 3*u^20 + 4*u^18 + 5*u^16 + 2*u^14 + 2*u^12 + u^10, 0*u^0, 
#          u^23 + u^21 + 4*u^19 + 3*u^17 + u^15 + u^13, 
#          u^22 + 2*u^20 + 4*u^18 + 5*u^16 + 2*u^14 + 2*u^12 + u^10, 
#          u^21 + 5*u^19 + 4*u^17 + 2*u^15 + u^13, 
#          2*u^20 + 4*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10, 
#          u^21 + 4*u^19 + 5*u^17 + 2*u^15 + u^13, u^17, 
#          2*u^18 + 4*u^16 + 2*u^14 + 2*u^12 + u^10, u^19 + u^17 + u^15, 
#          0*u^0, u^18 + 4*u^16 + 2*u^14 + 2*u^12 + u^10, 
#          2*u^18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^10, 
#          u^16 + u^14 + u^12 + u^10, u^18 + 4*u^16 + 4*u^14 + 4*u^12 + 2*u^10,
#          u^16 + 3*u^14 + 2*u^12 + 2*u^10, u^15, 
#          u^18 + 2*u^16 + 4*u^14 + u^12 + u^10, 
#          2*u^16 + 2*u^14 + 3*u^12 + 2*u^10, u^12 + u^10, 0*u^0, 0*u^0, 
#          0*u^0, u^14 + u^12 + u^10, u^12, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, 0*u^0, 2*u^12 + 2*u^10, 
#          u^12 + 2*u^10, 0*u^0, u^14 + u^12 + u^10, u^12 + 2*u^10, u^10, 
#          u^10, 2*u^10, u^11 + u^9, 0*u^0, u^10, 0*u^0, u^9, 0*u^0, u^10, 
#          u^9, u^8, 0*u^0, 0*u^0, 0*u^0, u^9, 0*u^0, u^8, 0*u^0, u^8, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^86 + u^84 + 3*u^82 + 3*u^80 + 6*u^78 + 8*u^76 + 11*u^74 + 14*u^72 + 
#            20*u^70 + 21*u^68 + 29*u^66 + 32*u^64 + 37*u^62 + 43*u^60 + 48*u^
#            58 + 48*u^56 + 57*u^54 + 55*u^52 + 58*u^50 + 60*u^48 + 58*u^46 + 
#            55*u^44 + 57*u^42 + 48*u^40 + 48*u^38 + 43*u^36 + 37*u^34 + 32*u^
#            32 + 29*u^30 + 21*u^28 + 20*u^26 + 14*u^24 + 11*u^22 + 8*u^20 + 
#            6*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^66 + 2*u^64 + 4*u^62 + 8*u^60 + 12*u^58 + 16*u^56 + 24*u^54 + 
#            28*u^52 + 35*u^50 + 41*u^48 + 44*u^46 + 46*u^44 + 50*u^42 + 45*u^
#            40 + 46*u^38 + 42*u^36 + 37*u^34 + 32*u^32 + 29*u^30 + 21*u^28 + 
#            20*u^26 + 14*u^24 + 11*u^22 + 8*u^20 + 6*u^18 + 3*u^16 + 3*u^
#            14 + u^12 + u^10, 2*u^54 + 3*u^52 + 8*u^50 + 12*u^48 + 18*u^46 + 
#            23*u^44 + 31*u^42 + 31*u^40 + 37*u^38 + 35*u^36 + 34*u^34 + 30*u^
#            32 + 28*u^30 + 21*u^28 + 20*u^26 + 14*u^24 + 11*u^22 + 8*u^20 + 
#            6*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^46 + 3*u^44 + 6*u^42 + 10*u^40 + 15*u^38 + 19*u^36 + 23*u^34 + 
#            22*u^32 + 25*u^30 + 19*u^28 + 19*u^26 + 14*u^24 + 11*u^22 + 8*u^
#            20 + 6*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^47 + 2*u^45 + 2*u^43 + 4*u^41 + 2*u^39 + 4*u^37 + 2*u^35 + 2*u^
#            33 + u^31, u^45 + 2*u^43 + 4*u^41 + 6*u^39 + 10*u^37 + 9*u^35 + 
#            13*u^33 + 10*u^31 + 8*u^29 + 7*u^27 + 3*u^25 + 2*u^23 + u^21, 
#          2*u^38 + 5*u^36 + 10*u^34 + 13*u^32 + 18*u^30 + 16*u^28 + 17*u^26 + 
#            13*u^24 + 11*u^22 + 8*u^20 + 6*u^18 + 3*u^16 + 3*u^14 + u^12 + u^
#            10, u^39 + 3*u^37 + 5*u^35 + 9*u^33 + 9*u^31 + 8*u^29 + 7*u^27 + 
#            3*u^25 + 2*u^23 + u^21, u^38 + u^36 + 3*u^34 + 2*u^32 + 2*u^
#            30 + u^28, u^36 + 4*u^34 + 9*u^32 + 14*u^30 + 15*u^28 + 17*u^26 + 
#            13*u^24 + 11*u^22 + 8*u^20 + 6*u^18 + 3*u^16 + 3*u^14 + u^12 + u^
#            10, 2*u^34 + 4*u^32 + 9*u^30 + 13*u^28 + 17*u^26 + 16*u^24 + 16*u^
#            22 + 11*u^20 + 8*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^32 + 4*u^30 + 6*u^28 + 9*u^26 + 10*u^24 + 9*u^22 + 7*u^20 + 6*u^
#            18 + 3*u^16 + 3*u^14 + u^12 + u^10, u^29 + u^27 + u^25, 
#          u^32 + 2*u^30 + 5*u^28 + 4*u^26 + 6*u^24 + 3*u^22 + 2*u^20 + u^18, 
#          u^28 + 3*u^26 + 6*u^24 + 7*u^22 + 7*u^20 + 6*u^18 + 3*u^16 + 3*u^
#            14 + u^12 + u^10, u^29 + 3*u^27 + 3*u^25 + 3*u^23 + 2*u^21, 
#          2*u^26 + 3*u^24 + 5*u^22 + 3*u^20 + 2*u^18 + u^16, 0*u^0, u^24, 
#          u^26 + u^24 + 3*u^22 + 2*u^20 + 4*u^18 + 2*u^16 + 3*u^14 + 2*u^
#            12 + u^10, u^24 + 4*u^22 + 5*u^20 + 5*u^18 + 3*u^16 + 3*u^14 + u^
#            12 + u^10, u^24 + 4*u^22 + 3*u^20 + 2*u^18 + u^16, 
#          2*u^24 + 2*u^22 + 2*u^20, 0*u^0, u^24, 
#          2*u^24 + 7*u^22 + 8*u^20 + 7*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          3*u^22 + 5*u^20 + 6*u^18 + 4*u^16 + 3*u^14 + u^12 + u^10, 
#          u^24 + u^22 + u^20, u^23 + 3*u^21 + 4*u^19 + 3*u^17 + u^15, 
#          u^22 + 3*u^20 + 5*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^21 + u^19, 2*u^21 + 2*u^19 + 2*u^17 + u^15, 
#          2*u^20 + 4*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^21 + 3*u^19 + 3*u^17 + u^15, 
#          2*u^18 + 2*u^16 + 3*u^14 + 2*u^12 + u^10, u^19 + u^17 + u^15, 
#          u^19 + u^17 + u^15, 3*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^19 + u^17, u^18 + u^16, 3*u^18 + 3*u^16 + 3*u^14 + u^12 + u^10, 
#          u^18 + 2*u^16 + 3*u^14 + 2*u^12 + u^10, 
#          u^18 + u^16 + 2*u^14 + u^12 + u^10, u^16 + 2*u^14 + 2*u^12 + u^10, 
#          u^16 + 2*u^14 + 3*u^12 + u^10, 2*u^17 + 2*u^15 + 2*u^13, u^12, 
#          u^14 + u^12 + u^10, u^12 + u^10, u^13, 0*u^0, u^12, u^12, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^14, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^12 + u^10, 3*u^12 + u^10, 0*u^0, u^12, u^12 + u^10, 
#          u^12 + u^10, u^10, 0*u^0, 0*u^0, u^11, u^10, 0*u^0, u^11, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^9, u^8, 0*u^0, 0*u^0, u^8, 0*u^0, 
#          0*u^0, u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^84 + u^80 + u^78 + u^76 + 2*u^74 + 3*u^72 + 2*u^70 + 5*u^68 + 4*u^
#            66 + 5*u^64 + 7*u^62 + 7*u^60 + 7*u^58 + 10*u^56 + 7*u^54 + 10*u^
#            52 + 10*u^50 + 9*u^48 + 10*u^46 + 10*u^44 + 7*u^42 + 10*u^40 + 
#            7*u^38 + 7*u^36 + 7*u^34 + 5*u^32 + 4*u^30 + 5*u^28 + 2*u^26 + 
#            3*u^24 + 2*u^22 + u^20 + u^18 + u^16 + u^12, 
#          u^62 + u^60 + 2*u^58 + 3*u^56 + 3*u^54 + 5*u^52 + 6*u^50 + 6*u^48 + 
#            8*u^46 + 8*u^44 + 7*u^42 + 9*u^40 + 7*u^38 + 7*u^36 + 7*u^34 + 
#            5*u^32 + 4*u^30 + 5*u^28 + 2*u^26 + 3*u^24 + 2*u^22 + u^20 + u^
#            18 + u^16 + u^12, u^52 + u^50 + 2*u^48 + 3*u^46 + 4*u^44 + 4*u^
#            42 + 7*u^40 + 5*u^38 + 7*u^36 + 6*u^34 + 5*u^32 + 4*u^30 + 5*u^
#            28 + 2*u^26 + 3*u^24 + 2*u^22 + u^20 + u^18 + u^16 + u^12, 
#          u^42 + u^40 + 2*u^38 + 4*u^36 + 3*u^34 + 5*u^32 + 3*u^30 + 5*u^28 + 
#            2*u^26 + 3*u^24 + 2*u^22 + u^20 + u^18 + u^16 + u^12, 
#          u^45 + u^43 + 2*u^39 + u^35 + u^33, u^39 + 3*u^35 + u^33 + 2*u^31 + 
#            2*u^29 + u^25, u^36 + u^34 + 3*u^32 + 3*u^30 + 4*u^28 + 2*u^26 + 
#            3*u^24 + 2*u^22 + u^20 + u^18 + u^16 + u^12, 
#          u^35 + u^33 + 2*u^31 + 2*u^29 + u^25, u^36 + u^32 + u^30, 
#          u^34 + u^32 + 3*u^30 + 4*u^28 + 2*u^26 + 3*u^24 + 2*u^22 + u^20 + u^
#            18 + u^16 + u^12, u^30 + 2*u^28 + 2*u^26 + 3*u^24 + 2*u^22 + 2*u^
#            20 + u^18 + u^16 + u^12, u^30 + u^28 + 2*u^26 + 2*u^24 + 2*u^
#            22 + u^20 + u^18 + u^16 + u^12, u^27, u^30 + 2*u^26 + u^22, 
#          u^26 + 2*u^22 + u^20 + u^18 + u^16 + u^12, u^25 + 2*u^23, u^20, 
#          u^28, 0*u^0, u^12, u^22 + u^20 + u^18 + u^16 + u^12, u^20, 2*u^22, 
#          0*u^0, 0*u^0, 2*u^22 + 2*u^20 + u^18 + u^16 + u^12, 
#          u^20 + u^18 + u^16 + u^12, u^22, u^19, u^20 + u^18 + u^16 + u^12, 
#          u^21, u^19, u^18 + u^16 + u^12, u^19, u^12, 0*u^0, 0*u^0, 
#          u^18 + u^16 + u^12, 0*u^0, u^18, u^18 + u^16 + u^12, u^12, 
#          u^16 + u^12, u^12, u^12, u^15, 0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 
#          0*u^0, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^9, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^82 + 2*u^78 + 3*u^76 + 2*u^74 + 5*u^72 + 7*u^70 + 5*u^68 + 11*u^
#            66 + 11*u^64 + 10*u^62 + 17*u^60 + 16*u^58 + 14*u^56 + 23*u^54 + 
#            18*u^52 + 18*u^50 + 24*u^48 + 18*u^46 + 18*u^44 + 23*u^42 + 14*u^
#            40 + 16*u^38 + 17*u^36 + 10*u^34 + 11*u^32 + 11*u^30 + 5*u^28 + 
#            7*u^26 + 5*u^24 + 2*u^22 + 3*u^20 + 2*u^18 + u^14, 
#          u^66 + u^64 + u^62 + 4*u^60 + 4*u^58 + 5*u^56 + 10*u^54 + 9*u^52 + 
#            11*u^50 + 16*u^48 + 13*u^46 + 15*u^44 + 19*u^42 + 13*u^40 + 15*u^
#            38 + 16*u^36 + 10*u^34 + 11*u^32 + 11*u^30 + 5*u^28 + 7*u^26 + 
#            5*u^24 + 2*u^22 + 3*u^20 + 2*u^18 + u^14, 
#          u^54 + u^52 + 3*u^50 + 5*u^48 + 5*u^46 + 8*u^44 + 11*u^42 + 8*u^
#            40 + 12*u^38 + 12*u^36 + 9*u^34 + 10*u^32 + 10*u^30 + 5*u^28 + 
#            7*u^26 + 5*u^24 + 2*u^22 + 3*u^20 + 2*u^18 + u^14, 
#          u^46 + u^44 + 4*u^42 + 3*u^40 + 5*u^38 + 7*u^36 + 6*u^34 + 6*u^32 + 
#            9*u^30 + 4*u^28 + 6*u^26 + 5*u^24 + 2*u^22 + 3*u^20 + 2*u^18 + u^
#            14, 0*u^0, u^49 + 2*u^45 + 3*u^43 + 2*u^41 + 5*u^39 + 5*u^37 + 
#            3*u^35 + 6*u^33 + 4*u^31 + 2*u^29 + 4*u^27 + u^25 + u^23 + u^21, 
#          u^38 + 2*u^36 + 2*u^34 + 3*u^32 + 5*u^30 + 3*u^28 + 5*u^26 + 4*u^
#            24 + 2*u^22 + 3*u^20 + 2*u^18 + u^14, 
#          2*u^39 + 2*u^37 + 2*u^35 + 4*u^33 + 3*u^31 + 2*u^29 + 4*u^27 + u^
#            25 + u^23 + u^21, 0*u^0, u^36 + u^34 + 2*u^32 + 4*u^30 + 2*u^28 + 
#            5*u^26 + 4*u^24 + 2*u^22 + 3*u^20 + 2*u^18 + u^14, 
#          u^36 + 2*u^34 + 3*u^32 + 5*u^30 + 6*u^28 + 7*u^26 + 6*u^24 + 6*u^
#            22 + 4*u^20 + 3*u^18 + u^16 + u^14, 
#          u^30 + u^28 + u^26 + 3*u^24 + u^22 + 2*u^20 + 2*u^18 + u^14, 0*u^0, 
#          u^30 + 2*u^28 + 3*u^24 + u^22 + u^20 + u^18, 
#          u^28 + 2*u^24 + u^22 + 2*u^20 + 2*u^18 + u^14, 0*u^0, 
#          2*u^28 + 2*u^26 + u^24 + 4*u^22 + u^20 + u^18 + u^16, 0*u^0, 
#          u^30 + u^28 + u^26 + 2*u^24 + u^20, 
#          u^28 + 2*u^24 + 2*u^22 + 2*u^20 + 4*u^18 + u^16 + 2*u^14 + u^12, 
#          u^20 + u^18 + u^14, u^26 + 2*u^22 + u^20 + u^18 + u^16, 0*u^0, 
#          0*u^0, u^26 + u^24 + u^20, u^22 + u^20 + 2*u^18 + u^16 + u^14, 
#          u^22 + u^20 + 2*u^18 + u^16 + u^14, 0*u^0, 
#          u^25 + u^23 + 2*u^21 + 2*u^19 + 2*u^17 + u^15, u^18 + u^14, 0*u^0, 
#          u^21 + u^17 + u^15, u^18 + u^14, u^21 + 2*u^17 + u^15, 
#          2*u^18 + u^16 + 2*u^14 + u^12, u^21 + u^19 + 2*u^17 + 2*u^15, 
#          0*u^0, u^14, u^17, 0*u^0, u^14, u^18 + u^16 + 2*u^14 + u^12, 0*u^0, 
#          2*u^14 + u^12, u^12, 0*u^0, u^18 + 2*u^16 + u^14 + 2*u^12, 
#          2*u^14 + u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, u^15, 0*u^0, 
#          u^12, 0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^11 + u^9, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^9, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^9, u^9, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^8, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^83 + u^81 + 2*u^79 + 5*u^77 + 6*u^75 + 9*u^73 + 15*u^71 + 17*u^69 + 
#            23*u^67 + 32*u^65 + 34*u^63 + 43*u^61 + 53*u^59 + 54*u^57 + 64*u^
#            55 + 72*u^53 + 70*u^51 + 79*u^49 + 82*u^47 + 76*u^45 + 82*u^43 + 
#            79*u^41 + 70*u^39 + 72*u^37 + 64*u^35 + 54*u^33 + 53*u^31 + 43*u^
#            29 + 34*u^27 + 32*u^25 + 23*u^23 + 17*u^21 + 15*u^19 + 9*u^17 + 
#            6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          u^65 + u^63 + 4*u^61 + 8*u^59 + 11*u^57 + 19*u^55 + 27*u^53 + 32*u^
#            51 + 44*u^49 + 52*u^47 + 55*u^45 + 65*u^43 + 67*u^41 + 64*u^39 + 
#            68*u^37 + 62*u^35 + 54*u^33 + 53*u^31 + 43*u^29 + 34*u^27 + 32*u^
#            25 + 23*u^23 + 17*u^21 + 15*u^19 + 9*u^17 + 6*u^15 + 5*u^13 + 2*u^
#            11 + u^9 + u^7, u^53 + 3*u^51 + 8*u^49 + 13*u^47 + 20*u^45 + 31*u^
#            43 + 37*u^41 + 43*u^39 + 51*u^37 + 50*u^35 + 48*u^33 + 49*u^31 + 
#            41*u^29 + 34*u^27 + 32*u^25 + 23*u^23 + 17*u^21 + 15*u^19 + 9*u^
#            17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          u^47 + u^45 + 5*u^43 + 10*u^41 + 14*u^39 + 22*u^37 + 29*u^35 + 30*u^
#            33 + 36*u^31 + 35*u^29 + 30*u^27 + 30*u^25 + 23*u^23 + 17*u^21 + 
#            15*u^19 + 9*u^17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          u^42 + u^40 + u^38 + 2*u^36 + u^34 + u^32 + u^30, 
#          2*u^44 + 3*u^42 + 6*u^40 + 11*u^38 + 12*u^36 + 15*u^34 + 18*u^32 + 
#            14*u^30 + 14*u^28 + 11*u^26 + 6*u^24 + 4*u^22 + 2*u^20, 
#          3*u^37 + 8*u^35 + 13*u^33 + 20*u^31 + 23*u^29 + 24*u^27 + 26*u^25 + 
#            21*u^23 + 17*u^21 + 15*u^19 + 9*u^17 + 6*u^15 + 5*u^13 + 2*u^
#            11 + u^9 + u^7, u^38 + 3*u^36 + 8*u^34 + 12*u^32 + 11*u^30 + 13*u^
#            28 + 11*u^26 + 6*u^24 + 4*u^22 + 2*u^20, 
#          u^35 + 2*u^33 + u^31 + u^29 + u^27, 
#          u^35 + 5*u^33 + 13*u^31 + 17*u^29 + 21*u^27 + 25*u^25 + 21*u^23 + 
#            17*u^21 + 15*u^19 + 9*u^17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7,
#          u^35 + 4*u^33 + 10*u^31 + 16*u^29 + 24*u^27 + 29*u^25 + 29*u^23 + 
#            26*u^21 + 21*u^19 + 13*u^17 + 8*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7
#            , u^31 + 5*u^29 + 7*u^27 + 13*u^25 + 15*u^23 + 13*u^21 + 13*u^
#            19 + 9*u^17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          u^30 + u^28 + u^26 + u^24, 3*u^29 + 4*u^27 + 7*u^25 + 8*u^23 + 5*u^
#            21 + 4*u^19 + 2*u^17, u^29 + 2*u^27 + 8*u^25 + 11*u^23 + 11*u^
#            21 + 13*u^19 + 9*u^17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          u^26 + u^24 + u^22 + u^20, 3*u^27 + 5*u^25 + 8*u^23 + 9*u^21 + 6*u^
#            19 + 4*u^17 + 2*u^15, 0*u^0, u^25 + u^23, 
#          u^29 + u^27 + 4*u^25 + 7*u^23 + 7*u^21 + 10*u^19 + 10*u^17 + 7*u^
#            15 + 7*u^13 + 4*u^11 + u^9 + u^7, 
#          u^25 + 3*u^23 + 6*u^21 + 9*u^19 + 7*u^17 + 6*u^15 + 5*u^13 + 2*u^
#            11 + u^9 + u^7, 4*u^23 + 6*u^21 + 5*u^19 + 4*u^17 + 2*u^15, 
#          u^23 + u^21 + u^19, 0*u^0, u^23, 3*u^23 + 8*u^21 + 11*u^19 + 10*u^
#            17 + 8*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          2*u^23 + 6*u^21 + 10*u^19 + 10*u^17 + 8*u^15 + 5*u^13 + 2*u^11 + u^
#            9 + u^7, 0*u^0, 3*u^22 + 6*u^20 + 7*u^18 + 6*u^16 + 3*u^14, 
#          u^21 + 5*u^19 + 6*u^17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          u^20 + u^18, 2*u^20 + 2*u^18 + 3*u^16 + 2*u^14, 
#          5*u^19 + 6*u^17 + 6*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          2*u^20 + 4*u^18 + 5*u^16 + 3*u^14, 
#          4*u^19 + 7*u^17 + 7*u^15 + 7*u^13 + 4*u^11 + u^9 + u^7, 
#          u^20 + 3*u^18 + 4*u^16 + 3*u^14, u^18 + u^16, 
#          u^19 + 3*u^17 + 5*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          2*u^18 + 2*u^16 + u^14, 0*u^0, 
#          3*u^17 + 5*u^15 + 5*u^13 + 2*u^11 + u^9 + u^7, 
#          4*u^17 + 6*u^15 + 7*u^13 + 4*u^11 + u^9 + u^7, 
#          u^17 + 2*u^15 + 3*u^13 + 2*u^11 + u^9 + u^7, 
#          2*u^17 + 4*u^15 + 7*u^13 + 5*u^11 + u^9 + u^7, 
#          u^17 + 2*u^15 + 4*u^13 + 4*u^11 + u^9 + u^7, u^12, 
#          u^17 + u^15 + 2*u^13 + 2*u^11, u^17 + 3*u^15 + 5*u^13 + 4*u^11 + u^
#            9 + u^7, u^13 + 2*u^11 + u^9 + u^7, u^12, 0*u^0, 0*u^0, 
#          u^13 + 2*u^11, u^13 + u^11, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13 + 4*u^11 + u^9 + u^7,
#          3*u^11 + u^9 + u^7, u^12, 2*u^11, 3*u^11 + u^9 + u^7, 
#          u^11 + u^9 + u^7, u^11 + u^9 + u^7, 2*u^11 + u^9 + u^7, 2*u^10, 
#          0*u^0, u^9 + u^7, u^11 + u^7, u^10, 0*u^0, u^9 + u^7, u^10, 
#          u^9 + u^7, u^7, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^7, u^7, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^7, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0 ], 
#      [ u^82 + 2*u^80 + 3*u^78 + 5*u^76 + 8*u^74 + 11*u^72 + 15*u^70 + 21*u^
#            68 + 25*u^66 + 31*u^64 + 40*u^62 + 44*u^60 + 51*u^58 + 60*u^56 + 
#            61*u^54 + 69*u^52 + 75*u^50 + 72*u^48 + 78*u^46 + 78*u^44 + 72*u^
#            42 + 75*u^40 + 69*u^38 + 61*u^36 + 60*u^34 + 51*u^32 + 44*u^30 + 
#            40*u^28 + 31*u^26 + 25*u^24 + 21*u^22 + 15*u^20 + 11*u^18 + 8*u^
#            16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^64 + 3*u^62 + 5*u^60 + 9*u^58 + 15*u^56 + 20*u^54 + 29*u^52 + 
#            38*u^50 + 43*u^48 + 53*u^46 + 59*u^44 + 60*u^42 + 66*u^40 + 64*u^
#            38 + 59*u^36 + 59*u^34 + 51*u^32 + 44*u^30 + 40*u^28 + 31*u^26 + 
#            25*u^24 + 21*u^22 + 15*u^20 + 11*u^18 + 8*u^16 + 5*u^14 + 3*u^
#            12 + 2*u^10 + u^8, 2*u^52 + 5*u^50 + 9*u^48 + 17*u^46 + 24*u^44 + 
#            31*u^42 + 41*u^40 + 45*u^38 + 47*u^36 + 50*u^34 + 46*u^32 + 42*u^
#            30 + 39*u^28 + 31*u^26 + 25*u^24 + 21*u^22 + 15*u^20 + 11*u^18 + 
#            8*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^46 + 3*u^44 + 6*u^42 + 11*u^40 + 18*u^38 + 22*u^36 + 29*u^34 + 
#            33*u^32 + 32*u^30 + 34*u^28 + 29*u^26 + 24*u^24 + 21*u^22 + 15*u^
#            20 + 11*u^18 + 8*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^43 + u^41 + 2*u^39 + 2*u^37 + 2*u^35 + 2*u^33 + u^31 + u^29, 
#          u^45 + 2*u^43 + 5*u^41 + 7*u^39 + 10*u^37 + 15*u^35 + 14*u^33 + 
#            16*u^31 + 15*u^29 + 10*u^27 + 9*u^25 + 5*u^23 + 2*u^21 + u^19, 
#          u^38 + 4*u^36 + 10*u^34 + 16*u^32 + 20*u^30 + 25*u^28 + 24*u^26 + 
#            22*u^24 + 20*u^22 + 15*u^20 + 11*u^18 + 8*u^16 + 5*u^14 + 3*u^
#            12 + 2*u^10 + u^8, 2*u^37 + 6*u^35 + 8*u^33 + 12*u^31 + 13*u^29 + 
#            10*u^27 + 9*u^25 + 5*u^23 + 2*u^21 + u^19, 
#          u^36 + 2*u^34 + 2*u^32 + 2*u^30 + u^28 + u^26, 
#          3*u^34 + 8*u^32 + 14*u^30 + 21*u^28 + 22*u^26 + 22*u^24 + 20*u^22 + 
#            15*u^20 + 11*u^18 + 8*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          2*u^34 + 6*u^32 + 12*u^30 + 19*u^28 + 25*u^26 + 28*u^24 + 27*u^22 + 
#            23*u^20 + 16*u^18 + 10*u^16 + 6*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          2*u^30 + 6*u^28 + 11*u^26 + 12*u^24 + 15*u^22 + 13*u^20 + 10*u^18 + 
#            8*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^31 + u^29 + 2*u^27 + u^25 + u^23, 
#          u^30 + 3*u^28 + 7*u^26 + 6*u^24 + 7*u^22 + 5*u^20 + 2*u^18 + u^16, 
#          u^28 + 5*u^26 + 7*u^24 + 12*u^22 + 12*u^20 + 10*u^18 + 8*u^16 + 5*u^
#            14 + 3*u^12 + 2*u^10 + u^8, u^27 + 2*u^25 + 2*u^23 + u^21 + u^19, 
#          u^28 + 3*u^26 + 7*u^24 + 7*u^22 + 8*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          0*u^0, u^26 + u^22, u^28 + 3*u^26 + 4*u^24 + 7*u^22 + 9*u^20 + 8*u^
#            18 + 9*u^16 + 7*u^14 + 4*u^12 + 3*u^10 + u^8, 
#          u^24 + 5*u^22 + 7*u^20 + 8*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^
#            10 + u^8, 2*u^24 + 4*u^22 + 6*u^20 + 5*u^18 + 2*u^16 + u^14, 
#          u^24 + 2*u^22 + u^20 + u^18, 0*u^0, u^22, 
#          u^24 + 5*u^22 + 10*u^20 + 11*u^18 + 9*u^16 + 6*u^14 + 3*u^12 + 2*u^
#            10 + u^8, 3*u^22 + 8*u^20 + 10*u^18 + 9*u^16 + 6*u^14 + 3*u^12 + 
#            2*u^10 + u^8, 0*u^0, u^23 + 4*u^21 + 7*u^19 + 7*u^17 + 4*u^15 + u^
#            13, 3*u^20 + 5*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^21 + u^19 + u^17, 3*u^19 + 3*u^17 + 2*u^15 + u^13, 
#          2*u^20 + 5*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          4*u^19 + 5*u^17 + 4*u^15 + u^13, u^20 + 4*u^18 + 8*u^16 + 7*u^14 + 
#            4*u^12 + 3*u^10 + u^8, 2*u^19 + 4*u^17 + 3*u^15 + u^13, 
#          u^19 + 2*u^17 + u^15, 2*u^18 + 5*u^16 + 5*u^14 + 3*u^12 + 2*u^
#            10 + u^8, u^19 + 2*u^17 + 2*u^15, 0*u^0, 
#          u^18 + 5*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^18 + 6*u^16 + 7*u^14 + 4*u^12 + 3*u^10 + u^8, 
#          2*u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8, 
#          3*u^16 + 5*u^14 + 5*u^12 + 3*u^10 + u^8, 
#          2*u^16 + 4*u^14 + 3*u^12 + 3*u^10 + u^8, u^15 + u^13 + u^11, 
#          u^16 + 2*u^14 + u^12 + u^10, 2*u^16 + 4*u^14 + 4*u^12 + 3*u^10 + u^8
#            , u^12 + 2*u^10 + u^8, u^13 + u^11, 0*u^0, 0*u^0, u^12 + u^10, 
#          u^12, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 2*u^12 + 3*u^10 + u^8, u^12 + 3*u^10 + u^8, 
#          u^13 + u^11, u^12 + u^10, u^12 + 3*u^10 + u^8, 2*u^10 + u^8, 
#          2*u^10 + u^8, 2*u^10 + u^8, u^11 + u^9, 0*u^0, u^10 + u^8, 
#          u^10 + u^8, u^9, 0*u^0, u^10 + u^8, u^9, u^8, u^8, 0*u^0, u^8, 
#          0*u^0, 0*u^0, u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^7, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^81 + 2*u^79 + 2*u^77 + 5*u^75 + 7*u^73 + 8*u^71 + 14*u^69 + 17*u^
#            67 + 19*u^65 + 28*u^63 + 31*u^61 + 34*u^59 + 44*u^57 + 45*u^55 + 
#            48*u^53 + 57*u^51 + 55*u^49 + 56*u^47 + 62*u^45 + 56*u^43 + 55*u^
#            41 + 57*u^39 + 48*u^37 + 45*u^35 + 44*u^33 + 34*u^31 + 31*u^29 + 
#            28*u^27 + 19*u^25 + 17*u^23 + 14*u^21 + 8*u^19 + 7*u^17 + 5*u^
#            15 + 2*u^13 + 2*u^11 + u^9, 2*u^63 + 3*u^61 + 5*u^59 + 10*u^57 + 
#            13*u^55 + 18*u^53 + 27*u^51 + 30*u^49 + 36*u^47 + 44*u^45 + 44*u^
#            43 + 47*u^41 + 51*u^39 + 45*u^37 + 44*u^35 + 43*u^33 + 34*u^31 + 
#            31*u^29 + 28*u^27 + 19*u^25 + 17*u^23 + 14*u^21 + 8*u^19 + 7*u^
#            17 + 5*u^15 + 2*u^13 + 2*u^11 + u^9, 
#          u^53 + 3*u^51 + 5*u^49 + 10*u^47 + 16*u^45 + 20*u^43 + 27*u^41 + 
#            33*u^39 + 33*u^37 + 36*u^35 + 37*u^33 + 31*u^31 + 30*u^29 + 27*u^
#            27 + 19*u^25 + 17*u^23 + 14*u^21 + 8*u^19 + 7*u^17 + 5*u^15 + 2*u^
#            13 + 2*u^11 + u^9, 2*u^45 + 4*u^43 + 6*u^41 + 13*u^39 + 15*u^37 + 
#            19*u^35 + 25*u^33 + 23*u^31 + 24*u^29 + 24*u^27 + 18*u^25 + 16*u^
#            23 + 14*u^21 + 8*u^19 + 7*u^17 + 5*u^15 + 2*u^13 + 2*u^11 + u^9, 
#          0*u^0, u^46 + u^44 + 4*u^42 + 6*u^40 + 7*u^38 + 12*u^36 + 12*u^34 + 
#            12*u^32 + 13*u^30 + 10*u^28 + 7*u^26 + 6*u^24 + 3*u^22 + u^20 + u^
#            18, u^39 + 2*u^37 + 5*u^35 + 10*u^33 + 12*u^31 + 16*u^29 + 18*u^
#            27 + 15*u^25 + 15*u^23 + 13*u^21 + 8*u^19 + 7*u^17 + 5*u^15 + 2*u^
#            13 + 2*u^11 + u^9, u^38 + 4*u^36 + 6*u^34 + 8*u^32 + 11*u^30 + 
#            9*u^28 + 7*u^26 + 6*u^24 + 3*u^22 + u^20 + u^18, 0*u^0, 
#          u^35 + 4*u^33 + 7*u^31 + 12*u^29 + 16*u^27 + 14*u^25 + 15*u^23 + 
#            13*u^21 + 8*u^19 + 7*u^17 + 5*u^15 + 2*u^13 + 2*u^11 + u^9, 
#          u^35 + 4*u^33 + 8*u^31 + 14*u^29 + 19*u^27 + 21*u^25 + 22*u^23 + 
#            19*u^21 + 14*u^19 + 10*u^17 + 6*u^15 + 3*u^13 + 2*u^11 + u^9, 
#          u^31 + 2*u^29 + 7*u^27 + 7*u^25 + 9*u^23 + 10*u^21 + 7*u^19 + 6*u^
#            17 + 5*u^15 + 2*u^13 + 2*u^11 + u^9, 0*u^0, 
#          u^31 + u^29 + 5*u^27 + 5*u^25 + 5*u^23 + 5*u^21 + 3*u^19 + u^17 + u^
#            15, 3*u^27 + 4*u^25 + 7*u^23 + 9*u^21 + 7*u^19 + 6*u^17 + 5*u^
#            15 + 2*u^13 + 2*u^11 + u^9, 0*u^0, 
#          u^29 + 2*u^27 + 6*u^25 + 7*u^23 + 6*u^21 + 6*u^19 + 3*u^17 + u^
#            15 + u^13, 0*u^0, u^27 + u^25 + u^23 + u^21, 
#          2*u^27 + 3*u^25 + 4*u^23 + 8*u^21 + 7*u^19 + 7*u^17 + 8*u^15 + 4*u^
#            13 + 3*u^11 + 2*u^9, 2*u^23 + 4*u^21 + 4*u^19 + 5*u^17 + 4*u^15 + 
#            2*u^13 + 2*u^11 + u^9, u^25 + 3*u^23 + 4*u^21 + 5*u^19 + 3*u^
#            17 + u^15 + u^13, 0*u^0, 0*u^0, u^23 + u^21, 
#          2*u^23 + 5*u^21 + 7*u^19 + 7*u^17 + 5*u^15 + 3*u^13 + 2*u^11 + u^9, 
#          2*u^23 + 5*u^21 + 7*u^19 + 7*u^17 + 5*u^15 + 3*u^13 + 2*u^11 + u^9, 
#          0*u^0, u^24 + 3*u^22 + 5*u^20 + 6*u^18 + 4*u^16 + 2*u^14 + u^12, 
#          u^21 + 2*u^19 + 4*u^17 + 4*u^15 + 2*u^13 + 2*u^11 + u^9, 0*u^0, 
#          u^20 + 3*u^18 + 2*u^16 + u^14 + u^12, 
#          u^21 + 2*u^19 + 4*u^17 + 4*u^15 + 2*u^13 + 2*u^11 + u^9, 
#          u^20 + 4*u^18 + 3*u^16 + 2*u^14 + u^12, 
#          u^21 + 2*u^19 + 5*u^17 + 7*u^15 + 4*u^13 + 3*u^11 + 2*u^9, 
#          u^20 + 4*u^18 + 4*u^16 + 2*u^14 + u^12, 0*u^0, 
#          2*u^17 + 3*u^15 + 2*u^13 + 2*u^11 + u^9, u^18 + u^16 + u^14, 0*u^0, 
#          u^17 + 3*u^15 + 2*u^13 + 2*u^11 + u^9, 
#          2*u^17 + 6*u^15 + 4*u^13 + 3*u^11 + 2*u^9, u^15 + u^13 + u^11 + u^9,
#          u^17 + 5*u^15 + 4*u^13 + 4*u^11 + 2*u^9, 
#          2*u^15 + 2*u^13 + 2*u^11 + 2*u^9, 0*u^0, 
#          u^17 + 3*u^15 + 3*u^13 + u^11 + u^9, 
#          u^17 + 3*u^15 + 3*u^13 + 3*u^11 + 2*u^9, u^11 + u^9, 0*u^0, 0*u^0, 
#          0*u^0, u^13 + u^11 + u^9, u^13 + u^11, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^13, u^13, 0*u^0, 0*u^0, 0*u^0, 
#          u^13 + 2*u^11 + 2*u^9, u^11 + 2*u^9, 0*u^0, u^13 + u^11 + u^9, 
#          u^11 + 2*u^9, u^9, u^9, u^11 + 2*u^9, u^12 + u^10 + u^8, 0*u^0, 
#          u^9, u^9, u^8, 0*u^0, u^9, u^8, u^9 + u^7, 0*u^0, 0*u^0, 0*u^0, 
#          u^8, 0*u^0, u^7, 0*u^0, u^7, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^7, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^81 + u^79 + u^77 + 3*u^75 + 3*u^73 + 4*u^71 + 7*u^69 + 7*u^67 + 9*u^
#            65 + 13*u^63 + 13*u^61 + 16*u^59 + 19*u^57 + 19*u^55 + 22*u^53 + 
#            24*u^51 + 24*u^49 + 25*u^47 + 26*u^45 + 25*u^43 + 24*u^41 + 24*u^
#            39 + 22*u^37 + 19*u^35 + 19*u^33 + 16*u^31 + 13*u^29 + 13*u^27 + 
#            9*u^25 + 7*u^23 + 7*u^21 + 4*u^19 + 3*u^17 + 3*u^15 + u^13 + u^
#            11 + u^9, u^63 + u^61 + 2*u^59 + 4*u^57 + 5*u^55 + 8*u^53 + 11*u^
#            51 + 13*u^49 + 16*u^47 + 19*u^45 + 20*u^43 + 21*u^41 + 22*u^39 + 
#            21*u^37 + 19*u^35 + 19*u^33 + 16*u^31 + 13*u^29 + 13*u^27 + 9*u^
#            25 + 7*u^23 + 7*u^21 + 4*u^19 + 3*u^17 + 3*u^15 + u^13 + u^11 + u^
#            9, u^51 + 2*u^49 + 4*u^47 + 7*u^45 + 9*u^43 + 12*u^41 + 15*u^39 + 
#            16*u^37 + 16*u^35 + 17*u^33 + 15*u^31 + 13*u^29 + 13*u^27 + 9*u^
#            25 + 7*u^23 + 7*u^21 + 4*u^19 + 3*u^17 + 3*u^15 + u^13 + u^11 + u^
#            9, u^45 + u^43 + 2*u^41 + 5*u^39 + 6*u^37 + 8*u^35 + 11*u^33 + 
#            11*u^31 + 11*u^29 + 12*u^27 + 9*u^25 + 7*u^23 + 7*u^21 + 4*u^19 + 
#            3*u^17 + 3*u^15 + u^13 + u^11 + u^9, 
#          u^42 + u^40 + u^38 + 2*u^36 + u^34 + u^32 + u^30, 
#          u^42 + u^40 + 2*u^38 + 4*u^36 + 4*u^34 + 5*u^32 + 5*u^30 + 4*u^28 + 
#            3*u^26 + 2*u^24 + u^22, 2*u^35 + 5*u^33 + 6*u^31 + 8*u^29 + 10*u^
#            27 + 8*u^25 + 7*u^23 + 7*u^21 + 4*u^19 + 3*u^17 + 3*u^15 + u^
#            13 + u^11 + u^9, u^36 + 2*u^34 + 3*u^32 + 4*u^30 + 4*u^28 + 3*u^
#            26 + 2*u^24 + u^22, u^35 + 2*u^33 + u^31 + u^29 + u^27, 
#          2*u^33 + 4*u^31 + 6*u^29 + 9*u^27 + 8*u^25 + 7*u^23 + 7*u^21 + 4*u^
#            19 + 3*u^17 + 3*u^15 + u^13 + u^11 + u^9, 
#          u^33 + 2*u^31 + 4*u^29 + 7*u^27 + 8*u^25 + 9*u^23 + 9*u^21 + 6*u^
#            19 + 4*u^17 + 3*u^15 + u^13 + u^11 + u^9, 
#          u^29 + 4*u^27 + 4*u^25 + 5*u^23 + 6*u^21 + 4*u^19 + 3*u^17 + 3*u^
#            15 + u^13 + u^11 + u^9, u^30 + u^28 + u^26 + u^24, 
#          2*u^27 + 2*u^25 + 2*u^23 + 2*u^21 + u^19, 
#          u^27 + 2*u^25 + 3*u^23 + 5*u^21 + 4*u^19 + 3*u^17 + 3*u^15 + u^
#            13 + u^11 + u^9, u^26 + u^24 + u^22 + u^20, 
#          u^25 + 2*u^23 + 2*u^21 + 2*u^19 + u^17, u^25, 0*u^0, 
#          u^27 + u^25 + u^23 + 3*u^21 + 2*u^19 + 2*u^17 + 3*u^15 + u^13 + u^
#            11 + u^9, u^23 + 3*u^21 + 3*u^19 + 3*u^17 + 3*u^15 + u^13 + u^
#            11 + u^9, u^23 + u^21 + 2*u^19 + u^17, u^23 + u^21 + u^19, u^23, 
#          0*u^0, u^23 + 3*u^21 + 4*u^19 + 4*u^17 + 3*u^15 + u^13 + u^11 + u^9,
#          u^21 + 3*u^19 + 4*u^17 + 3*u^15 + u^13 + u^11 + u^9, 0*u^0, 
#          u^20 + 3*u^18 + 2*u^16, u^21 + u^19 + 3*u^17 + 3*u^15 + u^13 + u^
#            11 + u^9, u^20 + u^18, u^18 + u^16, 
#          u^19 + 3*u^17 + 3*u^15 + u^13 + u^11 + u^9, 2*u^18 + 2*u^16, 
#          2*u^17 + 3*u^15 + u^13 + u^11 + u^9, u^18 + u^16, u^18 + u^16, 
#          2*u^17 + 3*u^15 + u^13 + u^11 + u^9, u^18 + u^16, 0*u^0, 
#          2*u^17 + 3*u^15 + u^13 + u^11 + u^9, 
#          u^17 + 3*u^15 + u^13 + u^11 + u^9, u^17 + 2*u^15 + u^13 + u^11 + u^9
#            , 2*u^15 + u^13 + u^11 + u^9, 2*u^15 + u^13 + u^11 + u^9, u^12, 
#          0*u^0, u^15 + u^13 + u^11 + u^9, u^11 + u^9, u^12, u^13, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^11 + u^9, u^11 + u^9, u^12, 0*u^0, 
#          u^11 + u^9, u^11 + u^9, u^11 + u^9, u^9, 0*u^0, 0*u^0, u^9, u^9, 
#          0*u^0, 0*u^0, u^9, 0*u^0, 0*u^0, u^9, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^7, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^71 + u^69 + 2*u^65 + u^63 + u^61 + 4*u^59 + u^57 + 2*u^55 + 5*u^
#            53 + u^51 + 4*u^49 + 5*u^47 + 5*u^43 + 4*u^41 + u^39 + 5*u^37 + 
#            2*u^35 + u^33 + 4*u^31 + u^29 + u^27 + 2*u^25 + u^21 + u^19, 
#          u^59 + u^55 + 2*u^53 + u^51 + 2*u^49 + 3*u^47 + 4*u^43 + 3*u^41 + u^
#            39 + 4*u^37 + 2*u^35 + u^33 + 4*u^31 + u^29 + u^27 + 2*u^25 + u^
#            21 + u^19, u^49 + u^47 + 2*u^43 + u^41 + u^39 + 3*u^37 + u^35 + u^
#            33 + 3*u^31 + u^29 + u^27 + 2*u^25 + u^21 + u^19, 
#          u^41 + u^37 + u^35 + 2*u^31 + u^29 + 2*u^25 + u^21 + u^19, 0*u^0, 
#          u^44 + u^42 + 2*u^38 + u^34 + 2*u^32 + u^28 + u^26 + u^22, 
#          u^31 + u^25 + u^21 + u^19, u^38 + u^34 + u^32 + u^28 + u^26 + u^22, 
#          0*u^0, u^31 + u^25 + u^21 + u^19, 
#          u^33 + u^31 + u^29 + u^27 + u^25 + u^23 + 2*u^21 + u^19 + u^17, 
#          u^19, 0*u^0, u^29 + u^23 + u^19, u^19, 0*u^0, u^27 + u^21 + u^17, 
#          0*u^0, u^29 + u^25 + u^19, u^23 + u^19 + u^17 + u^13, 0*u^0, u^17, 
#          0*u^0, 0*u^0, u^25 + u^19, u^17, u^17, 0*u^0, u^24 + u^18 + u^16, 
#          0*u^0, 0*u^0, u^16, 0*u^0, u^16, u^13, u^20 + u^16 + u^14, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^13, 0*u^0, u^13, 0*u^0, 0*u^0, 
#          u^17 + u^13 + u^11, u^13, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^11, 0*u^0, 
#          0*u^0, u^15, 0*u^0, 0*u^0, 0*u^0, u^11, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^7, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^7, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^78 + u^76 + 2*u^74 + 3*u^72 + 5*u^70 + 6*u^68 + 10*u^66 + 11*u^64 + 
#            15*u^62 + 18*u^60 + 22*u^58 + 24*u^56 + 30*u^54 + 31*u^52 + 35*u^
#            50 + 37*u^48 + 39*u^46 + 39*u^44 + 42*u^42 + 39*u^40 + 39*u^38 + 
#            37*u^36 + 35*u^34 + 31*u^32 + 30*u^30 + 24*u^28 + 22*u^26 + 18*u^
#            24 + 15*u^22 + 11*u^20 + 10*u^18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^
#            10 + u^8 + u^6, u^60 + 2*u^58 + 4*u^56 + 7*u^54 + 10*u^52 + 14*u^
#            50 + 19*u^48 + 23*u^46 + 27*u^44 + 32*u^42 + 33*u^40 + 35*u^38 + 
#            35*u^36 + 34*u^34 + 31*u^32 + 30*u^30 + 24*u^28 + 22*u^26 + 18*u^
#            24 + 15*u^22 + 11*u^20 + 10*u^18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^
#            10 + u^8 + u^6, u^50 + 2*u^48 + 5*u^46 + 8*u^44 + 14*u^42 + 17*u^
#            40 + 23*u^38 + 25*u^36 + 28*u^34 + 27*u^32 + 28*u^30 + 23*u^28 + 
#            22*u^26 + 18*u^24 + 15*u^22 + 11*u^20 + 10*u^18 + 6*u^16 + 5*u^
#            14 + 3*u^12 + 2*u^10 + u^8 + u^6, 2*u^42 + 3*u^40 + 6*u^38 + 10*u^
#            36 + 14*u^34 + 16*u^32 + 21*u^30 + 19*u^28 + 20*u^26 + 17*u^24 + 
#            15*u^22 + 11*u^20 + 10*u^18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^
#            10 + u^8 + u^6, u^41 + u^39 + 2*u^37 + u^35 + 2*u^33 + u^31 + u^29
#            , u^41 + 2*u^39 + 4*u^37 + 5*u^35 + 8*u^33 + 8*u^31 + 9*u^29 + 
#            8*u^27 + 6*u^25 + 4*u^23 + 2*u^21 + u^19, 
#          u^36 + 3*u^34 + 6*u^32 + 11*u^30 + 13*u^28 + 16*u^26 + 15*u^24 + 
#            14*u^22 + 11*u^20 + 10*u^18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^
#            10 + u^8 + u^6, u^35 + 3*u^33 + 5*u^31 + 7*u^29 + 8*u^27 + 6*u^
#            25 + 4*u^23 + 2*u^21 + u^19, u^34 + u^32 + 2*u^30 + u^28 + u^26, 
#          2*u^32 + 6*u^30 + 10*u^28 + 14*u^26 + 15*u^24 + 14*u^22 + 11*u^20 + 
#            10*u^18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^32 + 4*u^30 + 7*u^28 + 13*u^26 + 16*u^24 + 18*u^22 + 16*u^20 + 
#            14*u^18 + 8*u^16 + 6*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^30 + 2*u^28 + 5*u^26 + 8*u^24 + 10*u^22 + 9*u^20 + 9*u^18 + 6*u^
#            16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, u^27 + u^25 + u^23, 
#          u^28 + 2*u^26 + 5*u^24 + 4*u^22 + 4*u^20 + 2*u^18 + u^16, 
#          u^26 + 4*u^24 + 6*u^22 + 8*u^20 + 9*u^18 + 6*u^16 + 5*u^14 + 3*u^
#            12 + 2*u^10 + u^8 + u^6, u^25 + 2*u^23 + 2*u^21 + u^19, 
#          u^26 + 2*u^24 + 4*u^22 + 5*u^20 + 4*u^18 + 2*u^16 + u^14, 0*u^0, 
#          0*u^0, 2*u^24 + 3*u^22 + 4*u^20 + 7*u^18 + 5*u^16 + 6*u^14 + 4*u^
#            12 + 3*u^10 + u^8 + u^6, 2*u^22 + 4*u^20 + 7*u^18 + 5*u^16 + 5*u^
#            14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          2*u^22 + 3*u^20 + 4*u^18 + 2*u^16 + u^14, u^22 + 2*u^20 + u^18, 
#          0*u^0, 0*u^0, 2*u^22 + 6*u^20 + 9*u^18 + 7*u^16 + 6*u^14 + 3*u^12 + 
#            2*u^10 + u^8 + u^6, u^22 + 3*u^20 + 7*u^18 + 7*u^16 + 6*u^14 + 
#            3*u^12 + 2*u^10 + u^8 + u^6, u^20, 
#          u^21 + 3*u^19 + 5*u^17 + 4*u^15 + u^13, 
#          u^20 + 4*u^18 + 5*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^19 + u^17, u^19 + 2*u^17 + 2*u^15 + u^13, 
#          3*u^18 + 5*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^19 + 3*u^17 + 4*u^15 + u^13, 2*u^18 + 4*u^16 + 6*u^14 + 4*u^12 + 
#            3*u^10 + u^8 + u^6, u^17 + 3*u^15 + u^13, u^17 + 2*u^15, 
#          u^18 + 3*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^17 + 2*u^15, u^16, 3*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          2*u^16 + 6*u^14 + 4*u^12 + 3*u^10 + u^8 + u^6, 
#          u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^6, 
#          u^16 + 3*u^14 + 5*u^12 + 3*u^10 + u^8 + u^6, 
#          u^16 + 3*u^14 + 4*u^12 + 3*u^10 + u^8 + u^6, u^15 + u^13 + u^11, 
#          u^14 + u^12 + u^10, 2*u^14 + 3*u^12 + 3*u^10 + u^8 + u^6, 
#          u^12 + 2*u^10 + u^8 + u^6, u^11, 0*u^0, u^12, u^12 + u^10, u^12, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^12 + 3*u^10 + u^8 + u^6, u^12 + 3*u^10 + u^8 + u^6, 
#          u^11, u^10, 3*u^10 + u^8 + u^6, 2*u^10 + u^8 + u^6, 
#          2*u^10 + u^8 + u^6, u^10 + u^8 + u^6, u^9, 0*u^0, u^10 + u^8 + u^6, 
#          u^6, u^9, 0*u^0, u^8 + u^6, u^9, u^8 + u^6, u^8 + u^6, 0*u^0, u^8, 
#          0*u^0, 0*u^0, u^8 + u^6, u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^6, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^74 + u^72 + 2*u^70 + 3*u^68 + 5*u^66 + 7*u^64 + 9*u^62 + 11*u^60 + 
#            15*u^58 + 17*u^56 + 21*u^54 + 23*u^52 + 26*u^50 + 29*u^48 + 31*u^
#            46 + 32*u^44 + 34*u^42 + 33*u^40 + 34*u^38 + 32*u^36 + 31*u^34 + 
#            29*u^32 + 26*u^30 + 23*u^28 + 21*u^26 + 17*u^24 + 15*u^22 + 11*u^
#            20 + 9*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^58 + 2*u^56 + 4*u^54 + 6*u^52 + 9*u^50 + 13*u^48 + 17*u^46 + 20*u^
#            44 + 24*u^42 + 26*u^40 + 29*u^38 + 29*u^36 + 29*u^34 + 28*u^32 + 
#            26*u^30 + 23*u^28 + 21*u^26 + 17*u^24 + 15*u^22 + 11*u^20 + 9*u^
#            18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^48 + 3*u^46 + 5*u^44 + 9*u^42 + 12*u^40 + 17*u^38 + 19*u^36 + 
#            22*u^34 + 23*u^32 + 23*u^30 + 21*u^28 + 20*u^26 + 17*u^24 + 15*u^
#            22 + 11*u^20 + 9*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^
#            8 + u^6, u^42 + 2*u^40 + 5*u^38 + 7*u^36 + 11*u^34 + 13*u^32 + 
#            16*u^30 + 16*u^28 + 17*u^26 + 15*u^24 + 14*u^22 + 11*u^20 + 9*u^
#            18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 0*u^0, 
#          u^41 + 2*u^39 + 4*u^37 + 5*u^35 + 7*u^33 + 8*u^31 + 8*u^29 + 8*u^
#            27 + 6*u^25 + 5*u^23 + 3*u^21 + 2*u^19 + u^17, 
#          2*u^34 + 4*u^32 + 7*u^30 + 9*u^28 + 12*u^26 + 12*u^24 + 12*u^22 + 
#            10*u^20 + 9*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^35 + 3*u^33 + 5*u^31 + 6*u^29 + 7*u^27 + 6*u^25 + 5*u^23 + 3*u^
#            21 + 2*u^19 + u^17, 0*u^0, u^32 + 4*u^30 + 6*u^28 + 10*u^26 + 
#            11*u^24 + 12*u^22 + 10*u^20 + 9*u^18 + 7*u^16 + 5*u^14 + 3*u^12 + 
#            2*u^10 + u^8 + u^6, u^32 + 4*u^30 + 7*u^28 + 12*u^26 + 15*u^24 + 
#            17*u^22 + 15*u^20 + 14*u^18 + 10*u^16 + 7*u^14 + 4*u^12 + 2*u^
#            10 + u^8 + u^6, u^28 + 3*u^26 + 5*u^24 + 7*u^22 + 7*u^20 + 7*u^
#            18 + 6*u^16 + 5*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 0*u^0, 
#          u^28 + 2*u^26 + 4*u^24 + 4*u^22 + 4*u^20 + 3*u^18 + 2*u^16 + u^14, 
#          u^26 + 3*u^24 + 5*u^22 + 6*u^20 + 7*u^18 + 6*u^16 + 5*u^14 + 3*u^
#            12 + 2*u^10 + u^8 + u^6, 0*u^0, u^26 + 3*u^24 + 5*u^22 + 5*u^20 + 
#            5*u^18 + 3*u^16 + 2*u^14 + u^12, 0*u^0, u^24 + u^22 + u^20, 
#          u^26 + 2*u^24 + 4*u^22 + 5*u^20 + 7*u^18 + 7*u^16 + 7*u^14 + 5*u^
#            12 + 4*u^10 + 2*u^8 + u^6, u^22 + 2*u^20 + 4*u^18 + 4*u^16 + 4*u^
#            14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          2*u^22 + 3*u^20 + 4*u^18 + 3*u^16 + 2*u^14 + u^12, 0*u^0, 0*u^0, 
#          u^22 + u^20, u^22 + 3*u^20 + 6*u^18 + 6*u^16 + 6*u^14 + 4*u^12 + 
#            2*u^10 + u^8 + u^6, u^22 + 3*u^20 + 6*u^18 + 6*u^16 + 6*u^14 + 
#            4*u^12 + 2*u^10 + u^8 + u^6, 0*u^0, 
#          2*u^21 + 4*u^19 + 5*u^17 + 4*u^15 + 3*u^13 + u^11, 
#          2*u^18 + 3*u^16 + 4*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 0*u^0, 
#          u^19 + 2*u^17 + 2*u^15 + 2*u^13 + u^11, 
#          2*u^18 + 3*u^16 + 4*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          u^19 + 3*u^17 + 3*u^15 + 3*u^13 + u^11, 
#          2*u^18 + 4*u^16 + 6*u^14 + 5*u^12 + 4*u^10 + 2*u^8 + u^6, 
#          u^19 + 3*u^17 + 4*u^15 + 3*u^13 + u^11, 0*u^0, 
#          u^16 + 3*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, u^17 + u^15 + u^13, 
#          0*u^0, u^16 + 3*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6, 
#          2*u^16 + 5*u^14 + 5*u^12 + 4*u^10 + 2*u^8 + u^6, 
#          u^14 + u^12 + u^10 + u^8 + u^6, u^16 + 4*u^14 + 5*u^12 + 5*u^10 + 
#            2*u^8 + u^6, u^14 + 2*u^12 + 3*u^10 + 2*u^8 + u^6, 0*u^0, 
#          u^16 + 2*u^14 + 3*u^12 + 2*u^10 + u^8, 
#          u^16 + 3*u^14 + 4*u^12 + 4*u^10 + 2*u^8 + u^6, u^10 + u^8 + u^6, 
#          0*u^0, 0*u^0, 0*u^0, u^12 + 2*u^10 + u^8, u^12 + u^10, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^12, 0*u^0, 0*u^0, 
#          0*u^0, u^12 + 3*u^10 + 2*u^8 + u^6, 2*u^10 + 2*u^8 + u^6, 0*u^0, 
#          u^12 + 2*u^10 + u^8, 2*u^10 + 2*u^8 + u^6, u^8 + u^6, u^8 + u^6, 
#          2*u^10 + 2*u^8 + u^6, u^11 + 2*u^9 + u^7, 0*u^0, u^8 + u^6, 
#          u^10 + u^8 + u^6, u^9 + u^7, 0*u^0, u^8 + u^6, u^9 + u^7, 
#          u^8 + 2*u^6, u^6, 0*u^0, 0*u^0, u^9 + u^7, 0*u^0, 2*u^6, u^6, u^6, 
#          0*u^0, 0*u^0, 0*u^0, u^6, 0*u^0, u^6, 0*u^0, 0*u^0, 0*u^0, u^6, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^77 + u^75 + u^73 + 3*u^71 + 3*u^69 + 4*u^67 + 7*u^65 + 7*u^63 + 9*u^
#            61 + 13*u^59 + 12*u^57 + 15*u^55 + 19*u^53 + 17*u^51 + 21*u^49 + 
#            23*u^47 + 20*u^45 + 24*u^43 + 24*u^41 + 20*u^39 + 23*u^37 + 21*u^
#            35 + 17*u^33 + 19*u^31 + 15*u^29 + 12*u^27 + 13*u^25 + 9*u^23 + 
#            7*u^21 + 7*u^19 + 4*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^59 + u^57 + 3*u^55 + 5*u^53 + 6*u^51 + 10*u^49 + 13*u^47 + 13*u^
#            45 + 18*u^43 + 19*u^41 + 18*u^39 + 21*u^37 + 20*u^35 + 17*u^33 + 
#            19*u^31 + 15*u^29 + 12*u^27 + 13*u^25 + 9*u^23 + 7*u^21 + 7*u^
#            19 + 4*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^49 + 2*u^47 + 3*u^45 + 7*u^43 + 9*u^41 + 11*u^39 + 15*u^37 + 15*u^
#            35 + 15*u^33 + 17*u^31 + 14*u^29 + 12*u^27 + 13*u^25 + 9*u^23 + 
#            7*u^21 + 7*u^19 + 4*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^41 + 2*u^39 + 4*u^37 + 7*u^35 + 8*u^33 + 11*u^31 + 12*u^29 + 10*u^
#            27 + 12*u^25 + 9*u^23 + 7*u^21 + 7*u^19 + 4*u^17 + 3*u^15 + 3*u^
#            13 + u^11 + u^9 + u^7, u^42 + u^40 + u^38 + 2*u^36 + u^34 + u^
#            32 + u^30, 2*u^38 + 2*u^36 + 3*u^34 + 6*u^32 + 4*u^30 + 5*u^28 + 
#            5*u^26 + 2*u^24 + 2*u^22 + u^20, 
#          u^35 + 2*u^33 + 5*u^31 + 7*u^29 + 8*u^27 + 10*u^25 + 8*u^23 + 7*u^
#            21 + 7*u^19 + 4*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^34 + 3*u^32 + 3*u^30 + 5*u^28 + 5*u^26 + 2*u^24 + 2*u^22 + u^20, 
#          u^33 + u^31 + u^29 + u^27, 3*u^31 + 4*u^29 + 7*u^27 + 10*u^25 + 8*u^
#            23 + 7*u^21 + 7*u^19 + 4*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7
#            , u^31 + 2*u^29 + 5*u^27 + 8*u^25 + 9*u^23 + 10*u^21 + 9*u^19 + 
#            6*u^17 + 4*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^29 + 2*u^27 + 4*u^25 + 6*u^23 + 5*u^21 + 6*u^19 + 4*u^17 + 3*u^
#            15 + 3*u^13 + u^11 + u^9 + u^7, u^24, 
#          u^29 + u^27 + 2*u^25 + 4*u^23 + 2*u^21 + 2*u^19 + u^17, 
#          u^25 + 3*u^23 + 3*u^21 + 6*u^19 + 4*u^17 + 3*u^15 + 3*u^13 + u^
#            11 + u^9 + u^7, u^26 + u^24 + 2*u^22 + 2*u^20, 
#          u^23 + 3*u^21 + 2*u^19 + 2*u^17 + u^15, 0*u^0, 0*u^0, 
#          u^23 + u^21 + 2*u^19 + 3*u^17 + 2*u^15 + 3*u^13 + 2*u^11 + u^9 + u^7
#            , u^21 + 4*u^19 + 3*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          2*u^21 + 2*u^19 + 2*u^17 + u^15, u^21 + 2*u^19, 0*u^0, 0*u^0, 
#          3*u^21 + 6*u^19 + 5*u^17 + 4*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^21 + 3*u^19 + 4*u^17 + 4*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^21 + u^19, u^20 + 2*u^18 + 3*u^16 + u^14, 
#          2*u^19 + 3*u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, u^18, 
#          u^20 + u^18 + 2*u^16 + u^14, u^19 + 2*u^17 + 3*u^15 + 3*u^13 + u^
#            11 + u^9 + u^7, u^18 + 3*u^16 + u^14, 
#          u^17 + 2*u^15 + 3*u^13 + 2*u^11 + u^9 + u^7, u^16 + u^14, 
#          u^16 + u^14, u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, u^16, 
#          u^17 + u^15, u^17 + 3*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          2*u^15 + 3*u^13 + 2*u^11 + u^9 + u^7, 
#          u^15 + 2*u^13 + u^11 + u^9 + u^7, u^15 + 2*u^13 + 2*u^11 + u^9 + u^7
#            , u^15 + 2*u^13 + 3*u^11 + u^9 + u^7, u^16 + u^14 + 2*u^12, u^11, 
#          u^13 + u^11 + u^9 + u^7, u^11 + u^9 + u^7, u^12, 0*u^0, u^11, u^11, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^11 + u^9 + u^7, 3*u^11 + u^9 + u^7, 0*u^0, 
#          u^11, u^11 + u^9 + u^7, u^11 + u^9 + u^7, u^9 + u^7, u^7, 0*u^0, 
#          u^10, u^9 + u^7, 0*u^0, u^10, 0*u^0, u^7, 0*u^0, u^7, u^7, u^8, 
#          u^7, 0*u^0, 0*u^0, 2*u^7, u^7, 0*u^0, u^7, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^6, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^76 + 2*u^72 + 2*u^70 + 2*u^68 + 4*u^66 + 6*u^64 + 4*u^62 + 10*u^
#            60 + 9*u^58 + 9*u^56 + 14*u^54 + 14*u^52 + 12*u^50 + 20*u^48 + 
#            15*u^46 + 16*u^44 + 20*u^42 + 16*u^40 + 15*u^38 + 20*u^36 + 12*u^
#            34 + 14*u^32 + 14*u^30 + 9*u^28 + 9*u^26 + 10*u^24 + 4*u^22 + 6*u^
#            20 + 4*u^18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          u^60 + u^58 + u^56 + 4*u^54 + 4*u^52 + 5*u^50 + 10*u^48 + 9*u^46 + 
#            11*u^44 + 15*u^42 + 13*u^40 + 14*u^38 + 18*u^36 + 12*u^34 + 14*u^
#            32 + 14*u^30 + 9*u^28 + 9*u^26 + 10*u^24 + 4*u^22 + 6*u^20 + 4*u^
#            18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          2*u^48 + u^46 + 4*u^44 + 6*u^42 + 7*u^40 + 9*u^38 + 13*u^36 + 9*u^
#            34 + 13*u^32 + 12*u^30 + 9*u^28 + 9*u^26 + 10*u^24 + 4*u^22 + 6*u^
#            20 + 4*u^18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          u^42 + 2*u^40 + 2*u^38 + 6*u^36 + 5*u^34 + 7*u^32 + 9*u^30 + 8*u^
#            28 + 7*u^26 + 10*u^24 + 4*u^22 + 6*u^20 + 4*u^18 + 2*u^16 + 2*u^
#            14 + 2*u^12 + u^8, u^35 + u^31, u^43 + 2*u^39 + 2*u^37 + 2*u^35 + 
#            4*u^33 + 4*u^31 + 2*u^29 + 5*u^27 + 2*u^25 + u^23 + 2*u^21, 
#          u^36 + u^34 + 3*u^32 + 5*u^30 + 5*u^28 + 6*u^26 + 8*u^24 + 4*u^22 + 
#            6*u^20 + 4*u^18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          2*u^33 + 2*u^31 + 2*u^29 + 4*u^27 + 2*u^25 + u^23 + 2*u^21, 
#          u^32 + u^28, u^32 + 3*u^30 + 3*u^28 + 6*u^26 + 7*u^24 + 4*u^22 + 
#            6*u^20 + 4*u^18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          u^32 + 2*u^30 + 4*u^28 + 6*u^26 + 8*u^24 + 7*u^22 + 8*u^20 + 5*u^
#            18 + 4*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          u^28 + u^26 + 5*u^24 + 3*u^22 + 4*u^20 + 4*u^18 + 2*u^16 + 2*u^14 + 
#            2*u^12 + u^8, u^29 + u^25, 2*u^24 + 2*u^22 + 2*u^18, 
#          3*u^24 + 2*u^22 + 3*u^20 + 4*u^18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          u^21, u^26 + u^24 + 3*u^22 + 2*u^20 + u^18 + 2*u^16, 0*u^0, u^24, 
#          2*u^24 + 2*u^22 + 2*u^20 + 4*u^18 + 3*u^16 + 2*u^14 + 4*u^12 + u^8, 
#          u^24 + 3*u^20 + 2*u^18 + 2*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          u^22 + 2*u^20 + 2*u^16, u^20, 0*u^0, 0*u^0, 
#          3*u^20 + 2*u^18 + 3*u^16 + 2*u^14 + 2*u^12 + u^8, 
#          2*u^20 + 2*u^18 + 3*u^16 + 2*u^14 + 2*u^12 + u^8, 0*u^0, 
#          u^21 + 2*u^19 + 2*u^17 + 2*u^15 + u^13, 
#          2*u^18 + u^16 + 2*u^14 + 2*u^12 + u^8, u^19, u^15, 
#          2*u^18 + u^16 + 2*u^14 + 2*u^12 + u^8, 2*u^17 + u^15 + u^13, 
#          2*u^18 + 2*u^16 + 2*u^14 + 4*u^12 + u^8, u^17 + u^15 + u^13, u^17, 
#          u^16 + u^14 + 2*u^12 + u^8, 2*u^17 + u^13, 0*u^0, 
#          u^16 + u^14 + 2*u^12 + u^8, 2*u^16 + u^14 + 4*u^12 + u^8, 
#          u^16 + 2*u^12 + u^8, u^14 + 3*u^12 + u^10 + u^8, 
#          u^16 + 3*u^12 + u^8, 0*u^0, u^16 + 2*u^12, 
#          u^14 + 3*u^12 + u^10 + u^8, u^12 + u^8, 0*u^0, 0*u^0, 0*u^0, u^12, 
#          u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^12 + u^10 + u^8, u^8, u^11, 0*u^0, u^8, u^8, 
#          u^10 + u^8, u^8, 2*u^9, 0*u^0, u^8, u^6, 0*u^0, u^10, u^8, u^9, 
#          u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^9, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^73 + u^71 + 2*u^69 + 4*u^67 + 4*u^65 + 7*u^63 + 10*u^61 + 10*u^59 + 
#            15*u^57 + 18*u^55 + 18*u^53 + 25*u^51 + 26*u^49 + 26*u^47 + 33*u^
#            45 + 31*u^43 + 31*u^41 + 36*u^39 + 31*u^37 + 31*u^35 + 33*u^33 + 
#            26*u^31 + 26*u^29 + 25*u^27 + 18*u^25 + 18*u^23 + 15*u^21 + 10*u^
#            19 + 10*u^17 + 7*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          u^57 + 2*u^55 + 3*u^53 + 7*u^51 + 9*u^49 + 12*u^47 + 18*u^45 + 20*u^
#            43 + 23*u^41 + 29*u^39 + 27*u^37 + 29*u^35 + 31*u^33 + 26*u^31 + 
#            26*u^29 + 25*u^27 + 18*u^25 + 18*u^23 + 15*u^21 + 10*u^19 + 10*u^
#            17 + 7*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          u^47 + 3*u^45 + 5*u^43 + 9*u^41 + 14*u^39 + 16*u^37 + 21*u^35 + 
#            24*u^33 + 22*u^31 + 24*u^29 + 23*u^27 + 18*u^25 + 18*u^23 + 15*u^
#            21 + 10*u^19 + 10*u^17 + 7*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^
#            7 + u^5, u^41 + 2*u^39 + 5*u^37 + 7*u^35 + 12*u^33 + 14*u^31 + 
#            16*u^29 + 19*u^27 + 16*u^25 + 16*u^23 + 15*u^21 + 10*u^19 + 10*u^
#            17 + 7*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          u^38 + u^34 + u^32 + u^28, u^40 + u^38 + 3*u^36 + 6*u^34 + 5*u^32 + 
#            9*u^30 + 8*u^28 + 6*u^26 + 7*u^24 + 4*u^22 + 2*u^20 + 2*u^18, 
#          2*u^33 + 5*u^31 + 8*u^29 + 12*u^27 + 12*u^25 + 14*u^23 + 13*u^21 + 
#            10*u^19 + 10*u^17 + 7*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          u^34 + 2*u^32 + 6*u^30 + 6*u^28 + 6*u^26 + 7*u^24 + 4*u^22 + 2*u^
#            20 + 2*u^18, u^31 + u^29 + u^25, 
#          u^31 + 5*u^29 + 9*u^27 + 10*u^25 + 14*u^23 + 13*u^21 + 10*u^19 + 
#            10*u^17 + 7*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          u^31 + 4*u^29 + 7*u^27 + 12*u^25 + 16*u^23 + 17*u^21 + 16*u^19 + 
#            14*u^17 + 9*u^15 + 6*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          2*u^27 + 4*u^25 + 6*u^23 + 9*u^21 + 8*u^19 + 8*u^17 + 7*u^15 + 4*u^
#            13 + 4*u^11 + 2*u^9 + u^7 + u^5, u^26 + u^22, 
#          u^27 + 3*u^25 + 3*u^23 + 5*u^21 + 4*u^19 + 2*u^17 + 2*u^15, 
#          u^25 + 3*u^23 + 7*u^21 + 7*u^19 + 8*u^17 + 7*u^15 + 4*u^13 + 4*u^
#            11 + 2*u^9 + u^7 + u^5, u^24 + u^22 + u^18, 
#          u^25 + 3*u^23 + 4*u^21 + 6*u^19 + 4*u^17 + 2*u^15 + 2*u^13, 0*u^0, 
#          u^21, u^25 + 2*u^23 + 4*u^21 + 6*u^19 + 6*u^17 + 8*u^15 + 6*u^13 + 
#            5*u^11 + 4*u^9 + u^7 + u^5, 2*u^21 + 3*u^19 + 6*u^17 + 5*u^15 + 
#            4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          2*u^21 + 4*u^19 + 4*u^17 + 2*u^15 + 2*u^13, u^21 + u^17, 0*u^0, 
#          u^21, 2*u^21 + 5*u^19 + 8*u^17 + 7*u^15 + 6*u^13 + 4*u^11 + 2*u^
#            9 + u^7 + u^5, u^21 + 4*u^19 + 7*u^17 + 7*u^15 + 6*u^13 + 4*u^
#            11 + 2*u^9 + u^7 + u^5, 0*u^0, 
#          2*u^20 + 4*u^18 + 5*u^16 + 4*u^14 + 2*u^12, 
#          u^19 + 3*u^17 + 5*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, u^16, 
#          2*u^18 + 2*u^16 + 2*u^14 + 2*u^12, 
#          3*u^17 + 5*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          2*u^18 + 3*u^16 + 4*u^14 + 2*u^12, 
#          2*u^17 + 6*u^15 + 6*u^13 + 5*u^11 + 4*u^9 + u^7 + u^5, 
#          u^18 + 3*u^16 + 3*u^14 + 2*u^12, u^16 + u^14, 
#          u^17 + 3*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, u^16 + 2*u^14, 
#          0*u^0, 3*u^15 + 4*u^13 + 4*u^11 + 2*u^9 + u^7 + u^5, 
#          4*u^15 + 6*u^13 + 5*u^11 + 4*u^9 + u^7 + u^5, 
#          u^15 + 2*u^13 + 2*u^11 + 2*u^9 + u^7 + u^5, 
#          2*u^15 + 4*u^13 + 6*u^11 + 4*u^9 + u^7 + u^5, 
#          u^15 + 3*u^13 + 3*u^11 + 4*u^9 + u^7 + u^5, u^14 + u^10, 
#          u^15 + 2*u^13 + u^11 + 2*u^9, u^15 + 3*u^13 + 5*u^11 + 4*u^9 + u^
#            7 + u^5, u^11 + 2*u^9 + u^7 + u^5, u^10, 0*u^0, 0*u^0, 
#          u^11 + 2*u^9, u^11, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 2*u^11 + 4*u^9 + u^7 + u^5, 
#          u^11 + 4*u^9 + u^7 + u^5, u^10, u^11 + 2*u^9, 
#          u^11 + 4*u^9 + u^7 + u^5, 2*u^9 + u^7 + u^5, 2*u^9 + u^7 + u^5, 
#          3*u^9 + u^7 + u^5, u^10 + 2*u^8, 0*u^0, u^9 + u^7 + u^5, 
#          u^9 + u^7 + u^5, 2*u^8, 0*u^0, u^9 + u^7 + u^5, 2*u^8, 2*u^7 + u^5, 
#          u^7 + u^5, 0*u^0, u^7, u^8, 0*u^0, 2*u^7 + u^5, u^5, u^7, 0*u^0, 
#          0*u^0, 0*u^0, u^5, u^6, 0*u^0, 0*u^0, 0*u^0, u^5, u^5, 0*u^0, 
#          0*u^0, u^5, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^68 + u^64 + 2*u^62 + 2*u^60 + 3*u^58 + 5*u^56 + 4*u^54 + 7*u^52 + 
#            8*u^50 + 8*u^48 + 10*u^46 + 12*u^44 + 10*u^42 + 13*u^40 + 13*u^
#            38 + 12*u^36 + 13*u^34 + 13*u^32 + 10*u^30 + 12*u^28 + 10*u^26 + 
#            8*u^24 + 8*u^22 + 7*u^20 + 4*u^18 + 5*u^16 + 3*u^14 + 2*u^12 + 
#            2*u^10 + u^8 + u^4, u^52 + 2*u^50 + 2*u^48 + 4*u^46 + 6*u^44 + 
#            6*u^42 + 9*u^40 + 10*u^38 + 10*u^36 + 12*u^34 + 12*u^32 + 10*u^
#            30 + 12*u^28 + 10*u^26 + 8*u^24 + 8*u^22 + 7*u^20 + 4*u^18 + 5*u^
#            16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 
#          u^44 + u^42 + 3*u^40 + 4*u^38 + 6*u^36 + 8*u^34 + 9*u^32 + 8*u^30 + 
#            11*u^28 + 9*u^26 + 8*u^24 + 8*u^22 + 7*u^20 + 4*u^18 + 5*u^16 + 
#            3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 
#          u^38 + u^36 + 3*u^34 + 5*u^32 + 4*u^30 + 8*u^28 + 7*u^26 + 7*u^24 + 
#            7*u^22 + 7*u^20 + 4*u^18 + 5*u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^
#            8 + u^4, 0*u^0, 2*u^35 + u^33 + 3*u^31 + 3*u^29 + 3*u^27 + 3*u^
#            25 + 3*u^23 + 2*u^21 + u^19 + u^17, 
#          u^32 + u^30 + 4*u^28 + 4*u^26 + 5*u^24 + 6*u^22 + 6*u^20 + 4*u^18 + 
#            5*u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 
#          u^31 + 2*u^29 + 2*u^27 + 3*u^25 + 3*u^23 + 2*u^21 + u^19 + u^17, 
#          0*u^0, 2*u^28 + 3*u^26 + 4*u^24 + 6*u^22 + 6*u^20 + 4*u^18 + 5*u^
#            16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 
#          2*u^28 + 3*u^26 + 5*u^24 + 7*u^22 + 8*u^20 + 7*u^18 + 7*u^16 + 4*u^
#            14 + 3*u^12 + 2*u^10 + u^8 + u^4, u^26 + u^24 + 3*u^22 + 4*u^20 + 
#            3*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 0*u^0, 
#          u^26 + 2*u^22 + 2*u^20 + 2*u^18 + u^16 + u^14, 
#          2*u^22 + 3*u^20 + 3*u^18 + 4*u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^
#            8 + u^4, 0*u^0, u^24 + u^22 + 2*u^20 + 3*u^18 + 2*u^16 + u^14 + u^
#            12, 0*u^0, 0*u^0, u^22 + 3*u^20 + 2*u^18 + 4*u^16 + 4*u^14 + 3*u^
#            12 + 3*u^10 + 2*u^8 + u^4, u^20 + u^18 + 3*u^16 + 2*u^14 + 2*u^
#            12 + 2*u^10 + u^8 + u^4, u^20 + 2*u^18 + 2*u^16 + u^14 + u^12, 
#          0*u^0, 0*u^0, 0*u^0, u^20 + 2*u^18 + 4*u^16 + 3*u^14 + 3*u^12 + 2*u^
#            10 + u^8 + u^4, u^20 + 2*u^18 + 4*u^16 + 3*u^14 + 3*u^12 + 2*u^
#            10 + u^8 + u^4, 0*u^0, u^19 + 2*u^17 + 2*u^15 + 2*u^13 + u^11, 
#          2*u^16 + 2*u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 0*u^0, 
#          u^17 + u^15 + u^13 + u^11, 2*u^16 + 2*u^14 + 2*u^12 + 2*u^10 + u^
#            8 + u^4, u^17 + u^15 + 2*u^13 + u^11, 
#          2*u^16 + 3*u^14 + 3*u^12 + 3*u^10 + 2*u^8 + u^4, 
#          u^17 + u^15 + 2*u^13 + u^11, 0*u^0, 
#          u^16 + u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, u^13, 0*u^0, 
#          u^14 + 2*u^12 + 2*u^10 + u^8 + u^4, 
#          2*u^14 + 3*u^12 + 3*u^10 + 2*u^8 + u^4, u^12 + u^10 + u^8 + u^4, 
#          2*u^14 + 2*u^12 + 4*u^10 + 2*u^8 + u^4, 
#          u^14 + u^12 + 2*u^10 + 2*u^8 + u^4, 0*u^0, u^14 + u^12 + u^10 + u^8,
#          u^14 + 2*u^12 + 3*u^10 + 2*u^8 + u^4, u^10 + u^8 + u^4, 0*u^0, 
#          0*u^0, 0*u^0, u^10 + u^8, u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^12, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 2*u^10 + 2*u^8 + u^4, 
#          u^10 + 2*u^8 + u^4, 0*u^0, u^10 + u^8, u^10 + 2*u^8 + u^4, 
#          u^8 + u^4, u^8 + u^4, u^10 + 2*u^8 + u^4, u^9 + u^7, 0*u^0, 
#          u^8 + u^4, u^8 + u^4, u^7, 0*u^0, u^8 + u^4, u^7, u^8 + u^6 + u^4, 
#          u^4, 0*u^0, 0*u^0, u^7, 0*u^0, u^6 + u^4, u^4, u^6, 0*u^0, 0*u^0, 
#          0*u^0, u^4, 0*u^0, u^6, 0*u^0, 0*u^0, u^4, u^4, 0*u^0, 0*u^0, u^4, 
#          u^4, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^65 + u^63 + u^61 + 3*u^59 + 2*u^57 + 3*u^55 + 6*u^53 + 4*u^51 + 6*u^
#            49 + 9*u^47 + 6*u^45 + 9*u^43 + 11*u^41 + 7*u^39 + 11*u^37 + 11*u^
#            35 + 7*u^33 + 11*u^31 + 9*u^29 + 6*u^27 + 9*u^25 + 6*u^23 + 4*u^
#            21 + 6*u^19 + 3*u^17 + 2*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^53 + 2*u^49 + 3*u^47 + 3*u^45 + 5*u^43 + 7*u^41 + 5*u^39 + 9*u^
#            37 + 9*u^35 + 7*u^33 + 10*u^31 + 9*u^29 + 6*u^27 + 9*u^25 + 6*u^
#            23 + 4*u^21 + 6*u^19 + 3*u^17 + 2*u^15 + 3*u^13 + u^11 + u^9 + u^7
#            , u^43 + 2*u^41 + 2*u^39 + 5*u^37 + 5*u^35 + 5*u^33 + 8*u^31 + 
#            7*u^29 + 6*u^27 + 8*u^25 + 6*u^23 + 4*u^21 + 6*u^19 + 3*u^17 + 
#            2*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^37 + 2*u^35 + 2*u^33 + 4*u^31 + 5*u^29 + 4*u^27 + 6*u^25 + 6*u^
#            23 + 3*u^21 + 6*u^19 + 3*u^17 + 2*u^15 + 3*u^13 + u^11 + u^9 + u^7
#            , 0*u^0, u^38 + u^36 + u^34 + 3*u^32 + 2*u^30 + 2*u^28 + 4*u^
#            26 + u^24 + 2*u^22 + 2*u^20 + u^16, 
#          u^31 + 2*u^29 + 2*u^27 + 4*u^25 + 4*u^23 + 3*u^21 + 5*u^19 + 3*u^
#            17 + 2*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^32 + u^30 + 2*u^28 + 3*u^26 + u^24 + 2*u^22 + 2*u^20 + u^16, 
#          0*u^0, u^29 + u^27 + 4*u^25 + 3*u^23 + 3*u^21 + 5*u^19 + 3*u^17 + 
#            2*u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          u^29 + 2*u^27 + 4*u^25 + 5*u^23 + 6*u^21 + 6*u^19 + 5*u^17 + 4*u^
#            15 + 3*u^13 + 2*u^11 + u^9 + u^7, u^25 + 2*u^23 + u^21 + 3*u^19 + 
#            3*u^17 + u^15 + 3*u^13 + u^11 + u^9 + u^7, 0*u^0, 
#          2*u^23 + u^21 + u^19 + 2*u^17 + u^13, 
#          u^23 + u^21 + 2*u^19 + 3*u^17 + u^15 + 3*u^13 + u^11 + u^9 + u^7, 
#          0*u^0, u^23 + 3*u^21 + u^19 + 2*u^17 + 2*u^15 + u^11, 0*u^0, 
#          u^23 + u^19, u^23 + u^21 + 2*u^19 + 4*u^17 + 2*u^15 + 4*u^13 + 3*u^
#            11 + u^9 + 2*u^7, u^19 + u^17 + u^15 + 2*u^13 + u^11 + u^9 + u^7, 
#          u^21 + u^19 + u^17 + 2*u^15 + u^11, 0*u^0, 0*u^0, u^19, 
#          u^19 + 2*u^17 + 2*u^15 + 2*u^13 + 2*u^11 + u^9 + u^7, 
#          u^19 + 2*u^17 + 2*u^15 + 2*u^13 + 2*u^11 + u^9 + u^7, 0*u^0, 
#          u^20 + 2*u^18 + 2*u^16 + 2*u^14 + u^12 + u^10, 
#          u^17 + 2*u^13 + u^11 + u^9 + u^7, 0*u^0, u^16 + u^14 + u^10, 
#          u^17 + 2*u^13 + u^11 + u^9 + u^7, 2*u^16 + u^14 + u^12 + u^10, 
#          u^17 + u^15 + 3*u^13 + 3*u^11 + u^9 + 2*u^7, 
#          2*u^16 + 2*u^14 + u^12 + u^10, 0*u^0, u^13 + u^11 + u^9 + u^7, 
#          u^16 + u^12, 0*u^0, u^13 + u^11 + u^9 + u^7, 
#          u^15 + 2*u^13 + 3*u^11 + u^9 + 2*u^7, u^11 + u^7, 
#          2*u^13 + 2*u^11 + 2*u^9 + 2*u^7, 2*u^11 + 2*u^7, 0*u^0, 
#          u^15 + u^13 + 3*u^11 + u^7, 2*u^13 + 2*u^11 + 2*u^9 + 2*u^7, u^7, 
#          0*u^0, 0*u^0, 0*u^0, u^11 + u^7, u^9, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^11, 0*u^0, 0*u^0, 0*u^0, u^11 + u^9 + 2*u^7, 
#          2*u^7, 0*u^0, u^11 + u^7, 2*u^7, u^7, u^7, 2*u^7, 
#          u^10 + 2*u^8 + u^6, 0*u^0, u^7, u^7 + u^5, u^6, 0*u^0, u^7, 
#          u^8 + u^6, u^7 + u^5, 0*u^0, 0*u^0, 0*u^0, u^8 + u^6, u^8, u^5, 
#          0*u^0, u^5, 0*u^0, 0*u^0, u^7, u^5, 0*u^0, u^5, 0*u^0, 0*u^0, 
#          0*u^0, u^5, 0*u^0, u^5, 0*u^0, 0*u^0, u^4, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^64 + u^60 + u^58 + u^56 + u^54 + 2*u^52 + u^50 + 3*u^48 + 2*u^46 + 
#            2*u^44 + 3*u^42 + 3*u^40 + 2*u^38 + 4*u^36 + 2*u^34 + 3*u^32 + 
#            3*u^30 + 2*u^28 + 2*u^26 + 3*u^24 + u^22 + 2*u^20 + u^18 + u^
#            16 + u^14 + u^12 + u^8, u^48 + u^46 + u^44 + 2*u^42 + 2*u^40 + 
#            2*u^38 + 3*u^36 + 2*u^34 + 3*u^32 + 3*u^30 + 2*u^28 + 2*u^26 + 
#            3*u^24 + u^22 + 2*u^20 + u^18 + u^16 + u^14 + u^12 + u^8, 
#          u^40 + u^38 + 2*u^36 + u^34 + 3*u^32 + 2*u^30 + 2*u^28 + 2*u^26 + 
#            3*u^24 + u^22 + 2*u^20 + u^18 + u^16 + u^14 + u^12 + u^8, 
#          u^32 + u^30 + 2*u^28 + u^26 + 3*u^24 + u^22 + 2*u^20 + u^18 + u^
#            16 + u^14 + u^12 + u^8, u^35 + u^31, u^31 + u^27 + u^25 + u^21, 
#          u^28 + u^26 + 2*u^24 + u^22 + 2*u^20 + u^18 + u^16 + u^14 + u^
#            12 + u^8, u^27 + u^25 + u^21, u^28, 
#          u^26 + 2*u^24 + u^22 + 2*u^20 + u^18 + u^16 + u^14 + u^12 + u^8, 
#          u^24 + u^22 + 2*u^20 + u^18 + 2*u^16 + u^14 + u^12 + u^8, 
#          u^24 + u^22 + u^20 + u^18 + u^16 + u^14 + u^12 + u^8, 0*u^0, 
#          u^22 + u^18, u^18 + u^16 + u^14 + u^12 + u^8, u^21 + u^19, u^16, 
#          0*u^0, 0*u^0, u^16 + u^12 + u^8, u^16 + u^14 + u^12 + u^8, u^16, 
#          u^18, 0*u^0, 0*u^0, u^18 + 2*u^16 + u^14 + u^12 + u^8, 
#          u^16 + u^14 + u^12 + u^8, u^18, u^15, u^16 + u^14 + u^12 + u^8, 
#          0*u^0, u^15, u^14 + u^12 + u^8, u^15, u^12 + u^8, 0*u^0, u^13, 
#          u^14 + u^12 + u^8, 0*u^0, u^14, u^14 + u^12 + u^8, u^12 + u^8, 
#          u^12 + u^8, u^8, u^12 + u^10 + u^8, u^13 + u^11, 0*u^0, u^8, u^8, 
#          0*u^0, 0*u^0, u^10, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^8, u^10 + u^8, 
#          0*u^0, 0*u^0, u^8, u^8, u^8, 0*u^0, 0*u^0, u^9, u^8, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^7, u^6, 0*u^0, 0*u^0, u^6, 
#          0*u^0, 0*u^0, u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^5, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^4, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^57 + u^55 + u^53 + 2*u^51 + 2*u^49 + 3*u^47 + 4*u^45 + 4*u^43 + 5*u^
#            41 + 6*u^39 + 6*u^37 + 7*u^35 + 7*u^33 + 7*u^31 + 7*u^29 + 7*u^
#            27 + 7*u^25 + 6*u^23 + 6*u^21 + 5*u^19 + 4*u^17 + 4*u^15 + 3*u^
#            13 + 2*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          u^45 + u^43 + 2*u^41 + 3*u^39 + 4*u^37 + 5*u^35 + 6*u^33 + 6*u^31 + 
#            7*u^29 + 7*u^27 + 7*u^25 + 6*u^23 + 6*u^21 + 5*u^19 + 4*u^17 + 
#            4*u^15 + 3*u^13 + 2*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          u^37 + 2*u^35 + 3*u^33 + 4*u^31 + 5*u^29 + 6*u^27 + 6*u^25 + 6*u^
#            23 + 6*u^21 + 5*u^19 + 4*u^17 + 4*u^15 + 3*u^13 + 2*u^11 + 2*u^
#            9 + u^7 + u^5 + u^3, u^33 + u^31 + 2*u^29 + 4*u^27 + 4*u^25 + 5*u^
#            23 + 5*u^21 + 5*u^19 + 4*u^17 + 4*u^15 + 3*u^13 + 2*u^11 + 2*u^
#            9 + u^7 + u^5 + u^3, 0*u^0, u^30 + u^28 + 2*u^26 + 2*u^24 + 2*u^
#            22 + 2*u^20 + u^18 + u^16, u^27 + 2*u^25 + 3*u^23 + 4*u^21 + 4*u^
#            19 + 4*u^17 + 4*u^15 + 3*u^13 + 2*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          u^26 + 2*u^24 + 2*u^22 + 2*u^20 + u^18 + u^16, 0*u^0, 
#          u^25 + 2*u^23 + 4*u^21 + 4*u^19 + 4*u^17 + 4*u^15 + 3*u^13 + 2*u^
#            11 + 2*u^9 + u^7 + u^5 + u^3, u^25 + 2*u^23 + 4*u^21 + 5*u^19 + 
#            6*u^17 + 6*u^15 + 4*u^13 + 3*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          2*u^21 + 2*u^19 + 3*u^17 + 3*u^15 + 3*u^13 + 2*u^11 + 2*u^9 + u^
#            7 + u^5 + u^3, 0*u^0, u^21 + u^19 + 2*u^17 + u^15 + u^13, 
#          u^21 + u^19 + 3*u^17 + 3*u^15 + 3*u^13 + 2*u^11 + 2*u^9 + u^7 + u^
#            5 + u^3, 0*u^0, u^19 + 2*u^17 + 2*u^15 + u^13 + u^11, 0*u^0, 
#          0*u^0, u^21 + u^19 + 2*u^17 + 3*u^15 + 3*u^13 + 3*u^11 + 3*u^9 + 
#            2*u^7 + u^5 + u^3, u^17 + 2*u^15 + 2*u^13 + 2*u^11 + 2*u^9 + u^
#            7 + u^5 + u^3, u^17 + 2*u^15 + u^13 + u^11, 0*u^0, 0*u^0, 0*u^0, 
#          u^17 + 3*u^15 + 3*u^13 + 3*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          u^17 + 3*u^15 + 3*u^13 + 3*u^11 + 2*u^9 + u^7 + u^5 + u^3, 0*u^0, 
#          u^16 + 2*u^14 + 2*u^12 + u^10, u^15 + 2*u^13 + 2*u^11 + 2*u^9 + u^
#            7 + u^5 + u^3, 0*u^0, u^14 + u^12 + u^10, 
#          u^15 + 2*u^13 + 2*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          u^14 + 2*u^12 + u^10, u^15 + 2*u^13 + 3*u^11 + 3*u^9 + 2*u^7 + u^
#            5 + u^3, u^14 + 2*u^12 + u^10, 0*u^0, 
#          u^13 + 2*u^11 + 2*u^9 + u^7 + u^5 + u^3, u^12, 0*u^0, 
#          u^13 + 2*u^11 + 2*u^9 + u^7 + u^5 + u^3, 
#          u^13 + 3*u^11 + 3*u^9 + 2*u^7 + u^5 + u^3, 
#          u^11 + u^9 + u^7 + u^5 + u^3, u^13 + 2*u^11 + 4*u^9 + 2*u^7 + u^
#            5 + u^3, u^11 + 2*u^9 + 2*u^7 + u^5 + u^3, 0*u^0, 
#          u^11 + u^9 + u^7, u^13 + 2*u^11 + 3*u^9 + 2*u^7 + u^5 + u^3, 
#          u^9 + u^7 + u^5 + u^3, 0*u^0, 0*u^0, 0*u^0, u^9 + u^7, u^9, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 2*u^9 + 2*u^7 + u^5 + u^3, u^9 + 2*u^7 + u^5 + u^3, 0*u^0, 
#          u^9 + u^7, u^9 + 2*u^7 + u^5 + u^3, u^7 + u^5 + u^3, 
#          u^7 + u^5 + u^3, u^9 + 2*u^7 + u^5 + u^3, u^8 + u^6, 0*u^0, 
#          u^7 + u^5 + u^3, u^9 + u^7 + u^5 + u^3, u^6, 0*u^0, u^7 + u^5 + u^3,
#          u^6, u^7 + 2*u^5 + u^3, u^5 + u^3, 0*u^0, 0*u^0, u^6, 0*u^0, 
#          2*u^5 + u^3, u^5 + u^3, u^5, 0*u^0, 0*u^0, 0*u^0, u^5 + u^3, 0*u^0, 
#          u^5, 0*u^0, 0*u^0, u^3, u^5 + u^3, 0*u^0, 0*u^0, u^3, u^3, 0*u^0, 
#          0*u^0, u^3, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^56 + u^54 + u^52 + 2*u^50 + 2*u^48 + 2*u^46 + 4*u^44 + 3*u^42 + 4*u^
#            40 + 5*u^38 + 4*u^36 + 5*u^34 + 6*u^32 + 4*u^30 + 6*u^28 + 5*u^
#            26 + 4*u^24 + 5*u^22 + 4*u^20 + 3*u^18 + 4*u^16 + 2*u^14 + 2*u^
#            12 + 2*u^10 + u^8 + u^6 + u^4, u^44 + u^42 + 2*u^40 + 3*u^38 + 
#            3*u^36 + 4*u^34 + 5*u^32 + 4*u^30 + 6*u^28 + 5*u^26 + 4*u^24 + 
#            5*u^22 + 4*u^20 + 3*u^18 + 4*u^16 + 2*u^14 + 2*u^12 + 2*u^10 + u^
#            8 + u^6 + u^4, u^36 + 2*u^34 + 3*u^32 + 3*u^30 + 5*u^28 + 4*u^
#            26 + 4*u^24 + 5*u^22 + 4*u^20 + 3*u^18 + 4*u^16 + 2*u^14 + 2*u^
#            12 + 2*u^10 + u^8 + u^6 + u^4, u^32 + u^30 + 2*u^28 + 3*u^26 + 
#            3*u^24 + 4*u^22 + 4*u^20 + 3*u^18 + 4*u^16 + 2*u^14 + 2*u^12 + 
#            2*u^10 + u^8 + u^6 + u^4, u^27, u^29 + u^27 + u^25 + 2*u^23 + u^
#            21 + u^19 + u^17, u^26 + 2*u^24 + 3*u^22 + 3*u^20 + 3*u^18 + 4*u^
#            16 + 2*u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          u^25 + 2*u^23 + u^21 + u^19 + u^17, u^24, 
#          u^24 + 3*u^22 + 3*u^20 + 3*u^18 + 4*u^16 + 2*u^14 + 2*u^12 + 2*u^
#            10 + u^8 + u^6 + u^4, u^24 + 2*u^22 + 3*u^20 + 4*u^18 + 5*u^16 + 
#            3*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          2*u^20 + 2*u^18 + 3*u^16 + 2*u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^
#            4, u^21, u^20 + u^18 + u^16 + u^14, 
#          u^20 + u^18 + 3*u^16 + 2*u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          u^17, u^18 + u^16 + u^14 + u^12, 0*u^0, 0*u^0, 
#          u^20 + u^18 + 2*u^16 + 2*u^14 + 2*u^12 + 2*u^10 + 2*u^8 + u^6 + u^4,
#          2*u^16 + u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          u^16 + u^14 + u^12, u^16, 0*u^0, 0*u^0, 
#          2*u^16 + 2*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          u^16 + 2*u^14 + 3*u^12 + 2*u^10 + u^8 + u^6 + u^4, 0*u^0, 
#          u^15 + 2*u^13 + u^11, u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          u^15, u^13 + u^11, u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          2*u^13 + u^11, u^14 + 2*u^12 + 2*u^10 + 2*u^8 + u^6 + u^4, 
#          u^13 + u^11, u^13, 2*u^12 + 2*u^10 + u^8 + u^6 + u^4, u^13, 0*u^0, 
#          2*u^12 + 2*u^10 + u^8 + u^6 + u^4, 
#          2*u^12 + 2*u^10 + 2*u^8 + u^6 + u^4, u^12 + u^10 + u^8 + u^6 + u^4, 
#          u^12 + 2*u^10 + 2*u^8 + u^6 + u^4, u^12 + u^10 + 2*u^8 + u^6 + u^4, 
#          u^9, u^8, u^12 + 2*u^10 + 2*u^8 + u^6 + u^4, u^8 + u^6 + u^4, u^9, 
#          0*u^0, 0*u^0, u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 2*u^8 + u^6 + u^4, 
#          2*u^8 + u^6 + u^4, u^9, u^8, 2*u^8 + u^6 + u^4, u^8 + u^6 + u^4, 
#          u^8 + u^6 + u^4, u^8 + u^6 + u^4, u^7, 0*u^0, u^6 + u^4, 
#          u^8 + u^6 + u^4, u^7, 0*u^0, u^6 + u^4, u^7, u^6 + u^4, u^6 + u^4, 
#          0*u^0, u^6, 0*u^0, 0*u^0, u^6 + u^4, u^4, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, u^4, u^5, 0*u^0, 0*u^0, 0*u^0, u^4, u^4, 0*u^0, 0*u^0, u^4, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^3, 0*u^0, 0*u^0, 0*u^0, 0*u^0 ], 
#      [ u^52 + u^48 + u^46 + 2*u^42 + 2*u^40 + 3*u^36 + u^34 + u^32 + 4*u^
#            30 + u^28 + u^26 + 3*u^24 + 2*u^20 + 2*u^18 + u^14 + u^12 + u^8, 
#          u^42 + u^40 + 2*u^36 + u^34 + u^32 + 3*u^30 + u^28 + u^26 + 3*u^
#            24 + 2*u^20 + 2*u^18 + u^14 + u^12 + u^8, 
#          u^36 + u^32 + 2*u^30 + u^28 + u^26 + 2*u^24 + 2*u^20 + 2*u^18 + u^
#            14 + u^12 + u^8, u^30 + u^28 + 2*u^24 + u^20 + 2*u^18 + u^14 + u^
#            12 + u^8, 0*u^0, u^31 + u^27 + u^25 + u^21 + u^15, 
#          u^24 + u^20 + u^18 + u^14 + u^12 + u^8, u^27 + u^25 + u^21 + u^15, 
#          0*u^0, u^24 + u^20 + u^18 + u^14 + u^12 + u^8, 
#          u^24 + u^22 + 2*u^20 + u^18 + u^16 + u^14 + u^12 + u^10 + u^8, 
#          u^18 + u^12 + u^8, 0*u^0, u^22 + u^18 + u^12, u^18 + u^12 + u^8, 
#          0*u^0, u^20 + u^16 + u^10, 0*u^0, u^18, 
#          u^18 + u^16 + 2*u^12 + u^8 + u^6, u^8, u^16 + u^10, 0*u^0, 0*u^0, 
#          u^18, u^16 + u^10 + u^8, u^16 + u^10 + u^8, 0*u^0, 
#          u^17 + u^15 + u^9, u^8, 0*u^0, u^15 + u^9, u^8, u^15 + u^9, 
#          u^12 + u^8 + u^6, u^15 + u^13 + u^9, 0*u^0, u^8, 0*u^0, 0*u^0, u^8, 
#          u^12 + u^8 + u^6, 0*u^0, u^12 + u^8 + u^6, u^6, 0*u^0, 
#          u^12 + u^10 + u^6, u^12 + u^8 + u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^6, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^10, 
#          0*u^0, 0*u^0, 0*u^0, u^6, u^6, 0*u^0, u^10 + u^6, u^6, 0*u^0, 
#          0*u^0, u^6, u^9 + u^5, 0*u^0, 0*u^0, u^6, u^5, 0*u^0, 0*u^0, u^5, 
#          u^4, 0*u^0, 0*u^0, 0*u^0, u^7 + u^5, 0*u^0, u^4, 0*u^0, u^4, 0*u^0, 
#          0*u^0, u^6, 0*u^0, 0*u^0, u^4, 0*u^0, u^6, 0*u^0, u^4, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^3, 0*u^0, 0*u^0, 0*u^0, u^3, 0*u^0, 0*u^0, 
#          0*u^0 ], 
#      [ u^46 + u^42 + u^40 + u^38 + 2*u^36 + 2*u^34 + u^32 + 3*u^30 + 2*u^
#            28 + 2*u^26 + 3*u^24 + 2*u^22 + 2*u^20 + 3*u^18 + u^16 + 2*u^14 + 
#            2*u^12 + u^10 + u^8 + u^6 + u^2, 
#          u^36 + u^34 + u^32 + 2*u^30 + 2*u^28 + 2*u^26 + 3*u^24 + 2*u^22 + 
#            2*u^20 + 3*u^18 + u^16 + 2*u^14 + 2*u^12 + u^10 + u^8 + u^6 + u^2,
#          u^30 + u^28 + 2*u^26 + 2*u^24 + 2*u^22 + 2*u^20 + 3*u^18 + u^16 + 
#            2*u^14 + 2*u^12 + u^10 + u^8 + u^6 + u^2, 
#          u^26 + u^24 + 2*u^22 + u^20 + 3*u^18 + u^16 + 2*u^14 + 2*u^12 + u^
#            10 + u^8 + u^6 + u^2, 0*u^0, u^25 + u^21 + u^19 + u^15, 
#          u^22 + u^20 + 2*u^18 + u^16 + 2*u^14 + 2*u^12 + u^10 + u^8 + u^
#            6 + u^2, u^21 + u^19 + u^15, 0*u^0, 
#          u^20 + 2*u^18 + u^16 + 2*u^14 + 2*u^12 + u^10 + u^8 + u^6 + u^2, 
#          u^20 + 2*u^18 + 2*u^16 + 3*u^14 + 2*u^12 + 2*u^10 + u^8 + u^6 + u^2,
#          u^18 + u^16 + u^14 + 2*u^12 + u^10 + u^8 + u^6 + u^2, 0*u^0, 
#          u^16 + u^12, u^16 + u^14 + 2*u^12 + u^10 + u^8 + u^6 + u^2, 0*u^0, 
#          u^16 + u^14 + u^10, 0*u^0, 0*u^0, u^16 + u^14 + 2*u^12 + 2*u^10 + u^
#            8 + 2*u^6 + u^2, u^14 + u^12 + u^10 + u^8 + u^6 + u^2, 
#          u^14 + u^10, 0*u^0, 0*u^0, 0*u^0, 
#          u^14 + u^12 + 2*u^10 + u^8 + u^6 + u^2, 
#          u^14 + u^12 + 2*u^10 + u^8 + u^6 + u^2, 0*u^0, u^13 + u^11 + u^9, 
#          u^12 + u^10 + u^8 + u^6 + u^2, 0*u^0, u^9, 
#          u^12 + u^10 + u^8 + u^6 + u^2, u^11 + u^9, 
#          u^12 + 2*u^10 + u^8 + 2*u^6 + u^2, u^11 + u^9, 0*u^0, 
#          u^10 + u^8 + u^6 + u^2, u^11, 0*u^0, u^10 + u^8 + u^6 + u^2, 
#          2*u^10 + u^8 + 2*u^6 + u^2, u^10 + u^6 + u^2, 
#          u^10 + 2*u^8 + 2*u^6 + u^2, u^10 + 2*u^6 + u^2, 0*u^0, u^10 + u^6, 
#          u^10 + 2*u^8 + 2*u^6 + u^2, u^6 + u^2, 0*u^0, 0*u^0, 0*u^0, u^6, 
#          u^8, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, u^8 + 2*u^6 + u^2, 2*u^6 + u^2, 0*u^0, u^6, 
#          2*u^6 + u^2, u^6 + u^2, u^6 + u^2, 2*u^6 + u^2, u^7 + u^5, 0*u^0, 
#          u^6 + u^2, u^6 + u^4 + u^2, u^5, 0*u^0, u^6 + u^2, u^5, 
#          u^6 + u^4 + u^2, u^2, 0*u^0, 0*u^0, u^5, 0*u^0, u^4 + u^2, u^2, 
#          u^4, 0*u^0, 0*u^0, 0*u^0, u^4 + u^2, 0*u^0, u^4, 0*u^0, 0*u^0, u^2, 
#          u^4 + u^2, 0*u^0, u^4, u^2, u^2, u^3, 0*u^0, u^2, 0*u^0, 0*u^0, 
#          u^2, 0*u^0, 0*u^0 ], 
#      [ u^29 + u^23 + u^19 + u^17 + u^13 + u^11 + u^7 + u, 
#          u^23 + u^19 + u^17 + u^13 + u^11 + u^7 + u, 
#          u^19 + u^17 + u^13 + u^11 + u^7 + u, u^17 + u^13 + u^11 + u^7 + u, 
#          0*u^0, u^14, u^13 + u^11 + u^7 + u, u^14, 0*u^0, 
#          u^13 + u^11 + u^7 + u, u^13 + u^11 + u^9 + u^7 + u, u^11 + u^7 + u, 
#          0*u^0, u^11, u^11 + u^7 + u, 0*u^0, u^9, 0*u^0, 0*u^0, 
#          u^11 + u^7 + u^5 + u, u^7 + u, u^9, 0*u^0, 0*u^0, 0*u^0, 
#          u^9 + u^7 + u, u^9 + u^7 + u, 0*u^0, u^8, u^7 + u, 0*u^0, u^8, 
#          u^7 + u, u^8, u^7 + u^5 + u, u^8, 0*u^0, u^7 + u, 0*u^0, 0*u^0, 
#          u^7 + u, u^7 + u^5 + u, u, u^7 + u^5 + u, u^5 + u, 0*u^0, u^5, 
#          u^7 + u^5 + u, u, 0*u^0, 0*u^0, 0*u^0, u^5, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          u^5 + u, u^5 + u, 0*u^0, u^5, u^5 + u, u, u, u^5 + u, u^4, 0*u^0, 
#          u, u^5 + u, u^4, 0*u^0, u, u^4, u^3 + u, u, 0*u^0, 0*u^0, u^4, 
#          0*u^0, u^3 + u, u, u^3, 0*u^0, 0*u^0, 0*u^0, u, 0*u^0, u^3, 0*u^0, 
#          0*u^0, u, u^3 + u, 0*u^0, 0*u^0, u, u, u^2, 0*u^0, u, 0*u^0, u^2, 
#          u, u, 0*u^0 ], 
#      [ u^0, u^0, u^0, u^0, 0*u^0, 0*u^0, u^0, 0*u^0, 0*u^0, u^0, u^0, u^0, 
#          0*u^0, 0*u^0, u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^0, u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^0, u^0, 0*u^0, 0*u^0, u^0, 0*u^0, 0*u^0, 
#          u^0, 0*u^0, u^0, 0*u^0, 0*u^0, u^0, 0*u^0, 0*u^0, u^0, u^0, u^0, 
#          u^0, u^0, 0*u^0, 0*u^0, u^0, u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^0, u^0, 0*u^0, 0*u^0, u^0, u^0, u^0, u^0, 
#          0*u^0, 0*u^0, u^0, u^0, 0*u^0, 0*u^0, u^0, 0*u^0, u^0, u^0, 0*u^0, 
#          0*u^0, 0*u^0, 0*u^0, u^0, u^0, 0*u^0, 0*u^0, 0*u^0, 0*u^0, u^0, 
#          0*u^0, 0*u^0, 0*u^0, 0*u^0, u^0, u^0, 0*u^0, 0*u^0, u^0, u^0, 
#          0*u^0, 0*u^0, u^0, 0*u^0, 0*u^0, u^0, u^0, u^0 ] ] ];
