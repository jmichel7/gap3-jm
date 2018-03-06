#############################################################################
##
#A  hecke.g                CHEVIE library           Meinolf Geck, Jean Michel 
##
#A  $Id: hecke.g,v 1.3 1997/03/21 15:24:50 werner Exp $
##
#Y  Copyright (C) 1992 - 1996  Lehrstuhl D fur Mathematik, RWTH Aachen, and
#Y  Universite Paris VII.
##
##  This file contains  GAP functions for working with Hecke algebras and
##  their character tables.
##

############################################################################
##
#F  PoincarePolynomial( <H> ) . . . Poincare Polynomial of Hecke Algebra <H>
##
## returns  the Poincare Polynomial  of <H> (equal  to the Schur element
## associated  to the 'identity' character).  The point of separate data
## is  that the  Poincare polynomial  is often  more rational than other
## Schur  elements. It coincides with  the generic order for 1-parameter
## algebras.

PoincarePolynomial:=function(H)
  if not IsRec(H) or not IsBound(H.parameter) then
    Error("<H> should be a Hecke algebra");
  fi;
  return Product(ReflectionType(Group(H)),
     t->ReflTypeOps.PoincarePolynomial(t,H.parameter{t.indices}));
end;

############################################################################
##
#F  JInductionTable( <subgroup>, <group> ) . . . J-induction of characters
##  
##  This    function  works   like  'InductionTable'   but   computes  the
##  J-induction, as defined in [Lusztig-book, p.78].
##  
JInductionTable:=function(u,g)local it,i,j,au,ag;
  au:=LowestPowerGenericDegrees(u);
  ag:=LowestPowerGenericDegrees(g);
  if IsHeckeAlgebra(u) then u:=Group(u);fi;
  if IsHeckeAlgebra(g) then g:=Group(g);fi;
  it:=InductionTable(u,g);
  it.what:="JInductionTable";
  it.head:=function(t,option)
    if IsBound(option.TeX) then return SPrint("$J$-Induction from $",
      ReflectionName(t.u,option),"$ to $",ReflectionName(t.g,option),"$");
    else return SPrint("J-Induction from ",
      ReflectionName(t.u,option)," to ",ReflectionName(t.g,option));fi;end;
  for i in [1..Length(au)] do
    for j in [1..Length(ag)] do
      if ag[j]<>au[i] then it.scalar[j][i]:=0;fi;
    od;
  od;
  return it;
end;

if false then
############################################################################
##
#F  HeckeScalarProducts( <ti>, <char1>, <char2> )  . . . . . . . . . . . . .
#F   . . . . . . . . . . . . . . . . . .  scalar products between characters
##
##  'HeckeScalarProducts' specializes the parameters to  1  and  returns the
##  matrix of ordinary scalar products between the  specialized  characters.
##  This only works if the parameters are all equal to 1 or are polynomials
##  with value 1 at 1.
##
HeckeScalarProducts:=function(ti,char1,char2)local specialize;
  specialize:=function(p)
    if IsPolynomial(p) then return Value(p,1);
    else return p;
    fi;
  end;
  if ForAny(ti.parameter,x->specialize(-x[1]/x[2])<>1)
  then Error("parameters must be equal to 1 or have value 1 at 1\n");
  fi;
  char1:=List(char1,x->List(x,specialize));
  char2:=List(char2,x->List(x,specialize));
  return List(char1,i->List(char2,j->ScalarProduct(ti,i,j)));
end;
fi;

#############################################################################
##
#F  HeckeReflectionRepresentation( <H> )  . . . . . . . . . . . . . . . . . . 
#F  the reflection representation of Hecke algebra <H> of a Coxeter group,
#F  returned as a set of matrices for the generators. 
##
HeckeReflectionRepresentation:=function(H)local t,i,j,n,m,a,C,q,W;
  W:=Group(H);
  if Length(Set(H.parameter))>1 or not IsCoxeterGroup(W) then
    Error("Reflexion representation of Cyclotomic Hecke algebras or\n\
           Hecke algebras with unequal parameters not implemented");
  fi;
  q:=-H.parameter[1][1]/H.parameter[1][2];
  C:=[];
  for i in W.generatingReflections do
    C[i]:=[];
    for j in [1..i-1] do
      m:=CoxeterMatrix(W)[i][j];
      if m<>0 then m:=E(m)+E(m)^-1;else m:=2;fi;
      C[i][j]:=(2+m)*q^0;
      if m=-2 then C[j][i]:=0; else C[j][i]:=q;fi;
    od;
    C[i][i]:=q+q^0;
  od;
  t:=[];
  for i in W.generatingReflections do
    a:=[];
    for j in W.generatingReflections do
      a[j]:=W.generatingReflections*0*q;
      a[j][j]:=q;
      a[j][i]:=a[j][i]-C[i][j];
    od;
    Add(t,-H.parameter[1][2]*a);
  od;
  return t;
end;
  
###########################################################################
##
#F HeckeCharValuesGood( <H>, <w> ) . . . character values of T_w^d for good w
##
## 'HeckeCharValuesGood'  computes the values of the irreducible characters
## of the  Iwahori-Hecke algebra <H> on T_w^d, the $d$-th power of the basis
## element corresponding to the good element <w> of the Coxeter group, where 
## $d$ is the order of <w>. The point is that the character table of the Hecke
## algebra is not used, and that all the eigenvalues of T_w^d are monomials
## in H.parameters, so this can be used to find the absolute value of the
## eigenvalues of T_w, a step towards computing the character table of <H>.
#
HeckeCharValuesGood:=function(H,w)local HI,HJ,J,eig,wd,o,W;
  W:=Group(H);
  wd:=GoodCoxeterWord(W,w);if wd=false then 
   Error("Conjugacy Class representatives of ",W," not good");
  fi;
  wd:=Reversed(wd);
  if Length(wd)=0 or wd[Length(wd)][1]<>[1..W.semisimpleRank] then
         Add(wd,[[1..W.semisimpleRank],0]);fi;
  HJ:=HeckeSubAlgebra(H,[]);eig:=[1];
  for J in wd do
    HI:=HeckeSubAlgebra(H,J[1]);
    eig:=InductionTable(Group(HJ),Group(HI)).scalar*eig;
    o:=List(HeckeCentralMonomials(HI),x->x^(J[2]/2));
    eig:=List([1..Length(eig)],i->eig[i]*o[i]);
    HJ:=HI;
  od;
  return eig;
end;
