###########################################################################
##                                                                       ##
#A chevlie.g                                               Meinolf Geck  ##
##                                                                       ##
## This file contains functions for constructing the Lie algebra and the ##
## corresponding  Chevalley groups  associated with a given root system. ##
## It works both in GAP3 and in GAP4.                                    ##
##                                                                       ##
#Y Copyright (C) 2016    Lehrstuhl fuer Algebra, Universitaet Stuttgart  ##
##                                                                       ##
## This program is free software: you can redistribute it and/or modify  ##
## it under the terms of the GNU General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or     ##
## (at your option) any later version.                                   ##
##                                                                       ##
## This program is distributed in the hope that it will be useful,       ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
## GNU General Public License for more details.                          ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program. If not, see <http://www.gnu.org/licenses/>.  ##
##
###########################################################################
Print("                                                                \n");
Print(" ###############################################################\n");
Print(" ## chevlie - Constructing Lie algebras and Chevalley groups  ##\n");
Print(" ##     (by Meinolf Geck, version 1.0, 10 July 2016)          ##\n");
Print(" ##                                                           ##\n");
Print(" ## For first help type \"HelpChevLie();\"");                 ##\n");
Print("                      ##\n");
Print(" ## Updates at www.math.rwth-aachen.de/~CHEVIE/contrib.html   ##\n");
Print(" ## Works in GAP3 and GAP4 -- all comments welcome !          ##\n");
Print(" ###############################################################\n");
Print("                                                                \n");

##########################################################################
##
#F HelpChevLie() . . . . . . . . . . . . . . . . . . . . . . .  first help 
##
## 'HelpChevLie'  returns the list of the most important functions in this
## package. 
##
HelpChevLie:=function()
  Print(" WeylRecord\n"); 
  Print(" LieAdjointRepresentation\n");
  Print(" LieScRepresentation\n");
  Print(" CanonicalChevalleyBasis\n");
  Print(" OrbitWeight\n");
  Print(" MinusculeWeight\n");
  Print(" ChevalleyGroupGenerators\n");
  Print(" ChevalleyGroupAdj\n");
  Print(" ChevalleyGroupSc\n");
  Print(" ChevalleyRootElement\n");
  Print(" MonomialSubgroup\n");
  Print(" ChevalleyCommRels\n");
end;

GAP4:=VERSION[1]='4'; 
if GAP4 then
  Read("chevlie_gap4.g");
else
WeylRecord:=function(arg)local W,d,c,rng,i,j;
  W:=ApplyFunc(CoxeterGroup,arg);
  d:=BipartiteDecomposition(W);
  W.epsilon:=List(W.matgens,x->1);
  W.epsilon{d[2]}:=d[2]*0-1;
  W.longestword:=LongestCoxeterWord(W);
  W.longestperm:=LongestCoxeterElement(W);
  c:=W.cartan;
  rng:=DiagonalMat(W.rootLengths);
  for i in [1..Length(c)] do 
    for j in [1..Length(c)] do 
      if i<>j then
        if c[i][j]<>0 then
          if W.epsilon[i]=1 then
            rng[i][j]:=c[i][j]*W.rootLengths[i];
          else
            rng[j][i]:=c[j][i]*W.rootLengths[j];
          fi;
        fi;
      fi;
    od;
  od;
  W.ringel:=rng;
  return W;
end;

MinusculeWeights:=W->W.roots{Flat(WeightInfo(W).minusculeWeights)};

myisdiag:=IsDiagonalMat;
DirSumMat:=DiagonalMat;

ScalRootCoRoot:=function(cmat,alpha,i)
  return Sum(List([1..Length(alpha)],j->alpha[j]*cmat[i][j]));
end;
fi;
##########################################################################
##
#F LieAdjointRepresentation( <W> ) . . . . . . generators for adjoint rep.
##
## 'LieAdjointRepresentation'   returns matrices  for the generators  e_i, 
## f_i,h_i under the adjoint representation.
##
## See also 'LieScRepresentation' and 'LieMinusculeRepresentation'.

# check relations for Chevalley generators
ChevalleyCheck:=function(W,reps)
  local i,j,l,h,e,f,res;
  res:=true;
  l:=Length(reps[1]);
  e:=reps[1];
  f:=reps[2];
  h:=reps[3];
  for i in [1..l] do
    for j in [i+1..l] do
      if h[i]*h[j]<>h[j]*h[i] then
        Print("#W mist rels1\n");
        res:=false;
      fi;
    od;
  od;
  for i in [1..l] do
    for j in [1..l] do 
      if i=j then
        if e[i]*f[i]-f[i]*e[i]<>h[i] then 
          Print("#W mist rels2\n");
          res:=false;
        fi;
      else
        if e[i]*f[j]<>f[j]*e[i] then 
          Print("#W mist rels3\n");
          res:=false;
        fi;
      fi;
    od;
  od;
  for i in [1..l] do
    for j in [1..l] do 
      if h[j]*e[i]-e[i]*h[j]<>W.cartan[j][i]*e[i] then 
        Print("#W mist rels4\n");
        res:=false;
      fi;
      if h[j]*f[i]-f[i]*h[j]<>-W.cartan[j][i]*f[i] then 
        Print("#W mist rels5\n");
        res:=false;
      fi;
    od;
  od;
  if res=false then
    Print("# NOT OK !!!\n");
  fi;
  return res;
end;

LieAdjointRepresentation:=function(W)
  local R,Rp,Ei,Fi,Hi,l,s,a,b,n,i,j,me,mf,mh,omega;
  l:=Length(W.cartan);
  Rp:=W.roots{[1..W.N]};
  R:=Concatenation(Reversed(Rp),-Rp);
  Ei:=[]; Fi:=[]; Hi:=[];
  for s in [1..l] do
    me:=0*IdentityMat(2*W.N+l,Rationals);
    for a in [1..2*W.N] do
      if R[a]=-Rp[s] then
        me[W.N+s][a+l]:=1;
      elif R[a]+Rp[s] in R then 
        b:=Position(R,R[a]+Rp[s]);
        n:=1;
        while R[a]-n*Rp[s] in R do
          n:=n+1;
        od;
        if a<=W.N then 
          me[b][a]:=n;
        else
          me[b+l][a+l]:=n;
        fi;
      fi;
    od;
    for j in [1..l] do
      if W.cartan[j][s]>=0 then 
        me[W.N-s+1][W.N+j]:=W.cartan[j][s];
      else
        me[W.N-s+1][W.N+j]:=-W.cartan[j][s];
      fi;
    od;
    Add(Ei,me);
    mf:=0*IdentityMat(2*W.N+l,Rationals);
    for a in [1..2*W.N] do
      if R[a]=Rp[s] then
        mf[W.N+s][a]:=1;
      elif R[a]-Rp[s] in R then 
        b:=Position(R,R[a]-Rp[s]);
        n:=1;
        while R[a]+n*Rp[s] in R do
          n:=n+1;
        od;
        if a<=W.N then 
          mf[b][a]:=n;
        else
          mf[b+l][a+l]:=n;
        fi;
      fi;
    od;
    for j in [1..l] do
      if W.cartan[j][s]>=0 then 
        mf[W.N+l+s][W.N+j]:=W.cartan[j][s];
      else
        mf[W.N+l+s][W.N+j]:=-W.cartan[j][s];
      fi;
    od;
    Add(Fi,mf);
    mh:=0*IdentityMat(2*W.N+l,Rationals);
    for a in [1..W.N] do
      mh[a][a]:=ScalRootCoRoot(W.cartan,R[a],s);
      mh[W.N+l+a][W.N+l+a]:=ScalRootCoRoot(W.cartan,R[W.N+a],s);
    od;
    Add(Hi,mh);
  od;
  Print("#I dim = ",2*W.N+l,", Chevalley relations \c");
  omega:=0*IdentityMat(2*W.N+l,Rationals);
  for i in [1..Length(omega)] do 
    omega[i][Length(omega)+1-i]:=1;
  od;
  for i in [1..l] do 
    for j in [1..l] do 
      if i=j then 
        omega[W.N+i][W.N+j]:=1;
      else
        omega[W.N+i][W.N+j]:=0;
      fi;
    od;
  od;
  for i in [1..l] do 
    if Ei[i]*omega<>omega*Fi[i] or Hi[i]*omega<>-omega*Hi[i] then
      Error("mist!");
    else
      Print(".\c"); 
    fi;
  od;
  Print(" ",ChevalleyCheck(W,[Ei,Fi,Hi]),"\n");
  return [Ei,Fi,Hi,omega];
end;

##########################################################################
##
#F OrbitWeights( <W>, <l> ) . . . . . . . . . . . . . .  orbit of a weight
##
## 'OrbitWeight' returns the orbit of a weight <l> under the action of the
## Weyl group <W>.
##
## See also 'MinusculeWeight'.
##
OrbitMats:=function(mats,l)
  local orb,o,o1,g;
  orb:=[l];
  for o in orb do 
    for g in mats do
      o1:=o*g;
      if not o1 in orb then
        Add(orb,o1);
      fi;
    od;
  od;
  return orb;
end;

OrbitWeight:=function(W,l)
  local A,A1;
  A:=TransposedMat(W.cartan);
  A1:=A^(-1);
  return OrbitMats(List(W.matgens,g->A1*g*A),l);
end;

OrbitWeight1:=function(W,l)
  local A1,o,o1,w,w1;
  A1:=TransposedMat(W.cartan^(-1));
  o:=List(OrbitWeight(W,l),w->w-l);
  o1:=[];
  for w in o do
    Add(o1,w*A1);
  od;
  return o1;
end;


ScalWeightCoRoot:=function(cmat,alpha,i)
  return Sum(List([1..Length(alpha)],j->alpha[j]*cmat[i][j]));
end;

##########################################################################
##
#F LieMinusculeRepresentation( <W>, <lambda> ) . . . . . .  generators for 
## . . . . . . . . . . . . . .  representation with given miniscule weight
##
## 'LieMinusculeRepresentation'  returns matrices for the generators  e_i, 
## f_i,h_i under the representation with a given minuscule weight.
##
## See also 'MinusculeWeight' and 'OrbitWeight'.
##
LieMinusculeRepresentation:=function(W,lambda)
  local R,Rp,Ei,Fi,Hi,l,s,a,b,me,mf,mh;
  l:=Length(W.cartan);
  Rp:=TransposedMat(W.cartan);
  R:=OrbitWeight(W,lambda);
  Ei:=[]; Fi:=[]; Hi:=[];
  for s in [1..l] do
    me:=0*IdentityMat(Length(R),Rationals);
    mf:=0*IdentityMat(Length(R),Rationals);
    mh:=0*IdentityMat(Length(R),Rationals);
    for a in [1..Length(R)] do
      if R[a]+Rp[s] in R then 
        me[Position(R,R[a]+Rp[s])][a]:=1;
      fi;
      if R[a]-Rp[s] in R then 
        mf[Position(R,R[a]-Rp[s])][a]:=1;
      fi;
      #mh[a][a]:=ScalWeightCoRoot(W.cartan,R[a],s);
      mh[a][a]:=R[a][s];
    od;
    Add(Ei,me);
    Add(Fi,mf);
    Add(Hi,mh);
  od;
  Print("#I dim = ",Length(R),", Chevalley relations ",
                          ChevalleyCheck(W,[Ei,Fi,Hi]),"\n");
  return [Ei,Fi,Hi];
end;

##########################################################################
##
#F LieScRepresentation( <W> ) . . . . generators for representation giving
## . . . . . . . . . . . . . . . . . . .  simply-connected Checalley group
##
## 'LieScRepresentation' returns matrices for the generators e_i, f_i, h_i 
## such that corresponding Chevalley group is of simply-connected type.
##
## See also 'LieAdjointRepresentation'.
##
LieScRepresentation:=function(W)
  local l,vec1,rep1,rep2;
  l:=Length(W.cartan);
  if l>=1 and W.cartan=CartanMat("A",l) then
    vec1:=0*[1..l];
    vec1[1]:=1;
    return LieMinusculeRepresentation(W,vec1);
  elif l>=2 and W.cartan=CartanMat("B",l) then
    vec1:=0*[1..l];
    vec1[1]:=1;
    return LieMinusculeRepresentation(W,vec1);
  elif l>=2 and W.cartan=CartanMat("C",l) then
    vec1:=0*[1..l];
    vec1[l]:=1;
    return LieMinusculeRepresentation(W,vec1);
  elif l>=3 and W.cartan=CartanMat("D",l) then
    vec1:=0*[1..l];
    vec1[1]:=1;
    rep1:=LieMinusculeRepresentation(W,vec1);
    vec1:=0*[1..l];
    vec1[2]:=1;
    rep2:=LieMinusculeRepresentation(W,vec1);
    rep1:=[List([1..l],x->DirSumMat(rep1[1][x],rep2[1][x])),
             List([1..l],x->DirSumMat(rep1[2][x],rep2[2][x])),
               List([1..l],x->DirSumMat(rep1[3][x],rep2[3][x]))];
    Print("#I dim = ",Length(rep1[1][1]),
             ", Chevalley relations ",ChevalleyCheck(W,rep1),"\n");
    return rep1;
  elif l=6 and W.cartan=CartanMat("E",6) then
    vec1:=0*[1..l];
    vec1[1]:=1;
    return LieMinusculeRepresentation(W,vec1);
  elif l=7 and W.cartan=CartanMat("E",7) then
    vec1:=0*[1..l];
    vec1[l]:=1;
    return LieMinusculeRepresentation(W,vec1);
  else
    return LieAdjointRepresentation(W);
  fi;
end;

##########################################################################
##
#F G2Representation( <W> ) 
##
## 'G2Representation'  returns  the 6-dimensional matrix representation of 
## the Lie algebra of type G_2.
## 
G2Representation:=function(W)
  local W1,Ei,Fi,Hi,i,r;
  if W.cartan<>CartanMat("G",2) then
    Error("wrong Cartan matrix");
  fi;
  Ei:=[[[0,0,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],
       [0,0,0,0,0,1,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]],[[0,1,0,0,0,0,0],
       [0,0,0,0,0,0,0],[0,0,0,2,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,0,0],
       [0,0,0,0,0,0,1],[0,0,0,0,0,0,0]]];
  Fi:=[[[0,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,0,0,0,0,0],
       [0,0,0,0,0,0,0],[0,0,0,0,1,0,0],[0,0,0,0,0,0,0]],[[0,0,0,0,0,0,0],
       [1,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,2,0,0,0],
       [0,0,0,0,0,0,0],[0,0,0,0,0,1,0]]];
  Hi:=List([1..2],i->Ei[i]*Fi[i]-Fi[i]*Ei[i]);
  Print("#I dim = ",Length(Ei[1][1]),", Chevalley relations ",
                              ChevalleyCheck(W,[Ei,Fi,Hi]),"\n");
  return [Ei,Fi,Hi];
end;

##########################################################################
##
#F F4Representation( <W> ) 
##
## 'F4Representation' returns the 27-dimensional matrix representation of 
## the Lie algebra of type F_4.
## 
F4Representation:=function(W)
  local Ei,Fi,Hi,i,r;
  if W.cartan<>CartanMat("F",4) then
    Error("wrong Cartan matrix");
  fi;
  r:=LieScRepresentation(WeylRecord("E",6)); 
  Ei:=[r[1][2],r[1][4],r[1][3]+r[1][5],r[1][1]+r[1][6]];
  Fi:=[r[2][2],r[2][4],r[2][3]+r[2][5],r[2][1]+r[2][6]];
  Hi:=List([1..4],i->Ei[i]*Fi[i]-Fi[i]*Ei[i]);
  Print("#I dim = ",Length(Ei[1][1]),", Chevalley relations ",
                              ChevalleyCheck(W,[Ei,Fi,Hi]),"\n");
  return [Ei,Fi,Hi];
end;

TestChevalleyGenerators:=function(W,P,Q)
  local R,Rp,Ei,Fi,Hi,l,s,a,b,n,p,j,j1,me,mf,mh;
  #if P*TransposedMat(Q)<>W.cartan then
  #  Print("#W not a valid root datum\n");
  #  return false;
  #fi;
  # P=coroots, Q=roots; cartan of size l times l; P,Q of size l times n
  l:=Length(W.cartan);
  n:=Length(P[1]);
  Rp:=W.roots{[1..W.N]};
  R:=Concatenation(Reversed(Rp),-Rp);
  Ei:=[]; Fi:=[]; Hi:=[];
  for s in [1..l] do
    me:=0*IdentityMat(2*W.N+n,Rationals);
    for a in [1..2*W.N] do
      if R[a]=-Rp[s] then
        for j in [1..n] do 
          me[a+n][W.N+j]:=Q[s][j];
        od;
      elif R[a]+Rp[s] in R then 
        b:=Position(R,R[a]+Rp[s]);
        p:=1;
        while R[a]-p*Rp[s] in R do
          p:=p+1;
        od;
        if a<=W.N then 
          me[a][b]:=p;
        else
          me[a+n][b+n]:=p;
        fi;
      fi;
    od;
    for j in [1..n] do
      if P[s][j]>=0 then 
        me[W.N+j][W.N-s+1]:=P[s][j];
      else
        me[W.N+j][W.N-s+1]:=-P[s][j];
      fi;
    od;
    Add(Ei,TransposedMat(me));
    mf:=0*IdentityMat(2*W.N+n,Rationals);
    for a in [1..2*W.N] do
      if R[a]=Rp[s] then
        for j in [1..n] do 
          mf[a][W.N+j]:=Q[s][j];
        od;
      elif R[a]-Rp[s] in R then 
        b:=Position(R,R[a]-Rp[s]);
        p:=1;
        while R[a]+p*Rp[s] in R do
          p:=p+1;
        od;
        if a<=W.N then 
          mf[a][b]:=p;
        else
          mf[a+n][b+n]:=p;
        fi;
      fi;
    od;
    for j in [1..n] do
      if P[s][j]>=0 then 
        mf[W.N+j][W.N+n+s]:=P[s][j];
      else
        mf[W.N+j][W.N+n+s]:=-P[s][j];
      fi;
    od;
    Add(Fi,TransposedMat(mf));
    mh:=0*IdentityMat(2*W.N+n,Rationals);
    for a in [1..W.N] do
      mh[a][a]:=ScalRootCoRoot(W.cartan,R[a],s);
      mh[W.N+n+a][W.N+n+a]:=ScalRootCoRoot(W.cartan,R[W.N+a],s);
    od;
    Add(Hi,TransposedMat(mh));
  od;
  Print("#I Relations ",ChevalleyCheck(W,[Ei,Fi,Hi]),"\n");
  return [Ei,Fi,Hi];
end;

##########################################################################
##
#F ChevalleyGroupGenerators( <W>, <reps>, <t> ) . . . . . . generators for 
## . . . . . . . . . . . . . . . Chevalley group in a given representation
##
## 'ChevalleyGroupGenerators'  returns  the  generators  x_i(t) and y_i(t)
## of a Chevalley group with root system given by <W>.
##
## See also 'ChevalleyGroupAdj' and 'ChevalleyGroupSc':
##
ChevalleyGroupGenerators:=function(W,repr,t)
  local weiter,i,j,a,a1,g,p,pow,gens;
  gens:=[];
  for j in repr{[1,2]} do 
    for a in j do
      weiter:=true;
      a1:=a;
      i:=1;
      pow:=[ShallowCopy(a)];
      while weiter do
        a1:=a1*a;
        if a1=0*a1 then
          weiter:=false;
        else
          i:=i+1;
          a1:=i^(-1)*a1;
          Add(pow,a1);
          if ForAny(Set(Flat(a1)),x->not IsInt(x)) then
            Print(Set(Flat(a1)),"\n");
            return false;
          fi;
        fi;
      od;
      g:=IdentityMat(Length(a1),t^0);
      for p in [1..Length(pow)] do 
        g:=g+t^p*pow[p];
      od;
      Add(gens,g);
    od;
  od;
  return gens;
  #return Group(gens,IdentityMat(Length(gens[1]),K1[1]^0));
end;

##########################################################################
##
#F ChevalleyGroupRepr( <W>, <repr>, <K> ) . . Chevalley group over a field
##
## 'ChevalleyGroupRepr' returns the Chevalley group of type <W> in a given
## representation <repr> over the field <K>.
##
ChevalleyGroupRepr:=function(W,repr,K)
  local x,y,su,f,gens;
  gens:=[];
  f:=[];
  for x in Elements(K) do
    if x<>0*x and not x in f then 
      Append(gens,ChevalleyGroupGenerators(W,repr,x));
      su:=[x];
      for y in su do
        if not (x+y in su) then
          Add(su,x+y);
          AddSet(f,x+y);
        fi;
      od;
    fi;
  od;
  return Group(gens,gens[1]^0);
end;

ChevalleyGroupAdj:=function(W,K)
  return ChevalleyGroupRepr(W,LieAdjointRepresentation(W),K);
end;

ChevalleyGroupSc:=function(W,K)
  return ChevalleyGroupRepr(W,LieScRepresentation(W),K);
end;

ChevalleyGroupMin:=function(W,l,K)
  return ChevalleyGroupRepr(W,LieMinusculeRepresentation(W,l),K);
end;

ChevalleyGroupG2dim6:=function(K)
  local W;
  W:=WeylRecord("G",2);
  return ChevalleyGroupRepr(W,G2Representation(W),K);
end;

ChevalleyGroupF4dim27:=function(K)
  local W;
  W:=WeylRecord("F",4);
  return ChevalleyGroupRepr(W,F4Representation(W),K);
end;

SplitTorus:=function(W,repr,t)
  local gens1,gens2,s,l,gen,gent;
  l:=Length(W.cartan);
  gens1:=ChevalleyGroupGenerators(W,repr,t);
  gens2:=ChevalleyGroupGenerators(W,repr,-t^(-1));
  gent:=List([1..l],s->gens1[s]*gens2[l+s]*gens1[s]);
  gens1:=ChevalleyGroupGenerators(W,repr,-t^0);
  gens2:=ChevalleyGroupGenerators(W,repr,t^0);
  gen:=List([1..l],s->gent[s]*gens1[s]*gens2[l+s]*gens1[s]);
  for s in [1..l] do
    if myisdiag(gen[s])=false then
      Error("mist diag");
    fi;
  od;
  return gen;
end;

##########################################################################
##
#F MonomialSubgroup( <W>, <reps> ) . . . . . . . . . . . . .  lifts of s_i
##
## 'MononimalSubgroup'  returns the lifts of the simple reflections of <W>
## to elements in the given representations <reps> of the Lie algebra.
## (These satisfy the braid relations, they may have order 2 or 4.)
##
MonomialSubgroup:=function(W,repr)
  local omat,br,y,d,i,j,gens1,gens2,s,l,gen;
  l:=Length(W.cartan);
  gens1:=ChevalleyGroupGenerators(W,repr,1);
  gens2:=ChevalleyGroupGenerators(W,repr,-1);
  gen:=List([1..l],s->gens1[s]*gens2[l+s]*gens1[s]);
  omat:=IdentityMat(l,Rationals);
  for i in [1..l] do 
    for j in [1..i] do
      y:=gen[i]*gen[j];
      d:=1;
      while y<>y^0 do
        y:=y*gen[i]*gen[j];
        d:=d+1;
      od;
      omat[i][j]:=d;
      if i<>j then
        omat[j][i]:=d;
      fi;
    od;
  od;
  # check braid relations
  Print("#I Braid relations \c");
  br:=true;
  for i in [1..l] do
    Print(".\c");
    for j in [1..i-1] do
      if W.cartan[i][j]*W.cartan[j][i]=0 and 
                          gen[i]*gen[j]<>gen[j]*gen[i] then 
        Print("Mist: ",[i,j],"\n");
        br:=false;
      elif W.cartan[i][j]*W.cartan[j][i]=1 and 
        gen[i]*gen[j]*gen[i]<>gen[j]*gen[i]*gen[j] then 
          Print("Mist: ",[i,j],"\n");
          br:=false;
      elif W.cartan[i][j]*W.cartan[j][i]=2 and 
        (gen[i]*gen[j])^2<>(gen[j]*gen[i])^2 then 
          Print("Mist: ",[i,j],"\n");
          br:=false;
      elif W.cartan[i][j]*W.cartan[j][i]=3 and 
        (gen[i]*gen[j])^3<>(gen[j]*gen[i])^3 then 
          Print("Mist: ",[i,j],"\n");
          br:=false;
      fi;
    od;
  od;
  Print(" ",br,"\n");
  PrintArray(omat);
  return gen;
end;

CheckAdjMonomial:=function(W)
  local i,j,l,pj,adj,r,pr,z,mi,ms,mat;
  adj:=LieAdjointRepresentation(W);
  ms:=MonomialSubgroup(W,adj);
  for i in [1..Length(ms)] do 
    if ms[i]*adj[4]<>adj[4]*ms[i]^(-1) then
      Error("Mist 0a!");
    fi;
  od;
  Print("I +\c");
  for i in [1..Length(ms)] do 
    if ms[i]*adj[1][i]*ms[i]^(-1)<>-adj[2][i] or 
                   ms[i]*adj[2][i]*ms[i]^(-1)<>-adj[1][i] then
      Error("Mist 0b!");
    fi;
  od;
  Print("+\c");
  for i in [1..Length(ms)] do 
    for j in [1..Length(ms)] do 
      if ms[i]*adj[3][j]*ms[i]^(-1)<>adj[3][j]-W.cartan[j][i]*adj[3][i] then
        Error("Mist 0c!");
      fi;
    od;
  od; 
  Print("+\c");
  for i in [1..Length(ms)] do
    for j in [1..Length(ms)] do
      if i=j then
        if ms[i][W.N+j][W.N+j]<>-1 then 
          Error("Mist 1!");
        fi;
      else
        if ms[i][W.N+i][W.N+j]<>W.cartan[j][i] then 
          Error("Mist 2!");
        fi;
      fi;
    od;
    l:=Concatenation([1..W.N],[(W.N+Length(ms)+1)..Length(ms[1])]);
    mat:=ms[i]{l}{l};
    #PrintArray(mat);
    for r in [1..Length(W.rootsadjoint)] do 
      pr:=Position(W.rootsadjoint,
           W.roots[Position(W.roots,W.rootsadjoint[r])^W.permgens[i]]);
      z:=0;
      for j in [1..Length(mat)] do
        if mat[j][r]<>0 then
          z:=z+1;
          pj:=j;
        fi;
      od;
      if z<>1 then 
        Error("Mist 3a!");
      fi;
      if pj<>pr then 
    #    Print("\n #W ",[mat,pj,pr],"\n");
        Error("Mist 3b!");
      fi;
      mi:=1;
      while W.rootsadjoint[r]-mi*W.roots[i] in W.roots do
        mi:=mi+1;
      od;
      if mat[pr][r]<>-(-1)^mi then
        Error("Mist 4!");
      fi;
    od;
  od;
  Print(" Formulae for n_i's OK !\n");
  return true;
end;

RingelProd:=function(W,a,b)
  local i,j,r;
  r:=0;
  for i in [1..Length(a)] do
    for j in [1..Length(b)] do
      r:=r+a[i]*W.ringel[i][j]*b[j];
    od;
  od;
  return r;
end;

RingelRule:=function(W,z,x)
  local zx,xz;
  xz:=RingelProd(W,x,z);
  zx:=RingelProd(W,z,x);
  if zx<0 and xz>=0 then
    return [1,xz];
  elif zx>=0 and xz<0 then
    return [-1,zx];
  else
    Error("mist (ringelrule)");
  fi;
end;

#########################################################################
##
#F CanonicalChevalleyBasis( <W> ) . . . . . . . canonical Chevalley basis 
#F CanonicalChevalleyBasisRep( <W>, <repr> )  . . . . . . . . . . . . . . 
##
## 'CanonicalChevalleyBasis'  returns the canonical Chevalley basis  of a
## Lie algebra of type <W> (explicitly as list of matrices in the adjoint 
## representation  as returned by  'LieAdjointRepresentation').   It uses 
## the pm 1 function 'epsilon' defined in <W>.  If the user wishes to use  
## the opposite function, it can be given as an optional second argument.
## A component  .structureconstants is added to  <W>  which  contains the 
## structure constants N_rs with respect to the canonical basis. 
##
## The function 'CanonicalChevalleyBasisRep' is  similar  but returns the 
## matrices with respect to a  given  representation  <repr>   (e.g., see 
## 'LieScRepresentation' etc.)
##
CanonicalChevalleyBasisRep:=function(W,repr)
  local all,struct,i,j,i1,k,ind,r,s,a,x,q;
  all:=[];
  # first positive roots
  for r in [1..W.N] do
    if r<=Length(W.cartan) then
      Add(all,W.epsilon[r]*repr[1][r]);
    else
      i:=1; 
      while not ((W.roots[r]-W.roots[i]) in W.roots) do
        i:=i+1;
      od;
      s:=Position(W.roots,W.roots[r]-W.roots[i]);
      q:=1;
      while W.roots[s]-q*W.roots[i] in W.roots do 
        q:=q+1;
      od;
      a:=W.epsilon[i]*q^(-1)*(all[i]*all[s]-all[s]*all[i]);
      if ForAny(Set(Flat(a)),x->not IsInt(x)) then  
        Error("mist!");
      else
        Add(all,a);
      fi;
    fi;
  od;
  # now negative roots
  for r in [1..W.N] do
    if r<=Length(W.cartan) then 
      Add(all,-W.epsilon[r]*repr[2][r]);
    else
      i:=1;
      while not ((W.roots[W.N+r]+W.roots[i]) in W.roots) do
        i:=i+1;
      od;
      s:=Position(W.roots,W.roots[W.N+r]+W.roots[i]);
      q:=1;
      while W.roots[s]+q*W.roots[i] in W.roots do 
        q:=q+1;
      od;
      a:=-W.epsilon[i]*q^(-1)*(all[W.N+i]*all[s]-all[s]*all[W.N+i]);
      if ForAny(Set(Flat(a)),x->not IsInt(x)) then
        Error("mist!");
      else
        Add(all,a);
      fi;
    fi;
  od;
  if not IsBound(W.structureconstants) then 
    Print("#I Structure constants \c");
    #for s in W.structureconstants do 
    #  if all[s[1]]*all[s[2]]-all[s[2]]*all[s[1]]<>s[3]*all[s[4]] then
    #    Error("nochmal mist!");
    #  fi;
    #od;
    #Print("OK !\n");
    # else 
    struct:=[];
    ind:=[];
    for i in [1..Length(all[1])] do
      for j in [1..Length(all[1])] do
        Add(ind,[i,j]);
      od;
    od;
    for r in [1..Length(W.roots)] do 
      for s in [r+1..Length(W.roots)] do 
        if W.roots[r]+W.roots[s] in W.roots then
          k:=Position(W.roots,W.roots[r]+W.roots[s]);
          i1:=1;
          while all[k][ind[i1][1]][ind[i1][2]]=0 do 
            i1:=i1+1;
          od;
          a:=0;
          for i in [1..Length(all[1][1])] do 
            a:=a+all[r][ind[i1][1]][i]*all[s][i][ind[i1][2]]-
                      all[s][ind[i1][1]][i]*all[r][i][ind[i1][2]];
          od;
          q:=a/all[k][ind[i1][1]][ind[i1][2]];
          #a:=all[r]*all[s]-all[s]*all[r];
          #q:=a[ind[i1][1]][ind[i1][2]]/all[k][ind[i1][1]][ind[i1][2]];
          Add(struct,[r,s,q,k]);
          if Length(struct) mod 1000=0 then 
            Print(Length(struct)," \c");
          fi;
        fi;
      od;
    od;
    W.structureconstants:=struct;
    Print("-> ",Length(struct),"\n");
  fi;
  return all;
end;

CanonicalChevalleyBasis:=function(W)
  return CanonicalChevalleyBasisRep(W,LieAdjointRepresentation(W));
end;

CanonicalChevalleyBasisOld:=function(arg)
  local W,ht0,eps,w0,struct,pr,m,v,x,repr,r,r1,w,pw,i,j,i1,ind,k,l,n,n1,mypr;
  W:=arg[1];
  if Length(arg)=1 then
    eps:=W.epsilon;
  else
    eps:=arg[2];
  fi;
  mypr:=true;
  if Length(arg)=3 then 
    mypr:=arg[3];
  fi;
  repr:=LieAdjointRepresentation(W);
  m:=MonomialSubgroup(W,repr);
  w0:=W.longestword;
  w:=m[1]^0;
  l:=[];
  r:=[];
  v:=[];
  Print("#I Computing root vectors \c");
  for i in [1..W.N] do
    j:=1;
    while w[W.N+1-j][W.N+1-w0[i]]=0 do
      j:=j+1;
    od;
    Add(r,j);
    if w[W.N+1-j][W.N+1-w0[i]]=1 then
      Add(v,eps[w0[i]]);
    elif w[W.N+1-j][W.N+1-w0[i]]=-1 then
      Add(v,-eps[w0[i]]);
    else
      Error("mist!");
    fi;
    Add(l,v[i]*w*repr[1][w0[i]]*w^(-1));
    w:=w*m[w0[i]];
  od;
  # put roots in same order as in W.roots and check
  l:=List([1..W.N],i->l[Position(r,i)]);
  for i in [1..W.N] do
    Add(l,-repr[4]*l[i]*repr[4]);
  od;
  for j in [1..Length(W.cartan)] do 
    for i in [1..Length(l)] do 
      if repr[3][j]*l[i]-l[i]*repr[3][j]<>
                  ScalRootCoRoot(W.cartan,W.roots[i],j)*l[i] then
        Error("mist!");
      fi;
    od;
    Print(".\c");
  od;
  Print(" OK\n");
  if repr[4]^2<>repr[4]^0 then
    Error("mist!");
  fi;
  struct:=[];
  ind:=[];
  for i in [1..Length(repr[4])] do
    for j in [1..Length(repr[4])] do
      Add(ind,[i,j]);
    od;
  od;
  Print("#I epsilon-function : ",eps,"\n");
  ht0:=Sum(W.roots[W.N]);
  for i in [1..Length(l)] do
    for j in [i+1..Length(l)] do
      #if (i<=W.N and j>i and j<=W.N) or (i<=W.N and j<W.N+i and j>W.N) then 
      k:=W.roots[i]+W.roots[j];
      if AbsInt(Sum(k))<=ht0 and k in W.roots then
        k:=Position(W.roots,k);
        pr:=l[i]*l[j]-l[j]*l[i];
        i1:=1;
        while pr[ind[i1][1]][ind[i1][2]]=0 do 
          i1:=i1+1;
        od;
        n:=pr[ind[i1][1]][ind[i1][2]]/l[k][ind[i1][1]][ind[i1][2]];
        if (not IsInt(n)) or pr<>n*l[k] then 
          Error("supermist!");
        fi;
        if mypr=true and i<=W.N and j<=W.N then 
          Print("#I [ ");
          for i1 in W.roots[i] do
            Print(i1,"\c");
          od;
          Print(", ");
          if Sum(W.roots[j])>0 then
            Print(" ");
            for i1 in W.roots[j] do
              Print(i1,"\c");
            od;
          else 
            Print("-");
            for i1 in W.roots[j] do
              Print(-i1,"\c");
            od;
          fi;
          Print(" ] = ");
          if n>0 then
            Print("+",n);
          else
            Print(n);
          fi;
          r1:=RingelRule(W,W.roots[i],W.roots[j]);
          n1:=n;
          if Sum(W.roots[i]) mod 2=1 and Sum(W.roots[j]) mod 2=1 then
            n1:=-n;
          fi;
          if r1[1]=1 then
            Print(" (",r1[2],":",n1,")\n");
          else
            Print(" (",r1[2],":",-n1,")\n");
          fi;
          #Print("\n");
        fi;
        Add(struct,[i,j,n,k]);
      fi;
    od;
  od;
  W.canonicalbasis:=l;
  W.structureconstants:=struct;
  Print("#I ",Length(struct)," structure constants\n");
  return [struct,repr,l];
end;
  
# checks if the basis elements computed in CanonicalChevalleyBasis 
# satisfy relations for canonical basis of adjoint module
CheckCanChevBasis:=function(W,canbase)
  local i,j,ea,a,s,x,f,f1,Ei,Fi,Hi,p,k;
  Ei:=canbase[2][1];
  Fi:=canbase[2][2];
  Hi:=canbase[2][3];
  ea:=canbase[3];
  Print("#I \c");
  for i in [1..Length(W.cartan)] do 
    for j in [1..Length(W.cartan)] do 
      x:=-W.epsilon[j]*Ei[i]*Hi[j]+W.epsilon[j]*Hi[j]*Ei[i];
      if i=j and x<>2*ea[i] then
        Error("mist 1aa!");
      elif i<j and x<>-W.cartan[j][i]*ea[i] then 
        Error("mist 1ab!");
      fi;
      x:=-W.epsilon[j]*Fi[i]*Hi[j]+W.epsilon[j]*Hi[j]*Fi[i];
      if i=j and x<>2*ea[W.N+i] then
        Error("mist 1ba!");
      elif i<j and x<>-W.cartan[j][i]*ea[W.N+i] then 
        Error("mist 1bb!");
      fi;
    od;
    f1:=[];
    for a in [1..2*W.N] do 
      f:=0;
      x:=Ei[i]*ea[a]-ea[a]*Ei[i];
      k:=W.roots[a]+W.roots[i];
      if k in W.roots then  
        f:=1;
        k:=Position(W.roots,k);
        p:=1; 
        while W.roots[a]-p*W.roots[i] in W.roots do
          p:=p+1;
        od;
        if x<>p*ea[k] then 
          Error("mist 2a!");
        fi;
      elif a=W.N+i then 
        f:=2;
        if x<>-W.epsilon[i]*Hi[i] then 
          Error("mist 2b!");
        fi;
      else 
        f:=3;
        if x<>0*x then 
          Error("mist 2c!");
        fi;
      fi;
      if f=0 then
        Error("mist!");
      else
        Add(f1,f);
      fi;
      x:=Fi[i]*ea[a]-ea[a]*Fi[i];
      k:=W.roots[a]-W.roots[i];
      f:=0;
      if k in W.roots then  
        f:=1;
        k:=Position(W.roots,k);
        p:=1; 
        while W.roots[a]+p*W.roots[i] in W.roots do
          p:=p+1;
        od;
        if x<>p*ea[k] then 
          Error("mist 3a!");
        fi;
      elif a=i then 
        f:=2;
        if x<>-W.epsilon[i]*Hi[i] then 
          Error("mist 3b!");
        fi;
      else 
        f:=3;
        if x<>0*x then 
          Error("mist 2c!");
        fi;
      fi;
      if f=0 then
        Error("mist!");
      else
        Add(f1,f);
      fi;
    od;
    for j in Set(f1) do 
      Print(j);
    od;
    Print(" \c");
  od;
  Print(": \c");
  for s in W.structureconstants do 
    if ea[s[1]]*ea[s[2]]-ea[s[2]]*ea[s[1]]<>s[3]*ea[s[4]] then
      Error("nochmal mist!");
    fi;
  od;
  Print("\n#I Canonical basis relations true\n");
  return true;
end;

#########################################################################
##
#F LieNrs( <W>, <r>, <s> ) . . . . . . . . single structure constant with
## . . . . . . . . . . . . . . . . . respect to canonical Chevalley basis
##
## 'LieNrs'  returns the structure constant <nrs> such that 
##                        [e_r,e_s]=nrs*e_{r+s} if r,s and r+s are roots.
##
LieNrs:=function(W,r,s)
  local b,str,nrs,r1,s1,l;
  if not IsBound(W.structureconstants) then 
    b:=CanonicalChevalleyBasis(W);
  fi;
  str:=W.structureconstants;
  if not (r+s in W.roots) then 
    nrs:=0;
  else
    r1:=Position(W.roots,r);
    s1:=Position(W.roots,s);
    if r1<s1 then 
      l:=1;
      while str[l][1]<>r1 or str[l][2]<>s1 do 
        l:=l+1;
      od;
      nrs:=str[l][3];
    else
      l:=1;
      while str[l][1]<>s1 or str[l][2]<>r1 do 
        l:=l+1;
      od;
      nrs:=-str[l][3];
    fi;
    if nrs=0 then
      Error("mist!");
    fi;
  fi;
  return nrs;
end;

#########################################################################
##
#F ChevalleyCommRels( <W>, <r>, <s> ) . . . . . . . . structure constants
## . . . . . . . . . . . . . . . . . .  in Chevalley commutator relations
##
## 'ChevCommRelats'  returns  the structure constants  C(i,j,r,s) in  the 
## Chevalley commutator  relations  for  two  linearly independent  roots 
## <r>, <s> with respect to the canonical Chevalley basis. The formula is
## as follows (see Carter, Simple groups of Lie type, p.76):  For  t,u in 
## the base ring, we have 
##
##     x_s(u)^(-1)x_r(t)^(-1)x_s(u)x_r(t) 
##                       = prod_{i,j>0} x_{ir+js}(C(i,j,r,s)t^iu^j)
## 
## (Note that a minus sign in the Carter formula is absorbed into C.)
##

## Carter, p.76: C(i,1,r,s)=M(r,s,i), C(1,j,r,s)=(-1)^j*M(s,r,j)
## C(3,2,r,s)=M(r+s,r,2)/3, C(2,3,r,s)=-2*M(s+r,s,2)/3 where 
## M(r,s,i)=N(r,s)N(r,r+s)...N(r,(i-1)r+s)/i! First, helper function
LieMrsi:=function(W,r,s,i)
  local r1,s1;
  r1:=W.roots[r];
  s1:=W.roots[s];
  if i=1 then 
    return LieNrs(W,r1,s1);
  elif i=2 then 
    return LieNrs(W,r1,s1)*LieNrs(W,r1,r1+s1)/2;
  elif i=3 then 
    return LieNrs(W,r1,s1)*LieNrs(W,r1,r1+s1)*LieNrs(W,r1,2*r1+s1)/6;
  else 
    Error("mist!");
  fi;
end;

ChevalleyCommRels:=function(W,r,s)
  local cr,m,i,j,t;
  if not (W.roots[r]+W.roots[s] in W.roots) then 
    return [];
  fi;
  cr:=[];
  m:=LieMrsi(W,r,s,1);
  if m<>0 then 
    Add(cr,[[1,1],-m]);
  fi;
  m:=LieMrsi(W,s,r,2);
  if m<>0 then 
    Add(cr,[[1,2],-m]);
  fi;
  m:=LieMrsi(W,r,s,2);
  if m<>0 then 
    Add(cr,[[2,1],m]);
  fi;
  m:=LieMrsi(W,s,r,3);
  if m<>0 then 
    Add(cr,[[1,3],m]);
  fi;
  m:=LieMrsi(W,r,s,3);
  if m<>0 then 
    Add(cr,[[3,1],-m]);
  fi;
  t:=Position(W.roots,W.roots[r]+W.roots[s]);
  if IsInt(t) then 
    m:=LieMrsi(W,t,s,2);
    if m<>0 then 
      Add(cr,[[2,3],-2*m/3]);
    fi;
    m:=LieMrsi(W,t,r,2);
    if m<>0 then 
      Add(cr,[[3,2],-m/3]);
    fi;
  fi;
  return cr;
end;

NilExponential:=function(a,t)
  local i,x,a1;
  i:=1;
  x:=a^0;
  a1:=t*a;
  while a1<>0*a1 do 
    x:=x+a1;
    i:=i+1;
    a1:=i^(-1)*t*a1*a;
  od;
  return x;
end;

# check commutator relations with explicit matrix realisation
CheckCommRels:=function(W,r,s)
  local str,xr,xs,l,a1,a2,t,u,tt;
  if not IsBound(W.canonicalbasis) then 
    str:=CanonicalChevalleyBasisOld(W)[3];
  else
    str:=W.canonicalbasis;
  fi;
  if W.roots[r]+W.roots[s] in W.roots then 
    for tt in [[2,19],[1009,211],[23,10007]] do 
      t:=tt[1];
      u:=tt[2];
      xr:=NilExponential(str[r],t);
      xs:=NilExponential(str[s],u);
      a1:=xs^(-1)*xr^(-1)*xs*xr;
      a2:=a1^0;
      for l in ChevalleyCommRels(W,r,s) do 
        a2:=a2*NilExponential(str[Position(W.roots,
          l[1][1]*W.roots[r]+l[1][2]*W.roots[s])],l[2]*t^l[1][1]*u^l[1][2]);
      od;
     # Print([a1,a2],"\n");
      if a1<>a2 then 
        Error([r,s],"\n");
      fi;
    od;
    return a1=a2;
  else
   return true;
  fi;
end;

#########################################################################
##
#F ChevalleyRootElement( <W>, <repr>, <r>, <t> ) . . . . . . root element
##
## 'ChevalleyRootElement' returns  the  root  group element x_r(t) in the 
## Chevalley group defined with respect to the representation <repr>.
##
ChevalleyRootElement:=function(W,repr,r,t)
  local weiter,i,a,a1,pow,g,x,p,rp;
  rp:=Position(W.roots,r);
  a:=repr[Position(W.roots,r)];
  a1:=a;
  pow:=[];
  weiter:=true;
  i:=1;
  pow:=[ShallowCopy(a)];
  while weiter do
    a1:=a1*a;
    if a1=0*a1 then
      weiter:=false;
    else
      i:=i+1;
      a1:=i^(-1)*a1;
      Add(pow,a1);
      if ForAny(Set(Flat(a1)),x->not IsInt(x)) then
        Print(Set(Flat(a1)),"\n");
        return false;
      fi;
    fi;
  od;
  g:=IdentityMat(Length(a1),t^0);
  for p in [1..Length(pow)] do
    g:=g+t^p*pow[p];
  od;
  return g;
end;

ChevalleyRootElementAdj:=function(W,r,t)
  local str,e,weiter,i,a,a1,pow,g,x,p;
  if not IsBound(W.canonicalbasis) then 
    str:=CanonicalChevalleyBasisOld(W)[3];
  else
    str:=W.canonicalbasis;
  fi;
  a:=str[Position(W.roots,r)];
  a1:=a;
  pow:=[];
  weiter:=true;
  i:=1;
  pow:=[ShallowCopy(a)];
  while weiter do
    a1:=a1*a;
    if a1=0*a1 then
      weiter:=false;
    else
      i:=i+1;
      a1:=i^(-1)*a1;
      Add(pow,a1);
      if ForAny(Set(Flat(a1)),x->not IsInt(x)) then
        Print(Set(Flat(a1)),"\n");
        return false;
      fi;
    fi;
  od;
  g:=IdentityMat(Length(a1),t^0);
  for p in [1..Length(pow)] do
    g:=g+t^p*pow[p];
  od;
  return g;
end;

# Jordan blocks of a unipotent matrix
JordanBlocks:=function(mat)
  local a,b,l,i,j,r,r1;
  l:=[];
  r:=[];
  a:=mat^0;
  b:=IdentityMat(Length(mat),1);
  for i in [1..Length(mat)] do 
    for j in [i+1..Length(mat)] do 
      b[i][j]:=j-i+1;
    od;
    r[i]:=RankMat(a);
    a:=a*(mat-mat^0);
  od;
  r:=r*TransposedMat(b)^-1;
  r1:=[];
  for i in [1..Length(r)] do 
    for j in [1..r[i]] do 
      Add(r1,i);
    od;
  od;
  return r1;
end;
    
# check all  commutator relations
testcr:=function(W)
  local cr,l,i,j,x;
  cr:=CanonicalChevalleyBasisOld(W);
  l:=[];
  for i in [1..2*W.N] do 
    Print(".\c");
    for j in [1..2*W.N] do 
      x:=CheckCommRels(W,i,j);
      if x=false then
        Print("\n",[i,j],"\c");
      fi;
      Add(l,x);
    od;
  od;
 Print("\n");
  return Set(l);
end;
