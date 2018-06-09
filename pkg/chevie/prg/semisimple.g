#############################################################################
##
#A  semisimple.g           CHEVIE library      Cedric Bonnafe and Jean Michel
##
#Y  Copyright (C) 2004 - 2010  Lehrstuhl D fur Mathematik, RWTH Aachen,
#Y  University  Paris VII and Universite de Franche-Comte
##
##  This file contains functions dealing with semisimple elements of algebraic
##  groups.
##
## A reductive group G of rank n is represented in CHEVIE by
## CoxeterGroup(R,C)  where R  is the  inregral matrix  whose lines are the
## roots expressed in the canonical basis of X(T) and C the integral matrix
## whose lines are the coroots expressed in the canonical basis of Y(T).
##
## A finite order semisimple element is represented as an element of
## Y(T)\otimes Q/Z=(Q/Z)^n, that is a list of length n of rationals r
## such that 0<=r<1.
##
##  A semisimple element record contains the following fields:
##  .v     The defining list of rationals
##  .group The parent of the defining group

MvpOps.Mod1:=p->p; # horrible -- fix it
#############################################################################
##
#F  SemisimpleElementOps . . . . .  operations record for semisimple elements
##  
##  
SemisimpleElementOps:=OperationsRecord("SemisimpleElementOps");

# semisimple elements have a domain to be usable as group elements
# SemisimpleElement(W,v[,true or false])
SemisimpleElement:=function(arg)local W,v,res;
  W:=arg[1];v:=arg[2];
  if IsRec(W) then W:=Parent(W);fi;
  res:=rec(group:=W, operations:=SemisimpleElementOps,
    domain:=rec(operations:=rec(
      \in:=function(a,b)return a.operations=SemisimpleElementOps;end,
      Group:=function(D,gens,id) return rec(isDomain:=true,isGroup:=true,
	identity:=id,generators:=gens,operations:=GroupOps,isFinite:=true);
      end)));
  if Length(arg)=3 then res.additive:=arg[3];
  else res.additive:=ForAll(v,IsRat);fi;
  if res.additive then v:=Mod1(v);fi;
  res.v:=v;
  return res;
end;

IsSemisimpleElement:=x->IsRec(x) and IsBound(x.operations) and x.operations=
  SemisimpleElementOps;

SemisimpleElementOps.\*:=function(a,b)
  if IsList(a) then return List(a,x->x*b);
  elif IsList(b) then return List(b,x->a*x);
  elif a.additive then return SemisimpleElement(a.group,a.v+b.v,true);
  else return SemisimpleElement(a.group,
                 Zip(a.v,b.v,function(x,y)return x*y;end),false);
  fi;
end;

SemisimpleElementOps.\/:=function(a,b)return a*b^-1;end;

SemisimpleElementOps.\=:=function(a,b)return a.v=b.v;end;

SemisimpleElementOps.\<:=function(a,b)return a.v<b.v;end;

# n can be:
# Integer --- raise to that power
# Permutation --- act by elt of Coxeter group or coset
# Matrix --- act by element of GL(X(T))
# List --- it defines an element of X(T) in the basis of roots;
#          return the pairing
# SemisimpleElement --- return s if s also; else make n act on s
SemisimpleElementOps.\^:=function(s,n)
  if not IsSemisimpleElement(s) then return s.operations.\^(s,n);fi;
  if IsInt(n) then # raise to a power
    if s.additive then return SemisimpleElement(s.group,n*s.v,true);
          else return SemisimpleElement(s.group,List(s.v,x->x^n),false);
    fi;
  elif IsPerm(n) then # act by an element of W or WF
    return s^MatYPerm(s.group,n^-1);
  elif IsMat(n) then # act by an element of GL(X(T))
    if s.additive then return SemisimpleElement(s.group,s.v*n,true);
          else return SemisimpleElement(s.group,List(TransposedMat(n),
             v->Product(Zip(v,s.v,function(a,b)return b^a;end))),false);
    fi;
  elif IsSemisimpleElement(n) then return s;
  elif IsList(n) then # act by an element of the root lattice
    if s.additive then return Mod1(n*s.group.simpleRoots*s.v);
          else return Product(Zip(n*s.group.simpleRoots,s.v,
                            function(a,b)return b^a;end));
    fi;
  else Error("does not know how to make ",n," act on a semisimple element");
  fi;
end;

SemisimpleElementOps.Comm:=function(a,b)return a/a;end;

SemisimpleElementOps.Frobenius:=function(W,s,i)return s^(W.phi^i);end;

SemisimpleElementOps.String:=x->SPrint("<",Join(x.v),">");

SemisimpleElementOps.Print:=function(x)Print(String(x));end;

SemisimpleElementOps.Value:=function(s,x)
  s:=ShallowCopy(s);s.v:=Value(s.v,x);return s;end;

SubTorusOps:=OperationsRecord("SubTorusOps");

SubTorus:=function(arg)local W,V; W:=arg[1];
  if Length(arg)=1 then V:=IdentityMat(Rank(W));else V:=arg[2];fi;
  V:=ComplementIntMat(IdentityMat(Rank(W)),V);
  if ForAny(V.moduli,x->x<>1) then 
    Error("not a pure sublattice");return false;
  fi;
  return rec(generators:=V.sub,complement:=V.complement,group:=W,
    operations:=SubTorusOps);
end;

SubTorusOps.String:=T->SPrint("SubTorus(",T.group,",",
  T.generators,")");

SubTorusOps.Print:=function(T)Print(String(T));end;

SubTorusOps.Rank:=T->Length(T.generators);

# element in subtorus
SubTorusOps.\in:=function(s,T)local n,i,v,r,V; 
  n:=Lcm(List(s.v,Denominator)); s:=s.v*n;
  V:=List(T.generators,x->List(x,y->y mod n));
  i:=1;
  for v in Filtered(V,x->x<>0*x) do
    while v[i]=0 do if s[i]<>0 then return false;else i:=i+1;fi;od;
    r:=Gcdex(n,v[i]);
    v:=List(r.coeff2*v,x->x mod n);
    if s[i] mod v[i]<>0 then return false;
    else s:=s-s[i]/v[i]*v;s:=List(s,x->x mod n);
    fi;
  od;
  return s=0*s;
end;

# returns (Tso,s-stable representatives of T/Tso) for automorphism s of T
# here m is the matrix of s on Y(T)
# use ss 1.2(1): Ker(1+m+m^2+...)/Im(m-Id)
FixedPoints:=function(T,m)local n,fix,Y1,o;
  n:=List(T.generators*m,z->SolutionMat(T.generators,z)); #action on subtorus
  if false in n or not ForAll(Flat(n),IsInt) then 
    Error(m," does not stabilize ",T);
  fi;
  fix:=NullspaceIntMat(n-n^0); # pure sublattice Y(Tso)
  o:=OrderMat(n);
  Y1:=NullspaceIntMat(Sum([0..o-1],i->n^i));# pure sublattice Y(T1) where 
  # T=T1.Tso almost direct product, thus spaces Y.(1-s) and Y1.(1-s) coincide
  n:=BaseIntMat(n-n^0); # basis of Im(1-s)
  m:=List(Y1,v->SolutionMat(n,v));# basis of Im[(1-s)^{-1} restricted to Y1]
  # generates elements y of Y1\otimes\BQ such that (1-s)y\in Y1
  return [SubTorus(T.group,fix*T.generators),AbelianGenerators(List(m,v->
    SemisimpleElement(T.group,v*Y1*T.generators)))];
end;

#############################################################################
##
#F  AlgebraicCentre( <W> )  . . . centre of algebraic group W
##  
##  <W>  should be a Weyl group record  (or an extended Weyl group record).
##  The  function returns information  about the centre  Z of the algebraic
##  group defined by <W> as a record with fields:
##   Z0:         subtorus Z^0
##   complement: S=complement torus of Z0 in T
##   AZ:         representatives of Z/Z^0 given as a group of ss elts
##   [implemented only for connected groups 18/1/2010]
##   [I added something hopefully correct in general. JM 22/3/2010]
##   [introduced subtori JM 2017 and corrected AZ computation]
##   descAZ:  describes AZ as a quotient  of the fundamental group Pi (seen
##   as  the centre of  the simply connected  goup with same isogeny type).
##   Returns words in the generators of Pi which generate the kernel of the
##   map Pi->AZ
##
AlgebraicCentre:=function(W)local Z0,res,F0s,w,toAZ,hom,AZ,m,id;
  if IsExtendedGroup(W) then F0s:=List(W.F0s,TransposedMat); W:=W.group;fi;
  if W.simpleRoots=[] then Z0:=IdentityMat(W.rank);
  else Z0:=NullspaceIntMat(TransposedMat(W.simpleRoots));
  fi;
  Z0:=SubTorus(W,Z0);
  if Length(Z0.complement)=0 then AZ:=[];
  else AZ:=(Z0.complement*TransposedMat(W.simpleRoots))^-1*Z0.complement;
  fi;
  AZ:=List(AZ,x->SemisimpleElement(W,x));
  if IsBound(F0s) then # compute fixed space of F0s in Y(T)
    for m in F0s do
      AZ:=Filtered(AZ,s->(s/SemisimpleElement(W,s.v*m)) in Z0);
      if Rank(Z0)>0 then Z0:=FixedPoints(Z0,m);Append(AZ,Z0[2]);Z0:=Z0[1];fi;
    od;
  fi;
  id:=SemisimpleElement(W,[1..W.rank]*0);
  res:=rec(Z0:=Z0,AZ:=Group(AbelianGenerators(AZ),id));
  if IsBound(F0s) and Length(F0s)>0 then return res;fi;
  AZ:=List(WeightInfo(W).CenterSimplyConnected,x->SemisimpleElement(W,x));
  if Length(AZ)=0 then res.descAZ:=AZ;return res;fi;
  AZ:=ApplyFunc(Group,AZ);
  toAZ:=function(s)
    s:=s.v*W.simpleCoroots;
    s:=SolutionMat(Concatenation(res.Z0.complement,res.Z0.generators),
       Concatenation(s,[1..Length(res.Z0.generators)]*0));
    return SemisimpleElement(W,s{[1..W.semisimpleRank]}*res.Z0.complement);
  end;
  # map of root data Y(Wsc)->Y(W)
  hom:=GroupHomomorphismByImages(AZ,res.AZ,AZ.generators,
    List(AZ.generators,toAZ));
  res.descAZ:=List(Kernel(hom).generators,x->GetWord(AZ,x));
  return res;
end;

#############################################################################
##
#F  SemisimpleSubgroup( <W>, <V>, <n> )  . . . n-torsion of SubTorus(W,V)
#F  SemisimpleSubgroup( <T>, <n> )  . . . n-torsion of torus T
##
##  Returns the subgroup of semisimple elements of order dividing n in the 
##  subtorus of T represented by V, an integral basis of a sublattice of Y(T).
#
SemisimpleSubgroup:=function(arg)local T,n;
  T:=arg[1];n:=arg[Length(arg)];
  if IsGroup(T) then return
    ApplyFunc(Group,List(arg[2],v->SemisimpleElement(T,v/n)));
  else return 
    ApplyFunc(Group,List(T.generators,v->SemisimpleElement(T.group,v/n)));
  fi;
end;

#############################################################################
##
#F  Centralizer( <W>, <s>)  . . . Stabilizer of <s> in <W>
##  
##  Returns  the stabilizer  of the  semisimple element  <s> in  <W>, which
##  describes  also C_G(s),  if G  is the  algebraic group described by the
##  Weyl  group record W.  The stabilizer is  an extended reflection group,
##  with the reflection group part equal to the Weyl group of C_G^0(s), and
##  the diagram automorphism part being those induced by C_G(s)/C_G^0(s) on
##  C_G^0(s). It is accepted that <W> itself is an extended group.
##
CoxeterGroupOps.Centralizer:=function(W,s)local p,W0s,N,totalW;
  if s in Parent(W) or IsGroup(s) then return PermGroupOps.Centralizer(W,s);fi;
  if  not IsSemisimpleElement(s) then
    Error(s," must be an element of Parent(W) or a semisimple element");
  fi;
  if IsExtendedGroup(W) then 
    totalW:=Subgroup(Parent(W.group),Concatenation(W.group.generators,W.phis));
    W:=W.group;
  else totalW:=W;
  fi;
  if s.additive then 
       p:=Filtered(W.rootInclusion{[1..W.N]},i->s^Parent(W).roots[i]=0);
  else p:=Filtered(W.rootInclusion{[1..W.N]},i->s^Parent(W).roots[i]=s.v[1]^0);
  fi;
  W0s:=ReflectionSubgroup(W,p);
  N:=Normalizer(totalW,W0s);
  N:=Filtered(List(LeftCosets(N,W0s),Representative),w->s=s^w);
  N:=List(N,x->ReducedInRightCoset(W0s,x));
  N:=Subgroup(W,AbelianGenerators(N));
  if W.rank<>W.semisimpleRank then
    if Length(N.generators)=0 then N:=Group(W.matgens[1]^0);
    else N:=ApplyFunc(Group,List(N.generators,x->MatXPerm(W,x)));
    fi;
  fi;
  return ExtendedReflectionGroup(W0s,N.generators);
end;

#############################################################################
##
#F  IsQuasiIsolated( <W>, <s>)  . . . . whether <s> is quasi-isolated in <W>
##
IsQuasiIsolated:=function(W,s)
  return Rank(AlgebraicCentre(Centralizer(W,s)).Z0)=
    Rank(W)-SemisimpleRank(W);
end;

#############################################################################
##
#F  IsIsolated( <W>, <s>)  . . . . whether <s> is isolated in <W>
CoxeterGroupOps.IsIsolated:=function(W,s)
  return Rank(AlgebraicCentre(Centralizer(W,s).group).Z0)=
    Rank(W)-SemisimpleRank(W);
end;

#############################################################################
##
#F  FundamentalGroup(<W>)  fundamental group of algebraic group <W>
##
## This  function  returns  the  fundamental  group  of the algebraic group
## corresponding  to the Weyl group record <W> as the diagram automorphisms
## of  the  corresponding  affine  Weyl  group  induced  by  <W>, thus as a
## permutation group on the simple roots and the negative longest roots
## of each irreducible component.
## The fundamental group is defined as (P^\vee\cap Y(T))/Q^\vee
#
CoxeterGroupOps.FundamentalGroup:=function(W)local n,l,e,omega,r,moved,iszero;
  if W.cartan=[] then return Group(());fi;
  omega:=Mod1(W.cartan^-1*W.simpleCoroots);# simple coweights in basis of Y(T)
  Add(omega,omega[1]*0); # add a"zero" weight
  iszero:=x->x=x*0;
  l:=List(W.type,function(t)local n,r;n:=t.indices;
    # next line  uses that negative roots are listed by decreasing height!
    r:=First([2*W.N,2*W.N-1..1],i->Sum(W.roots[i]{n})<>0);
    return [n,r,
       Concatenation(Filtered(n,i->W.roots[r][i]=-1),[Length(omega)])];end);
  e:=Cartesian(List(l,x->x[3]));
  e:=Filtered(e,x->iszero(Mod1(Sum(omega{x}))));
  e:=List(e,x->Product([1..Length(x)],i->
    LongestCoxeterElement(W,W.rootInclusion{l[i][1]})*
       LongestCoxeterElement(W,W.rootInclusion{Difference(l[i][1],[x[i]])})));
  moved:=Concatenation(List(l,x->x[1]));Append(moved,List(l,x->x[2]));
  e:=List(e,p->RestrictedPerm(p,W.rootInclusion{moved}));
  return Group(AbelianGenerators(e),());
end;

CoxeterGroupOps.Dual:=W->CoxeterGroup(W.simpleCoroots,W.simpleRoots);

ExtendedGroupOps.FundamentalGroup:=W->FundamentalGroup(W.group);

#############################################################################
##
#F  QuasiIsolatedRepresentatives(<W>[,<p>]) . . representatives of W-orbits of 
##       quasi-isolated semisimple elements.
##
##  This function follows Theorem 4.6 in 
##  C.Bonnafe, ``Quasi-Isolated Elements in Reductive Groups''
##  Comm. in Algebra 33 (2005), 2315--2337
##  after one fixes the following bug: at the beginning of section 4.B
##  ``the stabilizer of $\Omega\cap\tilde\Delta_i$ in $\cal A_G$ acts
##    transitively on $\Omega\cap\tilde\Delta_i$''
##  should be
##  ``the stabilizer of $\Omega$ in $\cal A_G$ acts
##    transitively on $\Omega\cap\tilde\Delta_i$''
##
CoxeterGroupOps.QuasiIsolatedRepresentatives:=function(arg)
  local W,p,H,res,iso,w,ind,Z;
  W:=arg[1];
  if W.semisimpleRank=0 then return [SemisimpleElement(W,[1..W.rank]*0)];fi;
  H:=FundamentalGroup(W);
  if Length(arg)=1 then p:=0;else p:=arg[2];fi;
  iso:=W.cartan^-1*W.simpleCoroots; # coweights
  w:=[]; ind:=[];
  res:=List(Cartesian(List(W.type,
    function(t)local n,r,d,p; n:=t.indices; # n is \Delta_t
    # next line  uses that negative roots are listed by decreasing height!
    r:=First([2*W.N,2*W.N-1..1],i->Sum(W.roots[i]{n})<>0);
    d:=W.rootInclusion{Concatenation(n,[r])}; # d is \tilde\Delta_t
    Add(ind,d);
    Add(w,Concatenation(Zip(iso{n},-W.roots[r]{n},function(x,y)return
         x/y;end),[0*iso[1]]));
    p:=Concatenation(List([1..Size(H)],i->Combinations(d,i)));
    return Filtered(p,P->Length(Orbits(Stabilizer(H,P,OnSets),P))=1); 
      # possible sets \Omega_t
    end)),Concatenation);
  res:=Filtered(res, function(P)local S;S:=Stabilizer(H,P,OnSets);
    return ForAll(ind,I->Length(Orbits(S,Intersection(P,I)))=1);end);
  res:=List(Orbits(H,List(res,Set),OnSets),x->x[1]);# possible sets \Omega
  if p<>0 then 
    res:=Filtered(res,P->ForAll(Zip(ind,w,
    function(I,W)local J;J:=Intersection(P,I); 
      return Length(J) mod p<>0 and ForAll(W{List(J,x->Position(I,x))},
        v->Lcm(List(v,Denominator))mod p<>0);end),x->x));
  fi;
  res:=List(res,P->Sum(Zip(ind,w,
    function(I,p)local J;J:=Intersection(P,I); 
      return Sum(p{List(J,x->Position(I,x))})/Length(J);end)));
  res:=Set(List(res,s->SemisimpleElement(W,Mod1(s))));
  Z:=AlgebraicCentre(W).Z0;
  if Rank(Z)>0 then 
    res:=res{Filtered([1..Length(res)],i->not ForAny([1..i-1],j->
     res[i]/res[j] in Z))};
  fi;
  return res;
end;

#############################################################################
##
#F  StructureRationalPointsConnectedCenter(<MF>,q)  . . Structure of Z^0(M)^F
#
#   MF is a spets. Gives the abelian invariants of Z^0(M)^F(q)
#
StructureRationalPointsConnectedCentre:=function(MF,q)local M,W,Z0,Phi,Z0F;
  if IsSpets(MF) then M:=Group(MF);
  else M:=MF;MF:=Spets(M);
  fi;
  W:=Parent(M);
  Z0:=AlgebraicCentre(M).Z0;
  Phi:=MatYPerm(W,MF.phi);
  Z0F:=Z0.generators*(Phi*q-Phi^0);
  Z0F:=List(Z0F,x->SolutionIntMat(Z0.generators,x));
  Z0F:=DiagonalOfMat(DiagonalizeIntMat(Z0F).normal);
  return Filtered(Z0F,x->x<>1);
end;

#############################################################################
##  IntermediateGroup(W,I)
##  returns an intermediate semisimple group between the adjoint and simply 
##  connected group with same Cartan matrix as W whose X is generated by the
##  root lattice and the given subset I of minuscule weights
IntermediateGroup:=function(W,I)local C,w,d,R,v;
  C:=CartanMat(W);w:=TransposedMat(C)^-1; # w = weights in terms of roots
  R:=w^0;
  for v in I do
    if IsInt(v) then Add(R,w[v]); else Add(R,Sum(w{v})); fi;
  od;
  d:=Lcm(List(Flat(R),Denominator));R:=BaseIntMat(d*R)/d;
  return CoxeterGroup(R^-1,C*TransposedMat(R));
end;

#############################################################################
##
##   RootDatum(...)
##   return various known CoxeterGroups and CoxeterCosets
##   Thanks to Jay Taylor for csp and gpin.
##
RootDatum:=function(arg)local type,data,res; type:=arg[1];
  data:=rec();
  data.gl:=function(dim)local R;
    if dim=1 then return Torus(1);fi;
    R:=IdentityMat(dim);R:=R{[1..dim-1]}-R{[2..dim]};
    return CoxeterGroup(R,R);
  end;
  data.pgl:=dim->CoxeterGroup("A",dim-1);
  data.sl:=dim->CoxeterGroup("A",dim-1,"sc");
  data.u:=dim->CoxeterCoset(data.gl(dim),List(-IdentityMat(dim),Reversed));
  data.su:=function(dim)if dim=2 then return CoxeterCoset(data.sl(dim));
    else return CoxeterCoset(data.sl(dim),
         Product([1..QuoInt(dim-1,2)],i->(i,dim-i)));fi;end;
  data.psu:=function(dim)if dim=2 then return CoxeterCoset(data.pgl(dim));
    else return CoxeterCoset(data.pgl(dim),
         Product([1..QuoInt(dim-1,2)],i->(i,dim-i)));fi;end;
  data.sp:=function(dim)local R,R1,i;
    R:=IdentityMat(dim/2); for i in [2..dim/2] do R[i][i-1]:=-1;od;
    R1:=Copy(R);R1[1]:=2*R1[1];
    return CoxeterGroup(R1,R);
  end;
  data.csp:=function(dim)local R,cR;dim:=dim/2;
    R:=IdentityMat(dim+1);R:=R{[1..dim]}-R{[2..dim+1]}; cR:=Copy(R); 
    R[1]{[1,2]}:=[0,-1]; cR[1]{[1,2]}:=[1,-2];
    return CoxeterGroup(cR,R);
  end;
  data.gpin:=function(dim)local R,cR,d;
    d:=QuoInt(dim,2);
    R:=IdentityMat(d+1);R:=R{[1..d]}-R{[2..d+1]}; cR:=Copy(R);
    if dim mod 2=1 then
      R[1]{[1,2]}:=[0,-1]; cR[1]{[1,2]}:=[1,-2];
    else
      R[1]{[1..3]}:=[1,-1,-1]; cR[1]{[1..3]}:=[0,-1,-1];
      R:=List(R,x->Concatenation([0],x));
      cR:=List(cR,x->Concatenation([0],x));
      cR[1][1]:=1;
    fi;
    return CoxeterGroup(R,cR);
  end;
  data.("3gpin8"):=CoxeterCoset(data.gpin(8),[[1,1,1,0,0,0],
   [-2,0,-1,-1,-1,-1],[-1,0,-1,0,0,-1],[-1,0,-1,0,-1,0],
   [-1,0,-1,-1,0,0],[-1,-1,0,0,0,0]]);
  data.("gpin-"):=function(dim)local F,d;
    d:=dim/2;F:=IdentityMat(d+2);
    F[1]{[1..3]}:=[1,-1,1];
    F{[1..d+2]}[2]:=[1..d+2]*0-1;
    F{[2,3]}{[2,3]}:=-IdentityMat(2);
    return CoxeterCoset(data.gpin(dim),F);
  end;
  data.so:=function(dim)local R,R1,i;R:=IdentityMat(QuoInt(dim,2)); 
    for i in [2..Length(R)] do R[i][i-1]:=-1;od;
    if dim mod 2=1 then R1:=Copy(R);R1[1][1]:=2; return CoxeterGroup(R,R1);
    else R[1][2]:=1; return CoxeterGroup(R,R);
    fi;
  end;
  data.("so-"):=dim->CoxeterCoset(data.so(dim),(1,2));
  data.psp:=dim->CoxeterGroup("C",dim/2);
  data.pso:=function(dim)
    if dim mod 2=1 then return CoxeterGroup("B",QuoInt(dim,2));
    else return CoxeterGroup("D",dim/2);
    fi;
  end;
  data.("pso-"):=dim->CoxeterCoset(data.pso(dim),(1,2));
  data.spin:=function(dim)
    if dim mod 2=1 then return CoxeterGroup("B",QuoInt(dim,2),"sc");
    else return CoxeterGroup("D",dim/2,"sc");
    fi;
  end;
  data.("spin-"):=dim->CoxeterCoset(data.spin(dim),(1,2));
  data.halfspin:=function(dim)local R;
    if dim mod 4<>0 then 
      Error("Half-Spin groups only exist in dimension multiple of 4");
    fi;
    dim:=dim/2;R:=IdentityMat(dim);
    R[dim]:=Concatenation([-dim/2,1-dim/2],[2-dim,3-dim..-2],[2]);
    return CoxeterGroup(R,CartanMat("D",dim)*TransposedMat(R^-1));
  end;
  data.2I:=e->CoxeterCoset(CoxeterGroup("Isym",2,e),(1,2));
  data.suzuki:=CoxeterCoset(CoxeterGroup("Bsym",2),(1,2));
  data.2B2:=data.suzuki;
  data.G2:=CoxeterGroup("G",2);
  data.ree:=CoxeterCoset(CoxeterGroup("Gsym",2),(1,2));
  data.2G2:=data.ree;
  data.triality:=CoxeterCoset(CoxeterGroup("D",4),(1,2,4));
  data.3D4:=data.triality;
  data.3D4sc:=CoxeterCoset(CoxeterGroup("D",4,"sc"),(1,2,4));
  data.2E6:=CoxeterCoset(CoxeterGroup("E",6),(1,6)(3,5));
  data.2E6sc:=CoxeterCoset(CoxeterGroup("E",6,"sc"),(1,6)(3,5));
  data.2F4:=CoxeterCoset(CoxeterGroup("Fsym",4),(1,4)(2,3));
  data.F4:=CoxeterGroup("F",4);
# the following is Galois-stable for Aut(H3,[[1,2,1,2,3,2,1,2,1],[3],[2]])
  data.H3stable:=CoxeterGroup([[(5-ER(5))/2,-ER(5),-1+ER(5)],[-2*ER(5),1,0],
    [2*ER(5),1,0]],[[(5+ER(5))/40,(-1-ER(5))/4,(-3+ER(5))/16],
    [(-3*ER(5))/20,1/2,(2-ER(5))/8],[(3*ER(5))/20,1/2,(2+ER(5))/8]]);
# the following model is Galois-stable for Aut(H4,[[1,2,1,2,3,2,1,2,1],
#  [3],[2],[1,2,3,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,3,1,2,1,2,3,4,3,2,1]]);
  data.H4stable:=CoxeterGroup([[(5-ER(5))/2,-ER(5),-1+ER(5),0],[-2*ER(5),1,0,0],
    [2*ER(5),1,0,0],[(5-ER(5))/4,(-1+ER(5))/4,(-3-ER(5))/2,(60-12*ER(5))/25]],
    [[(5+ER(5))/40,(-1-ER(5))/4,(-3+ER(5))/16,0],[(-3*ER(5))/20,1/2,
      (2-ER(5))/8,0],[(3*ER(5))/20,1/2,(2+ER(5))/8,0],[-ER(5)/20,-1/2,
      (-4-ER(5))/8,(25-5*ER(5))/96]]);
  if IsBound(data.(type)) then res:=data.(type);
    if IsFunc(res) then res:=ApplyFunc(res,arg{[2..Length(arg)]});fi;
    res.name:=SPrint("RootDatum(",Join(List(arg,FormatGAP)),")");
    return res;
  fi;
  Print("Unknown type \"",type,"\". Known types are:\n");
  Cut(Join(List(Set(RecFields(data)),FormatGAP)));
  Error("\n");
end;

#############################################################################
##
#F  SemisimpleCentralizerRepresentatives(W[,p])
## Representatives of G-classes of C_G(s)^0.
## Same as W-orbits of subsets of \Pi\cup\{-\alpha_0\}
SemisimpleCentralizerRepresentatives:=function(arg)
  local W,p,cent,ED,J,R,indices,h;
  W:=arg[1];
  if Length(arg)=1 then p:=0;else p:=arg[2];fi;
  return List(Cartesian(List(ReflectionType(W),function(t)local r;
  indices:=W->W.rootInclusion{W.generatingReflections};
  cent:=[];
  r:=Filtered([1..W.N],i->Sum(W.roots[i])=Sum(W.roots[i]{t.indices}));
  h:=List(W.roots{r},Sum);
  ED:=Concatenation(t.indices,[r[Position(h,Maximum(h))]]);
  for J in Combinations(ED) do
    R:=ReflectionSubgroup(W,J);
    if ForAll(cent,G->IsomorphismType(R)<>IsomorphismType(G) or
      RepresentativeOperation(W,indices(R),indices(G),OnSets)=false) then
      Add(cent,R);
    fi;
  od;
  cent:=List(cent,indices);
  if p=0 then return cent;fi;
  return Filtered(cent,I->ForAll(Concatenation(SmithNormalFormMat(W.roots{I})),
     x->x=0 or x mod p<>0));end)),Concatenation);
end;
