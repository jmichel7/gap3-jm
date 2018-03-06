###########################################################################
# Contribution to the Chevie Package
#
#  This file contains some supplementary programs for working with braids
#    (C) Jean MICHEL 1995-2007
###########################################################################
##
#F  LeftDivisors(b)   returns all left divisors of Graside element b
##
#   Example:
#   gap> B:=Braid(CoxeterGroup("A",2)); 
#   function ( arg ) ... end
#   gap> LeftDivisors(B(1,2));
#   [ ., 1, 12 ]
##
LeftDivisors:=function(b)local s,M,w;
  M:=b.monoid;if b=M.Elt([]) then return [b];fi;
  w:=GarsideAlpha(b);
  s:=Filtered([1..M.nrAtoms],i->M.IsLeftDescending(w,i));
  s:=List(s,x->M.Elt([M.atoms[x]]));
  return Union([M.Elt([])],s,Union(List(s,x->x*LeftDivisors(x^-1*b))));
end;

RightDivisors:=b->List(LeftDivisors(b.monoid.Reverse(b)),b.monoid.Reverse);

########################################################################
##
#F  PiRoots(wF,d) this function returns all reduced d-th roots of pi in the
##   braid monoid of the Coxeter coset WF using a very bestial algorithm
##
#   Example:
#   gap> WF:=RootDatum("3D4");
#   3D4
#   gap> PiRoots(WF,12);
#   [ 31, 43, 23, 32, 34, 13 ]
##
PiRoots:=function(WF,d)local W,e,pi,F,p,phiorder;
  phiorder:=function(w,p)
    return First(DivisorsInt(OrderPerm(w*p)),i->(w*p)^i=p^i);
  end;
  if IsBound(WF.phi) then W:=CoxeterGroup(WF);p:=WF.phi;F:=Frobenius(WF);
		     else W:=WF;p:=();F:=x->x;
  fi;
  e:=CoxeterElements(W,2*W.N/d);
  e:=List(Filtered(e,x->phiorder(x,p)=d),Braid(W));
  pi:=Braid(W)(LongestCoxeterElement(W))^2;
  e:=Filtered(e,x->TwistedPower(d,x,F)=pi);
  return e;
end;

########################################################################
#
#F  CheckGeckKimPfeiffer(WF) Check that all roots of pi are cyclically
##   conjugate in the Braid monoid of the Coxeter coset WF
##
CheckGeckKimPfeiffer:=function(WF)local F,r,d,rr,n,pi,b,r,W,red,bad,rbad;
  if IsBound(WF.phi) then W:=CoxeterGroup(WF);F:=Frobenius(WF);
  else W:=WF;F:=x->x;fi;
  r:=RegularEigenvalues(WF);
  r:=r{[2..Length(r)]};
  pi:=Braid(W)(LongestCoxeterElement(W))^2;
  for d in r do
    rr:=PiRoots(WF,d);
    Print("d=",d," reduced roots:",Length(rr));
    rr:=ConjugacySet(rr[1],F,"Cyc");
    Print(" cycclass:", Length(rr));
    bad:=0;
    rbad:=0;
    red:=0;
    for b in rr do
      if Length(b.elm)<2 then red:=red+1;fi;
      if TwistedPower(d,b,F)<>pi then 
          if Length(b.elm)<2 then rbad:=rbad+1;fi;
	  bad:=bad+1;
      fi;
    od;
    Print(" reduced:",red," bad:",bad," rbad:",rbad,"\n");
  od;
end;

########################################################################
#
#F  AlphaI(b,I) find the longest prefix of Garside element b using only
##    b.monoid.atoms{I}
#   Example:
#   gap> W:=CoxeterGroup("A",4);
#   CoxeterGroup("A",4)
#   gap> w0:=Braid(W)(LongestCoxeterElement(W));
#   w0
#   gap> AlphaI(w0,[1,2,3]);
#   121321
##
AlphaI:=function(b,I)local M,res,i,s;
  M:=b.monoid;
  res:=M.Elt([]);
  i:=1;
  while i<=Length(I)  do
    if M.IsLeftDescending(GarsideAlpha(b),I[i]) then
      s:=M.Elt([M.atoms[I[i]]]);res:=res*s; b:=s^-1*b; i:=1;
    else i:=i+1;
    fi;
  od;
  return res;
end;

########################################################################
#
#F  AllWords(b)  find all the decompositions in atoms of Garside element b
#   Example
#   gap> W:=CoxeterGroup("A",2);                  
#   CoxeterGroup("A",2)
#   gap> pi:=Braid(W)(LongestCoxeterElement(W))^2;
#   w0.w0
#   gap> AllWords(pi);
#   [ [ 1, 1, 1, 2, 2, 1 ], [ 1, 2, 1, 1, 1, 2 ], [ 1, 2, 1, 2, 2, 1 ], 
#     [ 1, 2, 2, 2, 1, 1 ], [ 2, 1, 1, 1, 2, 2 ], [ 2, 1, 2, 1, 1, 2 ], 
#     [ 2, 1, 2, 2, 2, 1 ], [ 2, 2, 2, 1, 1, 2 ] ]
AllWords:=function(b)local M,s;
  M:=b.monoid;
  if b=M.Elt([]) then return [[]];fi;
  s:=GarsideAlpha(b);
  return Concatenation(List(Filtered([1..M.nrAtoms],i->M.IsLeftDescending(s,i)),
    i->List(AllWords(M.Elt([M.LeftQuotient(s,M.atoms[i])])*
    GarsideOmega(b)),d->Concatenation([i],d))));
end;

#############################################################################
# Hurwitz action of an element b in B_{Length(l)} on the sequence l of
#  group elements.
# b is a list of integers representing an element of B_{Length(l)}, or an
# integer i representating the list [i].
# Example:
#   gap> HurwitzAction([1,2,3,4,5],[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)]);
#   [ (2,3), (3,4), (4,5), (5,6), (6,7), (1,7) ]
#
HurwitzAction:=function(b,l)local i;
  if IsInt(b) then b:=[b];fi;
  l:=ShallowCopy(l);
  for i in b do
    if i>0 then l{[i,i+1]}:=[l[i+1],l[i]^l[i+1]];
    else i:=-i; l{[i,i+1]}:=[l[i]*l[i+1]*l[i]^-1,l[i]];
    fi;
  od;
  return l;
end;

#############################################################################
# Let l be a list of group elements.
# Returns [set of all group elements which appear as an item in one of
# the lists in the Hurwitz orbit of the list l, Hurwitz orbit of l]
# to find the dual atoms apply it to an expression for a Coxeter element
# Example:
#   gap> HurwitzOrbitItems([(1,2),(2,3)]);      
#   #I Hurwitz orbit: 1 items=3
#   #I Hurwitz orbit: 2 items=3
#   #I Hurwitz orbit: 3 items=3
#   [ [ (2,3), (1,2), (1,3) ], 
#     [ [ (2,3), (1,3) ], [ (1,2), (2,3) ], [ (1,3), (1,2) ] ] ]
HurwitzOrbitItems:=function(l)local orbit,new,old,r,refs;
  orbit:=[]; old:=[l];r:=Length(l); refs:=Set(l);
  repeat
    new:=Set(Concatenation(List([1..r-1],i->List(old,e->HurwitzAction(i,e)))));
    UniteSet(orbit,old);old:=Difference(new,orbit);
    UniteSet(refs,Set(Concatenation(old)));
    InfoChevie("#I Hurwitz orbit: ",Length(orbit)," items=",Length(refs),"\n");
  until Length(old)=0;
  return [refs,orbit];
end;

#############################################################################
##
#F  DualRelations
##  The relations of the dual braid monoid M
##
DualRelations:=function(M)local res,i,j,w,p,l,rels;
  res:=[];
  for i in [1..M.nrAtoms] do
    w:=M.atoms[i]^-1*M.delta;
    for j in [1..M.nrAtoms] do
      if M.IsLeftDescending(w,j) then Add(res,[i,j]);fi;
    od;
  od;
  p:=List(res,x->Product(M.atoms{x}));
  l:=Set(p);
  rels:=List(l,x->[]);
  for i in [1..Length(res)] do Add(rels[Position(l,p[i])],res[i]);od;
  return rels;
end;

PoincareDualSimples:=function(W,w)local d,M;
  d:=ReflectionDegrees(W);
  InfoChevie("#  predicted=",Product(d+Maximum(d))/Product(d),"\n");
  M:=DualBraidMonoid(W,w);
  d:=LeftDivisorsSimple(M,M.delta);
  return Sum([1..Length(d)],i->X(Rationals)^(i-1)*Length(d[i]));
end;

#############################################################################
#   Jean Michel and François Digne 2012.
#  Additional functions to study conjugation in braid monoids.

# returns [cycling(b),simple r] such that b^r=cycling(b)
Cycling:=function(b) local M,w,l,res;
  M:=b.monoid;l:=Length(b.elm);
  if l=0 then return [b,M.identity];fi;
  w:=M.DeltaAction(b.elm[1],-b.pd);
  if Length(b.elm)=1 then return [M.Elt([w],b.pd),w];fi;
  return
    [GarsideEltOps.Normalize(M.Elt(M.AddToNormal(b.elm{[2..l]},w),b.pd)),w];
end;

# returns [decycling(b),simple r] such that b^r=decycling(b)
DeCycling:=function(b) local M,w;
  M:=b.monoid;
  if Length(b.elm)=0 then return [b,M.Elt([])];fi;
  w:=M.Elt([b.elm[Length(b.elm)]])^-1;
  return [b^w,w];
end;

# returns [sliding(b),r] such that b^r=sliding(b)
Sliding:=function(b)local r;
  r:=PreferredPrefix(b);
  return [PositiveSimpleConjugation(b,r),r];
end;

# random braid of length l in garside monoid M
RandomBraid:=function(M,l)
  return M.B(List([1..l],x->Random([1..M.nrAtoms])));end;
