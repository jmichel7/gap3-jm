#############################################################################
##
#A  coxeter.g            CHEVIE library          Meinolf Geck, Frank Luebeck, 
#A                                            Jean Michel and Goetz Pfeiffer
##
#Y  Copyright (C) 1992 - 2018  Lehrstuhl D fuer Mathematik, RWTH Aachen, 
#Y  and   University Paris VII.
##
##  This file contains general functions that deal with finite Coxeter groups 
##  represented as permutations of the roots.

#############################################################################
##
#F  CoxeterTypeFromArg( <type>, <rk>, ... )  
# returns a reflectiontype + unparsable part of arg
#  E.g., 'CoxeterTypeFromArg( "F", 4 );'
#  For type I2(n), a third argument is needed: ( "I", 2, n ).
# For types B?, G?, F? a third argument is expected the ratio of the roots
CoxeterTypeFromArg:=function(arg)local res,t,shift; res:=[];
  shift:=function()local t;t:=arg[1];arg:=arg{[2..Length(arg)]};return t;end;
  while Length(arg)>0 and IsString(arg[1]) and arg[1]<>[] do
    t:=rec(series:=shift(),rank:=shift(),operations:=ReflTypeOps); 
    if t.rank>0 then
      if t.rank=1 then t.series:="A";fi;
      if t.series="B" then t.cartanType:=2;
      elif t.series="C" then t.series:="B";t.cartanType:=1;
      elif t.series="Bsym" then t.series:="B";t.cartanType:=ER(2);
      elif t.series="B?" then t.series:="B";t.cartanType:=shift();
      elif t.series="D" and t.rank=3 then t.series:="A";
      elif t.series="G" then t.series:="G";t.cartanType:=1;
      elif t.series="Gsym" then t.series:="G";t.cartanType:=ER(3);
      elif t.series="G?" then t.series:="G";t.cartanType:=shift();
      elif t.series="F" then t.series:="F";t.cartanType:=1;
      elif t.series="Fsym" then t.series:="F";t.cartanType:=ER(2);
      elif t.series="F?" then t.series:="F";t.cartanType:=shift();
      elif t.series="H" and t.rank=2 then t.series:="I";t.bond:=5;
      elif t.series="I" then t.bond:=shift();
	if t.bond mod 2=0 then t.cartanType:=1;fi;
      elif t.series="Isym" then t.bond:=shift();t.series:="I";
	t.cartanType:=E(2*t.bond)+E(2*t.bond)^-1;
      elif t.series="I?" then t.bond:=shift();t.series:="I";
	t.cartanType:=shift();
      fi;
      Add(res,t);
    fi;
  od;
  return [res,arg];
end;

#############################################################################
##
#F  CartanMat(W) or  CartanMat( <type>, <rk>, ... )  
#F    returns Cartan matrix of W or of given type(s) and rank(s)
##
##  If several couples <type>, <rk> are given, the direct sum of the
##  corresponding Cartan matrices is returned.
##
CartanMat:=function (arg) local t, res;
  if IsRec(arg[1]) then return ApplyFunc(Dispatcher("CartanMat"),arg);fi;
  return ApplyFunc(DiagonalMat,List(ApplyFunc(CoxeterTypeFromArg,arg)[1],
    CartanMat));
end;

#############################################################################
##
#F  FiniteCoxeterTypeFromCartanMat( <C> )  . . . . type of a Cartan matrix
##  
##  Returns a ReflectionType (a list of records (series:=,indices:=,rank:=)
##  which  describe irreducible components  of the Cartan  matrix <C>, such
##  that    C{indices}{indices}=CartanMat(series,rank)
##  A field .bond is added for type "I".
##  For series B,G,I,F a field .cartanType is added containing the ratio
##  of one conjugacy class of roots to the other
##  
##  If <C> is not the Cartan matrix of a finite Coxeter group this function
##  returns 'false'. So it can be used to check a matrix for this property.
  
FiniteCoxeterTypeFromCartanMat:=C->List(DecomposedMat(C),
  function(s)local co,m,r,x,l,m1,t,tmp,rev;
    m:=C{s}{s}; r:=Length(m); l:=TwoTree(m);
    if l=false then return false;fi;
    t:=rec(operations:=ReflTypeOps,rank:=r);
    if not IsList(l[Length(l)]) then # line: types A,B,C,F,G,H,I
      m1:=m{l}{l};
      co:=i->m1[i][i+1]*m1[i+1][i];
      rev:=function()l:=Reversed(l);m1:=m{l}{l};end;
      if r=1 then t.series:="A";
      elif r=2 then
	if co(1)=1 then t.series:="A";
	elif co(1)=2 then t.series:="B"; 
          if m1[1][2]=-1 or (m1[1][2]<>-2 and m1[2][1]>m1[1][2]) then rev();fi;
             #B2 preferred to C2
	  t.cartanType:=-m1[1][2];
	elif co(1)=3 then t.series:="G"; 
          if m1[1][2]<>-1 and m1[2][1]>m1[1][2] then rev();fi;
          t.cartanType:=-m1[1][2];
	else x:=NofCyc(co(1));
          if m1[2][1]=-1 or (m1[1][2]<>-1 and m1[2][1]>m1[1][2])then rev();fi;
	  if co(1)=2+E(x)+E(x)^-1 then t.bond:=x; else t.bond:=2*x; fi;
	  if t.bond=1 then return false;fi;
	  t.series:="I";
	  if t.bond mod 2=0 then t.cartanType:=-m1[1][2];fi;
	fi;
      else
        if co(r-1)<>1 then rev();fi;
	if co(1)=1 then
	  if co(2)=1 then t.series:="A";
	  elif r<>4 then return false;
	  else
	    if m1[2][3]<m1[3][2] then rev();fi;
	    t.series:="F"; t.cartanType:=-m1[2][3];
	  fi;
	elif co(1)=2 then t.series:="B"; t.cartanType:=-m1[1][2];
	elif not r in [3,4] then return false;
	elif co(1)=(3+ER(5))/2 then t.series:="H";
	else return false;
	fi; 
      fi;
      t.indices:=l;
    else tmp:=l[1];l:=l{[2..Length(l)]};
      if Length(l[2])=1 then t.series:="D";
       t.indices:=Concatenation(l[1],l[2],[tmp],l[3]);
      elif r in [6..8] then t.series:="E";
       t.indices:=Concatenation([l[2][2],l[1][1],l[2][1],tmp],l[3]);
      else return false;
      fi;
    fi;
    if CartanMat(t)<>m{t.indices}{t.indices} then return false;fi; 
       # countercheck if we get C from type
    t.indices:=s{t.indices};return t;
end);

# JM   10-2-2001  :  I   suppressed  the  sorting   on  type  and  rank  in
# ReflectionType, which was done to make easier to test isomorphism between
# Coxeter groups. This sorting should be confined to 'IsomorphismType'. The
# reason  is  that  it  affected  the  way character tables of products are
# constructed,  which caused all  kinds of trouble  in trying to compute HC
# and Lusztig induction of unipotent characters.

# IsomorphismType(W[, options])
# returns a string describing the isomorphism type of the reflection group
# or coset W. Options can be rec(TeX:=true).
IsomorphismType:=function(arg)local W,t,opt;
  W:=arg[1];if Length(arg)=1 then opt:=rec();else opt:=arg[2];fi;
  t:=Reversed(Collected(List(ReflectionType(W),x->ReflectionName(x,opt))));
  t:=Join(List(t,function(x)
    if x[2]=1 then return x[1];else return SPrint(x[2],x[1]);fi;end),"+");
  if IsSpets(W) then W:=Group(W);fi;
  if W.rank>W.semisimpleRank and IsBound(opt.torus) then
    if t<>"" then PrintToString(t,"+");fi;
    if IsBound(opt.TeX) then
         PrintToString(t,"T_{",W.rank-W.semisimpleRank,"}");
    else PrintToString(t,"T",W.rank-W.semisimpleRank);
    fi;
  fi;
  return String(t);
end;
    
#############################################################################
##
#F  RootsCartan( <C> ) . . the positive root system of the Cartan Matrix <C>
##
RootsCartan := function(C)local R,a,i,v,rootsIrredWeyl;
  # rootsIrredWeyl(<C>) root system for irreducible integral Cartan matrix,
  # in the basis of the simple roots.
  rootsIrredWeyl:=function(C)local R, I, l, r, p, old, new, j;
    l:=Length(C); I:=C^0; R:=[]; old:=I;
    while old<>[] do Add(R,old); new:=[];
      for r in old do for j in [1..l] do p:=C[j]*r;
	if p<0 or Length(R)-p-1>0 and r[j]>p and 
            r-(p+1)*I[j] in R[Length(R)-p-1] then 
          if not r+I[j] in new then Add(new,r+I[j]);fi;
	fi;
        old:=new;
      od; od;
    od;
    return Concatenation(R);# R[i] holds positive roots of height i.  
  end;
  if Length(C)=0 then return [];fi;
  if ForAll(Concatenation(C),IsInt) then
    R:=Concatenation(List(DecomposedMat(C),a->List(rootsIrredWeyl(C{a}{a}),
      function(x)local v;v:=0*C[1]; v{a}:=x; return v;end)));
    SortBy(R,x->[Sum(x),-x]);
  else R:=C^0;
    for a  in R  do for i  in [1..Length(C)] do
      if a <> R[i] then v:=a-R[i]*(a*C[i]);
	if not v in R then Add(R,v);fi;
      fi;
    od; od;
  fi;
  return R;
end;

########################################################################
##
#F  HighestShortRoot( <W> ) . . . . . . . . . . . . . highest short root 
##
##  <W>  should be  irreducible. Returns  the unique  short root of maximal
##  height;  if all roots have the same length then this is the unique root
##  of maximal height, which can be obtained in general by W.roots[W.N].
##
HighestShortRoot:=function(W)
  if Length(ReflectionType(W))>1 then Error(W," should be irreducible\n");fi;
  return First([W.N,W.N-1..1],i->W.rootLengths[
    W.rootRestriction[W.orbitRepresentative[i]]]=1);
end;

########################################################################
##
#F  BadPrimes( <W> ) . . . . . . . . . . . . . bad primes
##
##
BadPrimes:=function(W)local l;
  if IsList(W) then 
    l:=Set(Flat(W)); l:=Filtered(l,x->x>1);
    l:=Set(Concatenation(List(l,Factors)));
    return l;
  else return BadPrimes(W.roots);
  fi;
end;

#############################################################################
##
#F  SimpleRootsSubsystem( <W>, <l> ) . . . . . . . simple roots for subsystem
#F  of reflection subgroup
##  
##  <l>  must be  a subset  of [1..2*W.parentN].  Returns the subset of l 
#   indexing simple positive roots of ReflectionSubgroup(Parent(<W>),<l>).
##  
SimpleRootsSubsystem:=function(W,l)local orb, gen, tmp, i, j, n, s;
  l:=Set(l);
  if IsBound(W.parent) then W:=Parent(W); fi;
  
  # trivial case (fast handling of parabolic subgroups):
  if IsSubset([1..W.semisimpleRank],l) then return l; fi;
  
  gen:=List(l,i->Reflection(W,W.rootRestriction[i]));
  
  # compute orbit:
  tmp:=Set(Concatenation(PermGroupOps.Orbits(rec(generators:=gen),l,OnPoints)));
  orb:=Filtered(tmp,x->x<=W.N);
  
  tmp:=Set([]);
  
  # first check input (get the result faster, if input is set of simple roots):
  for i in [1..Length(l)] do
    if l[i]<=W.N then
      n:=0;
      for j in orb do if n<2 and j^gen[i]>W.N then n:=n+1; fi; od;
      if n=1 then AddSet(tmp,l[i]); fi;
    fi;  
  od;
  if tmp=l then return tmp; fi;
  
  for i in Difference(orb,l) do
    n:=0; s:=Reflection(W,W.rootRestriction[i]);
    for j in orb do if n<2 and j^s>W.N then n:=n+1; fi; od;
    if n=1 then AddSet(tmp,i); fi;
  od;
  
  return tmp;
end;

#############################################################################
##
#V  CoxeterGroupOps  . . . . . .  operation record for finite Coxeter groups
##
CoxeterGroupOps := OperationsRecord("CoxeterGroupOps",ComplexGroupOps);
Inherit(CoxeterGroupOps,AbsCoxOps); # this overrides from HasTypeOps 
# BraidRelations and ParabolicRepresentatives, from PermRootOps
# CartanMat, Elements, EltWord, ReducedInRightCoset, ReflectionSubgroup and
# ReflectionType, from PermRootOps CartanMat and from ComplexGroupOps Hecke.
# We want to keep:
CoxeterGroupOps.name:="CoxeterGroupOps";
CoxeterGroupOps.CartanMat:=PermRootOps.CartanMat;
CoxeterGroupOps.BraidRelations:=HasTypeOps.BraidRelations;
CoxeterGroupOps.ParabolicRepresentatives:=HasTypeOps.ParabolicRepresentatives;

CoxeterGroupOps.Dimension:=W->Rank(W)+2*W.N; # dimension of algebraic group

#############################################################################
##
#F  CoxeterLength( <W> , <w> )  . . . . . length of a permutation w in W
##
##  Counts the number of positive roots sent by w to negative roots.
##  
CoxeterGroupOps.CoxeterLength:=function(W,w)local res,N,r;
# return Number(OnTuples(W.rootInclusion{[1..W.N]},w),r->r>W.parentN);
  res:=0;N:=W.parentN; # Code below is 4 times faster than above line
  for r in OnTuples(W.rootInclusion{[1..W.N]},w) do
    if r>N then res:=res+1;fi;
  od;
  return res;
end;

# w is a permutation or a reduced Coxeter word. Returns inversions of w=
# list of roots of Parent(W) sent by w to negative roots
# if w is a word this list is ordered so that a+b is between a and b
CoxeterGroupOps.Inversions:=function(W,w) 
  if IsPerm(w) then return Filtered(W.rootInclusion{[1..W.N]},r->r^w>W.parentN);
  else return List([1..Length(w)],
                    i->W.rootInclusion[w[i]]^EltWord(W,w{[i-1,i-2..1]}));
  fi;
end;

# returns a basis of X(T)\otimes\Q formed of roots and a basis of X(T)^W
CoxeterGroupOps.BaseX:=function(W)
  if not IsBound(W.forMatX) then
    W.forMatX:=W.simpleRoots{[1..W.semisimpleRank]};
    if W.semisimpleRank<W.rank then
      if Length(W.forMatX)=0 then W.forMatX:=IdentityMat(W.rank);
      else Append(W.forMatX,NullspaceMat(TransposedMat(W.simpleCoroots)));
      fi;
    fi;
  fi;
  return W.forMatX;
end;

#############################################################################
##
#F  MatXPerm( <W>, <w> )  . . . . . . . . . . . . convert a permutation to
#F  corresponding matrix operating on X
##  
##  Let <w> be a permutation with the following property: The images of the
##  simple roots of <W> ( i.e., <W>.rootInclusion{<W>.generatingReflections})
##  must be roots of  <W> ( i.e., in the set <W>.rootInclusion ).
##  
##  Let X_1 be the sublattice of X consisting of the elements orthogonal to
##  all  coroots of  <W>. 'MatXPerm'  returns the  unique invertible matrix
##  which, as matrix acting on X, maps the simple roots of <W> as indicated
##  by <w> and which induces the identity map on X_1.
##  
##  If in particular <w> is  an element of  <W> then 'MatXPerm( <W>, <w> )
##  is  the matrix   representing <w> as  linear map   on X. This  follows
##  immediately from the formula for the generating reflections of <W>.
##  
CoxeterGroupOps.MatXPerm:=function(W, w)local X, M, I;
  if W.semisimpleRank=0 then return IdentityMat(W.rank);fi;
  X:=CoxeterGroupOps.BaseX(W);
  I:=[W.semisimpleRank+1..W.rank];
  M:=List(W.roots{W.rootRestriction{OnTuples(W.rootInclusion{
            [1..W.semisimpleRank]},w)}},i->Concatenation(i,I*0));
  return Concatenation(M,IdentityMat(W.rank){I})^X;
end;

# The permutations of <roots> effected by the reflections determined by cartan
PermutationsSimpleReflections:=function(roots,cartan)local p,N;
  p:=SortingPerm(roots);N:=[1..Length(roots)];
  return List([1..Length(cartan)],function(i)local v;v:=Copy(roots);
    v{N}[i]:=v{N}[i]-roots*cartan[i];return SortingPerm(v)/p;end);
end;

#############################################################################
##
#F  AddComponentsCoxeterGroup(W)  only used internally by 'CoxeterGroup'
#F  'ReflectionSubgroup' and \*
##  
AddComponentsCoxeterGroup:=function(W)local parent, tmp, i, p, f,C;
  if IsBound(W.parent) then parent:=W.parent; else parent:=W; fi;
  W.operations:=CoxeterGroupOps;
  W.isCoxeterGroup := true;
  W.semisimpleRank:=Length(W.cartan);
  W.generatingReflections:=[1..W.semisimpleRank];
  W.nbGeneratingReflections:=W.semisimpleRank;
  W.reflectionsLabels:=W.rootInclusion{[1..W.semisimpleRank]};
  W.N:=Length(W.roots)/2;
  W.parentN:=Length(parent.roots)/2;

  if IsBound(parent.rank) then W.rank:=parent.rank;
  elif W.semisimpleRank>0 then W.rank:=Length(W.simpleRoots[1]);
  else W.rank:=0;
  fi;

  # for each root in parent group we determine a simple root in the same orbit
  tmp:=PointsAndRepresentativesOrbits(W,2*W.parentN);
  W.orbitRepresentative:=[];
  W.orbitRepresentativeElements:=[];
  for i in [1..Length(tmp[1])] do 
    if tmp[1][i][1] in W.rootInclusion then
      W.orbitRepresentative{W.rootRestriction{tmp[1][i]}}
                                                :=0*tmp[1][i]+tmp[1][i][1];
      W.orbitRepresentativeElements{W.rootRestriction{tmp[1][i]}}
                                                               :=tmp[2][i];
    fi;
  od;

  # root lengths from parent group, if already known:
  if IsBound(parent.rootLengths) then
    W.rootLengths:=parent.rootLengths{parent.orbitRepresentative{
     W.rootInclusion{W.generatingReflections}}};
  else
  f:=function(t)local i,j,k,c;
    c:=CartanMat(t);
    for j in [1..Length(t.indices)] do for k in [1..j-1] do
      if c[j][k]<>c[k][j] then
        for i in t.indices do 
	  if W.orbitRepresentative[i]=W.orbitRepresentative[t.indices[j]]
	  then W.rootLengths[i]:=-c[k][j];
	  else W.rootLengths[i]:=-c[j][k];
	  fi;
	od;
	return;
      fi;   
    od;od;
    W.rootLengths{t.indices}:=List(t.indices,i->1);
  end;
  # root lengths (for parent group) using the classification:
    W.rootLengths:=[];
    for tmp in ReflectionType(W) do f(tmp);od;
  fi;
  
  # the set of coroots: i-th entry is coroot to roots[i]
  C:=TransposedMat(W.cartan);
  if W.cartan=C then W.coroots:=W.roots;
  elif not IsBound(parent.coroots) then
    # the dual group as permutation group on the coroots:
    tmp:=rec(roots:=RootsCartan(C));Append(tmp.roots,-tmp.roots);
    tmp.reflections:=PermutationsSimpleReflections(tmp.roots,C);
    W.coroots:=[];
    for i in [1..W.N] do
      p:=W.orbitRepresentative[i];
      p:=p^Product(tmp.reflections{CoxeterWord(W,
         W.orbitRepresentativeElements[i])});
      Add(W.coroots,tmp.roots[p]);
    od;
    Append(W.coroots,-W.coroots);
  else
    W.coroots:=[];
    for p in RootsCartan(C) do
      W.coroots[Position(parent.coroots{W.rootInclusion},p*
        parent.coroots{W.rootInclusion{[1..W.semisimpleRank]}})]:=p;
    od;
    Append(W.coroots,-W.coroots);
  fi;

  # Order of the Coxeter group 
  # (This doesn't use a stabilizer chain but allows a much faster 
  # construction of a stabilizer chain)
  # here also the 'ReflectionType' is computed.
  if ReflectionType(W)=false then Error("Unknown Cartan Matrix\n");
  else Size(W);
  fi;
  if not IsBound(W.name) then W.name:=ReflectionName(W);fi;
  W.matgens:=List(W.reflections{W.generatingReflections},x->MatXPerm(W,x));
end;

#############################################################################
##  
#F  CoxeterGroup( <cartan>)
#F  CoxeterGroup( <simpleRoots>, <simpleCoroots>)
#F  CoxeterGroup( <type>, <rank>[, <type>, <rank> ...][, "sc"])
#F  CoxeterGroup( <rec> )
##  
CoxeterGroup:=function(arg)local l, sc, W, typ, a, tmp;

  l:=Length(arg);
  
  if l=1 and IsRec(arg[1]) then # get back the .coxeter entry of a record:
    if IsBound(arg[1].coxeter) then return arg[1].coxeter;
    elif IsBound(arg[1].operations) and IsBound(arg[1].operations.CoxeterGroup)
    then return arg[1].operations.CoxeterGroup(arg[1]);
    fi;
  fi;
  
  if l>0 and arg[l]="sc" then sc:=true; l:=l-1;
  else sc:=false;
  fi;

  W:=rec();
  if l=2 and IsMat(arg[1]) then
    if sc then 
      Error("CoxeterGroup: argument \"sc\" not allowed ",
            "together with 2 matrices\n");
    fi;
    W.simpleRoots:=arg[1]; W.simpleCoroots:=arg[2];
    W.cartan:=W.simpleCoroots*TransposedMat(W.simpleRoots);
  else
    if l=0 then W.cartan:=[];
    elif IsString(arg[1]) or IsRec(arg[1]) then
	W.cartan:=ApplyFunc(CartanMat,arg{[1..l]});
    else W.cartan:=arg[1];
    fi;
    if Length(W.cartan)=0 then W.simpleRoots:=[];W.simpleCoroots:=[];
    elif sc then
      W.simpleRoots:=TransposedMat(W.cartan); W.simpleCoroots:=W.cartan^0;
    else
      W.simpleRoots:=W.cartan^0; W.simpleCoroots:=W.cartan;
    fi;
  fi;

  W.roots:=RootsCartan(W.cartan);Append(W.roots,-W.roots);
  tmp:=Group(PermutationsSimpleReflections(W.roots,W.cartan),());
  Inherit(tmp,W);W:=tmp;
  W.reflections:=ShallowCopy(W.generators);
  W.rootInclusion := [1..Length(W.roots)];
  W.rootRestriction:=[1..Length(W.roots)];
  AddComponentsCoxeterGroup(W);

  ###### for printing simplify matrix arguments: 
  # check if one matrix in argument is superfluous:
  if l=2 and W.rank=W.semisimpleRank and IsMat(arg[1]) then
    if arg[1]=arg[1]^0 then arg:=arg{[2..Length(arg)]}; l:=1;
    elif arg[2]=arg[2]^0 then arg[1]:=TransposedMat(arg[1]);arg[2]:="sc";l:=1;
    fi;
  fi;
  # check in case of one matrix if it can be created by 'CartanMat':
  if l=1 and W.cartan=ApplyFunc(DiagonalMat,List(ReflectionType(W),CartanMat)) 
  then W.name:=SPrint("CoxeterGroup(",Join(Concatenation(
    List(ReflectionType(W),t->ReflectionName(t,rec(arg:=true))),
    arg{[2..Length(arg)]})),")");
  else W.name:=SPrint("CoxeterGroup(",Join(List(arg,FormatGAP)),")");
  fi;
  ###### end of printing section

  AbsCoxOps.CompleteCoxeterGroupRecord(W);
  Inherit(W.operations,CoxeterGroupOps);
  return W;
end;


#############################################################################
##
#F  ReflectionSubgroup( <W>, <J>)
##  

CoxeterGroupOps.ReflectionSubgroup:=function(W,J)
  local simpleroots, tmp, res, i, n, r, pos;
  
  if IsBound(W.parent) then W:=W.parent;fi;
  
# generators must be given as indices of reflections of the parent group
  if IsSubset(W.reflectionsLabels,J) then
    J:=List(J,x->W.rootInclusion[Position(W.reflectionsLabels,x)]);
  elif not IsSubset(W.rootInclusion,J) then
    Error("argument should be labels for simple reflexions or a subset of W.rootInclusion");
  fi;
  
  res:=CHEVIE.GetCached(W,"ReflectionSubgroups",rec(callarg:=[J]),
      x->x.callarg);#check if subgroup in cache
  if IsBound(res.isGroup) then return res;fi;
  
  # determine a set of simple roots:
  simpleroots:=SimpleRootsSubsystem(W,J);
  if Set(J)=simpleroots and Length(J)=Length(Set(J)) then simpleroots:=J;fi;

  # The command 'Subgroup' returns W, if Set(J)=Set(W.generators).
  # To avoid this we substitute it by some more lines.
  Inherit(res,Group(()));
  res.parent     := W;
  res.generators := List(simpleroots,j->Reflection(W,j));
  for i in [1..Length(res.generators)] do res.(i):=res.generators[i]; od;
 
  if simpleroots=[] then # case of empty simpleroots:
    res.simpleRoots:=[]; res.simpleCoroots:=[];
    res.cartan:=[]; res.roots:=[]; r:=[];
  else
    res.simpleRoots:=W.roots{simpleroots}*W.simpleRoots;
    res.simpleCoroots:=W.coroots{simpleroots}*W.simpleCoroots;
    res.cartan:=W.coroots{simpleroots}*
                        W.cartan*TransposedMat(W.roots{simpleroots});
    res.roots:=RootsCartan(res.cartan);Append(res.roots,-res.roots);
    r:=res.roots*W.roots{simpleroots};
  fi;
  res.rootInclusion:=ListBlist([1..2*W.N],BlistList(W.roots,r));
  res.rootInclusion:=Permuted(res.rootInclusion,
                            PermListList(r,W.roots{res.rootInclusion}));
  res.rootRestriction:=[];
  res.rootRestriction{res.rootInclusion}:=[1..Length(res.roots)];
  res.reflections:=List(res.rootInclusion,i->Reflection(W,i));
  AddComponentsCoxeterGroup(res);
  AbsCoxOps.CompleteCoxeterGroupRecord(res);
  Inherit(res.operations,CoxeterGroupOps);
  if res.semisimpleRank=0 then
       res.name:=SPrint("ReflectionSubgroup(",W,", [])");
  else res.name:=SPrint("ReflectionSubgroup(",W,", ",
			res.rootInclusion{res.generatingReflections},")");
  fi;
  return res;
end;

#############################################################################
##
#F  LeftDescentSet( W, x)  the set of generators s such that sx < x
##  
##  the generators are numbered by the corresponding reflections of the parent
##  group
##  
CoxeterGroupOps.LeftDescentSet:=function(W,w) 
 return Filtered(W.rootInclusion{W.generatingReflections},i->i^w>W.parentN);
end;

#############################################################################
##
#F  CoxeterGroupOps.ReflectionType( <W> ) . . . . . type for CoxeterGroup record
##  
##  This is essentially the same as ReflectionType(<W>.cartan), but with
##  two exceptions: 
##  A subsystem  of type "A" or "D" is said to  be of type "~A" or "~D",
##  iff the following conditions are fulfilled:
##  
##   - <W> is not a parent group
##   - the subsystem is contained in an irreducible subsystem of 
##     Parent(<W>), which has long and short roots
##   - the subsystem consists of short roots
##  
CoxeterGroupOps.ReflectionType:= function(W)local p, a, b, res;
  res:=ReflectionType(W.cartan);
  
  if false in res then InfoChevie("# Unknown Cartan matrix\n");return false;fi;

  # parent or all roots of parent have same length:
  if not IsBound(W.parent) or Length(W.parent.rootLengths)=0 or
    Maximum(W.parent.rootLengths)=1 then return res;
  fi;
  
  for a in res do
    p:=Parent(W).orbitRepresentative[W.rootInclusion[a.indices[1]]];
    if a.series in ["A","D"] and Parent(W).rootLengths[p]=1 then
      for b in ReflectionType(Parent(W)) do
        if p in b.indices and b.series in ["B","C","F","G"] then
          a.tilde:=true;
        fi;
      od;
    fi;
  od;
  
  return res;
end;
        
#############################################################################
##
#F  CoxeterGroupOps.\=( <W1>, <W2> )  . . . equality test for Coxeter groups
##  
##  Here  two CoxeterGroup records  are defined to  be equal if  they are 
##  equal as permutation  groups, and if the simple  roots and  the simple 
##  coroots are equal.
##  
if false then
CoxeterGroupOps.\= := function(W1,W2)
  return IsCoxeterGroup(W1) and IsCoxeterGroup(W2) and
    W1.simpleRoots=W2.simpleRoots and W1.simpleCoroots=W2.simpleCoroots 
    and PermGroupOps.\=(W1,W2);
end;
fi;

CoxeterGroupOps.ReflectionCoDegrees:=W->ReflectionDegrees(W)-2;

#############################################################################
##
#F  ReflectionCharValue( <W>, <w> ) ... The reflection character on w
##
##  'ReflectionCharValue' returns  the value of the reflection character of W 
##  at the element w.
##  
CoxeterGroupOps.ReflectionCharValue:=function(W,w)
  return W.rank-W.semisimpleRank+Sum([1..W.semisimpleRank],
     i->W.roots[W.rootRestriction[W.rootInclusion[i]^w]][i]);
end;

CoxeterGroupOps.IsLeftDescending:=function(W,w,i)
   return W.rootInclusion[i]^w>W.parentN; end;

CoxeterGroupOps.FirstLeftDescending:=function(W,w)local i;
  for i in W.generatingReflections do
    if W.rootInclusion[i]^w>W.parentN then return i;fi;
  od;
  return false;
end;

#############################################################################
##
#F  PermMatX( <W> , <M> )  . . . . . . . . .  convert a matrix in X which
#F  preserves the roots to permutation of the roots of Parent(<W>)
#F  PermMatXCoxeterElement( <W> , <M> )  . . . . convert a matrix in X which
#F  represents an element of Parent(<W>) to corresponding permutation
#F  PermMatY( <W> , <M> )  . . . . . . . . .  convert a matrix in Y which
#F  preserves the coroots to permutation of the roots of Parent(<W>)
##  

# this is faster than 'PermMatX', but works for elements of the Coxeter 
# group only:
PermMatXCoxeterElement:=function(W, M) local tmp;
  
  # we always have to look at the parent group:
  if IsBound(W.parent) then W:=W.parent; fi;

  # the images of the points Base(W) determine the element uniquely:
  tmp:=W.roots*W.simpleRoots;
  tmp:=List(Base(W),p->Position(tmp,tmp[p]*M));
  return RepresentativeOperation(W,Base(W),tmp,OnTuples);
end;

CoxeterGroupOps.PermMatX:=function(W, M) local tmp;
  # we always have to look at the parent group:
  if IsBound(W.parent) and W<>W.parent then
    return W.parent.operations.PermMatX(W.parent,M);
  fi;
  if W.semisimpleRank=0 then return ();fi;
  tmp:=W.roots*W.simpleRoots;
  return PermListList(tmp,tmp*M);
end;

PermMatY:=function(W, M)
  return PermMatX(W, TransposedMat(M))^-1;
end;

CoxeterGroupOps.ProductRootEmbed:=function(W1,W2)local p;
  p:=SortingPerm(List(Concatenation(W1.roots,W2.roots),
       function(x)local s;s:=Sum(x);if not IsInt(s) then s:=evalf(s);fi;
         if s>0 then return s; else return W1.N+W2.N-s;fi;end));
  return [OnTuples([1..2*W1.N],p),OnTuples(2*W1.N+[1..2*W2.N],p)];
end;

CoxeterGroupOps.\*:=function(W1,W2)local g;
  if IsGroup(W2) then
    g:=ComplexGroupOps.\*(W1,W2);
    if not IsCoxeterGroup(W1) then return g;fi;
    g.cartan:=DiagonalMat(W1.cartan,W2.cartan);
    g.simpleRoots:=Concatenation(
       List(W1.simpleRoots,x->Concatenation(x,[1..W2.rank]*0)),
       List(W2.simpleRoots,x->Concatenation([1..W1.rank]*0,x)));
    g.simpleCoroots:=Concatenation(
       List(W1.simpleCoroots,x->Concatenation(x,[1..W2.rank]*0)),
       List(W2.simpleCoroots,x->Concatenation([1..W1.rank]*0,x)));
    AddComponentsCoxeterGroup(g);
    AbsCoxOps.CompleteCoxeterGroupRecord(g);
    Inherit(g.operations,CoxeterGroupOps);
    return g;
  else return CoxeterCoset(W1,W2);
  fi;
end;

#############################################################################
##
#F  ElementWithInversions(<W>,<N>)
#
# given the set N of positive roots of W negated by an element w, find w.
# N is a subset of [1..W.N] (not W.parentN!)
#
ElementWithInversions:=function(W,N)local r,p,w,n;
  w:=W.identity;n:=N;
  while Length(n)>0 do
    p:=Intersection(n,W.generatingReflections);
    if Length(p)=0 then return false;fi;
    r:=Reflection(W,p[1]);
    n:=W.rootRestriction{OnTuples(W.rootInclusion{Difference(n,p{[1]})},r)};
    w:=r*w;
  od;
  return w^-1;
end;

#############################################################################
##
#F  DescribeInvolution(<W>,<w>)
#
# Given  an involution w  there is a  unique H=ReflectionSubgroup(W,I) such
# that  w=LongestCoxeterElement(H)  and  w  is  central  in H. The function
# returns I. For now does not work for abscox groups.
#
DescribeInvolution:=function(W,w)
  return SimpleRootsSubsystem(W,Filtered(W.rootInclusion{[1..W.N]},
                                     i->i^w=i+W.parentN));
end;

#############################################################################
##
#F  StandardParabolic(<W>,<H>) H reflection subgroup or its simple roots.
#
# Given parabolic H returns w such that H^w is a standard parabolic subgroup
# of W
#
CoxeterGroupOps.StandardParabolic:=function(W,hr)local b,N,w;
  if not IsList(hr) then hr:=hr.rootInclusion{hr.generatingReflections};fi;
  if hr=[] then return ();fi;
  b:=W.roots{W.rootRestriction{hr}};
  b:=Concatenation(ListBlist(W.roots{W.generatingReflections},
    List(SemiEchelonMat(b).ishead,x->not x)),b);
# complete basis of I with part of S to make basis
  N:=ListBlist([1..W.N],List(W.roots{[1..W.N]}*b^-1,function(v)local x; 
   for x in v do if (IsRat(x) and x<0) or evalf(x)<0 then return true;
           elif (IsRat(x) and x>0) or evalf(x)>0 then return false;fi;od;end));
# find negative roots for associated order and make order standard
  w:=ElementWithInversions(W,N);
  if IsSubset(W.rootInclusion{W.generatingReflections},OnTuples(hr,w)) then
    return w;
  else return false;
  fi;
end;

#############################################################################
##
#F  RelativeGroup(<W>,<J>)
#
#   W  is a  Coxeter group;  if S=W.rootInclusion{W.generatingReflections},
#   then J should be a *distinguished* subset of S, that is if for s\in S-J
#   we  set v(s,J)=w_0^{J\cup s}w_0^J then J  is stable by all v(s,J). Then
#   R=N_W(W_J)/W_J  is a Coxeter group with  Coxeter system the v(s,J). The
#   program  return  R  in  its  reflection  representation  on X(ZL_J/ZG).
#   (according to Lusztig's "Coxeter Orbits...", the images of the roots of
#   W in X(ZL_J/ZG) form a root system).
#
#   R  is has the fields:
#     .relativeIndices:=Filtered(S,x->not x in J)
#     .parentMap:=the list of v(s,J)
#     .MappingFromNormalizer maps J-reduced elements of N_W(W_J) to
#         elements of R
#
CoxeterGroupOps.RelativeGroup:=function(W,J)local res,qr,S,I,vI,r;
  S:=W.rootInclusion{W.generatingReflections};
  if not IsSubset(S,J) then
    Error("implemented only for standard parabolic subgroups");
  fi;
  I:=Filtered(S,x->not x in J); # not Difference(S,J): keep the order
  vI:=List(I,i->LongestCoxeterElement(W,Concatenation([i],J))*
     LongestCoxeterElement(W,J));
  if ForAny(vI,g->OnSets(Set(J),g)<>Set(J)) then 
    Error(J," is not distinguished in ",W,"\n");
  fi;
  qr:=i->W.roots[W.rootRestriction[i]]{W.rootRestriction{I}};
  if J=[] then res:=W;
  else res:=CoxeterGroup(List(I,i->List(I,j->
       ProportionalityCoefficient(qr(j)-qr(j^vI[Position(I,i)]),qr(i)))));
  fi;
  res.relativeIndices:=I;
  res.parentMap:=vI;
  res.MappingFromNormalizer:= # maps reduced elts of N_W(W_J) to elts of res
    function(w)local c,r; c:=();
      while true do
	r:=PositionProperty(I,x->x^w>W.parentN);
        if r=false then return c;fi;
 	w:=vI[r]*w;
	c:=c*res.generators[r];
      od;
    end;
  return res;
end;

CoxeterGroupOps.Hecke:=function(arg)local H;
  H:=ApplyFunc(PermRootOps.Hecke,arg);
  H.operations:=HeckeAlgebraOps;
  Inherit(H.operations,CoxeterHeckeAlgebraOps);
  H.operations.CoxeterGroup:=x->x.reflectionGroup;
  return H;
end;

# Parabolic subgroups are represented by the indices of their generating roots.
ParabolicSubgroups:=function(W) local l;
  if IsBound(W.parabolicSubgroups) then return W.parabolicSubgroups;fi;
  l:=ParabolicRepresentatives(W);
  W.parabolicSubgroups:=Concatenation(List(l,function(x) local r; 
    r:=ReducedRightCosetRepresentatives(W,ReflectionSubgroup(W,x));
    return Set(List(r,y->Set(OnTuples(x,y))));end));
  return W.parabolicSubgroups;
end;

#############################################################################
#  Extended reflection groups by diagram automorphisms
#
# Fields phis: group of automorphisms as permutations
#         F0s: group of automorphisms as matrices
#
ExtendedGroupOps:=OperationsRecord("ExtendedGroupOps");

ExtendedGroupOps.ReflectionName:=function(W,opt)local g,ff;
  g:=W.F0s;
  if g=[] then return ReflectionName(W.group,opt);fi;
  if Length(g)=1 then return ReflectionName(Spets(W.group,g[1]),opt);fi;
  if ForAll(W.F0s,x->x^2=x^0) and Size(ApplyFunc(Group,W.F0s))=6 then
    return SPrint(ReflectionName(W.group,opt),".S3");
  fi;  
  ff:=List(W.phis,x->RestrictedPerm(x,
     W.group.rootInclusion{W.group.generatingReflections}));
  if ForAll(ff,x->x<>()) or Rank(W.group)=0 or ForAll(g,x->x=x^0) then
    return SPrint("Extended(",ReflectionName(W.group,opt),",",Join(ff),")");
  fi;
  g:=List(g,x->ReflectionName(Spets(W.group,x),opt));
  return SPrint("Extended(",Join(g),")");
end;

ExtendedGroupOps.Format:=ExtendedGroupOps.ReflectionName;

ExtendedGroupOps.String:=W->ReflectionName(W,rec());

ExtendedGroupOps.Print:=function(W)Print(String(W));end;

IsExtendedGroup:=x->IsRec(x) and IsBound(x.operations) and x.operations=
  ExtendedGroupOps;

ExtendedGroupOps.\*:=function(a,b)
  return ExtendedReflectionGroup(a.group*b.group,Concatenation(
    List(a.F0s,m->DiagonalMat(m,IdentityMat(b.group.rank))),
    List(b.F0s,m->DiagonalMat(IdentityMat(a.group.rank),m))));
end;

#############################################################################
#F  ExtendedReflectionGroup( <W>[, <automorphisms>])
##  
ExtendedReflectionGroup:=function(arg)local res,W,om,tmp;
  res:=rec(group:=arg[1],operations:=ExtendedGroupOps);
  if Length(arg)=1 then 
    if IsCoxeterCoset(arg[1]) then res.group:=Group(arg[1]);
         res.phis:=[arg[1].phi];res.F0s:=[arg[1].F0Mat];
    else res.phis:=[];res.F0s:=[];
    fi;
  elif IsPerm(arg[2]) then res.phis:=[arg[2]];
  elif IsGroup(arg[2]) then 
    if IsPermGroup(arg[2]) then res.phis:=arg[2].generators;
    elif IsMatGroup(arg[2]) then res.F0s:=arg[2].generators;
    fi;
  elif IsList(arg[2]) then
    if arg[2]=[] then res.phis:=[];res.F0s:=[];
    elif IsPerm(arg[2][1]) then res.phis:=arg[2];
    elif IsMat(arg[2][1]) or arg[2][1]=[] then res.F0s:=arg[2];
    elif IsMat(arg[2]) then res.F0s:=[arg[2]];
    fi;
  else Error("second argument should be a group of diagram automorphisms");
  fi;
  
  W:=Parent(res.group);
  if not IsBound(res.F0s) then
    if W.semisimpleRank<> W.rank then
      Error("#I automorphisms can be permutations ",
	    "only for semisimple parent datum.\n");
    fi;
    res.F0s:=List(res.phis,x->MatXPerm(W,x));
  fi;

  if W.roots=[] then res.phis:=List(res.F0s,M->());return res;fi;
  if not ForAll(res.F0s,M->IsNormalizing(W.roots*W.simpleRoots,M) and
    IsNormalizing(W.coroots*W.simpleCoroots,TransposedMat(M))) then
    Error("#I automorphisms should normalize roots and coroots.\n");
  fi;
  tmp:=W.roots*W.simpleRoots;
  res.phis:=List(res.F0s,M->PermListList(tmp,tmp*M));
  return res;
end;

ReadChv("prg/wclsinv"); # for CoxeterGroupOps.ClassInvariants
ReadChv("prg/semisimple"); # for CoxeterGroupOps.FundamentalGroup,
 # QuasiIsolatedRepresentatives, SemisimpleCentralizer, IsIsolated
