#############################################################################
##
#A  abscox.g              CHEVIE library       Frank Luebeck and  Jean Michel
##
#Y  Copyright (C) 1992 - 2010  Lehrstuhl D fuer Mathematik, RWTH Aachen, and
#Y  Universite Paris VII.
##
##  This file contains  GAP functions for working with abstract Coxeter groups.
##
# They are represented by groups W which should have the following fields
# and methods defined:
#
#   .nbGeneratingReflections
#
#   .reflections      such that W.reflections{[1..W.nbGeneratingReflections]}
#                     are the Coxeter generators
#
#   .operations.IsLeftDescending(W,w,i)
#      for w in W and i in [1..W.nbGeneratingReflections], true iff
#      W.reflections[i]*w has a shorter length in the generators than w.
#
# plus the methods contained in

AbsCoxOps:=OperationsRecord("AbsCoxOps");

#############################################################################
# Once the above fields are defined, the function
#     AbsCoxOps.CompleteCoxeterGroupRecord
#
# Can be called to fill in some other fields/methods (the value provided by
# AbsCoxOps.CompleteCoxeterGroupRecord can be overridden by filling first
# the value when it is marked as 'default' below):
#
#   .nbGeneratingReflections
#       default:   Length(W.reflections)
#
#   .reflectionsLabels
#       default:   [1..W.nbGeneratingReflections]
#
#   ReflectionFromName
#       default:   x->Position(W.reflectionsLabels,x);
#
#   isCoxeterGroup
#       true
#
#   .coxeterMat
#       should be pre-filled when there is an infinite bond, otherwise
#       there will be an infinite loop. An infinite bond is marked by 0
#
#
# CompleteCoxeterGroupRecord also fills the following fields:
#
#   .EigenvaluesGeneratingReflections=List([1..nbGeneratingReflections],x->1/2)
#   .OrdersGeneratingReflections=List([1..nbGeneratingReflections],x->2)
#     [The eigenvalues and orders of the generators, filled in for upward 
#      compatibility with complex reflection groups]
#
#   .orbitRepresentative=List of length .nbGeneratingReflections indicating
#     for each reflection the smallest index of a reflection in Parent(W)
#     to which it is conjugate
#
#    In addition, if W.isFinite is true,
#      W.longestElm
#      W.longestCoxeterWord
#      W.N  (number of reflections)
#
#    are computed by calling appropriate methods.

IsCoxeterGroup:=W->IsRec(W) and IsBound(W.isCoxeterGroup) and
                                W.isCoxeterGroup=true;

AbsCoxOps.CompleteCoxeterGroupRecord:=function(W)local i,j,ord;
  W.operations:=ShallowCopy(W.operations);Inherit(W.operations,AbsCoxOps);
  if not IsBound(W.nbGeneratingReflections) then
    W.nbGeneratingReflections:=Length(W.reflections);
  fi;
  if not IsBound(W.generatingReflections) then
    W.generatingReflections:=[1..W.nbGeneratingReflections];
  fi;
  W.EigenvaluesGeneratingReflections:=List(W.generatingReflections,x->1/2);
  W.isCoxeterGroup:=true;
  if not IsBound(W.rootInclusion) then
     W.rootInclusion:=W.generatingReflections;
     W.rootRestriction:=W.generatingReflections;
  fi;
  if not IsBound(W.reflectionsLabels) then
    W.reflectionsLabels:=W.rootInclusion{W.generatingReflections};
  fi;
  if not IsBound(W.operations.ReflectionFromName) then
    W.operations.ReflectionFromName:=function(W,x)
      return Position(W.reflectionsLabels,x);end;
  fi;
  if not IsBound(W.semisimpleRank) then
    W.semisimpleRank:=W.nbGeneratingReflections;
  fi;
  if not IsBound(W.rank) then
    if IsBound(W.matgens) then W.rank:=Length(W.matgens[1][1]);
    else W.rank:=W.semisimpleRank;fi;
  fi;
  if not IsBound(W.identity) and W.semisimpleRank>0 then
    W.identity:=W.reflections[1]^0;
  fi;
  if not IsBound(W.orbitRepresentative) then
    ord:=DecomposedMat(List(CoxeterMatrix(W),x->List(x,y->y mod 2)));
    W.orbitRepresentative:=[];
    for i in [1..Length(ord)] do
      W.orbitRepresentative{ord[i]}:=ord[i]*0+W.rootInclusion[Minimum(ord[i])];
    od;
  fi;
  if not IsBound(W.OrdersGeneratingReflections) then
     W.OrdersGeneratingReflections:=List(W.generatingReflections,x->2);
  fi;
  if not IsBound(W.type) then
    i:=ReflectionType(W);
    if i<>false then W.isFinite:=true;W.type:=i;
    else W.isFinite:=false;
    fi;
  fi;
  if IsFinite(W) then W.N:=Length(LongestCoxeterWord(W));fi;
end;

###########################################################################
#
# N.B. a difference with older versions of Chevie is that IsLeftDescending and
# FirstLeftDescending  return  an  index  in  the  list  of  generators, not a
# reflectionName,  for  efficiency.  LeftDescentSet  still  returns  names  of
# reflections.
#
###########################################################################

##  generic reduction of FirstLeftDescending to IsLeftDescending
##  in practical implementations of Coxeter groups this routine can often
##  be overriden by a more efficient one.

AbsCoxOps.FirstLeftDescending:=function(W, x)local i,ILD;
  ILD:=W.operations.IsLeftDescending; # avoid dispatching overhead
  for i in W.generatingReflections do if ILD(W,x,i) then return i; fi; od;
  return false;
end;

#############################################################################
##
#F  LeftDescentSet( W, x)  the set of generators s such that sx < x
##
##  the generators returned as their reflectionName
##

AbsCoxOps.LeftDescentSet:=function(W, x)
  return W.reflectionsLabels{Filtered(W.generatingReflections,
        i->IsLeftDescending(W,x,i))};
end;

#############################################################################
##
#F  RightDescentSet( W, x)  the set of generators s such that xs < x
##

RightDescentSet:=function(W,w) return LeftDescentSet(W,w^-1); end;

#############################################################################
##
#F  CoxeterWord( <W> , <w> )  . . . . . . . . for Coxeter groups
##
##  'CoxeterWord' returns a reduced word in  the standard generators of W
##  determined by the element w of W. The generators are denoted by the
##  corresponding element of W.reflectionsLabels.
##

AbsCoxOps.CoxeterWord:=function(W, w)local  l, i;
  l := [];
  while true do
    i := W.operations.FirstLeftDescending(W, w);
    if i=false then return W.reflectionsLabels{l}; fi;
    Add(l,i);
    w := W.reflections[i]*w;
  od;
end;

#############################################################################
##
#F  EltWord ( <W> , <w> ) . .  convert a word to a group element.
##
##  Returns the element of the Coxeter group <W> corresponding to the
##  Coxeter word <w>.

AbsCoxOps.EltWord:=function(W,w)
  if Length(w)=0 then return W.identity;fi;
  return Product(W.reflections{List(w,a->W.operations.ReflectionFromName(W,a))});
end;

#############################################################################
##
#F  ReducedCoxeterWord( <W> , w )  . . . . . . . . a reduced word for w in W
##
##  'ReducedCoxeterWord' returns a reduced expression for the element w, given
##  as  an  arbitrary  list  of  W.reflectionsLabels
##
ReducedCoxeterWord:=function( W, w )return CoxeterWord(W, EltWord (W, w)); end;

#############################################################################
##
#F  CoxeterLength( <W> , <w> )  . . . . . length of element w of W
##
##  'CoxeterLength'  returns  the length of w in the generating reflections.
##
AbsCoxOps.CoxeterLength:=function(W, w)return Length(CoxeterWord(W, w));end;

AbsCoxOps.Inversions:=function(W,w)local l;
  if IsList(w) then l:=List([Length(w),Length(w)-1..1],i->EltWord(W,
    Concatenation(w{[Length(w),Length(w)-1..i]},w{[i+1..Length(w)]})));
    return List(l,x->PositionProperty([1..W.N],j->Reflection(W,j)=x));
  else return Inversions(W,CoxeterWord(W,w));
  fi;
end;

#############################################################################
##
#F  CoxeterElements( <W>[, <l>]) . . the set of elements of W of length l
##
##  Returns as a list of elements of W the set of elements of length l,
##  sorted (according to the default order on elements of W).
##  <l> can be an integer or a list of integers.
##

AbsCoxOps.CoxeterElements:=function(arg)local W,l,e,i,H,x,ll;
  W:=arg[1];
  if Length(arg)=1 then
    if not IsBound(W.N) then Error("W should be finite");
    else return Concatenation(CoxeterElements(W,[0..W.N]));
    fi;
  fi;
  if IsInt(arg[2]) then ll:=[arg[2]]; else ll:=arg[2];Sort(ll);fi;
  if not IsBound(W.elts) then W.elts:=[[W.identity]];fi;
  for l in ll do
  if not IsBound(W.elts[l+1]) then
   if IsFinite(W) then
     if l>W.N then 
       InfoChevie("#W no elements of that length\n");
       return [];
     fi;
     H:=ReflectionSubgroup(W,
                        W.reflectionsLabels{[1..W.nbGeneratingReflections-1]});
     W.elts[l+1]:=[];
     if not IsBound(W.rc) then 
       W.rc:=List([0..W.N],x->[]);
       for x in ReducedRightCosetRepresentatives(W,H) do
          Add(W.rc[1+CoxeterLength(W,x)],x);
       od;
     fi;
     for i in [0..Minimum(l,H.N)] do
       e:=CoxeterElements(H,i);
       for x in W.rc[1+l-i] do Append(W.elts[l+1],e*x);od;
     od;
     W.elts[l+1]:=Set(W.elts[l+1]);
     InfoChevie2("#I Number of elements of length ",l,": ",
                                      Length(W.elts[l+1]),"\n");
     if W.N-l<>l then
       W.elts[1+W.N-l]:=Set(W.elts[l+1]*LongestCoxeterElement(W));
     fi;
   else
     W.elts[l]:=CoxeterElements(W,l-1);W.elts[l+1]:=[];
     for e in W.elts[l] do 
       for i in W.generatingReflections do
         if not IsLeftDescending(W,e,i) then 
           Add(W.elts[l+1],W.reflections[i]*e);
         fi;
       od;
     od;
     W.elts[l+1]:=Set(W.elts[l+1]);
   fi;
  fi;
  od;
  if IsInt(arg[2]) then return W.elts[arg[2]+1];
  else return W.elts{arg[2]+1};
  fi;
end;

#############################################################################
##
## The elements are returned in increasing length order...
##

AbsCoxOps.Elements:=W->AbsCoxOps.CoxeterElements(W);

#############################################################################
##
#F  CoxeterWords( <W> [, <l>] )  . . . . . .  words of W [ of length <l>]
##
##  'CoxeterWords' returns the list words in the Coxeter group W of length l.
##  If l is omitted then all words are returned, sorted by increasing length.
##

AbsCoxOps.CoxeterWords := function(arg)local W;
  W:=arg[1];
  if Length(arg)=1 then
    if not IsBound(W.N) then Error("W should be finite");
    else return Concatenation(List([0..W.N],i->CoxeterWords(W,i)));
    fi;
  else return List(CoxeterElements(W,arg[2]),w->CoxeterWord(W,w));
  fi;
end;

#############################################################################
##
#F  ReducedInRightCoset( <W> , <w> )  . . . . . reduced element in coset W.w
##
##  w is an automorphism of a parent Coxeter group of W.
##  'ReducedInRightCoset' returns the unique element in the right
##  coset W.w which is W-reduced.
##

AbsCoxOps.ReducedInRightCoset:=function(W,w)local i;
  while true  do
    i := FirstLeftDescending(W, w);
    if i=false then return w; fi;
    w := W.reflections[i] * w;
  od;
end;

#############################################################################
##
#F  ReducedRightCosetRepresentatives( <W>, <H>[, l])  . . . . . . . . . . .
#F  . . . . . . . . . . . . . . .  distinguished right coset representatives
##
##  ReducedRightCosetRepresentatives returns a list of reduced words in the
##  Coxeter  group <W>  that are  distinguished right coset representatives
##  for  the right cosets Hw where H is a reflection subgroup of W. if l is
##  given  (a length or a list of length) only representatives of length in
##  l are given. The result is sorted by increasing CoxeterLength.
##

ReducedRightCosetRepresentatives:=function(arg)local W,H,i,res,new,l,range;
  W:=arg[1];H:=arg[2];
  if Length(arg)=3 then range:=arg[3];
                        if IsInt(range) then range:=[range];fi;
  fi;
  res:=[[W.identity]];
  repeat
    new:=Concatenation(List(res[Length(res)],function(w)
      l:=W.generatingReflections;
      l:=Filtered(l,i->not IsLeftDescending(W,w^-1,i));
      l:=w*W.reflections{l};
      l:=Filtered(l,x->x=ReducedInRightCoset(H, x));
      return l;end));
    Add(res,Set(new));
  until Length(new)=0 or (IsBound(range) and Length(res)>Maximum(range));
  if IsBound(range) then 
    range:=Filtered(range,x->x<=1+Length(res));res:=res{1+range};
  fi;
  res:=Concatenation(res);
  InfoChevie2("#I nb. of cosets: ",Length(res),"\n");
  return res;
end;

#############################################################################
##
#F  PermsCosetsSubgroup ( <W>, <H> ) . . . . . . . . . . . . . .
#F  . . . .  permutation representation on the cosets of a reflection subgroup
##
##  'PermsCosetsSubgroup'  returns  the  list of permutations induced
##  by the standard generators of the Group  <W>  on the cosets  of  the
##  subgroup <H>. The cosets are in the order determined by the result of
##  the  function 'ReducedRightCosetRepresentatives( <W>, <H>)'.
##
PermCosetsSubgroup :=function(W,H)local D;
  D:=ReducedRightCosetRepresentatives(W, H);
  return List(W.reflections{W.generatingReflections},
    s->PermListList(List(D*s,x->ReducedInRightCoset(H,x)),D));
end;

#############################################################################
##
#F  Bruhat( <W>, <x>, <y> ) . . . . . . . Bruhat partial order
##
##  'Bruhat'  returns true, if the element  <x>  is less than or equal to the
##  element <y> of the Coxeter group <W>, and false otherwise. Both <x> and
##  <y>  must be given as group elements.
##
Bruhat:= function(W,x,y)local  s,i,d;
  if x=W.identity then return true;fi;
  # note: suppressing dispatch by calling W.operations.IsLeftDescending, etc..
  # gains a factor of 2 in time
  d:=W.operations.CoxeterLength(W,y)-W.operations.CoxeterLength(W,x);
  while d>0 do
    i:=W.operations.FirstLeftDescending(W,y);
    s := W.reflections[i];
    if W.operations.IsLeftDescending(W,x,i) then
      if x=s then return true;fi;
      x:=s*x; 
    else d:=d-1;
    fi;
    y:=s*y;
  od;
  return x=y;
end;

#############################################################################
##
#F  BruhatSmaller( <W>, <w>) . . List of elements Bruhat smaller than w
##
##  The result is a list whose i-th element is the list of elements of W
##  Bruhat smaller than w and of length i-1 (so the first element of the
##  list holds only W.identity and the CoxeterLength(W,w)-th element of
##  the lists holds only w).
##
BruhatSmaller:=function(W,w)local res,s,i,j;
  if w=W.identity then return [[w]];fi;
  i:=FirstLeftDescending(W,w);
  s:=W.reflections[i];
  res:=BruhatSmaller(W,s*w);Add(res,[]);
  for j in [1..Length(res)-1] do
    UniteSet(res[j+1],s*Filtered(res[j],x->not IsLeftDescending(W,x,i)));
  od;
  return res;
end;

#############################################################################
##
#F  BruhatPoset( <W>[, <w>]) . . Bruhat Poset of W [of elts. smaller than w]
##
##  The result is a poset for the required Bruhat interval with extra fields
##  .elts: elements smaller than w
##  .action: action[i] is action of multiplication on the left by 
##          W.reflections[i] on the interval
##  .W: the group W
BruhatPoset:=function(arg)local W,w,s,p,l,new,h,i,j,k;
  W:=arg[1];
  if Length(arg)=1 then w:=LongestCoxeterElement(W);else w:=arg[2];fi;
  if w=W.identity then return rec(elts:=[w],hasse:=[[]],
    action:=List(W.generatingReflections,x->[]),operations:=PosetOps,size:=1,
    label:=function(p,n,opt)return IntListToString(CoxeterWord(p.W,p.elts[n]));
           end,
    W:=W);
  fi;
  s:=FirstLeftDescending(W,w);p:=BruhatPoset(W,W.reflections[s]*w);l:=Size(p);
  new:=Filtered([1..l],k->not IsBound(p.action[s][k]));
  Append(p.elts,W.reflections[s]*p.elts{new});Append(p.hasse,List(new,x->[]));
  p.action[s]{new}:=l+[1..Length(new)];p.action[s]{l+[1..Length(new)]}:=new;
  for i in [1..Length(new)] do Add(p.hasse[new[i]],l+i);od;
  for i in [1..l] do j:=p.action[s][i];
    if j>i then
      for h in p.action[s]{p.hasse[i]} do
	if h>l then Add(p.hasse[j],h);
        k:=Position(W.reflections{W.generatingReflections},p.elts[h]/p.elts[j]);
	  if k<>false then p.action[k][j]:=h;p.action[k][h]:=j;fi;
	fi;
      od;
    fi;
  od;
  p.size:=Length(p.hasse);
  return p;
end;

#############################################################################
##
#F  ReducedExpressions( <W>, <w>) . . List of all redexp of w
##
ReducedExpressions:=function(W,w)local l;
  l:=LeftDescentSet(W,w);
  if Length(l)=0 then return [[]];fi;
  return Concatenation(List(l,x->List(ReducedExpressions(W,
    W.reflections[W.operations.ReflectionFromName(W,x)]*w),
      e->Concatenation([x],e))));
end;

# The following two methods will never finish if called for an infinite group

#############################################################################
##
#F  LongestCoxeterElement( <W> [, <I>] ) . . Longest element of W
##                  [ of the standard parabolic subgroup generated by I]
##
##  LongestCoxeterElement returns a group element
##

LongestCoxeterElement:=function(arg)local W,I,i,ILD,w;
  W:=arg[1];
  if Length(arg)=1 then
    if IsBound(W.longestElm) then return W.longestElm; 
    else I:=W.generatingReflections;
    fi;
  else I:=List(arg[2],x->W.operations.ReflectionFromName(W,x));
  fi;
  w:=W.identity;i:=1;
  ILD:=W.operations.IsLeftDescending; # supress dispatching overhead
  while i<=Length(I) do
    if ILD(W,w,I[i]) then i:=i+1;
    else w:=W.reflections[I[i]]*w;i:=1;fi;
  od;
  if Length(arg)=1 then W.longestElm:=w;fi;
  return w;
end;

#############################################################################
##
#F  LongestCoxeterWord( <W> )  . . . . . . . . . .  the longest element in W
##
## 'LongestCoxeterWord' returns a reduced expression in the  standard
##  generators of the unique longest element of the Coxeter group W.
##

LongestCoxeterWord:=function(W)local save;
  if not IsBound(W.longestCoxeterWord) then
    W.longestCoxeterWord:=CoxeterWord( W, LongestCoxeterElement(W));
  fi;
  return W.longestCoxeterWord;
end;

CoxeterMatrixFromCartanMat:=function(m)local res,i,j,find,c;
  find:=function(c)local x;
    x:=[2,3,4,6,0];if c in [0..4] then return x[c+1];fi;
    x:=NofCyc(c);
    if c=2+E(x)+E(x)^-1 then return x; 
    elif c=2+E(2*x)+E(2*x)^-1 then return 2*x;
    else Error("not a Cartan matrix of a Coxeter group");
    fi;
  end;
  res:=IdentityMat(Length(m));
  for i in [2..Length(m)] do for j in [1..i-1] do
    c:=find(m[i][j]*m[j][i]);res[i][j]:=c;res[j][i]:=c;
  od;od;
  return res;
end;

#############################################################################
##
#F CoxeterMatrix(<W>) . . . . . . . . . . .  the Coxeter  matrix <W>
##
## The result is the matrix whose entry  'm[i][j]' contains the  order of 
## $g_i\*g_j$ where  $g_i$ is the $i$-th Coxeter generator of <W>. 
## An infinite order is represented by the entry 0.
##

CoxeterMatrix:=function(W)local ord,i,j,m;
  if not IsBound(W.coxeterMat) then
    if IsBound(W.cartan) then
      W.coxeterMat:=CoxeterMatrixFromCartanMat(W.cartan);
    else
      W.coxeterMat := IdentityMat(W.nbGeneratingReflections);
      for i in W.generatingReflections do
	for j in [1..i-1] do
	  m:=W.reflections[i]*W.reflections[j];
	  if IsMat(m) and not IsFinite(W) # try to detect cheaply an infinite bond
	     and TraceMat(m)=Length(m) then ord:=0;
	  else ord:=Order(W,W.reflections[i]*W.reflections[j]);
	  fi;
	  W.coxeterMat[i][j]:=ord;
	  W.coxeterMat[j][i]:=ord;
	od;
      od;
    fi;
  fi;
  return W.coxeterMat;
end;

#############################################################################
##
#F CartanMatFromCoxeterMatrix( <m> )
##
## The  argument is a CoxeterMatrix for a  finite Coxeter group <W> and the
## result  is a Cartan Matrix for the standard reflection representation of
## <W>. Its diagonal terms are 2 and the coefficient between two generating
## reflections   s  and  t  is   -2\cos(\pi/m[s,t]),  where  by  convention
## \pi/m[s,t]=0  if  $m[s,t]=\infty$,  which  is  represented  in CHEVIE by
## setting m[s,t]:=0).

CartanMatFromCoxeterMatrix:=m->List(m,l->List(l,
  function(c)if c=0 then return -2;else return -E(2*c)^-1-E(2*c);fi;end));

AbsCoxOps.CartanMat:=function(W)
  if not IsBound(W.cartan) then
    W.cartan:=CartanMatFromCoxeterMatrix(CoxeterMatrix(W));
  fi;
  return W.cartan;
end;

AbsCoxOps.ReflectionType:=function(W)local t;
  t:=ReflectionType(CartanMat(W));
  if ForAny(t,x->x=false) then return false;fi;
  W.operations:=ShallowCopy(W.operations);
  Inherit(W.operations,HasTypeOps);
  return t;
end;

AbsCoxOps.BraidRelations:=function(W)local m,p;
  p:=function(i,j,b)
    return W.reflectionsLabels{List([1..b],k->i*(k mod 2)+j*((1-k)mod 2))};
  end;
  m:=CoxeterMatrix(W);
  return Concatenation(List([1..Length(m)],i->List([1..i-1],
    j->[p(i,j,m[i][j]),p(j,i,m[i][j])])));
end;

#############################################################################
##
#F  ReflectionSubgoup( <W>, <J> ) . . . subgroup gen. by reflections in J
##
## <J> is given as a list of reflection names for W
##

AbsCoxOps.ReflectionSubgroup:=function(W,J)local refs,P,inc,l;
  # checking if subgroup in cache:
  P:=CHEVIE.GetCached(W,"ReflectionSubgroups",rec(callarg:=[J]),x->x.callarg);
  if IsBound(P.isGroup) then return P;fi;
  if not IsSubset(W.reflectionsLabels,J) then # Use Dyer's method
    refs:=List(J,i->Reflection(W,i));
    Inherit(P,Subgroup(W,refs));
    P.reflections:=Union(Orbits(P,refs));
    P.rootInclusion:=List(P.reflections,function(x)local i;i:=0;
      repeat i:=i+1; until Reflection(W,i)=x;return i;end);
    SortParallel(P.rootInclusion,P.reflections);
    P.generatingReflections:=Filtered([1..Length(P.reflections)],
      t->Number(P.reflections,s->CoxeterLength(W,
        P.reflections[t]*s)<CoxeterLength(W,P.reflections[t]))=1);
    l:=Concatenation(P.generatingReflections,Difference(
      [1..Length(P.reflections)],P.generatingReflections));
    P.reflections:=P.reflections{l};P.rootInclusion:=P.rootInclusion{l};
    P.generatingReflections:=[1..Length(P.generatingReflections)];
    P.nbGeneratingReflections:=Length(P.generatingReflections);
    P.semisimpleRank:=Length(P.generatingReflections);
    P.reflectionsLabels:=P.rootInclusion{P.generatingReflections};
    P.operations.IsLeftDescending:=function(P,w,i)
      return CoxeterLength(P.reflectionParent,P.reflections[i]*w)<
             CoxeterLength(P.reflectionParent,w);end;
  elif Set(J)=Set(W.reflectionsLabels) then return W;
  else
    inc:=List(J,x->W.operations.ReflectionFromName(W,x));
    refs:=W.reflections{inc};
    Inherit(P,Subgroup(W,refs));
    P.reflections:=refs;
    P.semisimpleRank:=Length(P.reflections);
    P.rootInclusion:=W.rootInclusion{inc};
    P.reflectionsLabels:=W.reflectionsLabels{inc};
    P.operations.IsLeftDescending:=function(P,w,i)
      return IsLeftDescending(P.reflectionParent,w,P.rootInclusion[i]);end;
#   P.coxeterMat:=CoxeterMatrix(W){inc}{inc};
    P.cartan:=CartanMat(W){inc}{inc};
  fi;
  P.rootRestriction:=[];
  P.rootRestriction{P.rootInclusion}:=[1..Length(P.rootInclusion)];
  if IsBound(W.rank) then P.rank:=W.rank;fi;
  if IsBound(W.reflectionParent) then P.reflectionParent:=W.reflectionParent;
  else P.reflectionParent:=W;
  fi;
  AbsCoxOps.CompleteCoxeterGroupRecord(P);
  if P.semisimpleRank=0 then
       P.name:=SPrint("ReflectionSubgroup(",W,", [])");
  else P.name:=SPrint("ReflectionSubgroup(",W,", ",
			P.rootInclusion{P.generatingReflections},")");
  fi;
  return P;
end;

AbsCoxOps.BrieskornNormalForm:=function(W,w)local l,I,found,IL,i;
  l:=[];IL:=W.operations.IsLeftDescending;
  while true do
    I:=Filtered(W.generatingReflections,i->IL(W,w,i));
    if I=[] then return l;fi;
    Add(l,W.rootInclusion{I});
    repeat
      found:=false;
      for i in I do if IL(W,w,i) then found:=true;w:=W.reflections[i]*w; fi; od;
    until not found;
  od;
end;

# Calls f on each element of W
ForEachElement:=function(W,f)local l,g;
  if not IsFinite(W) then Error("only for finite Coxeter groups");
  elif W.nbGeneratingReflections=0 then f(W.identity);return;
  fi;
  l:=List([0..W.nbGeneratingReflections],i->
        ReflectionSubgroup(W,W.reflectionsLabels{[1..i]}));
  l:=List([1..Length(l)-1],i->ReducedRightCosetRepresentatives(l[i+1],l[i]));
  g:=function(x,v)local y;
    if Length(v)=0 then f(x);
    else for y in x*v[1] do g(y,v{[2..Length(v)]});od;
    fi;
  end;
  g(W.identity,l);
end;

# Calls f on each coxeter word for W
ForEachCoxeterWord:=function(W,f)local l,g;
  if not IsFinite(W) then Error("only for finite Coxeter groups");
  elif W.nbGeneratingReflections=0 then f(W.identity);return;
  fi;
  l:=List([0..W.nbGeneratingReflections],i->
        ReflectionSubgroup(W,W.reflectionsLabels{[1..i]}));
  l:=List([1..Length(l)-1],i->
    List(ReducedRightCosetRepresentatives(l[i+1],l[i]),w->CoxeterWord(W,w)));
  g:=function(x,v)local y;
    if Length(v)=0 then f(x);
    else for y in v[1] do g(Concatenation(x,y),v{[2..Length(v)]});od;
    fi;
  end;
  g([],l);
end;

# I is subset of W.rootInclusion{W.generatingReflections}
# return all W-conjugate subsets
StandardParabolicClass:=function(W,I)local res,n,new;
  res:=[]; new:=[I];
  repeat
    n:=Set(Concatenation(List(new,function(I)local J,rI;
      rI:=W.reflections{List(I,x->W.operations.ReflectionFromName(W,x))};
      J:=List(Difference(W.reflectionsLabels,I),
         function(i)local I1;I1:=Concatenation(I,[i]);
	   if IsFinite(W) or IsFinite(ReflectionSubgroup(W,I1)) 
	   then return LongestCoxeterElement(W,I1);
	   else return W.identity;fi;end);
      return Set(List(J,w->Set(List(OnTuples(rI,w),
         r->W.reflectionsLabels[Position(W.reflections,r)]))));end)));
    UniteSet(res,new); new:=Difference(n,res);
  until Length(new)=0;
  return res;
end;

AbsCoxOps.ParabolicRepresentatives:=function(W,s)local orbits,o,l;
  l:=Combinations(W.reflectionsLabels,s);
  orbits:=[];
  while Length(l)>0 do
    o:=StandardParabolicClass(W,l[1]);
    Add(orbits,o);
    l:=Difference(l,o);
  od;
  return List(orbits,x->x[1]);
end;

##############################################################################
##
## CoxeterHeckeAlgebraOps: operations for Hecke algebras for Coxeter groups.
##
CoxeterHeckeAlgebraOps:=OperationsRecord("CoxeterHeckeAlgebraOps",AbsHeckeOps);

AbsCoxOps.Hecke:=function(arg)local H;
  H:=ApplyFunc(PermRootOps.Hecke,arg);
  H.operations:=CoxeterHeckeAlgebraOps;
  return H;
end;

ReadChv("prg/kl"); # this will actually add new methods to basis T ...
ReadChv("prg/heckemod");
