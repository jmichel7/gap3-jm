#############################################################################
##
#A  spets.g                  CHEVIE library        Jean Michel
##
#Y  Copyright (C) 1992 - 2006  Lehrstuhl D fur Mathematik, RWTH Aachen
#Y  and   University  Paris VII.
##
##  This file contains general functions that deal with spets. The code
##  has been mostly yanked from coset.g
##
##  A Spets W\phi contains the following fields related to its torsion:
##  .phi     A standard representative of W\phi (which is always the same
##           for all Spets defining the same coset)
##  .F0Mat   The torsion as a matrix (may give torsion on a part of the
##           space not spanned by the roots)

#############################################################################
##
#F  SpetsOps . . . . . . .  operations record for Spets
##  
##  We  first copy  the  basic  functions   included in 'DomainOps'    and
##  overwrite some of them by more efficient ones.
##  
#if not IsBound(CHEVIE.PrintSpets) then CHEVIE.PrintSpets:=rec(GAP:=true);fi;
if not IsBound(CHEVIE.PrintSpets) then CHEVIE.PrintSpets:=rec();fi;

SpetsOps:=OperationsRecord("SpetsOps");

Inherit(SpetsOps,PermGroupOps,["Parent","ClassInvariants","PositionClass"]);

SpetsOps.TorusOrder:=PermRootOps.TorusOrder;

SpetsOps.IsFinite:=function(WF)return true;end;

SpetsOps.Group:= x->x.reflectionGroup;

SpetsOps.Random:=function(WF)return Random(Group(WF))*WF.phi;end;

SpetsOps.Rank:= function(WF)return Rank(Group(WF));end;

SpetsOps.Representative:=function(WF)return WF.phi;end;

SpetsOps.SemisimpleRank:= function(WF)return SemisimpleRank(Group(WF));end;

SpetsOps.Size:= function(WF)return Size(Group(WF));end;

SpetsOps.Elements:=function(WF)
  if WF.phi in Group(WF) then return Elements(Group(WF));
  else return Set(Elements(Group(WF))*WF.phi);
  fi;
end;

SpetsOps.\^:=function(WF,g)local W;
  W:=Group(WF);
  return Spets(ReflectionSubgroup(Parent(W),
    OnTuples(W.rootInclusion{W.generatingReflections},g)),WF.phi^g);
end;

SpetsOps.\in:=function(x,WF)return x/WF.phi in Group(WF);end;

SpetsOps.EltWord:=function(W,w)return EltWord(Group(W),w)*W.phi; end;

SpetsOps.\=:=function(WF1,WF2)
  return IsRec(WF1) and IsRec(WF2) and
         (WF1.F0Mat=WF2.F0Mat and Group(WF1)=Group(WF2));
end;

SpetsOps.Format:=function(WF,option)local res;
  if IsBound(option.GAP) then 
    res:="";
    if IsBound(WF.parent) and WF<>WF.parent then
      if IsCoxeterCoset(WF) then PrintToString(res,"CoxeterSubCoset(");
                            else PrintToString(res,"SubSpets(");fi;
      PrintToString(res,WF.parent,", ",
	     Group(WF).rootInclusion{Group(WF).generatingReflections});
      if WF.phi/WF.parent.phi<>() then
	PrintToString(res,", ",WF.phi/WF.parent.phi);
      fi;
    else
      if IsCoxeterCoset(WF) then PrintToString(res,"CoxeterCoset(");
                            else PrintToString(res,"Spets(");fi;
      PrintToString(res,Group(WF));
      if WF.F0Mat<>WF.F0Mat^0 then PrintToString(res,", ");
        if WF.F0Mat=MatXPerm(Group(WF),WF.phi) then PrintToString(res,WF.phi);
	else PrintToString(res,WF.F0Mat);
	fi;
      fi;
    fi;
    Append(res,")");return res;
  else return ReflectionName(WF,option);
  fi;
end;

#############################################################################
#F  SpetsOps.String( <WF> ) . . . . . . . . . . 'String' function for Spets
SpetsOps.String:=function(WF)return Format(WF,CHEVIE.PrintSpets);end;

SpetsOps.Display:=function(WF,opt)Print(Format(WF,opt),"\n");end;

SpetsOps.Print:=function(WF)Print(Format(WF,CHEVIE.PrintSpets));end;

SpetsOps.MatXPerm:=function(WF,w)
  if WF.F0Mat=[] then return [];fi; #case A0
  return WF.F0Mat*MatXPerm(Group(WF),WF.phi^-1*w);
end;

#############################################################################
##
#F  SpetsOps.FusionConjugacyClasses( <W1>, <W2> ) . . . . . . . . . . 
#F  'FusionConjugacyClasses' for Spets
##
SpetsOps.FusionConjugacyClasses:=function(u,g)
  if not IsSpets(g) then Error(u," is a coset but ",g," is not");fi;
  if Parent(u)=Parent(g) then
    return List(List(ConjugacyClasses(u),Representative),
                x->PositionClass(g,x));
  else
    Error("Not same parent: don't know how to compute fusion.\n");
  fi;
end;

SpetsOps.ParabolicRepresentatives:=function(WF,s)local W,l;
  W:=Group(WF);
  l:=ParabolicRepresentatives(W,s);
  if WF.phi=() then return l;fi;
  l:=List(l,function(I)local p,o;I:=Spets(ReflectionSubgroup(W,I));
    p:=RepresentativeOperation(W,I,I^WF.phi);
    if p=false then return false;fi;
    o:=Orbit(W,I);
    p:=PositionProperty(o,x->x^WF.phi=x);
    if p=false then return p;fi;
    return List(Group(o[p]).generators,x->Position(Reflections(W),x));
    end);
  return Filtered(l,x->x<>false);
end;

SpetsOps.ReflectionType:=function(WF)
##  The code  below computes the reflection  type of the Spets.  It is a
##  list of  records, one for  each orbit  of WF.phi on  the irreducible
##  components of W. The record has three components. If we denote by
##  phi the effect of WF.phi on [1..Length(W.roots)] (obtained via
##  W.rootInclusion)
##  
##   .orbit  The ReflectionType of the orbit (see HasTypeOps.ReflectionType)
##   .twist  RestrictedPerm(phi^Length(orbit),.orbit[1].indices)
##   .scalar An "Ennola-scalar" not accounted for by the twist
##  
##  The indices in .orbit are given such that for all i
##      .orbit([i+1].indices)=OnTuples(.orbit[i].indices,phi)
##  
##  If the  components in .orbit  are of type D_4  and phi^Length(orbit)
##  has order  2 on  .orbit[1].indices then .orbit[1].indices  is sorted
##  such that phi^Length(orbit) permutes the first two entries.
##  
##  If the  components in .orbit  are of type D_4  and phi^Length(orbit)
##  has order  3 on  .orbit[1].indices then .orbit[1].indices  is sorted
##  such that phi^Length(orbit) permutes the entries 1->2->4->1.
##  
##  This function may modify WF.phi and WF.F0Mat so should be called before
##  computing Spets data which depend on them (it is called in Spets and
##  SubSpets).
  local W,a,type,c,l,i,j,t,scal,z,zg,v,m,next,PW,w0,scals,fixG333;
  W:=WF.reflectionGroup;PW:=Parent(W);
  type:=List(ReflectionType(W),ShallowCopy);
  fixG333:=function()local G333,a,c;
    for a in type do
      if a.series="ST" and IsBound(a.p) and [a.p,a.q,a.rank]=[3,3,3] then
	c:=CHEVIE.R("ReducedInRightCoset","timp")(a.subgroup,WF.phi);
	if c<>false and WF.phi<>c.phi or c.gen<>a.subgroup.rootInclusion{c.gen}
        then 
	  c.gen:=a.subgroup.rootInclusion{c.gen};
	  a.indices:=W.rootRestriction{c.gen};
	  a.subgroup:=ReflectionSubgroup(W,c.gen); 
	  WF.phi:=c.phi;
	  WF.reflectionGroup:=ReflectionSubgroup(W,
	    Concatenation(List(type,x->W.rootInclusion{x.indices})));
	  W:=WF.reflectionGroup;
	fi;
	if Comm(a.subgroup.3,WF.phi)=() 
	then G333:=[1,2,3,44];
	else G333:=[1,50,3,12];
	fi;
	a.subgens:=List(G333,i->Reflection(a.subgroup,i));
	a.indices:=W.rootRestriction{a.subgroup.rootInclusion{G333}};
      fi;
    od;
  end;
  if WF.phi=() then
    WF.type:=List(type,x->rec(orbit:=[x],operations:=ReflTypeOps,twist:=()));
  else

  for a in type do
    a.subgroup:=ReflectionSubgroup(W,W.rootInclusion{a.indices});
    a.subgens:=ShallowCopy(a.subgroup.generators);
  od;
  if WF.phi<>() and ForAny(type,a->
      a.series="ST" and IsBound(a.p) and [a.p,a.q,a.rank]=[3,3,3]) then
      fixG333();
  fi;
  c:=List(type,a->Set(a.subgens));
  c:=PermListList(List(c,x->OnSets(x,WF.phi)),c);
  if c=false then 
    WF.operations.Print:=function(x)Print("***********");end;
    Error("fixG333() failed");
  fi;
# if c=false then fixG333();
#   c:=List(type,a->Set(a.subgens));
#   c:=PermListList(List(c,x->OnSets(x,WF.phi)),c);
#   if c=false then
#     ChevieErr("WF.phi does not preserves generators");
#     return false;
#   fi;
# fi;
  c:=Cycles(c^-1,[1..Length(type)]);

  WF.type:=List(c,x->rec(orbit:=type{x},operations:=ReflTypeOps));
  scals:=function(roots,images,phi)
    roots:=PW.roots{OnTuples(W.rootInclusion{roots},phi)};
    images:=PW.roots{W.rootInclusion{images}};
    return List([1..Length(roots)],
      i->ProportionalityCoefficient(roots[i],images[i]));
  end;
  for c in WF.type do
    l:=Length(c.orbit);
    c.scalar:=[];
    for i in [1..l] do
      a:=c.orbit[i];
      if i=l then next:=c.orbit[1]; else next:=c.orbit[i+1]; fi;
      t:=PermListList(next.subgens,OnTuples(a.subgens,WF.phi));
      if i<>l then
	next.indices:=Permuted(next.indices,t);
	next.subgroup:=ReflectionSubgroup(W,W.rootInclusion{next.indices});
        next.subgens:=ShallowCopy(next.subgroup.generators);
	scal:=scals(a.indices,next.indices,WF.phi);
      else 
	c.twist:=t;
	scal:=scals(a.indices,Permuted(next.indices,t^-1),WF.phi);
      fi;
      if Length(Set(scal))>1 or ForAny(scal,x->x=false) then 
        if a.series="B" and a.rank=2 and not false in scal then
	  scal:=[GetRoot(Product(scal),2)];
	else ChevieErr("ReflectionType failed: no element of ",
	  ReflectionName(W),".",FormatGAP(WF.F0Mat),
          " acts as scalar on orbits: scalars=",scal," twist=",t,"\n");
        #ReflTypeOps.Print:=function(x)Print("");end;
        return false;
	fi;
      fi;
      scal:=scal[1];
# now simplify scalar as much as possible using centre
      zg:=Centre(a.subgroup);z:=Size(zg);
      if z>1 then
	zg:=First(Elements(zg),x->OrderPerm(x)=z);
	v:=AsRootOfUnity(scals(a.indices{[1]},a.indices{[1]},zg)[1]);
	zg:=zg^PowerMod(Numerator(v),-1,Denominator(v));# distinguished generator
	v:=AsRootOfUnity(scal)+[0..z-1]/z;
	m:=Minimum(List(v,Denominator));
	m:=PositionProperty(v,n->Denominator(n)=m);
	scal:=E(Denominator(v[m]))^Numerator(v[m]);
	WF.phi:=zg^(m-1)*WF.phi;
      else
	v:=[AsRootOfUnity(scal)];m:=1;
      fi;
# simplify again by -1 in types 2A(>1), 2D(odd), 2E6
      if Denominator(v[m]) mod 2 <> Denominator(2*v[m]) mod 2 and 
        (a.series in ["A","D"] or (a.series="E" and a.rank=6))then 
	w0:=Product(a.subgroup.generators{LongestCoxeterWord(CoxeterGroup(a))});
        WF.phi:=w0*WF.phi;
        t:=PermListList(next.subgens,OnTuples(a.subgens,WF.phi));
#Print("l=",l," i=",i," scal before:",scal);
	if i<>l then
	  next.indices:=Permuted(next.indices,t);
	  next.subgroup:=ReflectionSubgroup(W,W.rootInclusion{next.indices});
          next.subgens:=ShallowCopy(next.subgroup.generators);
	  scal:=scals(a.indices,next.indices,WF.phi);
	else c.twist:=t;
	  scal:=scals(a.indices,Permuted(next.indices,t^-1),WF.phi);
	fi;
#Print(" scal after:",scal,"\n");
	if Length(Set(scal))>1 then Error();
	else scal:=scal[1];
	fi;
      fi;
      Add(c.scalar,scal);
    od;
  od;

  # some adjustment such that in type ^2D_4 the first two simple
  # roots are permuted by res.twist, and in type ^3D_4 the permutation 
  # of the simple roots is 1 -> 2 -> 4 -> 1:
  # in type 4G333 and 3G333 restore the indices to be in 1,2,3
  WF.type:=List(WF.type,function(a) local b, rf, j,i;    
    i:=a.orbit[1].indices;b:=a.orbit[1];
    if b.series="D" and b.rank=4 then
      if OrderPerm(a.twist)=2 then # ^2D_4
	rf:=Filtered([1..4],x->Reflection(W,i[x])<>Reflection(W,i[x])^WF.phi);
	rf:=Concatenation(rf,[3],Difference([1,2,4],rf));
	for j in a.orbit do j.indices:=j.indices{rf};od;
      elif OrderPerm(a.twist)=3 then # triality group ^3D_4
	if Reflection(W,i[1])^WF.phi<>Reflection(W,i[2]) then 
	  for j in a.orbit do j.indices:=j.indices{[1,4,3,2]};od;
	fi;
      fi;
    elif b.series="ST" and IsBound(b.p) and b.p=3 and b.q=3 and b.rank=3 then
      if OrderPerm(a.twist)=4 then
	for j in a.orbit do 
           j.indices[2]:=W.rootRestriction[j.subgroup.rootInclusion[2]];od;
      fi;
    fi;
    return a;
  end);
  fi;
  if IsBound(WF.scalars) then # add extra scalars
    for a in WF.type do
      if Length(Set(WF.scalars))=1 then scal:=Set(WF.scalars);
      else
      scal:=WF.scalars{a.orbit[1].indices};
      if Length(Set(scal))>1 then 
        if a.orbit[1].series="B" and a.orbit[1].rank=2 then
	  scal:=[GetRoot(Product(scal),2)];
	else
	ChevieErr("ReflectionType failed: no element of ",
	  ReflectionName(W),".",FormatGAP(WF.F0Mat),
          " acts as scalar on orbits: scalars=",scal," twist=",t,"\n");
        return false;
	fi;
      fi;
      fi;
      if IsBound(a.scalar) then a.scalar:=a.scalar*scal[1];
      else a.scalar:=[scal[1]];
      fi;
    od;
  fi;
  WF.operations:=ShallowCopy(WF.operations);
  Inherit(WF.operations,HasTypeOps);
  return WF.type;
end;

#############################################################################
##
#F  Spets( <W>[, <F0Mat> or  <perm>] ) . . . .  create a Spets
#F  Spets( <rec> ) . . . . . . . . . return component <rec>.spets
##  
##  In the first form <W> must be a complex reflection group.
##  <F0Mat> must be an invertible square  matrix of rank <W>.rank such that:
##    - <F0Mat> has finite order
##    - <F0Mat>: X->X lets invariant the set of roots of <W> and of Parent(<W>)
##  
##  If <W>.semisimpleRank  = <W>.rank it   is allowed to give an  argument
##  <perm> (which  is a permutation of at least IndependentRoots(W) --- 
##  perhaps of [1..W.nbGeneratingReflections]) instead  of <F0Mat>. Then
##  <F0Mat> is  computed as the  unique matrix which  maps the
##  independent roots of <W> on the roots given by i^<perm>. 
##  
##  If only  <W> is  given then the  default for  <F0Mat> is  the identity
##  matrix.
##  
##  'Spets' returns a record with the following components:
##    .isDomain     set to true
##    .isFinite     set to true
##    .isCoset      set to true
##    .reflectionGroup        <W> from the argument
##    .F0Mat        as described above
##    .phi          'canonical' element in the coset W*perm
##    .operations   SpetsOps
##    .callArgs     information used by 'Print'
##  
##  The form Spets(<rec>) returns the component <rec>.spets if present.
##   
CHEVIE.Cache.Spetss:=true;
Spets:=function(arg)local W, WF,perm,t,s,i,j;
  # just return  component .spets:
  if Length(arg)=1 then
    if IsRec(arg[1]) and IsBound(arg[1].spets) then return arg[1].spets;
    elif arg[1]="3G422" then return Spets(PermRootGroup(
[[2,(-1+ER(3))*E(3)], [2,(-1+ER(3))*E(3)^2], [2,(-1+ER(3))]],
[[(3+ER(3))/2,ER(3)*E(3)^2],[(3+ER(3))/2,ER(3)*E(3)],[(3+ER(3))/2,ER(3)]]/3),
(1,2,3));
    elif arg[1]="2G5" then return Spets(ComplexReflectionGroup(5),(1,2));
    elif arg[1]="3G333" then return Spets(ComplexReflectionGroup(3,3,3),(1,2,44));
    elif arg[1]="3pG333" then return Spets(ComplexReflectionGroup(3,3,3),(1,44,2));
    elif arg[1]="4G333" then return
   Spets(ComplexReflectionGroup(3,3,3),(2,12,22,38)(3,35,14,32)(5,26,8,33)
 (6,39,13,20)(7,42,24,27)(9,51,29,54)(10,15,21,46)(16,50,30,44)(17,53,23,36)
 (18,47,28,49)(19,41,40,34)(31,45,37,43));
    fi;
  fi;
  W:=arg[1];
  if IsCoxeterGroup(W) then return ApplyFunc(CoxeterCoset,arg);fi;
  arg:=arg{[2..Length(arg)]};
  WF:=rec(isDomain:=true, isFinite:=true, isCoset:=true,
   reflectionGroup:=W);

  if Length(arg)=0 then WF.F0Mat:=IdentityMat(W.rank);perm:=();WF.phi:=();
  elif IsPerm(arg[1]) then 
    perm:=arg[1];WF.F0Mat:=MatXPerm(Parent(W),perm);
    perm:=PermMatX(Parent(W),WF.F0Mat);
  else WF.F0Mat:=arg[1];
    if WF.F0Mat=[] then perm:=();
    else perm:=PermMatX(W,WF.F0Mat); # checking if roots are normalized
    fi;
  fi;

  if perm=false then
# check if there exists a permutation perm and for each W-orbit of roots O 
# a scalar l_O such that W.roots{O}*WF.F0Mat=l_O*W.roots{OnTuples(O,perm)}
    t:=CollectBy(W.generatingReflections,j->W.orbitRepresentative[j]);
    s:=List(t,function(inds)local scal,l;
      scal:=List(inds,function(y)local r;
        r:=W.roots[y]*WF.F0Mat;
        r:=List(W.roots,x->ProportionalityCoefficient(r,x));
        l:=Filtered([1..Length(W.roots)],i->r[i]<>false);
        return [l,r{l}];end);
      if Length(Set(List(scal,x->Set(x[2]))))>1 then Error("theory");fi;
      l:=Intersection(List(scal,x->x[2]));
      if Length(l)=0 then return false;fi;
      l:=List(l,AsRootOfUnity);
      if false in l then 
        Error("only reflection cosets of finite order implemented");
        return false;
      fi;
      # choose simplest scal
      l:=Minimum(List(l,x->[Denominator(x),Numerator(x)]));
      l:=E(l[1])^l[2];
      return [l,List(scal,x->x[1][Position(x[2],l)])];
    end);
    if false in s then
      ChevieErr("Spets(",ReflectionName(W),",F=",FormatGAP(WF.F0Mat),
      " must normalize set of roots of parent up to scalars.\n");
      return false;
    fi;
    WF.scalars:=List(s,x->x[1]){W.orbitRepresentative{W.generatingReflections}};
    perm:=[];
    for i in [1..Length(t)] do perm{t[i]}:=s[i][2];od;
    repeat
      i:=Filtered([1..Length(perm)],j->IsBound(perm[j]));
      for j in W.generatingReflections do
        perm{OnTuples(i,W.reflections[j])}:=
          OnTuples(perm{i},Reflection(W,perm[j]));
      od;
    until Length(i)=Length(W.roots);
    if Length(Set(perm))<Length(perm) then return false;fi;
    perm:=PermList(perm);
  fi;

  if not IsBound(WF.phi) then
    if not perm in W then WF.phi:=ReducedInRightCoset(W,perm);
    else WF.phi:=();fi;
    if IsRec(WF.phi) then # for type G3,3,3 change of generators may be needed
	WF.reflectionGroup:=WF.phi.reflectionGroup;
	WF.phi:=WF.phi.phi;
    fi;
    if WF.phi=false then 
      ChevieErr("Spets(",ReflectionName(W),",F of order ",
	OrderPerm(perm),"): no system of generators is stabilized\n");
      return false;
    fi;
    if WF.phi<>perm then WF.F0Mat:=MatXPerm(Parent(W),WF.phi/perm)*WF.F0Mat;fi;
  fi;

  WF.operations:=SpetsOps;WF.callArgs:=[];
  perm:=WF.phi;
  if false=ReflectionType(WF) then return false;fi; 
  # needed call --- may change WF.phi
  if WF.phi<>perm then WF.F0Mat:=MatXPerm(Parent(W),WF.phi/perm)*WF.F0Mat;fi;
  
  # for printing:
  if Length(arg)=0 or WF.F0Mat=IdentityMat(W.rank) then WF.callArgs:=[];
  else WF.callArgs:=WF.F0Mat;
  fi;
  # check for cached cosets. Especially useful for the trivial coset
  return CHEVIE.GetCached(W,"Spetss",WF,x->x.callArgs);
end;

ReflectionCoset:=Spets;

#############################################################################
##
#F  SubSpets(  <WF>, <I>[, <w> ] ) . . . . . . sub-Spets of a Spets
##   
##  This has basically the same effect as
##      Spets(ReflectionSubgroup(Group(WF),I), w*WF.phi)
##  except that WF is recorded as a parent (when <WF> itself has a parent,
##  this in turn is recorded as a parent).
##  The default for <w> is ().
##  
SubSpets:=function(arg)local WF, W, w, tmp, I, res, i,s;
  if not IsSpets(arg[1]) then Error("first argument must be a Spets\n");fi;
  WF:=Parent(arg[1]);w:=arg[1].phi/WF.phi;W:=Group(WF);
  I:=arg[2];
  if Length(arg)>2 then
    if arg[3] in W then w:=arg[3]*w;
    else Error("must give w in Group(WF).\n");
    fi;
  fi;
  res:=ReflectionSubgroup(W,I); # includes ReflectionType may change generators
  res:=Spets(res,w*WF.phi);
  if not IsRec(res) then return res;fi;
# Print("SubSpets<",arg,"> parent=",WF," W=",W,"\n");
  res.parent:=WF;
  return CHEVIE.GetCached(Group(res),"Cosets",res,x->[x.reflectionGroup,x.F0Mat]);
end;

ReflectionSubCoset:=SubSpets;
SpetsOps.ReflectionSubgroup:=SubSpets;

#############################################################################
##
#F  ReflectionEigenvalues(W [,c]) ..... Eigenvalues [of cth] class in
#F        reflection representation
##
## The eigenvalue E(n)^i is represented by i/n
##
SpetsOps.ReflectionEigenvalues:=function(arg)local W;
  W:=arg[1];
  if not IsBound(W.eigenvalues) then
    W.eigenvalues:=List(ConjugacyClasses(W),function(c)local p;
      p:=MatXPerm(W,Representative(c));
      if p=[] then return [];fi;
      p:=CycPol(CharacteristicPolynomial(p));
      return Concatenation(List(p.vcyc,x->List([1..x[2]],i->x[1])));end);
  fi;
  if Length(arg)>1 then return W.eigenvalues[arg[2]];
  else return W.eigenvalues;
  fi;
end;

##   Intended only for unclassified Spets (overriden by HasTypeOps)
SpetsOps.ConjugacyClasses:=function(WF)local W,g,c;
  W:=Group(WF);
  g:=ApplyFunc(Group,Concatenation(W.generators,[WF.phi]));
  c:=Filtered(ConjugacyClasses(g),x->x.representative*WF.phi^-1 in W);
  for g in c do g.group:=W;od;
  return c;
end;

# the constant by which Phi acts on the discriminant (see Spets 1.5)
PhiOnDiscriminant:=function(WF)
  return Product(ReflectionType(WF),function(t)
   if IsBound(t.scalar) then return Product(t.scalar)^
     Sum(ReflectionDegrees(t.orbit[1])+ReflectionCoDegrees(t.orbit[1]));
   else return 1;fi;end);
end;

# warning! The code has been changed  since the sign is not always equal
# to  (-1)^rank(WF)*DeterminantMat(WF.F0Mat).  This last  formula  works
# only for tori  and spetses containing 1-regular  elements; The correct
# formula  is (-1)^rank(WF)*Product(PhiFactors(WF))*xi  where xi  is the
# constant by which phi acts on the discriminant (if phi is zeta-regular
# then xi=\zeta^(N+Nhyp)) (see Spets 1.5)
GenericSign:=WF->(-1)^Group(WF).rank*Product(ReflectionDegrees(WF),x->x[2])*
  PhiOnDiscriminant(WF)^-1;

## If [d1,e1],..,[dn,en] are the generalized degrees the generic
## order  of the  corresponding Spets  is q^N  theta product(q^di-ei^-1)
## where theta is the constant by which phi acts on the discriminant (if
## phi is zeta-regular then theta=\zeta^(N+Nhyp))
SpetsOps.GenericOrder:=function(WF,q)
  if Group(WF).rank=0 then return q^0;
  else return GenericSign(WF)*q^Sum(ReflectionCoDegrees(WF),x->x[1]+1)*
    Product(ReflectionDegrees(WF),p->1-q^p[1]*p[2]^-1);
  fi;
end;

SpetsOps.FakeDegrees:=function(WF,q)
  return PermRootOps.FakeDegrees(WF,q)*(-1)^Group(WF).rank/GenericSign(WF);
end;

#############################################################################
##
#F  SpetsOps.CharTable( <WF> ) . .  character table for Spets
##  
##  The character  table   of a coset   W*phi is  defined as follows:  The
##  characters  of <W, phi>  with non  zero values on  the classes  in the
##  coset W*phi   are   precisely  those   whose restriction   to   W   is
##  irreducible. To each  character of W  which is extendable to  <W, phi>
##  there are |phi| extensions.
##  
##  Following Lusztig for rational Spets we choose one of these extensions in
##  each case, called the  *preferred*  extension.  This is defined  in
##  [Lusztig, Character sheaves IV, Chapter 17.2].
##  
##  The character table is  the (square) table which  gives the  values of
##  the *preferred* extensions on the  conjugacy classes contained in  the
##  coset W*F.
##  
##  It is  a GAP  character  table record, but if   <WF>.FOMat is  not the
##  identity matrix then there are no  components .powermap and .orders.
##  
##  Superseded by HasTypeOps.CharTable but can recompute data
#  **** To be written now for general spets ****

SpetsOps.CharTable:=function(WF)
  local W,WFt,cl,ct,good,index,i,p,pos,unique,f,ct1;
  InfoChevie("# using SpetsOps.CharTable \n");
  W:=Group(WF);
  WFt:=Group(Concatenation(W.generators,[WF.phi]),());
  f:=List(ConjugacyClasses(W),x->PositionClass(WFt,Representative(x)));#fusion
  ct1:=CharTable(W);
  cl:=List(ConjugacyClasses(WF),Representative);
  good:=List(cl,x->PositionClass(WFt,x));
  ct:=CharTable(WFt);
  ct.classes:=ct.classes{good};
  index:=ct.size/Size(WF);
  ct.size:=Size(WF);
  ct.order:=Size(WF);
  ct.centralizers:=ct.centralizers{good}/index;
  ct.orders:=ct.orders{good};
  ct.irreducibles:=Filtered(ct.irreducibles,x->ForAny(x{good},y->y<>0));
  ct.irredinfo:=List(ct.irreducibles,x->ct1.irredinfo[Position(ct1.irreducibles,
    x{f})]);
  ct.irreducibles:=List(ct.irreducibles,x->x{good});
  ct.group:=WF;
  Unbind(ct.powermap);
  unique:=[];
  pos:=[];
  for i in [1..Length(ct.irreducibles)] do
    p:=PositionProperty(unique,j->ProportionalityCoefficient(
      ct.irreducibles[j],ct.irreducibles[i])<>false);
    if p=false then Add(unique,i);pos[Length(unique)]:=[i];
    else Add(pos[p],i);
    fi;
  od;
  p:=List(pos,x->List(x,y->Lcm(List(ct.irreducibles[y],NofCyc))));
  pos:=List([1..Length(pos)],i->pos[i][Position(p[i],Minimum(p[i]))]);
  ct.irreducibles:=ct.irreducibles{pos};
  ct.irredinfo:=ct.irredinfo{pos};
  ct.operations.StringEntry := function(x)
    if x=0*x then return ".";else return Format(x);fi;end;
  if ChevieClassInfo(WF)<>false then
    ct.classnames:=ChevieClassInfo(WF).classnames;fi;
  return ct;
end;

#############################################################################
##
##  Here we use the   CoxeterClassParam function to  distinguish classes
##  which are not distinguished by cycle  type (and cycle type of elements
##  multiplied by center elements and not lying in small classes).
##  
SpetsOps.ClassInvariants:=function(arg)local WF, W, WxF;
  if IsList(arg[1]) then arg:=arg[1]; fi;
  WF:=arg[1]; W:=Group(WF);
  if WF.phi=() then return ClassInvariants(W);fi;
  WxF:=Group(Concatenation(W.generators,[WF.phi]),());
  # now cheat
  WxF.conjugacyClasses:=ConjugacyClasses(WF);
  WxF.size:=Size(WF);
  return PermGroupOps.ClassInvariants(WxF,ShortClassListFunction);
end;

SpetsOps.RelativeCoset:=function(arg)local WF,W,R,res,t;
# Print("SpetsOpsRelativeCoset ",arg," called \n");
  WF:=arg[1];W:=Group(WF);
  R:=ApplyFunc(RelativeGroup,Concatenation([W],arg{[2..Length(arg)]}));
  res:=Spets(R,GetRelativeAction(W,ReflectionSubgroup(W,arg[2]),WF.phi));
  if arg[2]=[] then 
     Group(res).MappingFromNormalizer:=R.MappingFromNormalizer;
     return res;
# else Group(res).MappingFromNormalizer:=function(x)Error("MappingFromNormalizer failed");return false;end;
  fi;
  return res;
end;

#############################################################################
##
#F  Frobenius( <W> )  . . . . . . . .  returns function which gives the action
#F  of W.phi on a permutation or on the set of roots.
##  
##  
##  The function f returned by 'Frobenius( <W> )' can be applied to:
##    a permutation x:           then f(x) = x^(<W>.phi^-1)
##    a list l of integers:      then f(l) = List(l, i->i^(<W>.phi^-1))
##    a record r with component   
##        .operations.Frobenius: then f(r) = r.operations.Frobenius(<W>,r)
##  
##  If the simple roots are   invariant under <W>.phi then the   second
##  form gives  the action of the Frobenius  on elements of the Coxeter group
##  written as words in the generators.
##  
SpetsOps.Frobenius:=function(W)local finv; finv:=W.phi^-1;
  return function(arg)local i,x;x:=arg[1];
    if Length(arg)=1 then i:=1;else i:=arg[2];fi;
    if IsPerm(x) or IsInt(x) then return x^(finv^i);
    elif IsList(x) then return OnTuples(x,finv^i);
    elif IsRec(x) then
      if IsBound(x.operations) and  IsBound(x.operations.Frobenius) then 
	return x.operations.Frobenius(W,x,i);
      else Error("no method Frobenius for ",arg[1],".\n");
      fi;
   fi;
  end;
end;

IsSpets:=W->IsRec(W) and IsBound(W.isCoset) and W.isCoset;
