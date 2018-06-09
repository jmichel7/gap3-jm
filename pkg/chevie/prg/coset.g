#############################################################################
##
#A  coset.g                  CHEVIE library        Frank Luebeck, Jean Michel
##
#Y  Copyright (C) 1992 - 2016  Lehrstuhl D fur Mathematik, RWTH Aachen, IWR
#Y  der Universitat Heidelberg, and   University Paris VII.
##
##  This file contains  functions for Coxeter cosets.
##
#############################################################################
##
#F  CoxeterCosetOps . . . . . . .  operations record for Coxeter group cosets
##  
##  We  first copy  the  basic  functions   included in 'DomainOps'    and
##  overwrite some of them by more efficient ones.
##  
CoxeterCosetOps:=OperationsRecord("CoxeterCosetOps",SpetsOps);
# we can always classify them
Inherit(CoxeterCosetOps,HasTypeOps); 

CoxeterCosetOps.CoxeterGroup:=CoxeterCosetOps.Group;
# for compatibility with previous versions of CHEVIE

IsCoxeterCoset:=WF->IsRec(WF) and IsBound(WF.isCoxeterCoset)
  and WF.isCoxeterCoset;

#############################################################################
##
#F CoxeterCosetOps.ReflectionType(  <WF> ) .  . .  . . .  the reflection
#F                                                 type of a CoxeterCoset
##  
##  Gives  a list which describes  WF.  The list  has one element for each
##  orbit of <WF>.phi on the irreducible components of  W. This element is
##  a record with two components:
##  
##    .orbit  The ReflectionType of the orbit (see HasTypeOps.ReflectionType)
##    .twist  A permutation  of  .orbit[1].indices  which describes  the
##            effect of <WF>.phi^Length(orbit) on it.
##  
##  The indices in .orbit are given such that
##      .orbit([2].indices)=OnTuples(.orbit[1].indices,<WF>.phi), etc...
##  
##  If   the    components   in   .orbit    are   of   type    D_4   and
##  <WF>.phi^Length(orbit)  has   order  2  on   .orbit[1].indices  then
##  .orbit[1].indices   is  sorted   such  that   <WF>.phi^Length(orbit)
##  permutes the first two entries.
##  
##  If   the    components   in   .orbit    are   of   type    D_4   and
##  <WF>.phi^Length(orbit)  has   order  3  on   .orbit[1].indices  then
##  .orbit[1].indices   is  sorted   such  that   <WF>.phi^Length(orbit)
##  permutes the entries 1 -> 2 -> 4 -> 1.
##  
CoxeterCosetOps.ReflectionType:=function(WF)local W, c, type, phifus;
  W:=Group(WF);
  type:=ReflectionType(W);
  c:=List(type,a->Set(W.rootInclusion{a.indices}));
  type:=List(Cycles(Permutation(WF.phi,c,OnSets),[1..Length(type)]),x->type{x});
  phifus:=WF.phi^MappingPermListList(W.rootInclusion,[1..Length(W.roots)]);
  type:=List(type,function(c)local i;
    for i in [2..Length(c)] do
      c[i]:=ShallowCopy(c[i]);c[i].indices:=OnTuples(c[i-1].indices,phifus);
    od;
    return rec(orbit:=c,
     twist:=RestrictedPerm(phifus^Length(c),c[1].indices),
     operations:=ReflTypeOps);
    end);

  # some adjustment such that in type ^2D_4 the first two simple
  # roots are permuted by res.twist, and in type ^3D_4 the permutation 
  # of the simple roots is 1 -> 2 -> 4 -> 1:

  type:=List(type,function(a) local rf,o,i;    
    i:=a.orbit[1].indices;o:=OrderPerm(a.twist);
    if o=1 or a.orbit[1].series<>"D" or Length(i)<>4 then return a;fi;
    if o=2 then # ^2D_4
      rf:=Filtered([1..4],x->i[x]<>i[x]^a.twist);
      rf:=Concatenation(rf,[3],Difference([1,2,4],rf));
      a.orbit:=List(a.orbit,function(j)
        j:=ShallowCopy(j); j.indices:=j.indices{rf};return j;end);
      return a;
    fi;
    # triality group ^3D_4
    if i[1]^a.twist<>i[2] then 
      a.orbit:=List(a.orbit,function(j)
        j:=ShallowCopy(j);j.indices:=j.indices{[1,4,3,2]};return j;end);
    fi;
    return a;
  end);

##    some sorting: don't do it, it is for isomorphismtype
##    hi:=List(type,a->[a.orbit[1].rank,a.orbit[1].series]);
##    SortParallel(hi,type,function(a,b) 
##                      return a[1]>b[1] or (a[1]=b[1] and a[2]<b[2]); end);  
  return type;
end;

CoxeterCosetOps.CoxeterWord:=function(W,w)return CoxeterWord(Group(W),w); end;

CoxeterCosetOps.CoxeterLength:=function(W,w)return CoxeterLength(Group(W),w);
end;

CoxeterCosetOps.LeftDescentSet:=function(W,w)return LeftDescentSet(Group(W),w);
end;

CoxeterCosetOps.FirstLeftDescending:=function(W,w)return
  FirstLeftDescending(Group(W),w);
end;

CoxeterCosetOps.IsLeftDescending:=function(W,w,i)return
  IsLeftDescending(Group(W),w,i);
end;

CoxeterCosetOps.CoxeterElements:=function(arg)local W;
  W:=arg[1];arg[1]:=Group(W);
  return ApplyFunc(CoxeterElements,arg)*W.phi;
end;

CoxeterCosetOps.ReflectionCharValue:=function(W,w)
  return CoxeterGroupOps.ReflectionCharValue(Group(W),w);
end;

CoxeterCosetOps.CoxeterWords:=function(arg)local a;
  a:=ShallowCopy(arg);a[1]:=Group(a[1]);
  return ApplyFunc(CoxeterWords,a);
end;

CoxeterCosetOps.ParabolicRepresentatives:=function(WF,s)local W,l;
  W:=Group(WF);
  l:=Combinations(W.rootInclusion{W.generatingReflections},s);
  l:=Filtered(l,x->OnSets(Set(x),WF.phi)=Set(x));
  l:=CollectBy(l,x->IsomorphismType(ReflectionSubgroup(W,x)));
  return Concatenation(List(l,v->v{Filtered([1..Length(v)],i->
    not ForAny([1..i-1],
      j->RepresentativeOperation(W,Set(v[i]),Set(v[j]),OnSets)<>false))}));
end;

#############################################################################
##
#F  CoxeterCoset( <W>[, <F0Mat> or  <perm>] ) . . . .  create a CoxeterGroup
#F  coset record
#F  CoxeterCoset( <rec> ) . . . . . . . . . return component <rec>.coxeterCoset
##  
##  In the first form <W> must be a CoxeterGroup record.
##  <F0Mat> must be a square  matrix of rank <W>.rank such that:
##    - <F0Mat> is invertible
##    - <F0Mat> has finite order
##    - <F0Mat>: X->X lets invariant the set of roots of <W> and of Parent(<W>)
##    - TransposedMat(<F0Mat>): Y->Y lets invariant the set of coroots
##        of <W> and of Parent(<W>)
##  
##  If <W>.semisimpleRank  = <W>.rank it   is allowed to give an  argument
##  <perm> (which  is a permutation) instead  of <F0Mat>. Then  <F0Mat> is
##  computed as the  unique matrix which  maps the simple  roots of <W> on
##  the roots given by i^<perm>, i in [1..<W>.rank]. 
##  
##  If only  <W> is  given then the  default for  <F0Mat> is  the identity
##  matrix.
##  
##  'CoxeterCoset' returns a record with the following components:
##    .isDomain     set to true
##    .isFinite     set to true
##    .isCoset      set to true
##    .isCoxeterCoset  set to true
##    .reflectionGroup      <W> from the argument
##    .F0Mat        as described above
##    .phi          unique element in W.<perm> stabilizing simple roots of <W>
##    .operations   CoxeterCosetOps
##    .callArgs     information used by 'Print'
##  
##  In the second form 'CoxeterCoset' returns the component <rec>.coxeterCoset,
##  if present.
##   
CHEVIE.Cache.CoxeterCosets:=true;
CoxeterCoset:=function(arg)local W, PW, perm, tmp, WF, i;

  # just return  component .coxeterCoset:
  if Length(arg)=1 and IsRec(arg[1]) and IsBound(arg[1].coxeterCoset) then
    return arg[1].coxeterCoset;
  fi;

  W:=arg[1];arg:=arg{[2..Length(arg)]};
  WF:=rec(isDomain:=true,
           isFinite:=true,
           isCoxeterCoset:=true,
           isCoset:=true,
           reflectionGroup:=W);

  if Length(arg)>0 then
    if IsPerm(arg[1]) then 
      if Set(W.rootInclusion)<>OnSets(Set(W.rootInclusion),arg[1]) then
	InfoChevie("#I permutation for F0 must normalize set of roots.\n");
	return false;
      fi;
      if IsSubset(W.rootInclusion,MovedPoints(arg[1])) then
           WF.F0Mat:=MatXPerm(W,arg[1]);
      else WF.F0Mat:=MatXPerm(Parent(W),arg[1]);
      fi;
      if W.rank=0 then perm:=();
      else perm:=PermMatX(W,WF.F0Mat);
      fi;
    else WF.F0Mat:=arg[1];
      # checking finite order of F0Mat:
      if W.rank>W.semisimpleRank then OrderMat(WF.F0Mat); fi;

      if WF.F0Mat<>[] then perm:=PermMatX(W,WF.F0Mat);else perm:=();fi;
      if perm=false then
	InfoChevie("#I matrix for F0 must normalize set of roots.\n");
	return false;
      fi;
    fi;
  else WF.F0Mat:=IdentityMat(W.rank);perm:=();
  fi;
  
  PW:=Parent(W);
  # checking if coroots of parent are normalized:
  if IsCoxeterGroup(PW) and PW.semisimpleRank>0 and  # parent may be Spets
    not IsNormalizing(PW.coroots*PW.simpleCoroots,TransposedMat(WF.F0Mat)) then
    InfoChevie("#I transposed of matrix for F0 must normalize set of",
      " coroots of parent.\n");
    return false;
  fi;

  WF.phi:=ReducedInRightCoset(W,perm);
  if WF.phi<>perm then WF.F0Mat:=MatXPerm(Parent(W),WF.phi/perm)*WF.F0Mat;fi;

  WF.operations:=CoxeterCosetOps;
  
  # for printing:
  if Length(arg)=0 or W.rank=0 or WF.F0Mat=WF.F0Mat^0 then WF.callArgs:=[];
  else
    # check for simplification:
    tmp:=W.rootInclusion{[1..W.semisimpleRank]};
    if W.rank=W.semisimpleRank and Set(OnTuples(tmp,WF.phi))=Set(tmp) then
      WF.callArgs:=MappingPermListList(tmp, OnTuples(tmp,WF.phi));
    else
      WF.callArgs:=arg[1];
    fi;
  fi;

  WF.name:=ReflectionName(WF);

  # check for cached cosets. Especially useful for the trivial coset
  return CHEVIE.GetCached(W,"CoxeterCosets",WF,x->x.callArgs);
end;

#############################################################################
##
#F  CoxeterSubCoset( <WF>, <I>[, <w> ] ) . . . . subcoset of a Coxeter coset
##   
##  <I> must be as in 'R:=ReflectionSubgroup(Group(<WF>), <I>)'.
##  <w>  must  be an  element of  Group(<WF>) such that the root system  of
##  'R' is invariant under <w>*<WF>.phi. The default for <w> is ().
##  
CHEVIE.Cache.SubCosets:=false;
CoxeterSubCoset:=function(arg)local WF, W, w, tmp, I, res, i;
  if Length(arg)=1 then arg:=arg[1]; fi;
  if arg[1]=false then return false;fi;
  WF:=Parent(arg[1]);w:=arg[1].phi/WF.phi;W:=Group(WF);
  I:=arg[2];
  if Length(arg)>2 then
    if arg[3] in W then w:=arg[3]*w;
    else Error("must give w in Group(WF).\n");
    fi;
  fi;
  res:=rec(isDomain:=true, isFinite:=true,
    isCoxeterCoset:=true, isCoset:=true, parent:=WF);
  
  res.reflectionGroup:=ReflectionSubgroup(W,I);

  # checking, if w*WF.phi normalizes subroot system:
  tmp:=Set(res.reflectionGroup.rootInclusion);
  if not OnSets(tmp,w*WF.phi)=tmp then
    Error("must give w, such that w * WF.phi normalizes subroot system.\n");
  fi;
  
  res.phi:=ReducedInRightCoset(res.reflectionGroup,w*WF.phi);
  res.F0Mat:=MatXPerm(W,res.phi/WF.phi)*WF.F0Mat;
  res.operations:=CoxeterCosetOps;

  return CHEVIE.GetCached(W,"SubCosets",res,x->[x.reflectionGroup,x.F0Mat]);
end;

CoxeterCosetOps.ReflectionSubgroup:=CoxeterSubCoset;

# Hecke(Coset,H) or Hecke(Coset,para[,rootpara])
# make a hecke coset from Coxeter coset WF and Hecke algebra H or parameters
CHEVIE.Cache.HeckeCosets:=true;
CoxeterCosetOps.Hecke:=function(arg)local W,i,WF,H,res;
  WF:=arg[1];W:=Group(WF);
  if Length(arg)=2 and IsRec(arg[2]) and IsBound(arg[2].reflectionGroup) then
    H:=arg[2];
  else
    arg[1]:=W;H:=ApplyFunc(Hecke,arg);
  fi;
  if W.semisimpleRank>0 then
  for i in Cycles(WF.phi,W.rootInclusion{[1..W.semisimpleRank]}) do
    if Length(Set(H.parameter{W.rootRestriction{i}}))>1 then
      Error("Hecke algebra parameters should be equal for", 
            " phi-conjugate reflections");
    fi;
  od;
  fi;
  res:=rec(hecke:=H,spets:=WF,isDomain:=true,operations:=HeckeCosetOps);
  return CHEVIE.GetCached(WF,"HeckeCosets",res,x->[x.hecke.parameter,
       x.hecke.rootParameter]);
end;

CoxeterCosetOps.RelativeCoset:=function(WF,J)local res,p;
# Print("CoxeterCosetOps.RelativeCoset ",WF,J," called \n");
  res:=RelativeGroup(Group(WF),J);
  if J=[] then p:=WF.phi;
  else p:=PermListList(res.parentMap,OnTuples(res.parentMap,WF.phi));
  fi;
  return CoxeterCoset(res,p);
end;

########################################################################
##
#F  TwistedPower(n,x,F) let x be an element of a monoid on which F acts.
##  this function returns  (xF)^n 
##
TwistedPower:=function(n,x,F)  
  if n=0 then return x^0;
  else return x* F(TwistedPower(n-1,x,F));
  fi;
end;

########################################################################
# next functions apply to any spets

TwistingElements:=function(WF,J)local L,e,W,W_L,H,h,N,WF_L,gens,i;
  if IsSpets(WF) then W:=Group(WF); else W:=WF;WF:=Spets(W);fi;
  if J=[] then return List(ConjugacyClasses(WF),Representative)*WF.phi^-1;fi;
  if IsCoxeterGroup(W) 
#   and IsSubset(W.rootInclusion{W.generatingReflections},J) 
  then
    h:=RepresentativeOperation(W,OnSets(Set(J),WF.phi),Set(J),OnSets);
    if h=false then 
      InfoChevie("\n   # no subspets for ",J);
      return [];
    fi;
    W_L:=Stabilizer(W,Set(J),OnSets);
    e:=List(ConjugacyClasses(Group(Concatenation(W_L.generators,
       [WF.phi*h]),W_L.identity)),Representative);
    e:=Filtered(e,x->WF.phi*h*x^-1 in W_L);
    return e*WF.phi^-1;
  fi;
  L:=ReflectionSubgroup(W,J); N:=Normalizer(W,L); W_L:=FactorGroup(N,L);
  if WF.phi<>() then Error("not implemented for twisted parent Spets");fi;
  if Size(W_L)>=10 then
    H:=Group(List(W_L.generators,x->GetRelativeAction(W,L,
         Representative(x.element))),IdentityMat(W.rank-L.semisimpleRank));
    if Size(H)=Size(W_L) then 
      h:=GroupHomomorphismByImages(H,W_L,H.generators,W_L.generators);
      h.isMapping:=true;
      e:=List(ConjugacyClasses(H),Representative);
      e:=List(e,x->Representative(Image(h,x).element));
      return e;
    fi;
  fi;
  return List(ConjugacyClasses(W_L),x->Representative(Representative(x).element));
end;

########################################################################
##
## Twistings( <W>, <L> )'
##
## <L>  should be a reflection subgroup of group or coset <W>, or a sublist
## of the generating reflections.
## Returns  a list of  representatives, up to  <W>-conjugacy, of reflection
## sub-cosets whose reflection group is <L>.
Twistings:=function(arg)local WF,J,tt,t,i,gens,W;
  WF:=arg[1];
  if Length(arg)=1 then
    W:=WF;
    tt:=CollectBy(ReflectionType(W),ReflectionName);
    gens:=[];
    for t in tt do
      for i in [1..Length(t)-1] do
        Add(gens,Product(Zip(W.rootInclusion{t[i].indices},
                             W.rootInclusion{t[i+1].indices},
          function(i,j)return (i,j);end)));
      od;
      J:=W.rootInclusion{t[1].indices};
      if t[1].series="A" then 
        if t[1].rank>1 then Add(gens,
          Product([1..QuoInt(t[1].rank,2)],i->(J[i],J[t[1].rank+1-i])));
        fi;
      elif t[1].series="D" then Add(gens,(J[1],J[2]));
        if t[1].rank=4 then Add(gens,(J[1],J[4]));fi;
      elif t[1].series="E" and t[1].rank=6 then
        Add(gens,(J[1],J[6])(J[3],J[5]));
      fi;
    od;
    if W<>Parent(W) then
      Error(W," must not be a proper subgroup of another reflection group");
    fi;
    return List(Elements(Group(gens,())),x->CoxeterCoset(W,x));
  fi;
  J:=arg[2];
  if IsGroup(J) then J:=J.rootInclusion{J.generatingReflections};fi;
  if not IsSpets(WF) then WF:=Spets(WF);fi;
  return Filtered(List(TwistingElements(WF,J),
       x->SubSpets(WF,J,x)),x->not IsIdentical(x,false));
end;

Torus:=function(arg)local W,res,i;
  W:=arg[1];
  if IsInt(W) then 
    res:=Copy(CoxeterGroup());Unbind(res.degrees);res.rank:=W;
    res.name:=SPrint("Torus(",W,")");return res;
  elif IsMat(W) then 
    res:=Copy(CoxeterGroup());Unbind(res.degrees);res.rank:=Length(W);
    return Spets(res,W);
  fi;
  i:=arg[2];
  if IsSpets(W) then 
    return SubSpets(W,[],Representative(ConjugacyClasses(W)[i])/W.phi);
  fi;
  return Spets(ReflectionSubgroup(W,[]),Representative(ConjugacyClasses(W)[i]));
end;

ReadChv("prg/wclsinv"); #for CoxeterCosetOps.ClassInvariants
ReadChv("prg/sscoset"); #for CoxeterCosetOps.Centralizer,
#  QuasiIsolatedRepresentatives, IsIsolated
