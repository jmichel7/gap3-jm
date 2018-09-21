#############################################################################
##
#A  garside.g                 CHEVIE library                      Jean Michel.
##
#Y  Copyright (C) 2002-   The CHEVIE Team
##
##  This file contains functions to work with elements of (locally) Garside
##  monoids
#  Maybe the best illustration is the session log below
#  
#  gap> w:=DualBraid(CoxeterGroup("A",4))(1,2,3,4,3);
#  1234.3
#  gap> w*w;
#  1234.1249.34
#  gap> w^4;
#  1234.1249.1249.1249.267.4.
#
#  Elements are stored  as delta^pd w_1.w_2...w_n where pd is positive or
#  negative and w_1.w_2...w_n is the  left-greedy canonical form: w_1 is
#  the  longest simple  which is  a left  prefix, then  the rest  of the
#  element is  recursively decomposed  in the same  way. This  is called
#  "left-greedy normal  form". There is no  delta if the monoid  is only
#  locally Garside
#
#  Elements are represented as records with fields:
#
#   .pd       the power of delta (positive or negative -- may be absent)
#             called inf by Birman et al.
#   .elm      the canonical left-greedy decomposition (a list of simples)
#             Length(elm) is called 'canonical length' by Birman et al.
#   .monoid   a pointer to the monoid info
#   .operations=GarsideEltOps  see below
#
#   Most examples are 'generated groups' where simples are identified to
#   elements of some finite group.
#############################################################################

GarsideEltOps:=OperationsRecord("GarsideEltOps");

# suppress instances of delta and identity in .elm
# not its business no check that .elm is a normal form otherwise.
GarsideEltOps.Normalize:=function(b)local res,s,delta; 
  if IsBound(b.monoid.delta) then delta:=b.monoid.delta;else delta:=0;fi;
  res:=[];
  for s in b.elm do
    if s=delta then b.pd:=b.pd+1;
    elif s<>b.monoid.identity then Add(res,s);
    fi;
  od;
  b.elm:=res;return b;
end;

GarsideEltOps.IsLeftDescending:=function(b,s)
  if IsBound(b.pd) then
    if b.pd<0 then Error("only for elements of the Garside monoid");
    elif b.pd>0 then return true;
    fi;
  fi;
  if Length(b.elm)=0 then return false;fi;
  return b.monoid.IsLeftDescending(b.elm[1],s);
end;

GarsideEltOps.FirstLeftDescending:=function(b)
  if IsBound(b.pd) then
    if b.pd<0 then Error("only for elements of the Garside monoid");
    elif b.pd>0 then return 1;
    fi;
  fi;
  if Length(b.elm)=0 then return false;fi;
  return PositionProperty([1..b.monoid.nrAtoms],i->
    b.monoid.IsLeftDescending(b.elm[1],i));
end;

GarsideEltOps.LeftDescentSet:=function(b)
  if IsBound(b.pd) then
    if b.pd<0 then Error("only for elements of the Garside monoid");
    elif b.pd>0 then return [1..b.monoid.nrAtoms];
    fi;
  fi;
  if Length(b.elm)=0 then return [];fi;
  return Filtered([1..b.monoid.nrAtoms],
    i->b.monoid.IsLeftDescending(b.elm[1],i));
end;

GarsideEltOps.Format:=function(b,option)local M,res,p,seqCollected,W;
  M:=b.monoid;
  seqCollected:=function(l)local res,v,i;
    if l=[] then return [];fi;
    res:=[];v:=[l[1],1];
    for i in [2..Length(l)] do
      if l[i]=v[1] then v[2]:=v[2]+1; else Add(res,v);v:=[l[i],1]; fi;
    od;
    Add(res,v);
    return res;
  end;
  p:=function(b)local l;
    if IsBound(option.Brieskorn) then
      return Join(List(seqCollected(BrieskornNormalForm(b)),function(x)local s;
        W:=M.group;
        if x[1]=W.rootInclusion{W.generatingReflections} then s:=M.stringDelta;
	elif Length(x[1])=1 then s:=String(x[1][1]);
	else s:=SPrint("w_",IntListToString(x[1]));
	fi;
	return FormatMonomial(s,x[2],option);end),".");
    fi;
    l:=List(b.elm,v->M.FormatSimple(v,option));
    if IsBound(option.GAP) then l:=Join(l,","); else l:=Join(l,"."); fi;
    if not IsBound(b.pd) or b.pd=0 then return l;fi;
    if IsBound(option.GAP) then return SPrint("[",l,"],",b.pd);
    else
      if l<>"" then l:=Concatenation(".",l);fi;
      if IsBound(option.Greedy) then return SPrint(M.stringDelta,"^",b.pd,l);
      else return SPrint(Join(List([1..b.pd],i->M.stringDelta),"."),l);
      fi;
    fi;
  end;
  if IsBound(option.Greedy) or not IsBound(b.pd) then b:=["",p(b)];
  else b:=List(AsFraction(b),p);
  fi;
  if IsBound(option.GAP) then res:=List(b,x->SPrint("B(",x,")"));
       res:=[SPrint(res[1],"^-1"),"*",res[2]];
  else res:=[SPrint("(",b[1],")^-1"),".",b[2]];
  fi;
  if b[1]<>"" then 
    if b[2]<>"" then return Concatenation(res);
    else return res[1];
    fi;
  elif b[2]<>"" or IsBound(option.GAP) then return res[3];
  else return ".";
  fi;
end;

GarsideEltOps.String:=b->GarsideEltOps.Format(b,rec());

#F Print(b) print b (using Format)

GarsideEltOps.Print:=function(b)Print(Format(b,CHEVIE.PrintGarside));end;

if not IsBound(CHEVIE.PrintGarside) then CHEVIE.PrintGarside:=rec();fi;

GarsideEltOps.\*:=function(a,b)local M,x,res;
  if IsList(a) then return List(a,x->x*b);
  elif IsList(b) then return List(b,x->a*x);
  fi;
  M:=a.monoid;
  if M<>b.monoid then 
    Error("not elements of the same Garside monoid or group");
  fi;
  res:=ShallowCopy(a);res.elm:=ShallowCopy(res.elm);
  if IsBound(M.delta) and b.pd<>0 then
    res.pd:=a.pd+b.pd;res.elm:=List(res.elm,x->M.DeltaAction(x,b.pd));
  fi;
  for x in b.elm do res.elm:=M.AddToNormal(res.elm,x); od;
#   res:=M.Elt(ShallowCopy(b.elm)); # variant
#   if IsBound(M.delta) then
#     res.pd:=a.pd+b.pd;a:=List(a.elm,x->M.DeltaAction(x,b.pd));
#   else a:=a.elm;
#   fi;
#   for x in Reversed(a) do res.elm:=M.PrefixToNormal(x,res.elm);od;
  return GarsideEltOps.Normalize(res);
end;

GarsideEltOps.\/:=function(a,b)return a*b^-1;end;

GarsideEltOps.\^:=function(b,n)local M,res;
  if not IsInt(n) then return n^-1*b*n;fi;
  M:=b.monoid;
  if n<0 then
    if not IsBound(M.delta) then 
      Error("negative power in a monoid which is not Garside");
    fi;
    b:=M.Inverse(b); n:=-n;
  fi;
  if n=1 then return b;fi;
  res:=M.Elt([]);
  while true do
    if n mod 2 <> 0 then res:=res*b;fi;
    n:=QuoInt(n,2);
    if n=0 then return res;fi;
    b:=b*b;
  od;
end;

# Algorithm following Hao Zheng "A new approach to extracting roots in Garside
# groups" Comm. Algebra 34 (2006) 1793--1802
GarsideEltOps.GetRoot:=function(b0,n)local b,conj,M,tM,inner,a,l,k;
  if not IsBound(b0.pd) or b0.pd<0 then Error("Only for Garside monoids");fi;
  M:=b0.monoid;
  tM:=TwistedPowerMonoid(M,n);
  b:=tM.Elt(List(Concatenation(List([1..b0.pd],x->M.delta),b0.elm),
    function(s)local v;v:=List([1..n],i->M.identity);
      v[1]:=s;return tM.atom(v,false);end));b.elm[1].t:=true;
  inner:=function(b)local cst,sc,class,a,m,e;
    cst:=b->ForAll(b.elm,x->Length(Set(x.v))=1);
    sc:=RepresentativeSC(b); conj:=[sc.conj]; class:=[sc.circuit[1]];
    if cst(class[1]) then return conj[1];fi;
    for a in class do for m in AtomicMaps(a) do
      if not m.tgt in class then
        e:=conj[Position(class,a)]*m.map;
#       Print(".\c");
        if cst(m.tgt) then return e;fi;
        Add(class,m.tgt); Add(conj,e);
      fi;
    od; od;
    return false;
  end;
  conj:=inner(b);if conj=false then return false;fi;
  a:=b^conj;
  a:=M.Elt(List(a.elm,x->x.v[1]),a.pd);
  l:=Length(conj.elm);
  k:=Number(conj.elm,x->x.t);
  conj:=M.Elt(Concatenation(List([1..k],i->conj.elm[i].v[1+(i mod n)]),
                List([k+1..l],i->conj.elm[i].v[1+(k mod n)])),conj.pd);
  return conj*a*conj^-1;
end;

########################################################################
##
#F  GarsideAlpha(b) returns the largest simple (for divisiblity) dividing b
##
GarsideAlpha:=function(b)
  if IsBound(b.pd) then
    if b.pd<0 then Error("GarsideAlpha: element should be positive");
    elif b.pd>0 then return b.monoid.delta;
    fi;
  fi;
  if Length(b.elm)>0 then return b.elm[1];
  else return b.monoid.identity;
  fi;
end;

########################################################################
##
#F  GarsideOmega(b) returns Alpha(b)^-1 * b
##
GarsideOmega:=function(b)
  if IsBound(b.pd) then
    if b.pd<0 then Error("GarsideAlpha: element should be positive");
    elif b.pd>0 then b:=ShallowCopy(b);b.pd:=b.pd-1;return b;
    fi;
  fi;
  if Length(b.elm)>0 then 
    b:=ShallowCopy(b);b.elm:=b.elm{[2..Length(b.elm)]};return b;
  else return b.monoid.identity;
  fi;
end;

#############################################################################
##
#F LeftGcd(a1,..,an) returns [d,d^-1a1,..,d^-1a_n] where d=left gcd(a1,..,an)
##
LeftGcd:=function(arg)local x,gcd,M,m;
  if Length(arg)=0 then Error("LeftGcd needs an argument");
  elif Length(arg)=1 then return [arg[1],arg[1]^0];
  fi;
  M:=arg[1].monoid;
  if IsBound(arg[1].pd) then 
    m:=Minimum(List(arg,x->x.pd));
    arg:=List(arg,b->Concatenation(List([1..b.pd-m],i->M.delta),b.elm));
  else arg:=List(arg,b->b.elm);
  fi;
  gcd:=[];
  while true do
    if ForAny(arg,b->Length(b)=0) then x:=[M.identity];
    else x:=ApplyFunc(M.LeftGcdSimples,List(arg,x->x[1]));
    fi;
    if x[1]=M.identity then
      if IsBound(M.delta) then gcd:=M.Elt(gcd,m);else gcd:=M.Elt(gcd);fi;
      return List(Concatenation([gcd],List(arg,M.Elt)),GarsideEltOps.Normalize);
    else gcd:=M.AddToNormal(gcd,x[1]);
      arg:=List([1..Length(arg)],
        i->M.PrefixToNormal(x[i+1],arg[i]{[2..Length(arg[i])]}));
    fi;
  od;
end;

#############################################################################
##
#F  LeftLcm(a1,..,an) returns [m,ma1^-1,..,man^-1] where m=Left lcm(a1,..,an)
##
LeftLcm:=function(arg)local M,x;
  if Length(arg)=0 then Error("LeftLcm needs an argument");
  elif Length(arg)=1 then return [arg[1],arg[1]^0];
  fi;
  M:=arg[1].monoid; 
  if not IsBound(M.Inverse) then 
    Error("LeftLcm: should be elements of a Garside group");
  fi;
  x:=ApplyFunc(LeftGcd,List(arg,M.Inverse));
  x[1]:=M.Inverse(x[1]);
  return x;
end;

#########################################################################
##
#F ReversedWord
##
ReversedWord:=function(b)
  if not IsBound(b.monoid.Reverse) then
    Error(b.monoid," does not have a Reverse operation");
  fi;
  return b.monoid.Reverse(b);
end;

########################################################################
##
#F  RightLcm(a,b) returns triple [m,a^-1m,b^-1m] where m=right lcm(a,b)
##
RightLcm:=function(arg)
  return List(ApplyFunc(LeftLcm,List(arg,ReversedWord)),ReversedWord);
end;

########################################################################
##
#F  RightGcd(a,b) returns triple [d,ad^-1,bd^-1] where d=right gcd(a,b)
##
RightGcd:=function(arg)
  return List(ApplyFunc(LeftGcd,List(arg,ReversedWord)),ReversedWord);
end;

#############################################################################
##
#F  AsFraction(b) . . . . . . . . . . . .  Fraction normal form. Returns
##  a pair of elements x,y of monoid with trivial left gcd and b=x^-1 y
##
AsFraction:=function(b)local M; M:=b.monoid;
  if not IsBound(M.Inverse) then 
    Error("AsFraction: should be element of a Garside group");fi;
  if b.pd>=0 then return [M.Elt([]),b];fi;
  if -b.pd>Length(b.elm) then return [M.Inverse(b),M.Elt([])];fi;
  return [M.Elt(M.Inverse(b).elm{[Length(b.elm)+b.pd+1..Length(b.elm)]}),
     M.Elt(b.elm{[1-b.pd..Length(b.elm)]})];
end;

############################################################################
##
#F  AsWord(w) . . . . . . . . . . . . . . . .  sequence of atoms
##
##  w  is  an  element  of  a  Garside  monoid.  AsWord(w)  returns  the
##  expression of w as  a list of atoms given as  a sequence of integers
##  (the index in the list of  atoms). As a special case (symmetric with
##  the input methods of 'Braid') an  element which is not in the monoid
##  is returned in fraction normal form,  the inverse of the atoms coded
##  as negative integers.

AsWord:=function(b)local res,i,M;
  M:=b.monoid;
  res:=[];
  if IsBound(b.pd) then
    if b.pd<0 then 
      b:=AsFraction(b);
      return Concatenation(-Reversed(AsWord(b[1])),AsWord(b[2]));
    fi;
    for i in [1..b.pd] do Append(res,M.AtomListSimple(M.delta));od;
  fi;
  for i in b.elm do Append(res,M.AtomListSimple(i));od;
  return res;
end;

# A (locally) Garside monoid M is a record with at least the following fields
# defined:
#   .IsLeftDescending(s,i), .IsRightAscending(s,i), 
#   .IsLeftAscending(s,i), .IsRightDescending(s,i),
# where s is a simple and i an atom index
#   .atoms: list of atoms
#   .identity: the trivial simple
#   .delta, .DeltaAction :     should be bound if M is Garside
#   .LeftQuotient(a,b) .RightQuotient(a,b) .Product(a,b) (on simples when legal)
#
# The next function ensures that the following methods are available in the
# Garside monoid M:
#  nrAtoms, Elt, AtomListSimple, LeftGcdSimples, RightGcdSimples, alpha2,
#  AddToNormal, FormatSimple
#
#  and if M is Garside:
#  RightLcmSimples, LeftLcmSimples
#
#  The second argument is a record specifying additional properties:
#  .interval  if bound define automatically LeftQuotient, RightQuotient, Product
#  .eltOps    if bound add these methods to GarsideEltOps for  elements of M.
CompleteGarsideRecord:=function(M,opt)local eltops;
  if IsBound(opt.eltOps)then
       eltops:=ShallowCopy(GarsideEltOps);Inherit(eltops,opt.eltOps);
  else eltops:=GarsideEltOps;
  fi;
  if not IsBound(M.Elt) then
    M.Elt:=function(arg)local res; # M.elt(elm[,pd])
       res:=rec(elm:=arg[1],operations:=eltops,monoid:=M);
       if Length(arg)>1 then res.pd:=arg[2];
       elif IsBound(M.delta) then res.pd:=0;
       fi;
       return res;
    end;
  fi;
  if IsBound(opt.interval) then
    M.interval:=true;
    M.LeftQuotient:=function(a,b)return a^-1*b;end;
    M.RightQuotient:=function(a,b)return a/b;end;
    M.Product:=function(a,b)return a*b;end;
    if IsBound(M.delta) then
      M.identity:=M.delta^0;
      M.RightComplementToDelta:=x->x^-1*M.delta;
      M.LeftComplementToDelta:=x->M.delta/x;
      if not IsBound(M.DeltaAction) then
        if IsBound(M.orderDelta) then
          M.DeltaAction:=function(s,i)
            return s^(M.delta^(i mod M.orderDelta));end;
        else M.DeltaAction:=function(s,i)return s^(M.delta^i);end;
        fi;
      fi;
      if M.delta^2=M.identity then
        M.ReverseSimple:=x->x^-1;
      fi;
    else M.identity:=M.atoms[1]^0;
    fi;
  fi;
  if IsBound(M.atoms) then M.nrAtoms:=Length(M.atoms);fi;
  if not IsBound(M.AtomListSimple) then
    M.AtomListSimple:=function(s)local i,res;
      res:=[]; i:=1;
      while i<=M.nrAtoms do
        if M.IsLeftDescending(s,i) then s:=M.LeftQuotient(M.atoms[i],s);
             Add(res,i);i:=1;
        else i:=i+1;
        fi;
      od;
      if s<>M.identity then Error("monoid seems not Garside");fi;
      return res;
    end;
  fi;
  if not IsBound(M.FormatSimple) then
    M.FormatSimple:=function(s,opt)local l;l:=M.AtomListSimple(s);
      if IsBound(opt.GAP) then return Join(l);
      else return IntListToString(l,"[]");
      fi;
    end;
  fi;
  if IsBound(M.delta) then
    if not IsBound(M.stringDelta) then
      M.stringDelta:="D";
    fi;
    if not IsBound(M.RightComplementToDelta) then
      M.RightComplementToDelta:=x->M.LeftQuotient(x,M.delta);
    fi;
    if not IsBound(M.LeftComplementToDelta) then
      M.LeftComplementToDelta:=x->M.RightQuotient(M.delta,x);
    fi;
    if not IsBound(M.RightLcmSimples) then
      M.RightLcmSimples:=function(arg)local x;
	x:=ApplyFunc(M.RightGcdSimples,List(arg,M.RightComplementToDelta));
        x[1]:=M.Product(arg[1],x[2]);return x;
      end;
    fi;
    if not IsBound(M.LeftLcmSimples) then
      M.LeftLcmSimples:=function(arg)local x;
	x:=ApplyFunc(M.LeftGcdSimples,List(arg,M.LeftComplementToDelta));
        x[1]:=M.Product(x[2],arg[1]);return x;
      end;
    fi;
    if not IsBound(M.DeltaAction) then
      M.DeltaAction:=function(x,i)
        while i>0 do
          x:=M.LeftQuotient(M.RightComplementToDelta(x),M.delta);i:=i-1;od;
        while i<0 do
          x:=M.RightQuotient(M.delta,M.LeftComplementToDelta(x));i:=i+1;od;
        return x;
      end;
    fi;
    if not IsBound(M.Inverse) then
      M.Inverse:=function(b)local k;k:=Length(b.elm);
       return M.Elt(List([k,k-1..1],i->M.DeltaAction(
         M.RightComplementToDelta(b.elm[i]),-i-b.pd)),-b.pd-k);
      end;
    fi;
  fi;
  if not IsBound(M.LeftGcdSimples) then
  # M.LeftGcdSimples(a1,..,an) returns [gcd,gcd^-1*a1,...,gcd^-1*an]
    M.LeftGcdSimples:=function(arg)local i,s,x;x:=M.identity;
      repeat
        s:=M.identity;
        for i in [1..M.nrAtoms] do
	  if ForAll(arg,b->M.IsLeftDescending(b,i)) then
	    s:=M.Product(s,M.atoms[i]);arg:=M.LeftQuotient(M.atoms[i],arg);
          fi;
        od;
        x:=M.Product(x,s);
      until s=M.identity;
      return Concatenation([x],arg);
    end;
  fi;
  if not IsBound(M.RightGcdSimples) then
    if IsBound(M.ReverseSimple) then
    M.RightGcdSimples:=function(arg)return List(ApplyFunc(M.LeftGcdSimples,
      List(arg,M.ReverseSimple)),M.ReverseSimple);end;
    else
    M.RightGcdSimples:=function(arg)local i,x;
      i:=1;x:=M.identity;
      while i<=M.nrAtoms do
	if ForAll(arg,b->M.IsRightDescending(b,i)) then
          x:=M.Product(M.atoms[i],x);arg:=M.RightQuotient(arg,M.atoms[i]);
	  i:=1;
	else i:=i+1;
	fi;
      od;
      return Concatenation([x],arg);
    end;
    fi;
  fi;
  if not IsBound(M.alpha2) then
    # Given x,v simples returns x', v' such that x' * v' = x * v and where
    #  x' = x * (longest subword z of v such that xz is simple)
    M.alpha2:=function(x,v)local i;
      if IsBound(M.delta) then
        i:=M.LeftGcdSimples(M.RightComplementToDelta(x),v);
        return [M.Product(x,i[1]),i[3]];
      else
	i:=1;
	repeat
	  if M.IsLeftDescending(v,i) and M.IsRightAscending(x,i)
	  then x:=M.Product(x,M.atoms[i]);v:=M.LeftQuotient(M.atoms[i],v);i:=1;
	  else i:=i+1;
	  fi;
	until i>M.nrAtoms;
	return [x,v];
      fi;
    end;
  fi;
  if not IsBound(M.AddToNormal) then
    # given a normal form v and a simple x return the normal form of v.x,
    # discarding a final 1 but not an initial delta in the result.
    # Warning! This modifies v. Make a ShallowCopy if needed.
    M.AddToNormal:=function(v,x)local i,l;
      Add(v,x);
      for i in [Length(v),Length(v)-1..2] do
        x:=M.alpha2(v[i-1],v[i]);
        if x[2]=M.identity then v:=v{[1..i-1]};v[i-1]:=x[1]; # i=Length(v)
	elif x[1]=v[i-1] then return v;
 	elif IsBound(M.delta) and x[1]=M.delta then
 	  v[i]:=x[2];v{[2..i-1]}:=List(v{[1..i-2]},x->M.DeltaAction(x,1));
 	  v[1]:=x[1];
 	  return v;
        else v{[i-1,i]}:=x;
	fi;
      od;
      return v;
    end;
  fi;
  if not IsBound(M.PrefixToNormal) then
    # given a normal form v and a simple x return the normal form of x.v ;
    # the result may have one initial delta.
    M.PrefixToNormal:=function(x,v)local res,i;
      res:=[];
      for i in [1..Length(v)] do
        x:=M.alpha2(x,v[i]);Add(res,x[1]);x:=x[2];
	if x=v[i] then Append(res,v{[i..Length(v)]});return res;
	elif x=M.identity then Append(res,v{[i+1..Length(v)]});return res;
	fi;
      od;
      Add(res,x);
      return res;
    end;
  fi;
  if not IsBound(M.Reverse) and IsBound(M.ReverseSimple) then
    if IsBound(M.delta) then
      M.Reverse:=function(b)local res,s;
        if Length(b.elm)=0 then return b;fi;
        res:=[];
        for s in Reversed(b.elm) do
          res:=M.AddToNormal(res,M.DeltaAction(M.ReverseSimple(s),-b.pd));
        od;
        return GarsideEltOps.Normalize(M.Elt(res,b.pd));
      end;
    else
      M.Reverse:=function(b)local res,s;
        if Length(b.elm)=0 then return b;fi;
        res:=[];
        for s in Reversed(b.elm) do
          res:=M.AddToNormal(res,M.ReverseSimple(s));
        od;
        return GarsideEltOps.Normalize(M.Elt(res));
      end;
    fi;
  fi;
  if not IsBound(M.B) then
##  M.B: a function which can be called in the forms
#F  or B([s1,...,sn] [, pd])  
#F  or B(s1,...sn)  
#
#  and  returns  an  element  of  M,  the  product  of the simples of index
#  s1,...,sn.  As a special, if any of the s_i is a negative integer, it is
#  taken  to  represent  the  inverse  of  the corresponding atom. If pd (a
#  positive  or  negative  integer)  is  given,  the  resulting  element is
#  multiplied by M.delta^pd.
  M.B:=function(arg)local res,s,p;
    res:=M.Elt([]);
    if Length(arg)<>0 then
      if IsList(arg[1]) then 
        if Length(arg)=2 then res.pd:=arg[2];fi;
        arg:=arg[1];
      fi;
      for s in Reversed(arg) do # faster in reversed order (see \*)
        if IsInt(s) and s<0 then p:=-1;s:=-s;else p:=1;fi;
        res:=M.Elt([M.atoms[s]])^p*res; 
      od;
    fi;
    return res;
  end;
  fi;
end;

#############################################################################
##
#F  GarsideWords(M,l)  . . . . . . . . . . .  monoid elements of length l
##
##  M is a (locally) Garside monoid with atoms of length 1.
##  GarsideWords(W,l) returns the set of elements of M of length l.

GarsideWords:=function(M,l)local i,bb,w,b,B,s,r,rr,RA,lb;
  if l=0 then return [M.Elt([])];fi;
  if not IsBound(M.words) then 
    M.words:=[Set(List(M.atoms,x->M.Elt([x])))];
  fi;
  RA:=w->Filtered([1..M.nrAtoms],i->M.IsRightAscending(w,i));
  if not IsBound(M.words[l]) then 
    M.words[l]:=[];
    for b in GarsideWords(M,l-1) do
      lb:=Length(b.elm);
      if lb=0 then r:=[];
      else
        r:=RA(b.elm[lb]);if lb=1 then rr:=[]; else rr:=RA(b.elm[lb-1]);fi;
        for s in r do
          w:=M.Product(b.elm[lb],M.atoms[s]);
          if not ForAny(rr,i->M.IsLeftDescending(w,i)) then
            bb:=ShallowCopy(b);bb.elm:=ShallowCopy(bb.elm);
            bb.elm[Length(bb.elm)]:=w;
            Add(M.words[l],GarsideEltOps.Normalize(bb));
          fi;
        od;
      fi;
      for s in Drop(M.atoms,r) do
        bb:=ShallowCopy(b);bb.elm:=Concatenation(bb.elm,[s]);
        Add(M.words[l],GarsideEltOps.Normalize(bb));
      od;
    od;
    M.words[l]:=Set(M.words[l]);
  fi;
  return M.words[l]; 
end;

# LeftDivisorsSimple(<M>,<s>[,<i>]) 
# all divisors of the simple s of monoid M [of length i]
LeftDivisorsSimple:=function(arg)local M,s,res,i,rest,x,new;
  M:=arg[1];s:=arg[2];
  rest:=[[M.identity,s]];
  res:=[];
  while Length(rest)>0 do
    if Length(arg)=3 and Length(res)=arg[3] then return List(rest,y->y[1]);fi;
    Add(res,rest);
    new:=[];
    for x in rest do
      for i in [1..M.nrAtoms] do
        if M.IsLeftDescending(x[2],i) then
          Add(new,[M.Product(x[1],M.atoms[i]),M.LeftQuotient(M.atoms[i],x[2])]);
        fi;
      od;
    od;
    rest:=Set(new);
#   InfoChevie("new:",Length(new)," different:",Length(rest),"\n");
  od;
# return List(res,x->GarsideEltOps.Normalize(M.Elt([x[1]])));
  return List(res,x->List(x,y->y[1]));
end;

#############################################################################
##                                                                         ##
##  First example of Garside monoid: the Artin-Tits braid monoid           ##
##                                                                         ##
#############################################################################
#  Maybe the best illustration is the session log below
#  
#  gap> w:=Braid(CoxeterGroup("A",4))(1,2,3,4,3);
#  12343
#  gap> w*w;
#  12132432.43
#  gap> w^4;
#  w0^2.
#############################################################################

BraidMonoid:=function(W)local M,ILD,l;
  if not IsCoxeterGroup(W) then Error("W should be a Coxeter Group");fi;
  if IsBound(W.braidMonoid) then  return W.braidMonoid;fi;
  ILD:=W.operations.IsLeftDescending; # avoid dispatching overhead 
  M:=rec(group:=W,
    operations:=rec(Print:=function(M)Print("BraidMonoid(",W,")");end),
    IsLeftDescending:=function(w,i)return ILD(W,w,i);end,
    IsLeftAscending:=function(w,i)return not ILD(W,w,i);end,
    IsRightAscending:=function(w,i)return not ILD(W,w^-1,i);end,
    IsRightDescending:=function(w,i)return ILD(W,w^-1,i);end,
    AtomListSimple:=v->W.operations.CoxeterWord(W,v),
    atoms:=W.reflections{W.generatingReflections}
  );
  if IsFinite(W) then 
    M.delta:=LongestCoxeterElement(W);
    M.orderDelta:=OrderPerm(PermListList(M.atoms,OnTuples(M.atoms,M.delta)));
    M.stringDelta:="w0";
  fi;
  CompleteGarsideRecord(M,rec(interval:=true,eltOps:=rec(
    Frobenius:=function(W,x,i)
         x:=ShallowCopy(x);x.elm:=List(x.elm,y->Frobenius(W)(y,i));
         return x;
    end,
    BrieskornNormalForm:=function(b)local W,res,v;
	 W:=b.monoid.group;
	 res:=[];
	 if IsBound(b.pd) then
	 if b.pd<0 then Error(b," should be in the Braid Monoid");fi;
	 for v in [1..b.pd] do Add(res,W.reflectionsLabels);od;
       fi;
       for v in b.elm do Append(res,BrieskornNormalForm(W,v));od;
       return res;
     end)));
  if IsBound(W.operations.LeftGcdSimples) then
    M.LeftGcdSimples:=W.operations.LeftGcdSimples;
  fi;
  W.braidMonoid:=M;
##  M.B: a function which can be called in the forms
#F     B(<element of W> [, pw0])  
#F  or B([s1,...,sn] [, pw0])  
#F  or B(s1,...sn)  
##
##  returns  an element  of the  Artin braid  monoid of  W: either  the
##  reduced braid  corresponding to an  element of the Coxeter  group W,
##  or  (non  necessarily  reduced)  braid  product  of  the  generators
##  corresponding to the  labels for Coxeter generators  s1,...,sn. As a
##  special case, if any  of the s_i is a negative  integer, it is taken
##  to represent the inverse of the corresponding generator.
##  If  pw0 (a  positive or  negative integer)  is given,  the resulting
##  element is multiplied by w0^pw0.
#############################################################################
  M.B:=function(arg)local res,s,r,pd;
    res:=M.Elt([]);
    r:=x->M.Elt([W.reflections[W.operations.ReflectionFromName(W,x)]]);
    if Length(arg)<>0 then
      if not IsWordFor(W,arg[1]) and not arg[1] in W.reflectionsLabels 
          and not (IsInt(arg[1]) and -arg[1] in W.reflectionsLabels)
      then
      # assume h is a group element
        if arg[1]<>W.identity then res.elm:=[arg[1]];fi;
        if Length(arg)=2 then res.pd:=arg[2];fi;
        return GarsideEltOps.Normalize(res);
      else
        pd:=0;
        if IsList(arg[1]) then 
          if Length(arg)=2 then pd:=arg[2];fi;
          arg:=arg[1];
        fi;
        for s in Reversed(arg) do # faster in reversed order (see \*)
          if IsInt(s) and s<0 then res:=r(-s)^-1*res; else res:=r(s)*res;fi;
        od;
        if IsBound(res.pd) then res.pd:=res.pd+pd;fi;
      fi;
    fi;
    return res;
  end;
  return M;
end;

#############################################################################
##
#F  Braid(W) . . . . . . . . . . . . . returns braid-element-making function
##
#############################################################################
Braid:=W->BraidMonoid(W).B;

#############################################################################
##
#F  EltBraid(w)  . . . . . . . . . . . . . . .  Image in the Coxeter group.
##
##  returns the image in w.monoid.group of w (works for interval monoids)
##
EltBraid:=function(x)local res,M;
  M:=x.monoid;
  if not IsBound(M.interval) then 
    Error(x," should be an element of an interval monoid"); 
  fi;
  if IsBound(x.pd) then res:=M.delta^x.pd;
  else res:=M.identity;
  fi;
  if Length(x.elm)=0 then return res;
  else return res*Product(x.elm);
  fi;
end;

###########################################################################
##
#F GoodCoxeterWord( <W>, <w> ) . . . . . . . . . . . . check if w is a good 
## . . . . . . . . element in the braid group in the sense of [Geck-Michel]
##
## checks if w^d (where d is the order of the image of w in W) is a product
## of  longest elements in a decreasing  chain of parabolic subgroups of W.
## If  this is  true, then  the list  of the  corresponding subsets  of the
## generators is returned. Otherwise, false is returned.

GoodCoxeterWord:=function(W,w)local res,i,bw,d;
  bw:=Braid(W)(w);
  d:=Order(W,EltBraid(bw));
  res:=BrieskornNormalForm(bw^d);
  if not ForAll([1..Length(res)-1],i->IsSubset(res[i],res[i+1])) then
     InfoChevie("#I ",bw," is not good\n");return false;
  fi;
  res:=Collected(res);SortBy(res,x->-Length(x[1]));
  InfoChevie("#I ",bw,"^",d,"=",Join(List(res,i->FormatMonomial(
    SPrint("w_",IntListToString(i[1],"{}")),i[2],rec())),"."),"\n");
  return res;
end;

VeryGoodCoxeterWord:=function(W,w)local res,i,bw,d;
  bw:=Braid(W)(w);
  d:=Order(W,EltBraid(bw));
  if d mod 2<>0 then return GoodCoxeterWord(W,w);fi;
  d:=d/2;
  res:=BrieskornNormalForm(bw^d);
  if not ForAll([1..Length(res)-1],i->IsSubset(res[i],res[i+1])) then
     InfoChevie("#I ",bw," is not very good\n");return false;
  fi;
  res:=Collected(res);SortBy(res,x->-Length(x[1]));
  InfoChevie("#I ",bw,"^",d,"=",Join(List(res,i->FormatMonomial(
    SPrint("w_",IntListToString(i[1],"{}")),i[2],rec())),"."),"\n");
  return res;
end;
# BipartiteDecomposition a decomposition of the generating reflections
# of W in two sets L, R such that the subgroups W_L and W_R are abelian
BipartiteDecomposition:=function(W)local L,R,comm,r,rest;
  L:=[];R:=[];rest:=Set(W.generatingReflections);
  comm:=function(x,y)x:=W.reflections[x];y:=W.reflections[y];return x*y=y*x;end;
  while Length(rest)>0 do
    r:=PositionProperty(rest,x->ForAny(L,y->not comm(x,y)));
    if r<>false then 
      if ForAny(R,y->not comm(rest[r],y)) then
        Error(W," has no bipartite decomposition");
      fi;
      Add(R,rest[r]); SubtractSet(rest,[rest[r]]);
    else r:=PositionProperty(rest,x->ForAny(R,y->not comm(x,y)));
      if r<>false then Add(L,rest[r]); SubtractSet(rest,[rest[r]]);
      else Add(L,rest[1]); SubtractSet(rest,[rest[1]]);
      fi;
    fi;
  od;
  return [W.rootInclusion{L},W.rootInclusion{R}];
end;

#############################################################################
##
##  Functions for the dual braid monoid
##
#  Maybe the best illustration is the session log below
#  
#  gap> w:=DualBraid(CoxeterGroup("A",4))(1,2,3,4,3);
#  124.34
#  gap> w*w;
#  c.[5,6,10].34.3
#  gap> w^4;
#  c.c.c.c.c
#
#############################################################################

PermRootOpsDualBraidMonoid:=function(arg)local W,new,M,n,r;
  W:=arg[1];
  M:=rec(group:=W,
    operations:=rec(Print:=function(M)Print("DualBraidMonoid(",W,",",
      M.deltaword,")");end),
    IsLeftDescending:=function(simp,i)
      return ReflectionLength(W,M.atoms[i]^-1*simp)<ReflectionLength(W,simp);
    end,
    IsRightDescending:=function(simp,i)
      return ReflectionLength(W,simp/M.atoms[i])<ReflectionLength(W,simp);
    end,
    IsRightAscending:=function(simp,i) 
      return M.IsLeftDescending(simp^-1*M.delta,i);
    end,
    IsLeftAscending:=function(simp,i) 
      return M.IsRightDescending(M.delta*simp^-1,i);
    end,
    identity:=W.identity
  );
  if not IsFinite(W) then Error("W should be finite");fi;
  if Length(arg)=1 then
    if IsCoxeterGroup(W) then 
      M.deltaword:=Concatenation(BipartiteDecomposition(W));
    elif Length(W.generatingReflections)>W.semisimpleRank or Length(W.type)>1 
    then Error("W should be irreducible and well-generated");
    else r:=PositionRegularClass(W,Maximum(ReflectionDegrees(W)));
      M.deltaword:=ChevieClassInfo(W).classtext[r];
    fi;
    M.delta:=EltWord(W,M.deltaword);
  elif IsWordFor(W,arg[2]) then
    M.deltaword:=arg[2];M.delta:=EltWord(W,M.deltaword);
  else M.delta:=arg[2];M.deltaword:=CoxeterWord(W,M.delta);
  fi;
  M.atoms:=Reflections(W);
  if not IsCoxeterGroup(W) then 
    n:=ReflectionLength(W,M.delta);
    M.atoms:=Filtered(M.atoms,r->ReflectionLength(W,M.delta/r)<n);
  fi;
  M.AtomListSimple:=function(v)
    if IsCoxeterGroup(W) then return ReflectionWord(W,v);
    else return ReflectionWord(W,v,M.atoms);
    fi;
  end;
  M.orderDelta:=Order(W,M.delta);
  M.stringDelta:="c";
  M.Reverse:=function(b)local res,s;
    if Length(b.elm)=0 then return M.revMonoid.Elt([],b.pd);fi;
    res:=[];
    for s in List(Reversed(b.elm),y->M.revMonoid.DeltaAction(y^-1,b.pd))
    do res:=M.revMonoid.AddToNormal(res,s);od;
    return GarsideEltOps.Normalize(M.revMonoid.Elt(res,b.pd));
  end;
  CompleteGarsideRecord(M,rec(interval:=true,eltOps:=rec(
    Frobenius:=function(W,x,i)local F,c; F:=Frobenius(W);c:=x.monoid.delta;
     if F(c)<>c then Error("delta should be invariant by Frobenius");fi;
     x:=ShallowCopy(x);x.elm:=List(x.elm,y->Frobenius(W)(y,i));
     return x;
    end)));
##  Add to M.B the possibility of a call B(<element of W>)  
  M.oldB:=M.B;
  M.B:=function(arg)local res,s,p;
    if Length(arg)=1 and not IsInt(arg[1]) and 
      (not IsList(arg[1]) or IsMat(arg[1])) then # assume got a group element
      arg:=ReflectionWord(W,arg[1],M.atoms);
    fi;
    return ApplyFunc(M.oldB,arg);
  end;
  if IsCoxeterGroup(W) then
  # give expression in usual monoid of simple s
  M.ToOrdinary:=function(s)local W,B,h,w;
    W:=M.group;B:=Braid(W);
    if not IsBound(M.OrdinaryAtoms) then
      W:=M.group;
      h:=OrderPerm(M.delta);
      if h mod 2=0 then w:=Concatenation(List([1..h/2],x->M.deltaword));
      else w:=Concatenation(List([1..QuoInt(h,2)],x->M.deltaword));
	Append(w,M.deltaword{[1..Length(M.deltaword)/2]});
      fi;
      if CoxeterLength(W,EltWord(W,w))<>Length(w) then 
	w:=Concatenation(List([1..h],x->M.deltaword));
      fi;
      w:=List([1..Length(w)],i->Concatenation(w{[1..i]},-w{[i-1,i-2..1]}));
      w:=Set(List(w,x->ApplyFunc(B,x)));
      SortBy(w,x->Position(Reflections(W),EltBraid(x)));
      M.OrdinaryAtoms:=w;
    fi;
    s:=M.AtomListSimple(s);
    if s=[] then return B();else return Product(M.OrdinaryAtoms{s});fi;
  end;
  fi;
  if Length(arg)=3 then M.revMonoid:=arg[3];
  else M.revMonoid:=DualBraidMonoid(W,Reversed(M.deltaword),M);
  fi;
  return M;
end;

#############################################################################
##
#F  DualBraid(W [, coxeter]) . . . . . . . returns Dual-Braid making function
##
#############################################################################
DualBraid:=function(arg)return ApplyFunc(DualBraidMonoid,arg).B;end;

#############################################################################
# Empirical reduction of the complexity of a generating set <bs> of
# a given subgroup of a Garside group
#############################################################################
ShrinkGarsideGeneratingSet:=function(bs)local b,a,ai,aib,bai,try,simplify,j,i,f;
  try:=function(p)local pos;
    p:=f(p); #   Print("p=",p,"\n");
    pos:=Position(bs,p);
    if pos<>false and pos<>j then
      InfoChevie(" eliminated");pos:=Maximum(pos,j);
      bs:=Drop(bs,pos);return true;
    fi;
    if p{[1,2]}>=bs[j]{[1,2]} then return false;fi;
    InfoChevie(" reduced to <",p[1],"-",p[2],"> ");
    bs[j]:=p;
    simplify:=true;
    return false;
  end;
  f:=function(b)local l;l:=List(AsFraction(b),x->Length(AsWord(x)));
    if l[1]>l[2] then return [l[1]+l[2],l[2],b^-1];
                 else return [l[1]+l[2],l[1],b];
    fi;
  end;
  bs:=List(bs,f);
  repeat
    simplify:=false;
    bs:=Set(bs);
    InfoChevie("#I total length ",Sum(bs,x->x[1]),
               "   maximal length ",bs[Length(bs)][1],"\n");
    for j in [Length(bs),Length(bs)-1..2] do
      InfoChevie("#I <",j,":",bs[j][1],"-",bs[j][2],">");
      i:=0;
      repeat i:=i+1;InfoChevie(".\c");
        a:=bs[i][3];ai:=a^-1;b:=bs[j][3];aib:=ai*b;bai:=b*ai;
#       Print("i=",i,"bs[i]=",bs[i]," bs[j]=",bs[j],"\n");
      until bs[i]{[1,2]}>=bs[j]{[1,2]} or try(aib) or try(bai) or 
        try(aib*a) or try(a*bai) or try(a*b) or try(b*a);# or try(aiba/b);
      InfoChevie("\n");
    od;
  until not simplify;
  return List(bs,x->x[3]);
end;

#############################################################################
# returns a presentation of the Garside group corresponding to the monoid M 
# (by theorem 4.1 of Dehornoy-Paris)
Presentation:=function(M)local F,rels,i,j,v;
  F:=FreeGroup(Length(M.atoms));
  rels:=[];
  for i in [1..Length(M.atoms)] do
    for j in [1..i-1] do
      v:=M.LeftLcmSimples(M.atoms[i],M.atoms[j]);
      v:=List(v{[2,3]},M.AtomListSimple);
      Add(v[1],i);Add(v[2],j);
      Add(rels,Product(v[1],i->F.(i))*Product(v[2],i->F.(i))^-1);
    od;
  od;
  return PresentationFpGroup(F/rels);
end;

TwistedPowerMonoid:=function(M1,n)local M,m,a,rot,lrot,atom,atomops,dec;
  rot:=function(i,t)if t then return 1+(i mod n);else return i;fi;end;
  lrot:=function(i,t)if t then return 1+((i-2)mod n);else return i;fi;end;
  m:=M1.nrAtoms;
  dec:=function(i)i:=i-2;return 1+[QuoInt(i,m),i mod m];end;
  atomops:=rec(Print:=function(r)
      if r.t then Print("t");fi;
      Print("(",Join(List(r.v,a->IntListToString(M1.AtomListSimple(a)))),")");
      end
  );
  M:=rec(operations:=rec(
    Print:=function(M)Print("twisted ",Ordinal(n)," power of ",M1);end));
  M.atom:=function(v,t)return rec(v:=v,t:=t,operations:=atomops);end;
  M.identity:=M.atom(List([1..n],x->M1.identity),false);
  M.atoms:=[M.atom(M.identity.v,true)];
  Append(M.atoms,Concatenation(List([1..n],i->List(M1.atoms,function(a)local r;
      r:=ShallowCopy(M.identity);r.v:=ShallowCopy(r.v);
      r.v[i]:=a;return r;end))));
  M.IsLeftDescending:=function(s,i)if i=1 then return s.t;fi;
    i:=dec(i); return M1.IsLeftDescending(s.v[lrot(i[1],s.t)],i[2]); end;
  M.IsRightDescending:=function(s,i)if i=1 then return s.t;fi;
    i:=dec(i); return M1.IsRightDescending(s.v[i[1]],i[2]);end;
  M.IsRightAscending:=function(s,i)if i=1 then return not s.t;fi;
    i:=dec(i); return M1.IsRightAscending(s.v[i[1]],i[2]);end;
  M.delta:=M.atom(List([1..n],i->M1.delta),true);
  M.Product:=function(a,b)if IsList(b) then return List(b,x->M.Product(a,x));fi;
    return M.atom(List([1..n],i->M1.Product(a.v[rot(i,b.t)],b.v[i])),a.t or b.t);
  end;
  M.LeftQuotient:=function(a,b)
    if IsList(b) then return List(b,x->M.LeftQuotient(a,x));fi;
    return M.atom(List([1..n],i->M1.LeftQuotient(a.v[rot(i,b.t<>a.t)],b.v[i])),
                b.t<>a.t);
  end;
  M.RightQuotient:=function(a,b)
    if IsList(a) then return List(a,x->M.RightQuotient(x,b));fi;
    return M.atom(List([1..n],
      i->M1.RightQuotient(a.v[lrot(i,b.t)],b.v[lrot(i,b.t)])),a.t<>b.t);
  end;
  CompleteGarsideRecord(M,rec());
  return M;
end;
