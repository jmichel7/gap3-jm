###########################################################################
#
# affa.g                             (C) François Digne (Amiens university)
#
# This file contains:
#   -functions for computing with periodic permutations of the integers
#     The main functions are:
#     PPerm to create a periodic permutation.
#     Operations *, /, ^
#     functions Cycles and CycleType
#
#   -the function CoxeterGroupAtildeGroup which creates the Coxeter group
#    of type ~A_n as a group of periodic permutations.
#
#   -the function DualBraidMonoid for the above kind of groups.
# It is based on the papers
# [Digne] Presentations duales pour les groupes de tresses de type affine A
#         Comment. Math. Helv. 81 (2006) 23--47
# [Shi] The Kazhdan-Lusztig cells in certain affine Weyl groups 
#       Springer LNM 1179 (1986) 
#--------------------------------------------------------------------------
# Functions for computing with periodic permutations f of the integers
# - of period n, which means f(i+n)=f(i)+n
# - with no shift, which means Sum([1..n],f)=Sum([1..n])
# Such permutations are stored as 
#    rec(perm:=[f(1),...,f(n)],operations:=PPermOps)
CycleType:=Dispatcher("CycleType");

PPermOps:=OperationsRecord("PPermOps"); # operations for periodic permutations
PPermOps.PrintOptions:=rec();

PPermOps.New:=perm->rec(perm:=perm,operations:=PPermOps);

# Function to define periodic permutations of period n. Forms:
# PPerm(f(1),...,f(n))
# PPerm([f(1),...,f(n)])
# PPerm(cycle_1,...,cycle_m,n) returns the product of the permutations
# represented by cycle_1*..*cycle_m where cycle_i is a list [i_1,..,i_k,[d]] 
# where i_1,..,i_k differ mod n representing the permutation 
# i_1->i_2->..->i-k->i_1+d*n         if d=0 then [d] can be omitted
PPerm:=function(arg)local v,res,n;
  if Length(arg)=0 then return PPerm([]);
  elif IsList(arg[1]) and IsInt(arg[Length(arg)]) then 
    n:=arg[Length(arg)];
    res:=Product(arg{[1..Length(arg)-1]},function(cyc)local i,perm,k,d;
      if not IsList(cyc[Length(cyc)]) then d:=0;
      else d:=cyc[Length(cyc)][1];cyc:=cyc{[1..Length(cyc)-1]};fi;
      if Length(Set(List(cyc,x->x mod n)))<>Length(cyc) then
        Error(cyc," : the images must be distinct mod ",n);
      fi;
      perm:=[1..n];cyc:=cyc-1;
      for i in [1..Length(cyc)] do
	k:=1+cyc[i] mod n;
	perm[k]:=cyc[1+i mod Length(cyc)]+k-cyc[i];
      od;  
      perm[k]:=perm[k]+d*n;
      return PPermOps.New(perm);
    end);
  else 
    if Length(arg)=1 and IsList(arg[1]) then v:=arg[1];else v:=arg;fi;
    n:=Length(v);
    res:=PPermOps.New(v);
  fi;  
  if Set(List(res.perm,x->x mod n))<>[0..n-1] then
    Error("images must be distinct mod ",n);
  fi;
  if Sum(res.perm)<>Sum([1..n]) then Error("sum of shifts must be 0");fi;
  return res;
end;

PPermOps.\*:=function(x,y)local ymodn,n;
  if IsList(x) then return List(x,z->PPermOps.\*(z,y));
  elif IsList(y) then return List(y,z->PPermOps.\*(x,z));
  fi;
  n:=Length(x.perm);ymodn:=1+List(y.perm-1,p->p mod n);
  return PPermOps.New(x.perm{ymodn}+y.perm-ymodn);
end;

PPermOps.\^:=function(x,i)local y,res,n,l,ll;
  if IsInt(x) then 
    y:=1+(x-1)mod Length(i.perm);return x+i.perm[y]-y;
  fi;
  if not IsInt(i) then return i^-1*x*i;fi;
  if i=0 then return PPermOps.New([1..Length(x.perm)]);
  elif i=1 then return x;
  elif i>0 then res:=x^0;y:=x;
    while i>0 do
     if i mod 2<>0 then res:=res*y;fi;
     y:=y*y;
     i:=QuoInt(i,2);
    od;
    return res;
  elif i=-1 then n:=Length(x.perm);
    l:=List(x.perm-1,i->i-(i mod n));
    ll:=OnTuples([1..n],PermList(x.perm-l)^-1);
    return PPerm(ll-l{ll});
  else return (x^-1)^-i;
  fi;
end;

PPermOps.\/:=function(x,y)return x*y^-1;end;

#PPermOps.\=:=function(a,b)return a.perm=b.perm;end;

# Non-trivial cycles of a PPerm; each cycle i_1,..,i_k,[d] is normalized
# such that i_1 mod n is the smallest of i_j mod n and i_1 is in [1..n]
PPermOps.Cycles:=function(a)local l,x,cyc,n,res,j;
  res:=[];n:=Length(a.perm);l:=[1..n];
  while Length(l)>0 do
    x:=l[1];cyc:=[];
    repeat Add(cyc,x);x:=x-1;j:=x mod n;x:=a.perm[1+j]+x-j;
    until (x-cyc[1])mod n=0;
    SubtractSet(l,1+List(cyc-1,x->x mod n));
    Add(cyc,[(x-cyc[1])/n]);
    if Length(cyc)>2 or cyc[2]<>[0] then Add(res,cyc);fi;
  od;
  return res;
end;    

PPermOps.Format:=function(a,opt)local n,stringdec,c; n:=Length(a.perm);
  stringdec:=function(d)
    if IsBound(opt.sgn) then
      if d>0 then return List([1..d],x->'+');
      else return List([1..-d],x->'-');
      fi;
    else
      if d=0 then return "";else return SPrint("[",d,"]");fi;
    fi;
  end;
  if IsBound(opt.GAP) then 
    return SPrint("PPerm(",FormatGAP(a.perm),")");
  fi;
  c:=PPermOps.Cycles(a);
  if c=[] then return "()";fi;
  return Concatenation(List(c,function(cyc)local res,d;
    d:=cyc[Length(cyc)][1];cyc:=cyc{[1..Length(cyc)-1]};
    res:=SPrint("(",Join(List(cyc,function(y)local x;
     y:=y-1;x:=y mod n;return SPrint(1+x,stringdec((y-x)/n));end)),")");
     PrintToString(res,stringdec(d));
     return res;end));
end;    

PPermOps.String:=a->PPermOps.Format(a,PPermOps.PrintOptions);

PPermOps.Print:=function(a)Print(String(a));end;

#The following formula is from [Shi] Lemma 4.2.2
PPermOps.CoxeterLength:=function(w)local n;
  n:=Length(w.perm);
  return Sum([1..n],j->Sum([1..j-1],i->AbsInt(Floor((j^w-i^w)/n))));
end;

PPermOps.IsRightDescending:=function(w,i)
  if i=Length(w.perm) then return w.perm[i]>w.perm[1]+Length(w.perm);
  else return w.perm[i]>w.perm[1+i];fi;
end;

PPermOps.IsLeftDescending:=function(w,i)
  return w.operations.IsRightDescending(w^-1,i);
end;

PPermOps.FirstLeftDescending:=function(x)local i,IRD;
  x:=x^-1;IRD:=x.operations.IsRightDescending;
  for i in [1..Length(x.perm)] do if IRD(x,i) then return i;fi;od;
  return false;
end;

# for this function see [Digne],2.8
PPermOps.ReflectionLength:=function(w)local res,d,n,p0,pos,neg,m,vec;
  vec:=function(pp,i)
    pp:=List(PartitionsSet([1..Length(pp)],i),x->List(x,y->Sum(pp{y})));
    if pp[1][1]<0 then pp:=-pp;fi;
    return Set(List(pp,Collected));
  end;
  d:=List(Cycles(w),x->x[Length(x)][1]);n:=Length(w.perm);
  p0:=Number(d,x->x=0);
  res:=n+Length(d)-2*p0-Number([1..n],i->i=i^w);
  if p0=Length(d) then return res;fi;
  pos:=Filtered(d,x->x>0);neg:=Filtered(d,x->x<0);
  m:=Minimum(Length(pos),Length(neg));
  return res-2*First([m,m-1..1],
     i->Length(Intersection(vec(pos,i),vec(neg,i)))>0);
end;

PPermOps.ReflectionWord:=function(w)local n,ff,res,s;
  n:=Length(w.perm);
  ff:=function(w)local l,d,i,s;
    l:=ReflectionLength(w); d:=1;
    while true do
      for i in [1..n] do
        s:=PPerm([i,i+d],n);
        if ReflectionLength(s*w)<l then return s; fi;
      od;
      d:=d+1;
      if d mod n=0 then d:=d+1;fi;
    od;
  end;
  res:=[];
  while w<>w^0 do
    s:=ff(w);Add(res,s);w:=s*w;
  od;
  return res;
end;

PPermOps.CycleType:=function(a)local res;
  res:=List(PPermOps.Cycles(a),cyc->[Length(cyc)-1,cyc[Length(cyc)]]);
  Sort(res);
  return res;
end;
#--------------------------------------------------------------------------
# The next function constructs W(~A_{n-1}) as a group of periodic
# permutations of period n.
 
AtildeGroupOps:=OperationsRecord("AtildeGroupOps",GroupOps);

AtildeGroupOps.Print:=function(W)Print(W.name);end;

AtildeGroupOps.PrintDiagram:=function(W)
    CHEVIE.R("PrintDiagram","AffineA")(W.reflectionsLabels);end;
 
AtildeGroupOps.IsRightDescending:=function(W,w,i)
  return PPermOps.IsRightDescending(w,i);
end;

AtildeGroupOps.FirstLeftDescending:=function(W,x)
  return PPermOps.FirstLeftDescending(x);
end;

AtildeGroupOps.IsLeftDescending:=function(W,w,i)
  return PPermOps.IsLeftDescending(w,i);
end;

# Reflections (a,b[i]) are enumerated by lexicographical order of [i,a,b-a]
# with i positive --- recall that when a>b this reflection is printed (b,a[-i])
AtildeGroupOps.Reflection:=function(W,i) local n,p,r,ecart,pos;
  n:=W.semisimpleRank;
  p:=QuoInt(i-1,n*(n-1)); r:=(i-1) mod (n*(n-1));
  ecart:= QuoInt(r,n); pos:=r mod n;
  return PPerm([1+pos,2+pos+ecart+p*n],n);
end;

AtildeGroupOps.CoxeterLength:=function(W,w)
  return PPermOps.CoxeterLength(w);
end;

AtildeGroupOps.ReflectionLength:=function(W,w)
  return PPermOps.ReflectionLength(w);
end;
  
CoxeterGroupAtildeGroup:=function(n)local W,i,refs;
  W:=rec(isGroup:=true, isDomain:=true,
    name:=SPrint("CoxeterGroupAtildeGroup(",n,")"),
    operations:=AtildeGroupOps);
  W.generators:=List([1..n-1],i->PPermOps.New(Permuted([1..n],(i,i+1))));
  W.generators[n]:=PPermOps.New(Concatenation([0],[2..n-1],[n+1]));
  for i in [1..n] do W.(i):=W.generators[i];od;
  W.reflections:=W.generators;
  W.operations.\in:=function(e,W)return Length(e.perm)=W.rank;end;
  AbsCoxOps.CompleteCoxeterGroupRecord(W);
  return W;
end;

# DualMonoid(W[,M])
# constructs dual monoid for tilde A_{n-1}
#
# If a second argument is given, constructs the reversed monoid
# [which allows to fill the field M.revMonoid after building the dual monoid]

AtildeGroupOps.DualBraidMonoid:=function(arg)local M,n,W,delta;
  W:=arg[1];n:=W.rank;
  if Length(arg)=1 then
       delta:=rec(perm:=Concatenation([1-n],[3..n],[2+n]),operations:=PPermOps);
  else delta:=arg[2].delta^-1;
  fi;
  M:=rec(rank:=n,
    delta:=delta,
    stringDelta:="c",
    group:=W,
    identity:=W.identity,
    FormatSimple:=function(a,opt)local option;
      option:=ShallowCopy(PPermOps.PrintOptions);Inherit(option,opt);
      if IsBound(option.word) then return IntListToString(CoxeterWord(W,a));
      else return Format(a,option);fi;end,
    Reverse:=function(b)local res,s;
      if Length(b.elm)=0 then return M.revMonoid.Elt([],b.pd);fi;
      res:=[];
      for s in List(Reversed(b.elm),y->M.revMonoid.DeltaAction(y^-1,b.pd))
      do res:=M.revMonoid.AddToNormal(res,s);od;
      return GarsideEltOps.Normalize(M.revMonoid.Elt(res,b.pd));
    end,
    RightComplementToDelta:=a->a^-1*M.delta,
# descent sets are encoded as a pair: a list of atoms, and a list
# of atoms of the form [1,u] representing all atoms [1,u+in]
# This uses lemma 2.20 of [Digne] and is valid only if there are
# 0 or 2 cycles with a non-zero shift
    LeftDescentSet:=function(a)local d,j,k,x,res;
      res:=[[],[]];
      for x in PPermOps.Cycles(a) do
        d:=x[Length(x)][1];x:=x{[1..Length(x)-1]};
	if x[1]<>1 or Length(x)<>1 then
	  for j in [1..Length(x)] do
	    for k in [j+1..Length(x)] do 
	      Add(res[1],[x[j],x[k]]);
	      if d<>0 then Add(res[1],[x[j]+n,x[k]]);fi;  
	    od;
	    if d<>0 then Add(res[2],[1,1+((x[j]-1) mod n)]);fi;  
	  od;
	fi;
      od;
      return List(res,i->List(i,x->PPerm(x,n)));
    end,
    RightAscentSet:=a->M.LeftDescentSet(M.RightComplementToDelta(a)),
    operations:=rec(Print:=function(M) Print("DualBraidMonoid(",W,")");end)
  );
  M.DeltaAction:=function(s,i)return s^(M.delta^i);end;
  M.FirstIntersectionLDS:=function(a,b)local t;
    for t in a[1] do 
      if t in b[1] or 
        List(t.perm,u->1+((u-1) mod n)) in List(b[2],x->x.perm) then return t;
      fi;
    od;
    for t in b[1] do 
      if List(t.perm,u->1+((u-1) mod n)) in List(a[2],x->x.perm) then return t;
      fi;
    od;
    for t in a[2] do if t in b[2] then return t;fi; od;
    return false;
  end;
  M.LeftGcdSimples:=function(a,b)local t,x;
    x:=M.identity;
    while true do
      t:=M.FirstIntersectionLDS(M.LeftDescentSet(a),M.LeftDescentSet(b));
      if t=false then return [x,a,b];fi;
      x:=x*t;a:=t^-1*a;b:=t^-1*b;
    od;
  end;
  M.Elt:=function(arg)local res;
    res:=rec(elm:=arg[1],operations:=GarsideEltOps,monoid:=M);
    if Length(arg)>1 then res.pd:=arg[2]; else res.pd:=0; fi;
    return GarsideEltOps.Normalize(res);
  end;
  M.B:=function(arg)local x,p;
    if IsList(arg[1]) then x:=[arg[1]]; else x:=arg; fi;
    p:=PositionProperty(x,y->not Number(Cycles(y),c->c[Length(c)][1]<>0)in [0,2]);
    if p<>false then Error(x[p]," is not a dual simple");fi;
    x:=Concatenation(List(x,M.AtomListSimple));
    if not ForAll(x,M.IsDualAtom) then
      Error("not atom of dual monoid: ",First(x,y->not M.IsDualAtom(y)));fi;
    return GarsideEltOps.Normalize(Product(List(x,y->M.Elt([y]))));
  end;
  M.AtomListSimple:=function(w) local s,res,v;res:=[];v:=w;
    while v<>M.group.identity do
     s:=Flat(M.LeftDescentSet(v))[1];
     v:=s^-1*v;
     Add(res,s);
    od;
    return res;
  end;
  M.IsDualAtom:=function(a)local c;
   c:=Cycles(a);
   return Length(c)=1 and Length(c[1])=3 and c[1][3]=[0] and
    (c[1][1]=1 or AbsInt(c[1][1]-c[1][2])<M.rank);
  end;
  CompleteGarsideRecord(M,rec(interval:=true));
  if Length(arg)=1 then M.revMonoid:=DualBraidMonoid(W,M);
  else M.revMonoid:=arg[2];
  fi;
  return M;
end;

AtildeBraid:=n->DualBraidMonoid(CoxeterGroupAtildeGroup(n)).B;
Atilde:=PPerm;
