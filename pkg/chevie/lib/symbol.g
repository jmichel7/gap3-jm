#############################################################################
##
#A  symbol.g                CHEVIE library       Meinolf Geck, Frank Luebeck, 
#A                                                Jean Michel, Goetz Pfeiffer
##
#Y  Copyright (C) 1992 -  2011 Lehrstuhl D fuer  Mathematik, RWTH Aachen
#Y  and the universities of Kaiserslautern, Galway, and Paris VII.
##
##  This file contains a few functions dealing with combinatorics of
##  multiple partitions and symbols.
## 
##  A basic reference is G. Malle, "Unipotente Grade...", J.Algebra 177 (1995)
##  
#############################################################################

#############################################################################
#F  PartitionTupleToString( <mpart>[,opt] ) . . . . . . . . . as the name says
##  
##  'IntListToString' is applied to the  partitions in <mpart> and
##  the results are concatenated by a separating '.':
##  
##     [[],[2,1,1]]      --->    ".211"
##     [[4,3],[2,1,1]]   --->    "43.211"
##     [[14,3],[2,1,1]]  --->    "(14,3).211"
##  
##  Knows how to handle partitions rotationally invariant of G(de,e,r)
##  
PartitionTupleToString:=function(arg)local s,v,opt;s:=arg[1];
  if Length(arg)=1 then opt:=rec();else opt:=arg[2];fi;
  if IsInt(s[Length(s)]) then v:=E(s[Length(s)-1])^s[Length(s)];
    if v=1 then v:="+";elif v=-1 then v:="-";else v:=Format(v,opt);fi;
    s:=Join(List(s{[1..Length(s)-2]},IntListToString),".");
    if IsBound(opt.TeX) then Append(s,"\\cdot ");fi;
    return SPrint(s,v);
  else return String(Join(List(s,IntListToString),"."));
  fi;
end;

# same as above but the result is separated by ',' and bracketed
StringSymbol:=function(s)local v,e;
  e:=Length(s);
  if e=0 then return "()";
  elif IsInt(s[e]) then v:=E(s[e-1])^s[e];
    if v=1 then v:="+";elif v=-1 then v:="-";else v:=Format(v);fi;
    return SPrint("(",Join(List(s{[1..e-2]},IntListToString)),v,")");
  else return SPrint("(",Join(List(s,IntListToString)),")");
  fi;
end;

# shift a list of beta-numbers
ShiftBeta:=function(beta,n)
  if n>=0 then return Concatenation([0..n-1],beta+n);
  elif beta{[1..-n]}<>[0..-n-1] then Error("Cannot shift ",beta," by ",n,"\n");
  else return beta{[1-n..Length(beta)]}+n;
  fi;
end;

# gives the partition corresponding to a list b of beta-numbers
PartBeta:=function(b)
  if Length(b)=0 then return []; # needed since []+[] is an error!!!
  else return Reversed(Filtered(b-[0..Length(b)-1],x->x<>0));
  fi;
end;

# SymbolPartitionTuple(p,s) symbol of shape s associated to multipartition p.
# Special cases: 
#  -  when s  is a  positive integer  it is  interpreted as [s,0,0,...] and
#  negative  as [0,-s,-s,....] so when p is a double partition one gets the
#  symbol of 'defect' s associated to p;
#  also the principal series of G(e,1,r) is SymbolPartitionTuple(p,1) and
#  that of G(e,e,r) is  SymbolPartitionTuple(p,0);
#  - for parameters for G(e,e,r) invariant by some cyclic shift thus
#  truncated, works provided s=0.
SymbolPartitionTuple:=function(p,s)local e,l;
  if IsInt(p[Length(p)]) then l:=Length(p)-2;e:=l*p[Length(p)-1];
  else e:=Length(p);l:=Length(p);
  fi;
  if IsInt(s) then if s<0 then s:=Concatenation([0],[1..l-1]*0-s);
                   else s:=Concatenation([s],[1..l-1]*0);fi;
  else s:=s{[1..l]};fi;
  s:=List(p{[1..l]},Length)-s; s:=Maximum(s)-s;
  p:=Copy(p);p{[1..l]}:=List([1..l],i->ShiftBeta(BetaSet(p[i]),s[i]));
  return p;
end;

FullSymbol:=function(S)local m;
  if Length(S)=0 or IsList(S[Length(S)]) then return S;
  else m:=S[Length(S)-1]; return Concatenation(List([1..m],
     i->List(S{[1..Length(S)-2 ]},ShallowCopy)));
  fi;
end;

# defect of a symbol with 2 lines
DefectSymbol:=function(s)s:=FullSymbol(s);return Length(s[1])-Length(s[2]);end;

# rank of a symbol
RankSymbol:=function(s)local ss,e;s:=FullSymbol(s);
  ss:=Sum(s,Length);e:=Length(s);
  if e=0 then return 0;fi;
  return Sum(s,Sum)-QuoInt((ss-1)*(ss-e+1),2*e);
end;

# returns all the 2-line symbols of defect d and rank n
Symbols:=function(n,d)
  n:=n-QuoInt(d*d,4);
  if n<0 then return [];fi;
  if d>0 then return List(PartitionTuples(n,2),x->SymbolPartitionTuple(x,d));fi;
  return List(CHEVIE.R("CharInfo","imp")(2,2,n).charparams,
	      CHEVIE.R("symbolcharparam","D"));
end;

# Returns  the  fake  degree  of  the  unipotent character corresponding to
# symbol s as a 'CycPol' (see cycpol.g).
# Works for symbols of G(e,1,r) (d=1) and G(e,e,r) (d=0).
# When given a second argument p dividing e, works for the coset G(e,e,r).s_1^p
# See Malle, "Unipotente Grade", 2.11 and 5.7
#
CycPolFakeDegreeSymbol:=function(arg)local s,p,res,delta,theta,r,d,e,q,rot;
  s:=FullSymbol(arg[1]); q:=X(Cyclotomics); e:=Length(s); r:=RankSymbol(s);
  if e=0 then return CycPol(1);fi;
  if Length(arg)=2 then p:=E(e)^arg[2];else p:=1;fi;

  delta:=S->Product(Combinations(S,2),x->CycPol(q^(e*x[2])-q^(e*x[1])));
  theta:=S->Product(Filtered(S,x->x>0),l->Product([1..l],h->CycPol(q^(e*h)-1)));

  if Length(s)=1 then d:=0; else d:=DefectSymbol(s); fi;

  if   d=1 then res:=theta([r]);
  elif d=0 then res:=theta([r-1])*CycPol(q^r-p);
  else res:=CycPol(0*q); # non-principal series
  fi;

  res:=res*Product(s,S->delta(S)/theta(S))/
    CycPol(q^Sum(e*[1..QuoInt(Sum(s,Length),e)-1]+(d mod 2),x->x*(x-1)/2));

  if d=1 then res:=res*CycPol(q^([0..e-1]*List(s,Sum))); 
  else rot:=Rotations(s);
    res:=res*CycPol(List([0..e-1],j->p^j)*List(rot,s->q^([0..e-1]*List(s,Sum))))
           /Number(rot,i->i=s);
    if e=2 and p=-1 then res:=-res;fi; # fix type 2D
  fi;

  if r=2 and e>2 and p=E(e) then 
    res:=CycPol(Value(res,E(2*e)*q)/E(2*e)^Degree(res));
    # not CycPolOps.EnnolaTwist since not a cycpol!
  fi;
  return res;
end;

# works for same cases as CycPolFakeDegreeSymbol
LowestPowerFakeDegreeSymbol:=function(s)local res,gamma,d,e;
  s:=FullSymbol(s);
  if Length(s)=1 then d:=0; else d:=DefectSymbol(s); fi;
  if not d in [0,1] then return -1;fi;
  e:=Length(s);
  res:=e*Sum(Filtered(s,x->x<>[]),S->S*[Length(S)-1,Length(S)-2..0]);
  gamma:=i->List(i+[0..e-1],x->x mod e)*List(s,Sum);
  if d=1 then res:=res+gamma(0);
  else res:=res+Minimum(List([0..e-1],gamma));
  fi;
  return res-Sum(e*[1..QuoInt(Sum(s,Length),e)-1]+(d mod 2),x->x*(x-1)/2);
end;

# works for same cases as CycPolFakeDegreeSymbol
HighestPowerFakeDegreeSymbol:=function(arg)local s,p,res,gamma,r,d,e;
  s:=FullSymbol(arg[1]); e:=Length(s);
  if Length(s)=1 then d:=0; else d:=DefectSymbol(s); fi;
  if Length(arg)=2 then p:=arg[2];else p:=e;fi;
  if not d in [0,1] then return -1;fi;
  r:=RankSymbol(s); 
  if d=1 then res:=e*r*(r+1)/2;
  else        res:=e*r*(r-1)/2+r;
  fi;
  res:=res+e*Sum(Filtered(s,x->x<>[]),S->S*[0..Length(S)-1]-Sum(S,l->l*(l+1)/2));
  gamma:=i->List(i+[0..e-1],x->x mod e)*List(s,Sum);
  if d=1 then res:=res+gamma(0);
  else        res:=res+Maximum(List([0..e-1],gamma));
  fi;
  return res-Sum(e*[1..QuoInt(Sum(s,Length),e)-1]+(d mod 2),x->x*(x-1)/2);
end;

# Generic degree of the unipotent character corresponding to symbol S
# as a 'CycPol' (see cycpol.g). Works for symbols for
# G(e,1,r) (d=1, defect=0), for G(e,e,r) (d=0, defect=0) 
# and for {}^2 G(e,e,2) (d=0, defect=1) (this includes 2Dn, 2B2, 2G2)
# here d=Inhalt mod. e
# See Malle's "Unipotente Grade", 3.9 and 6.4
CycPolGenericDegreeSymbol:=function(S)local sh,res,theta,m,r,d,e,defect,q;
  S:=FullSymbol(S);r:=RankSymbol(S); e:=Length(S); sh:=List(S,Length); 
  if e=0 then return CycPol(1);fi;
  m:=QuoInt(Sum(sh),e); d:=Sum(sh)mod e; q:=X(Cyclotomics);
  defect:=(Binomial(e,2)*m-sh*[0..e-1]) mod e; 
  theta:=S->Product(S,l->Product([1..l]*e,h->CycPol(q^h-1)));

  # initialize with the q'-part of the group order
  if d=1 then res:=theta([r]);
  elif d=0 then res:=theta([r-1])*CycPol(q^r-E(e)^defect);
  fi;

  res:=res*(-1)^([0..e-1]*List(sh,x->Binomial(x,2)))*
   Product([0..e-1],i->Product([i..e-1],j->Product(List(S[i+1],l->q^l*E(e)^i),
    l->Product(Filtered(List(S[j+1],l->q^l*E(e)^j),
    m->i<j or Degree(m)<Degree(l)),m->CycPol(l-m)))))/
     (Product(S,theta)*(E(4)^Binomial(e-1,2)*ER(e)^e)^m
     *CycPol(q^Sum(e*[1..m-1]+d,x->Binomial(x,2))));
  
  if r=1 then m:=Position(S,[]);if m<>false then res:=(-1)^m*res;fi;
  elif r=2 and e=3 and S in 
    [[[1],[0,1,2],[0,1,2]],[[],[0,2],[0,1]],[[],[0,1],[0,2]]] then res:=-res;
  elif r=3 and e=3 and d=0 and S in 
    [[[0,1,2],[],[]],[[0,1,2],[0,1,2],[]]] then res:=-res;
  elif r=4 and e=3 and d=0 and S in 
    [[[0,1,3],[],[]],[[0,1,2,3],[1],[0]],[[0,1,2,3],[0],[1]],
    [[0,1,3],[0,1,2],[]],[[0,1,2],[0,1,3],[]],[[0,1,2,3],[0,1,2,3],[1]]]
    then res:=-res;
  elif r=5 and e=3 and d=0 and S in 
  [[[0,2,3],[],[]],[[0,1,2,4],[1],[0]],[[0,1,2,4],[0],[1]],[[0,1,2,3,4],[1,2],
  [0,1]],[[0,1,2,3],[1],[1]],[[0,1,2,3,4],[0,1],[1,2]],[[0,1,4],[],[]],[[0,1,2,
  3],[2],[0]],[[0,1,2,3],[0],[2]],[[0,2,3],[0,1,2],[]],[[0,1,3],[0,1,3],[]],[[0,
  1,2,4],[0,1,2,3],[1]],[[0,1,2],[0,2,3],[]],[[0,1,2,3],[0,1,2,4],[1]],[[0,1,2,
  3,4],[0,1,2,3,4],[1,2]],[[0,1,4],[0,1,2],[]],[[0,1,2],[0,1,4],[]],[[0,1,2,3],
  [0,1,2,3],[2]]] then res:=-res;
  fi;
  # Dudas' sign change
  
  if d=0 then res:=res*First([1..e],i->S{1+List(i+[0..e-1],j->j mod e)}=S);fi;
  if defect=0 or r<>2 or e<=2 then return res;
  else return E(e)^-1*CycPolOps.EnnolaTwist(res,E(2*e)); # 2I(e)
  fi;
end;

# valuation of the generic degree of the unipotent character of symbol p 
# Works for symbols of 2Dn, and for symbols for G(e,1,r) and G(e,e,r)
LowestPowerGenericDegreeSymbol:=function(p)local m,e;
  p:=FullSymbol(p);e:=Length(p);p:=Concatenation(p); Sort(p); m:=Length(p);
  return p*[m-1,m-2..0]-QuoInt(m*(m-e)*(2*m-3-e),12*e);
end;

# degree of the generic degree of the unipotent character of symbol p 
# Works for symbols of 2Dn, and for symbols for G(e,1,r) and G(e,e,r)
HighestPowerGenericDegreeSymbol:=function(p)local m,r,e;
  p:=FullSymbol(p);r:=RankSymbol(p);e:=Length(p);
  p:=Concatenation(p); Sort(p); m:=Length(p);
  if (m mod e)=1 then r:=e*r*(r+1)/2;
  else r:=e*r*(r-1)/2+r;
  fi;
  return r+p*[0..m-1]-Sum(p,x->e*x*(x+1)/2)-QuoInt(m*(m-e)*(2*m-3-e),12*e);
end;

# list of standard tableaux of shape the partition-tuple S (that is, a filling
# of S with the numbers [1..Sum(S,Sum)] such that the numbers increase across
# the rows and down the columns).
# If S is a single partition return single tableaux.
Tableaux:=function(S)local res,w,single;
  if S=[] then return [];fi;
  single:=not IsList(S[1]);if single then S:=[S];fi;
  w:=Sum(S,Sum);
  if w=0 then return [List(S,x->List(x,y->[]))];fi;
  res:=Concatenation(List([1..Length(S)],function(i)local rim,l;
    l:=Length(S[i]);rim:=Filtered([1..l-1],j->S[i][j+1]<S[i][j]);
    if l<>0 and S[i][l]<>0 then Add(rim,l);fi;
    return Concatenation(List(rim,function(p)local n,t;
      n:=Copy(S);n[i][p]:=n[i][p]-1;n:=Tableaux(n);
      for t in n do Add(t[i][p],w);od;
      return n;
    end));end));
  if single then return List(res,x->x[1]);else return res;fi;
end;

Compositions:=function(arg)local n,k;
  n:=arg[1];
  if Length(arg)=1 then if n=0 then return [[]];fi;
    return Concatenation(List([1..n],i->List(Compositions(n-i),
      c->Concatenation(c,[i]))));
  fi;
  k:=arg[2];if k=1 then return [[n]];fi;
    return Concatenation(List([1..n-1],i->List(Compositions(n-i,k-1),
      c->Concatenation(c,[i]))));
end;

#############################################################################
##
#F  DifferencePartitions( <gamma>, <alpha> ) . . .  difference of partitions.
##
##  If <gamma>-<alpha> doesn't exist or  contains a 2x2 box then 'false'
##  is returned.  Otherwise this  difference is a  union of  skew hooks.
##  The  function returns  a  record with  fields  cc:=number of  hooks,
##  ll:=combined leg length and d:=the  contents of a box underneath the
##  last hook.
##  This function is used for the computation of single character values
##  for Hecke algebras of types A, B or D.
##  
DifferencePartitions:=function(gamma, alpha)local i, dp, old, new, int, inhook;
   dp:= rec(cc:= 0, ll:= 0);
   if Length(alpha) > Length(gamma) then return false; fi;
   old:= []; inhook:= false;  
   alpha:=Concatenation(alpha,[Length(alpha)+1..Length(gamma)]*0);
   for i in [1..Length(gamma)] do
      if alpha[i] > gamma[i] then return false; fi;
      new:= [alpha[i]+1..gamma[i]];
      int:= Intersection(old, new);
      if Length(int) > 1 then return false;
      elif Length(int) = 1 then dp.ll:= dp.ll + 1;
      elif inhook then dp.cc:= dp.cc + 1; dp.d:= old[1] - i; inhook:= false;
      fi;
      if new <> [] then inhook:= true; fi;
      old:= new;
   od;

   if inhook then dp.cc:= dp.cc + 1; dp.d:= old[1] - Length(gamma) - 1; fi;
   return dp;
end;

#############################################################################
##
#F  LessSymbols( <x>, <y> ) . . . . . . . . the '<' function for symbols
##  
##  A symbol is smaller than another if the shape is lexicographically smaller,
##  or the shape is the same and the the symbol is lexicographically bigger
##                                               
LessSymbols:=function(x,y)return 
  List(x,Length)<List(y,Length) or (List(x,Length)=List(y,Length) and x>y);
end;

# e-symbols of rank r, Malle-defect def and content=c (mod e)
# The list is returned sorted by HC series (principal series first)
# SymbolsDefect(d,r,0,1) gives symbols of unipotent characters of G(d,1,r)
# SymbolsDefect(e,r,0,0) gives symbols of unipotent characters of G(e,e,r)
SymbolsDefect:=function(e,r,def,c)
  local defShape, shapesSymbols, IsReducedSymbol, S;

  # Malle-defect of symbol of shape s
  defShape:=function(s)local e; e:=Length(s);
   return (Binomial(e,2)*QuoInt(Sum(s),e)-s*[0..e-1]) mod e;
  end;

  # possible shapes for cuspidal symbols of rank<=r, length e, content=c mod e
  shapesSymbols:=function(r,e,c)local f,res,m,new;
    f:=function(lim2,sum,nb,max)local res,a;
      if nb=1 then 
	if sum=0 then return [[sum]];
	else return [];
	fi;
      fi;
      res:=[];a:=QuoInt(sum,nb-1);
      while a<=max and Binomial(a,2)<=lim2 and a<=sum do
	Append(res,List(f(lim2-Binomial(a,2),sum-a,nb-1,a),
	     x->Concatenation([a],x)));
	a:=a+1;
      od;
      return res;
    end;
    res:=[];m:=0;
    repeat new:=f(r+QuoInt((m*e+c-1)*(m*e+c-e+1),2*e),c+m*e,e,c+m*e);
	   Append(res,new); m:=m+1;
    until Length(new)=0;
    res:=Concatenation(List(res,x->Arrangements(x,e)));
    return Filtered(res,s->defShape(s)=def and 
      ForAll(Rotations(s){[2..Length(s)]},x->defShape(x)<>def or x<=s));
  end;

  IsReducedSymbol:=s->ForAll(Rotations(s){[2..Length(s)]},
    x->s=x or LessSymbols(x,s));

  S:=Concatenation(List(shapesSymbols(r,e,c),s->
    List(PartitionTuples(r-RankSymbol(List(s,x->[0..x-1])),e),
                                     x->SymbolPartitionTuple(x,s))));
  if c<>0 then return S;else return Filtered(S,IsReducedSymbol);fi;
end;

# XSP(rho,s,n[,d]) returns the union of the Lusztig-Spaltenstein 
# \tilde X^{\rho-s,s}_{n,d} for all d even when the 4th argument is present
#                               all d odd otherwise
# In "Character sheaves on disconnected groups II, 13.2" the notation is
#  {}^\rho X^s_{n,d}.
# The result is a list of lists, each one corresponding to a similarity class.
# If s = 0, only positive defects are considered.
# XSP(2,1,n) LS symbols for Sp_2n
# XSP(4,2,n) LS symbols for Sp_2n in char.2
# XSP(2,0,n) LS symbols for SO_{2n+1} [defect odd]
# XSP(2,0,n,true) LS symbols for SO_{2n} [defect even]
# XSP(4,0,n,true) LS symbols for SO_{2n} in char 2
# returns records with fields:
# .symbol  
# .dimBu
# .sp  parameter (double partition) of the generalized Springer correspondent
#                              character of the relative Weyl group
# .Au  describes a character of A(u) as a list: true->sgn, false->Id
#      representing the local system of the Springer correspondent
XSP:=function(arg)local rho,s,n,d,res,S,xsp;
  rho:=arg[1];s:=arg[2];n:=arg[3];
# Lusztig-Spaltenstein \tilde X^{\rho-s,s}_{n,d} 
  xsp:=function(rho,s,n,d)local nrsd;
    nrsd:=rho*QuoInt(d^2,4)-s*(d-d mod 2)/2;
    if n<nrsd then return [];fi;
    return List(PartitionTuples(n-nrsd,2),function(S)
      S:=SymbolPartitionTuple(S,d);
      S:=List(S,function(x)
        if x=[] then return x;else return x+[0..Length(x)-1]*(rho-1);fi;end);
      return [S[1],S[2]+s];end);
  end;

  if Length(arg)=4 then d:=0;else d:=1;fi;
  res:=[];
  repeat S:=xsp(rho,s,n,d);
    if d=0 then S:=List(S,function(x)Sort(x);return x;end);S:=Set(S);fi;
    Append(res,S);d:=d+2;until S=[];
  if s<> 0 and d mod 2<>0 then d:=-1;
    repeat S:=xsp(rho,s,n,d);Append(res,S);d:=d-2;until S=[];
  fi;
  return List(CollectBy(res,x->[ApplyFunc(Union,x),ApplyFunc(Intersection,x)]),
   # here f is a similarity class of symbols
   function(f)local ii,d,i,j,dist,r,n,m;
    ii:=[];d:=ApplyFunc(SymmetricDifference,f[1]);
    if Length(d)>0 then i:=[d[1]];
    for j in d{[2..Length(d)]} do
      if j-i[Length(i)]<rho then Add(i,j);
      else Add(ii,i);i:=[j];
      fi;
    od;
    Add(ii,i); ii:=Filtered(ii,x->x[1]>=s); # intervals
    fi;
    # dist will be the distinguished symbol in its similarity class
    dist:=First(f,function(x)local l;
      if not DefectSymbol(x) in [0,1] then return false;fi;
      l:=[];l{[1,3..2*Length(x[1])-1]}:=x[1];l{[2,4..2*Length(x[2])]}:=x[2];
      return SortingPerm(l)=();end);
    n:=Sum(dist,Length);d:=DefectSymbol(dist);m:=QuoInt(n,2);
    i:=Concatenation(dist);Sort(i);i:=Reversed(i);
    n:=i*[0..n-1]-rho*m*(m-1)*(4*m-5+6*d)/6-s*m*(m+d-1);# geck-malle 2.22
    return List(f,function(S)local r;
      r:=function(x,s)if x=[] then return x;
       else return Reversed(Filtered(x-[0..Length(x)-1]*rho-s,y->y<>0));fi;end;
      r:=rec(symbol:=S,sp:=[r(S[1],0),r(S[2],s)],dimBu:=n);
      if DefectSymbol(S)=0 then
        if r.sp>Reversed(r.sp) then r.sp:=Reversed(r.sp);fi;
	if r.sp[1]=r.sp[2] then r.sp:=[r.sp[1],2,0];fi;
      elif DefectSymbol(S)<0 then
# See Shoji "Unipotent characters of finite classical groups" 5.7.3 and 5.8
        r.sp:=Reversed(r.sp);
      fi;
      r.Au:=List(ii,i->Intersection(S[1],i)<>Intersection(dist[1],i));
      if s=0 and not S[1]=S[2] then
	if r.Au[Length(r.Au)] then r.Au:=List(r.Au,x->not x);fi;
	r.Au:=r.Au{[1..Length(r.Au)-1]};
      fi;
      r.operations := rec(
        Display := function(r, opt)
      	 Print(IntListToString(r.symbol[1]), "\n", IntListToString(r.symbol[2]), "\n");
      	 end,
      	String := r -> SPrint("(symbol:=", PartitionTupleToString(r.symbol),
      	 	", sp:=", PartitionTupleToString(r.sp),
      		", dimBu:=", r.dimBu,
      		", Au:=", Format(r.Au),
      		")"),
      	Print := function(r) Print(String(r)); end);
      return r;end);
    end);
end;
