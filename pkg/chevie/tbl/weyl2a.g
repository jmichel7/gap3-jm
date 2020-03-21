#############################################################################
##
#A  tbl/weyl2a.g             CHEVIE library                     Frank Luebeck
##
#Y  Copyright (C) 1992 - 2001  The CHEVIE Team
##
##  This file contains  the data for the Coxeter coset  W(A_l).F0 where F0
##  induces the non trivial graph automophism on the diagram of of type A.
##
#############################################################################

#############################################################################
##
#F  WordsClassRepresentatives( <n> ) . . . . . . representatives of conjugacy
##  classes (of minimal length)
##
# args: n[, partitions]   to avoid 2 calls of Partitions in ClassInfo
CHEVIE.AddData("WordsClassRepresentatives", "2A",
function(arg)local n, part, guesslongest, redw,l,w0,p;
  n := arg[1];
  if Length(arg) > 1 then part := arg[2]; else part := Partitions(n+1); fi;
  # make word from permutation:
  redw:=function(n,w) local l, i; l:=[];
    while true  do
      i:=PositionProperty([1..n],j->j^w>(j+1)^w);
      if i=false then return l;fi;
      Add( l, i ); w := (i,i+1) * w;
    od;
    return l;
  end;
  # returns longest element in class p  see Geck/Kim ...
  guesslongest := function(p) local x, off, i;
    p:=Concatenation(Filtered(p,i->i mod 2=0),
	  	     Filtered(p,i->i<>1 and i mod 2=1));
    x:=(); off:=0;
    for i in p do
      x:=x*Product(List([2..i],j->(l[off+1],l[off+j])));
      off:=off+i;
    od;
    return x;
  end;
  l:=[]; w0:=();
  for p in [1..QuoInt(n+1,2)] do Append(l,[p,n-p+2]); w0:=w0*(p,n-p+2); od;
  if n mod 2 =0 then Add(l,QuoInt(n,2)+1); fi;
  return List(part,p->redw(n,guesslongest(p)*w0));
end);

##  _checkshortest:=function(W,x)local l, p, els;
##    l:=Length(x);
##    p:=EltWord(W,x)/LongestCoxeterElement(W);
##    els:=Elements(ConjugacyClass(W,p))*LongestCoxeterElement(W);
##    return ForAll(els,p->CoxeterLength(W,p)>=l);
##  end;

#############################################################################
##
#F  CHEVIE.R("ClassInfo","2A")( <n> ) . . conjugacy classes for type 2A
##
##  'CHEVIE.R("ClassInfo","2A")' returns a record with three components:
##    classtext:   representatives of minimal length in  the  conjugacy
##                 classes, as words in generators in standard order
##    classparams:  partitions, parameterizing the classes
##    classnames:  strings for partitions
##    classes:  size of classes
##
##  The  ordering  corresponds  to the  order of  the columns of the ordinary
##  character table of the symmetric  group $S_{n+1}$, as returned by the GAP
##  function 'CharTable("Symmetric", <n+1>)'.
##
CHEVIE.AddData("ClassInfo","2A",function(n)local res;
  res :=CHEVIE.R("ClassInfo", "A")(n);
  res.classtext :=
         CHEVIE.R("WordsClassRepresentatives", "2A")(n, res.classparams);
  Unbind(res.orders);
  return res;
end);

CHEVIE.AddData("NrConjugacyClasses", "2A", n->NrPartitions(n+1));

#############################################################################
##
#F  ClassParam2A( <n>, <w> )  . . . . . . . . . . . . . class parameter of w
##
##  given an element w of a Coxeter group W of type A_n as word in standard
##  generators,  'ClassParam2A' returns the  classparams of its F-conjugacy
##  class under the nontrivial F action permuting the generators.
##
CHEVIE.AddData("ClassParameter","2A",function(n,w)local i,j,x,res,mark,cyc;
  x:=(); for i in w do x:=x*(i,i+1); od;
  # shifting from w_0 class to conjugacy class:
  for i in [1..QuoInt(n+1,2)] do x:=x*(i,n+2-i); od;
  res:=[];
  mark:=[1..n+1];
  for i in [1..n+1] do
    if mark[i]<>0 then
      cyc:=CyclePermInt(x,i);
      Add(res,Length(cyc));
      for j in cyc do mark[j]:=0; od;
    fi;
  od;
  Sort(res);
  return Reversed(res);
end);

# parameters for characters
CHEVIE.AddData("CharParams","2A",n->Partitions(n+1));
CHEVIE.AddData("CharName","2A",function(arg)return IntListToString(arg[2]);end);

CHEVIE.AddData("CharInfo","2A",n->CHEVIE.R("CharInfo","A")(n));

#############################################################################
##
#F  CHEVIE.R("CharTable","2A")( <l> ) outer character table of W(A_l)
##
##  This   function  returns   the   part  of   the   character  table   of
##  CoxeterGroup("A",l).2  on  the outer  classes.  Each character  of
##  CoxeterGroup("A",l)  has  two  extensions to  the  whole  group.
##  'CharTable2A' gives the values of the *preferred* extensions defined in
##  [CS,17.2, case A_l].
##
CHEVIE.AddData("CharTable","2A",function(r)local i, tbl;
  tbl := CHEVIE.R("CharTable","A")(r);
  tbl.identifier := SPrint("W(^2A",r,")");
  for i in [1..Length(tbl.irreducibles)] do
    # Lusztig [Character Sheaves] 17.2:
    # Preferred extension: \sigma acts on \tilde E by (-1)^a_E.w_0
    tbl.irreducibles[i]:=(-1)^CHEVIE.R("LowestPowerFakeDegree","A")
        (tbl.irredinfo[i].charparam)*tbl.irreducibles[i];
  od;
  Inherit(tbl,CHEVIE.R("ClassInfo","2A")(r));
  return tbl;
end);

CHEVIE.AddData("HeckeCharTable","2A", function(r,param,rootparam)
  local q, v, W, T, qE, H, tbl, cl, i;
  q:=-param[1][1]/param[1][2];
  if not IsBound(rootparam[1]) then v:=GetRoot(q,2,"CharTable(Hecke(2A))");
  else v:=rootparam[1];
  fi;
  W:=CoxeterGroup("A",r);
# If q_E is the square root which deforms to 1 of the eigenvalue of T_{w_0}
# on E which deforms to 1, then we have:
#  E~(T_w\phi)=\overline(E(T_{w^-1w_0}))q_E (trivial extension)
#  E~(T_w\phi)=(-1)^a_E\overline(E(T_{w^-1w_0}))q_E (preferred extension)
# where \overline means q->q^-1
  qE:=HeckeCentralMonomials(Hecke(W,v));
  H:=Hecke(W,v^-2);T:=Basis(H,"T");
  tbl:=ShallowCopy(CharTable(H));
  Inherit(tbl,CHEVIE.R("ClassInfo","2A")(r));
  tbl.identifier:=SPrint("H(^2A",r,")");
  cl:=List(tbl.classtext,x->T(EltWord(W,x)*LongestCoxeterElement(W)));
  tbl.irreducibles:=TransposedMat(List(cl,HeckeCharValues));
  for i in [1..Length(tbl.irreducibles)] do
    tbl.irreducibles[i]:=(-1)^CHEVIE.R("LowestPowerFakeDegree","A")
      (tbl.irredinfo[i].charparam[1])*qE[i]*tbl.irreducibles[i];
  od;
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end);

CHEVIE.AddData("PhiFactors","2A",n->List([2..n+1],x->(-1)^x));

CHEVIE.AddData("HeckeRepresentation","2A",function(n,param,sqrtparam,i)
  local H,res,W,p; W:=CoxeterGroup("A",n);
  H:=Hecke(W,-param[1][1]/param[1][2]);p:=Partitions(n+1)[i];
  res:=rec(gens:=SpechtModel(H,p));
  res.F:=Product(res.gens{LongestCoxeterWord(W)})/
    GetRoot(HeckeCentralMonomials(H)[i])*(-1)^CHEVIE.R("LowestPowerFakeDegree","A")(p);
  return res;
end);

CHEVIE.AddData("Representation","2A",function(n,i)
  return CHEVIE.R("HeckeRepresentation","2A")
       (n,List([1..n],x->[1,-1]),[1..n*0+1],i);end);

CHEVIE.AddData("FakeDegree","2A",function(n,p,q)local res;
  res:=CHEVIE.R("FakeDegree","A")(n,p,Indeterminate(Rationals));
  return (-1)^Valuation(res)*Value(res,-q);end);

#############################################################################
#F  PartitionTwoCoreQuotient( <d>, <pp> )
##  returns partition associated to 2-core <d> and pair of partitions <p>.
##
PartitionTwoCoreQuotient:=function(d,p)local x;
# if QuoInt(d,2) mod 2=0 then  # changed jm 9-2-2001 to fit experimental
# evidence -- to be checked in F-S, Lusztig ?
     x:=SymbolPartitionTuple(Reversed(p),-d);
# else x:=SymbolPartitionTuple(p,-d);
# fi;
  return PartBeta(Set(Concatenation(2*x[1],2*x[2]+1)));
end;

# [LuB, 4.4, 4.16, 4.19]
CHEVIE.AddData("UnipotentCharacters","2A",function(l)local  uc, d, k, s, i, r;
  uc:=CHEVIE.R("UnipotentCharacters","A")(l);
  uc.charSymbols:=List(CHEVIE.R("CharParams","A")(l),i->[i]);
  uc.almostHarishChandra:=uc.harishChandra;
  uc.almostHarishChandra[1].relativeType:=rec(orbit:=[rec(series:="A",
     indices:=[1..l],rank:=l)],twist:=Product([1..QuoInt(l,2)],i->(i,l+1-i)),
     rank:=l);
  uc.harishChandra:=[];
  d:=0;
  while d*(d+1)/2 <= l+1 do
    k:=l+1-QuoInt(d*(d+1),2);
    if k mod 2 = 0 then
      r:=QuoInt(k,2);
      s:=rec(levi:=[r+1..l-r],
             relativeType:=rec(series:="B",indices:=[r,r-1..1],rank:=r),
             eigenvalue:=(-1)^(Product(d+[-1..2])/8));
      if d=0 then s.relativeType.cartanType:=1;fi; # type C
      # for the eigenvalue see Lusztig CBMS proof of 3.34 (ii)
      if r<>0 then s.parameterExponents:=Concatenation([2*d+1],0*[2..r]+2);
      else s.parameterExponents:=[];
      fi;
      if k<l then 
        if l-k<10 then s.cuspidalName:=SPrint("{}^2A_",l-k,"");
	else s.cuspidalName:=SPrint("{}^2A_{",l-k,"}");
	fi;
      else s.cuspidalName:="";
      fi;

      # see Fong/Srinivasan for this map
      s.charNumbers:=List(CHEVIE.R("CharParams","B")(r),
        a->Position(uc.charSymbols,[PartitionTwoCoreQuotient(d,a)]));

      FixRelativeType(s);
      Add(uc.harishChandra,s);
    fi;
    d:=d+1;
  od;

  # for delta see Lusztig's book page 124 line 7
  for i in [1..Length(uc.families)] do 
    if 0<> (uc.a[i]+uc.A[i]) mod 2 then 
       uc.families[i]:=Family("C'1",uc.families[i].charNumbers);fi;
  od;
  return uc;
end);

CHEVIE.AddData("UnipotentClasses","2A",function(r,p)local uc,c,t,WF,m;
  uc:=Copy(CHEVIE.R("UnipotentClasses","A")(r,p));
  for c in uc.classes do
    t:=Parent(c.red);
    if Length(t.type)>1 then Error();fi;
    if Length(t.type)=0 or Rank(t)=1 then
      WF:=CoxeterCoset(Parent(c.red));
    else WF:=CoxeterCoset(Parent(c.red),Product([1..QuoInt(t.rank,2)],
      i->(i,t.rank+1-i)));
    fi;
    t:=Twistings(WF,c.red.rootInclusion{c.red.generatingReflections});
    m:=List(t,x->ReflectionEigenvalues(x,PositionClass(x,x.phi)));
    m:=List(m,x->Number(x,y->y=1/2));
    p:=Position(m,Maximum(m));
    c.F:=t[p].phi;
  od;
  return uc;end);
