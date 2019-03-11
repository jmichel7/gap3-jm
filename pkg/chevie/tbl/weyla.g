#############################################################################
##
#A  tbl/weyla.g             CHEVIE library        Goetz Pfeiffer, Jean Michel
##
#Y  Copyright (C) 1994 - 2001  The CHEVIE Team
##
##  This   file  contains   special  functions   for  Coxeter   groups  and
##  Iwahori-Hecke algebras of type A.
##
CHEVIE.AddData("CartanMat","A",function ( n ) local  a, i;
  a := IdentityMat(n);
  for i  in [ 1 .. n ]  do
    a[i][i] := 2;
    if i<n then a[i][i+1] := -1; fi;
    if i>1 then a[i][i-1] := -1; fi;
  od;
  return a;
end);

CHEVIE.AddData("ReflectionDegrees", "A", n->[2..n+1]);

CHEVIE.AddData("PrintDiagram","A",function(r,indices,title)local i;
  Print(title," ",Join(indices," - "),"\n");
end);

# simple roots in R^{l+1} as in Bourbaki
CHEVIE.AddData("GeneratingRoots", "A", function(l)local r, i;
  r := List([1..l],i->0*[1..l+1]);
  for i in [1..l] do r[i]{[i,i+1]} := [1, -1]; od;
  return r;
end);

CHEVIE.AddData("ParabolicRepresentatives", "A", function(l,s)
  return CHEVIE.R("ParabolicRepresentatives","imp")(1,1,l,s);end);

#############################################################################
##
#F  ClassInfo( <n> )  . . . . . . . . . . . conjugacy classes for type A
##
##  ClassInfo returns a record with five components:
##    classtext:   representatives of minimal length in  the  conjugacy
##                 classes, as words in generators in standard order
##    classparams:  partitions, parameterizing the classes
##    classnames:  strings for partitions
##    classes:  size of classes
##    orders:  orders of classes
##
##  The  ordering  corresponds  to the  order of  the columns of the ordinary
##  character table of the symmetric  group $S_{n+1}$, as returned by the GAP
##  function 'Char(acter)Table("Symmetric", <n+1>)'.
##

# very good representatives in the sense of [Geck-Michel]
CHEVIE.AddData("WordClass", "A", function(pi)local w,i,l,r;
  w:= []; i:= 0;
  for l in pi do #  build alternating cycle for each part.
    r:=l mod 2;
    Append(w,i+Concatenation([1,3..l-1-r],[2,4..l+r-2]));
    i:= i+l;
  od;
  return w;
end);

CHEVIE.AddData("ClassInfo","A",function(n)local res;
  res := rec(classparams := Partitions(n+1));
  res.classnames := List(res.classparams, IntListToString);
  res.classtext :=List(res.classparams, CHEVIE.R("WordClass", "A"));
  res.classes:=List(res.classparams,pi->Factorial(n+1)/
    CharTableSymmetric.centralizers[1](n,pi));
  res.orders:=List(res.classparams,Lcm);
  return res;
end);

CHEVIE.AddData("NrConjugacyClasses", "A", n-> NrPartitions(n+1));

CHEVIE.AddData("WeightInfo","A",n->rec(minusculeWeights:=[1..n],
  decompositions:=List([1..n],i->[i]),
  moduli:=[n+1]));

#############################################################################
##
#V  ClassParameter( <W>, <w> )  . . . . . .  class parameter of w
##
##  given  an element w  of a Coxeter  group W of  type A, given as word in
##  standard generators, returns the classparam of its conjugacy class.
##
CHEVIE.AddData("ClassParameter","A",function(n, w) local i, x, res, mark, cyc;
  x:=(); for i in w do x:=x*(i,i+1); od; res:=[]; mark:=[1..n+1];
  for i in [1..n+1] do
    if mark[i]<>0 then cyc:=CyclePermInt(x,i); Add(res,Length(cyc));
      mark{cyc}:=cyc*0;
    fi;
  od;
  Sort(res);
  return Reversed(res);
end);

#############################################################################
##
#F  characters for type A
##
CHEVIE.AddData("CharParams","A",n->Partitions(n+1));
CHEVIE.AddData("LowestPowerFakeDegree","A",p->p*[0..Length(p)-1]);
CHEVIE.AddData("HighestPowerFakeDegree","A",p->Sum(p)*(Sum(p)-1)/2
   -CHEVIE.R("LowestPowerFakeDegree","A")(AssociatedPartition(p)));

# A.CharName(rank,param[,opt])
CHEVIE.AddData("CharName","A",function(arg)return IntListToString(arg[2]);end);

CHEVIE.AddData("CharInfo","A",function(n)local res;
  res:=rec(charparams:=CHEVIE.R("CharParams","A")(n));
  res.extRefl:=List([0..n],
     i->Position(res.charparams,Concatenation([n+1-i],[1..i]*0+1)));
  res.b:=List(res.charparams,CHEVIE.R("LowestPowerFakeDegree","A"));
  res.B:=List(res.charparams,CHEVIE.R("HighestPowerFakeDegree","A"));
  res.a:=res.b; res.A:=res.B;
  return res;
end);

#############################################################################
##  Char(acter)Table is essentially taken from GAP library, but labeling for
##  classes and characters is added. Unfortunately the code for this is
##  quite different for GAP 3 and GAP 4. The GAP 4 code cannot even be read
##  by GAP 3.
CHEVIE.AddData("CharTable", "A", CHEVIE.compat.CharTableA);

#############################################################################
##
#F  HeckeCharTableA( <n>, <q> )  . . . . . . . . character table of $H(A_n)$.
##
##  'HeckeCharTableA'  returns  the character  table  of the Hecke algebra of
##  type $A_n$ with parameter <q>.
##
CHEVIE.tmp:=ShallowCopy(CharTableSymmetric);

CHEVIE.tmp.identifier:= "HeckeA";

CHEVIE.tmp.specializedname := nq-> SPrint("H(A", nq[1]-1, ")");

CHEVIE.tmp.order:=nq->Factorial(nq[1]);
CHEVIE.tmp.size:=nq->Factorial(nq[1]);

CHEVIE.tmp.domain:=nq->IsList(nq) and Length(nq)=2 and IsInt(nq[1]) and nq[1]>0;

CHEVIE.tmp.text:= "generic character table of Hecke algebras of type A";

CHEVIE.tmp.classparam:=[nq->Partitions(nq[1])];
CHEVIE.tmp.charparam:=[nq->Partitions(nq[1])];

# overwrite whole entry because of the ShallowCopy above.
CHEVIE.tmp.irreducibles := [[function(nq, gamma, pi)
   local n, q, k, val, alpha, dif, AHk;

   n:= nq[1]; q:= nq[2]; #  for the sake of clearness.

   #  termination condition.
   if n = 0 then return q^0; fi;

   k:= pi[1]; #  get length of longest cycle.

   val:= 0 * q;

   AHk:=CHEVIE.R("Hk","A").irreducibles[1][1];
   #  loop over the partitions of n-k.
   for alpha in Partitions(n-k) do
      dif:= DifferencePartitions(gamma, alpha);
      if dif <> false then
         val:= val + (q-1)^(dif.cc-1) * (-1)^dif.ll * q^(k-dif.ll-dif.cc)
          *AHk([n-k, q], alpha, pi{[2..Length(pi)]});
      fi;
   od;

   # return the result.
   return val;
end]];

##  Info.. statements commented since not compatible with GAP 4 (FL)
CHEVIE.tmp.matrix:= function(nq)
  local scheme, beta, pm, i, m, n, q, k, t, res, charCol, hooks;
  n:= nq[1]; q:= nq[2]; pm:= []; scheme:= [];
  #  how to encode all hooks.
  hooks:= function(beta, m)
    local i, j, hk, pr, cbs, prs, leg, hks, ll, lg, lh, gamma, new;
    prs:= [];
    #  find all hooks.
    for i in beta do leg:= 0;
      for j in [i-1,i-2..0] do
	if  j in beta then leg:= leg + 1;
	else Add(prs, [rec(from:= i, to:= j, leg:= leg)]);
	fi;
      od;
    od;
    #  construct combinations.
    cbs:= ShallowCopy(prs); hks:= List([1..m], x->[]);
    for hk in cbs do 
      #  extend.
      for pr in prs do
	if pr[1].to>hk[Length(hk)].from then Add(cbs,Concatenation(hk,pr));fi;
      od;
      #  encode.
      ll:=Sum(hk,x->x.from-x.to);lg:=Sum(hk,x->x.leg);lh:= Length(hk);
      new:= rec(wgt:= (-1)^lg * q^(ll-lg-lh) * (q-1)^(lh-1), adr:= 1);
      #  recalculate address.
      if ll < m-1 then
	gamma:= Difference(beta, List(hk, x-> x.from));
	UniteSet(gamma, List(hk, x-> x.to));
	if 0 in gamma then
	  j:= 0; while gamma[j+1] = j do j:= j+1; od;
	  gamma:= gamma{[j+1..Length(gamma)]} - j;
	fi;
	new.adr:= Position(pm[m-ll], gamma);
      fi;
      #  insert.
      Add(hks[ll], new);
    od;
    return hks;
  end;
  #  collect hook encodings.
##  InfoCharTable2("#I  Scheme: \c");
  for i in [1..n] do
##    InfoCharTable2(i, " \c");
    pm[i]:= List(Partitions(i), BetaSet);
    scheme[i]:= [];
    for beta in pm[i] do Add(scheme[i], hooks(beta, i)); od;
  od;
##  InfoCharTable2("done.\n");
  #  how to construct a new column.
  charCol:= function(n, t, k) local col, pi, hk, val;
    col:= [];
    for pi in scheme[n] do val:= 0*q;
      for hk in pi[k] do val:= val + hk.wgt * t[hk.adr]; od;
      Add(col, val);
    od;
    return col;
  end;
  pm:= List([1..n-1], x-> []);
  #  construct the columns.
##  InfoCharTable2("#I  Cycles: \c");
  for m in [1..QuoInt(n,2)] do
##     InfoCharTable2(m, " \c");
     Add(pm[m], charCol(m, [1], m));
     for k in [m+1..n-m] do
       for t in pm[k-m] do Add(pm[k], charCol(k, t, m)); od;
     od;
  od;
##  InfoCharTable2("done.\n");
##  InfoCharTable2("#I  Tables: \c");
  res:= [];
  for k in [1..n-1] do
##    InfoCharTable2(k, " \c");
    for t in pm[n-k] do Add(res, charCol(n, t, k)); od;
  od;
##  InfoCharTable2("done.\n");
  Add(res, charCol(n, [1], n));
  return TransposedMat(res);
end;

CHEVIE.AddData("Hk","A",ShallowCopy(CHEVIE.tmp));

CHEVIE.AddData("HeckeCharTable", "A", CHEVIE.compat.HeckeCharTableA);

#############################################################################
##
#F  PoincarePolynomial( <n>, <q> ) . . Poincare polynomial for type A.
##
##  PoincarePolynomial returns the Poincare  polynomial of the
##  Coxeter group $W$ of type $A_n$, ie.  the sum of $q^l(w)$ over all
##  elements $w$ of $W$.
##
CHEVIE.AddData("PoincarePolynomial","A",function(n,param)
   return Product([1..n],i->Sum([0..i],k->(-param[1][1]/param[1][2])^k));
end);

#############################################################################
##
#F SchurElement( <n>, <q> ) . . . . . . . Schur elements for type A.
##
##  SchurElement returns the list of Schur elements for the character
##  table of the Hecke algebra of type $A_n$ with parameter $q$.
##
##  [Reference: Carter II, 13.5, p. 446.]
##
CHEVIE.AddData("SchurElement","A",function(n,alpha,param,sqrtparam)
   local i, j, lambda, res,q;
   q:=-param[1][1]/param[1][2];
   lambda:= BetaSet(alpha);
   res:= q^Binomial(Length(lambda), 3);
   for i in lambda do
      for j in [0..i-1] do
         if j in lambda then res:= res/q^j;
         else res:= res * Sum([0..i-j-1], e->q^e);
         fi;
      od;
   od;
   return res;
end);

CHEVIE.AddData("FactorizedSchurElement","A",function(arg)
  return CHEVIE.R("FactorizedSchurElement","imp")
    (1,1,arg[1]+1,[arg[2]],arg[3],[]);
end);

CHEVIE.AddData("HeckeRepresentation","A",function(n,param,sqrtparam,i)local H;
  H:=Hecke(CoxeterGroup("A",n),-param[1][1]/param[1][2]);
  return SpechtModel(H,Partitions(n+1)[i]);
end);

CHEVIE.AddData("Representation","A",function(n,i)
  return CHEVIE.R("Representation","imp")(1,1,n+1,i){[2..n+1]};end);

CHEVIE.AddData("FakeDegree","A",function(n,p,q)return
  CHEVIE.R("PoincarePolynomial","A")(Sum(p)-1,[[q,-1]])/
  CHEVIE.R("SchurElement","A")(Sum(p)-1,p,[[q,-1]],[]);end);

CHEVIE.AddData("DecompositionMatrix","A",function(l,p)
  return [[[1..NrPartitions(l+1)],
    MatrixDecompositionMatrix(DecompositionMatrix(Specht(p,p),l+1))]];
end);
 
CHEVIE.AddData("UnipotentCharacters","A",function(l)local ci;
  ci:=CHEVIE.R("CharInfo","A")(l);
  return rec(harishChandra:=[rec(
    levi:=[],
    relativeType:=rec(series:="A",indices:=[1..l],rank:=l),
    parameterExponents:=[1..l]*0+1,
    cuspidalName:="",
    eigenvalue:=1,
    charNumbers:=[1..Length(ci.charparams)])],
    families:=List([1..Length(ci.charparams)],i->Family("C1",[i])),
    charParams:=ci.charparams,
    charSymbols:=List(ci.charparams,x->[BetaSet(x)]),
    a:=ci.a,
    A:=ci.A);
end);

CHEVIE.AddData("Invariants","A",function(n)local m;
  m:=CHEVIE.A.GeneratingRoots(n);Add(m,[1..n+1]*0+1);
  return List([2..n+1],i->function(arg)local v;
    v:=ShallowCopy(arg);Add(v,0*v[1]);v:=v*m;
    return Sum(Arrangements([1..n+1],i),a->Product(v{a}));end);
  end);

CHEVIE.AddData("UnipotentClasses","A",function(n,p)local uc,i,j,cl,d,ss;
  uc:=rec(classes:=List(Partitions(n+1),p->rec(parameter:=p)),
    springerSeries:=Concatenation(List(DivisorsInt(n+1),d->
   List(PrimeResidues(d),i->rec(relgroup:=CoxeterGroup("A",(n+1)/d-1),
        Z:=[E(d)^i],levi:=Filtered([1..n+1],i->i mod d<>0),locsys:=[])))));
  ss:=z->First(uc.springerSeries,x->x.Z=[z]);
  for i in [1..Length(uc.classes)] do
    cl:=uc.classes[i];p:=cl.parameter;d:=Gcd(p);
    cl.name:=IntListToString(p);
    cl.Au:=ComplexReflectionGroup(d,1,1);
    cl.balacarter:=Concatenation(List([1..Length(p)],i->
      Sum(p{[1..i-1]})+[1..p[i]-1]));
    p:=Concatenation(List(p,x->[1-x,3-x..x-1]));Sort(p);
    cl.dynkin:=List([1..Length(p)-1],i->p[i+1]-p[i]);
    cl.red:=[];p:=1;for j in Collected(cl.parameter) do
      Append(cl.red,[p..p+j[2]-2]);p:=p+j[2];od;
    cl.red:=ReflectionSubgroup(CoxeterGroup("A",p-2),cl.red);
    cl.AuAction:=ExtendedReflectionGroup(cl.red,[IdentityMat(cl.red.rank)]);
    if d=2 then Add(ss(1).locsys,[i,2]);Add(ss(-1).locsys,[i,1]);
    else for j in [0..d-1] do Add(ss(E(d)^j).locsys,[i,j+1]);od;
    fi;
  od;
  uc.orderClasses:=Hasse(Poset(List(uc.classes,x->List(uc.classes,y->
            Dominates(y.parameter,x.parameter)))));
  return uc;
end);

CHEVIE.AddData("KLeftCellRepresentatives","A",
function(n)local W,l,f;
  f:=function(i)
    if i<>() then i:=Product(CoxeterWord(W,i),j->(j,j+1));fi;
    i:=List(RobinsonSchenstedCorrespondent(n+1,i).P,Length);
    return Position(CharParams(W),[i]);
  end;
  W:=CoxeterGroup("A",n);
  l:=Filtered(Elements(W),x->x^2=());
  return List(l,x->rec(duflo:=OnTuples([1..n],x),reps:=[],character:=[f(x)]));
  end);
