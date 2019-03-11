#############################################################################
##
#A  tbl/weylbc.g                CHEVIE library    Goetz Pfeiffer, Jean Michel
##
#Y  Copyright (C) 1994 - 2012  The CHEVIE Team
##
##  This   file  contains   special  functions   for  Coxeter   groups  and
##  Iwahori-Hecke algebras of type B and C.
##  The  functions which behave  differently depending on  the isogeny type
##  take an argument after the rank, the cartanType.
##   cartanType=1 type C
##   cartanType=2 type B
##   cartanType=ER(2) and rank=2 Suzuki groups.
##  other  values  happen  for  root  systems  which  occur  inside complex
##  reflection groups.
##
CHEVIE.AddData("CartanMat","B",function (arg) local  n,type,a;
  n:=arg[1];
  if Length(arg)=2 then type:=arg[2];else type:=2;fi;
  a:=CHEVIE.R("CartanMat","A")(n);
  a[1][2]:=-type;a[2][1]:=2/a[1][2];
  return a;
end);

CHEVIE.AddData("PrintDiagram","B",function(r,indices,title,type)local i;
  Print(title," ");
  if type=1 then Print(indices[1]," >=> ",indices[2]);
  elif type=2 then Print(indices[1]," <=< ",indices[2]);
  elif type=ER(2) then Print(indices[1]," = ",indices[2]);
  else Print(indices[1]," ?=? ",indices[2]);
  fi;
  for i in [3..r] do Print(" - ",indices[i]);od;
  Print("\n");
end);

CHEVIE.AddData("ReflectionName","B",function(arg)local i,r,type,option;
  r:=arg[1];option:=arg[2];
  if Length(arg)=3 then type:=arg[3];else type:=2;fi;
  if type=2 then 
    if IsBound(option.TeX) then return SPrint("B_",TeXBracket(r));
    elif IsBound(option.arg) then return SPrint("\"B\",",r);
    else return SPrint("B",r);fi;
  elif type=1 then 
    if IsBound(option.TeX) then return SPrint("C_",TeXBracket(r));
    elif IsBound(option.arg) then return SPrint("\"C\",",r);
    else return SPrint("C",r); fi;
  elif type=ER(2) then 
    if IsBound(option.TeX) then return SPrint("B^{\\hbox{sym}}_",TeXBracket(r));
    elif IsBound(option.arg) then return SPrint("\"Bsym\",",r);
    else return SPrint("Bsym",r); fi;
  elif IsBound(option.TeX) then 
         return SPrint("B^?_",TeXBracket(r),"(",Format(type,option),")");
  elif IsBound(option.arg) then return SPrint("\"B?\",",r,",",type);
  else return SPrint("B?",r,"(",Format(type),")");
  fi;
end);

CHEVIE.AddData("GeneratingRoots", "B", function(l,type)local rts, i;
  rts := List([1..l],i->0*[1..l]);
  for i in [1..l-1] do rts[i]{[i,i+1]} := [1,-1]; od;
  rts[l][l] := 2/type;
  return rts{[l,l-1..1]}; # CHEVIE ordering
end);

CHEVIE.AddData("ParabolicRepresentatives", "B", function(l,s)
  return CHEVIE.R("ParabolicRepresentatives","imp")(2,1,l,s);end);

CHEVIE.AddData("ReflectionDegrees","B", n->2*[1..n]);

CHEVIE.AddData("Size", "B",function(arg)return 2^arg[1]*Factorial(arg[1]);end);

CHEVIE.AddData("NrConjugacyClasses", "B", function(arg)
  return NrPartitionTuples(arg[1], 2);
end);

CHEVIE.AddData("WeightInfo","B",function(n,type)
  if type=2 then return rec(minusculeWeights:=[1],minusculeCoweights:=[n],
    decompositions:=[[1]],moduli:=[2]);
  else           return rec(minusculeWeights:=[n],minusculeCoweights:=[1],
    decompositions:=[[1]],moduli:=[2]);
  fi;
end);

#############################################################################
##
#F  WordClass( <pi> )  . . . very good representative in the sense of
##  [Geck-Michel] of the conjugacy class parametrized by pi 
##
CHEVIE.AddData("WordClass", "B", function(pi) local  w, i, l,r;
  w:= [];  i:= 1;
  for l in Reversed(pi[2]) do #  handle signed cycles in reversed order.
    Append(w, [i, i-1 .. 2]); Append(w, [1 .. i+l-1]); i:= i+l;
  od;
  for l in pi[1] do #  the unsigned cycles.
    r:=l mod 2;
    Append(w,i+Concatenation([1,3..l-1-r],[2,4..l+r-2]));
    i:= i+l;
  od;
  return w;
end);

#############################################################################
##
#F  ClassInfo( <n> ) . . . . . . . . conjugacy classes for type B
##
##  ClassInfo returns a record with three components:
##    classtext:   representatives of minimal length in  the  conjugacy
##                 classes, as words in generators in standard order
##    classparams:  double partitions, parameterizing the classes
##    classnames:  strings representing the double partitions
##
##  The  ordering  corresponds  to the  order of  the columns of the ordinary
##  character table of the Coxeter group of type $B_n$, as returned by the GAP
##  function 'CharTable("WeylB", <n>)'.
##

CHEVIE.AddData("ClassInfo", "B", function(n) local res;
  res:=CHEVIE.R("ClassInfo","imp")(2,1,n);
  res.classtext:=List(res.classparams,CHEVIE.R("WordClass", "B"));
  res.classes:=List(res.centralizers,x->res.centralizers[1]/x);
  return res;
end);

#############################################################################
##
#V  ClassParameter( n, <w> )  . . . . . .  class parameter of w in type Bn.
##
##  returns the classparam of an element w of a Coxeter group of type Bn.
##
CHEVIE.AddData("ClassParameter", "B", function(n,w)local x,i,res,mark,cyc,j;
  x:=();
  for i in w do
    if i=1 then x:=x*(1,n+1); else x:=x*(i-1,i)(i-1+n,i+n);fi;
  od;
  res:=[[],[]];
  mark:=[1..n];
  for i in [1..n] do
    if mark[i]<>0 then
      cyc:=CyclePermInt(x,i);
      if i+n in cyc then Add(res[2],Length(cyc)/2);
      else Add(res[1],Length(cyc));
      fi;
      for j in cyc do
        if j>n then mark[j-n]:=0; else mark[j]:=0; fi;
      od;
    fi;
  od;

  Sort(res[1]); Sort(res[2]);
  return [Reversed(res[1]),Reversed(res[2])];
end);

#############################################################################
##
#F  CharParams( <n> ) . . . . . . . . . . . . . . . characters for type B
##
CHEVIE.AddData("CharParams", "B", n->PartitionTuples(n,2));

#how to make charname from charparam
CHEVIE.AddData("CharName", "B",
  function(arg) return PartitionTupleToString(arg[2]); end);

CHEVIE.AddData("LowestPowerFakeDegree","B", function(p) local pp, m, res;
  pp:=SymbolPartitionTuple(p,1); m:=Length(pp[2]);
  res:=pp[1]*[m,m-1..0];
  if pp[2]<>[] then res:=res+pp[2]*[m-1,m-2..0]; fi;
  return 2*res+Sum(pp[2])-m*(m-1)*(4*m+1)/6;
end);

CHEVIE.AddData("CharInfo","B",function(n)local res;
  res:=rec(charparams:=CHEVIE.R("CharParams","B")(n));
  res.extRefl:=Concatenation(
   List([0..n-1],i->Position(res.charparams,[[n-i],[1..i]*0+1])),
   [Position(res.charparams,[[],[1..n]*0+1])]);
  res.a:=List(res.charparams,
    p->LowestPowerGenericDegreeSymbol(SymbolPartitionTuple(p,1)));
  res.A:=List(res.charparams,
    p->HighestPowerGenericDegreeSymbol(SymbolPartitionTuple(p,1)));
  res.b:=List(res.charparams,CHEVIE.R("LowestPowerFakeDegree","B"));
  res.B:=res.a+res.A-res.b;
  return res;
end);

##  essentially the library call to table "WeylB", l, but different
##  in GAP 3 and GAP 4
CHEVIE.IndirectAddData("CharTable", ["B"], CHEVIE.compat.CharTableB);

#############################################################################
##
#V  HeckeCharTable( <n>, <para> ). . .  character table of $H(B_n)$.
##
##  HeckeCharTable returns the character  table of the Hecke
##  algebra  $H(B_n)$ with parameters <para>.
##
CHEVIE.tmp := ShallowCopy(CharTableWeylB);

CHEVIE.tmp.identifier := "HeckeB";

CHEVIE.tmp.specializedname:=nq->SPrint("H(B",nq[1],")");

CHEVIE.tmp.order:=nq-> 2^nq[1]*Factorial(nq[1]);
CHEVIE.tmp.size:=nq-> 2^nq[1]*Factorial(nq[1]);

CHEVIE.tmp.domain:= function(nq)
   return IsList(nq) and Length(nq) = 3 and IsInt(nq[1]) and nq[1] > 0;
end;

CHEVIE.tmp.text:= "generic character table of Hecke algebras of type B";

CHEVIE.tmp.classparam:=[nq -> PartitionTuples(nq[1],2)];
CHEVIE.tmp.charparam:=[nq -> PartitionTuples(nq[1],2)];

CHEVIE.tmp.irreducibles := [[function(nq, gamma, pi)
   local n, q, Q, k, val, t, nn, alpha, dif, BHk;

   n:= nq[1];  q:= nq[3];  Q:= nq[2]; #  for the sake of clearness.

   #  termination condition.
   if n = 0 then
      return q^0;
   fi;

   val:= 0*q; #  initialize character value.

   BHk:=CHEVIE.R("Hk","B").irreducibles[1][1];

   #  positive cycles first.
   if pi[1] <> [] then

      k:= pi[1][1]; #  get length of the longest cycle.

      #  loop over double paritions of n-k.
      for alpha in PartitionTuples(n-k, 2) do
         dif:= [];
         dif[1]:= DifferencePartitions(gamma[1], alpha[1]);
         dif[2]:= DifferencePartitions(gamma[2], alpha[2]);
         if dif[1] <> false and dif[2] <> false then
            dif:= rec(cc:= dif[1].cc + dif[2].cc, ll:= dif[1].ll + dif[2].ll);
            val:= val + (q-1)^(dif.cc-1) * (-1)^dif.ll * q^(k-dif.ll-dif.cc)
                * BHk([n-k, Q, q], alpha,[pi[1]{[2..Length(pi[1])]}, pi[2]]);
         fi;
      od;

   else # pi[2] = []

      #  get length of the longest cycle.
      k:= pi[2][1];

      #  loop over k-hooks in gamma[1].
      nn:= Sum(gamma[1]);
      if nn >= k then
         for alpha in Partitions(nn - k) do
            dif:= DifferencePartitions(gamma[1], alpha);
            if dif <> false and dif.cc = 1 then
               val:= val + Q * (-1)^dif.ll * q^(n+dif.d)
                 * BHk([n-k, Q, q], [alpha, gamma[2]],
                                    [pi[1], pi[2]{[2..Length(pi[2])]}]);
            fi;
         od;
      fi;

      #  loop over hooks in gamma[2].
      nn:= Sum(gamma[2]);
      if nn >= k then
         for alpha in Partitions(nn - k) do
            dif:= DifferencePartitions(gamma[2], alpha);
            if dif <> false and dif.cc = 1 then
               val:= val + (-1)^(dif.ll+1) * q^(n+dif.d)
                 * BHk([n-k, Q, q], [gamma[1], alpha],
                                    [pi[1], pi[2]{[2..Length(pi[2])]}]);
            fi;
         od;
      fi;

   fi;

   #  return the result.
   return val;

end]];

##  Info... statements commented out to make it GAP4 compatible
CHEVIE.tmp.matrix:= function(nq)
  local scheme,beta,pm,i,m,k,t,n,x,y,np,col,res,charCol,hooks,DoublePartitions;
  ############################################################################
  #F  DoublePartitions( <n> ) . . . . . . . . . . . . .  pairs of partitions.
  ##
  ##  JM:  seems  to  be  used  just  for  the  Sortex  at the end. Should be
  ##  suppressed by using a different enumeration...
  DoublePartitions:= function(n) local m, k, pm, t, s, res;
    if n = 0 then return [[[], []]]; fi;
    pm:= List([1..n], x->[]);
    #  second position.
    for m in [1..n] do
       Add(pm[m], [[], [m]]); #  add the m-cycle.
       for k in [m+1..n] do
	  for t in pm[k-m] do
	     s:= [[], [m]]; Append(s[2], t[2]); Add(pm[k], s);
	  od;
       od;
    od;
    #  first position.
    for m in [1..QuoInt(n,2)] do
       Add(pm[m], [[m], []]); #  add the m-cycle.
       for k in [m+1..n-m] do
	  for t in pm[k-m] do
	     s:= [[m], t[2]]; Append(s[1], t[1]); Add(pm[k], s);
	  od;
       od;
    od;
    #  collect.
    res:= [];
    for k in [1..n-1] do
       for t in pm[n-k] do
	  s:= [[k], t[2]]; Append(s[1], t[1]); Add(res, s);
       od;
    od;
    Add(res, [[n], []]); Append(res, pm[n]);
    return res;
  end;
  n:= nq[1]; x:= nq[3]; y:= nq[2]; pm:= []; scheme:= [];
  hooks:=function(beta,m) #  how to encode all hooks.
    local i,j,k,hk,pr,cbs,prs,leg,hks,lb,ll,lg,lh,gamma,new;
    hks:= List([1..m], x->[]); prs:= [];
    lb:= [Length(beta[1]), Length(beta[2])];

    #  find all hooks.
    for i in [1, 2] do
      prs[i]:= [];
      for j in beta[i] do leg:= 0;
	for k in Reversed([0..j-1]) do
	  if  k in beta[i] then leg:= leg + 1;
	  else Add(prs[i],rec(from:=j, to:=k, leg:=leg, pow:= m+k-lb[i]));
	  fi;
	od;
      od;
    od;

    #  construct combinations.
    cbs:= List(prs[1], x-> [[x], []]);
    Append(cbs, List(prs[2], x-> [[], [x]]));
    for hk in cbs do

       #  extend.
       for pr in prs[1] do
	  if hk[2] = [] and pr.to > hk[1][Length(hk[1])].from then
	     new:= List(hk, ShallowCopy); Add(new[1], pr); Add(cbs, new);
	  fi;
       od;
       for pr in prs[2] do
	  if hk[2] = [] or pr.to > hk[2][Length(hk[2])].from then
	     new:= List(hk, ShallowCopy); Add(new[2], pr); Add(cbs, new);
	  fi;
       od;
       #  encode.
       ll:= Sum(hk[1], x-> x.from - x.to) + Sum(hk[2], x-> x.from - x.to);
       lg:= Sum(hk[1], x-> x.leg) + Sum(hk[2], x-> x.leg);
       lh:= Length(hk[1]) + Length(hk[2]);
       new:= rec(wgt:= [(-1)^lg * x^(ll-lg-lh) * (x-1)^(lh-1), 0], adr:= 1);
       if lh = 1 then
	  if IsBound(hk[1][1]) then new.wgt[2]:= (-1)^lg * y * x^hk[1][1].pow;
	  else                      new.wgt[2]:= (-1)^(lg+1) * x^hk[2][1].pow;
	  fi;
       fi;
       #  recalculate address.
       if ll < m then
	  gamma:= [];
	  for i in [1, 2] do
	    gamma[i]:= Difference(beta[i], List(hk[i], x-> x.from));
	    UniteSet(gamma[i], List(hk[i], x-> x.to));
	    if 0 in gamma[i] then
	      j:= 0;
	      while j < Length(gamma[i]) and gamma[i][j+1] = j do j:= j+1; od;
	      gamma[i]:= gamma[i]{[j+1..Length(gamma[i])]}-j;
	    fi;
	  od;
	  new.adr:= Position(pm[m-ll], gamma);
       fi;
       Add(hks[ll], new); #  insert.
    od;
    return hks;
  end;

  #  collect hook encodings.
##     InfoCharTable2("#I  Scheme: \c");
  for i in [1..n] do
##        InfoCharTable2(i, " \c");
     pm[i]:= List(PartitionTuples(i, 2), p-> List(p, BetaSet));
     scheme[i]:= [];
     for beta in pm[i] do
	Add(scheme[i], hooks(beta, i));
     od;
  od;
##     InfoCharTable2("done.\n");

  #  how to construct a new column.
  charCol:= function(n, t, k, p) local col, pi, hk, val; col:= [];
     for pi in scheme[n] do
	val:= 0*y;
	for hk in pi[k] do val:= val + hk.wgt[p] * t[hk.adr]; od;
	Add(col, val);
     od;
     return col;
  end;

  #  construct the columns.
##     InfoCharTable2("#I  Cycles: \c");
  pm:= List([1..n], x->[]);

  #  second position.
  for m in [1..n] do
##        InfoCharTable2(m, " \c");
     Add(pm[m], charCol(m, [1], m, 2)); #  add the m-cycle.
     for k in [m+1..n] do
	for t in pm[k-m] do
	   Add(pm[k], charCol(k, t, m, 2));
	od;
     od;
  od;

  #  first position.
  for m in [1..QuoInt(n,2)] do
##        InfoCharTable2(m, " \c");
     Add(pm[m], charCol(m, [1], m, 1)); #  add the m-cycle.
     for k in [m+1..n-m] do
	for t in pm[k-m] do
	   Add(pm[k], charCol(k, t, m, 1));
	od;
     od;
  od;
##     InfoCharTable2("done.\n");

  #  collect.
##     InfoCharTable2("#I  Tables: \c");
  res:= [];
  for k in [1..n-1] do
##        InfoCharTable2(k, " \c");
     for t in pm[n-k] do Add(res, charCol(n, t, k, 1)); od;
  od;
  Add(res, charCol(n, [1], n, 1)); Append(res, pm[n]);
##     InfoCharTable2("done.\n");
  res:=Permuted(res,Sortex(DoublePartitions(n))/Sortex(PartitionTuples(n,2)));
  return TransposedMat(res);
end;

CHEVIE.AddData("Hk","B",CHEVIE.tmp);
Unbind(CHEVIE.tmp);

CHEVIE.IndirectAddData("HeckeCharTable",["B"],CHEVIE.compat.HeckeCharTableB);

#############################################################################
##
#F  PoincarePolynomial( <n>, <para> )  Poincare polynomial for type B.
##
##  PoincarePolynomial returns the Poincare polynomial of the
##  Coxeter group $W$ of type $B_n$, ie.  the sum of $q^l(w)$ over all
##  elements $w$ of $W$.
##
CHEVIE.AddData("PoincarePolynomial","B",function(n, para)local q1,q2;
   q1:=-para[1][1]/para[1][2]; q2:=-para[2][1]/para[2][2];
   return Product([0..n-1], i-> (q2^i*q1 + 1) * Sum([0..i], k-> q2^k));
end);

#############################################################################
##
#F  SchurElement( <n> , <char>, <param>, <sqrtparam> ) . . . . Schur element
#F  for type B.
##
##  JM 2/2011: With the new division-free formula by Chlouveraki, it is as
##  fast to call the general routine for G(2,1,n)
##
CHEVIE.AddData("SchurElement","B",function(arg)
  return CHEVIE.R("SchurElement","imp")(2,1,arg[1],arg[2],arg[3],[]);
end);

CHEVIE.AddData("FactorizedSchurElement","B",function(arg)
  return CHEVIE.R("FactorizedSchurElement","imp")(2,1,arg[1],arg[2],arg[3],[]);
end);

CHEVIE.AddData("HeckeRepresentation","B",function(arg)
return CHEVIE.R("HeckeRepresentation","imp")(2,1,arg[1],arg[2],[],arg[4]);
end);

CHEVIE.AddData("Representation","B",function(n,i)
  return CHEVIE.R("Representation","imp")(2,1,n,i);
end);

#############################################################################
#F  FakeDegree( n, <c>,q ) Fake Degree of char. with charparam <c>
##
CHEVIE.AddData("FakeDegree","B",function(n,c,q)
  return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,1)),q);
end);

CHEVIE.AddData("DecompositionMatrix","B",function(l,p)local pp,dd,pt,decS;
  decS:=i->MatrixDecompositionMatrix(DecompositionMatrix(Specht(p,p),i));
  pp:=List([0..l],Partitions); pt:=PartitionTuples(l,2);
  if p=2 then # restriction matrix from  Bn to Sn * dec(Sn)
    return [[[1..Length(pt)],List(pt,function(p)
      p:=LittlewoodRichardsonRule(p[1],p[2]);
      return List(pp[l+1],function(x)
         if x in p then return 1; else return 0;fi;end);end)*decS(l)]];
  else dd:=Concatenation([[[1]],[[1]]],List([2..l],decS));
    return List([0..l],i->
      [List(Cartesian(pp[i+1],pp[l+1-i]),x->Position(pt,x)),
       List(Cartesian(dd[i+1],dd[l+1-i]),x->List(Cartesian(x),Product))]);
  fi;
end);
 
## some cleanup for "B" <-> "C" ? XXX
CHEVIE.AddData("UnipotentCharacters","B",function(arg)
  local uc,symbols,r,d,s,rank;
  rank:=arg[1];
  uc:=rec(harishChandra:=[],charSymbols:=[]);
  for d in 1+2*[0..QuoInt(-1+RootInt(1+4*rank,2),2)] do
    s:=(d^2-1)/4;
    s:=rec(relativeType:=rec(series:="B",indices:=[1+s..rank],rank:=rank-s),
	   levi:=[1..s],
	   eigenvalue:=(-1)^QuoInt(d+1,4),
	   parameterExponents:=Concatenation([d],[2+s..rank]*0+1),
	   cuspidalName:=SPrint("B_",TeXBracket(s)));
    Add(uc.harishChandra,s);
    symbols:=Symbols(rank,d);
    s.charNumbers:=[1..Length(symbols)]+Length(uc.charSymbols);
    FixRelativeType(s);
    Append(uc.charSymbols,symbols);
  od;
  uc.harishChandra[1].cuspidalName:="";
  uc.a:=List(uc.charSymbols,LowestPowerGenericDegreeSymbol);
  uc.A:=List(uc.charSymbols,HighestPowerGenericDegreeSymbol);
  uc.families:=FamiliesClassical(uc.charSymbols);
  if Length(arg)=2 and arg[2]=1 then
    uc.harishChandra[1].relativeType.cartanType:=1;
  fi;
  return uc;
end);

# References:
# [Lu]   G.Lusztig,   Character   sheaves   on   disconnected   groups,  II
# Representation Theory 8 (2004) 72--124
#
# [GM]  M.Geck and G.Malle, On the existence of a unipotent support for the
# irreducible  characters of  a finite  group of  Lie type,  Trans. AMS 352
# (1999) 429--456
# 
# [S]   N.Spaltenstein,  Classes  unipotentes  et  sous-groupes  de  Borel,
# Springer LNM 946 (1982)
#
#  The function is called with type=1 for type C and type=2 for type B
#
CHEVIE.AddData("UnipotentClasses","B",function(r,type,char)local cl,uc,
  i,l,s,cc,ss,symbol2para,part2dynkin,addSpringer,d,LuSpin,trspringer,j;
  part2dynkin:=function(part)local p,res; #partition->Dynkin-Richardson diagram
    p:=Concatenation(List(part,d->[1-d,3-d..d-1]));
    Sort(p);p:=p{[QuoInt(3+Length(p),2)..Length(p)]};
    if type=1 then res:=[2*p[1]];else res:=[p[1]];fi;
    Append(res,p{[2..Length(p)]}-p{[1..Length(p)-1]});
    return res;
  end;
  addSpringer:=function(s)local ss,p;
    ss:=First(uc.springerSeries,x->x.defect=DefectSymbol(s.symbol));
    if s.sp=[[],[]] then p:=1;
    elif s.sp=[[1],[]] then p:=2;
    elif s.sp=[[],[1]] then p:=1;
    else p:=Position(CharParams(ss.relgroup),[s.sp]);
    fi;
    ss.locsys[p]:=[Length(uc.classes),Position(CharParams(cc.Au),
      List(s.Au,function(x)if x then return [1,1];else return [2];fi;end))];
  end;
  if type=ER(2) then type:=2;char:=2;fi; #treat 2B2 as B2; make sure char=2
  if char=2 then ss:=XSP(4,2,r);
  elif type=1 then ss:=XSP(2,1,r);
  else ss:=XSP(2,0,r);
  fi;
  l:=Union(List(ss,c->List(c,x->[DefectSymbol(x.symbol),Sum(x.sp,Sum)])));
  SortBy(l,x->[AbsInt(x[1]),-SignInt(x[1])]);
  uc:=rec(classes:=[],springerSeries:=List(l,function(d)local res;
   res:=rec(relgroup:=CoxeterGroup("C",d[2]),defect:=d[1],locsys:=[],
     levi:=[1..r-d[2]]);
   if char=2 then res.Z:=[1];
   elif type=1 then res.Z:=[(-1)^(r-d[2])];
   elif IsInt(ER(2*(r-d[2])+1)) then res.Z:=[1];else res.Z:=[-1];fi;
   return res;end));
  if char<>2 then # Invert [GM] 2.6 and 2.10
   symbol2para:=function(S)local c,i,l,part,d; #Lusztig symbol->Partition
     c:=Concatenation(S);Sort(c); i:=1;part:=[];
     d:=type mod 2;
     while i<=Length(c) do
       if i=Length(c) or c[i+1]-c[i]>0 then Add(part,2*(c[i]-(i-1))+1-d);i:=i+1;
       else l:=2*(c[i]-(i-1))-d;Append(part,[l,l]);i:=i+2;
       fi;
     od;
     Sort(part);part:=Filtered(part,y->y<>0);
     return Reversed(part);
   end;
  else # Invert [GM] 2.7
   symbol2para:=function(S)local c,i,l,part,ex;
     c:=Concatenation(S);Sort(c);i:=1;part:=[];ex:=[];
     while i<=Length(c) do
       if i=Length(c) or c[i+1]-c[i]>1 then Add(part,2*(c[i]-2*(i-1)));i:=i+1;
       elif c[i]=c[i+1] then
	 l:=2*(c[i]-2*(i-1))-2;Append(part,[l,l]);Add(ex,l);i:=i+2;
       elif c[i]+1=c[i+1] then
	 l:=2*(c[i]-2*(i-1))-1;Append(part,[l,l]);i:=i+2;
       fi;
     od;
     Sort(part);part:=Filtered(part,y->y<>0);
     return [Reversed(part),ex];
   end;
  fi;
  if char=2 then type:=1;fi;
  for cl in ss do
    cc:=rec(parameter:=symbol2para(cl[1].symbol));
    cc.Au:=ApplyFunc(CoxeterGroup,Concatenation(List(cl[1].Au,x->["A",1])));
    if char<>2 then
      cc.dynkin:=part2dynkin(cc.parameter);
      cc.name:=IntListToString(cc.parameter);
    else type:=1; 
      cc.dimBu:=cl[1].dimBu;
      cc.name:=Join(List(Reversed(Collected(cc.parameter[1])),
	function(x)local res;res:=IntListToString([1..x[2]]*0+x[1],"[]");
	  if x[1] in cc.parameter[2] then return SPrint("(",res,")");fi;
	  return res;end),"");
    fi;
    cc.red:=CoxeterGroup();
    if char=2 then j:=cc.parameter[1];else j:=cc.parameter;fi;
    for j in Collected(j) do
      if j[1]mod 2=type mod 2 then cc.red:=cc.red*CoxeterGroup("C",j[2]/2);
      elif j[2]mod 2<>0 then 
        if j[2]>1 then cc.red:=cc.red*CoxeterGroup("B",(j[2]-1)/2);fi;
      elif j[2]>2 then cc.red:=cc.red*CoxeterGroup("D",j[2]/2);
      else cc.red:=cc.red*Torus(1);
      fi;
    od;
    Add(uc.classes,cc); for s in cl do addSpringer(s);od;
  od;
  uc.orderClasses:=Hasse(Poset(List(uc.classes,x->List(uc.classes,
    function(y)local m,f,fx,fy,i;
     if char<>2 then return Dominates(y.parameter,x.parameter);fi;
     # cf. [S] 2.10 page 24
     m:=Maximum(x.parameter[1][1],y.parameter[1][1]);
     f:=x->List([1..m],i->Sum(Filtered(x,z->z<i))+i*Number(x,z->z>=i));
     fx:=f(x.parameter[1]);fy:=f(y.parameter[1]);
     for i in [1..m] do 
       if fx[i]<fy[i] then return false;
       elif fx[i]=fy[i] and i in y.parameter[2] then
	 if i in Difference(x.parameter[1],x.parameter[2])then return false;fi;
	 if i<m and (fx[i+1]-fy[i+1]) mod 2=1 then return false;fi;
       fi;
     od;
     return true;end))));
  if char<>2 and type=2 then
    LuSpin:=function(p)local t,a,b,i,j,l,d; # cf [Lu] 14.2
      Sort(p);a:=[];b:=[];d:=[0,1,0,-1];d:=d{List(p,x->1+x mod 4)};
      i:=1;
      while i<=Length(p) do l:=p[i];t:=Sum(d{[1..i-1]});
	if 1=l mod 4 then Add(a,(l-1)/4-t);i:=i+1;
	elif 3=l mod 4 then Add(b,(l-3)/4+t);i:=i+1;
	else j:=i;while i<=Length(p) and p[i]=l do i:=i+1;od;j:=[1..(i-j)/2]*0;
	  Append(a,j+(l+l mod 4)/4-t);Append(b,j+(l-l mod 4)/4+t);
	fi;
      od;
      a:=Filtered(a,x->x<>0);a:=Reversed(a);
      b:=Filtered(b,x->x<>0);b:=Reversed(b);
      if Sum(d)>=1 then return [a,b];else return [b,a];fi;
    end;
    addSpringer:=function(f,i,s,k)local ss,p;
      ss:=First(uc.springerSeries,f);
      if s in [[[],[1]],[[],[]]] then p:=1;
      elif s=[[1],[]] then p:=2;
      else p:=Position(CharParams(ss.relgroup),[s]);fi;
      ss.locsys[p]:=[i,k];
    end;
    trspringer:=function(i,old,new)local ss,c,p;
      for ss in uc.springerSeries do for c in ss.locsys do
	if c[1]=i then p:=Position(old,c[2]);
	  if p<>false then c[2]:=new[p];fi;fi;od;od;
    end;
    d:=0;
    while 4*d^2-3*d<=r do i:=4*d^2-3*d;
      if (r-d) mod 2=0 then
	l:=Concatenation([1..i],[i+2,i+4..r]);
	Add(uc.springerSeries,rec(relgroup:=CoxeterGroup("B",(r-i)/2),
	  levi:=l,Z:=[-1],locsys:=[]));
	i:=4*d^2+3*d;
	if i<=r and d<>0 then 
	  l:=Concatenation([1..i],[i+2,i+4..r]);
	  Add(uc.springerSeries,
	  rec(relgroup:=CoxeterGroup("B",(r-i)/2),levi:=l,Z:=[-1],locsys:=[]));
	fi;
      fi;
      d:=d+1;
    od;
    l:=Filtered([1..Length(uc.classes)],i->
      ForAll(Collected(uc.classes[i].parameter),c->c[1] mod 2=0 or c[2]=1));
    for i in l do
      cl:=uc.classes[i];
      s:=LuSpin(cl.parameter);
      if Size(cl.Au)=1 then cl.Au:=CoxeterGroup("A",1);trspringer(i,[1],[2]);
        d:=1;
      elif Size(cl.Au)=4 then cl.Au:=CoxeterGroup("B",2);
        trspringer(i,[1,2,3,4],[1,3,5,4]);d:=2;
      else Error("Au non-commutative of order ",Size(cl.Au)*2," not implemented");
      fi;
      addSpringer(ss->ss.Z=[-1] and ss.relgroup.rank=Sum(s,Sum),i,s,d);
    od;
  fi;
  return uc;
end);

CHEVIE.AddData("Invariants","B",function (n,type)local  m;
  m:=[1..n]*0+1;m[1]:=2/type;
  m:=DiagonalMat(m)*CHEVIE.R("GeneratingRoots","imp")(2,1,n);
  return List(CHEVIE.imp.Invariants(2,1,n),f->function(arg)
    return ApplyFunc(f,arg*m);end);
end);
