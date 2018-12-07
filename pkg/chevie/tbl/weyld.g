#############################################################################
##
#A  tbl/weyld.g                CHEVIE library     Goetz Pfeiffer, Jean Michel
##
#Y  Copyright (C) 1994 - 2001  The CHEVIE Team
##
##  This file contains special functions for Coxeter groups and Iwahori-Hecke 
##  algebras  of type  D.
##
CHEVIE.AddData("CartanMat","D",function ( n ) local  a,m;
  if n<3 then m:=3;else m:=n;fi;
  a:=CHEVIE.R("CartanMat","A")(m);
  a{[1..3]}{[1..3]}:=[[2,0,-1],[0,2,-1],[-1,-1,2]];
  return a{[1..n]}{[1..n]};
end);

CHEVIE.AddData("Size", "D", 
        function(arg) return 2^(arg[1]-1)*Factorial(arg[1]); end);

CHEVIE.AddData("PrintDiagram","D",function(r,indices,title)local i,s;
  Print(title," ",indices[1],"\n");s:=String("",Length(title)+1);
  Print(s," \\\n",s,"  ",indices[3]);
  for i in [4..r] do Print(" - ",indices[i]);od;Print("\n");
  Print(s," /\n",s,indices[2],"\n");
end);

CHEVIE.AddData("GeneratingRoots", "D", function(l)local r, rts, i;
  rts := [];
  for i in [1..l-1] do r:=0*[1..l];r{[i,i+1]}:=[1,-1];Add(rts, r); od;
  r:=0*[1..l];r{[l-1,l]}:=[1,1];Add(rts,r);
  return Reversed(rts); # CHEVIE ordering
end);

CHEVIE.AddData("WeightInfo","D",function(n)local res;
  if n mod 2=1 then return rec(minusculeWeights:=[1,2,n],
    decompositions:=[[1],[3],[2]],moduli:=[4]);
  else return rec(minusculeWeights:=[1,2,n],decompositions:=[[1,0],[0,1],[1,1]],
    moduli:=[2,2]);
  fi;
end);

CHEVIE.AddData("ParabolicRepresentatives", "D", function(l,s)
  return CHEVIE.R("ParabolicRepresentatives","imp")(2,2,l,s);end);

# arg: n[, classparams]  (to avoid recomputation in ClassInfo)
CHEVIE.AddData("WordsClassRepresentatives","D", function(arg) 
  local n, param, res, w, i, pi, l,r;
  n := arg[1];
  if Length(arg)=2 then param := List(arg[2], a-> List(a, ShallowCopy));
  else param := PartitionTuples(n, 2);
  fi;
  res:= [];

  #  loop over the labels for type B.
  for pi in param do
    #  take those with an even number of signed parts.
    if pi[2] = '+' then pi[2] := []; fi;
    if IsList(pi[2]) and  Length(pi[2]) mod 2 = 0 then
         w:= []; i:= 1;

         #  handle signed parts in reversed order.
         for l in Reversed(pi[2]) do

            #  the first sign is empty.
	    if i = 1 then Append(w, [2 .. i+l-1]);
	    else Append(w, [i, i-1 .. 3]); Append(w, [1 .. i+l-1]);
            fi;
            i:= i+l;
         od;
   
         #  the unsigned cycles.
         for l in pi[1] do 
	   r:=l mod 2;
	   Append(w,i+Concatenation([1,3..l-1-r],[2,4..l+r-2]));
	   i:= i+l; 
         od;

         #  cosmetics for lexicographics.
         if w <> [] and w[1] = 2 then w[1]:= 1; fi;
   
         # classes are labelled with '+', if they have representatives
	 # in parabolic subgroup of type A_{l-1}, given by {1,3,4,..}
	 if pi[2] = [] and ForAll(pi[1], x->x mod 2=0)  then
           Add(res, w); w := ShallowCopy(w); w[1] := 2;
	 fi;
         Add(res, w);
      fi;
   od;

   #  return the list.
   return res;
end);

#############################################################################
##
#F  ClassInfo( <n> ) . . . . . . .  conjugacy classes for type D.
##
##  ClassInfo returns a record with components:
##    classtext:   representatives of minimal length in  the  conjugacy  
##                 classes, as words in generators in standard order
##    classparams:  double partitions or [partition,sign],
##                 parameterizing the classes
##    classnames:  strings for double partitions or [partition,sign].
##    classes:  cardinality of the classes
##    centralizers:  cardinality of the classes
##  
##  The ordering corresponds to the order of the columns of  the  ordinary
##  character table of the Coxeter group of type $D_n$, as returned by the GAP
##  function 'CharTable("WeylD", <n>)'.
##
CHEVIE.AddData("ClassInfo", "D", function(n) local res;
  res := CHEVIE.R("ClassInfo", "imp")(2,2,n);
  res.classparams:= List(res.classparams,function(x)
    if Length(x)=2 then return x;fi;
    if x[3]=0 then return [x[1],'+'];else return [x[1],'-'];fi;end);
  res.classtext:=CHEVIE.R("WordsClassRepresentatives","D")
     (n,res.classparams);
  return res;
end);

CHEVIE.AddData("NrConjugacyClasses", "D", function(n)
  if n mod 2 = 1 then return NrPartitionTuples(n, 2) / 2;
  else return (NrPartitionTuples(n, 2) + 3*NrPartitions(n/2)) /2;
  fi;
end);

#############################################################################
##
#F  CharInfo( <n> )  . . . . . . . . . . . characters for type D
##  
CHEVIE.AddData("CharInfo","D",n->CHEVIE.R("CharInfo","imp")(2,2,n));

# how to make a .charname from a .charparam
CHEVIE.AddData("CharName","D", 
  function(arg)return PartitionTupleToString(arg[2]);end);

#############################################################################
##
#F  ClassParameter( <n>, <w> )  . . . . . . . . . class parameter of w
##  
##  given an element w  of a Coxeter group W of type A  as word in  standard
##  generators, ClassParameter returns the classparam of its conjugacy class.
##  

# Used by 'ClassParamD' for distinguishing classes with '+' or '-' in label:
# (precomputed for D_n with n=4,6,8)
CHEVIE.AddData("gensMODA","D",
    [,,,
    [[(1,2)(7,8),(3,4)(5,6),(2,3)(6,7),(3,5)(4,6)],[[4],[,,2]],
    [[2],[1,,1]]],,
    [[(1,2)(8,11)(12,14)(15,17)(16,18)(19,21)(22,25)(31,32),
    (3,4)(5,6)(7,9)(10,13)(20,23)(24,26)(27,28)(29,30),
    (2,3)(6,8)(9,12)(13,16)(17,20)(21,24)(25,27)(30,31),
    (3,5)(4,6)(12,15)(14,17)(16,19)(18,21)(27,29)(28,30),
    (5,7)(6,9)(8,12)(11,14)(19,22)(21,25)(24,27)(26,28),
    (7,10)(9,13)(12,16)(14,18)(15,19)(17,21)(20,24)(23,26)],
    [[16],[4,,6],[1,,,,5]],[[12],[2,,6],[,2,,,4]]
    ],,
    [[(1,2)(8,11)(12,15)(16,20)(17,21)(22,26)(23,27)(28,33)
    (29,34)(30,35)(36,41)(37,42)(43,50)(44,51)(52,59)
    (60,68)(61,69)(70,77)(78,85)(79,86)(87,92)(88,93)
    (94,99)(95,100)(96,101)(102,106)(103,107)(108,112)(109,113)
    (114,117)(118,121)(127,128),(3,4)(5,6)(7,9)(10,13)
    (14,18)(19,24)(25,31)(32,38)(39,45)(40,46)(47,53)
    (48,54)(49,55)(56,62)(57,63)(58,64)(65,71)(66,72)
    (67,73)(74,80)(75,81)(76,82)(83,89)(84,90)(91,97)
    (98,104)(105,110)(111,115)(116,119)(120,122)(123,124)(125,126),
    (2,3)(6,8)(9,12)(13,17)(18,23)(20,25)(24,30)
    (26,32)(33,39)(34,40)(41,48)(42,49)(50,57)(51,58)
    (53,61)(59,67)(62,70)(68,76)(71,78)(72,79)(80,87)
    (81,88)(89,95)(90,96)(97,103)(99,105)(104,109)(106,111)
    (112,116)(117,120)(121,123)(126,127),
    (3,5)(4,6)(12,16)(15,20)(17,22)(21,26)(23,29)
    (27,34)(30,37)(35,42)(39,47)(45,53)(48,56)(54,62)
    (57,65)(58,66)(63,71)(64,72)(67,75)(73,81)(76,84)
    (82,90)(87,94)(92,99)(95,102)(100,106)(103,108)(107,112)
    (109,114)(113,117)(123,125)(124,126),
    (5,7)(6,9)(8,12)(11,15)(22,28)(26,33)(29,36)
    (32,39)(34,41)(37,44)(38,45)(40,48)(42,51)(46,54)
    (49,58)(55,64)(65,74)(71,80)(75,83)(78,87)(81,89)
    (84,91)(85,92)(88,95)(90,97)(93,100)(96,103)(101,107)
    (114,118)(117,121)(120,123)(122,124),
    (7,10)(9,13)(12,17)(15,21)(16,22)(20,26)(25,32)
    (31,38)(36,43)(41,50)(44,52)(48,57)(51,59)(54,63)
    (56,65)(58,67)(62,71)(64,73)(66,75)(70,78)(72,81)
    (77,85)(79,88)(86,93)(91,98)(97,104)(103,109)(107,113)
    (108,114)(112,117)(116,120)(119,122),
    (10,14)(13,18)(17,23)(21,27)(22,29)(26,34)(28,36)
    (32,40)(33,41)(38,46)(39,48)(45,54)(47,56)(52,60)
    (53,62)(59,68)(61,70)(67,76)(69,77)(73,82)(75,84)
    (81,90)(83,91)(88,96)(89,97)(93,101)(95,103)(100,107)
    (102,108)(106,112)(111,116)(115,119),
    (14,19)(18,24)(23,30)(27,35)(29,37)(34,42)(36,44)
    (40,49)(41,51)(43,52)(46,55)(48,58)(50,59)(54,64)
    (56,66)(57,67)(62,72)(63,73)(65,75)(70,79)(71,81)
    (74,83)(77,86)(78,88)(80,89)(85,93)(87,95)(92,100)
    (94,102)(99,106)(105,111)(110,115)],
    [[64],[16,,24],[,,32],[4,,,,20],[,,,,,,16]],
    [[56],[12,,24],[6,,28],[2,4,,,18],[1,,3,,,,14]]]]);


CHEVIE.AddData("ClassParameter","D",function(n,w)
  local x, i, res, mark, cyc, j, tmp, gens;
  
  x:=();
  for i in w do
    if i=1 then x:=x*(1,n+2)(2,n+1); else x:=x*(i-1,i)(i-1+n,i+n);fi;
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
        if j>n then mark[j-n]:=0;
        else mark[j]:=0;
        fi;
      od;
    fi;
  od;
  
  if res[2]=[] and ForAll(res[1],i->i mod 2 = 0 ) then
    # for classes with '+' or '-' we use the cycle type for the
    # permutation representation on the cosets of the parabolic
    # subgroup [2..n]: (the generators for this representation are 
    # stored in the global variable 'CHEVIE.R("gensMODA","D")')
    
    if not IsBound(CHEVIE.R("gensMODA","D")[n]) then
      tmp:=CoxeterGroup("D",n);
      gens:=PermCosetsSubgroup (tmp,
                               ReflectionSubgroup(tmp,[2..n]));
      tmp:=CHEVIE.R("ClassInfo","D")(n);
      tmp:=tmp.classtext{Filtered([1..Length(tmp.classnames)],i->
                   '+' in tmp.classnames[i] or '-' in tmp.classnames[i])};
      tmp:=List(tmp,a->CycleStructurePerm(Product(gens{a})));
      CHEVIE.R("gensMODA","D")[n]:=[gens,tmp{2*[1..Length(tmp)/2]-1},
                             tmp{2*[1..Length(tmp)/2]}];
    fi;
    
    tmp:=CycleStructurePerm(Product(CHEVIE.R("gensMODA","D")[n][1]{w}));
    if tmp in CHEVIE.R("gensMODA","D")[n][2] and 
                         not tmp in CHEVIE.R("gensMODA","D")[n][3] then
      res[2]:='+';
    elif not tmp in CHEVIE.R("gensMODA","D")[n][2] 
                         and tmp in CHEVIE.R("gensMODA","D")[n][3] then
      res[2]:='-';
    fi;
  fi;  
  
  Sort(res[1]);
  if IsList(res[2]) then
    Sort(res[2]);
    return [Reversed(res[1]),Reversed(res[2])];
  else
    return [Reversed(res[1]),res[2]];
  fi;
end);

#############################################################################
##
#F  CharTable( <l> ) . . .  character table of CoxeterGroup("D",l)
##  
##  This function returns  the part of the character table the Coxeter group
##  of  type B_l on classes  inside  a reflection subgroup  of  type D_l.
##  
##  If l is even then some of the classes  and restrictions split into two
##  classes or  characters, respectively. Their  values  are given by  the
##  character  values    for W(B_l) and   those  for  the  symmetric group
##  S_(l/2). This is described in [Pfeiffer, G., Character Tables of  Weyl
##  Groups in GAP]. 
##  
CHEVIE.AddData("CharTable","D", CHEVIE.compat.CharTableD);

###########################################################################
##
#F  HeckeCharTable( <n>, <u> ) . . . .  character table of $H(D_n)$.
##
##  HeckeCharTable returns the character table of the Hecke algebra associated
##  with the finite  Coxeter group of  type  $D_n$   with parameter <u>.
##
CHEVIE.tmp:=ShallowCopy(CharTableWeylD);

CHEVIE.tmp.identifier:= "HeckeD"; 

CHEVIE.tmp.specializedname:=nq->SPrint("H(D",nq[1],")");

CHEVIE.tmp.size:=nq->2^(nq[1]-1)*Factorial(nq[1]);  
CHEVIE.tmp.order:=nq->2^(nq[1]-1)*Factorial(nq[1]);  

CHEVIE.tmp.domain:= function(nq)
   return IsList(nq) and Length(nq)=2 and IsInt(nq[1]) and nq[1]>1;
end;

CHEVIE.tmp.text:= "generic character table of Hecke algebras of type D";

CHEVIE.tmp.classparam:= [nq-> CharTableWeylD.classparam[1](nq[1])];
CHEVIE.tmp.charparam:= [nq-> CHEVIE.R("CharInfo","D")(nq[1]).charparams];

CHEVIE.tmp.irreducibles:=[[function(nq, alpha, pi)
   local delta, va,vb, val,n, q,AHk,BHk,s;
   n:=nq[1];q:=nq[2];s:="+-";
   if q=1 then return CharTableWeylD.irreducibles[1][1](n, alpha, pi); fi;
   
   AHk:=CHEVIE.R("Hk","A").irreducibles[1][1];
   BHk:=CHEVIE.R("Hk","B").irreducibles[1][1];

   if not IsList(alpha[2]) then
      delta:= [alpha[1], alpha[1]];
      if not IsList(pi[2]) then
         vb:= BHk([n, 1, q], delta, [pi[1], []])/2;
	 va:=(q+1)^Length(pi[1])/2*AHk([n/2,q^2],alpha[1],pi[1]/2);
         if s[alpha[3]+1]=pi[2] then val:=vb+va;
         else val:=vb-va;
         fi;
      else val:= BHk([n, 1, q], delta, pi)/2;
      fi;
   else
      if not IsList(pi[2]) then val:= BHk([n, 1, q], alpha, [pi[1], []]);
      else val:= BHk([n, 1, q], alpha, pi);
      fi;
   fi;
   return val;
end]];

CHEVIE.AddData("Hk","D",ShallowCopy(CHEVIE.tmp));

CHEVIE.AddData("HeckeCharTable","D",
  CHEVIE.compat.HeckeCharTableD);

CHEVIE.AddData("FactorizedSchurElement","D",function(arg)local p,i,n;
  p:=arg[2];n:=arg[1];
  if p[2] in "+-" then p:=[p[1],p[1]]; fi;
  return CHEVIE.R("FactorizedSchurElement","imp")(2,2,n,p,arg[3],[]);
end);

CHEVIE.AddData("HeckeRepresentation","D",function(arg)local p,i,n;
  i:=arg[4];n:=arg[1];
  p:=CHEVIE.R("CharInfo","D")(n).charparams[i];
  if p[Length(p)]=0 then i:=i+1; elif p[Length(p)]=1 then i:=i-1; fi;
  return CHEVIE.R("HeckeRepresentation","imp")(2,2,n,arg[2],[],i);
end);

CHEVIE.AddData("Representation","D",function(n,i)local p;
  p:=CHEVIE.R("CharInfo","D")(n).charparams[i];
  if p[Length(p)]=0 then i:=i+1; elif p[Length(p)]=1 then i:=i-1; fi;
  return CHEVIE.R("Representation","imp")(2,2,n,i);
end);

#############################################################################
##
#F  PoincarePolynomialD . . . . . . . . . . . Poincare polynomial for type D.
##
CHEVIE.AddData("PoincarePolynomial","D",function(n, para)local q;
  q:=-para[1][1]/para[1][2];
  return Sum([0..n-1],k->q^k)*Product([1..n-1],i->(q^i+1)*Sum([0..i-1],k->q^k));
end);


CHEVIE.AddData("symbolcharparam","D",c->SymbolPartitionTuple(c,0));

CHEVIE.AddData("Invariants","D",function(n)local  m;
  m:=CHEVIE.R("GeneratingRoots","imp")(2,2,n);
  return List(CHEVIE.imp.Invariants(2,2,n),f->function(arg)
    return ApplyFunc(f,arg*m);end);
end);

#############################################################################
##
#F  CycPolGenericDegree( <para> ) .  . . . Generic Degree  for type D.
##
##  CycPolGenericDegree  returns the  generic degree of the character
##  with parameter <para> as a CycPol (see CycPol.g).
##
##  [Reference: Lusztig, 'Irreducible Repr. of classical groups']
##
CHEVIE.AddData("CycPolGenericDegree","D",
       c->CycPolGenericDegreeSymbol(SymbolPartitionTuple(c,0)));

CHEVIE.AddData("SchurElement","D", function(n,phi,q,sqrtparam)
  return CHEVIE.R("PoincarePolynomial","D")(n,q)/
      Value(CHEVIE.R("CycPolGenericDegree","D")(phi),-q[1][1]/q[1][2]);
end);

#############################################################################
##
#F  FakeDegree( n,<c>,q ) Fake Degree of char. with charparam <c>
##
CHEVIE.AddData("FakeDegree","D",function(n,c,q)
  return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,0)),q);
end);

CHEVIE.AddData("UnipotentCharacters","D",function(rank)local uc,symbols,r,d,s;
  uc:=rec(harishChandra:=[],charSymbols:=[]);
  for d in 4*[0..RootInt(QuoInt(rank,4),2)] do
    r:=d^2/4;
    s:=rec(relativeType:=rec(series:="B",indices:=[1+r..rank],rank:=rank-r),
      levi:=[1..r], eigenvalue:=(-1)^QuoInt(d+1,4),
      parameterExponents:=Concatenation([d],[2+r..rank]*0+1));
    if r<10 then s.cuspidalName:=SPrint("D_",r,"");
    else s.cuspidalName:=SPrint("D_{",r,"}");
    fi;
    if d=0 then
      s.relativeType.series:="D";
      s.cuspidalName:="";
      s.parameterExponents[1]:=1;
    fi;
    Add(uc.harishChandra,s);
    symbols:=Symbols(rank,d);
    s.charNumbers:=[1..Length(symbols)]+Length(uc.charSymbols);
    FixRelativeType(s);
    Append(uc.charSymbols,symbols);
  od;
  uc.a:=List(uc.charSymbols,LowestPowerGenericDegreeSymbol);
  uc.A:=List(uc.charSymbols,HighestPowerGenericDegreeSymbol);
  uc.families:=FamiliesClassical(uc.charSymbols);
  return uc;
end);

CHEVIE.AddData("ReflectionDegrees","D",n->Concatenation(2*[1..n-1],[n]));

# References for unipotent classes:
# [Lu] G.Lusztig, Character sheaves on disconnected groups, II 
#   Representation Theory 8 (2004) 72--124
#
# [GM]  M.Geck and G.Malle, On the existence of a unipotent support for the
# irreducible  characters of  a finite  group of  Lie type,  Trans. AMS 352
# (1999) 429--456
# 
# [S]  N.Spaltenstein,  Classes  unipotentes  et  sous-groupes  de
# Borel, Springer LNM 946 (1982)
# 
CHEVIE.AddData("UnipotentClasses","D",function(n,char)local s,uc,cl,cc,l,ss,k,
  d,i,symbol2partition,addSpringer,partition2DR,trspringer,LuSpin,j;
  addSpringer:=function(s,i)local ss,p;
    ss:=First(uc.springerSeries,x->x.defect=DefectSymbol(s.symbol));
    if s.sp in [[[],[1]],[[],[]]] then p:=1;
    elif s.sp=[[1],[]] then p:=2;
    else p:=Position(CharParams(ss.relgroup),[s.sp]);fi;
    ss.locsys[p]:=[i,Position(CharParams(cc.Au),
      List(s.Au,function(x)if x then return [1,1];else return [2];fi;end))];
  end;
  partition2DR:=function(part)local p;
    p:=Concatenation(List(part,x->[1-x,3-x..x-1]));
    Sort(p);p:=p{[1+Length(p)/2..Length(p)]};
    return Concatenation([p[1]+p[2]],List([1..Length(p)-1],i->p[i+1]-p[i]));
  end;
  if char=2 then ss:=XSP(4,0,n,1);
    symbol2partition:=function(S)local c,i,l,part,ex; # see [GM] 2.17
      c:=Concatenation(S);Sort(c); i:=1;part:=[];ex:=[];
      while i<=Length(c) do
	if i=Length(c) or c[i+1]-c[i]>1 then 
	     Add(part,2*(c[i]-2*(i-1))+2);i:=i+1;
	elif c[i+1]-c[i]>0 then 
	     l:=2*(c[i]-2*(i-1))+1;Append(part,[l,l]);i:=i+2;
	else l:=2*(c[i]-2*(i-1));Append(part,[l,l]);i:=i+2;Add(ex,l);
	fi;
      od;
      Sort(part);part:=Filtered(part,y->y<>0);
      return [Reversed(part),ex];
    end;
  else ss:=XSP(2,0,n,1); # see [GM] 2.10
    symbol2partition:=function(S)local c,i,l,part;
      c:=Concatenation(S);Sort(c); i:=1;part:=[];
      while i<=Length(c) do
	if i=Length(c) or c[i+1]-c[i]>0 then Add(part,2*(c[i]-(i-1))+1);i:=i+1;
	else l:=2*(c[i]-(i-1));Append(part,[l,l]);i:=i+2;
	fi;
      od;
      Sort(part);part:=Filtered(part,y->y<>0);
      return Reversed(part);
    end;
  fi;
  l:=Union(List(ss,c->List(c,x->[DefectSymbol(x.symbol),
    Sum(FullSymbol(x.sp),Sum)])));
  SortBy(l,x->[AbsInt(x[1]),-SignInt(x[1])]);
  uc:=rec(classes:=[],springerSeries:=List(l,function(d)local res;
    res:=rec(defect:=d[1],locsys:=[],levi:=[1..n-d[2]]);
    if (n-d[2]) mod 4=0 or char=2 then 
       if n mod 2=0 then res.Z:=[1,1];else res.Z:=[1];fi;
    else 
       if n mod 2=0 then res.Z:=[-1,-1];else res.Z:=[-1];fi;
    fi;
    if d[1]=0 then res.relgroup:=CoxeterGroup("D",d[2]);
    else res.relgroup:=CoxeterGroup("B",d[2]);fi;
    return res;end));
  for cl in ss do
    cc:=rec(parameter:=symbol2partition(cl[1].symbol));
    if char=2 then
      cc.dimBu:=cl[1].dimBu;
      cc.name:=Join(List(Reversed(Collected(cc.parameter[1])),
	function(x)local res;res:=IntListToString([1..x[2]]*0+x[1],"[]");
	  if x[1] in cc.parameter[2] then return SPrint("(",res,")");fi;
	  return res;end),"");
    else
      cc.dynkin:=partition2DR(cc.parameter);
      cc.name:=IntListToString(cc.parameter);
    fi;
    cc.Au:=ApplyFunc(CoxeterGroup,Concatenation(List(cl[1].Au,x->["A",1])));
    CharNames(cc.Au);
    if char<>2 then
      cc.red:=CoxeterGroup();
      j:=cc.parameter;
      for j in Collected(j) do
	if j[1]mod 2=0 then cc.red:=cc.red*CoxeterGroup("C",j[2]/2);
	elif j[2]mod 2<>0 then 
	  if j[2]>1 then cc.red:=cc.red*CoxeterGroup("B",(j[2]-1)/2);fi;
	elif j[2]>2 then cc.red:=cc.red*CoxeterGroup("D",j[2]/2);
	else cc.red:=cc.red*Torus(1);
	fi;
      od;
#   else cc.red:=?????:
    fi;
    if not IsList(cl[1].sp[2]) then cl[1].sp[3]:=1-((n/2) mod 2);fi;
    Add(uc.classes,cc); for s in cl do addSpringer(s,Length(uc.classes));od;
    if not IsList(cl[1].sp[2]) then
      cl[1].sp[3]:=1-cl[1].sp[3];Add(cc.name,'+');
      cc:=Copy(cc);cc.name[Length(cc.name)]:='-';
      if IsBound(cc.dynkin) then cc.dynkin{[1,2]}:=cc.dynkin{[2,1]};fi;
      Add(uc.classes,cc); for s in cl do addSpringer(s,Length(uc.classes));od;
    fi;
  od;
  if char=2 then # cf. [S] 2.10 page 24
    uc.orderClasses:=Hasse(Poset(List(uc.classes,x->List(uc.classes,
     function(y)local m,f,fx,fy,i;
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
      if x.parameter=y.parameter and x<>y then return false;fi;
      return true;end))));
  else
    uc.orderClasses:=Hasse(Poset(List([1..Length(uc.classes)],
    i->List([1..Length(uc.classes)],
    j->Dominates(uc.classes[j].parameter,uc.classes[i].parameter) 
      and (uc.classes[j].parameter<>uc.classes[i].parameter or i=j)))));
  fi;
  
  if char<>2 then
    d:=0;
    while 4*d^2-d<=n do i:=4*d^2-d;
      if (n-d) mod 2=0 then
	l:=Concatenation([1..i],[i+2,i+4..n]);
	s:=rec(relgroup:=CoxeterGroup("B",(n-i)/2),levi:=l,locsys:=[]);
	if n mod 2=0 then s.Z:=[1,-1];else s.Z:=[E(4)];fi;
	Add(uc.springerSeries,s);
        if d=0 then l:=Concatenation([1],[4,6..n]);fi;
	s:=rec(relgroup:=CoxeterGroup("B",(n-i)/2),levi:=l,locsys:=[]);
	if n mod 2=0 then s.Z:=[-1,1];else s.Z:=[-E(4)];fi;
	Add(uc.springerSeries,s);
	i:=4*d^2+d;
	if d<>0 and i<=n then
	  l:=Concatenation([1..i],[i+2,i+4..n]);
	  s:=rec(relgroup:=CoxeterGroup("B",(n-i)/2),levi:=l,locsys:=[]);
	  if n mod 2=0 then s.Z:=[1,-1];else s.Z:=[E(4)];fi;
	  Add(uc.springerSeries,s);
	  s:=rec(relgroup:=CoxeterGroup("B",(n-i)/2),levi:=l,locsys:=[]);
	  if n mod 2=0 then s.Z:=[1,1];else s.Z:=[-E(4)];fi;
	  Add(uc.springerSeries,s);
	fi;
      fi;
      d:=d+1;
    od;
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
    trspringer:=function(i,new)local ss,c;
      for ss in uc.springerSeries do for c in ss.locsys do
	if c[1]=i then c[2]:=new[c[2]];fi;od;od;
    end;
    l:=Filtered([1..Length(uc.classes)],i->
      ForAll(Collected(uc.classes[i].parameter),c->c[1] mod 2=0 or c[2]=1));
    for i in l do
      cl:=uc.classes[i];
      s:=LuSpin(cl.parameter);
      if Size(cl.Au)=1 then cl.Au:=CoxeterGroup("A",1);
        trspringer(i,[2]);k:=[1,1];
      elif Size(cl.Au)=2 then cl.Au:=CoxeterGroup("A",1,"A",1);
        trspringer(i,[2,4]);k:=[1,3];
      elif Size(cl.Au)=8 then cl.Au:=CoxeterGroup("A",1,"B",2);
        trspringer(i,[1,6,8,5,10,3,4,9]);k:=[2,7]; # error?
      else Error("Au non-commutative of order ",Size(cl.Au)*2,
                 " not implemented");
      fi;
      if not '-' in cl.name then addSpringer(ss->ss.Z in [[1,-1],[E(4)]] and 
	  ss.relgroup.rank=Sum(s,Sum),i,s,k[1]);
      fi;
      if not '+' in cl.name then addSpringer(ss->ss.Z in [[-1,1],[-E(4)]] and 
          ss.relgroup.rank=Sum(s,Sum),i,s,k[2]);
      fi;
    od;
  fi;
  return uc;
end);
