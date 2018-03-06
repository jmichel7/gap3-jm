#############################################################################
##
#A  tbl/weyl2d.g               CHEVIE library                   Frank Luebeck
##
#Y  Copyright (C) 1994 - 2001  The CHEVIE Team
##
##  This  file contains  data  for  the Coxeter  coset  of  type ^2D.  One
##  convenient way of considering this is to look at the non-trivial coset
##  of W(D_l) inside W(B_l).
##

CHEVIE.AddData("ClassParams", "2D", function(n)local B;
  B := CHEVIE.R("ClassParams", "B")(n);
  return Filtered(B, a-> Length(a[2]) mod 2 = 1);
end);

CHEVIE.AddData("WordsClassRepresentatives", "2D", function(n)
  return CHEVIE.R("ClassInfo", "2D")(n).classtext;
end);

#############################################################################
##
#F  CHEVIE.R("ClassInfo","2D")( <n> ) . . conjugacy classes for type 2D
##
##  'CHEVIE.R("ClassInfo","2D")' returns a record with three components:
##    classtext:   representatives of minimal length in  the  conjugacy
##                 classes, as words in generators in standard order
##    classparams:  partitions, parameterizing the classes
##    classnames:  strings for partitions
##
##  We  parametrize the  F-conjugacy classes  by the  classes in  the coset
##  Bn-Dn.  If n is  odd, since F  is inner acting  as w0, it would also be
##  possible  to  parametrize  them  by  the  classes  in  Dn  (you get the
##  F-classes by translating with w0). This gives two possible labelings of
##  the  F-classes: one by the Dn-classes  and one by the outer Bn-classes.
##  They  correspond as follows: Let [a,b] the double partition for w in Dn
##  and let [c,d] be the double partition for w.w0Bn in Bn. Then c contains
##  the  even entries of a and the odd entries of b and d contains the even
##  entries of b and the odd entries of a.
##
CHEVIE.AddData("ClassInfo","2D",function(n)local  l, B;
  B:=CHEVIE.R("ClassInfo","B")(n);
  l:=Filtered([1..Length(B.classtext)],i->Length(B.classparams[i][2]) mod 2=1);
  return rec(classnames:=B.classnames{l},
             classparams:=B.classparams{l},
             classes:=B.classes{l},
             classtext:=List(B.classtext{l},function(l)local res,i,n;
  # 1 is the automorphism, and u=121 is the new generator of W(2D). We deal
  # with words with an odd number of 1. To gather one 1 at right we go left
  # to right doing substitutions 11->, 12->u1, and 1a->a1 if a<>2.
    res:=[];n:=1;
    for i in [1..Length(l)] do
      if l[i]=1 then n:=(n+1) mod 2;
      elif l[i]=2 then Add(res,2-n);
      else Add(res,l[i]);
      fi;
    od;
    return res;
  end));
end);

CHEVIE.AddData("NrConjugacyClasses", "2D", function(n)
  if n mod 2 = 1 then return NrPartitionTuples(n, 2) / 2;
  else return (NrPartitionTuples(n, 2) - NrPartitions(n/2)) /2;
  fi;
end);

#############################################################################
##
#F  CHEVIE.R("ClassParameter","2D")( <n>, <w> )  . class parameter of w
##
##  given  an element w of a Coxeter group  W of type D as word in standard
##  generators,    'CHEVIE.R("ClassParameter","2D")'    returns   the
##  classparam  of its F-conjugacy  class under the  F action permuting the
##  first two generators.
##
CHEVIE.AddData("ClassParameter","2D",function(n,w)local x, i, res, mark, cyc, j;
  x:=();
  for i in w do
    if i=1 then x:=x*(1,n+2)(2,n+1); else x:=x*(i-1,i)(i-1+n,i+n);fi;
  od;

  # now shifting in the coset:
  x:=x*(1,n+1);

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

  Sort(res[1]);Sort(res[2]);
  return [Reversed(res[1]),Reversed(res[2])];
end);

CHEVIE.AddData("IsPreferred","2D",
# test if a character of W(B) corresponds to the preferred extension
# for ^2D, see [CS,17.2] and [Lusztig-book,4.4,4.18]:
 function(pp)pp:=SymbolPartitionTuple(pp,0);return pp[1]>pp[2];end);

CHEVIE.AddData("IsGood","2D",pp->pp[1]>pp[2]);
# whether a character of W(B) corresponds to the "good" extension for 2D

CHEVIE.AddData("testchar","2D",CHEVIE.R("IsPreferred","2D"));

CHEVIE.AddData("CharParams","2D",n->Filtered(CHEVIE.R("CharParams","B")
         (n),CHEVIE.R("testchar","2D")));

CHEVIE.AddData("CharName","2D",
                function(arg) return PartitionTupleToString(arg[2]); end);

CHEVIE.AddData("CharInfo","2D",function(n)local res,resparams;
  res:=rec(charparams:=CHEVIE.R("CharParams","2D")(n));
  res.extRefl:=List([0..n-2],i->[[1..i]*0+1,[n-i]]);
  Append(res.extRefl,[[[1],[1..n-1]*0+1],[[],[1..n]*0+1]]);
  res.extRefl:=List(res.extRefl,x->PositionProperty(res.charparams,y->y=x
     or y=Reversed(x)));
  resparams:=CHEVIE.R("CharInfo","D")(n).charparams;
  res.charRestrictions:=List(res.charparams,x->PositionProperty(resparams,
    y->y=x or y=Reversed(x)));
  res.nrGroupClasses:=Length(resparams);
  return res;
end);

#############################################################################
##
#F  CHEVIE.R("CharTable","2D")( <l> ) outer character table of
##                                          CoxeterGroup("D",l).2
##
##  This function returns the part of the character table the Coxeter group
##  of  type B_l on classes  outside a reflection subgroup  of type D_l for
##  the characters which remain irreducible on restriction to this subgroup
##  and which correspond to the *preferred* extensions defined in [CS,17.2,
##  case D_l].
##
##  Alternatively  you can get the  *good* extension instead of *preferred*
##  extension by defining testchar appropriately.
##
CHEVIE.AddData("CharTable","2D", CHEVIE.compat.CharTable2D);

CHEVIE.AddData("HeckeCharTable", "2D", CHEVIE.compat.HeckeCharTable2D);

CHEVIE.AddData("FakeDegree","2D",function(n,c,q)
  return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,0),1),q);end);

CHEVIE.AddData("PhiFactors","2D",function(n)local res;
  res:=[1..n-1]*0+1; Add(res,-1); return res;
end);

CHEVIE.AddData("UnipotentCharacters","2D",function(rank)
  local symbols,uc,n,i,d,s,r,f,z,Defect0to2;
  uc:=rec(harishChandra:=[],charSymbols:=[],almostHarishChandra:=[]);
  for d in 4*[0..QuoInt(RootInt(rank)-1,2)]+2 do
    r:=d^2/4;
    s:=rec(relativeType:=rec(series:="B",indices:=[1+r..rank],rank:=rank-r),
           levi:=[1..r],
           eigenvalue:=1, # see Geck-malle
           parameterExponents:=Concatenation([d],[2+r..rank]*0+1));
    if r<10 then s.cuspidalName:=SPrint("{}^2D_",r,"");
    else s.cuspidalName:=SPrint("{}^2D_{",r,"}");
    fi;
    Add(uc.harishChandra,s);
    if d=2 then s.levi:=[]; s.cuspidalName:="";fi;
    symbols:=Symbols(rank,d);
    s.charNumbers:=[1..Length(symbols)]+Length(uc.charSymbols);
    FixRelativeType(s);
    Append(uc.charSymbols,symbols);
  od;
  uc.a:=List(uc.charSymbols,LowestPowerGenericDegreeSymbol);
  uc.A:=List(uc.charSymbols,HighestPowerGenericDegreeSymbol);
  uc.almostCharSymbols:=[];
  for d in 4*[0..RootInt(QuoInt(rank,4),2)] do
    r:=d^2/4;
    s:=rec(relativeType:=rec(series:="B",indices:=[1+r..rank],rank:=rank-r),
           levi:=[1..r], eigenvalue:=(-1)^QuoInt(d+1,4));
    if r<10 then s.cuspidalName:=SPrint("D_",r,"");
    else s.cuspidalName:=SPrint("D_{",r,"}");
    fi;
    r:=s.relativeType.rank;
    symbols:=Symbols(rank,d);
    if QuoInt(d+1,4) mod 2<>0 then symbols:=List(symbols,Reversed);fi;
    if d=0 then
      s.relativeType.series:="D";
      s.relativeType:=rec(orbit:=[s.relativeType],twist:=(1,2));
      s.cuspidalName:="";
      symbols:=List(CHEVIE.R("CharParams","2D")(rank),
                    x->SymbolPartitionTuple(x,0));
    fi;
#the map which goes from almost characters to unipotent characters for 2Dn
    Defect0to2:=function(ST)local a;
      a:=Minimum(SymmetricDifference(ST[1],ST[2]));
      ST:=[SymmetricDifference(ST[1],[a]),SymmetricDifference(ST[2],[a])];
      if Length(ST[1])>Length(ST[2]) then return ST;else return Reversed(ST);fi;
    end;
    s.charNumbers:=List(symbols,s->Position(uc.charSymbols,Defect0to2(s)));
    uc.almostCharSymbols{s.charNumbers}:=symbols;
    if d<>0 then FixRelativeType(s);fi;
    Add(uc.almostHarishChandra,s);
  od;
  # note: delta is always 1 since a+A is always even
  z:=x->rec(Z1:=SymmetricDifference(x[1],x[2]),Z2:=Intersection(x));
  uc.families:=List(Set(List(uc.charSymbols,z)),function(f)local sharp,res;
    sharp:=s->SymmetricDifference(Difference(s[2],f.Z2),
                                  f.Z1{[1,3..Length(f.Z1)-1]});
    res:=rec(charNumbers:=Filtered([1..Length(uc.charSymbols)],
         i->z(uc.charSymbols[i])=f));
    res.almostCharNumbers:=res.charNumbers;
    res.fourierMat:=List(uc.charSymbols{res.charNumbers},
       u->List(uc.almostCharSymbols{res.almostCharNumbers},
       a->2^(-QuoInt(Length(f.Z1)-1,2))*
         (-1)^Length(Intersection(sharp(u),sharp(a)))));
    if Length(res.fourierMat)=16 then # JM jan 2015: fix this horrible kludge
      res.fourierMat[16]:=-res.fourierMat[16];
      res.fourierMat{[1..16]}[16]:=-res.fourierMat{[1..16]}[16];
    fi;
    res.eigenvalues:=res.charNumbers*0+1; # see Geck-Malle
    res.sh:=List(res.charNumbers,y->1);  # is that correct for Geck-Malle?
    if Length(res.eigenvalues)=1 then res.charLabels:=[""];res.special:=1;
    else res.charLabels:=List(uc.charSymbols{res.charNumbers},
      function(M)local v,D,v1,v2,s;
      M:=SymmetricDifference(Difference(M[2],f.Z2),f.Z1{[3,5..Length(f.Z1)-1]});
      v:=List(f.Z1,z->Number(M,y->y>=z)mod 2);
      D:=Length(v);
      v1:=v{[2,4..D-(D mod 2)]};
      v2:=v{[3,5..D-1+(D mod 2)]};
      if D mod 2=1 then Add(v1,0);fi;
      # v1, v2 is coordinates in (e1,e3,e5,..) and in (e2,e4,..) basis
      v1:=List([1..Length(v2)],i->Sum(v1{[i,i+1]})mod 2);
      # coordinates in e1, e1+e3, e1+e3+e5, ...
      s:="+-";
      return ConcatenationString(s{v2+1},",",s{v1+1});end);
    fi;
    res.special:=PositionProperty(res.charLabels,x->ForAll(x,y->y in "+,"));
    res.name:=Concatenation(f.Z1,f.Z2,f.Z2);Sort(res.name);
    res.name:=IntListToString(res.name);
    res.explanation:="classical family";
    res.perm:=();
    res.size:=Length(res.charNumbers);
    res.operations:=FamilyOps;
    return res;
  end);
  return uc;
end);

CHEVIE.AddData("UnipotentClasses","2D",function(r,p)
  return CHEVIE.R("UnipotentClasses","D")(r,p);end);
