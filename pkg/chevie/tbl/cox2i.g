#############################################################################
##
#A  tbl/cox2i.g                CHEVIE library                    Meinolf Geck
##
#Y  Copyright (C) 1992 - 2001  The CHEVIE Team
##
##  This file contains information about the coset W.F where W is
##  the group I_2(m), m>3:   1-2 and  F is (1,2).
##
##  The HeckeCharTable was obtained by explicit computation using
##  representations of the Hecke algebra of type I_2(m).
##
CHEVIE.AddData("ReflectionName","2I",
  function(m,option)
    if IsBound(option.TeX) then return SPrint("{}^2I_2(",m,")");
    else return SPrint("2I2(",m,")");
    fi;end);

CHEVIE.AddData("NrConjugacyClasses","2I",m->QuoInt(m+3,2));

CHEVIE.AddData("WordsClassRepresentatives","2I", function(m)local r, x, i;
  r:=[[]]; x:=[1];
  for i in [1..QuoInt(m+1,2)] do Add(r,ShallowCopy(x)); Append(x,[2,1]); od;
  return r;
end);

CHEVIE.AddData("ClassInfo","2I",function(m)local res;
  res:=rec(classtext:=CHEVIE.R("WordsClassRepresentatives", "2I")(m));
  res.classnames:=List(res.classtext,IntListToString);
  res.classparams:=res.classnames;
  res.classes:=[m];Append(res.classes,[1..QuoInt(m,2)]*0+2);
  if m mod 2=1 then Add(res.classes,1);fi;
  res.orders:=List(res.classtext,i->2*m/Gcd(2*m,Length(i))); res.orders[1]:=2;
  return res;
end);

CHEVIE.AddData("PhiFactors","2I",m->[1,-1]);

#############################################################################
##
#F  CHEVIE.R("ClassParameter","2I")( <w> )   . class params
##
##  given an element w  of a Coxeter group W of type I_2 as word in standard
##  generators, returns the classparam of its F-conjugacy class under the
##  nontrivial F action permuting the generators  1<->2.
##
CHEVIE.AddData("ClassParameter","2I",function(m,w)local l;
  if IsInt(Length(w)/2) then
    return "";
  else
    l:=ShallowCopy(w); if l[1]=2 then l:=l{[2..Length(l)]};Add(l,1); fi;
    return IntListToString(l);
  fi;
end);

#how to make a .charname from a .charparam
CHEVIE.AddData("CharName","2I",function(m,x,option)local s;
  if IsList(x[1]) then return PartitionTupleToString(x);
  else
    if IsBound(option.TeX) then s:="phi";else s:="\\phi_";fi;
    s:=SPrint(s,"{",x[1],",",x[2],"}");
    if Length(x)=3 then Append(s,x[3]);fi;
    return String(s);
  fi;
end);

CHEVIE.AddData("CharInfo","2I",function(m)local res;
  res:=rec(extRefl:=[1,3,2]);
  if m=4 then res.charparams:= [[[2],[]], [[],[1,1]],[[1],[1]]];
  else res.charparams:=Concatenation([[1,0],[1,m]],
      List([1..QuoInt(m-1,2)],i->[2,i]));
  fi;
  return res;
end);

#########################################################################
##
#F  HeckeCharTable( <m>, <v> )  . . . . . . . . . . . . . . . .
##  . . . . . . . . character table of the Hecke algebra of dihedral type
##  The parameter of the algebra is  <v>^2; we assume that <m> > 3, since
##  the case <m>=3 is covered by 2A_2. But we also include the case where
##  <m>=4 (2B_2) or <m>=6 (2G_2).
##
CHEVIE.AddData("HeckeCharTable","2I",function(m,param,sqrtparam)
  local q,i,j,ct,cos,cl,l,ident,ord,v,tbl;
  q:=-param[1][1]/param[1][2];
  if m=4 then ident:="2B"; elif m=6 then ident:="2G"; else ident:="2I2"; fi;
  ident:=SPrint(ident,"(",m,")");
  if q<>1 then ident:=SPrint("H(",ident,")");fi;
  if not IsBound(sqrtparam[1]) then v:=GetRoot(q,2,"CharTable(",ident,")");
  else v:=sqrtparam[1];
  fi;
  cl:=CHEVIE.R("ClassInfo","2I")(m);
  cos:=i->E(2*m)^i+E(2*m)^-i;
  ct:=[List(cl.classtext,i->q^Length(i)),List(cl.classtext,i->(-1)^Length(i))];
  for i in [1..QuoInt(m-1,2)] do
    ct[i+2]:=[0];
    for j in [1..QuoInt(m+1,2)] do ct[i+2][j+1]:=-v^(2*j-1)*cos(i*(2*j-1));
    od;
#   ct[i+2]:=ct[i+2]*(-1)^i; # to make Ennola duality work
  od;
  tbl:=rec(identifier:=ident, name:=ident,
     cartan:=[[2,-cos(1)],[-cos(1),2]],  size:=2*m,
     parameter:=[q,q],  sqrtparameter:=[v,v],
     irreducibles:=ct*v^0,
     irredinfo:=List(CHEVIE.R("CharInfo","2I")(m).charparams,x->rec(
      charparam:=x,charname:=CHEVIE.R("CharName","2I")(m,x,rec(TeX:=true)))));
  Inherit(tbl,cl);
  tbl.centralizers:=List(cl.classes,i->tbl.size/i);
  tbl := CHEVIE.compat.MakeCharacterTable(tbl);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end);

CHEVIE.AddData("CharTable","2I",
 m->CHEVIE.R("HeckeCharTable","2I")(m,[[1,-1],[1,-1]],[1,1]));

CHEVIE.AddData("Representation","2I",function(m,i)
  return CHEVIE.R("HeckeRepresentation","2I")
 (m,[[1,-1],[1,-1]],[1,1],i);end);

CHEVIE.AddData("HeckeRepresentation","2I",
  function(m,param,sqrtparam,i)local v,q,e;
  q:=-param[1][1]/param[1][2];
  if not IsBound(sqrtparam[1]) then
       v:=GetRoot(q,2,"Representation(Hecke(2I2(",m,")),[",i,"])");
  else v:=sqrtparam[1];
  fi;
  e:=E(2*m);
  if i=1 then return rec(gens:=[[[v^2]],[[v^2]]], F:=[[1]]);
  elif i=2 then return rec(gens:=[[[-1]],[[-1]]]*v^0, F:=[[1]]);
  else i:=i-2;return rec(gens:=[[[-1,0],[v*(e^i+e^-i),v^2]],
                [[v^2,v*(e^i+e^-i)],[0,-1]]]*v^0, F:=-[[0,1],[1,0]]);
  fi;
end);

CHEVIE.AddData("UnipotentCharacters","2I",function(e)
  local nc,uc,i,ac,c,n,symUnp,untUnp,TeXpref,eig;
##############################################################################
##  Symbols, unipotent degrees, Fourier matrix, Frobenius eigenvalues       ##
##  as in: G. Malle, "Unipotente Grade...", J. Algebra 177 (1995)           ##
##  and: M. Geck, G.M. "Fourier transforms...", J. Algebra 260 (2003)       ##
##  for non-trivial family of twisted dihedral group {^e}G(e,e,2)           ##
##                                                               GM 05.12.01##
##############################################################################
# There are 3 families: Is, St and a big one M.
#  M has principal series almost characters R(0,j) where 0<j<e/2
#    with fake degree q^{e-j}-q^j.
#  M has cuspidal almost characters R(k,l) where 0<k<l<e-k
#  M has unipotent characters indexed by S(k,l) where 0<k<l<2e-k with k,l odd.
# Malle associates a symbol (s_0,...,s_{e-1}) where
#   s_i=[0] for i\ne 0,k,l,k+l, s_0=s_{k+l}=[0,1] and s_k=s_l=[]
# The principal series chars in M are S(0,l) where 0<l<e/2 with 
#   associated symbol s_i=[0] for i\ne 0,l and s_0=s_l=[1].
# The b is min(k+l,e-k-l) and the eigenvalue of F is E(e)^(k*l).
# Fourier is as given below.
  uc:=rec();
  n:=QuoInt(e-1,2);
  nc:=Concatenation(List([1..n],i->List([i+1..e-i-1],j->[i,j])));
  ac:=Concatenation(List([1..n],l->[0,l]),nc); # almost-char symbols
  symUnp:=Concatenation(List([1..n],i->List([i+1..e-i],j->2*[i,j]-1)));
  if e mod 2=1 then
  # each character is aligned with its Ennola-correspondent almost char.
    untUnp:=List(symUnp,function(s)local res;
    res:=List((e-Reversed(s))/2,x->x mod e);
    if res[1]>res[2] then res:=List(-res,x->x mod e);fi;
    return res;
    end);
    SortParallel(untUnp,symUnp);
  fi;
  if e=4 then TeXpref:="B_2"; elif e=6 then TeXpref:="G_2";
  else TeXpref:=SPrint("I_2(",e,")"); fi;
  uc.harishChandra:=Concatenation(
   [ rec(relativeType:=rec(series:="A",indices:=[1],rank:=1),
         parameterExponents:=[e], levi:=[], eigenvalue:=1,
	 cuspidalName:="", charNumbers:=[2,1])],
   List(symUnp,x->
     rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
     parameterExponents:=[], levi:=[1,2], eigenvalue:=E(2*e)^(Product(x)),
     cuspidalName:=SPrint("{}^2",TeXpref,"[",x[1],",",x[2],"]"),
     charNumbers:=[2+Position(symUnp,x)])));
  uc.almostHarishChandra:=Concatenation(
   [ rec(relativeType:=rec(
      orbit:=[rec(series:="I",indices:=[1,2],rank:=2,bond:=e)],
      twist:=(1,2)), parameterExponents:=[1,1], levi:=[], eigenvalue:=1,
	 cuspidalName:="",charNumbers:=[1..n+2])],
   List(nc,x->
     rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
         parameterExponents:=[], levi:=[1,2], eigenvalue:=E(e)^(-Product(x)),
	 cuspidalName:=SPrint(TeXpref,"[",x[1],",",x[2],"]"),
         charNumbers:=[n+2+Position(nc,x)])));
  if e=4 then uc.almostHarishChandra[2].cuspidalName:="B_2";
              uc.almostHarishChandra[1].relativeType.orbit[1]:=
                 rec(series:="B",indices:=[1,2],rank:=2,cartanType:=ER(2));
  elif e=6 then eig:=[E(3)^2,-1,E(3),1];
    uc.almostHarishChandra[1].relativeType.orbit[1]:=
       rec(series:="G",indices:=[1,2],rank:=2,cartanType:=ER(3));
    for i in [1..4] do 
      uc.almostHarishChandra[i+1].cuspidalName:=
        SPrint("G2[",Format(eig[i],rec(TeX:=1)),"]");
    od;
  fi;
  uc.charParams:=Concatenation(CHEVIE.R("CharInfo","I")(e).charparams{[1,2]},ac);
  uc.almostCharSymbols:=Concatenation([List([1..e],x->[0]),
   List([1..e],x->[0,1])],List(ac,function(s)local S,k,l;
      S:=List([1..e],i->[0]);k:=s[1];l:=s[2];
      S[k+1]:=[];S[l+1]:=[];Add(S[1],1);Add(S[k+l+1],1);return S; end));
  uc.almostCharSymbols[1][e]:=[2];uc.almostCharSymbols[2][e]:=[1,2];
  uc.charSymbols:=Concatenation([List([1..e],x->[0]),List([1..e],x->[0,1])],
    List(symUnp,function(s)local S,k,l;
    k:=(s[1]+1)/2;l:=(s[2]+1)/2;S:=List([1..e],function(i)
      if i=k+1 or i=l+1 then return [];else return [0];fi;end);
      Add(S[1],1);Add(S[k+l],1);return S; end));
  uc.charSymbols[1]{[1,2]}:=[[0,2],[]];uc.charSymbols[2]{[1,2]}:=[[0,1,2],[1]];
  c:=a->E(2*e)^a+E(2*e)^(-a);
  uc.families:=[Family("C1",[1]), Family("C1",[2]),
    Family(rec(eigenvalues:=List(symUnp,s->E(2*e)^(s[1]*s[2])),
      fourierMat:=List(ac,j->List(symUnp,i->(c(i*Reversed(j))-c(i*j))/e)),
      sh:=List(ac,s->E(e)^(-s[1]*s[2])),
      charNumbers:=2+[1..Length(ac)],
      special:=1))];
  uc.a:=Concatenation([0,e],List(ac,x->1));
  uc.A:=Concatenation([0,e],List(ac,x->e-1));
# if e=4 then uc.families[3].fourierMat:=ER(2)/2*[[1,1],[1,-1]];fi;
  if e=5 then 
# Modified 25-8-2004 to fit with H3, H4
# Asterisque, Geck-Malle, H4 in Duke are like the old version
# 'Unipotente Grade' and I2, imprimitive, current H3 H4 agree with new version
    uc.families[3]:=uc.families[3]^13; # "GaloisCyc(f,13)"
    for c in uc.harishChandra do c.eigenvalue:=GaloisCyc(c.eigenvalue,13);od;
  fi;
  return uc;
# Properties: S^-1=TransposedMat(S)
# f.fakdeg:=Concatenation(List(ac,i->x^(e-i[2])-x^i[2]),List(nc,i->0*x));
# f.unpdeg:=List(symUnp,i->
#      (c(i[1])-c(i[2]))/e*x*(x^2-1)*(x^e+1)/Product(i,a->(x-z^a)*(x-z^-a)));
# S*f.unpdeg=f.fakdeg
# (DiagonalMat(f.eigenvalues)*TransposedMat(S)*DiagonalMat(f.sh)^-1*S)^2=S^0
end);
