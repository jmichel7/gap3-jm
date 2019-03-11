#############################################################################
##
#A  tbl/compat4.g                 CHEVIE library                 Frank Lübeck
##
#Y  Copyright (C) 2001  The CHEVIE Team
##
##  This file  contains some functions  for GAP 4  which allow to  use the
##  remaining data files to  be shared by the GAP 3 and  GAP 4 versions of
##  CHEVIE.
##
# maybe change these names in GAP3 ? XXX
# (or define synonyms in GAP3 and use the new names in data files)
SymbolPartitionTuple := SymbolPartitionTupleShape;
Symbols := SymbolsRankDefect;
PermCosetsSubgroup := PermsCosetsReflectionSubgroup;
PartBeta := PartitionBetaSet;
Hecke := HeckeAlgebra;
Format := String;

#############################################################################
##
##  CHEVIE.compat  . . . . . . . . . . . . record to collect functions which
##  are different for GAP 3 and GAP 4
##
CHEVIE.compat := rec();
##  Don't use this for actual `Info' statements, only used in data files
##  which should be GAP 3 compatible.
CHEVIE.compat.InfoChevie := function(arg)
  if InfoLevel(InfoChevie) > 0 then
    CallFuncList(Print, arg);
  fi;
end;
##  For copying components of a record / component object
CHEVIE.compat.Inherit := function(c1, c2)
  local a;
  for a in NamesOfComponents(c2) do
    c1!.(a) := c2!.(a);
  od;
end;

##  just a new name
CHEVIE.compat.CharTable := CharacterTable;
##  the `!.' syntax is not readable by GAP 3
CHEVIE.compat.AdjustHeckeCharTable := function(tbl, param)
  local r, irr, i;
  r := [1..Length(tbl!.WordsClassRepresentatives)];
  irr := List(Irr(tbl), ShallowCopy);
  for i  in r  do
    irr{r}[i] := irr{r}[i]
           * Product( (-1 * param{tbl!.WordsClassRepresentatives[i]}[2]) );
  od;
  tbl!.Irr := List(irr, a-> Character(tbl, a));
end;

# in GAP 4 we may have to translate first from GAP 3 version
CHEVIE.compat.MakeCharacterTable := function(tbl)
  if IsRecord(tbl) then
    tbl := CHEVIE.TranslateCharTableGAP4GAP3(tbl);
  fi;
  return tbl;
end;

# again `!.' syntax
CHEVIE.compat.ChangeIdentifier := function(obj, name)
  obj!.Identifier := name;
end;

##  change name in GAP3 or GAP4?
PartitionTupleToString := IntListTupleToString;

##  The remaining entries are for character tables in classical types. Their
##  creation is very incompatible between GAP3/4.
CHEVIE.compat.CharTableA := function(rank)
  local   tbl,  p;
  tbl := CharacterTable("Symmetric", rank+1);
  p := Immutable(Partitions(rank+1));
  ##  here we overwrite existing attribute values!
  tbl!.CharacterParameters := p;
  tbl!.ClassParameters := p;
  ##  a bit smoother in GAP 4
  p := List(p, IntListToString);
  Setter(CharacterNames)(tbl, p);
  Setter(ClassNames)(tbl, p);
  return tbl;
end;
CHEVIE.compat.HeckeCharTableA := function(n, param, sqrtparam)
  local q, tbl, r, irr, i;
  q:=-param[1][1]/param[1][2];
  tbl:= CharacterTableSpecialized(CHEVIE.R("Hk","A"), [n+1, q]);
  Setter(CartanMatAttr)(tbl, CartanMat("A", n));
  Setter(Parameters)(tbl, List([1..n], x-> q));
  tbl!.WordsClassRepresentatives:=CHEVIE.R("ClassInfo","A")(n).classtext;
  tbl!.ClassParameters := List(tbl!.ClassParameters, x->x[2]);
  Setter(ClassNames)(tbl, List(tbl!.ClassParameters, IntListToString));
  tbl!.CharacterParameters := CHEVIE.Data("CharParams", "A", n);
  Setter(CharacterNames)(tbl, List(CharacterParameters(tbl),
                                   CHEVIE.R("CharName","A")));
  CHEVIE.compat.AdjustHeckeCharTable(tbl, param);
  return tbl;
end;
##  case ^2A
CHEVIE.compat.CharTable2A := function(r)
  local aphi, tbl, irr, cp, eps, i;

  # [LuB, 4.19] Preferred extension: -1 acts on \tilde E where E corresponds
  # to partition p by (-1)^aphi(p)
  aphi := function(p)
    return Sum(List(AssociatedPartition(p),i->i*(i-1)/2));
  end;

  tbl := CHEVIE.R("CharTable","A")(r);
  tbl!.Identifier := String(Concatenation("W(2A",String(r),")"));
  irr := ShallowCopy(Irr(tbl));
  cp := CharacterParameters(tbl);
  for i in [1..Length(irr)] do
    eps := (-1)^aphi(cp[i]);
    if eps = -1 then
      irr[i] := -irr[i];
    fi;
  od;
  tbl!.Irr := Immutable(irr);

  # delete computed power maps
  tbl!.ComputedPowerMaps := [];

  return tbl;
end;
CHEVIE.compat.HeckeCharTable2A := function(r,param,rootparam)
  local q, v, W, phi, WF, qE, H, tbl, w0, cl, as, i;
  q:=-param[1][1]/param[1][2];
  if not IsBound(rootparam[1]) then v:=GetRoot(q,2);
  else v:=rootparam[1];
  fi;
  W:=CoxeterGroupByReflectionDatum("A",r);
  if r=1 then phi:=();
  else phi:=Product([1..QuoInt(r,2)],i->(i,r+1-i));
  fi;
  WF:=CoxeterCoset(W,phi);
# If q_E is the square root which deforms to 1 of the eigenvalue of T_{w_0}
# on E which deforms to 1, then we have:
#  E~(T_w\phi)=\overline(E(T_{w^-1w_0}))q_E (trivial extension)
#  E~(T_w\phi)=(-1)^a_E\overline(E(T_{w^-1w_0}))q_E (preferred extension)
# where \overline means q->q^-1
  qE:=HeckeCentralMonomials(HeckeAlgebra(W,v));
  H:=HeckeAlgebra(W,v^-2);
  tbl:=CharacterTable(H);
  tbl!.WordsClassRepresentatives :=
                        CHEVIE.R("WordsClassRepresentatives","2A")(r);
  tbl!.Identifier := SPrint("H(^2A",r,")");
  w0 := LongestCoxeterElement(W);
  cl := List(tbl!.WordsClassRepresentatives, x->ElementCoxeterWord(W,x) * w0);
  tbl!.Irr := TransposedMat(List(cl, x->
                                   HeckeCharValues(HeckeElement(H,"T",x))));
  as := List(CharacterParameters(tbl), CHEVIE.R("LowestPowerFakeDegree","A"));
  for i in [1..Length(tbl!.Irr)] do
    tbl!.Irr[i] := (-1)^as[i] * qE[i] * tbl!.Irr[i];
  od;
  CHEVIE.compat.AdjustHeckeCharTable(tbl, param);
  return tbl;
end;

# for tname in ["B", "Bsym", "C"]
CHEVIE.compat.CharTableB := function(tname)
  return function(rank)
    local   tbl,  p;
    tbl := CharacterTable("WeylB", rank);
    # overwrite some existing attributes
    tbl!.Identifier := SPrint("W(",tname,rank,")");
    p := Immutable(CHEVIE.Data("CharParams", "B", rank));
    tbl!.CharacterParameters := p;
    tbl!.ClassParameters := p;
    p := List(p, CHEVIE.R("CharName", "B"));
    Setter(CharacterNames)(tbl, p);
    Setter(ClassNames)(tbl, p);
    return tbl;
  end;
end;
CHEVIE.compat.HeckeCharTableB := function(tname)
  return function(n, param, sqrtparam)
    local q, tbl, r, irr, i;
    q := [-param[1][1]/param[1][2], -param[2][1]/param[2][2]];
    tbl := CharacterTableSpecialized(CHEVIE.R("Hk","B"), [n, q[1], q[2]]);
    tbl!.Identifier := SPrint("H(", tname, n, ")");
    Setter(CartanMatAttr)(tbl, CartanMat(tname, n));
    Setter(Parameters)(tbl, List([1..n], x-> q));
    tbl!.WordsClassRepresentatives :=
                        CHEVIE.R("WordsClassRepresentatives", tname)(n);
    tbl!.ClassParameters := List(tbl!.ClassParameters, x->x[2]);
    Setter(ClassNames)(tbl, List(tbl!.ClassParameters,
                                      CHEVIE.R("ClassName", tname)));
    tbl!.CharacterParameters := CHEVIE.Data("CharParams", tname, n);
    Setter(CharacterNames)(tbl, List(CharacterParameters(tbl),
                                     CHEVIE.R("CharName",tname)));
    CHEVIE.compat.AdjustHeckeCharTable(tbl, param);
    return tbl;
  end;
end;
CHEVIE.compat.CharTableD := function(l)
  local hi, erg, i, j, lst, lst2, p, sy, syc, tmp, cli, pow,  S, Sval;

  hi:=CharacterTable("WeylB",l);
  cli:=CHEVIE.R("ClassInfo","D")(l);
  if l mod 2 = 0 then
    S:=CharacterTable("Symmetric",l/2);
    Sval:=function(aa,pp)
      return Irr(S)[Position(ClassParameters(S), [1,aa[1]])]
             [Position(ClassParameters(S), [1,pp[1]/2])];
    end;
  fi;

  # classes in subgroup:
  lst:=[];
  for i in [1..Length(ClassParameters(hi))] do
    if Length(ClassParameters(hi)[i][2][2]) mod 2=0 then
      # degenerate classes:
      if ClassParameters(hi)[i][2][2] = [] and
                  ForAll(ClassParameters(hi)[i][2][1]/2,IsInt) then
        Append(lst,[i,i]);
      else
        Add(lst,i);
      fi;
    fi;
  od;

  # similarly for characters:
  lst2:=[];
  for i in [1..NrConjugacyClasses(hi)] do
    sy:= CharacterParameters(hi)[i][2];
    if sy[1]<=sy[2] then
      # degenerate characters:
      if sy[1]=sy[2] then
        Append(lst2,[i,i]);
      else
        Add(lst2,i);
      fi;
    fi;
  od;

  # for power maps:
  pow:=function(i,p,lst)
    local r;
    r:=Position(lst,p[lst[i]]);
    if '-' in cli.classparams[i] and '+' in cli.classparams[r] then
      r:=r+1;
    fi;
    return r;
  end;

  erg:=rec(Identifier:=String(Concatenation("W(D",String(l),")")));
  erg.UnderlyingCharacteristic := 0;
  erg.Size:=Size(hi)/2;
  erg.SizesCentralizers:=SizesCentralizers(hi){lst}/2;
  erg.OrdersClassRepresntatives := OrdersClassRepresentatives(hi){lst};
  erg.ClassParameters:=cli.classparams;
  erg.ClassNames:=cli.classnames;
  erg.SizesConjugacyClasses:= SizesConjugacyClasses(hi){lst};
  erg.ComputedPowerMaps:=[];
  for i in [1..Length(ComputedPowerMaps(hi))] do
    if IsBound(ComputedPowerMaps(hi)[i]) then
      erg.ComputedPowerMaps[i]:=List([1..Length(lst)],j->
                                     pow(j,ComputedPowerMaps(hi)[i],lst));
    fi;
  od;
  erg.InfoText:="extracted from generic character table of type WeylB";
  erg.CharacterParameters := List(CharacterParameters(hi){lst2}, a->a[2]);
  erg.CharacterNames := CharacterNames(hi){lst2};
  erg.Irr:=Irr(hi){lst2}{lst};

  # corrections for degenerate classes and characters:
  for i in [1..Length(lst)] do
    if '+' in erg.ClassParameters[i] or '-' in erg.ClassParameters[i] then
      erg.SizesCentralizers[i]:=2*erg.SizesCentralizers[i];
      erg.SizesConjugacyClasses[i]:= erg.SizesConjugacyClasses[i]/2;
    fi;
  od;

  i:=1;
  while i<=Length(lst2) do
    sy:= erg.CharacterParameters[i];
    if sy[1]=sy[2] then
      erg.Irr[i]:= erg.Irr[i]/2;
      erg.Irr[i+1]:=erg.Irr[i+1]/2;
      erg.CharacterParameters[i] := [erg.CharacterParameters[i][1], '+'];
      erg.CharacterParameters[i+1] := [erg.CharacterParameters[i+1][1], '-'];
      erg.CharacterNames[i]:= CHEVIE.R("CharName","D")
                              (erg.CharacterParameters[i]);
      erg.CharacterNames[i+1]:= CHEVIE.R("CharName","D")
                              (erg.CharacterParameters[i+1]);
      for j in [1..Length(lst)] do
        syc:=erg.ClassParameters[j];
        if '+' in syc then
          tmp:=2^(Length(syc[1])-1)*Sval(sy,syc);
          erg!.Irr[i][j]:=erg!.Irr[i][j]+tmp;
          erg!.Irr[i+1][j+1]:=erg!.Irr[i+1][j+1]+tmp;
          erg!.Irr[i][j+1]:=erg!.Irr[i][j+1]-tmp;
          erg!.Irr[i+1][j]:=erg!.Irr[i+1][j]-tmp;
        fi;
      od;
      i:=i+2;
    else
      i:=i+1;
    fi;
  od;
  erg.CharacterNames := List(erg.CharacterParameters,
                             CHEVIE.R("CharName","D"));
  ConvertToCharacterTableNC(erg);
  SetFilterObj(erg, Tester(ClassNames) and Tester(CharacterNames));
  return erg;
end;
CHEVIE.compat.HeckeCharTableD := function(n, param, sqrtparam)
  local q, tbl, r, irr, i;
  q := -param[1][1]/param[1][2];
  tbl := CharacterTableSpecialized(CHEVIE.R("Hk","D"), [n, q]);
  tbl!.Identifier := SPrint("H(D", n, ")");
  Setter(CartanMatAttr)(tbl, CartanMat("D", n));
  Setter(Parameters)(tbl, List([1..n], x-> q));
  tbl!.WordsClassRepresentatives :=
                          CHEVIE.R("WordsClassRepresentatives", "D")(n);
  tbl!.ClassParameters := List(tbl!.ClassParameters, x->x[2]);
  Setter(ClassNames)(tbl, List(tbl!.ClassParameters,
                                    CHEVIE.R("ClassName", "D")));
  tbl!.CharacterParameters := CHEVIE.Data("CharParams", "D", n);
  Setter(CharacterNames)(tbl, List(CharacterParameters(tbl),
                                   CHEVIE.R("CharName","D")));
  CHEVIE.compat.AdjustHeckeCharTable(tbl, param);
  return tbl;
end;

CHEVIE.compat.CharTable2D := function(l)
  local hi, tbl, i, lst, lst2, p, sy, cli, test, cp;

  hi := CharacterTable("WeylB", l);
  cli := CHEVIE.R("ClassInfo","2D")(l);
  cp := ClassParameters(hi);
  lst := Filtered([1..Length(cp)], i->Length(cp[i][2][2]) mod 2=1);
  tbl := rec(identifier := String(Concatenation("W(2D",String(l),")")));
  tbl.size := Size(hi)/2;
  tbl.centralizers := SizesCentralizers(hi){lst}/2;
  tbl.orders := OrdersClassRepresentatives(hi){lst};
  tbl.classparam := cli.classparams;
  tbl.classtext := cli.classtext;
  tbl.classnames := cli.classnames;
  tbl.classes := SizesConjugacyClasses(hi){lst};
  tbl.text := "extracted from generic character table of type B";
  sy := List(CharacterParameters(hi), a->a[2]);
  lst2 := Filtered([1..Length(sy)],i-> CHEVIE.R("testchar","2D")(sy[i]));
  tbl.irredinfo := List(CharacterParameters(hi){lst2}, a->rec(charparam := a[2],
                 charname := CHEVIE.R("CharName","2D")(a[2])));
  tbl.irreducibles := Irr(hi){lst2}{lst};
  tbl := CHEVIE.TranslateCharTableGAP4GAP3(tbl);
  return tbl;
end;
CHEVIE.compat.HeckeCharTable2D := function(l,param,rootparam)
  local q, hi, clp, lst2, lst, tbl, cli, tchar, chp, H, Replace;
  q := -param[1][1] / param[1][2];
  q := Concatenation([[q^0,-1]], List([2..l],i->[q,-1]));
  hi := CHEVIE.R("HeckeCharTable","B")(l, q, []);
  clp := ClassParameters(hi);
  lst2 := [1..Length(clp)];
  lst := Filtered(lst2, i-> Length(clp[i][2]) mod 2 = 1);
  tbl := rec(Identifier := SPrint("H(^2D",l,")"),
           Size := Size(hi)/2,
           OrdersClassRepresentatives := OrdersClassRepresentatives(hi){lst},
           SizesCentralizers := SizesCentralizers(hi){lst}/2,
	   SizesConjugacyClasses := SizesConjugacyClasses(hi){lst},
           InfoText := "extracted from generic character table of HeckeB");
  cli := CHEVIE.R("ClassInfo","2D")(l);
  tbl.WordsClassRepresentatives := cli.classtext;
  tbl.ClassNames := cli.classnames;
  tbl.ClassParameters := cli.classparams;

  tchar := CHEVIE.R("testchar","2D");
  chp := CharacterParameters(hi);
  lst2 := Filtered(lst2, i-> tchar(chp[i]));
  tbl.CharacterParameters := chp{lst2};
  tbl.CharacterNames := List(tbl.CharacterParameters,
                                          CHEVIE.R("CharName","2D"));
  H := HeckeAlgebra(CoxeterGroupByReflectionDatum("B",l), q);
  Replace := SubstitutionSublist;
  tbl.Irr := TransposedMat(List(tbl.WordsClassRepresentatives,
    x->HeckeCharValues(Basis(H,"T")
      (Concatenation([1],Replace(x,[1],[1,2,1]))),
      List(Irr(hi), AsList){lst2})));
  tbl.UnderlyingCharacteristic := 0;
  ConvertToCharacterTableNC(tbl);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end;
