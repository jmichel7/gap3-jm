#############################################################################
##
#A  tbl/compat3.g                 CHEVIE library                 Frank Lübeck
##
#Y  Copyright (C) 2001  The CHEVIE Team
##
##  This file  contains some functions  for GAP 3  which allow to  use the
##  remaining data files to  be shared by the GAP 3 and  GAP 4 versions of
##  CHEVIE.
##
if not IsBound(CHEVIE.Revision) then
  CHEVIE.Revision := Ignore;
fi;

## just another name in GAP 4
#CallFuncList := ApplyFunc;

#############################################################################
##
##  CHEVIE.compat  . . . . . . . . . . . . record to collect functions which
##  are different for GAP 3 and GAP 4
##
CHEVIE.compat := rec();

CHEVIE.compat.InfoChevie := function(arg)
  # looks circumstantially, but works after resetting InfoChevie
  ApplyFunc(InfoChevie, arg);
end;

##  new name in GAP 4
CHEVIE.compat.CharTable := CharTable;

##  very different for GAP 4
CHEVIE.compat.AdjustHeckeCharTable:=function(tbl,param)local r,i;
  r:=[1..Length(tbl.classtext)];
  for i in r do tbl.irreducibles{r}[i]:=tbl.irreducibles{r}[i]*
         Product(-param{tbl.classtext[i]}[2]);
  od;
end;

##  need `!.' in GAP 4
CHEVIE.compat.ChangeIdentifier := function(obj, name)
  obj.identifier := name;
end;

##  since function is not in GAP 4
CHEVIE.compat.Inherit := Inherit;

##  trivial here
CHEVIE.compat.MakeCharacterTable := function(tbl)
  tbl.operations:=CharTableOps;
  return tbl;
end;

##  The remaining entries are for character tables in classical types. Their
##  creation is very incompatible between GAP3/4.
CHEVIE.compat.CharTableA := function(n)local tbl;
  tbl:=CharTable("Symmetric",n+1);
  tbl.identifier:=SPrint("W(A",n,")");
  tbl.cartan:=CartanMat("A",n);
  Inherit(tbl,CHEVIE.R("ClassInfo","A")(n));
  tbl.irredinfo:=List(CHEVIE.R("CharParams","A")(n),
  x->rec(charparam:=x, charname:=CHEVIE.R("CharName","A")(n,x,rec(TeX:=true))));
  return tbl;
end;
CHEVIE.compat.HeckeCharTableA := function(n, param,sqrtparam)local tbl,q;
  q:=-param[1][1]/param[1][2];
  tbl:= CharTableSpecialized(CHEVIE.R("Hk","A"), [n+1, q]);
  tbl.cartan:= CartanMat("A", n);
  tbl.parameter:= List([1..n], x-> q);
  Inherit(tbl,CHEVIE.R("ClassInfo","A")(n));
  tbl.irredinfo:=List(CHEVIE.R("CharParams","A")(n),x->
    rec(charparam:=x,charname:=CHEVIE.R("CharName","A")(n,x,rec(TeX:=true))));
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end;
##  case ^2A
CHEVIE.compat.CharTable2A := function(r)local i, tbl;
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
end;
CHEVIE.compat.HeckeCharTable2A := function(r,param,rootparam)
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
end;
##  for tname in ["B", "Bsym", "C"]
CHEVIE.compat.CharTableB := function(tname)
  return function(rank)local tbl;
    tbl:=CharTable("WeylB",rank);
    tbl.identifier:=SPrint("W(",tname,rank,")");
    tbl.cartan:=CartanMat(tname,rank);
    Inherit(tbl,CHEVIE.R("ClassInfo","B")(rank));
    tbl.irredinfo:=List(CHEVIE.R("CharParams","B")(rank),
      x->rec(charparam:=x,
charname:=CHEVIE.R("CharName","B")(rank,x,rec(TeX:=true))));
    return tbl;
  end;
end;
CHEVIE.compat.HeckeCharTableB := function(tname)
  return function(n, param,sqrtparam) local q, tbl;
    q:=List(param,x->-x[1]/x[2]);
    tbl:= CharTableSpecialized(CHEVIE.R("Hk","B"), [n, q[1], q[2]]);
    tbl.identifier := SPrint("H(", tname, n, ")");
    tbl.cartan:= CartanMat(tname, n);
    tbl.parameter:= q;
    tbl.classtext:= CHEVIE.R("ClassInfo",tname)(n).classtext;
    tbl.classparams:= List(tbl.classparam, x-> x[2]);
    tbl.classnames:= List(tbl.classparams, PartitionTupleToString);
    tbl.irredinfo:=List(CHEVIE.R("CharParams",tname)(n),x->
             rec(charparam:=x,charname:=CHEVIE.R("CharName","B")(n,x,rec(TeX:=true))));
    CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
    return tbl;
  end;
end;
CHEVIE.compat.CharTableD := function(n)
  local hi, res, i, j, k, cls, chr, p, sy, syc, cli, pow,  S, Sval;

  hi:=CharTable("WeylB",n);
  if not IsBound(hi.size) then hi.size:=2^n*Factorial(n); fi;
  if n mod 2 = 0 then
    S:=CharTable("Symmetric",n/2);
    Sval:=function(aa,pp)
      return S.irreducibles[Position(S.classparam,[1,aa[1]])]
                           [Position(S.classparam,[1,pp[1]/2])];
    end;
  fi;

  # classes in subgroup:
  cls:=[];
  for i in [1..Length(hi.classparam)] do
    if Length(hi.classparam[i][2][2]) mod 2=0 then
      # degenerate classes:
      if hi.classparam[i][2][2] = [] and ForAll(hi.classparam[i][2][1]/2,IsInt)
         then
        Append(cls,[i,i]);
      else
        Add(cls,i);
      fi;
    fi;
  od;

  # similarly for characters:
  chr:=[];
  for i in [1..Length(hi.irredinfo)] do
    sy:=hi.irredinfo[i].charparam[2];
    if sy[1]<=sy[2] then
      # degenerate characters:
      if sy[1]=sy[2] then Append(chr,[i,i]); else Add(chr,i); fi;
    fi;
  od;

  cli:=CHEVIE.R("ClassInfo","D")(n);
  # for power maps:
  pow:=function(i,p,cls)local r; r:=Position(cls,p[cls[i]]);
    if '-' in cli.classparams[i] and '+' in cli.classparams[r] then
        r:=r+1;
    fi;
    return r;
  end;

  res:=rec(identifier:=SPrint("W(D",n,")"));
  res.size:=hi.size/2;
  res.cartan:=CartanMat("D",n);
  Inherit(res,cli);
  res.powermap:=[];
  for i in [1..Length(hi.powermap)] do
    if IsBound(hi.powermap[i]) then
      res.powermap[i]:=List([1..Length(cls)],j->pow(j,hi.powermap[i],cls));
    fi;
  od;
  res.text:="extracted from generic character table of type WeylB";
  res.operations:=CharTableOps;
  res.irredinfo:=List(CHEVIE.R("CharInfo","D")(n).charparams,p->
      rec(charparam:=p,charname:=PartitionTupleToString(p)));
  res.irreducibles:=hi.irreducibles{chr}{cls};

  i:=1;
  while i<=Length(chr) do
    sy:=res.irredinfo[i].charparam;
    if IsInt(sy[2]) then
      for k in [0,1] do
        res.irreducibles[i]:=res.irreducibles[i]/2;
	for j in [1..Length(cls)] do
	  syc:=res.classparams[j];
	  if '+' in syc then
	    res.irreducibles[i]{[j,j+1]}:=res.irreducibles[i]{[j,j+1]}
	      +(-1)^k*2^(Length(syc[1])-1)*Sval(sy,syc)*[1,-1];
	  fi;
	od;
        i:=i+1;
      od;
    else i:=i+1;
    fi;
  od;

  return res;
end;

CHEVIE.compat.HeckeCharTableD:=function(n, param,sqrtparam)local r, tbl,u;
   u:=-param[1][1]/param[1][2];
   tbl:= CharTableSpecialized(CHEVIE.R("Hk","D"), [n, u]);
   tbl.cartan:= CartanMat("D", n);
   tbl.parameter:= List([1..n], x-> u);
   tbl.irredinfo:=List(CHEVIE.R("CharInfo","D")(n).charparams,p->
      rec(charparam:=p,charname:=PartitionTupleToString(p)));
   Inherit(tbl,CHEVIE.R("ClassInfo","D")(n));
   CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
   return tbl;
end;

CHEVIE.compat.CharTable2D := function(l)local hi, tbl,lst, chr;
  hi:=CHEVIE.R("CharTable","B")(l);
  chr:=[1..Length(hi.classparams)];
  lst:=Filtered(chr, i->Length(hi.classparams[i][2]) mod 2=1);
  chr:=Filtered(chr,i->CHEVIE.R("testchar","2D")
                     (hi.irredinfo[i].charparam));
  tbl:=rec(identifier:=SPrint("W(^2D",l,")"),
    size:=hi.size/2,
    centralizers:=hi.centralizers{lst}/2,
    orders:=hi.orders{lst},
    classes:=hi.classes{lst},
    text:="extracted from generic character table of type B",
    operations:=CharTableOps,
    irredinfo:=List(hi.irredinfo{chr},a->rec(charparam:=a.charparam,
		   charname:=CHEVIE.R("CharName","2D")(l,a.charparam))),
    irreducibles:=hi.irreducibles{chr}{lst});
   Inherit(tbl,CHEVIE.R("ClassInfo","2D")(l));
   return tbl;
end;

CHEVIE.compat.HeckeCharTable2D := function(l,param,rootparam)
  local hi,cli,lst,tbl,chr,q;
  q:=-param[1][1]/param[1][2];
  q:=Concatenation([[q^0,-1]],List([2..l],i->[q,-1]));
  hi:=CHEVIE.R("HeckeCharTable","B")(l,q,[]);
  chr:=[1..Length(hi.classparams)];
  lst:=Filtered(chr,i->Length(hi.classparams[i][2]) mod 2=1);
  tbl:=rec(identifier:=SPrint("H(^2D",l,")"),
	   size:=hi.size/2,
	   orders:=hi.orders{lst},
	   centralizers:=hi.centralizers{lst}/2,
       classes:=hi.classes{lst},
	   text:="extracted from generic character table of HeckeB",
	   operations:=CharTableOps);
  Inherit(tbl,CHEVIE.R("ClassInfo","2D")(l));
  chr:=Filtered(chr,
	i->CHEVIE.R("testchar","2D")(hi.irredinfo[i].charparam));
  tbl.irredinfo:=List(hi.irredinfo{chr},a->rec(charparam:=a.charparam,
      charname:=CHEVIE.R("CharName","2D")(l,a.charparam)));
  tbl.irreducibles:=TransposedMat(List(tbl.classtext,
    x->HeckeCharValues(Basis(Hecke(CoxeterGroup("B",l),q),"T")
      (Concatenation([1],Replace(x,[1],[1,2,1]))),hi.irreducibles{chr})));
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end;
