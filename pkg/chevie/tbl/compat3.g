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
CHEVIE.compat.CharTableA := function(n)local tbl,ci;
  tbl:=CharTable("Symmetric",n+1);
  tbl.identifier:=SPrint("W(A",n,")");
  tbl.cartan:=CartanMat("A",n);
  Inherit(tbl,CHEVIE.R("ClassInfo","A")(n));
  ci:=CHEVIE.R("CharInfo","A")(n);
  tbl.irredinfo:=Zip(ci.charparams,ci.charnames,function(x,y)return
    rec(charparam:=x,charname:=y);end);
  return tbl;
end;

CHEVIE.compat.HeckeCharTableA := function(n, param,sqrtparam)local tbl,q,ci;
  q:=-param[1][1]/param[1][2];
  tbl:= CharTableSpecialized(CHEVIE.R("Hk","A"), [n+1, q]);
  tbl.cartan:= CartanMat("A", n);
  tbl.parameter:= List([1..n], x-> q);
  Inherit(tbl,CHEVIE.R("ClassInfo","A")(n));
  ci:=CHEVIE.R("CharInfo","A")(n);
  tbl.irredinfo:=Zip(ci.charparams,ci.charnames,function(x,y)return
    rec(charparam:=x,charname:=y);end);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end;
CHEVIE.compat.CharTableB := function(rank)local tbl,ci;
  tbl:=CharTable("WeylB",rank);
  tbl.identifier:=SPrint("W(B",rank,")");
  tbl.cartan:=CartanMat("B",rank);
  Inherit(tbl,CHEVIE.R("ClassInfo","B")(rank));
  ci:=CHEVIE.R("CharInfo","B")(rank);
  tbl.irredinfo:=Zip(ci.charparams,ci.charnames,function(x,y)return
    rec(charparam:=x,charname:=y);end);
  return tbl;
end;
CHEVIE.compat.HeckeCharTableB := function(n, param,sqrtparam) local q, tbl,ci;
  q:=List(param,x->-x[1]/x[2]);
  tbl:= CharTableSpecialized(CHEVIE.R("Hk","B"), [n, q[1], q[2]]);
  tbl.identifier := SPrint("H(B", n, ")");
  tbl.cartan:= CartanMat("B", n);
  tbl.parameter:= q;
  tbl.classtext:= CHEVIE.R("ClassInfo","B")(n).classtext;
  tbl.classparams:= List(tbl.classparam, x-> x[2]);
  tbl.classnames:= List(tbl.classparams, PartitionTupleToString);
  ci:=CHEVIE.R("CharInfo","B")(n);
  tbl.irredinfo:=Zip(ci.charparams,ci.charnames,function(x,y)return
    rec(charparam:=x,charname:=y);end);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end;
CHEVIE.compat.CharTableD := function(n)
  local hi, res, i, j, k, cls, chr, p, sy, syc, cli, pow,  S, Sval;

  hi:=CharTable("WeylB",n);
  if not IsBound(hi.size) then hi.size:=2^n*Factorial(n); fi;
  if n mod 2 = 0 then
    S:=CharTable("Symmetric",QuoInt(n,2));
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
  res.size:=QuoInt(hi.size,2);
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
