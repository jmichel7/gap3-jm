#############################################################################
##
#A  refltype.g  CHEVIE library    Jean Michel  2-2000
##
#Y  Copyright (C) 1999 Lehrstuhl D fur Mathematik, RWTH Aachen,
#Y  and University Paris VII.
##
##  This file contains operations for reflexion groups which are classified
##  and irreducible (for a Spets, irreducible means only one orbit under the
##  twist -- this includes descent of scalars).
##  It fetches the answers from the CHEVIE tables.
##  
##  A group is classified when its field .type has been filled.
##  (this attribute is obtained by the function 'ReflectionType').
##  This field is a list of records describing each irreducible component
##  (or an orbit of such components under the twist).
##  The argument of the ReflTypeOps routines is one such record.

ReflTypeOps:=OperationsRecord("ReflTypeOps");

ReflTypeOps.DegreeDescent:=function(t)
  if IsBound(t.orbit) then return Length(t.orbit); else return 1; fi;
end;

ReflTypeOps.CartanMat:=function(t)
  if IsBound(t.cartan) then return t.cartan;
  elif IsBound(t.cartanType) then 
       return CHEVIE.Data("CartanMat",t,t.cartanType);
  else return CHEVIE.Data("CartanMat",t);
  fi;
end;

ReflTypeOps.ReflectionDegrees:=function(t)local f,d,a;
  if not IsBound(t.orbit) then return CHEVIE.Data("ReflectionDegrees",t);fi;
  d:=CHEVIE.Data("ReflectionDegrees",t.orbit[1]);
# Let  t.scalar=[s_1,..,s_r],  where  r=Length(t.orbit)  and  let  p be the
# PhiFactor   of  t.twist  associated  to  the  reflection  degree  d_i  of
# t.orbit[1].   If   G0   is   the   Spets  described  by  t.orbit[1],  and
# G1:=Ennola(Product(t.scalar),G0)  then G is isomorphic  to the descent of
# scalars  of G1. According to spets 1.5, a Phifactor of Ennola(zeta,G0) is
# \zeta^{d_i}  times  that  of  G0;  and  by  spets  1.5  or [Digne-Michel,
# parabolic A.1] those of an a-descent of scalars are
# \zeta_a^j\zeta_i^{1/a} (all the a-th roots of \zeta_i).
  if OrderPerm(t.twist)>1 then f:=CHEVIE.Data("PhiFactors",t);
    if f=false then return false;fi;
  else f:=d*0+1;fi;
  f:=TransposedMat([d,f]);
  if IsBound(t.scalar) then 
    d:=Product(t.scalar); f:=List(f,p->[p[1],p[2]*d^p[1]]);
  fi;
  a:=Length(t.orbit);
  return Concatenation(TransposedMat(List(f,function(p)local r;
    r:=GetRoot(p[2],a);return List([0..a-1],i->[p[1],r*E(a)^i]);end)));
end;

ReflTypeOps.ReflectionCoDegrees:=function(t)local f,d,a;
  if not IsBound(t.orbit) then d:=CHEVIE.Data("ReflectionCoDegrees",t);
    if d<>false then return d;
    else d:=CHEVIE.Data("ReflectionDegrees",t);
      return Reversed(Maximum(d)-d); # assume well-generated
    fi;
  fi;
  d:=CHEVIE.Data("ReflectionCoDegrees",t);
  if d=false then  # assume well-generated
    d:=CHEVIE.Data("ReflectionDegrees",t.orbit[1]);
    a:=Position(d,Maximum(d));d:=Reversed(d[a]-d);
    if t.twist<>() then 
	# since factor is well-generated.
   # then by [BonLehMi 6.3 we have \eps^*_i=\eps_r/\eps_i where
   # \eps_r is attached to the Coxeter degree
       f:=CHEVIE.Data("PhiFactors",t);if f=false then return false;fi;
       f:=Reversed(List(f,x->f[a]/x));
    else f:=d*0+1;fi;
    d:=TransposedMat([d,f]);
  elif t.twist=() then d:=TransposedMat([d,d*0+1]);
  fi;
  if IsBound(t.scalar) then 
    f:=Product(t.scalar); d:=List(d,p->[p[1],p[2]*f^p[1]]);
  fi;
  a:=Length(t.orbit);
  return Concatenation(TransposedMat(List(d,function(p)local r;
    r:=GetRoot(p[2],a);return List([0..a-1],i->[p[1],r*E(a)^i]);end)));
end;

ReflTypeOps.GeneratingRoots:=t->ShallowCopy(CHEVIE.Data("GeneratingRoots",t));
ReflTypeOps.GeneratingCoRoots:=t->ShallowCopy(CHEVIE.Data("GeneratingCoRoots",t));

ReflTypeOps.WeightInfo:=function(t)local s;
  s:=["WeightInfo",t];
  if IsBound(t.cartanType) then Add(s,t.cartanType);fi;
  s:=ApplyFunc(CHEVIE.Data,s);
  if s=false then return rec(minusculeWeights:=[],minusculeCoweights:=[],
    decompositions:=[],moduli:=[]);
  fi;
  if not IsBound(s.minusculeCoweights) then 
    s.minusculeCoweights:=s.minusculeWeights;
  fi;
  return s;
end;

ReflTypeOps.ParabolicRepresentatives:=function(t,s)
   return CHEVIE.Data("ParabolicRepresentatives",t,s);
end;

ReflTypeOps.Invariants:=function(t)
  if IsBound(t.cartanType) then return CHEVIE.Data("Invariants",t,t.cartanType);
  else return CHEVIE.Data("Invariants",t);
  fi;
end;

ReflTypeOps.Discriminant:=function(t)
# if IsBound(t.cartanType) then return CHEVIE.Data("Discriminant",t,t.cartanType);
# else 
    return CHEVIE.Data("Discriminant",t);
# fi;
end;

ReflTypeOps.PoincarePolynomial:=function(t,param)
  return CHEVIE.Data("PoincarePolynomial",t,param);
end;

ReflTypeOps.CharName:=function(arg)
  return ApplyFunc(CHEVIE.Data,Concatenation(["CharName"],arg));
end;

ReflTypeOps.FakeDegree:=function(t,p,q)
  if IsBound(t.scalar) then 
       return CHEVIE.Data("FakeDegree",t,p,Product(t.scalar,s->q*s^-1));
  else return CHEVIE.Data("FakeDegree",t,p,q^ReflTypeOps.DegreeDescent(t));
  fi;
end;

ReflTypeOps.EigenvaluesGeneratingReflections:=function(t)
  if t.series<>"ST" then return List(CartanMat(t),x->1/2);
  else return CHEVIE.Data("EigenvaluesGeneratingReflections",t);
  fi;
end;

ReflTypeOps.ChevieClassInfo:=function(t)local r;
  r:=ShallowCopy(CHEVIE.Data("ClassInfo",t));
  if IsBound(t.orbit) and Length(t.orbit)>1 then
    if IsBound(r.classes) then r.classes:=r.classes*
       Product(ReflectionDegrees(t.orbit[1]))^(Length(t.orbit)-1);fi;
    # for nontrivial F the orders are  not easy to determine (and cannot
    # be determined  from the  ReflectionType for  nonsemisimple Spets),
    Unbind(r.orders);
  fi;
  if IsBound(t.orbit) then t:=t.orbit[1];fi;
  r.classtext:=List(r.classtext,c->t.indices{c});
  return r;
end;

ReflTypeOps.WGraph:=function(arg)
  return ApplyFunc(CHEVIE.Data,Concatenation(["WGraph"],arg));
end;

ReflTypeOps.HeckeCharTable:=function(arg)local tbl,t;
  tbl:=ApplyFunc(CHEVIE.Data,Concatenation(["HeckeCharTable"],arg));
  if tbl=false then return tbl;fi;
  if not IsBound(tbl.name) then tbl.name:=tbl.identifier; fi;
  if not IsBound(tbl.order) then tbl.order:=tbl.size; fi;
  t:=arg[1];
# descent of scalars: essentially the table is that of 
# the group corresponding to t.orbit[1] with t.twist acting on it.
  if IsBound(t.orbit) and Length(t.orbit)>1 then
     tbl.name:="";tbl.identifier:=tbl.name;
     tbl.size:=tbl.size^Length(t.orbit);
     tbl.classes:=List(tbl.centralizers,x->tbl.size/x);
  fi;
  return tbl;
end;

#############################################################################
##
#F  ReflTypeOps.CharTable( <t> ) . . . . . . . . . . . 'CharTable', using
#F  the classification
##
ReflTypeOps.CharTable:=function(t)local tbl;
  tbl:=CHEVIE.Data("CharTable",t);
  if tbl=false then return false;fi;
  if not IsBound(tbl.name) then tbl.name:=tbl.identifier; fi;
  if not IsBound(tbl.order) then tbl.order:=tbl.size; fi;
  if IsBound(t.orbit) then
    # for nontrivial F the orders are  not easy to determine (and cannot
    # be determined  from the  ReflectionType for  nonsemisimple Spets),
    # and powermaps do not make sense:
    if t.twist<>() then
      Unbind(tbl.orders);
      Unbind(tbl.powermap);
    fi;
    if Length(t.orbit)>1 then
      tbl.size:=tbl.size^Length(t.orbit);
      tbl.order:=tbl.size;
      tbl.classes:=List(tbl.centralizers,a->tbl.size/a);
    fi;
 #  if IsBound(t.scalar) then  # scalar-twisted Spetses
 #    # is the arbitrary choice below reasonable?
 #    tbl.irreducibles:=tbl.irreducibles*Product(t.scalar);
 #  fi;
  fi;
  return tbl;
end;

CHEVIE.GenericRepresentation:=function(opname)
  return function(arg)local t,i,res,l,F,i,r,f,tup,j;
  t:=arg[1];i:=arg[Length(arg)];
  res:=ApplyFunc(CHEVIE.Data,Concatenation([opname],arg));
  if res=false then Print(opname," not implemented for ",
    ReflectionName(t),"\n");return false;
  fi;
  if not IsBound(t.orbit) then return res;fi;
  if IsBound(t.orbit) and not IsRec(res) then
    res:=rec(gens:=res,F:=res[1]^0);fi;
  l:=Length(t.orbit);if l=1 then return res;fi;
  r:=Length(res.gens[1]);
  f:=function(x)if x=0 then return 1;else return 0;fi;end;
  res.gens:=Concatenation(List([1..l],i->List(res.gens,
    m->ApplyFunc(KroneckerProduct,List([1..l],k->m^f(k-i))))));
  F:=res.F; res.F:=NullMat(r^l);
  for i in [1..r^l] do tup:=CartesianAt([1..l]*0+r,i);
    for j in [1..r] do res.F[i][PositionCartesian([1..l]*0+r,
       Concatenation([j],tup{[1..Length(tup-1)]}))]:=F[tup[Length(tup)]][j];
    od;
  od;
  return res;
end;
end;

ReflTypeOps.Representation:=CHEVIE.GenericRepresentation("Representation");
ReflTypeOps.HeckeRepresentation:=CHEVIE.GenericRepresentation("HeckeRepresentation");

ReflTypeOps.NrConjugacyClasses:=function(arg)
  return ApplyFunc(CHEVIE.Data,Concatenation(["NrConjugacyClasses"],arg));
end;

ReflTypeOps.PrintDiagram:=function(t)local s,t;
  if IsBound(t.orbit) then
     if Length(t.orbit)>1 then
       Print("phi permutes the next ",Length(t.orbit)," components\n");
     fi;
     if t.twist<>() then
       Print("phi");
       if Length(t.orbit)>1 then Print("^",Length(t.orbit));fi;
       Print(" acts as ",t.twist," on the component below\n");
     fi;
     PrintDiagram(t.orbit);
  else 
    s:=["PrintDiagram",t,t.indices,ReflTypeOps.ReflectionName(t,rec())];
    if IsBound(t.cartanType) then Add(s,t.cartanType);fi;
    ApplyFunc(CHEVIE.Data,s);
  fi;
end;

ReflTypeOps.BraidRelations:=function(t)local m,p,r;
  if IsBound(t.orbit) then Error("not for Spets");
  elif t.series<>"ST" then
   m:=CoxeterMatrixFromCartanMat(CartanMat(t));
   p:=function(i,j,b)return List([1..b],k->i*(k mod 2)+j*((1-k)mod 2));end;
   r:=Concatenation(List([1..Length(m)],i->List([1..i-1],
     j->[p(i,j,m[i][j]),p(j,i,m[i][j])])));
  else r:=CHEVIE.Data("BraidRelations",t);
  fi;
  if IsBound(t.indices) then return List(r,x->List(x,y->t.indices{y}));
  else return r;
  fi;
end;

#############################################################################
##
#F  ReflectionName(<type>) . . . . . . . . . . . . name of a ReflectionType
##  
##  Gives a string which describes the isomorphism type of W (like "A2B3" for
##  the direct product of a root system of type A2 by one of type B3, or like
##  "I2(7)"). 
##

ReflTypeOps.ReflectionName:=function(arg)local t,option,res;
  t:=arg[1];
  if Length(arg)=1 then option:=rec();else option:=arg[2];fi;
  res:="";
  if IsBound(t.orbit) then  # Coxeter cosets
    if OrderPerm(t.twist)<>1 then 
      if IsBound(option.TeX) then Append(res,"{}^");fi;
      Append(res,String(OrderPerm(t.twist)));
      if t.twist=(1,4,2) then Append(res,"'");fi;
      if IsBound(option.TeX) then Append(res,"\\kern-0.2em ");fi;
    fi;
    if Length(t.orbit)=1 then 
         Append(res,ReflTypeOps.ReflectionName(t.orbit[1],option));
    else PrintToString(res,"(",ReflectionName(t.orbit,option),")");
    fi;
    if IsBound(t.scalar) and ForAny(t.scalar,x->x<>1) then
      PrintToString(res,"[",Join(List(t.scalar,x->Format(x,option))),"]");
    fi;
  else 
    if IsBound(t.tilde) then 
      if IsBound(option.TeX) then Append(res,"\\tilde ");
      else Append(res,"~");fi;
    fi;
    if IsBound(t.cartanType) then 
         Append(res,CHEVIE.Data("ReflectionName",t,option,t.cartanType));
    else Append(res,CHEVIE.Data("ReflectionName",t,option));
    fi;
  fi;
  return String(res);
end;

ReflTypeOps.String:=function(t)local l,ff;
  ff:=Filtered(RecFields(t),x->x in
    [ "bond", "cartanType", "indices", "orbit", "p", "q", 
      "rank", "series", "ST", "twist", "scalar", "tilde"]);
  l:=Maximum(List(ff,Length));
  return SPrint("rec(",Join(List(ff,function(f)local res;
     res:=SPrint(String(f,-l)," := ");
     if IsString(t.(f)) then Append(res,FormatGAP(t.(f)));
     else Append(res,String(t.(f)));fi;
     return res;end),",\n  "),")");
end;

ReflTypeOps.Print:=function(t)Print(String(t));end;

ReflTypeOps.ChevieCharInfo:=function(t)local res,f,get,uc;
  get:=function(f,F)local tmp;
    if not IsBound(res.(f)) then 
      tmp:=CHEVIE.Data(F,t);if tmp<>false then res.(f):=tmp;fi;
    fi;
  end;
  res:=CHEVIE.Data("CharInfo",t);
# if not IsBound(res.charnames) then
#   res.charnames:=List(res.charparams,
#          x->ApplyFunc(CHEVIE.Data,["CharName",t,x,rec()]));
# fi;
  res.positionId:=res.extRefl[1];
  res.positionDet:=res.extRefl[Length(res.extRefl)];
  get("b","LowestPowerFakeDegrees");
  get("B","HighestPowerFakeDegrees");
  get("a","LowestPowerGenericDegrees");
  get("A","HighestPowerGenericDegrees");
  if not IsBound(res.a) then
    uc:=CHEVIE.Data("UnipotentCharacters",t);
    if uc<>false then
      res.a:=uc.a{uc.harishChandra[1].charNumbers};
      res.A:=uc.A{uc.harishChandra[1].charNumbers};
    fi;
  fi;
  if IsBound(t.orbit) and not IsBound(res.charRestrictions) then
    res.charRestrictions:=[1..Length(res.charparams)];
    res.nrGroupClasses:=Length(res.charparams);
    # assume ortit twist trivial
  fi;
  for f in ["a","A","b","B"] do
  if IsBound(res.(f)) then res.(f):=res.(f)*ReflTypeOps.DegreeDescent(t);fi;
  od;
  return res;
end;

ReflTypeOps.KLeftCellRepresentatives:=function(t)
  return CHEVIE.Data("KLeftCellRepresentatives",t);
end;

ReflTypeOps.SemisimpleRank:=function(t)
  return CHEVIE.Data("SemisimpleRank",t);
end;

ReflTypeOps.DecompositionMatrix:=function(arg)local m,res,n;
  m:=ApplyFunc(CHEVIE.Data,Concatenation(["DecompositionMatrix"],arg));
  if m=false then return false;fi;
  n:=ApplyFunc(CHEVIE.Data,Concatenation(["NrConjugacyClasses"],arg{[1..Length(arg)-1]}));
  Append(m,List(Difference([1..n],Union(List(m,x->x[1]))),i->[[i],[[1]]]));
  res:=ApplyFunc(DiagonalMat,List(m,x->x[2]));
  SortParallel(Concatenation(List(m,x->x[1])),res);
  return res;
end;
