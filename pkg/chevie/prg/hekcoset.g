#############################################################################
##
#A  hekcoset.g                     CHEVIE library                 Jean Michel
##
#Y  Copyright (C) 1996  Lehrstuhl D fuer Mathematik, RWTH Aachen, and
#Y  Universite Paris VII.
##
##  This file contains  GAP functions for working with Hecke cosets and
##  their character tables.
##############################################################################
HeckeCosetOps:=OperationsRecord("HeckeCosetOps",HeckeAlgebraOps);

HeckeCosetOps.Hecke:=HF->HF.hecke;
HeckeCosetOps.Group:=HF->Group(HF.hecke);

HeckeCosetOps.CharTable:=function(HF)local res,tbl,WF,H,tmp;
  WF:=ReflectionCoset(HF);
  H:=Hecke(HF);
  tmp:=List(ReflectionType(WF),function(o)local param,rootparam;
    param:=H.parameter{o.orbit[1].indices};
    rootparam:=Sublist(H.rootParameter,o.orbit[1].indices);
    tbl:=ReflTypeOps.HeckeCharTable(o,param,rootparam);
    return tbl;
  end);
  res:=tmp[1];
  for tbl in tmp{[2..Length(tmp)]} do
    res:=CharTableDirectProduct(res,tbl);
  od;
  res.classes:=List(res.centralizers,x->res.size/x);
  res.irredinfo:=List(CharParams(WF),x->rec(charparam:=x,
                                        charname:=CharName(WF,x,rec(TeX:=true))));
  Inherit(res,ChevieClassInfo(WF));
  res.name:=Concatenation("H(",ReflectionName(WF),")");
  res.identifier:=res.name;
  return res;
end;

HeckeCosetOps.Basis:=function(H,basis)
  if not IsBound(H.operations.(basis)) then 
      Error("basis ",basis," unknown");fi;
  return function(arg)local res;
   res:=Hecke(H).operations.(basis).MakeBasisElt(Hecke(H),basis,arg);
   res.coset:=H;res.elm:=res.elm*H.spets.phi;
   return res;
   end;
end;

HeckeCosetOps.HeckeSubAlgebra:=function(HF,subW)local WF,W;
  WF:=ReflectionCoset(HF);
  if IsList(subW) then subW:=CoxeterSubCoset(WF,subW);fi;
  return Hecke(subW,HeckeSubAlgebra(Hecke(HF),Group(subW)));
end;
