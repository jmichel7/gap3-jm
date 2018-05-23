#############################################################################
##
#A  dispatch.g               CHEVIE library                       Jean Michel
##
#Y  Copyright (C) 1996 - 2017  Lehrstuhl D fur Mathematik, RWTH Aachen,
#Y  University of St. Andrews, and  University Paris VII.
##
##  This file contains generic functions used in the package 'chevie'
##
##  Note: Inheritance diagram for Chevie reflection groups and cosets:
##
##          PermGroupOps
##              |          HasTypeOps
##          PermRootOps     /    . \   
##              |          /     . SpetsOps
##            ComplexGroupOps    .    |
##                     |         . CoxeterCosetOps
##                     |         .    
##                     |        AbsCoxOps  MatGroupOps
##                     |        /     \     |
##                 CoxeterGroupOps    GenCoxOps
##                                        |
##                                  AffineCoxeterGroupOps
##
##  AbsCox groups which happen to be finite inherit from HasType

AlphaInvolution:=Dispatcher("AlphaInvolution");
AltInvolution:=Dispatcher("AltInvolution");
BetaInvolution:=Dispatcher("BetaInvolution");
BraidRelations:=Dispatcher("BraidRelations");
BrieskornNormalForm:=Dispatcher("BrieskornNormalForm");

CharName:=function(arg)
  if Length(arg)=2 then Add(arg,rec());fi;
  return ApplyFunc(Dispatcher("CharName"),arg);end;

ClassName:=Dispatcher("ClassName");
Name:=Dispatcher("Name");
CoxeterElements:=Dispatcher("CoxeterElements");
CoxeterLength:=Dispatcher("CoxeterLength");
CoxeterWord:=Dispatcher("CoxeterWord");
CoxeterWords:=Dispatcher("CoxeterWords");
DualBraidMonoid:=Dispatcher("DualBraidMonoid");
EltWord:=Dispatcher("EltWord");
FactorizedSchurElement:=Dispatcher("FactorizedSchurElement");
FirstLeftDescending:=Dispatcher("FirstLeftDescending");
Fourier:=Dispatcher("Fourier");
Frobenius:=Dispatcher("Frobenius");
GenericOrder:=Dispatcher("GenericOrder");
Hecke:=Dispatcher("Hecke");
HeckeCharValues:=Dispatcher("HeckeCharValues");
HeckeClassPolynomials:=Dispatcher("HeckeClassPolynomials");
HeckeSubAlgebra:=Dispatcher("HeckeSubAlgebra");
HighestPowerGenericDegrees:=Dispatcher("HighestPowerGenericDegrees");
Invariants:=Dispatcher("Invariants");
Inversions:=Dispatcher("Inversions");
IsIsolated:=Dispatcher("IsIsolated");
IsLeftDescending:=Dispatcher("IsLeftDescending");
LeftDescentSet:=Dispatcher("LeftDescentSet");
LowestPowerGenericDegrees:=Dispatcher("LowestPowerGenericDegrees");
MatXPerm:=Dispatcher("MatXPerm");
NrConjugacyClasses:=Dispatcher("NrConjugacyClasses");
PermMatX:=Dispatcher("PermMatX");
QuasiIsolatedRepresentatives:=Dispatcher("QuasiIsolatedRepresentatives");
ReducedInRightCoset:=Dispatcher("ReducedInRightCoset");
ReflectionCharValue:=Dispatcher("ReflectionCharValue");
ReflectionLength:=Dispatcher("ReflectionLength");
ReflectionSubgroup:=Dispatcher("ReflectionSubgroup");
Reflections:=Dispatcher("Reflections");
ReflectionEigenvalues:=Dispatcher("ReflectionEigenvalues");
RelativeGroup:=Dispatcher("RelativeGroup");
RelativeCoset:=Dispatcher("RelativeCoset");
Representations:=Dispatcher("Representations");
SchurElement:=Dispatcher("SchurElement");
StandardParabolic:=Dispatcher("StandardParabolic");

Variables:=function(arg)
  if IsList(arg[1]) then return Union(List(arg[1],Variables));fi;
  return ApplyFunc(Dispatcher("Variables"),arg);
end;

TorusOrder:=Dispatcher("TorusOrder");
WeightInfo:=Dispatcher("WeightInfo");
WGraph:=Dispatcher("WGraph");

CharNumbers:=AttributeDispatcher("CharNumbers");
ChevieCharInfo:=AttributeDispatcher("ChevieCharInfo","charInfo");
ChevieClassInfo:=AttributeDispatcher("ChevieClassInfo","classInfo");
FactorizedSchurElements:=AttributeDispatcher("FactorizedSchurElements");
FieldOfDefinition:=AttributeDispatcher("FieldOfDefinition","field");
KLeftCellRepresentatives:=
  AttributeDispatcher("KLeftCellRepresentatives","cells0");
ReflectionDegrees:=AttributeDispatcher("ReflectionDegrees","degrees");
ReflectionCoDegrees:=AttributeDispatcher("ReflectionCoDegrees","codegrees");
SchurElements:=AttributeDispatcher("SchurElements");
SemisimpleRank:=AttributeDispatcher("SemisimpleRank");

FakeDegrees:=function(W,x)
  if IsRec(W) then 
    if not IsBound(W.fakeDegrees)  then
      if IsBound(W.operations) and IsBound(W.operations.FakeDegrees)
      then W.fakeDegrees:=W.operations.FakeDegrees(W,X(Cyclotomics));
      else Error( W, " has no method for FakeDegrees");
      fi;
    fi;
    return Value(W.fakeDegrees,x);
  else Error(W," should be a record"); fi;
end;

ReflectionName:=function(arg)local res;
   if Length(arg)=1 then Add(arg,rec());fi;
   if IsList(arg[1]) then 
     res:=List(arg[1],x->ReflectionName(x,arg[2]));
     if IsBound(arg[2].TeX) then return Join(res,"\\times ");
     else return Join(res,"x");
     fi;
   fi;
   return ApplyFunc(Dispatcher("ReflectionName"),arg);
end;

PositionId:=function(W)local ci;
  if IsCharTable(W) then
    if not IsBound(W.positionId) then
      W.positionId:=Position(W.irreducibles,W.irreducibles[1]*0+1);
    fi;
    return W.positionId;
  elif IsGroup(W) or IsSpets(W) then ci:=ChevieCharInfo(W);
    if not IsBound(ci.positionId) then 
      ci.positionId:=PositionId(CharTable(W));
    fi;
    return ci.positionId;
  else Error(W, " has no method for PositionId");
  fi;
end;

LowestPowerFakeDegrees:=function(W)local ci;
  ci:=ChevieCharInfo(W);
  if not IsBound(ci.b) then ci.b:=PermRootOps.LowestPowerFakeDegrees(W);fi;
  return ci.b;
end;

HighestPowerFakeDegrees:=function(W)local ci;
  ci:=ChevieCharInfo(W);
  if not IsBound(ci.B) then ci.B:=PermRootOps.HighestPowerFakeDegrees(W);fi;
  return ci.B;
end;

CharParams:=function(W)local ci;
  if IsRec(W) and IsBound(W.operations) and IsBound(W.operations.ChevieCharInfo)
  then ci:=ChevieCharInfo(W);
    if not IsBound(ci.charparams) then ci.charparams:=CharParams(CharTable(W));
    fi;
    return ci.charparams;
  elif IsCharTable(W) then
    if IsBound(W.irredinfo) and IsBound(W.irredinfo[1].charparam) then
      W.charparam:=List(W.irredinfo,x->x.charparam);
    else W.charparam:=[1..Length(W.irreducibles)];
    fi;
    return W.charparam;
  elif IsGroup(W) then return CharParams(CharTable(W));
  else Error(W, " has no method for CharParams");
  fi;
end;

CharNames:=function(arg)local W,option,res;
  W:=arg[1];
  if Length(arg)=2 then option:=arg[2];
  else option:=rec();
  fi;
  if option=rec() and IsRec(W) and IsBound(W.charNames) then
    return W.charNames;
  fi;
  if IsRec(W) and IsBound(W.operations) and IsBound(W.operations.CharNames) 
  then res:=W.operations.CharNames(W,option);
    if option=rec() then W.charNames:=res;fi;
    return res;
  elif IsGroup(W) then return CharNames(CharTable(W));
  else Error(W," has no method for CharNames");
  fi;
end;

AlmostCharNames:=function(arg)local W,option,res;
  W:=arg[1];
  if Length(arg)=2 then option:=arg[2];
  else option:=rec();
  fi;
  if IsRec(W) and IsBound(W.operations) and IsBound(W.operations.AlmostCharNames) then
   if not IsBound(option.TeX) and IsBound(W.almostCharNames) then return W.almostCharNames;fi;
   res:=W.operations.AlmostCharNames(W,option);
   if not IsBound(option.TeX) then W.almostCharNames:=res;fi;
   return res;
  else Error(W," has no method for AlmostCharNames");
  fi;
end;

ParabolicRepresentatives:=function(arg)
  if Length(arg)=1 then return Concatenation(List([0..SemisimpleRank(arg[1])],
    s->ParabolicRepresentatives(arg[1],s)));
  fi;
  return Dispatcher("ParabolicRepresentatives")(arg[1],arg[2]);
end;

#############################################################################
#F  ReflectionType( <W> )  . type of a Reflection group or Spets
##  
##  returns  a  list  of records (series:=,indices:=,rank:=) 
##  which  describe  irreducible components of W.
##
ReflectionType:=function(C)
  # dispatcher function, if not called with Cartan matrix:
  if IsRec(C) then return AttributeDispatcher("ReflectionType","type")(C);
  else return FiniteCoxeterTypeFromCartanMat(C);
  fi;
end;
#############################################################################
##
#F  PrintDiagram(<type>) . . . . . . . .Prints diagram for a ReflectionType
#F  PrintDiagram(<W>) . Prints diagram of Coxeter Datum or Hecke algebra, etc...
##  
##  The first form allows  PrintDiagram(ReflectionType(<cartan matrix>));
##
PrintDiagram:=function(W)local t;
  if IsRec(W) then Dispatcher("PrintDiagram")(W);return;fi;
  for t in W do PrintDiagram(t);od;
end;

Specialization:=function(H1,H2,f)
  if Group(H1)<>Group(H2) then 
    Error("Hecke algebras should be algebras for the same reflection group");fi;
  if ForAny([1..Length(H1.parameter)],i->ForAny([1,2],j->
      f(H1.parameter[i][j])<>H2.parameter[i][j])) then
    Error("f should send parameters of H1 to those of H2");fi;
  return t->t.operations.Specialization(t,H2,f);
end;
