#############################################################################
##
#A  init.g            CHEVIE library / init file of the CHEVIE package.
## Meinolf Geck, Frank Luebeck, Gunter Malle, Jean Michel and G\"otz Pfeiffer.
##
#Y  Copyright (C) 1992 - 2018 Lehrstuhl D f\"ur Mathematik, RWTH Aachen,
#Y  Universitat Kaiserslautern, Germany and University Paris VII,  France.
## 
if not IsBound(CHEVIE) then 
  CHEVIE:=rec(path:=LOADED_PACKAGES.chevie);
fi;

CHEVIE.name:="chevie";
CHEVIE.date:=[2019,03,11];
CHEVIE.homepage:="http://webusers.imj-prg.fr/~jean.michel/chevie";
CHEVIE.copyright:=
["If you use CHEVIE in your work please cite the authors as follows:",
"[Jean Michel] The development version of the CHEVIE package of GAP3",
"Journal of algebra 435 (2015) 308--336",
"[Meinolf Geck, Gerhard Hiss, Frank Luebeck, Gunter Malle, Goetz Pfeiffer]",
"CHEVIE -- a system for computing and processing generic character tables",
"Applicable Algebra in Engineering Comm. and Computing 7 (1996) 175--210"];

PrintPkgInit(CHEVIE);

ReadChv:=function(name)
  if not ReadPath(CHEVIE.path, name, ".g", "ReadChv") then
    Error("CHEVIE library file '", name, ".g' must exist and be readable");
  fi;
end;

#############################################################################
##     AUTO section
##    
AUTO(ReadChv("lib/complex"),
  Complex, ComplexConjugate, ComplexOps, Cyclotomic, IsComplex);

AUTO(ReadChv("lib/cycpol"),
  AsRootOfUnity, CycPol, CycPolOps,
  IsCycPol, LcmCycPol);

AUTO(ReadChv("lib/decimal"),
  DecimalOps, Exp, IsDecimal, Log, Pi, Rational, SetDecimalPrecision, evalf);

AUTO(ReadChv("lib/symbol"),
  Compositions,
  CycPolFakeDegreeSymbol, CycPolGenericDegreeSymbol, DefectSymbol,
  DifferencePartitions, FullSymbol, HighestPowerFakeDegreeSymbol,
  HighestPowerGenericDegreeSymbol, LessSymbols,
  LowestPowerFakeDegreeSymbol, LowestPowerGenericDegreeSymbol,
  PartBeta, PartitionTupleToString, RankSymbol, ShiftBeta, StringSymbol,
  SymbolPartitionTuple, Symbols, SymbolsDefect, Tableaux, XSP);

AUTO(ReadChv("lib/util"),
  AbelianGenerators, ApplyWord, BraidRelation, CartesianAt,
  CharRepresentationWords, CheckRelation, ChevieIndeterminate,
  CollectCoefficients, DecomposeTensor, DetPerm, Dictionary, DifferenceMultiSet,
  Drop, EvalWords, EvalPolRoot, FastValue, FOrbit, FOrbits, GetRoot, GetWord,
  InductionTable, InductionTableOps, Inherit, InverseListsMap,
  MinimalWordProperty, PointsAndRepresentativesOrbits, PositionCartesian,
  PositionDet, Replace, Rotations, SymmetricDifference, TwoTree, 
  Dtime, Stime, Elapsed,
  EvScheme);

AUTO(ReadChv("lib/po"),
  Chains, GcdPartitions, Hasse, Incidence, LcmPartitions, LinearExtension, 
  Partition, Poset, PosetOps, ReversedPoset, IsJoinLattice, IsMeetLattice);

AUTO(ReadChv("lib/matrix"),
  BigCellDecomposition, BlocksMat, CoFactors, DecomposedMat,
  DistHelpedRepresentativeOperation, EigenvaluesMat, ExteriorPower,
  IndependentLines, IsNormalizing, MatStab, OnMatrices, PermutedByCols,
  PermMatMat, ProportionalityCoefficient, RepresentativeDiagonalConjugation, 
  RepresentativeRowColPermutation, SchurFunctor, SymmetricPower, Transporter);

AUTO(ReadChv("lib/format"),
  BracketIfNeeded, Format, FormatCoefficient, FormatCyclotomic, FormatGAP,
  FormatLaTeX, FormatMonomial, FormatPolynomial, FormatQuadratic,
  FormatQuotient, FormatTable, FormatTeX, IntListToString, FormatMaple,
  TeXBracket, ListUnbnd, TeXStrip);

AUTO(ReadChv("lib/factschur"),
  FactorizedSchurElementsOps);

AUTO(ReadChv("lib/sperm"),
  CyclesSignedPerm, SignedPermListList, SignedPermutationMat, SignPermuted,
  SignedPerm);

AUTO(ReadChv("prg/abscox"),
  AbsCoxOps, Bruhat, BruhatSmaller, CartanMatFromCoxeterMatrix,
  CoxeterHeckeAlgebraOps, CoxeterMatrix, CoxeterMatrixFromCartanMat,
  ForEachCoxeterWord, ForEachElement, IsCoxeterGroup, LongestCoxeterElement,
  LongestCoxeterWord, PermCosetsSubgroup, ReducedCoxeterWord ,
  ReducedRightCosetRepresentatives, RightDescentSet,
  BruhatPoset, StandardParabolicClass, ReducedExpressions);

AUTO(ReadChv("prg/abshecke"), 
  AbsHeckeOps, CheckHeckeDefiningRelations, IsHeckeAlgebra); 

AUTO(ReadChv("prg/affine"), 
  Affine, AffineCoxeterGroupOps, AffineRootAction);

AUTO(ReadChv("prg/compatib"), 
  AlgebraicFundamentalGroup,
  BraidWords, CartanName, CartanType, CharHeckeRepresentation,
  CoxeterConjugacyClasses, CoxeterElementsLength, CoxeterWordReflections,
  DirectSumCartanMat, DirectSumMat, DoublePartitionToString,
  DoublePartitions, HeckeCharTable, HeckeCharTableDirectProduct,
  HeckeFusion, HeckePowermap, HeckeReducedChar, HeckeScalarProducts, 
  KLCoefficient, LeftCellRepresentation, 
  LongestWeylWord, ParametersCentralizers, PermBraid, PhiFactors,
  PermCoxeterWord, PermRepresentationRoots, PermWeylWord, PositionSgn,
  PrintDynkinDiagram, ReducedInCoxeterCoset, ReducedWeylWord, Rootsystem,
  SemisimpleCentralizer, SimpleLeftDivisors,
  SimpleReflexionMatrices, Weyl, WeylClassPolynomials, WeylConjugacyClasses,
  WeylCosetPermRepresentation, WeylElements, WeylLengthPerm, WeylMueMat,
  WeylReflections, WeylRightCosetRepresentatives, WeylWordPerm, WordBraid);

AUTO(ReadChv("prg/complexr"),
  ComplexGroupOps, ComplexReflectionGroup, HeckeAlgebraOps, CyclicHeckeOps,
  HeckeCentralMonomials);

AUTO(ReadChv("prg/conjbrai"),
  AtomicMaps, CentralizerGenerators, MinConjugating, CategoryOps,
  CategoryByAtoms, ConjugacyCategory, ConjugacySet, PositiveSimpleConjugation,
  PreferredPrefix, RepresentativeConjugation, RepresentativeSC, ShowMaps);

AUTO(ReadChv("prg/coset"),
  CoxeterCoset, CoxeterCosetOps, CoxeterSubCoset, IsCoxeterCoset,
  Torus,TwistedPower, TwistingElements, Twistings);

AUTO(ReadChv("prg/coxeter"),
  AddComponentsCoxeterGroup , BadPrimes, CartanMat, CoxeterGroup, 
  CoxeterGroupOps, DescribeInvolution, ElementWithInversions, 
  ExtendedReflectionGroup, ExtendedGroupOps, ParabolicSubgroups,
  FiniteCoxeterTypeFromCartanMat, HighestShortRoot, IsExtendedGroup,
  IsomorphismType, Inversions, ParseTypeFromArg,
  PermMatXCoxeterElement, PermMatY, RootsCartan, SimpleRootsSubsystem,
  CoxeterTypeFromArg, PermutationsSimpleReflections);

AUTO(ReadChv("prg/dispatch"),
  AlmostCharNames, AlphaInvolution, AltInvolution, 
  BetaInvolution, BraidRelations, BrieskornNormalForm,
  CharName, CharNames, CharNumbers, CharParams, ChevieCharInfo, ChevieClassInfo,
  ClassName, CoxeterElements, CoxeterLength, CoxeterWord, 
  CoxeterWords, DualBraidMonoid, EltWord, 
  FactorizedSchurElement, FactorizedSchurElements, 
  FakeDegrees, FieldOfDefinition,
  FirstLeftDescending, Fourier, Frobenius, GenericOrder, Hecke, HeckeCharValues,
  HeckeClassPolynomials, HeckeSubAlgebra, HighestPowerFakeDegrees,
  HighestPowerGenericDegrees, Invariants, IsIsolated, IsLeftDescending, 
  KLeftCellRepresentatives, LeftDescentSet,
  LowestPowerFakeDegrees, LowestPowerGenericDegrees, MatXPerm,Name,
  NrConjugacyClasses, QuasiIsolatedRepresentatives, PermMatX, PositionId, 
  PrintDiagram, ParabolicRepresentatives,
  ReducedInRightCoset, ReflectionCharValue, ReflectionCoDegrees,
  ReflectionDegrees, ReflectionEigenvalues, ReflectionLength, ReflectionName,
  ReflectionSubgroup, ReflectionType, Reflections, RelativeCoset,
  RelativeGroup, Representations, SchurElement, SchurElements, 
  SemisimpleRank, Specialization,
  StandardParabolic, TorusOrder, Variables, WeightInfo, WGraph);

AUTO(ReadChv("prg/eigenspaces"),
  GetRelativeAction, GetRelativeRoot, EigenspaceProjector,
  PositionRegularClass, RegularEigenvalues, RelativeDegrees, SplitLevis);

AUTO(ReadChv("prg/garside"), 
  AsFraction, AsWord, BipartiteDecomposition, Braid, BraidMonoid,
  CompleteGarsideRecord, DualBraid, PermRootOpsDualBraidMonoid,
  EltBraid, GarsideAlpha,
  GarsideEltOps, GarsideOmega, GarsideWords, GoodCoxeterWord, 
  LeftDivisorsSimple, LeftGcd, LeftLcm, Presentation, ReversedWord, 
  RightGcd, RightLcm, ShrinkGarsideGeneratingSet,TwistedPowerMonoid,
  VeryGoodCoxeterWord);

AUTO(ReadChv("prg/gt"),
  ClassTypes,RationalUnipotentClasses,ClosedSubsets);

AUTO(ReadChv("prg/gencox"),
  CoxeterGroupByCartanMatrix, CoxeterGroupByCoxeterMatrix,
  CoxeterGroupHyperoctaedralGroup, CoxeterGroupSymmetricGroup,
  GenCoxOps);

AUTO(ReadChv("prg/hastype"),
  FakeDegree, HasTypeOps, ReflectionGroup);

AUTO(ReadChv("prg/hecke"),
  HeckeAlgebraOps, HeckeReflectionRepresentation, 
  HeckeCharValuesGood, JInductionTable, PoincarePolynomial);

AUTO(ReadChv("prg/heckeelt"), 
  CreateHeckeBasis, HeckeEltOps,HeckeElt,IsHeckeElt); 

AUTO(ReadChv("prg/hekcoset"),
  HeckeCosetOps);

AUTO(ReadChv("prg/kl"),
  AsymptoticAlgebra, CriticalPair, DualWGraph, KLMueMat, 
  KazhdanLusztigCoefficient, KazhdanLusztigMue, KazhdanLusztigPolynomial,
  LeftCells, LeftCellOps, RootParameter, WGraphToRepresentation,
  LeftCell, LeftStar, RightStar, LeftStars,
  PrepareForPolynomials, PrepareForMvp, LeftStarNC, QXHalf, RootParameter,
  WGraph2Representation);

AUTO(ReadChv("prg/permroot"),
  AsReflection, CartanCoefficient, HyperplaneOrbits, IsParabolic, IsWordFor, 
  ParabolicClosure, PermRootGroup, PermRootGroupNC, PermRootOps, MatYPerm, 
  Reflection, ReflectionCharacter, ReflectionWord, jInductionTable);

AUTO(ReadChv("prg/refltype"), 
  ReflTypeOps);

AUTO(ReadChv("prg/semisimple"),
  AlgebraicCentre, FixedPoints, IntermediateGroup,
  IsSemisimpleElement, IsQuasiIsolated, RootDatum, 
  SemisimpleCentralizerRepresentatives,
  SemisimpleElement, SemisimpleElementOps, SemisimpleSubgroup,
  StructureRationalPointsConnectedCentre, SubTorus);

AUTO(ReadChv("prg/spets"),
  GenericSign, IsSpets, PhiOnDiscriminant, 
  ReflectionCoset, Spets, SpetsOps, SubSpets);

AUTO(ReadChv("prg/sscoset"),
  RelativeDatum);

AUTO(ReadChv("prg/wclsinv"),
  ComponentWordsPerm, ComponentWordsPermCoset, CoxeterClassParamCheckFunction,
  CoxeterCosetClassParamCheckFunction, CoxeterCosetOpsClassInvariants,
  CoxeterGroupOpsClassInvariants);

AUTO(ReadChv("tbl/cmplximp"),
  ImprimitiveCuspidalName);

AUTO(ReadChv("tbl/exceptio"),VcycSchurElement, VFactorSchurElement);

if not IsBound(InfoChevie) then InfoChevie:=Ignore;fi;
if not IsBound(InfoChevie2) then InfoChevie2:=Ignore;fi;
if not IsBound(ChevieErr) then ChevieErr:=Print;fi; # for non-fatal errors

RequirePackage("vkcurve"); # Mvp requires VKcurve
RequirePackage("algebra"); # AlgebraElement requires Algebra
RequirePackage("specht");  # MatrixDecompositionMatrix requires Specht

ReadChv("auto"); # Reads record CHEVIE.AUTO

############################################################################
##
#F  CHEVIE.R( <f>, <field>)   . . . read function or data
#F   CHEVIE.(field).(f) from file determined by CHEVIE.AUTO
##  
CHEVIE.R:=function(f,field)
  if IsBound(CHEVIE.(field)) and IsBound(CHEVIE.(field).(f))
  then return CHEVIE.(field).(f);
  fi;
  if not IsBound(CHEVIE.AUTO.(field)) then
    Error(field," not in CHEVIE tables");
  elif not IsBound(CHEVIE.AUTO.(field).(f)) then
    InfoChevie2("#I no stored data for ",field,".",f,"\n");
    return false;
  fi;
  ReadChv(CHEVIE.AUTO.(field).(f));
  if not IsBound(CHEVIE.(field)) or not IsBound(CHEVIE.(field).(f))
  then Error("error in chevie/auto.g ",field,".",f,"\n");
  fi;
  return CHEVIE.(field).(f);
end;

############################################################################
##
#F  CHEVIE.Field(<type>)   . . . returns the field <f> such that information
##  about <type> is stored in CHEVIE.(f), plus extra arguments needed
##  to specify <type>
##  example: 
##  gap> CHEVIE.Field(rec(series:="B",rank:=2));
##  ["B",2]
##
CHEVIE.Field:=function(t)local field,orderphi,phi;
  if IsBound(t.orbit) then 
    phi:=t.twist;orderphi:=OrderPerm(t.twist);t:=t.orbit[1];
  else orderphi:=1;
  fi;
  field:=t.series;
  if field="I" then 
    if orderphi=2 then field:=SPrint(orderphi,field);fi;
    return [field,t.bond];
  elif field="ST" then 
    if IsBound(t.ST) then 
      if orderphi<>1 then return [SPrint(orderphi,"G",t.ST)]; # 2G5 and 2G7
      elif t.ST in [4..22] then return ["G4_22",t.ST];
      else return [SPrint("G",t.ST)];
      fi;
    elif orderphi<>1 then return ["timp",t.p,t.q,t.rank,phi];
    else return ["imp",t.p,t.q,t.rank];
    fi;
  elif orderphi=1 then 
    if field in ["A","B","D"] then return [field,t.rank];
    elif field in ["H","E","F","G"] then return [SPrint(field,t.rank)];
    fi;
  else field:=SPrint(orderphi,field);
    if field in ["2A","2D"] then return [field,t.rank];
    elif field="2B" then return ["2I",4];
    elif field="2G" then return ["2I",6];
    else return [SPrint(field,t.rank)]; # 2F4, 2E6, 3D4
    fi;
  fi;
end;

############################################################################
##
#F  CHEVIE.Data( <fname>, <type>,<extraargs>)   . . . general function
#F   to access information for irreducible types from libraries
##  
##  <fname> is the name of the function whose value is required, e.g. 
##  "CharParams" or "HeckeCharTable"
##
##  <type> is an element of the result of ReflectionType for
##  a CoxeterGroup, ComplexReflectionGroup, CoxeterCoset or Spets
##  
##  example: 
##  gap> CHEVIE.Data("CharParams",rec(series:="B",rank:=2));
##  [[[1,1],[ ]], [[1],[1]], [[],[1,1]], [[2],[]],[[],[2]]]
##
CHEVIE.Data:=function(arg)local f,t;
  t:=CHEVIE.Field(arg[2]);
  f:=CHEVIE.R(arg[1],t[1]);
  if IsFunc(f) then 
    return ApplyFunc(f,Concatenation(t{[2..Length(t)]},arg{[3..Length(arg)]}));
  fi;
  return f; # a variable
end;

CHEVIE.AddData:=function(name,types,data)local type;
  if not IsString(types[1]) then types:=[types];fi;
  for type in types do
    if not IsBound(CHEVIE.(type)) then
      InfoChevie2("#I cheviedata: creating type ",type,"\n");
      CHEVIE.(type):=rec();
    fi;
    InfoChevie2("#I adding CHEVIE.",type,".",name,"\n");
    CHEVIE.(type).(name):=data;
  od;
end;

CHEVIE.IndirectAddData:=function(f,types,fun)local t;
  for t in types do CHEVIE.AddData(f,t,fun(t));od;
end;

############################################################################
#F  CHEVIE.GetCached(parent,name,obj,hash) .. Function to maintain a cache
##  
##  GetCached   tests   if   an   object   x   in  parent.(name)  satisfies
##    hash(x)=hash(obj) and if so  returns x else adds obj to parent.(name)
##    and returns obj.
##  The cache is activated only if CHEVIE.Cache.(name)=true
##
CHEVIE.Cache:=rec(ReflectionSubgroups:=true);

CHEVIE.GetCached:=function(parent,name,obj,hash)local i;
  if not IsBound(CHEVIE.Cache.(name)) or not CHEVIE.Cache.(name) 
  then return obj;
  fi;
  if not IsBound(parent.(name)) then parent.(name):=[];fi;
  i:=hash(obj);i:=PositionProperty(parent.(name),x->hash(x)=i);
  if i=false then Add(parent.(name),obj);return obj;
  else return parent.(name)[i];
  fi;
end;

ReadChv("unip/init");
ReadChv("contr/init");
ReadChv("tbl/compat3");
