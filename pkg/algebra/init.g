#############################################################################
##
#A  init.g       ALGEBRA package         C'edric Bonnaf'e
##
#Y  Copyright (C) 2005 - 2010  University of Franche-Comt'e
##
##  This is the init file of the ALGEBRA package.
## 
#############################################################################

if not IsBound(ALGEBRA) then 
  ALGEBRA:=rec(path:=LOADED_PACKAGES.algebra);
fi;
ALGEBRA.date:=[2010,10];
ALGEBRA.name:="algebra";
ALGEBRA.copyright:="(C) C'edric Bonnaf'e -- Finite dimensional algebras";

PrintPkgInit(ALGEBRA);

ReadAlgebra:= function(name)
  if not ReadPath(ALGEBRA.path, name, ".g", "ReadAlgebra") then
     Error("Algebra library file '", name, "' must exist and be readable");
  fi;
end;

InfoAlgebra:=Print;

AUTO(ReadAlgebra("algebra"), AdamsOperation, AlgebraElement, AlgebraEltOps,
AlgebraHomomorphismByLinearity, FDAlgebraOps, Augmentation, ByDigits,
CentralIdempotents, CentralizerAlgebra, CharacterDecomposition,
CyclotomicModP, Digits, GeneralizedSolomonAlgebra, GeneralizedSolomonAlgebraOps,
GrothendieckRing, GrothendieckRingOps, GroupAlgebra, GroupAlgebraCentre,
GroupAlgebraCentreOps, GroupAlgebraOps, IsCentralElement,
IsDominated, LeftIdeal, LeftIndecomposableProjectives,
LeftTraces, LoewyLength, PBlocks, PiComponent, PiPart,
PiPrimeSections, PiSections, PolynomialQuotientAlgebra, PRank,
ProjectionMatrix, QuaternionAlgebra, QuotientAlgebra, RadicalPower,
RestrictionHomomorphism, RightIdeal, RightTraces, SolomonAlgebra,
SolomonAlgebraOps, SolomonHomomorphism, SubAlgebra, TablePrint, TwoSidedIdeal,
VectorSpaceProjection, ZeroHeckeAlgebra, ZeroHeckeAlgebraOps);

LeftTrace:=Dispatcher("LeftTrace");
RightTrace:=Dispatcher("RightTrace");
CartanMatrix:=Dispatcher("CartanMatrix");
