#############################################################################
##
#W  flags.g	 	   	Matrix Packages                  Frank Celler
#W                                                           & Eamonn O'Brien
#W                                                              & Anthony Pye
##
#H  @(#)$Id: flags.g,v 1.1 1997/03/10 13:52:02 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
RevisionMatrix.flags_g :=
    "@(#)$Id: flags.g,v 1.1 1997/03/10 13:52:02 gap Exp $";


#############################################################################
##

#F  AbsolutelyReducibleFlag ( <module> )
##
AbsolutelyReducibleFlag := function( module )
    if IsBound(module.absolutelyReducible) = false then return "unknown"; fi;
    return module.absolutelyReducible;
end;


#############################################################################
##
#F  AbstractGeneratorsFlag( <next> )
##
AbstractGeneratorsFlag := function( next )
    if IsBound(next.abstractGenerators)=false then return "unknown"; fi;
    return next.abstractGenerators;
end;


#############################################################################
##
#F  AlgElCharPolFacFlag ( <module> )
##
AlgElCharPolFacFlag := function( module )
   if IsBound(module.algElCharPolFac) = false then return "unknown"; fi;
   return module.algElCharPolFac;
end;


#############################################################################
##
#F  AlgElCharPolFlag( <module> )
##
AlgElCharPolFlag := function( module )
    if IsBound(module.algElCharPol) = false then return "unknown"; fi;
    return module.algElCharPol;
end;


#############################################################################
##
#F  AlgElFlag( <module> )
##
AlgElFlag := function( module )
    if IsBound(module.algEl) = false then return "unknown"; fi;
    return module.algEl;
end;


#############################################################################
##
#F  AlgElMatFlag( <module> )
##
AlgElMatFlag := function( module )
    if IsBound(module.algElMat) = false then return "unknown"; fi;
    return module.algElMat;
end;


#############################################################################
##
#F  AlgElNullspaceDimensionFlag ( <module> )
##
AlgElNullspaceDimensionFlag := function( module )
   if IsBound(module.algElNullspaceDim) = false then return "unknown"; fi;
   return module.algElNullspaceDim;
end;


#############################################################################
##
#F  AlgElNullspaceVecFlag ( <module> )
##
AlgElNullspaceVecFlag := function( module )
    if IsBound(module.algElNullspaceVec) = false then return "unknown"; fi;
    return module.algElNullspaceVec;
end;


#############################################################################
##
#F  AssignLayersFlag(next,i1,i2,mat)
##
AssignLayersFlag := function(next,i1,i2,mat)
    local tmp;
    
    tmp := next.layers;
    tmp[i1][i2] := mat;
    next.layers := tmp;
    
end;

    
#############################################################################
##
#F  AssignLayersVecFlag( next, i1, i2, v )
##
AssignLayersVecFlag := function( next, i1, i2, v )
    local tmp;
    
    tmp := next.layersVec;
    tmp[i1][i2] := v;
    next.layersVec := tmp;
    
end;


#############################################################################
##
#F  BasisFlag( <next> )
##
BasisFlag := function( next )
    if IsBound(next.basis)=false then return "unknown"; fi;
    return next.basis;
end;


#############################################################################
##
#F  BasisSubmoduleFlag( <next> )
##
BasisSubmoduleFlag := function( next )
    if IsBound(next.basisSubmodule)=false then return "unknown"; fi;
    return next.basisSubmodule;
end;


#############################################################################
##
#F  BlockFlag( <BlockSystem> )
##
BlockFlag := function( BlockSystem )
    return BlockSystem.block;
end;


#############################################################################
##
#F  BlockSystemFlag ( <module> )
##
BlockSystemFlag := function( module )
   if IsBound(module.blockSystem) = false then return "unknown"; fi;
   return module.blockSystem;
end;


#############################################################################
##
#F  CentMatFlag ( <module> )
##
CentMatFlag := function( module )
    if IsBound(module.centMat) = false then return "unknown"; fi;
    return module.centMat;
end;


#############################################################################
##
#F  CentMatMinPolyFlag ( <module> )
##
CentMatMinPolyFlag := function( module )
    if IsBound(module.centMatMinPoly) = false then return "unknown"; fi;
    return module.centMatMinPoly;
end;


#############################################################################
##
#F  ClassicalTypeFlag( <classic> )
##
ClassicalTypeFlag := function( c )
    if not IsBound(c.type)  then
        return "unknown";
    else
        return c.type;
    fi;
end;


#############################################################################
##
#F  DegreeFieldExtFlag ( <module> )
##
DegreeFieldExtFlag := function( module)
    if IsBound(module.degreeFieldExt) = false then return "unknown"; fi;
    return module.degreeFieldExt;
end;


#############################################################################
##
#F  DimensionAboveFlag( <next> )
##
DimensionAboveFlag := function( next )
    if IsBound(next.dimAbove)=false then return "unknown"; fi;
    return next.dimAbove;
end;


#############################################################################
##
#F  DimensionFlag ( <module> )
##
DimensionFlag := function( module )
    if IsBound(module.dimension) = false then return "unknown"; fi;
    return module.dimension;
end;


#############################################################################
##
#F  DimensionFlag( <next> )
##
DimensionFlag := function( next )
    if IsBound(next.dimension)=false then return "unknown"; 
    fi;
    return next.dimension;
end;


#############################################################################
##
#F  DimensionQuotientFlag( <next> )
##
DimensionQuotientFlag := function( next )
    if IsBound(next.dimQuotient)=false then return "unknown"; fi;
    return next.dimQuotient;
end;


#############################################################################
##
#F  DualFormFlag( <classic> )
##
DualFormFlag := function( c )
    if IsBound(c.dualForm)  then
        return c.dualForm;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  ExtraSpecialFlag ( <module> )
##
ExtraSpecialFlag := function( module )
    if IsBound(module.extraSpecial) = false then return "unknown"; fi;
    return module.extraSpecial;
end;


#############################################################################
##
#F  ExtraSpecialGroupFlag ( <module> )
##
ExtraSpecialGroupFlag := function( module)
    if IsBound(module.extraSpecialGroup) = false then return "unknown"; fi;
    return module.extraSpecialGroup;
end;


#############################################################################
##
#F  ExtraSpecialPrimeFlag ( <module> )
##
ExtraSpecialPrimeFlag := function( module )
    if IsBound(module.extraSpecialPrime) = false then return "unknown"; fi;
    return module.extraSpecialPrime;
end;


#############################################################################
##
#F  FieldFlag ( <module> )
##
FieldFlag := function( module )
    if IsBound(module.field) = false then return "unknown"; fi;
    return module.field;
end;


#############################################################################
##
#F  FpGroupFlag( <next> )
##
FpGroupFlag := function( next )
    if IsBound(next.fpGroup)=false then return "unknown"; fi;
    return next.fpGroup;
end;


#############################################################################
##
#F  FpHomomorphismFlag( <next> )
##
FpHomomorphismFlag := function( next )
    if IsBound(next.fpHomomorphism)=false then return "unknown"; fi;
    return next.fpHomomorphism;
end;


#############################################################################
##
#F  FrobeniusAutomorphismsFlag ( <module> )
##
FrobeniusAutomorphismsFlag := function( module )
    if IsBound(module.frobeniusAutomorphisms) = false  then
        return "unknown";
    fi;
    return module.frobeniusAutomorphisms;
end;


#############################################################################
##
#F  GeneratorsFlag( <module> )
##
GeneratorsFlag := function( module )
    if IsBound(module.generators) = false then return "unknown"; fi;
    return module.generators;
end;


#############################################################################
##
#F  IdentityBlockFlag( <next> )
##
IdentityBlockFlag := function( next )
    if IsBound(next.identityBlock)=false then return "unknown"; fi;
    return next.identityBlock;
end;


#############################################################################
##
#F  IdentityFlag( <next> )
##
IdentityFlag := function( next )
    if IsBound(next.identity)=false then return "unknown"; fi;
    return next.identity;
end;


#############################################################################
##
#F  IdentityQuotientFlag( <next> )
##
IdentityQuotientFlag := function( next )
    if IsBound(next.identityQuotient)=false then return "unknown"; fi;
    return next.identityQuotient;
end;


#############################################################################
##
#F  ImprimitiveFlag ( <module> )
##
ImprimitiveFlag := function( module )
   if IsBound(module.imprimitive) = false then return "unknown"; fi;
   return module.imprimitive;
end;


#############################################################################
##
#F  InvariantFormFlag( <classic> )
##
InvariantFormFlag := function( c )
    if IsBound(c.invariantForm)  then
        return c.invariantForm;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  InvariantFormFlag( <classic> )
##
InvariantFormFlag := function( c )
    if not IsBound(c.invariantForm)  then
        return "unknown";
    else
        return c.invariantForm;
    fi;
end;


#############################################################################
##
#F  IsFaithfulFlag( <next> )
##
IsFaithfulFlag := function( next )
    if IsBound(next.isFaithful)=false then return "unknown"; fi;
    return next.isFaithful;
end;

#############################################################################
##
#F  IsGenericFlag ( <next> )
##
IsGenericFlag  := function( c )
    if IsBound(c.isGeneric)  then
        return c.isGeneric;
    else
        return "unknown";
    fi;
end;

#############################################################################
##
#F  IsOrthogonalGroupFlag( <classic> )
##
IsOrthogonalGroupFlag := function( c )
    if IsBound(c.containsSO)  then
        return c.containsSO;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsPossibleImprimitiveFlag( <classic> )
##
IsPossibleImprimitiveFlag := function( c )
    if IsBound(c.isPossibleImprimitive)  then
        return c.isPossibleImprimitive;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsPossibleNormalizerPGroupFlag( <classic> )
##
IsPossibleNormalizerPGroupFlag := function( c )
    if IsBound(c.isPossibleMysteriousP)  then
        return c.isPossibleMysteriousP;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsPossibleSemiLinearFlag( <classic> )
##
IsPossibleSemiLinearFlag := function( c )
    if IsBound(c.isPossibleLargerField)  then
        return c.isPossibleLargerField;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsPossibleSmallerFieldFlag( <classic> )
##
IsPossibleSmallerFieldFlag := function( c )
    if IsBound(c.isPossibleSmallerField)  then
        return c.isPossibleSmallerField;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsPossibleTensorPowerFlag( <classic> )
##
IsPossibleTensorPowerFlag := function( c )
    if IsBound(c.isPossibleTensorPower)  then
        return c.isPossibleTensorPower;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsPossibleTensorProductFlag( <classic> )
##
IsPossibleTensorProductFlag := function( c )
    if IsBound(c.isPossibleTensorProduct)  then
        return c.isPossibleTensorProduct;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsSLContainedFlag( <classic> )
##
IsSLContainedFlag := function( c )
    if IsBound(c.containsSL)  then
        return c.containsSL;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsSymplecticGroupFlag( <classic> )
##
IsSymplecticGroupFlag := function( c )
    if IsBound(c.containsSP)  then
        return c.containsSP;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  IsUnitaryGroupFlag( <classic> )
##
IsUnitaryGroupFlag := function( c )
    if IsBound(c.containsSU)  then
        return c.containsSU;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  KernelFlag( <next> )
##
KernelFlag := function( next )
    if IsBound(next.kernel)=false then return "unknown"; fi;
    return next.kernel;
end;


#############################################################################
##
#F  LayerDimensionsFlag( <next> )
##
LayerDimensionsFlag := function( next )
    if IsBound(next.layerDimensions) = false then return "unknown"; fi;
    return next.layerDimensions;
end;


#############################################################################
##
#F  LayerNumberFlag( <next> )
##
LayerNumberFlag := function( next )
    if IsBound(next.layerNumber)=false then return "unknown"; fi;
    return next.layerNumber;
end;


#############################################################################
##
#F  LayersFlag( <next> )
##
LayersFlag := function( next )
    if IsBound(next.layers)=false then return "unknown"; fi;
    return next.layers;
end;


#############################################################################
##
#F  LayersVecFlag( <next> )
##
LayersVecFlag := function( next )
    if IsBound(next.layersVec)=false then return "unknown"; fi;
    return next.layersVec;
end;


#############################################################################
##
#F  LinearPartFlag( <module> )
##
LinearPartFlag := function( module )
    if IsBound(module.linearPart) = false then return "unknown"; fi;
    return module.linearPart;
end;


#############################################################################
##
#F  MapsFlag( <BlockSystem> )
##
MapsFlag := function( BlockSystem )
    return BlockSystem.maps;
end;

#############################################################################
##
#F  MaximumStripFlag( <next> )
##
MaximumStripFlag := function( next )
    if IsBound(next.maximumStrip)= false then return "unknown"; fi;
    return next.maximumStrip;
end;


#############################################################################
##
#F  NumberBlocksFlag( <BlockSystem> )
##
NumberBlocksFlag := function( BlockSystem )
   return BlockSystem.numberBlocks;
end;


#############################################################################
##
#F  PGroupFlag( <next> )
##
PGroupFlag := function( next )
    if IsBound(next.pGroup)=false then return "unknown"; fi;
    return next.pGroup;
end;


#############################################################################
##
#F  PermDomainFlag( <next> )
##
PermDomainFlag := function( next )
    if IsBound(next.permDomain)=false then return "unknown"; fi;
    return next.permDomain;
end;


#############################################################################
##
#F  PermGroupFlag( <BlockSystem> )
##
PermGroupFlag := function( BlockSystem )
    return BlockSystem.permGroup;
end;

#############################################################################
##
#F  PermGroupPFlag( <next> )
##
PermGroupPFlag := function( next )
    if IsBound(next.permGroupP)=false then return "unknown"; fi;
    return next.permGroupP;
end;


#############################################################################
##
#F  PossibleAlmostSimpleFlag( <classic> )
##
PossibleAlmostSimpleFlag := function( c )
    if not IsBound(c.possibleAlmostSimple)  then
        return "unknown";
    else
        return c.possibleAlmostSimple;
    fi;
end;

#############################################################################
##
#F  PossibleAlternatingGroupsFlag( <classic> )
##
PossibleAlternatingGroupsFlag := function( c )
    if not IsBound(c.possibleAlternatingGroups)  then
        return "unknown";
    else
        return c.possibleAlternatingGroups;
    fi;
end;


#############################################################################
##
#F  PossibleChevalleyGroupsFlag( <classic> )
##
PossibleChevalleyGroupsFlag := function( c )
    if not IsBound(c.possibleChevalleyGroups)  then
        return "unknown";
    else
        return c.possibleChevalleyGroups;
    fi;
end;


#############################################################################
##
#F  PossibleImprimitiveDimensionsFlag( <next> )
##
PossibleImprimitiveDimensionsFlag := function( next )
    if IsBound (next.possibleImprimitiveDimensions) = false then
        return "unknown";
    fi;
    return next.possibleImprimitiveDimensions;
end;

#############################################################################
##
#F  PossibleNearlySimpleFlag( <classic> )
##
PossibleNearlySimpleFlag := function( c )
    if not IsBound(c.possibleNearlySimple)  then
        return "unknown";
    else
        return c.possibleNearlySimple;
    fi;
end;

#############################################################################
##
#F  PossibleOverLargerFieldFlag ( <next> )
##
PossibleOverLargerFieldFlag  := function( c )
    if not IsBound(c.possibleOverLargerField)  then
        return "unknown";
    else
        return c.possibleOverLargerField;
    fi;
end;


#############################################################################
##
#F  PossibleSmallerFieldFlag( <classic> )
##
PossibleSmallerFieldFlag := function( c )
    if IsBound(c.possibleSmallerField)  then
        return c.possibleSmallerField;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  PossibleSporadicGroupsFlag( <classic> )
##
PossibleSporadicGroupsFlag := function( c )
    if not IsBound(c.possibleSporadicGroups)  then
        return "unknown";
    else
        return c.possibleSporadicGroups;
    fi;
end;

#############################################################################
##
#F  PossibleTensorDimensionsFlag( <next> )
##
PossibleTensorDimensionsFlag := function( next )
   if IsBound (next.possibleTensorDimensions) = false then
       return "unknown";
   fi;
   return next.possibleTensorDimensions;
end;


#############################################################################
##
#F  PrimeFlag( <next> )
##
PrimeFlag := function( next )
    if IsBound(next.prime)=false then return "unknown"; fi;
    return next.prime;
end;


#############################################################################
##
#F  PrintLevelFlag( <next> )
##
PrintLevelFlag := function( next )
    if IsBound(next.printLevel)=false then return "unknown"; fi;
    return next.printLevel;
end;


#############################################################################
##
#F  QuadraticFormFlag( <classic> )
##
QuadraticFormFlag := function( c )
    if not IsBound(c.quadraticForm)  then
        return "unknown";
    else
        return c.quadraticForm;
    fi;
end;


#############################################################################
##
#F  QuotientFlag( <next> )
##
QuotientFlag := function( next )
    if IsBound(next.quotient)=false then return "unknown"; fi;
    return next.quotient;
end;


#############################################################################
##
#F  RecogniseFlag( <next> )
##
RecogniseFlag := function( next )
    if IsBound(next.recognise) = false then return "unknown"; fi;
    return next.recognise;
end;


#############################################################################
##
#F  ReducibleFlag ( <module> )
##
ReducibleFlag := function( module )
    if IsBound(module.reducible) = false then return "unknown"; fi;
    return module.reducible;
end;


#############################################################################
##
#F  SemiLinearFlag ( <module> )
##
SemiLinearFlag := function( module )
    if IsBound(module.semiLinear) = false then return "unknown"; fi;
    return module.semiLinear;
end;


#############################################################################
##
#F  SizeExtensionFlag( <classic> )
##
SizeExtensionFlag := function( c )
    if not IsBound(c.sizeExtension)  then
        return "unknown";
    else
        return c.sizeExtension;
    fi;
end;


#############################################################################
##
#F  SizeFlag( <classic> )
##
SizeFlag := function( c )
    if not IsBound(c.size)  then
        return "unknown";
    else
        return c.size;
    fi;
end;


#############################################################################
##
#F  SizeQuotientFlag( <next> )
##
SizeQuotientFlag := function( next )
    if IsBound(next.sizeQuotient)=false then return "unknown"; fi;
    return next.sizeQuotient;
end;


#############################################################################
##
#F  SubbasisFlag ( <module> )
##
SubbasisFlag := function( module )
    if IsBound(module.subbasis) = false then return "unknown"; fi;
    return module.subbasis;
end;


#############################################################################
##
#F  SuccessiveStripFlag( <next> )
##
SuccessiveStripFlag := function( next )
    if IsBound(next.successiveStrip)=false then return "unknown"; fi;
    return next.successiveStrip;
end;


#############################################################################
##
#F  SymTensorBasisFlag ( <module> )
##
SymTensorBasisFlag := function( module )
    if IsBound(module.symTensorBasis) = false then return "unknown"; fi;
    return module.symTensorBasis;
end;


#############################################################################
##
#F  SymTensorFactorsFlag ( <module> )
##
SymTensorFactorsFlag := function( module )
    if IsBound(module.symTensorFactors) = false then return "unknown"; fi;
    return module.symTensorFactors;
end;


#############################################################################
##
#F  SymTensorPermFlag ( <module> )
##
SymTensorPermFlag := function( module )
    if IsBound(module.symTensorPerm) = false then return "unknown"; fi;
    return module.symTensorPerm;
end;


#############################################################################
##
#F  SymTensorProductFlag ( <module> )
##
SymTensorProductFlag := function( module )
    if IsBound(module.symTensorProduct) = false then return "unknown"; fi;
    return module.symTensorProduct;
end;


#############################################################################
##
#F  TensorBasisFlag ( <module> )
##
TensorBasisFlag := function( module )
    if IsBound(module.tensorBasis) = false then return "unknown"; fi;
    return module.tensorBasis;
end;


#############################################################################
##
#F  TensorFactorsFlag ( <module> )
##
TensorFactorsFlag := function( module )
    if IsBound(module.tensorFactors) = false then return "unknown"; fi;
    return module.tensorFactors;
end;


#############################################################################
##
#F  TensorProductFlag ( <module> )
##
TensorProductFlag := function( module )
    if IsBound(module.tensorProduct) = false then return "unknown"; fi;
    return module.tensorProduct;
end;


#############################################################################
##
#F  TypeFlag( <next> )
##
TypeFlag := function( next )
    if IsBound(next.type)=false then return "unknown"; fi;
    return next.type;
end;


#############################################################################
##
#F  UndoAbsolutelyReducibleFlag( <module> )
##
UndoAbsolutelyReducibleFlag := function( module )
    Unbind (module.absolutelyReducible);
end;


#############################################################################
##
#F  UndoAbstractGeneratorsFlag( <next> )
##
UndoAbstractGeneratorsFlag := function( next )
    Unbind(next.abstractGenerators);
end;


#############################################################################
##
#F  UndoBasisFlag( <next> )
##
UndoBasisFlag := function( next )
    Unbind(next.basis);
end;


#############################################################################
##
#F  UndoBasisSubmoduleFlag( <next> )
##
UndoBasisSubmoduleFlag := function( next )
    Unbind(next.basisSubmodule);
end;


#############################################################################
##
#F  UndoCentMatFlag ( <module> )
##
UndoCentMatFlag := function( module )
    Unbind(module.centMat);
end;


#############################################################################
##
#F  UndoDegreeFieldExtFlag( <module> )
##
UndoDegreeFieldExtFlag := function( module)
    Unbind(module.degreeFieldExt);
end;


#############################################################################
##
#F  UndoDimensionAboveFlag( <next> )
##
UndoDimensionAboveFlag := function( next )
    Unbind(next.dimAbove);
end;


#############################################################################
##
#F  UndoDimensionQuotientFlag( <next> )
##
UndoDimensionQuotientFlag := function( next )
    Unbind(next.dimQuotient);
end;


#############################################################################
##
#F  UndoFpHomomorphismFlag( <next> )
##
UndoFpHomomorphismFlag := function( next )
    Unbind(next.fpHomomorphism);
end;


#############################################################################
##
#F  UndoGeneratorsFlag( <next> )
##
UndoGeneratorsFlag := function( next )
    Unbind(next.generators);
end;


#############################################################################
##
#F  UndoIdentityBlockFlag( <next> )
##
UndoIdentityBlockFlag := function( next )
    Unbind(next.identityBlock);
end;


#############################################################################
##
#F  UndoIdentityFlag( <next> )
##
UndoIdentityFlag := function( next )
    Unbind(next.identity);
end;


#############################################################################
##
#F  UndoIdentityQuotientFlag( <next> )
##
UndoIdentityQuotientFlag := function( next )
    Unbind(next.identityQuotient);
end;


#############################################################################
##
#F  UndoKernelFlag( <next> )
##
UndoKernelFlag := function( next )
    Unbind(next.kernel);
end;


#############################################################################
##
#F  UndoMaximumStripFlag( <next> )
##
UndoMaximumStripFlag := function( next )
    Unbind(next.maximumStrip);
end;


#############################################################################
##
#F  UndoPermDomainFlag( <next> )
##
UndoPermDomainFlag := function( next )
    Unbind(next.permDomain);
end;


#############################################################################
##
#F  UndoPermGroupPFlag( <next> )
##
UndoPermGroupPFlag := function( next )
    Unbind(next.permGroupP);
end;


#############################################################################
##
#F  UndoQuotientFlag( <next> )
##
UndoQuotientFlag := function( next )
    Unbind(next.quotient);
end;


#############################################################################
##
#F  UndoReducibleFlag ( <module> )
##
UndoReducibleFlag := function( module )
    Unbind (module.reducible);
end;


#############################################################################
##
#F  UndoSuccessiveStripFlag( <next> )
##
UndoSuccessiveStripFlag := function( next )
    Unbind(next.successiveStrip);
end;


#############################################################################
##
#F  UndoTensorBasisFlag( <module> )
##
UndoTensorBasisFlag := function( module )
   Unbind (module.TensorBasis);
end;


#############################################################################
##
#F  UndoTensorFactorsFlag ( <module> )
##
UndoTensorFactorsFlag := function( module )
    Unbind(module.tensorFactors);
end;


#############################################################################
##
#F  UndoTensorProductFlag ( <module> )
##
UndoTensorProductFlag := function( module )
    Unbind(module.tensorProduct);
end;


#############################################################################
##
#F  UnitaryFormFlag( <classic> )
##
UnitaryFormFlag := function( c )
    if IsBound(c.unitaryForm)  then
        return c.unitaryForm;
    else
        return "unknown";
    fi;
end;


#############################################################################
##

#F  SetAbsolutelyReducibleFlag( <module>, <absolutelyReducible> )
##
SetAbsolutelyReducibleFlag := function( module, absolutelyReducible )
    module.absolutelyReducible := absolutelyReducible;
end;


#############################################################################
##
#F  SetAbstractGeneratorsFlag( <next>, <gens> )
##
SetAbstractGeneratorsFlag := function( next, gens )
    next.abstractGenerators := gens;
end;


#############################################################################
##
#F  SetAlgElCharPolFacFlag( <module>, <p> )
##
SetAlgElCharPolFacFlag := function( module, p )
    module.algElCharPolFac := p;
end;


#############################################################################
##
#F  SetAlgElCharPolFlag( <module>, <p> )
##
SetAlgElCharPolFlag := function( module, p )
   module.algElCharPol := p;
end;


#############################################################################
##
#F  SetAlgElFlag( <module>, <algel> )
##
SetAlgElFlag := function( module, algel )
    module.algEl := algel;
end;


#############################################################################
##
#F  SetAlgElMatFlag( <module>, <M> )
##
SetAlgElMatFlag := function( module, M )
    module.algElMat := M;
end;


#############################################################################
##
#F  SetAlgElNullspaceDimensionFlag ( <module>, <ndim> )
##
SetAlgElNullspaceDimensionFlag := function( module, ndim )
    module.algElNullspaceDim := ndim;
end;


#############################################################################
##
#F  SetAlgElNullspaceVecFlag ( <module>, <v> )
##
SetAlgElNullspaceVecFlag := function( module, v )
    module.algElNullspaceVec := v;
end;


#############################################################################
##
#F  SetBasisFlag( <next>, <mat> )
##
SetBasisFlag := function( next, mat )
    next.basis := mat;
end;


#############################################################################
##
#F  SetBasisSubmoduleFlag( <next>, <mat> )
##
SetBasisSubmoduleFlag := function( next, mat )
    next.basisSubmodule := mat;
end;


#############################################################################
##
#F  SetBlockSystemFlag ( <module>, <SystemOfBlocks> )
##
SetBlockSystemFlag := function( module, SystemOfBlocks )
   module.blockSystem := SystemOfBlocks;
end;


#############################################################################
##
#F  SetCentMatFlag( <module>, <mat> )
##
SetCentMatFlag := function( module, mat )
    module.centMat := mat;
end;


#############################################################################
##
#F  SetCentMatMinPolyFlag ( <module>, <p> )
##
SetCentMatMinPolyFlag := function( module, p )
    module.centMatMinPoly := p;
end;


#############################################################################
##
#F  SetClassicalTypeFlag( <classic>, <type> )
##
SetClassicalTypeFlag := function( c, type )
    c.type := type;
end;


#############################################################################
##
#F  SetDegreeFieldExtFlag ( <module>, <deg> )
##
SetDegreeFieldExtFlag := function( module, deg )
    module.degreeFieldExt := deg;
end;


#############################################################################
##
#F  SetDimensionAboveFlag( <next>, <d> )
##
SetDimensionAboveFlag := function( next, d )
    next.dimAbove := d;
end;


#############################################################################
##
#F  SetDimensionFlag( <module>, <dimension> )
##
SetDimensionFlag := function( module, dimension )
    module.dimension := dimension;
end;


#############################################################################
##
#F  SetDimensionQuotientFlag( <next>, <d> )
##
SetDimensionQuotientFlag := function( next, d )
    next.dimQuotient := d;
end;


#############################################################################
##
#F  SetDualFormFlag( <classic>, <form> )
##
SetDualFormFlag := function( c, form )
    c.dualForm := form;
end;


#############################################################################
##
#F  SetExtraSpecialFlag ( <module>, <extraspecial> )
##
SetExtraSpecialFlag := function( module, extraspecial )
    module.extraSpecial := extraspecial;
end;


#############################################################################
##
#F  SetExtraSpecialGroupFlag ( <module>, <extraspecialgroup> )
##
SetExtraSpecialGroupFlag := function( module, extraspecialgroup )
    module.extraSpecialGroup := extraspecialgroup;
end;


#############################################################################
##
#F  SetExtraSpecialPrimeFlag ( <module>, <extraspecialprime> )
##
SetExtraSpecialPrimeFlag := function( module, extraspecialprime )
    module.extraSpecialPrime := extraspecialprime;
end;


#############################################################################
##
#F  SetFieldFlag( <module>, <F> )
##
SetFieldFlag := function( module, F )
    module.field := F;
end;


#############################################################################
##
#F  SetFpGroupFlag( <next>, <f> )
##
SetFpGroupFlag := function( next, f )
    next.fpGroup := f;
end;


#############################################################################
##
#F  SetFpHomomorphismFlag( <next>, <map> )
##
SetFpHomomorphismFlag := function( next, map )
    next.fpHomomorphism := map;
end;


#############################################################################
##
#F  SetFrobeniusAutomorphismsFlag( <module>, <frobeniusAutomorphisms> )
##
SetFrobeniusAutomorphismsFlag := function( module, frobeniusAutomorphisms )
    module.frobeniusAutomorphisms := frobeniusAutomorphisms;
end;


#############################################################################
##
#F  SetGeneratorsFlag( <module>, <generators> )
##
SetGeneratorsFlag := function( module, generators )
    module.generators := generators;
end;


#############################################################################
##
#F  SetIdentityBlockFlag( <next>, <mat> )
##
SetIdentityBlockFlag := function( next, mat )
    next.identityBlock := mat;
end;


#############################################################################
##
#F  SetIdentityFlag( <next>, <mat> )
##
SetIdentityFlag := function( next, mat )
    next.identity := mat;
end;


#############################################################################
##
#F  SetIdentityQuotientFlag( <next>, <mat> )
##
SetIdentityQuotientFlag := function( next, mat )
    next.identityQuotient := mat;
end;


#############################################################################
##
#F  SetImprimitiveFlag ( <module>, <imprimitive> )
##
SetImprimitiveFlag := function( module, imprimitive )
    module.imprimitive := imprimitive;
end;


#############################################################################
##
#F  SetInvariantFormFlag( <classic>, <form> )
##
InvariantFormFlag := function( c, form )
    c.invariantFormFlag := form;
end;


#############################################################################
##
#F  SetInvariantFormFlag( <classic>, <form> )
##
SetInvariantFormFlag := function( c, form )
    c.invariantForm := form;
end;


#############################################################################
##
#F  SetIsFaithfulFlag( <next>, <bool> )
##
SetIsFaithfulFlag := function( next, bool )
    next.isFaithful := bool;
end;


#############################################################################
##
#F  SetIsOrthogonalGroupFlag( <classic>, <flag> )
##
SetIsOrthogonalGroupFlag := function( c, flag )
    c.containsSO := flag;
end;


#############################################################################
##
#F  SetIsPossibleImprimitiveFlag( <classic>, <flag> )
##
SetIsPossibleImprimitiveFlag := function( c, flag )
    c.isPossibleImprimitive := flag;
end;


#############################################################################
##
#F  SetIsPossibleNormalizerPGroupFlag( <classic>, <field> )
##
SetIsPossibleNormalizerPGroupFlag := function( c, f )
    c.isPossibleMysteriousP := f;
end;


#############################################################################
##
#F  SetIsPossibleSemiLinearFlag( <classic>, <field> )
##
SetIsPossibleSemiLinearFlag := function( c, f )
    c.isPossibleLargerField := f;
end;


#############################################################################
##
#F  SetIsPossibleSmallerFieldFlag( <classic>, <field> )
##
SetIsPossibleSmallerFieldFlag := function( c, f )
    c.isPossibleSmallerField := f;
end;


#############################################################################
##
#F  SetIsPossibleTensorPowerFlag( <classic>, <flag> )
##
SetIsPossibleTensorPowerFlag := function( c, flag )
    c.isPossibleTensorPower := flag;
end;


#############################################################################
##
#F  SetIsPossibleTensorProductFlag( <classic>, <flag> )
##
SetIsPossibleTensorProductFlag := function( c, flag )
    c.isPossibleTensorProduct := flag;
end;


#############################################################################
##
#F  SetIsSLContainedFlag( <classic>, <flag> )
##
SetIsSLContainedFlag := function( c, flag )
    c.containsSL := flag;
end;


#############################################################################
##
#F  SetIsSymplecticGroupFlag( <classic>, <flag> )
##
SetIsSymplecticGroupFlag := function( c, flag )
    c.containsSP := flag;
end;


#############################################################################
##
#F  SetIsUnitaryGroupFlag( <classic>, <flag> )
##
SetIsUnitaryGroupFlag := function( c, flag )
    c.containsSU := flag;
end;


#############################################################################
##
#F  SetKernelFlag( <next>, <record> )
##
SetKernelFlag := function( next, record )
    next.kernel := record;
end;


#############################################################################
##
#F  SetLayerDimensionsFlag( <next>, <dims> )
##
SetLayerDimensionsFlag := function( next, dims )
    next.layerDimensions := dims;
end;


#############################################################################
##
#F  SetLayerNumberFlag( <next>, <l> )
##
SetLayerNumberFlag := function( next, l )
    next.layerNumber := l;
end;


#############################################################################
##
#F  SetLayersFlag( <next>, <layers> )
##
SetLayersFlag := function( next, layers )
    next.layers := layers;
end;


#############################################################################
##
#F  SetLayersVecFlag( <next>, <layersVec> )
##
SetLayersVecFlag := function( next, layersVec )
    next.layersVec := layersVec;
end;


#############################################################################
##
#F  SetLinearPartFlag ( <module>, <linearPart> )
##
SetLinearPartFlag := function( module, linearPart )
    module.linearPart := linearPart;
end;


#############################################################################
##
#F  SetMaximumStripFlag( <next>, <m> )
##
SetMaximumStripFlag := function( next, m )
    next.maximumStrip := m;
    if KernelFlag(next) = "unknown" then
        return;
    fi;
    SetMaximumStripFlag(KernelFlag(next),m);
end;


#############################################################################
##
#F  SetPGroupFlag( <next>, <igs> )
##
SetPGroupFlag := function( next, igs )
    next.pGroup := igs;
end;


#############################################################################
##
#F  SetPermDomainFlag( <next>, <permD> )
##
SetPermDomainFlag := function( next, permD )
    next.permDomain := permD;
end;


#############################################################################
##
#F  SetPermGroupPFlag( <next>, <permgp> )
##
SetPermGroupPFlag := function( next, permgp )
    next.permGroupP := permgp;
end;


#############################################################################
##
#F  SetPossibleAlmostSimpleFlag( <next>, <poss> )
##
SetPossibleAlmostSimpleFlag := function( next, poss )
    next.possibleAlmostSimple := poss;
end;


#############################################################################
##
#F  SetPossibleAlternatingGroupsFlag( <classic>, <degrees> )
##
SetPossibleAlternatingGroupsFlag := function( c, degrees )
    c.possibleAlternatingGroups := degrees;
end;


#############################################################################
##
#F  SetPossibleChevalleyGroupsFlag( <classic>, <groups> )
##
SetPossibleChevalleyGroupsFlag := function( c, groups )
    c.possibleChevalleyGroups := groups;
end;


#############################################################################
##
#F  SetPossibleImprimitiveDimensionsFlag( <next>, <possibleImpDims> )
##
SetPossibleImprimitiveDimensionsFlag := function( next, possibleImpDims )
    next.possibleImprimitiveDimensions := possibleImpDims;
end;


#############################################################################
##
#F  SetPossibleSmallerFieldFlag( <classic>, <field> )
##
SetPossibleSmallerFieldFlag := function( c, f )
    c.possibleSmallerField := f;
end;


#############################################################################
##
#F  SetPossibleSporadicGroupsFlag( <classic>, <names> )
##
SetPossibleSporadicGroupsFlag := function( c, names )
    c.possibleSporadicGroups := names;
end;


#############################################################################
##
#F  SetPossibleTensorDimensionsFlag( <next>, <possibleTensorDims> )
##
SetPossibleTensorDimensionsFlag := function( next, possibleTensorDims )
    next.possibleTensorDimensions := possibleTensorDims;
end;


#############################################################################
##
#F  SetPrimeFlag( <next>, <p> )
##
SetPrimeFlag := function( next, p )
    next.prime := p;
end;


#############################################################################
##
#F  SetPrintLevelFlag( <next>, <lev> )
##
SetPrintLevelFlag := function( next, lev )
    if lev < 1 or lev > 3 then 
        Error(" Print level out of range [1..3] ");
    fi;
    next.printLevel := lev;
    if KernelFlag(next) = "unknown" then
        return;
    fi;
    SetPrintLevelFlag(KernelFlag(next),lev);
end;


#############################################################################
##
#F  SetQuadraticFormFlag( <classic>, <form> )
##
SetQuadraticFormFlag := function( c, form )
    c.quadraticForm := form;
end;


#############################################################################
##
#F  SetQuotientFlag( <next>, <g> )
##
SetQuotientFlag := function( next, g )
    next.quotient := g;
end;


#############################################################################
##
#F  SetReducibleFlag ( <module>, <reducible> )
##
SetReducibleFlag := function( module, reducible )
    module.reducible := reducible;
end;


#############################################################################
##
#F  SetSemiLinearFlag ( <module>, <semiLinear> )
##
SetSemiLinearFlag := function( module, semiLinear )
    module.semiLinear := semiLinear;
end;


#############################################################################
##
#F  SetSizeExtensionFlag( <classic>, <size> )
##
SetSizeExtensionFlag := function( c, size )
    c.sizeExtension := size;
end;


#############################################################################
##
#F  SetSizeFlag( <classic>, <size> )
##
SetSizeFlag := function( c, size )
    c.size := size;
end;


#############################################################################
##
#F  SetSizeQuotientFlag( <next>, <o> )
##
SetSizeQuotientFlag := function( next, o )
    next.sizeQuotient := o;
end;


#############################################################################
##
#F  SetSubbasisFlag ( <module>, <subbasis> )
##
SetSubbasisFlag := function( module, subbasis )
    module.subbasis := subbasis;
end;


#############################################################################
##
#F  SetSuccessiveStripFlag( <next>, <n> )
##
SetSuccessiveStripFlag := function( next, n )
    next.successiveStrip := n;
    if KernelFlag(next) = "unknown" then
        return;
    fi;
    SetSuccessiveStripFlag(next,n);
end;


#############################################################################
##
#F  SetSymTensorBasisFlag ( <module>, <symtensorbasis> )
##
SetSymTensorBasisFlag := function( module, symtensorbasis )
    module.symTensorBasis := symtensorbasis;
end;


#############################################################################
##
#F  SetSymTensorFactorsFlag ( <module>, <symtensorfactors> )
##
SetSymTensorFactorsFlag := function( module, symtensorfactors )
    module.symTensorFactors := symtensorfactors;
end;


#############################################################################
##
#F  SetSymTensorPermFlag ( <module>, <symtensorperm> )
##
SetSymTensorPermFlag := function( module, symtensorperm )
    module.symTensorPerm := symtensorperm;
end;


#############################################################################
##
#F  SetSymTensorProductFlag( <module>, <symtensorproduct> )
##
SetSymTensorProductFlag := function( module, symtensorproduct )
    module.symTensorProduct := symtensorproduct;
end;


#############################################################################
##
#F  SetTensorBasisFlag ( <module>, <basis> )
##
SetTensorBasisFlag := function( module, basis )
    module.tensorBasis := basis;
end;


#############################################################################
##
#F  SetTensorFactorsFlag ( <module>, <factors> )
##
SetTensorFactorsFlag := function( module, factors )
    module.tensorFactors := factors;
end;


#############################################################################
##
#F  SetTensorProductFlag ( <module>, <tensorProduct> )
##
SetTensorProductFlag := function( module, tensorProduct )
    module.tensorProduct := tensorProduct;
end;


#############################################################################
##
#F  SetTypeFlag( <next>, <string> )
##
SetTypeFlag := function( next, string )
    next.type := string;
end;


#############################################################################
##
#F  SetUnitaryFormFlag( <classic>, <form> )
##
SetUnitaryFormFlag := function( c, form )
    c.unitaryForm := form;
end;


#############################################################################
##

#E  flags.g . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
