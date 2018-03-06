PrintPkgInit(rec(name:="matrix",version:="1.0",copyright:=
  ["(C) Frank Celler, Derek Holt, Charles Leedham-Green, Alice Niemeyer",
   "    Eamon O'Brien, Cheryl Praeger, Anthony Pye, Sarah Rees"]));

#F  set up a record for revision numbers ####################################
RevisionMatrix := rec();


#F  these are defined later #################################################
Smash          := rec();
MatrixMTX      := Smash;
IsGModule      := Error;
CRecSL         := rec();

#############################################################################
##

#F  fix problems in the GAP 3.4.3 library
##
#AUTO(ReadPkg( "matrix", "fixes-3.4", "dispatch" ), BlowupVec, Rewrite);
# dispatch also defines MatGroupOps.IsPrimitive
ReadPkg( "matrix", "fixes-3.4", "dispatch" );
AUTO(ReadPkg( "matrix", "fixes-3.4", "matrix"   ),ProjectiveOrderMat);

ReadDataPkg := function( arg )
    local   ind,  fln,  i;

    # store old indent value, add two spaces
    ind := ReadIndent;
    ReadIndent := ConcatenationString( ReadIndent, "  " );

    # construct complete path
    fln := Copy( LOADED_PACKAGES.(arg[1]) );
    for i  in [ 2 .. Length(arg)-1 ]  do
    	Append( fln, arg[i] );
    	Add( fln, '/' );
    od;
    Append( fln, arg[Length(arg)] );
    IsString(fln);
    InfoRead1( "#I", ReadIndent, "ReadPkg( \"", fln, "\" )\n" );

    # read in file
    if not READ(fln)  then
	Error("share library file \"",fln,"\" must exist and be readable");
    fi;

    # restore old indentation
    ReadIndent := ind;

end;


#############################################################################
##
#F  general library, flags
##
AUTO( ReadPkg( "matrix", "lib", "random" ),
      InitPseudoRandom,
      PseudoRandom, 
      RandomSeedFlag,
      SetRandomSeedFlag
);


AUTO( ReadPkg( "matrix", "lib", "orbhash" ),
      PermGroupRepresentation
);


AUTO( ReadPkg( "matrix", "lib", "matrix" ),
      DisplayMat
);


AUTO( ReadPkg( "matrix", "lib", "recclass" ),
      ClassicalForms,
      ConstructivelyRecogniseClassical,
      ConstructivelyRecognizeClassical,
      RecogniseClassical,
      RecogniseClassicalCLG,
      RecognizeClassical,
      RecognizeClassicalCLG
);


AUTO( ReadPkg( "matrix", "lib", "module" ),
      DualFrobeniusGModule,
      GModule,
      IsGModule
);


AUTO( ReadPkg( "matrix", "lib", "flags" ),
      AbsolutelyReducibleFlag,
      AbstractGeneratorsFlag,
      AlgElCharPolFacFlag,
      AlgElCharPolFlag,
      AlgElFlag,
      AlgElMatFlag,
      AlgElNullspaceDimensionFlag,
      AlgElNullspaceVecFlag,
      AssignLayersFlag,
      AssignLayersVecFlag,
      BasisFlag,
      BasisSubmoduleFlag,
      BlockFlag,
      BlockSystemFlag,
      CentMatFlag,
      CentMatMinPolyFlag,
      ClassicalTypeFlag,
      DegreeFieldExtFlag,
      DimensionAboveFlag,
      DimensionFlag,
      DimensionFlag,
      DimensionQuotientFlag,
      DuamFormFlag,
      ExtraSpecialFlag,
      ExtraSpecialGroupFlag,
      ExtraSpecialPrimeFlag,
      FieldFlag,
      FpGroupFlag,
      FpHomomorphismFlag,
      FrobeniusAutomorphismsFlag,
      GeneratorsFlag,
      IdentityBlockFlag,
      IdentityFlag,
      IdentityQuotientFlag,
      ImprimitiveFlag,
      InvariantFormFlag,
      InvariantFormFlag,
      IsFaithfulFlag,
      IsGenericFlag,
      IsOrthogonalGroupFlag,
      IsPossibleImprimitiveFlag,
      IsPossibleNormalizerPGroupFlag,
      IsPossibleSemiLinearFlag,
      IsPossibleSmallerFieldFlag,
      IsPossibleTensorPowerFlag,
      IsPossibleTensorProductFlag,
      IsSLContainedFlag,
      IsSymplecticGroupFlag,
      IsUnitaryGroupFlag,
      KernelFlag,
      LayerDimensionsFlag,
      LayerNumberFlag,
      LayersFlag,
      LayersVecFlag,
      LinearPartFlag,
      MapsFlag,
      MaximumStripFlag,
      NumberBlocksFlag,
      PGroupFlag,
      PermDomainFlag,
      PermGroupFlag,
      PermGroupPFlag,
      PossibleAlmostSimpleFlag,
      PossibleAlternatingGroupsFlag,
      PossibleChevalleyGroupsFlag,
      PossibleImprimitiveDimensionsFlag,
      PossibleNearlySimpleFlag,
      PossibleOverLargerFieldFlag,
      PossibleSmallerFieldFlag,
      PossibleSporadicGroupsFlag,
      PossibleTensorDimensionsFlag,
      PrimeFlag,
      PrintLevelFlag,
      QuadraticFormFlag,
      QuotientFlag,
      RecogniseFlag,
      ReducibleFlag,
      SemiLinearFlag,
      SetAbsolutelyReducibleFlag,
      SetAbstractGeneratorsFlag,
      SetAlgElCharPolFacFlag,
      SetAlgElCharPolFlag,
      SetAlgElFlag,
      SetAlgElMatFlag,
      SetAlgElNullspaceDimensionFlag,
      SetAlgElNullspaceVecFlag,
      SetBasisFlag,
      SetBasisSubmoduleFlag,
      SetBlockSystemFlag,
      SetCentMatFlag,
      SetCentMatMinPolyFlag,
      SetClassicalTypeFlag,
      SetDegreeFieldExtFlag,
      SetDimensionAboveFlag,
      SetDimensionFlag,
      SetDimensionQuotientFlag,
      SetDualFormFlag,
      SetExtraSpecialFlag,
      SetExtraSpecialGroupFlag,
      SetExtraSpecialPrimeFlag,
      SetFieldFlag,
      SetFpGroupFlag,
      SetFpHomomorphismFlag,
      SetFrobeniusAutomorphismsFlag,
      SetGeneratorsFlag,
      SetIdentityBlockFlag,
      SetIdentityFlag,
      SetIdentityQuotientFlag,
      SetImprimitiveFlag,
      SetInvariantFormFlag,
      SetInvariantFormFlag,
      SetIsFaithfulFlag,
      SetIsOrthogonalGroupFlag,
      SetIsPossibleImprimitiveFlag,
      SetIsPossibleNormalizerPGroupFlag,
      SetIsPossibleSemiLinearFlag,
      SetIsPossibleSmallerFieldFlag,
      SetIsPossibleTensorPowerFlag,
      SetIsPossibleTensorProductFlag,
      SetIsSLContainedFlag,
      SetIsSymplecticGroupFlag,
      SetIsUnitaryGroupFlag,
      SetKernelFlag,
      SetLayerDimensionsFlag,
      SetLayerNumberFlag,
      SetLayersFlag,
      SetLayersVecFlag,
      SetLinearPartFlag,
      SetMaximumStripFlag,
      SetPGroupFlag,
      SetPermDomainFlag,
      SetPermGroupPFlag,
      SetPossibleAlmostSimpleFlag,
      SetPossibleAlternatingGroupsFlag,
      SetPossibleChevalleyGroupsFlag,
      SetPossibleImprimitiveDimensionsFlag,
      SetPossibleSmallerFieldFlag,
      SetPossibleSporadicGroupsFlag,
      SetPossibleTensorDimensionsFlag,
      SetPrimeFlag,
      SetPrintLevelFlag,
      SetQuadraticFormFlag,
      SetQuotientFlag,
      SetReducibleFlag,
      SetSemiLinearFlag,
      SetSizeExtensionFlag,
      SetSizeFlag,
      SetSizeQuotientFlag,
      SetSubbasisFlag,
      SetSuccessiveStripFlag,
      SetSymTensorBasisFlag,
      SetSymTensorFactorsFlag,
      SetSymTensorPermFlag,
      SetSymTensorProductFlag,
      SetTensorBasisFlag,
      SetTensorFactorsFlag,
      SetTensorProductFlag,
      SetTypeFlag,
      SetUnitaryFormFlag,
      SizeExtensionFlag,
      SizeFlag,
      SizeQuotientFlag,
      SubbasisFlag,
      SuccessiveStripFlag,
      SymTensorBasisFlag,
      SymTensorFactorsFlag,
      SymTensorPermFlag,
      SymTensorProductFlag,
      TensorBasisFlag,
      TensorFactorsFlag,
      TensorProductFlag,
      TypeFlag,
      UndoAbsolutelyReducibleFlag,
      UndoAbstractGeneratorsFlag,
      UndoBasisFlag,
      UndoBasisSubmoduleFlag,
      UndoCentMatFlag,
      UndoDegreeFieldExtFlag,
      UndoDimensionAboveFlag,
      UndoDimensionQuotientFlag,
      UndoFpHomomorphismFlag,
      UndoGeneratorsFlag,
      UndoIdentityBlockFlag,
      UndoIdentityFlag,
      UndoIdentityQuotientFlag,
      UndoKernelFlag,
      UndoMaximumStripFlag,
      UndoPermDomainFlag,
      UndoPermGroupPFlag,
      UndoQuotientFlag,
      UndoReducibleFlag,
      UndoSuccessiveStripFlag,
      UndoTensorBasisFlag,
      UndoTensorFactorsFlag,
      UndoTensorProductFlag,
      UnitaryFormFlag
);


#############################################################################
##
#F  classical group recognition, Niemeyer/Praeger version
##
AUTO( ReadPkg( "matrix", "classic.np", "classical" ),
      GenericParameters,
      InitRecog,
      IsGenericNearlySimple,
      IsAlternating,
      IsExtensionField,
      IsGeneric,
      IsMatthieu,
      IsPSL,
      IsReducible,
      RecogniseClassicalNP,
      RecogniseClassicalNPCase,
      RecognizeClassicalNP,
      RuledOutExtFieldParameters,
      SetNotAbsIrredFlags,
      SetReturnNPFlags,
      TestRandomElement
);


#############################################################################
##
#F  special linear group recognition, Celler/Leedham-Green version
##
AUTO( ReadPkg( "matrix", "classic.clg", "recsl", "recsl" ),
      RecSL,
      RecogniseSL,
      RecognizeSL,
      SporadicGroupsInfo
);


AUTO( ReadPkg( "matrix", "classic.clg", "recsl", "crecsl" ),
      CRecSL,
      CRecogniseSL,
      CRecognizeSL
);


AUTO( ReadPkg( "matrix", "classic.clg", "recsl", "chevgrp" ),
      ChevA,
      ChevB,
      ChevC,
      ChevD,
      ChevE,
      ChevF,
      ChevG,
      Chev2A,
      Chev2B,
      Chev2D,
      Chev2E,
      Chev2F,
      Chev2G,
      Chev3D
);


AUTO( ReadPkg( "matrix", "classic.clg", "recsp", "recsp" ),
      RecSP,
      RecogniseSP,
      RecognizeSP
);


AUTO( ReadPkg( "matrix", "classic.clg", "recsp", "crecsp" ),
      CRecSP,
      CRecognizeSP
);


AUTO( ReadPkg( "matrix", "classic.clg", "recsu", "recsu" ),
      RecSU,
      RecogniseSU,
      RecognizeSU
);


AUTO( ReadPkg( "matrix", "classic.clg", "recso", "recso" ),
      RecSO,
      RecogniseSO,
      RecognizeSO
);


AUTO( ReadPkg( "matrix", "classic.clg", "recso", "recso0" ),
      RecSO0,
      RecognizeSO0
);


AUTO( ReadPkg( "matrix", "classic.clg", "recso", "recsom" ),
      RecSOm,
      RecognizeSOm
);


AUTO( ReadPkg( "matrix", "classic.clg", "recso", "recsop" ),
      RecSOp,
      RecognizeSOp
);


#############################################################################
##
#F  utility function for "classic.clg"
##
AUTO( ReadPkg( "matrix", "classic.clg", "util", "retree" ),
      PrintTree,
      RexpTree,
      SizeTree
);


AUTO( ReadPkg( "matrix", "classic.clg", "util", "semifacs" ),
      SemiFactorsInt,
      SemiPrimePowersInt
);


AUTO( ReadPkg( "matrix", "classic.clg", "util", "cycred" ),
      InitCRWord,
      NextCRWord
);


AUTO( ReadPkg( "matrix", "classic.clg", "util", "ppd" ),
      IsPpdElement
);


AUTO( ReadPkg( "matrix", "classic.clg", "util", "classic" ),
      BlockMat,
      DiagonalMat_mtx,
      EichlerTransformation,
      Gamma,
      GeneralOrthogonalGroup,
      LeftNullspaceMat,
      O,
      RightNullspaceMat,
      SPwithForm,
      SpinorNorm,
      WallForm
);


#############################################################################
##
#F  Smash
##
AUTO( ReadPkg( "matrix", "smash",  "absreduc" ),
      CompleteBasis,
      FieldGenCentMat,
      FrobeniusAction,
      UndoAbsolutelyIrreducibleFlags
);


AUTO( ReadPkg( "matrix", "smash",  "cdualmod" ),
      DualGModule
);


AUTO( ReadPkg( "matrix", "smash",  "cinduced" ),
      StrongGenImages,
      EltImage,
      StrongGenImagesCall,
      InducedGModule,
      PermGModule
);


AUTO( ReadPkg( "matrix", "smash",  "ctensorp" ),
      TensorProductGModule,
      WedgeGModule
);


AUTO( ReadPkg( "matrix", "smash",  "cwrthprd" ),
      ImprimitiveWreathProduct,
      PowerWreathProduct,
      PowerPerm
);


AUTO( ReadPkg( "matrix", "smash",  "choose" ),
      AddRandomTranslatingConjugate,
      ChooseRandomElements,
      Commutators,
      ElementOfOrder,
      ElementWithCharPol,
      LargestPrimeOrderElement,
      LargestPrimePowerOrderElement,
      RandomConjugate
);


AUTO( ReadPkg( "matrix", "smash",  "composit" ),
      Distinguish,
      MinimalSubGModule
);


AUTO( ReadPkg( "matrix", "smash",  "extraspl" ),
      ExtraSpecialDecomposition,
      Stripped
);


AUTO( ReadPkg( "matrix", "smash",  "homomorp" ),
      HomGModule,
      IsomorphismGModule,
      MatrixSum,
      MinimalSubGModules,
      SortHomGModule
);

AUTO( ReadPkg( "matrix", "smash",  "minblcks" ),
      AmalgamateBlocks,
      BlockImage,
      CountBlocks,
      EquateBlocks,
      FilterVector,
      InitialiseBlock,
      InverseColumn,
      MinBlocks,
      Rep,
      SetImage,
      SetRep,
      SetupBasis,
      SetupBlockPermGroup,
      StandardVector
);


AUTO( ReadPkg( "matrix", "smash",  "misc" ),
      FactorsToInt,
      IntPower,
      InverseMod,
      IsScalar
);


AUTO( ReadPkg( "matrix", "smash",  "polyfdeg" ),
      FactorsPolDeg,
      FactorsSquarefreePolDeg
);


AUTO( ReadPkg( "matrix", "smash",  "isprimit" ),
      AddSmashQueue,
      BasicReductionTests, 
      BlockNumbers,
      BlockSizes,
      CallSmashGModule,
      CharPolPrimeOrder,
      CharPolPrimePowerOrder,
      CharPolStructure,
      CompositeOrders,
      ConstructBlock,
      DeletePrimitivityComponents,
      ExamineSmashResult,
      ExponentGL,
      FindLargestPower,
      FinishComputation,
      FreeRank,
      GcdOrderGL,
      GcdSeq,
      IndexMinimumRankElement,
      InverseSet,
      IsBlockSystem,
      IsOrderValid,
      IsValidSymOrder,
      LcmSeq,
      OrderOfElement,
      PolynomialCoefficients,
      PolynomialQuotient,
      PolynomialRemainder,
      PrimitiveTest,
      ReportResult,
      ResolveTensor,
      SemiLinearTest,
      SetBlockNumbers,
      SetBlockSizes,
      SettleComputation,
      SmashElement,
      StartPrimitivityTest
);


AUTO( ReadPkg( "matrix", "smash",  "reducibl" ),
      EnlargeIrreducibleGModule,
      GoodElementGModule,
      InducedAction,
      OrthogonalVector,
      QuotientGModule,
      RandomIrreducibleSubGModule,
      SpinBasis,
      SubGModule,
      SubGModuleAction
);


AUTO( ReadPkg( "matrix", "smash",  "semilinr" ),
      SemiLinearDecomposition,
      PowerMaps,
      IsSemiLinear
);


AUTO( ReadPkg( "matrix", "smash",  "smash" ),
      SemiSimpleDecomposition,
      SmashGModule,
      TranslatesDirectSum,
      TranslatesIrreducible,
      TranslatesSGModules
);


AUTO( ReadPkg( "matrix", "smash",  "stabilis" ),
      BlockStabiliserTest,
      ChooseFirstElement,
      ExamineCompositionFactors,
      ExtractLargestPrimeOrderElement,
      ExtractLargestPrimePowerOrderElement,
      FindExpressions,
      FindMinimalSubGModules,
      FindSolutions,
      FixedPointFreeElement,
      InvariantsOfCF,
      ProcessLattice,
      ProcessSubGModule,
      SampleOfCommutators,
      SampleOfElements,
      SetupElements,
      SetupFactor,
      SomePowerFixesBlock,
      SortCFs,
      TestSolution
);


AUTO( ReadPkg( "matrix", "smash",  "system" ),
      BaseRingZero,
      LargestMovedPoint,
      OrderKnownDividend,
      PolCoefficients,
      PrimitiveElement, 
      Root
);


AUTO( ReadPkg( "matrix", "smash",  "tensor" ),
      KroneckerFactors,
      MultipleTensorProductDecomposition,
      SwapFactors,
      SymTensorProductDecomposition,
      TensorProductDecomposition,
      UndoTensorProductFlags
);


#T  1996/12/25 fceller we must read the following files in order to
#T                     set the entries in 'Smash'

ReadPkg( "matrix", "smash", "absreduc" );
ReadPkg( "matrix", "smash", "composit" );
ReadPkg( "matrix", "smash", "isprimit" );
ReadPkg( "matrix", "smash", "reducibl" );


#############################################################################
##
#F  Tensor
##

AUTO( ReadPkg( "matrix", "tensor",  "direct" ),
      DirectSumSpaces
);


AUTO( ReadPkg( "matrix", "tensor",  "facpoly" ),
      ApplySymmetry,
      ComputeTensorTable,
      DecideFactorisation,
      ExponentsOfFactors,
      FactorisePolynomials,
      FindFactorisation,
      IsSystemSolvable,
      ListFactors,
      NonNegativeSolution,
      PolynomialTensorProduct,
      PowerOfSmallOrder,
      ProcessVector,
      SetupMatrices
);


AUTO( ReadPkg( "matrix", "tensor",  "findpnt" ),
      AreProportional,
      BasisMatrix,
      BlocksOfMatrix,
      ComponentsOfSum,
      ConstructMatrices,
      ConstructNewFlat,
      FindIsom,
      FindPoint,
      FoundSingularElement,
      GeneralFindPoint,
      IdentifySubset,
      InvestigateMatrices,
      IsPoint,
      ProjectivitiesGenerateField,
      SearchForSingularElement,
      SetupBlocks,
      SetupMatrix
);


AUTO( ReadPkg( "matrix", "tensor",  "isprojec" ),
      FactorsFlag, 
      GFlag, 
      IsProjectivity, 
      ProcessElement,
      ProjectivityTest,
      SetFactorsFlag, 
      SetGFlag
);


AUTO( ReadPkg( "matrix", "tensor",  "istensor" ),
      ConstructTensorFactors,
      IsTensor,
      TensorTest
);


AUTO( ReadPkg( "matrix", "tensor",  "local" ),
      BestCompLength,
      CompositionSeriesLength,
      ComputeLattice,
      HasShortLattice,
      InvestigateLattice,
      LocalSubgroup,
      LocalSubgroups,
      LocalTest,
      SubmoduleLattice,
      SubmoduleLatticeAbort
);

AUTO( ReadPkg( "matrix", "tensor",  "misc" ),
      CentralisedMod,
      CentralisedSpace,
      CentralisedSpaceSet,
      CentralisedSpaceSetMod,
      ComplementSpace,
      ComplementSpaceSet,
      DeletePair,
      ElementCommutes,
      EliminateRepetitions,
      FirstNonZeroEntry,
      HashSet,
      HashSpace,
      IsPowerOfPolynomial,
      LogQ,
      NullSpace,
      ProjectiveOrder,
      Remove,
      Swap
);


AUTO( ReadPkg( "matrix", "tensor",  "numbers" ),
      CoPrimeFactorisations,
      DistinctPrimes,
      FactorList,
      InverseSet,
      IsAPower,
      IsPrimePower,
      PairOfFactors,
      ProperDivisors
);


AUTO( ReadPkg( "matrix", "tensor",  "order" ),
      ExistsFactorisation,
      FindBestPrime,
      LeastLinearSemiSimple,
      LeastProjective,
      LeastProjectiveSemiSimple,
      OrderTest,
      PossibleFactorisation,
      Score,
      SetPartitions
);


AUTO( ReadPkg( "matrix", "tensor",  "pdash" ),
      EmbedLargerField,
      RestrictSmallerField,
      pDash,
      pDashLocalElements,
      pDashLocals,
      pDashSubgroup
);


AUTO( ReadPkg( "matrix", "tensor",  "plocal" ),
      CommutingElement,
      ComputeSpaces,
      ExtractBlock,
      IspNormal,
      SubQuotAction,
      pLocal,
      pLocalElements,
      pLocalSubgroups
);


AUTO( ReadPkg( "matrix", "tensor",  "stabilis" ),
      GenerateElements,
      GetSeed,
      SetSeed,
      SpaceStabiliser,
      StabiliserOfSet
);


#############################################################################
##
#F  Reducible
##
AUTO( ReadPkg( "matrix", "reduce", "kernel" ),
      ApproximateKernel,
      ReduceLayerNumberByOne
);


AUTO( ReadPkg( "matrix", "reduce", "niceprnt" ),
      DisplayMatRecord,
      PrintBasicInfo,
      PrintLayer,
      PrintLayers,
      PrintNumberOfLayers
);


AUTO( ReadPkg( "matrix", "reduce", "split" ),
      ExponentsLayer,
      ExtractQuotientBlock,
      FinishEnlargingQuotient,
      GetBlocks,
      GoDownChain,
      GoDownGeneratorsRels,
      GoDownKernelPGroup,
      GoDownLargerQuotientSL,
      GoDownNewSubmodule,
      GoDownPGroup,
      GoDownPerm,
      GoDownRandomRels,
      GoDownSL,
      GoDownUnknown,
      InitSplit,
      InitialiseKernel,
      InitialiseQuotient,
      IsFaithfulMtx,    # jm 4/2016 renamed from IsFaithful (conflict with lib)
      LayerMat,
      ProcessClassicResult,
      ProcessPrimitiveResult,
      RandomRelsSL,
      RecogniseMatrixGroup,
      RecognizeMatrixGroup,
      RemoveRepeats,
      SetupPermRep,
      SplitMatGroup,
      StartEnlargingQuotient,
      TensorPowerTest,
      TidyMatRecord
);


AUTO( ReadPkg( "matrix", "reduce", "system" ),
      BaseRing, 
      GetFpHomomorphism,
      RandomRelsPerm
);


AUTO( ReadPkg( "matrix", "reduce", "wraps" ),
      EvaluateRelation,
      RandomRelation
);
