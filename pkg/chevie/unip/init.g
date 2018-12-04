AUTO(ReadChv("unip/uc"),
   CycPolUnipotentDegrees,
   FixRelativeType,
   ReflTypeOpsUnipotentCharacters,
   HasTypeOpsUnipotentCharacters,
   UnipotentCharacters,
   UnipotentCharactersOps,
   UnipotentDegrees);

AUTO(ReadChv("unip/ucl"),
   DistinguishedParabolicSubgroups, InducedLinearForm,
   BalaCarterLabels,
   HasTypeOpsUnipotentClasses, DimUnipotentClass, GreenTable,
   ICCTable, MellinValues, OrderClassesClassical, QuotientAu,
   ReflTypeOpsUnipotentClasses, SpecialPieces, SpringerSeries, 
   UnipotentClasses, UnipotentClassesOps, UnipotentValues);

AUTO(ReadChv("unip/families"),
   DrinfeldDouble,
   Family,
   FamiliesClassical,
   FamilyOps,
   FusionAlgebra,
   NrDrinfeldDouble,
   OnFamily,
   SubFamilyij);

AUTO(ReadChv("unip/unichars"),
   AlmostCharacter,
   UnipCharOps,
   UnipotentCharacter,
   IsUnipotentCharacter,
   DeligneLusztigLefschetz,
   DeligneLusztigLefschetzTable,
   DeligneLusztigCharacter,
   DeligneLusztigCharacterTable,
   HarishChandraInduction,
   LusztigInduction,
   LusztigRestriction);

AUTO(ReadChv("unip/lusztig"),
   FindCuspidalInLevi,
   FindIntSol,
   FindSeriesInParent,
   HarishChandraInductionTable,
   LusztigInductionTable);

AUTO(ReadChv("unip/unipotent"),
   IsUnipotentElement,
   UnipotentAbelianPart,
   UnipotentDecompose,
   UnipotentElement,
   UnipotentElementOps,
   UnipotentGroup,
   UnipotentGroupOps);

HasTypeOps.UnipotentCharacters:=
  function(arg)return ApplyFunc(HasTypeOpsUnipotentCharacters,arg);end;
ReflTypeOps.UnipotentCharacters:=
  function(arg)return ApplyFunc(ReflTypeOpsUnipotentCharacters,arg);end;
HasTypeOps.UnipotentClasses:=HasTypeOpsUnipotentClasses;
ReflTypeOps.UnipotentClasses:=ReflTypeOpsUnipotentClasses;
if not IsBound(CHEVIE.families) then CHEVIE.families:=rec();fi;
