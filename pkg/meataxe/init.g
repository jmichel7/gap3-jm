PrintPkgInit(rec(name:="meataxe",date:=[1996,12,11]));
#############################################################################
##
##  Notify the functions to be defined for {\MeatAxe}.
##
AUTO( ReadPkg( "meataxe", "gap", "meataxe" ),
      MeatAxe, GapObject, SplittingField );

AUTO( ReadPkg( "meataxe", "gap", "mapermut" ),
      MeatAxePermOps, IsMeatAxePerm, MeatAxePermutationsOps,
      MeatAxePermutations, OrderMeatAxePerm, MeatAxePerm );

AUTO( ReadPkg( "meataxe", "gap", "mamatrix" ),
      MeatAxeMatricesOps, MeatAxeMatOps, MeatAxeMatrices, PrintMeatAxeInput,
      IsMeatAxeMat, MeatAxeMat, OrderMeatAxeMat, ZEVPols, StringZEVPols,
      BrauerCharacterValue );

AUTO( ReadPkg( "meataxe", "gap", "mamodule" ),
      MeatAxeModuleOps, BasisMeatAxeModuleOps, 
      SemiEchelonBasisMeatAxeModuleOps, StandardBasisMeatAxeModuleOps, 
      LatticeSummands, GeneratorsSubmodule, GeneratorsSubmodules, 
      MeatAxeFactorModuleOps );

AUTO( ReadPkg( "meataxe", "gap", "mamatalg" ),
      MeatAxeMatAlgebraOps, MeatAxeMatGroupOps, RandomOrders );

#############################################################################
##
#F  InfoMeatAxe( ... ) . . . . . . . info function for the {\MeatAxe} package
##
if not IsBound( InfoMeatAxe ) then InfoMeatAxe:= Ignore; fi;

# rubbish...
if not IsBound( IsBlockMat ) then IsBlockMat:= obj -> false; fi;
