#############################################################################
##
#A  init.g                  CrystGap library                     Bettina Eick
#A                                                              Franz G"ahler
#A                                                              Werner Nickel
##
#Y  Copyright 1990-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  CrystGap - the crystallographic groups package for GAP (initialization)
##  
PrintPkgInit(rec(name:="cryst",date:=[1997]));
#############################################################################
AUTO(ReadPkg( "cryst", "lib", "cryst" ),
AddInternalBasis, AddTranslationsCrystGroup, AdjustStabilizer, AugmentedMatrix,
CheckTranslations, CheckWyckoffOrbits, CocycleInfo, CoefficientsMod,
ComplementsSG, ConjugatedCrystGroup, ConjugateInternalGenerators,
ConjugateWyckoffPositions, CrystGroup, CrystGroupOps, EliminateEquivExtensions,
FindTranslationsCrystGroup, FixedPointsModZ, FlattenMatMat, FractionModOne,
GroupExtEquations, HermiteNormalForm, IdentifyWyckoffOrbits, ImageEquivAffSpace,
InAffineSpaces, IsCrystGroup, IsSpaceGroup, IsStandardCrystGroup,
IsStandardSpaceGroup, IsSymmorphicSpaceGroup, IsWyckoffPosition,
ListOneCohomology, MakeSpaceGroup, MatJacobianMatrix,
MaximalSubgroupsRepresentatives, MaximalSubgroupsRepsKG, MaximalSubgroupsRepsSG,
MaximalSubgroupsRepsTG, NormalizeAffineSubspace, NoSmash, NullMatMat,
OneCoboundariesSG, OneCocyclesSG, OneCocyclesVector, OneCohomologySG,
OrbitAffineSpaces, PointGroup, ReducedTranslationBasis, RowEchelonForm,
RowEchelonFormT, RowEchelonFormVector, SimpleGenerators, SolutionHomEquations,
SolutionInhomEquations, SolveHomEquationsModZ, SolveInhomEquationsModZ,
SolveOneInhomEquationModZ, SpaceGroupsPointGroup, StandardTranslation,
TranslationsCrystGroup, TriangularizeMatVector, UseSmash, WyckoffBasis,
WyckoffOrbit, WyckoffPosClass, WyckoffPositionOps, WyckoffPositions,
WyckoffPositionsByStabilizer, WyckoffPositionsQClass, WyckoffSpaceGroup,
WyckoffStabilizer, WyckoffTranslation, WyPos);

AUTO(ReadPkg( "cryst", "lib", "cryst2" ),
AffineInequivalentSubgroups, AffineLift, AffineNormalizer, CentralizerElement,
CentralizerGL, CrystGroupOpHomOps, ElementsConjugacyClassInfiniteSubgroups, 
IntersectionModule, IntSolutionMat, NormalizerGL, PointGroupsBravaisClass, 
QuadFormEquations, QuadFormSpace, StabilizerInfiniteGroup, 
TranslationNormalizer, UnionModule, VectorModL);

#############################################################################
##
##  color group support
##
AUTO( ReadPkg( "cryst", "lib", "color" ),
      IsColorGroup, ColorSubgroup, ColorCosets, ColorOfElement, 
      ColorPermGroup, ColorHomomorphism, ColoredSubGroup,
      ColoredPointGroup, ColorGroup );

#############################################################################
##
##  Wyckoff lattice support
##
AUTO( ReadPkg( "cryst", "lib", "wylat" ),
      IsSubspaceAffineSubspace, WyckoffPosRelations, 
      WyckoffLatticeRecord, WyckoffLattice );
