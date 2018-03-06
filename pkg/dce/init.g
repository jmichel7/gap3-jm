#H##########################################################################
##A  dce.g      Steve Linton    1993--1995 
##
##  GAP double coset enumerator.
##
## Documentation of the data structures, invariant assumptions
## and the like is in a separate file, for sanity
##
PrintPkgInit(rec(name:="dce",date:=[1995]));

AUTO(ReadPkg("dce","dce"),
DCE, DCEAddExtraRels, DCEAddSchreier, DCECheckDK, DCECheckTab, DCEClearDeds,
DCEClearPDs, DCEColAdj, DCEColAdjCont, DCEColAdjContRegular, DCEColAdjSingle,
DCEColumn, DCEColumnRegular, DCECompBetas, DCEDedComp, DCEDefine, DCEDegree,
DCEDKToString, DCEDoMultiScan, DCEDoStacks, DCEDoubleCosetReps,
DCEDoubleCosetRepStabs, DCEDPush, DCEDPush1, DCEEmptyOne, DCEEqual,
DCEExtendBetas, DCEFollow, DCEHandlePEC, DCEHOrbits, DCEHOrbitsRegular,
DCEInfoPrint, DCEInitMBetas, DCEInitStacks, DCEIPerm, DCEIsFlippable,
DCELeftCosetReps, DCELin, DCEMultiFollow, DCEOPerm, DCEPack, DCEPerm,
DCEPerms, DCEPermString, DCEPrintRow, DCEPrintTab, DCEProcessA, DCEProcessB,
DCEProcessC, DCEProcessDed, DCEProcessRelator, DCEPush, DCEPushAll, DCERead,
DCERelGroups, DCERepWord, DCERightCosetReps, DCERin, DCERoot, DCERowOps,
DCERun, DCERunF, DCERunHavas, DCEScan, DCESetup, DCESetupCols, DCESetupEDPs,
DCESetupGens, DCESetupRels, DCESetupStrat, DCEStackA, DCEStackB, DCEStackC,
DCEStackDed, DCEStackPD, DCESubgroup, DCESuperPack, DCEUniverseOps, DCEWord,
DCEWordGroup, DCEWordOps, DCEWrite, DCEXonEl, DCEXonSG, IsDCEUniverse,
IsDCEWord, OldDCEInfoPrint, SetupSymmetricPresentation);
