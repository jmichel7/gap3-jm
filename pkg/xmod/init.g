##############################################################################
##
#A  init.g                    XMOD library                    version 13/ 1/97
##
##  This file is read by GAP on issuing the command RequirePackage( "xmod" );
##
#H  $Log: init.g,v $
#H  Revision 1.1  1997/03/27 13:33:46  gap
#H  Adding xmod
#H
#H  	SL
#H
##
PrintPkgInit(rec(name:="xmod",version:="1.3.1",date:=[1997,1],
  copyright:="(C) Chris Wensley, Murat Alp -- Crossed-modules, Cat1-groups"));

##############################################################################
##
##  xmod_dir . . . . . . . . . . . . . subdirectory of .../pkg  for xmod files
##
##  Change the initial part of the directory name in the following command
##  to point to the appropriate subdirectory on your own system:

xmod_dir := LOADED_PACKAGES.xmod;

###############################################################################
##  Crossed module record for  X  :-
##
##                         /--------->  X.source
##                         |                |
##     AUT(X.source) >=  X.aut              |
##                         ^                |  X.boundary
##                         |                |
##               X.action  |                V
##                         \----------- X.range
##
###############################################################################
##  cat1-group record for  C  :-
##                                      C.embedKernel
##                /--------->  C.source <------------ Kernel(C.tail) ~ C.kernel
##                |             |    |                                     |
##                |             |    |                          C.boundary |
##  C.embedRange  ^      C.tail |    | C.head              ( ~ restriction V
##  (= embedding  |             |    |                     (   of C.head   |
##  ( of C.range  |             |    |                     (   to C.kernel |
##  (in C.source  |             V    V                                     |
##                \----<-----  C.range  <---------------------------------/
##
##############################################################################
#H   HISTORY
#H   =======
#H
#H  13/ 1/97  Version 1.3.1 posted on the ftp server
#H   2/12/96  Extended Cat1List to include groups of size <= 47
#H  28/11/96  Rewrote PermGroupAutoGroup as AutomorphismPair
#H  27/11/96  Added: XModMorphismAutoPerm
#H  26/11/96  Added: intermediate dispatcher functions Derivations, Sections
#H  20/11/96  Major revision to derivation record: now stored as:
#H            XModDerivationByImages  with appropriate opertions record.
#H            All derivation, section and actor functions modified.
#H  13/11/96  Added: InducedXMod when P->Q is a surjection
#H   4/11/96  Added: InducedCat1Group using free product G *_R Q
#H  18/10/96  Moved all dispatcher functions to file:  dispatch.g
#H
#H  30/ 8/96  Version 1.2 posted on the ftp server
#H  29/ 8/96  GreedyCommonTransversal replaced by CommonTransversal
#H  28/ 8/96  Added:  DistinctRepresentatives, CommonRepresentatives
#H  23/ 8/96  Shifted  CommonTransversal  functions to  utils.g
#H            Added:  FactorXMod  and modified NormalSubXMods
#H  16/ 8/96  Added: InclusionMorphism & InverseMorphism to Cat1GroupOps
#H  25/ 7/96  Used modified EndomorphismClasses in AllCat1Groups
#H            Added: XModOps.Parent & .IsParent  and fields .parent, .isParent
#H  19/ 7/96  Add: DirectProductCat1Group, TrivialActionXMod, DirectProductXMod
#H            New: XModOps.AutomorphismGroup via direct product of autS,autR
#H  18/ 7/96  Added:  InnerDerivationGroup
#H  16/ 7/96  Moved utility functions from to file util.g
#H  13/ 7/96  Added new  Cat1GroupXMod  with a perm group as source
#H  10/ 7/96  Major rewrite of Derivations record
#H            Regular derivations now come first in the tables.
#H   9/ 7/96  Made use of the new  T.protected  Tietze field in InducedXMod
#H   8/ 7/96  Added:  CentralExtensionXMod  and  IsCentralExtension
#H   3/ 7/96  Renamed:  XModNormalInclusion  as  ConjugationXMod  etc.
#H            Added:  InnerAutomorphismGroup, InnerAutomorphismXMod
#H  23/ 5/96  Now calls TzPartition in InducedXMod before calling TzGo
#H  14/ 5/96  Added: ReverseCat1Group, ReverseIsomorphismCat1Group
#H   7/ 5/96  Added: RegularSectionsCat1Group, AllSectionsCat1Group
#H  30/ 4/96  Used PresentationViaCosetTable in InducedXMod
#H
#H  15/12/95  Version 1.1 posted on the ftp server
#H   7/12/95  Added:  XModMorphismOps.Kernel, IsRModule, XModRModule
#H  30/11/95  Existing stand-alone code for induced xmods converted to
#H            functions:  InducedXModData, InducedXMod, AllInducedXMods
#H  22/11/95  Added:  InnerDerivation  and  ListInnerDerivations
#H  21/11/95  Added:  IdentitySubXMod.  Modified some function names.
#H  17/11/95  Added to XModMorphismOps: .IsMonomorphism, .IsEndomorphism, etc.
#H  14/11/95  Added: DerivationSection, SectionDerivation, SectionsCat1Group
#H   3/11/95  Now 2 cat1-group embeddings: C.embedRange & C.embedKernel
#H  27/10/95  Added constructor  Cat1Group(S,t,h) & adjusted IsCat1Group
#H  20/10/95  Added XMod<-->Cat1Group functor & ActorCat1Group
#H  14/ 8/95  Attempted to ignore case X.action trivial in IsXMod
#H  11/ 8/95  Used Tuples in XModDerivations to select images for genrng
#H  14/ 7/95  Added:  XModInnerMorphism, XModCenter, XModInnerActor
#H  23/ 6/95  Added:  XModDerivationAction, XModActor
#H  15/ 6/95  Added:  XModOps and XModMorphismOps
#H  14/ 6/95  Added functions for automorphisms of xmods
#H            Added:  SourceEndomorphismDerivation, XModMorphismDerivation
#H  26/ 5/95  Renamed  XModPermNormalSubgroup  as  XModNormalInclusion
#H  25/ 5/95  Added SubXMod, NormalSubXMods, XModSize, XModElements, 
#H  11/ 5/95  First functions for derivations
#H   1/ 4/94  first GAP version for induced xmods converted from CAYLEY
#H
##############################################################################
##
##
#V  XModPrintLevel . . . . . . . . . variable for controlling detailed output
##
XModPrintLevel := 1;

##  Temporary commands which will not be required when the updated version
##  of  fptietze.g  is available:
ReadLib( "fpgrp" );
ReadLib( "fpsgpres" );
ReadPkg( "xmod", "lib/pact.g" );
ReadPkg( "xmod", "lib/mytz.g" );

## needed by  PresentationViaCosetTable:
MappingsOps.Permutation := GroupElementsOps.Permutation;

#F  AUTO( ReadPkg( "xmod", "lib/dispatch" )
AUTO( ReadPkg( "xmod", "lib/dispatch" ),
      InclusionMorphism, ZeroMorphism, IsAutomorphismGroup, 
      AutomorphismPermGroup, InnerAutomorphismGroup,
      IsConjugation, IsAspherical, IsSimplyConnected, IsCentralExtension,
      IsTrivialAction, IsAutomorphismXMod, IsZeroBoundary, IsRModule,
      ActorSquareRecord, WhiteheadPermGroup,
      Whitehead, Lue, Norrie, Actor, InnerActor, InnerMorphism );

#F  AUTO( ReadPkg( "xmod", "lib/util" )
AUTO( ReadPkg( "xmod", "lib/util" ),
      GroupOpsInclusionMorphism, GroupOpsZeroMorphism,
      AutomorphismPair, IsAutomorphismPair, GroupOpsIsAutomorphismGroup,
      GroupOpsAutomorphismPermGroup, GroupOpsInnerAutomorphismGroup,
      EndomorphismClasses, EndomorphismImages, IdempotentImages,
      IsFpPair, IsSemidirectPair, PairIsomorphism,
      RegularFpPair, MinTransitiveFpPair, MinBitransitiveFpPair,
      MinTransitiveCyclicFpPair, SmallFpPair, FpPairPermGroup,
      FpPair, SemidirectPair, PrintList, DistinctRepresentatives,
      CommonRepresentatives, CommonTransversal, IsCommonTransversal );

#F  AUTO( ReadPkg( "xmod", "lib/xmod" )
AUTO( ReadPkg( "xmod", "lib/xmod" ),
      XMod, IsXMod, XModPrint, XModName, ConjugationXMod, CentralExtensionXMod,
      AutomorphismXMod, InnerAutomorphismXMod, TrivialActionXMod,
      IsRModuleRecord, RModuleXMod, XModSelect, IdentitySubXMod,
      SubXMod, IsSubXMod, IsNormalSubXMod, NormalSubXMods, FactorXMod,
      XModOps, WhatTypeXMod );

#F  AUTO( ReadPkg( "xmod", "lib/xmodmor" )
AUTO( ReadPkg( "xmod", "lib/xmodmor" ),
      XModMorphism, IsXModMorphism, XModMorphismPrint, XModMorphismName,
      XModMorphismOps, ImageXModMorphism, SourceXModXPModMorphism,
      XModMorphismAutoPerm );

#F  AUTO( ReadPkg( "xmod", "lib/cat1" )
AUTO( ReadPkg( "xmod", "lib/cat1" ),
      Cat1, IsCat1, Cat1Print, Cat1Name, ConjugationCat1, ReverseCat1,
      Cat1Ops, XModCat1, Cat1XMod, SemidirectCat1XMod, Cat1Select,
      SubCat1, IdentitySubCat1, NormalSubCat1s, RepSubCat1s, AllCat1s );

#F  AUTO( ReadPkg, "xmod", "lib/cat1mor" )
AUTO( ReadPkg( "xmod", "lib/cat1mor" ),
      Cat1Morphism, IsCat1Morphism, Cat1MorphismName, Cat1MorphismPrint,
      Cat1MorphismOps, Cat1MorphismSourceHomomorphism, ReverseIsomorphismCat1,
      Cat1MorphismXModMorphism, XModMorphismCat1Morphism,
      Cat1InclusionMorphism, AreIsomorphicCat1s );

#F  AUTO( ReadPkg, "xmod", "lib/cat1list" )
AUTO( ReadPkg( "xmod", "lib/cat1list" ),
      Cat1List, Cat1ListMaxSize, NumbersOfIsomorphismClasses );

#F  AUTO( ReadPkg, "xmod", "lib/deriv" );
AUTO( ReadPkg( "xmod", "lib/deriv" ),
      Cat1SectionByImagesOps, Cat1SectionsOps, Cat1SectionByImages,
      XModDerivationByImagesOps, XModDerivationsOps, XModDerivationByImages,
      IsSection, IsDerivation, DerivationSection, SectionDerivation,
      CompositeSection, CompositeDerivation, AreSections, AreDerivations,
      GenerationOrder, CheckGenerationPairs,
      DerivationTable, DerivationImages, DerivationImage,
      BacktrackDerivationsJ, BacktrackDerivations, DerivationsSorted,
      BacktractSectionsJ, BacktrackSections, SectionsByEndoClasses,
      Sections, RegularSections, AllSections,
      Derivations, RegularDerivations, AllDerivations,
      WhiteheadGroupTable, WhiteheadMonoidTable,
      WhiteheadPermGroup, InverseDerivations, ListInverseDerivations,
      InnerDerivation, ListInnerDerivations, InnerDerivationGroup,
      SourceEndomorphismDerivation, TableSourceEndomorphismDerivations,
      RangeEndomorphismDerivation, TableRangeEndomorphismDerivations,
      XModEndomorphismDerivation, ImageAutomorphismDerivation,
      SourceEndomorphismSection, RangeEndomorphismSection,
      Cat1EndomorphismSection );

#F  AUTO( ReadPkg, "xmod", "lib/induce" )
AUTO( ReadPkg( "xmod", "lib/induce" ),
      InducedXModData, InducedXModByCopower, InducedXMod, AllInducedXMods,
      InducedCat1Data, InducedCat1ByFreeProduct, InducedCat1,
      AllInducedCat1s, IsomorphicXMod );
