#############################################################################
##
#A  init.g                  GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file is read by GAP upon startup. It displays the banner and then
##  defines all the functions in the Guava library as autoreadable
##
#H  $Log: init.g,v $
#H  Revision 1.1  1997/01/20 15:15:54  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.37  1994/11/10  11:11:24  jcramwin
#H  removed CreateBoundsTable
#H
#H  Revision 1.36  1994/11/10  11:03:54  jcramwin
#H  added CreateBoundsTable
#H
#H  Revision 1.35  1994/11/10  10:56:51  jcramwin
#H  removed ComplexTable
#H
#H  Revision 1.34  1994/11/03  08:56:20  jcramwin
#H  moved IsPerfectCode, IsMDSCode from codefun to bounds
#H
#H  Revision 1.33  1994/11/02  14:20:52  jcramwin
#H  removed the infolines in the function guava() and now guava can be
#H  called with a string or list of strings with names of files
#H
#H  Revision 1.32  1994/11/02  12:17:22  jcramwin
#H  added the function History in codeops.g
#H
#H  Revision 1.31  1994/11/02  10:33:32  jcramwin
#H  Canged DescriptionForCode in HistoryOfCode
#H
#H  Revision 1.30  1994/10/24  08:35:26  jcramwin
#H  LowerBoundMinimumDistance, UpperBoundMinimumDistance in codeops.g
#H
#H  Revision 1.29  1994/10/20  16:04:09  jcramwin
#H  DescriptionForCode added
#H
#H  Revision 1.28  1994/10/20  14:25:32  rbaart
#H  Added MinimumWeightWords
#H
#H  Revision 1.27  1994/10/20  14:11:02  jcramwin
#H  changed OptimalLinearCode in BestKnownLinearCode
#H
#H  Revision 1.26  1994/10/19  11:42:57  jcramwin
#H  minor layout change
#H
#H  Revision 1.25  1994/10/19  11:34:19  jcramwin
#H  SimplexCode Added
#H
#H  Revision 1.24  1994/10/17  10:35:46  rbaart
#H  Corrected some filenames ( codeword(s) and nordrob(.g) )
#H
#H  Revision 1.23  1994/10/14  08:57:49  rbaart
#H  Added BoundsMinimumDistance
#H
#H  Revision 1.22  1994/10/14  08:17:19  rbaart
#H  Added GUAVA_REF_LIST as a global variable
#H
#H  Revision 1.21  1994/10/13  15:30:13  jcramwin
#H  changed Distance to DistanceCodeword, Weight to WeightCodeword, Poly to
#H  PolyCodeword and Vector to VectorCodeword
#H
#H  Revision 1.20  1994/10/12  14:17:23  rbaart
#H  Changed banner
#H
#H  Revision 1.19  1994/10/12  14:01:23  rbaart
#H  Removed Guava's PR function
#H
#H  Revision 1.18  1994/10/12  11:57:15  rbaart
#H  Changed initial value of GUAVA_BOUNDS_TABLE
#H
#H  Revision 1.17  1994/10/12  08:24:41  rbaart
#H  Changed banner
#H
#H  Revision 1.16  1994/10/12  07:46:50  rbaart
#H  Added our names in the banner
#H
#H  Revision 1.15  1994/10/11  11:56:48  rbaart
#H  Added ConstructionBCode
#H
#H  Revision 1.14  1994/10/10  16:31:25  jcramwin
#H  changed SimpleTable in ComplexTable
#H
#H  Revision 1.13  1994/10/04  12:48:32  rbaart
#H  Added ResidueCode
#H
#H  Revision 1.12  1994/09/30  13:08:03  rbaart
#H  Changed OptimalCode to OptimalLinearCode
#H
#H  Revision 1.11  1994/09/30  12:27:39  rbaart
#H  Minor layout change
#H
#H  Revision 1.10  1994/09/30  11:12:17  rbaart
#H  Added OptimalCode
#H
#H  Revision 1.9  1994/09/29  11:25:57  rbaart
#H  Added ConcatenationCode
#H
#H  Revision 1.8  1994/09/29  10:36:35  rbaart
#H  Show progress of unit reading status
#H
#H  Revision 1.7  1994/09/29  10:29:00  rbaart
#H  Show progress of reading units
#H
#H  Revision 1.6  1994/09/29  10:26:34  rbaart
#H  *** empty log message ***
#H
#H  Revision 1.5  1994/09/29  10:09:02  rbaart
#H  Function guava() prints read units
#H
#H  Revision 1.4  1994/09/29  10:05:46  jcramwin
#H  SimpleTable added
#H
#H  Revision 1.3  1994/09/28  13:39:49  jcramwin
#H  IntFFE patch removed
#H
#H  Revision 1.2  1994/09/28  11:07:12  rbaart
#H  Inserted header
#H
#H  Revision 1.1  1994/09/28  09:44:04  jcramwin
#H  Initial revision
#H
##
if BANNER and not QUIET then
    Print("\n");
    # Use the PR function from GAP's init.g
    PR("   ___________                   |                  ");
    PR("  /            \\           /   --+--  Version 1.2   ");
    PR(" /      |    | |\\\\        //|    |                  ");
    PR("|    _  |    | | \\\\      // |                        ");
    PR("|     \\ |    | |--\\\\    //--|     Jasper Cramwinckel  ");
    PR(" \\     ||    | |   \\\\  //   |     Erik Roijackers     ");
    PR("  \\___/  \\___/ |    \\\\//    |     Reinald Baart       ");
fi;

#############################################################################
##
#V  GUAVA_TEMP_VAR . . . . . . . . variable for interfacing external programs
##
GUAVA_TEMP_VAR := 0;

#############################################################################
##
#V  GUAVA_BOUNDS_TABLE . . . . . . . . .  contains a list of tables of bounds
##
GUAVA_BOUNDS_TABLE := [ [], [] ];

#############################################################################
##
#V  GUAVA_REF_LIST  . . . contains a record of references for bounds on codes
##
GUAVA_REF_LIST := rec();

#############################################################################
##
#F  guava() . Temporary file that reads all files, for debugging purposes only
##
guava := function(arg)
    local i, files;
    if Length(arg) = 0 then
        files := ["util", "matrices", "bounds", "codeword", "codeops",
                  "codefun", "decoders", "codegen", "codeman", "nordrob" ];
    elif Length(arg) > 1 or IsString( arg[1] ) then
        files := arg;
    else
        files := arg[1];
    fi;
    for i in files do
        Read(Flat(Concatenation(PKGNAME,"guava/lib/",i,".g")));
    od;
end;
#F  AUTO( ReadPkg( "guava", "lib/codeops" )
AUTO( ReadPkg( "guava", "lib/bounds" ),
      UpperBoundHamming, UpperBoundSingleton, UpperBoundPlotkin,
      UpperBoundGriesmer, UpperBoundElias, UpperBoundJohnson, UpperBound,
      IsPerfectCode, IsMDSCode, OptimalityCode, OptimalityLinearCode,
      BoundsOps, BoundsMinimumDistance );

AUTO( ReadPkg( "guava", "lib/codecr" ),
      BoundsCoveringRadius, CoveringRadiusSearch, CoveringRadius );

AUTO( ReadPkg( "guava", "lib/codecstr" ),
      AmalgatedDirectSum, BlockwiseDirectSumCode, ExtendedDirectSumCode );

AUTO( ReadPkg( "guava", "lib/codefun" ),
      GuavaToLeon, WeightHistogram, MergeHistories );

AUTO( ReadPkg( "guava", "lib/codegen" ),
      ElementsCode, RandomCode, HadamardCode, ConferenceCode, MOLSCode,
      GeneratorMatCode, CheckMatCode, RandomLinearCode, HammingCode,
      SimplexCode, ReedMullerCode, LexiCode, GreedyCode, AlternantCode,
      GoppaCode, CordaroWagnerCode, GeneralizedSrivastavaCode, SrivastavaCode,
      ExtendedBinaryGolayCode, ExtendedTernaryGolayCode, BestKnownLinearCode,
      GeneratorPolCode, CheckPolCode, RepetitionCode, WholeSpaceCode,
      CyclicCodes, NrCyclicCodes, BCHCode, ReedSolomonCode, RootsCode, QRCode,
      NullCode, FireCode, BinaryGolayCode, TernaryGolayCode );

AUTO( ReadPkg( "guava", "lib/codeman" ),
      DualCode, AugmentedCode, EvenWeightSubcode, ConstantWeightSubcode,
      ExtendedCode, ShortenedCode, PuncturedCode, ExpurgatedCode,
      AddedElementsCode, RemovedElementsCode, LengthenedCode, ResidueCode,
      ConstructionBCode, PermutedCode, StandardFormCode, ConversionFieldCode,
      CosetCode, DirectSumCode, UnrestrictedDirectSumCode,
      LinearDirectSumCode, ConcatenationCode, DirectProductCode, UUVCode,
      UnionCode, IntersectionCode );

AUTO( ReadPkg( "guava", "lib/codemisc" ),
      CodeWeightEnumerator, CodeDistanceEnumerator, CodeMacWilliamsTransform,
      WeightVector, IsSelfComplementaryCode, AllOneVector, AllOneCodeword,
      IntCeiling, IntFloor, KroneckerDelta, SemiStandardForm );

AUTO( ReadPkg( "guava", "lib/codenorm" ),
      SubCoordinateCode, SubCode, CoordNorm, CodeNorm, CoordinateAcceptable,
      GeneralizedCodeNorm, IsNormalCode );

AUTO( ReadPkg( "guava", "lib/codeops" ),
      CodeOps, LinCodeOps, CycCodeOps, IsCode, WordLength, IsLinearCode,
      Redundancy, GeneratorMat, CheckMat, IsCyclicCode, GeneratorPol,
      CheckPol, MinimumDistance, LowerBoundMinimumDistance,
      UpperBoundMinimumDistance, MinimumWeightWords, WeightDistribution,
      InnerDistribution, OuterDistribution, Decode, IsSelfDualCode,
      SyndromeTable, StandardArray, AutomorphismGroup, IsSelfOrthogonalCode,
      CodeIsomorphism, RootsOfCode, DistancesDistribution, Syndrome,
      CodewordNr, History );

AUTO( ReadPkg( "guava", "lib/codeprob" ),
      ProbMinimumDistance, FindCoveringRadiusAtConstantWeight,
      FindCoveringRadius );

AUTO( ReadPkg( "guava", "lib/codeword" ),
      IsCodeword, Codeword, VectorCodeword, PolyCodeword, DistanceCodeword,
      WeightCodeword, Support, TreatAsVector, TreatAsPoly, CodewordOps,
      NullWord );

AUTO( ReadPkg( "guava", "lib/crbounds" ),
      LowerBoundCoveringRadiusVanWee, LowerBoundCoveringRadiusEmbedded1,
      LowerBoundCoveringRadiusEmbedded2, LowerBoundCoveringRadiusExcess1,
      LowerBoundCoveringRadiusExcess2, LowerBoundCoveringRadiusExcess3,
      LowerBoundCoveringRadiusInduction1, LowerBoundCoveringRadiusInduction2,
      LowerBoundCoveringRadiusInduction3,
      LowerBoundCoveringRadiusSphereCovering, GeneralLowerBoundCoveringRadius,
      UpperBoundCoveringRadiusDelsarte, GeneralUpperBoundCoveringRadius );

AUTO( ReadPkg( "guava", "lib/decoders" ),
      BCHDecoder, HammingDecoder );

AUTO( ReadPkg( "guava", "lib/matrices" ),
      KrawtchoukMat, GrayMat, SylvesterMat, HadamardMat, IsLatinSquare,
      AreMOLS, MOLS, VerticalConversionFieldMat, HorizontalConversionFieldMat,
      IsInStandardForm, PutStandardForm );

AUTO( ReadPkg( "guava", "lib/nordrob" ),
      NordstromRobinsonCode );

AUTO( ReadPkg( "guava", "lib/tblgener" ),
      CreateBoundsTable );

AUTO( ReadPkg( "guava", "lib/util" ),
      SphereContent, Krawtchouk, PermutedCols, ReciprocalPolynomial,
      CyclotomicCosets, PrimitiveUnityRoot, RemoveFiles, NullVector );

