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
#H  Revision 1.2  1997/01/20 15:06:13  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.41  1996/06/07  09:35:51  eminkes
#H  Removed UpperBoundCoveringRadiusJanwa, because it was functionally
#H  identical to UpperBoundCoveringRadiusGriesmerLike.
#H
#H  Revision 1.40  1996/05/23  12:01:18  eminkes
#H  Changed the banner to include my name...
#H
#H  Revision 1.39  1996/04/26  10:31:40  eminkes
#H  Changed banner, and startup variables.
#H
#H  Revision 1.38  1996/04/26  10:25:52  eminkes
#H  Made changes that were necessary for the new release of GUAVA.
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
if not IsBound( InfoCoveringRadius ) then InfoCoveringRadius := Print; fi;
if not IsBound( InfoMinimumDistance ) then InfoMinimumDistance := Ignore; fi;
if not IsBound( CRMemSize ) then CRMemSize := 2^15; fi;

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
                  "codefun", "decoders", "codegen", "codeman", "nordrob",
                  "codecr", "codecstr", "codenorm", "codemisc", "util2" ];
    elif Length(arg) > 1 or IsString( arg[1] ) then
        files := arg;
    else
        files := arg[1];
    fi;
    for i in files do
        ReadPkg( "guava", "lib", i );
    od;
end;
#F  AUTO( ReadPkg( "guava", "lib/codeops" )
AUTO( ReadPkg( "guava", "lib/codeops" ),
      IsCode, CodeOps, LinCodeOps, CycCodeOps, WordLength,
      Redundancy, IsLinearCode, GeneratorMat, CheckMat, IsCyclicCode,
      GeneratorPol, CheckPol, MinimumDistance, WeightDistribution,
      InnerDistribution, OuterDistribution, Decode,
      IsSelfDualCode, SyndromeTable, StandardArray,
      IsSelfOrthogonalCode, CodeIsomorphism, RootsOfCode, Syndrome,
      DistancesDistribution, CodewordNr, MinimumWeightWords, History,
      LowerBoundMinimumDistance, UpperBoundMinimumDistance );
#F  AUTO( ReadPkg( "guava", "lib/codeman" )
AUTO( ReadPkg( "guava", "lib/codeman" ),
      ExpurgatedCode, AddedElementsCode, PuncturedCode, ShortenedCode,
      RemovedElementsCode, LengthenedCode, MergeNames, DualCode,
      DirectSumCode, DirectProductCode, UUVCode, PermutedCode, ExtendedCode,
      ConstantWeightSubcode, AugmentedCode, EvenWeightSubcode, ResidueCode,
      StandardFormCode, UnionCode, IntersectionCode,  ConversionFieldCode,
      CosetCode, ConcatenationCode, ConstructionBCode );
#F  AUTO( ReadPkg( "guava", "lib/codefun" )
AUTO( ReadPkg( "guava", "lib/codefun" ),
      GuavaToLeon, WeightHistogram, MergeHistories );
#F  AUTO( ReadPkg( "guava", "lib/bounds" )
AUTO( ReadPkg( "guava", "lib/bounds" ),
      BoundsTableRangeQ, BoundsTableRangeN, BoundsTableRangeT, BoundsTable,
      UpperBoundHamming, UpperBoundSingleton, UpperBoundPlotkin,
      UpperBoundGriesmer, UpperBoundElias, UpperBoundJohnson, UpperBound,
      UpperMinimumDistanceBound, OptimalityCode, OptimalityLinearCode,
      BoundsMinimumDistance, IsPerfectCode, IsMDSCode);
#F  AUTO( ReadPkg( "guava", "lib/codeword" )
AUTO( ReadPkg( "guava", "lib/codeword" ),
      IsCodeword, Codeword, VectorCodeword, PolyCodeword, WeightCodeword, 
      Support, TreatAsVector, TreatAsPoly, CodewordOps, NullWord, 
      DistanceCodeword );
#F  AUTO( ReadPkg( "guava", "lib/decoders" ), BCHDecoder, HammingDecoder )
AUTO( ReadPkg( "guava", "lib/decoders" ), BCHDecoder, HammingDecoder );
#F  AUTO( ReadPkg( "guava", "lib/matrices" )
AUTO( ReadPkg( "guava", "lib/matrices" ), 
      KrawtchoukMat, GrayMat, SylvesterMat, HadamardMat, IsLatinSquare,
      AreMOLS, MOLS, AllWords, VerticalConversionFieldMat,
      HorizontalConversionFieldMat, IsInStandardForm, PutStandardForm );
#F  AUTO( ReadPkg( "guava", "lib/nordrob" ), NordstromRobinsonCode )
AUTO( ReadPkg( "guava", "lib/nordrob" ), NordstromRobinsonCode );
#F  AUTO( ReadPkg( "guava", "lib/util" )
AUTO( ReadPkg( "guava", "lib/util" ),
      SphereContent, Krawtchouk, PermutedCols, ReciprocalPolynomial,
      CyclotomicCosets, PrimitiveUnityRoot, RemoveFiles, NullVector,
      IntPPFFE );
#F  AUTO( ReadPkg( "guava", "lib/codegen" )
AUTO( ReadPkg( "guava", "lib/codegen" ),
      ElementsCode, RandomCode, HadamardCode, ConferenceCode, MOLSCode,
      GeneratorMatCode, CheckMatCode, RandomLinearCode, HammingCode,
      ReedMullerCode, LexiCode, GreedyCode, AlternantCode, GoppaCode, 
      CordaroWagnerCode, GeneralizedSrivastavaCode, SrivastavaCode, 
      ExtendedBinaryGolayCode, ExtendedTernaryGolayCode, BestKnownLinearCode,
      GeneratorPolCode, CheckPolCode, RepetitionCode, WholeSpaceCode,
      CyclicCodes, NrCyclicCodes, BCHCode, ReedSolomonCode, RootsCode,
      QRCode, NullCode, FireCode, BinaryGolayCode, TernaryGolayCode,
      SimplexCode );

AUTO( ReadPkg( "guava", "lib/codecr" ),
      CoveringRadius, BoundsCoveringRadius,
      IncreaseCoveringRadiusLowerBound, ExhaustiveSearchCoveringRadius,
      GeneralLowerBoundCoveringRadius, GeneralUpperBoundCoveringRadius,
      LowerBoundCoveringRadiusSphereCovering,
      LowerBoundCoveringRadiusVanWee1,
      LowerBoundCoveringRadiusVanWee2,
      LowerBoundCoveringRadiusEmbedded1,
      LowerBoundCoveringRadiusEmbedded2, 
      LowerBoundCoveringRadiusCountingExcess,
      LowerBoundCoveringRadiusInduction,
      UpperBoundCoveringRadiusRedundancy,
      UpperBoundCoveringRadiusDelsarte,
      UpperBoundCoveringRadiusGriesmerLike,
      UpperBoundCoveringRadiusCyclicCode,
      UpperBoundCoveringRadiusForGriesmerCodes,
      UpperBoundCoveringRadiusStrength );

AUTO( ReadPkg( "guava", "lib/codecstr" ),
      AmalgatedDirectSumCode, BlockwiseDirectSumCode,
      ExtendedDirectSumCode, PiecewiseConstantCode,
      GabidulinCode, EnlargedGabidulinCode, DavydovCode, TombakCode,
      EnlargedTombakCode );

AUTO( ReadPkg( "guava", "lib/codemisc" ),
      CodeWeightEnumerator, CodeDistanceEnumerator, CodeMacWilliamsTransform,
      WeightVector, RandomVector, IsSelfComplementaryCode,
      IsAffineCode, IsAlmostAffineCode, IsGriesmerCode, CodeDensity,
      AlternativeOuterDistribution );

AUTO( ReadPkg( "guava", "lib/codenorm" ),
      CoordinateSubCode, CoordinateNorm, CodeNorm, IsCoordinateAcceptable,
      GeneralizedCodeNorm, IsNormalCode );

AUTO( ReadPkg( "guava", "lib/util2" ),
      AllOneVector, AllOneCodeword, IntCeiling, IntFloor,
      KroneckerDelta, BinaryRepresentation,
      SortedGaloisFieldElements );
