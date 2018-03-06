PrintPkgInit(rec(name:="sisyphos",version:="0.8.1"));
#############################################################################
##
##  Notify the functions to be defined for {\SISYPHOS}.
##
AUTO( ReadPkg( "sisyphos", "gap", "sisgroup" ),
      OrderGL, IsCompatiblePCentralSeries, EstimateAmount,
      AgGroupNormalizedAgGroup, PrintSISYPHOSWord, PrintSisyphosInputPGroup,
      SisyphosAutomorphisms, OuterAutomorphisms,
      NormalizedAutomorphisms, NormalizedOuterAutomorphisms,
      PresentationAutomorphisms, AgNormalizedAutomorphisms,
      AgNormalizedOuterAutomorphisms, IsIsomorphic, Isomorphisms,
      CorrespondingAutomorphism, AutomorphismGroupElements,
      SisGModule, SisCohomology, SisTrivialModule, SisExtension,
      SisSplitExtension, SisDirectSumModule, SisTensorProductModule,
      SisDualModule);

AUTO( ReadPkg( "sisyphos", "gap", "sisgprin" ),
      NormalizedUnitsGroupRing );

#############################################################################
##
#V  SISYPHOS
#V  p
##
##  These are global variables.
##  A perhaps existing variable 'p' will be saved always (except if the
##  computation is interrupted during the run of one of the programs in this
##  file).
##
if not IsBound( SISYPHOS ) then
  SISYPHOS:= rec( SISISO  := false,
                  SISBOOL := false,
                  SISCODE := false,
                  SISTMEM := "1000000",
                  SISPMEM :=  "300000",
                  SISCALL := "bin/sis -q -b -s gap",
                  SISOps  := rec( Print:= function( r )
                                 Print( "rec(\n",
                                        "sizeAutG := ", r.sizeAutG, ",\n",
                                        "sizeInnG := ", r.sizeInnG, ",\n",
                                        "sizeOutG := ", r.sizeOutG, ",\n" );
                                 if IsBound( r.epimorphism ) then
                                   Print( "epimorphism := ",
                                           r.epimorphism, ",\n" );
                                 fi;
                                 Print( "generators := \n", r.generators,
                                         " )" );
                                  end ),
                  SISOpsGmodule  := rec( Print:= function( r )
                                 Print( "rec(\n",
                                        "dimension := ", r.dimension, ",\n",
                                        "matrices := ", r.matrices, ")\n" );
                                  end ),
                  SISOpsCohomology  := rec( Print:= function( r )
                                 Print( "rec(\n",
                                        "degree := ", r.degree, ",\n",
                                        "dimension := ", r.dimension, ",\n",
                                                   "cycleDimension := ", 
                                                   r.cycleDimension, ")\n" );
                                  end ),
                  inputfile  := "sisinp.xxx",
                  outputfile := "sisout.xxx"  );

fi;

if VERSYS{ [ 1 .. 5 ] } = "msdos" then

  # We do not want to swap the whole {\GAP} again just for removing files.
  # So under MS-DOS the input files for {\SISYPHOS} is always written
  # to the file 'SISYPHOS.inputfile', and the output is written to the
  # file 'SISYPHOS.outputfile'.
  # One should remove them by hand after one is sure that they are no longer
  # needed, for example after leaving {\GAP}.

  SISYPHOS.Exec     := Ignore;
  SISYPHOS.TmpName1 := function() return SISYPHOS.inputfile ; end;
  SISYPHOS.TmpName2 := function() return SISYPHOS.outputfile; end;

else

  # Under UNIX etc. we want to remove files as soon as possible,
  # and choose safe names for the temporary files.

  SISYPHOS.Exec     := Exec;
  SISYPHOS.TmpName1 := TmpName;
  SISYPHOS.TmpName2 := TmpName;

fi;

if not IsBound( p ) then p:= false; fi;
