###########################################################################
##
#A  init.g                   autag package                    Werner Nickel
##
##  Januar 1997
##
##  This file is part of the share package version (distributed with GAP
##  3.4.4) of the autag package, which contains GAP code for computing
##  automorphism groups of finite soluble groups given by special soluble
##  group presentations. 
## 
###########################################################################
PrintPkgInit(rec(name:="autag",date:=[1997,1]));

LOCALNAMEautag := Cat(LOADED_PACKAGES.autag, "autag/" );

ReadLocalDir := function ( dirname, name )
    local   readIndent;
    readIndent := ReadIndent;
    ReadIndent := ConcatenationString( ReadIndent, "  " );
    InfoRead1( "#I",ReadIndent,"ReadLib( \"", name, "\" )\n" );
    if not READ( ConcatenationString( dirname, name, ".g" ) )  then
        Error("the library file '",name,"' must exist and be readable");
    fi;
    ReadIndent := readIndent;
    if ReadIndent = ""  then
        InfoRead1( "#I  ReadLib( \"", name, "\" ) done\n" );
    fi;
end;

AUTO(ReadPkg("autag","autag/misc.g"),
IsSagGroup, Select, PrintVec, PrintMat, NiceMatList,
Nice, Plural, Seconds, Time, IsNonZeroElt, IsZeroElt,
RenamedGensSagGroup);

LOCALNAMEmod := Cat(LOCALNAMEautag,"mod/");

AUTO( ReadLocalDir( LOCALNAMEmod, "initial" ),
  GModOps, GModRecOps, InfoDecompose, InfoHom, InfoStack );

AUTO( ReadLocalDir( LOCALNAMEmod, "modaut" ),
  CentraliserMatGroup, CentralizerMatGroup );

AUTO( ReadLocalDir( LOCALNAMEmod, "modhom" ),
  SpinHom );

AUTO( ReadLocalDir( LOCALNAMEmod, "modrec" ),
  Gmodule );

LOCALNAMEaut := Cat(LOCALNAMEautag,"aut/");

AUTO( ReadLocalDir( LOCALNAMEaut, "initial" ),
  InfoOrbit, InfoOrbit2, InfoAutgroup, InfoAut, AutGroupOps, AutOps,
  AutGroupElements, OrbitOps, PairOps );

AUTO( ReadLocalDir( LOCALNAMEaut, "autgrp" ),
  AutGroupSagGroup, AutGroupFactors, AutGroupStructure, AutGroupSeries,
  AutGroupConverted );

AUTO( ReadLocalDir( LOCALNAMEaut, "autrec" ),
  IsAut, AutGroupElementsOps, Aut, InnerAut, IsInnerAut );

AUTO( ReadLocalDir( LOCALNAMEaut, "eqns" ),
  EqnOps, InfoEqns, InfoEqns2 );

AUTO( ReadLocalDir( LOCALNAMEaut, "nonsplit" ),
  InfoLift );

if IsBound(USERISMJS) then
    InfoOrbit := Print;
    InfoOrbit2 := Ignore;
    InfoAutgroup := Print;
    
    InfoIntertwine := Ignore;
    InfoDecompose := Ignore;
    InfoModuleHom := Ignore;
    InfoModuleHomTiming := Ignore;

    QUIET := function ()
        InfoOrbit := Ignore;
        InfoOrbit2 := Ignore;
        InfoAutgroup := Ignore;
    end;

    LOUD := function ()
        InfoOrbit := Print;
        InfoOrbit2 := Ignore;
        InfoAutgroup := Print;
    end;

fi;
