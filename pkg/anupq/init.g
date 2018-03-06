PrintPkgInit(rec(name:="anupq",version:="1.4"));

## setting 'ANUPQlog' to 'LogTo' will create a log file "debug.out" #########
ANUPQlog := function( str )
    InfoRead1 := Ignore;
    ANUPQlog  := Ignore;
end;
#ANUPQlog := LogTo;

## functions used from inside the ANU pga ###################################
AUTO( ReadPkg( "anupq", "gap", "anustab" ),
      ANUPQoutputResult, ANUPQstabilizer );

## functions to use the ANU pga from inside GAP #############################
AUTO( ReadPkg( "anupq", "gap", "anupga" ),
      ANUPQauto, ANUPQautoList, SavePqList, PqList, PqDescendants,
      DimensionFrattiniFactor );

## functions to use the ANU pq from inside GAP ##############################
AUTO( ReadPkg( "anupq", "gap", "anupq" ), Pq, PqHomomorphism );

## function to use the ANU standard presentation from inside GAP ############
AUTO( ReadPkg( "anupq", "gap", "anusp" ), StandardPresentation,
      IsIsomorphicPGroup, IsomorphismPcpStandardPcp, AutomorphismsPGroup );
