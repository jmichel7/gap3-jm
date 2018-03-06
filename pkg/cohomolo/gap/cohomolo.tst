SizeScreen( [ 70, 24 ] );
#>[ 70, 24 ]

TestDir:= "../gapdata/";
#>"../gapdata/"

InfoCohomology:= Print;
#>function (...) internal; end

RequirePackage( "cohomolo" );
Read( Concatenation( TestDir, "d8" ) );
chr:= CHR(G,2,F,m2);;
M:= SchurMultiplier( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>[ 2 ]

D:= CoveringGroup( chr );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
Size( D );
#>16

F:= FirstCohomologyDimension( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>2

F:= SecondCohomologyDimension( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>3

E:= SplitExtension( chr );;
Size(E);
#>32

E:= NonsplitExtension( chr );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
Size(E);
#>32

E:= NonsplitExtension( chr, [0,1,1] );;
Size(E);
#>32

Read( Concatenation( TestDir, "a6" ) );
chr:= CHR( G, 3, F, m3 );;
M:= SchurMultiplier( chr );
#>#Indices in the subgroup chain are:  10 4 
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>[ 3 ]

D:= CoveringGroup( chr );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
Size(D);
#>1080

F:= SecondCohomologyDimension( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>2

E:= SplitExtension( chr );;
SE:= Subgroup( E, [E.4,E.5,E.6,E.7,E.8] );;
Index( E, SE );
#>1080

E:= NonsplitExtension( chr, [1,2] );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
SE:= Subgroup( E, [E.4,E.5,E.6,E.7,E.8] );;
Index( E, SE );
#>1080

