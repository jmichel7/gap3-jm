#############################################################################
##
##  This file contains the test functions for the classical Recognition code.
##
##
RequirePackage("matrix");
ReadDataPkg( "matrix", "tst", "nsimple.g" );
RecogniseClassicalNP;
InfoRecog2 := Ignore;
InfoRecog1 := Print;


##
##  Linear groups
##
Print("#I Testing Linear Groups\n");
failures := 0;
for q in [ 2, 3, 5, 4, 9, 27 ] do
	for d in [ 10, 15 .. 60 ] do
	        Print( "d = ", d, " q = ", q, "\n");
		grp := GL( d,  q );
	        x := RecogniseClassicalNP( grp, "linear", 100 );
	        if x = false then failures := failures + 1; fi;
	od;
od;
if failures > 0 then 
    Print( "#I Didn't pass linear groups test, ", failures, " failures\n");
else
    Print( "#I Passed linear groups test\n ");
fi;


##
##  Symplectic groups
##

Print("\n#I Testing Symplectic Groups\n");
failures := 0;
for q in [ 2, 3, 5, 4, 9, 27 ] do
	for d in [ 10, 16 .. 40  ] do
	        Print( "d = ", d, " q = ", q, "\n");
		grp := SP( d,  q );
	        x := RecogniseClassicalNP( grp, "symplectic", 100 );
	        if x = false then failures := failures + 1; fi;
	od;
od;
if failures > 0 then 
    Print( "#I Didn't pass symplectic groups test, ", failures, " failures\n");
else
    Print( "#I Passed symplectic groups test\n ");
fi;


##
##  Unitary Groups
##
Print("\n#I Testing Unitary Groups\n");
failures := 0;
for q in [ 2, 3, 5, 4  ] do
	for d in [ 10, 15 .. 40  ] do
	        Print( "d = ", d, " q = ", q, "\n");
		grp := GU( d,  q );
	        x := RecogniseClassicalNP( grp, "unitary", 100 );
	        if x = false then failures := failures + 1; fi;
	od;
od;
if failures > 0 then 
    Print( "#I Didn't pass unitary groups test, ", failures, " failures\n");
else
    Print( "#I Passed unitary groups test\n ");
fi;
 
##
##  Orthogonal groups
##
Print("\n#I Testing Orthogonal+ Groups\n");
failures := 0;
for q in [ 2, 3, 5, 4  ] do
	for d in [ 10, 16 .. 28  ] do
 	      Print( "d = ", d, " q = ", q, "\n");
 	      grp := O( 1, d,  q );
              x := RecogniseClassicalNP( grp, "orthogonalplus", 100 );
 	      if x = false then failures := failures + 1; fi;
 	od;
od;
if failures > 0 then 
    Print( "#I Didn't pass O+ groups test, ", failures, " failures\n");
else
    Print( "#I Passed O+ groups test\n ");
fi;
 
Print("\n#I Testing Orthogonal- Groups\n");
failures := 0;
for q in [ 2, 3, 5, 4  ] do
        for d in [ 10, 16 .. 28  ] do
 	        Print( "d = ", d, " q = ", q, "\n");
 		grp := O( -1, d,  q );
 	        x := RecogniseClassicalNP( grp, "orthogonalminus", 100 );
 	        if x = false then failures := failures + 1; fi;
 	od;
od;
if failures > 0 then 
     Print( "#I Didn't pass O- groups test, ", failures, " failures\n");
else
     Print( "#I Passed O- groups test\n ");
fi;
 
Print("\n#I Testing Orthogonal o Groups\n");
failures := 0;
for q in [ 2, 3, 5, 4  ] do
 	for d in [ 11, 17 .. 29  ] do
 
 	        Print( "d = ", d, " q = ", q, "\n");
 		grp := O( 0, d,  q );
 	        x := RecogniseClassicalNP( grp, "orthogonalcircle", 100 );
 	        if x = false then failures := failures + 1; fi;
 	od;
od;
 
if failures > 0 then 
     Print( "#I Didn't pass Oo groups test, ", failures, " failures\n");
else
     Print( "#I Passed Oo groups test\n ");
fi;


InfoRecog2 := Print;

Print("\n#I Testing non-classical groups\n");
failures := 0;
for grp in [ a7d4q2, m11d5q3, 2m12d6q3, m23d11q2_1, m23d11q2_2, m24d11q2_1,
             m24d11q2_2, psl27d3q11_1, psl27d3q11_2 ] do
                Print("Calling RecogniseClassicalNP\n");
  	        x := RecogniseClassicalNP( grp, "linear", 30 );
  	        if x <> false then failures := failures + 1; fi;
                Print("Calling IsGenericNearlySimple\n");
                x := IsGenericNearlySimple (grp, "linear", 30 );
  	        if x <> true then failures := failures + 1; fi;

od;
if failures > 0 then 
    Print("#I Didn't pass nearly simple groups test, ",failures," failures\n");
else
    Print( "#I Passed nearly simple groups test\n ");
fi;
  

