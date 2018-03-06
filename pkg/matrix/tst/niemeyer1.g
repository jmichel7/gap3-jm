#############################################################################
##
##  This file contains the test functions for the classical Recognition code.
##  It requires a 2-3 minutes of real time.
##
##
RequirePackage ("matrix");
ReadDataPkg( "matrix", "tst", "nsimple.g" );
InfoRecog2 := Ignore;
InfoRecog1 := Print;

##
##  Linear groups
##
Print("#I Testing Linear Groups\n");
failures := 0;
for q in [  3, 5, 4,  ] do
	for d in [ 10, 30 ] do
	        Print( "d = ", d, " q = ", q, "\n");
		grp := GL( d,  q );
	        x := RecogniseClassicalNPCase( grp, "linear", 100 );
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
for q in [ 5, 4  ] do
	for d in [ 10, 16 .. 22  ] do
	        Print( "d = ", d, " q = ", q, "\n");
		grp := SP( d,  q );
	        x := RecogniseClassicalNPCase( grp, "symplectic", 100 );
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
for q in [  5, 4  ] do
	for d in [ 10, 15 .. 20  ] do

	        Print( "d = ", d, " q = ", q, "\n");
		grp := GU( d,  q );
	        x := RecogniseClassicalNPCase( grp, "unitary", 100 );
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
for q in [  3, 4  ] do
	for d in [ 10, 16 .. 22  ] do
 	      Print( "d = ", d, " q = ", q, "\n");
 	      grp := O( 1, d,  q );
              x := RecogniseClassicalNPCase( grp, "orthogonalplus", 100 );
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
for q in [ 2, 5  ] do
        for d in [ 10, 16 .. 22  ] do
 	        Print( "d = ", d, " q = ", q, "\n");
 		grp := O( -1, d,  q );
 	        x := RecogniseClassicalNPCase( grp, "orthogonalminus", 100 );
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
for q in [ 2, 3,  4  ] do
 	for d in [ 11, 17 .. 23  ] do
 
 	        Print( "d = ", d, " q = ", q, "\n");
 		grp := O( 0, d,  q );
 	        x := RecogniseClassicalNPCase( grp, "orthogonalcircle", 100 );
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
                Print("Calling RecogniseClassicalNPCase\n");
  	        x := RecogniseClassicalNPCase( grp, "linear", 30 );
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
  

