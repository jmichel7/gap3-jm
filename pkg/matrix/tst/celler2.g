# timings
RecogniseClassicalNP;
InfoRecog1 := Ignore;
InfoRecog2 := Ignore;

if not IsBound(NOS)  then NOS := 2;      fi;
if not IsBound(ALL)  then ALL := false;  fi;

Print( "#I  number of test runs: ", NOS, "\n" );

#############################################################################

# a7d4q2
G := Group( 
[[0,1,0,0],
 [1,1,0,0],
 [0,0,0,1],
 [0,0,1,1]] * GF(2).one,

[[0,0,1,0],
 [0,1,1,0],
 [1,0,1,1],
 [0,1,1,1]] * GF(2).one );

r := RecogniseClassical(G,"clg");
if not 7 in r.possibleAlternatingGroups  then
    Print("#W  failed: a7d4q2 (",r.possibleAlternatingGroups,")\n");
fi;


#############################################################################
tests := [
    "a5.gap",      5,
    "a5d4.gap",    5,
    "a6c6d12.gap", 6,
    "a6d16.gap",   6,
    "a7c6d24.gap", 7,
    "a7d20.gap",   7,
    "a8d20.gap",   8,
    "a8d64.gap",   8,
    "a9d28.gap",   9,
    "2a11d55.gap", 11,
    "2a11d56.gap", 11,
    "s22d21.gap",  22,
];
for i  in [ 1, 3 .. Length(tests)-1 ]  do
    file := tests[i];
    name := tests[i+1];

    # NP
    ReadDataPkg( "matrix", "data", file );
    tmp := ClassicalForms(G);
if ALL then
    if ForAny( tmp, x -> x[1] = "unknown" )  then
        tmp := First( tmp, x -> x[1] = "unknown" );
        if Length(tmp) > 1  then
            Print( "#I  no form, ", tmp[2], "\n" );
        else
            Print( "#I  no form found\n" );
        fi;
    else
        start := 0;
        for n  in [ 1 .. NOS ]  do
            ReadDataPkg( "matrix", "data", file );
            start := start - Runtime();
            r := RecogniseClassical(G,"np");
            start := start + Runtime();
        od;
        Print( "#I  NP  ", String(file,-13), ": ",
               String(QuoInt(start,NOS),11), " msec\n" );
    fi;
fi;

    # CLG
    start := 0;
    for n  in [ 1 .. NOS ]  do
        ReadDataPkg( "matrix", "data", file );
        start := start - Runtime();
        r := RecogniseClassical(G,"clg");
        if not name in r.possibleAlternatingGroups  then
            Print("#W  failed: ",name," (",r.possibleAlternatingGroups,")\n");
        fi;
        start := start + Runtime();
    od;
    Print( "#I  CLG ", String(file,-13), ": ",
           String(QuoInt(start,NOS),11), " msec\n" );
if ALL then
    Print( "\n" );
fi;

od;


#############################################################################
tests := [
    "2co1d24.gap",  "Co1",
    "co3d22.gap",   "Co3",
    "rud28f17.gap", "Ru",
    "m12d55.gap",   "M12",
    "m11d44.gap",   "M11",
    "3m22d12.gap",  "M22",
    "m22d54.gap",   "M22",
    "suz.gap",      "Suz",
    "suzd12a.gap",  "Suz",
    "suzd12b.gap",  "Suz",
    "he2d102.gap",  "He",
    "he2d50.gap",   "He",
];

for i  in [ 1, 3 .. Length(tests)-1 ]  do
    file := tests[i];
    name := tests[i+1];

    # NP
if ALL then
    start := 0;
    for n  in [ 1 .. NOS ]  do
        ReadDataPkg( "matrix", "data", file );
        start := start - Runtime();
        r := RecogniseClassical(G,"np");
        start := start + Runtime();
    od;
    Print( "#I  NP  ", String(file,-13), ": ",
           String(QuoInt(start,NOS),11), " msec\n" );
fi;

    # CLG
    start := 0;
    for n  in [ 1 .. NOS ]  do
        ReadDataPkg( "matrix", "data", file );
        start := start - Runtime();
        r := RecogniseClassical(G,"clg");
        if not name in r.possibleAlmostSimple  then
            Print("#W  failed: ", name, " (", r.possibleAlmostSimple, ")\n");
        fi;
        start := start + Runtime();
    od;
    Print( "#I  CLG ", String(file,-13), ": ",
           String(QuoInt(start,NOS),11), " msec\n" );
if ALL then
    Print( "\n" );
fi;

od;
