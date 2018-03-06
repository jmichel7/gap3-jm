#############################################################################
##
#A  codefun.g                GUAVA                              Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains non-dispatched functions to get info of codes
##
#H  $Log: codefun.g,v $
#H  Revision 1.1  1997/01/20 15:14:53  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.23  1994/11/08  17:39:37  jcramwin
#H  fixed the same bug
#H
#H  Revision 1.22  1994/11/08  17:35:19  jcramwin
#H  fixed a bug in MergeHistories
#H
#H  Revision 1.21  1994/11/08  16:39:32  jcramwin
#H  MergeHistories now writes U, V, W instead of C1, C2, C3
#H
#H  Revision 1.20  1994/11/03  09:27:55  jcramwin
#H  fixed another bug in MergeHistories (in the future I might just check the
#H  functions I write ;-)
#H
#H  Revision 1.19  1994/11/03  08:54:40  jcramwin
#H  moved IsMDSCode and IsPerfectCode to bounds.g
#H
#H  Revision 1.18  1994/11/03  08:46:17  jcramwin
#H  renamed HistoryOfCodes in MergeHistories and fixed a bug
#H
#H  Revision 1.17  1994/11/02  14:09:24  jcramwin
#H  fixed a bug in HistoryOfCodes
#H
#H  Revision 1.16  1994/10/31  15:03:23  jcramwin
#H  replaced DescriptionForCode in HistoryOfCode
#H
#H  Revision 1.15  1994/10/24  08:06:52  jcramwin
#H  Lower and UpperBounMinimumDistance removed
#H
#H  Revision 1.14  1994/10/21  15:00:05  jcramwin
#H  UpperBoundMinimumdistance now checks for the weight of the generators
#H
#H  Revision 1.13  1994/10/21  11:39:12  jcramwin
#H  DescriptionForCode now wants codes in stead of .names
#H
#H  Revision 1.12  1994/10/21  11:22:24  jcramwin
#H  some cosmetic changes to HistoryForCode
#H
#H  Revision 1.11  1994/10/21  11:13:30  jcramwin
#H  Wrote a history for codes. We decided not to use this function until
#H  we know more about the the GAP wants things to be printed
#H
#H  Revision 1.10  1994/10/20  16:31:59  jcramwin
#H  written DescriptionForCode
#H
#H  Revision 1.9  1994/10/14  15:12:15  rbaart
#H  UpperBoundMinimumDistance: fixed bug with uppererBound
#H
#H  Revision 1.8  1994/10/14  13:14:20  rbaart
#H  UpperBoundMinimumDistance: fixed bug
#H
#H  Revision 1.7  1994/10/14  10:34:53  jcramwin
#H  used 'ApplyFunc' in *BoundMinimumDistance
#H
#H  Revision 1.6  1994/10/14  08:59:29  jcramwin
#H  changed the functions lower and upperboundminimumdistance.
#H  they both call BoundsMunimumdistance now
#H
#H  Revision 1.5  1994/10/06  16:57:06  rbaart
#H  Refined WeightHistogram
#H
#H  Revision 1.4  1994/10/04  10:56:31  rbaart
#H  LowerBoundMinimumDistance: changed boolean argument
#H
#H  Revision 1.3  1994/09/30  17:04:41  jcramwin
#H  LowerBoundMinimumDistance now works with ( <n>, <k> [, <q> ] )
#H
#H  Revision 1.2  1994/09/28  11:52:33  jcramwin
#H  nur einige kommentaren
#H
#H  Revision 1.1  1994/09/28  09:40:14  jcramwin
#H  Initial revision
#H
##

#############################################################################
##
#F  GuavaToLeon( <C>, <file> )  .  converts a code to a form Leon can read it
##
##  converts a code in Guava format to a library in a format that is readable
##  by Leon's programs.
##
GuavaToLeon := function (arg)
    local C, file, w, vector, G, coord, n, k, TheMat;
    if Length(arg) = 2 then
        C := arg[1];
        file := arg[2];
    else
        Error("usage: GuavaToLeon( C, filename )");
    fi;
    G := GeneratorMat(C);
    k := Dimension(C);
    n := WordLength(C);
    PrintTo(file, "LIBRARY code;\n");
    AppendTo(file,"code=seq(",String(Size(Field(C))),",",String(k),
            ",",String(n),",seq(\n");
    for vector in [1..k] do
        for coord in [1..n-1] do
            AppendTo(file,IntFFE(G[vector][coord]),",");
        od;
        AppendTo(file, IntFFE(G[vector][n]));
        if vector < k then
            AppendTo(file,",\n");
        fi;
    od;
    AppendTo(file, "\n));\nFINISH;");
end;
  
#############################################################################
##
#F  WeightHistogram ( <C> [, <height>] )  . . . . .  plots the weights of <C>
##
##  The maximum length of the columns is <height>. Default height is one
##  third of the screen size.
##
WeightHistogram := function(arg)
    local C, wd, max, data, height, i, j, n, spaces, nr, Scr, char;
    Scr := SizeScreen();
    if Length(arg) = 1 then
        C := arg[1];
        height := Int(Scr[2]/3);
    elif Length(arg) = 2 then
        C := arg[1];
        height := arg[2];
    else
        Error("usage: Histo(C [, height])");
    fi;
    char := "*";
    n := WordLength(C);
    if n+2 >= Scr[1] then
        Error("histogram does not fit on screen");
    elif n+2 > Int(Scr[1] / 2) then
        spaces := "";
        nr := 0;
    elif n+2 > Int(Scr[1] / 4) then
        spaces := " ";
        nr := 1;
    else
        spaces := "  ";
        nr := 2;
    fi;
    wd := WeightDistribution(C);
    max := Maximum(wd);
    if max < height then
        height := max;
    fi;
    data := List(wd, w -> Int(w/max*height));
    Print(max);
    for i in [0..n*(nr+1)-Length(String(max))] do
        Print("-");
    od;
    Print("\n");
    for i in height - [0..height-1] do
        for j in data do
            if j >= i then
                Print(Concatenation(char,spaces));
            else
                Print(Concatenation(" ",spaces));
            fi;
        od;
        Print("\n");
    od;
    for i in [0..n] do
        if wd[i+1] = 0 then
            Print("-");
        else
            Print("+");
        fi;
        for j in [2..nr+1] do
            Print("-");
        od;
        #Print("-");
    od;
    Print("\n");
    for i in [0..n] do
        Print(i mod 10,spaces);
    od;
    Print("\n ",spaces);
    for i in [1..n] do
        if i mod 10 = 0 then
            Print(Int(i / 10),spaces);
        else
            Print(" ",spaces);
        fi;
    od;
    Print("\n");
end;

#############################################################################
##
#F  MergeHistories( <C>, <S> [, <C1> .. <Cn> ] ) . . . . . .  list of strings
##
##

MergeHistories := function( arg )
    local i, his, names;
    
    if Length( arg ) > 1 then
        names := "UVWXYZ";
        his := [];
        for i in [1..Length(arg)] do
            Add(his, Concatenation( [ names[i] ], ": ", arg[i][1]) );
            Append( his, List( arg[i]{[2..Length(arg[i])]}, line -> 
                    Concatenation("   ", line ) ) );
        od;
        return his;
    else
        Error("usage: MergeHistories( <C1> , <C2> [, <C3> .. <Cn> ] )");
    fi;
end;

#[26,10,d] linear, non-cyclic code over GF(3), with d in [1..5], is
#Hamming code with r = 3
#
#(26,10,2) code over GF(2) is MOLSCode
#
#[26,10,d] cyclic code over GF(3), with d in [1..5], is
#U | U + V construction of C1 and C2
#
#[26,10,8] linear, non-cyclic code over GF(2) is punctured code
#
# Problems that are still remaining:
#    - how should a code be printed if it doesn't fit on one line
#    - how should a history be printed if it doesn't fit on one line
#
#CodeOps.Name := function(C)
#    local line;
#    line := [ Concatenation( "(", String(WordLength(C)), ",",
#                    String(Size(C)),",") ];
#    if UpperBoundMinimumDistance(C) = LowerBoundMinimumDistance(C) then
#        Append(line[1], String(MinimumDistance(C)));
#        Append(line[1], ") ");
#    else
#        Append(line[1], "d) ");
#    fi;
#    Add(line, "code " );
#    Add(line, Concatenation("over GF(", String(Size(Field(C))), ")"));
#    if UpperBoundMinimumDistance(C) <> LowerBoundMinimumDistance(C) then
#        Append( line[ Length(line) ], ", ");
#        Add( line, Concatenation( "with d in [",
#                String( LowerBoundMinimumDistance(C) ), "..",
#                String( UpperBoundMinimumDistance(C) ), "], ") );
#    else
#        Append( line[ Length(line) ], " ");
#    fi;
#    Add( line, "is ");
#    if not IsBound(C.name) then
#        C.name := ["user ", "defined ", "code"];
#    fi;
#    return Concatenation( line,  C.name  );
#end;
#
#LinCodeOps.Name := function(C)
#    local line;
#    line := [ Concatenation( "[", String(WordLength(C)), ",",
#                    String(Dimension(C)),",") ];
#    if UpperBoundMinimumDistance(C) = LowerBoundMinimumDistance(C) then
#        Append(line[1], String(MinimumDistance(C)));
#        Append(line[1], "] ");
#    else
#        Append(line[1], "d] ");
#    fi;
#    Append(line, ["linear ","code "] );
#    Add(line, Concatenation("over GF(", String(Size(Field(C))), ")"));
#    if UpperBoundMinimumDistance(C) <> LowerBoundMinimumDistance(C) then
#        Append( line[ Length(line) ], ", ");
#        Add( line, Concatenation( "with d in [",
#                String( LowerBoundMinimumDistance(C) ), "..",
#                String( UpperBoundMinimumDistance(C) ), "], ") );
#    else
#        Append( line[ Length(line) ], " ");
#    fi;
#    Add( line, "is ");
#    if not IsBound(C.name) then
#        C.name := ["user ", "defined ","linear ", "code"];
#    fi;
#    return Concatenation( line, C.name );
#end;
#
#CycCodeOps.Name := function(C)
#    local line;
#    line := [ Concatenation( "[", String(WordLength(C)), ",",
#                    String(Dimension(C)),",") ];
#    if UpperBoundMinimumDistance(C) = LowerBoundMinimumDistance(C) then
#        Append(line[1], String(MinimumDistance(C)));
#        Append(line[1], "] ");
#    else
#        Append(line[1], "d] ");
#    fi;
#    Append(line, ["cyclic ", "code "] );
#    Add(line, Concatenation("over GF(", String(Size(Field(C))), ")"));
#    if UpperBoundMinimumDistance(C) <> LowerBoundMinimumDistance(C) then
#        Append( line[ Length(line) ], ", ");
#        Add( line, Concatenation( "with d in [",
#                String( LowerBoundMinimumDistance(C) ), "..",
#                String( UpperBoundMinimumDistance(C) ), "], ") );
#    else
#        Append( line[ Length(line) ], " ");
#    fi;
#    Add( line, "is ");
#    if not IsBound(C.name) then
#        C.name := ["user ","defined ","cyclic ","code"];
#    fi;
#    return Concatenation( line, C.name );
#end;
#
#############################################################
#
#CodeOps.String := function(C)
#    return Concatenation( C.operations.Name(C) );
#end;
#
#LinCodeOps.String := CodeOps.String;
#
#CycCodeOps.String := CodeOps.String;
#
#############################################################
#
#CodeOps.Print := function(C)
#    local i;
#    #I am not sure that it does what we want it to do
#    for i in C.operations.Name(C) do
#        Print(i);
#    od;
#end;
#
#LinCodeOps.Print := CodeOps.Print;
#
#CycCodeOps.Print := CodeOps.Print;
