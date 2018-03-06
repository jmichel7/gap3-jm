#############################################################################
##
#A  codegen.g               GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains functions for generating codes
##
#H  $Log: codegen.g,v $
#H  Revision 1.2  1997/01/20 15:05:14  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.85  1996/05/23  11:06:13  eminkes
#H  Fixed some bugs in HadamardCode.
#H
#H  Revision 1.84  1996/04/26  09:46:31  eminkes
#H  Made changes that were necessary for the new release of GUAVA.
#H
#H  Revision 1.83  1996/03/21  14:37:59  eminkes
#H  Fixed the same bug in CheckMatCode (ESM)
#H
#H  Revision 1.82  1996/03/21  14:35:39  eminkes
#H  Bug: GeneratorMatCode didn't accept Codewords anymore. Fixed. (ESM)
#H
#H  Revision 1.81  1995/11/14  14:43:23  jcramwin
#H  Changed ReedMullerCode(k, r) to ReedMullerCode(r, k). I hope it won't be
#H  too much work to change the manual!
#H
#H  Revision 1.80  1995/11/14  12:41:56  jcramwin
#H  RandomLinearCode(n, 0, F) now returns the NullCode(n, F)
#H
#H  Revision 1.79  1995/11/14  11:42:05  jcramwin
#H  Inserted '(' and ')' twice. I am not even sure if they were needed.
#H
#H  Revision 1.78  1995/11/09  14:22:53  jcramwin
#H  If I change the GeneratorMatCode, I must also change
#H  the CheckMatCode!
#H
#H  Revision 1.77  1995/11/09  13:34:29  jcramwin
#H  fixed a dont-forget-the-"fi"-bug
#H
#H  Revision 1.76  1995/11/09  13:06:28  jcramwin
#H  GeneratorMatCode now includes a wordlength an gives an error
#H  if a user tries to make a nullcode with this function
#H
#H  Revision 1.75  1995/02/15  09:27:25  jcramwin
#H  ElementsCode now checks the length of it's elements
#H
#H  Revision 1.74  1995/01/19  10:18:51  jcramwin
#H  Changed CyclicCodes(n,q) in CyclicCodes(n[,k],q)
#H
#H  Revision 1.73  1995/01/06  10:45:51  jcramwin
#H  changed the errormessage in ElementsCode
#H
#H  Revision 1.72  1994/11/15  10:26:33  jcramwin
#H  started a change in CyclicCodes
#H
#H  Revision 1.71  1994/11/10  10:04:24  jcramwin
#H  fixed a bug which I just introjuced with my last change
#H
#H  Revision 1.70  1994/11/10  09:43:57  jcramwin
#H  ElementsCode didn't check for MinDist = 0
#H
#H  Revision 1.69  1994/11/09  18:30:20  rbaart
#H  CyclicCodes and FireCode: changed way of printing
#H
#H  Revision 1.68  1994/11/09  17:12:54  rbaart
#H  LexiCode: changed n into k
#H
#H  Revision 1.67  1994/11/09  13:17:33  jcramwin
#H  changed some names
#H
#H  Revision 1.66  1994/11/09  13:00:14  rbaart
#H  Changed the way a BCH code is printed
#H
#H  Revision 1.65  1994/11/08  17:09:42  jcramwin
#H  SymmetricGroup(1000) in RepetitionCode(1000) took to long
#H
#H  Revision 1.64  1994/11/08  15:34:34  rbaart
#H  Fixed grammatical bugs in name of codes
#H
#H  Revision 1.63  1994/11/04  12:28:08  jcramwin
#H  changed the names (and histories) of the extended Golay codes
#H
#H  Revision 1.62  1994/11/04  11:56:12  jcramwin
#H  BestKnownLinearCode: lowerBoundMinimumDistance updated with table entry
#H
#H  Revision 1.61  1994/11/03  11:42:32  rbaart
#H  BCHCode: extended table to n=255
#H
#H  Revision 1.60  1994/11/03  11:04:05  jcramwin
#H  changed NullVector
#H
#H  Revision 1.59  1994/11/03  11:00:04  rbaart
#H  BCHCode: added a table with known minimum distances
#H
#H  Revision 1.58  1994/11/02  15:26:34  jcramwin
#H  fixed a bug in WholeSpaceCode
#H
#H  Revision 1.57  1994/11/02  13:24:34  jcramwin
#H  fixed a bug with .name
#H
#H  Revision 1.56  1994/10/31  15:05:14  jcramwin
#H  changed all the functions for the new way of printing
#H
#H  Revision 1.55  1994/10/28  11:24:57  jcramwin
#H  lexiCode rewritten
#H
#H  Revision 1.54  1994/10/27  14:15:43  rbaart
#H  ReedSolomonCode: fixed bug in DescriptionForCode
#H
#H  Revision 1.53  1994/10/27  10:35:44  rbaart
#H  Added automorphism groups for trivial codes
#H
#H  Revision 1.52  1994/10/26  12:13:26  rbaart
#H  Added upperbound for minimum distance of BCH codes
#H
#H  Revision 1.51  1994/10/25  11:10:11  rbaart
#H  BCHCode: improved construction and changed name
#H
#H  Revision 1.50  1994/10/21  08:20:35  rbaart
#H  Golay codes: fixed bug in return
#H
#H  Revision 1.49  1994/10/20  17:38:43  rbaart
#H  RepetitionCode: changed name to lower case
#H
#H  Revision 1.48  1994/10/20  17:19:31  rbaart
#H  Changed C.name := ... in DescriptionForCode
#H
#H  Revision 1.47  1994/10/20  16:30:10  jcramwin
#H  changed C.name for RepetitionCode
#H
#H  Revision 1.46  1994/10/20  14:22:04  jcramwin
#H  changed the name OptimalLinearCode in BestKnownLinearCode
#H
#H  Revision 1.45  1994/10/19  15:58:11  jcramwin
#H  Changed some accurences where functions call GeneratorMatCode didn't give
#H  a Field as argument
#H
#H  Revision 1.44  1994/10/19  15:07:52  rbaart
#H  ReedMullerCode: fixed bug with .wordLength
#H
#H  Revision 1.43  1994/10/19  14:39:23  jcramwin
#H  GeneratorMatCode and CheckMatCode no longer take the BaseMat of
#H  the first argument (if it doesn't need to)
#H
#H  Revision 1.42  1994/10/19  10:40:05  jcramwin
#H  fixed the bug I just created
#H
#H  Revision 1.41  1994/10/19  10:33:21  jcramwin
#H  wrote a comment with the definition of B.construction and introduced a bug
#H
#H  Revision 1.40  1994/10/18  15:56:31  jcramwin
#H  changed 'beautyful' in 'beautiful'
#H
#H  Revision 1.39  1994/10/18  15:53:50  jcramwin
#H  made the most beautiful function in Guava even more beautiful
#H
#H  Revision 1.38  1994/10/17  16:12:05  rbaart
#H  OptimalLinearCode: fixed small bug
#H
#H  Revision 1.37  1994/10/17  16:08:48  rbaart
#H  OptimalLinearCode: checks for lower and upper bounds
#H
#H  Revision 1.36  1994/10/17  15:55:49  jcramwin
#H  OptimalLinearCode does a better checking of arg
#H
#H  Revision 1.35  1994/10/17  15:18:21  rbaart
#H  OptimalLinearCode: changed error message
#H
#H  Revision 1.34  1994/10/17  15:02:56  rbaart
#H  OptimalLinearCode uses the new result of BoundsMinimumDistance
#H  All code generating functions now accept an integer for the baseField
#H
#H  Revision 1.33  1994/10/17  08:00:02  rbaart
#H  ReedMullerCode: fixed bug in weight distribution for r=2
#H
#H  Revision 1.32  1994/10/14  15:04:07  rbaart
#H  Very minor changes
#H
#H  Revision 1.31  1994/10/14  11:12:46  jcramwin
#H  all code generating functions should use the functions ElementsCode
#H  or Generator*Code or Check*Code now.
#H
#H  Revision 1.30  1994/10/14  09:04:05  rbaart
#H  Changed OptimalLinearCode to BoundsMinimumDistance and moved it to bounds.g
#H
#H  Revision 1.29  1994/10/14  08:05:41  jcramwin
#H  lower and upperboundstable are now in one file
#H  GUAVA_REF_LIST is implemented
#H
#H  Revision 1.28  1994/10/13  16:21:58  rbaart
#H  What have I changed?
#H
#H  Revision 1.27  1994/10/13  15:57:32  rbaart
#H  OptimalLinearCode: changed reference in RecurseBounds
#H
#H  Revision 1.26  1994/10/13  15:37:16  rbaart
#H  OptimalLinearCode: fixed small bug
#H
#H  Revision 1.25  1994/10/13  15:04:28  rbaart
#H  Changed codeword functions
#H
#H  Revision 1.24  1994/10/12  15:03:00  rbaart
#H  OptimalLinearCode: new conventions in numbering of operations
#H
#H  Revision 1.23  1994/10/12  12:43:30  rbaart
#H  Changes in OptimalLinearCode
#H
#H  Revision 1.22  1994/10/12  11:11:19  rbaart
#H  OptimalLinearCode: removed AddZeros
#H
#H  Revision 1.21  1994/10/12  08:08:30  rbaart
#H  OptimalLinearCode: should we check the actual mindist against the table?
#H
#H  Revision 1.20  1994/10/06  09:15:51  rbaart
#H  Added known minimum distance to Extended Golay codes
#H
#H  Revision 1.19  1994/10/06  09:09:28  jcramwin
#H  Added a line in OptimalLinearCode which would fill the upperbound
#H
#H  Revision 1.18  1994/10/06  08:54:52  rbaart
#H  OptimalLinearCode: added codes in library
#H
#H  Revision 1.17  1994/10/05  10:06:23  rbaart
#H  OptimalLinearCode: fixed wordlength of ResidueCode-part
#H
#H  Revision 1.16  1994/10/04  14:24:42  rbaart
#H  OptimalLinearCode: added ResidueCode
#H
#H  Revision 1.15  1994/10/04  11:09:44  jcramwin
#H  some minor bugfixes in OptimalLinearCode
#H
#H  Revision 1.14  1994/10/04  10:31:09  rbaart
#H  OptimalLinearCode: Added warning concerning LB
#H
#H  Revision 1.13  1994/10/04  09:56:54  rbaart
#H  New OptimalLinearCode (seperated calculation for bounds and codes)
#H
#H  Revision 1.12  1994/09/30  17:07:57  rbaart
#H  Added option to return distance to OptimalLinearCode
#H
#H  Revision 1.11  1994/09/30  14:09:32  rbaart
#H  Fixed bug in OptimalLinearCode (UUV construction)
#H
#H  Revision 1.10  1994/09/30  14:01:32  rbaart
#H  *** empty log message ***
#H
#H  Revision 1.9  1994/09/30  13:56:37  rbaart
#H  NullCode added to trivial cases of OptimalLinearCode
#H
#H  Revision 1.8  1994/09/30  13:50:28  rbaart
#H  Added comments to OptimalLinearCode
#H
#H  Revision 1.7  1994/09/30  13:10:35  rbaart
#H  Changed OptimalCode to OptimalLinearCode and adjusted function
#H  for different values of <q>
#H
#H  Revision 1.6  1994/09/30  12:58:32  rbaart
#H  OptimalLinearCode changed for GF(<>2)
#H
#H  Revision 1.5  1994/09/30  11:10:53  rbaart
#H  Added OptimalCode
#H
#H  Revision 1.4  1994/09/29  10:23:55  rbaart
#H  Reverted print messages
#H
#H  Revision 1.3  1994/09/29  10:22:24  rbaart
#H  Checked reading of file codegen.g
#H
#H  Revision 1.2  1994/09/29  10:18:53  rbaart
#H  Minor modifications in layout
#H
#H  Revision 1.1  1994/09/28  10:08:34  rbaart
#H  Initial revision
#H
##

#############################################################################
##
#F  ElementsCode( <L> [, <name> ], <F> )  . . . . . . code from list of words
##
ElementsCode := function (arg)
    local test, L, F, R;
    if not Length(arg) in [2..3] then
        Error("usage: ElementsCode(<E> [, <name> ], <F> )");
    fi;
    F := arg[Length(arg)];
    if IsInt(F) then
        F := GF(F);
    fi;
    L := Codeword( arg[1], F );
    test := WordLength( L[1] );
    if ForAny( L, i -> WordLength( i ) <> test ) then
        Error("All elements must have the same length");
    fi;
    test := Set( L );
    if Length( test ) < Length( L ) then
        L := test;
    fi;
    R := rec(
             isDomain := true,
             isCode := true,
             operations := CodeOps,
             baseField := F,
             elements := L
             );
    if Length(arg) = 3 then
        R.name :=  arg[2];
    else
        R.name :=  "user defined unrestricted code";
    fi;
    return R;
end;

#############################################################################
##
#F  RandomCode( <n>, <M> [, <F>] )  . . . . . . . .  random unrestricted code
##
RandomCode := function(arg)
    local n, M, F, L;
    if Length(arg) = 2 then
        n := arg[1];
        M := arg[2];
        F := GF(2);
    elif Length(arg) = 3 then
        n := arg[1];
        M := arg[2];
        if IsInt(arg[3]) then
            F := GF(arg[3]);
        else
            F := arg[3];
        fi;
    else
        Error("usage: RandomCode(<n>, <M> [, <F>])");
    fi;
    if Size(F)^n < M then
        Error(Size(F),"^",n," < ",M);
    fi;
    L := [];
    while Length(L) < M do
        AddSet(L, List([1..n], i -> Random(F)));
    od;
    return ElementsCode(L, "random unrestricted code", F);
end;

#############################################################################
##
#F  HadamardCode( <H | n> [, <t>] ) . Hadamard code of <t>'th kind, order <n>
##
HadamardCode := function(arg)
    local n, t, A, C;
    if Length(arg) = 1 then
        if IsMat(arg[1]) then
            A := arg[1];
            n := Length(A);
        else
            n := arg[1];
            A := HadamardMat(n);
        fi;
        t := 3;
    elif Length(arg) = 2 then
        if IsInt(arg[1]) then
            n := arg[1];
            A := HadamardMat(n);
        else
            A := arg[1];
            n := Length(A);
        fi;
        t := arg[2];
    else
        Error("usage: HadamardCode( <H | n> [, <t> ] )");
    fi;
    if A * TransposedMat(A) <> n * IdentityMat(n,n) then
        Error("The matrix is not a Hadamard matrix");
    fi;

    A := (A-1)/(-2);
    if t in [1, "1", "A", "a"] then
        A := TransposedMat(TransposedMat(A){[2..n]});
        C := ElementsCode(A, Concatenation("Hadamard code of order ",
                     String(n)), GF(2) );
        C.lowerBoundMinimumDistance := n/2;
        C.upperBoundMinimumDistance := n/2;
        C.innerDistribution := NullVector(n);
        # this was ... := NullVector(n+1);
        # but this seems to be wrong -- Eric Minkes
        C.innerDistribution[1] := 1;
        C.innerDistribution[n/2+1] := Size(C) - 1;
    elif t in [2, "2", "B", "b"] then
        A := TransposedMat(TransposedMat(A){[2..n]});
        Append(A, 1 - A);
        C := ElementsCode(A, Concatenation("Hadamard code of order ",
                     String(n)), GF(2) );
        C.lowerBoundMinimumDistance := n/2 - 1;
        C.upperBoundMinimumDistance := n/2 - 1;
        C.innerDistribution := NullVector(n);
        C.innerDistribution[1] := 1;
        C.innerDistribution[n] := 1;
        # this was ... [n+1]...
        # but this seems to be wrong -- Eric Minkes
        C.innerDistribution[n/2] := Size(C)/2 - 1;
        C.innerDistribution[n/2+1] := Size(C)/2 - 1;
    else
        Append(A, 1 - A);
        C := ElementsCode( A, Concatenation("Hadamard code of order ",
                     String(n)), GF(2) );
        C.lowerBoundMinimumDistance := n/2;
        C.upperBoundMinimumDistance := n/2;
        C.innerDistribution := NullVector(n);
        C.innerDistribution[1] := 1;
        C.innerDistribution[n+1] := 1;
        C.innerDistribution[n/2+1] := Size(C) - 2;
    fi;
    return(C);
end;

#############################################################################
##
#F  ConferenceCode( <n | M> ) . . . . . . . . . . code from conference matrix
##
ConferenceCode := function(arg)
    local n, C, S, I, J, E, zero, w, wd, F, LegendreSym, QRes, els;
    
    LegendreSym := function (i)
        if i = zero then
            return 0;
        elif i in QRes then
            return 1;
        else
            return -1;
        fi;
    end;
       
    if Length(arg) <> 1 then
        Error("usage: ConferenceCode( <n> | <M> )");
    fi;
    if IsInt(arg[1]) then
        n := arg[1];
        if (not IsPrimePowerInt(n)) or (n mod 4 <> 1) then
            Error("n must be a primepower and n mod 4 = 1");
        fi;
        E := Elements(GF(n));
        zero := E[1];
        QRes := [];
        for I in E do
            AddSet(QRes, I^2);
        od;
        S := List(E, i-> List(E, j->LegendreSym(j-i)));
    else
        S := Copy(arg[1]);
        n := Length(S);
        if S*TransposedMat(S) <> (n-1)*IdentityMat(n) or
           TransposedMat(S) <> S then
            Error("argument must be a symmetric conference matrix");
        fi;
        # Normalize S by multiplying rows and columns:
        for I in [2..n] do
            if S[I][1] <> 1 then
                for J in [1..n] do
                    S[I][J] := S[I][J] * -1;
                od;
            fi;
        od;
        for J in [2..n] do
            if S[1][J] <> 1 then
                for I in [1..n] do
                    S[I][J] := S[I][J] * -1;
                od;
            fi;
        od;
        # Strip first row and first column:
        S := List([2..n], i-> S[i]{[2..n]});
        n := n - 1;
    fi;
    F := GF(2);
    els := [NullWord(n, F)];
    I := IdentityMat(n);
    J := NullMat(n,n) + 1;
    Append(els, Codeword(1/2 * (S+I+J), F));
    Append(els, Codeword(1/2 * (-S+I+J), F));
    Add(   els, NullWord(n,F) + 1);
    C := ElementsCode( els, "conference code", F);
    w := WeightCodeword(C.elements[2]);
    wd := List([1..n+1], x -> 0);
    wd[1] := 1; wd[n+1] := 1;
    wd[w+1] := Size(C) - 2;
    C.weightDistribution := wd;
    C.lowerBoundMinimumDistance := (n-1) / 2;
    C.upperBoundMinimumDistance := (n-1) / 2;
    return C;
end;

#############################################################################
##
#F  MOLSCode( [ <n>, ] <q> )  . . . . . . . . . . . . . . . .  code from MOLS
##
##  MOLSCode([n, ] q) returns a (n, q^2, n-1) code over GF(q)
##  by creating n-2 mutually orthogonal latin squares of size q.
##  If n is omitted, a wordlength of 4 will be set.
##  If there are no n-2 MOLS known, the code will return an error
##
MOLSCode := function(arg)
    local q, n, M, i, j, k, C, els;
    if Length(arg) = 1 then
        n := 4;
        q := arg[1];
        M := MOLS(q);
    elif Length(arg) = 2 then
        n := arg[1];
        q := arg[2];
        M := MOLS(q, n-2);
    else
        Error("usage:= MOLSCode( [n, ] q )");
    fi;
    if not IsInt(q) then
        # Argument was a field
        q := Size(q);
    fi;
    if M = false then
        Error("No ",n-2," MOLS of order ",q," are known");
    else
        els:= [];
        for i in [1..q] do
            for j in [1..q] do
                els[(i-1)*q + j]:=[];
                els[(i-1)*q + j][1]:=i-1;
                els[(i-1)*q + j][2]:=j-1;
                for k in [3..n] do
                    els[(i-1)*q + j][k]:=M[k-2][i][j];
                od;
            od;
        od;
        C := ElementsCode( els, Concatenation("code generated by ",
                     String(n-2), " MOLS of order ", String(q)), GF(q) );
        C.lowerBoundMinimumDistance := n - 1;
        C.upperBoundMinimumDistance := n - 1;
        C.weightDistribution:=List( [1..n+1], x -> 0 );
	C.weightDistribution[1]:=1;
        C.weightDistribution[n]:=(q-1) * n;
        C.weightDistribution[n+1]:=(q-1) * (q + 1 - n);
        return C;
    fi;
end;

#############################################################################
##
#F  QuadraticCosetCode( <Q> ) . . . . . . . . . .  coset of RM(1,m) in R(2,m)
##
##  QuadraticCosetCode(Q) returns a coset of the ReedMullerCode of
##  order 1 (R(1,m)) in R(2,m) where m is the size of square matrix Q.
##  Q is the upper triangular matrix that defines the quadratic part of
##  the boolean functions that are in the coset.
##
#QuadraticCosetCode := function(arg)
#    local Q, m, V, R, RM1, k, f, C, h, wd;
#    if Length(arg) = 1 and IsMat(arg[1]) then
#        Q := arg[1];
#    else
#        Error("usage: QuadraticCosetCode( <mat> )");
#    fi;
#    m := Length(Q);
#    V := Tuples(Elements(GF(2)), m);
#    R := V*Q*TransposedMat(V);
#    f := List([1..2^m], i->R[i][i]);
#    RM1 := Concatenation(NullMat(1,2^m,GF(2))+GF(2).one, 
#                   TransposedMat(Tuples(Elements(GF(2)), m)));
#    k := Length(RM1);
#    C := rec(
#             isDomain := true,
#             isCode := true,
#             operations := CodeOps,            
#             baseField := GF(2),
#             wordLength := 2^m,
#    elements := Codeword(List(Tuples(Elements(GF(2)), k) * RM1, t-> t+f)),
#             lowerBoundMinimumDistance := 2^(m-1),
#             upperBoundMinimumDistance := 2^(m-1)
#            );
#    h := RankMat(Q + TransposedMat(Q))/2;
#    wd := NullMat(1, 2^m+1)[1];
#    wd[2^(m-1) - 2^(m-h-1) + 1] := 2^(2*h);
#    wd[2^(m-1) + 1] := 2^(m+1) - 2^(2*h + 1);
#    wd[2^(m-1) + 2^(m-h-1) + 1] := 2^(2*h);
#    C.weightDistribution := wd;
#    C.name := "quadratic coset code";
#    return C;
#end;

#############################################################################
##
#F  GeneratorMatCode( <G> [, <name> ], <F> )  . .  code from generator matrix
##
GeneratorMatCode := function (arg)
    local G, G2, F, C;

    if not Length(arg) in [2..3] then
        Error("usage: GeneratorMatCode( <G> [, <name> ], <F> )");
    fi;
    F := arg[ Length( arg ) ];
    if IsInt( F ) then
        F := GF( F );
    fi;
    if (Length(arg[1]) <> 0) and (IsCodeword(arg[1][1])) then
        arg[1] := VectorCodeword( arg[ 1 ], F );
    fi;
    if (Length(arg[1]) = 0) or (Length(arg[1][1]) = 0 ) then
        Error("use NullCode to generate a code with dimension 0");
    fi;
    G := VectorCodeword( arg[ 1 ], F );
    G2 := BaseMat( G );
    if Length( G2 ) < Length( G ) then
        G := G2;
    fi;
    C := rec(
             generatorMat := G,
             baseField := F,
             wordLength := Length(G[1]),
             isDomain := true,
             isCode := true,
             isLinearCode := true,
             operations := LinCodeOps );
    if Length( arg ) = 3 then
        C.name := arg[2];
    else
        C.name := "code defined by generator matrix";
    fi;
    return C;
end;

#############################################################################
##
#F  CheckMatCode( <H> [, <name> ], <F> )  . . . . . .  code from check matrix
##
CheckMatCode := function (arg)
    local H, H2, F, C;

    if not Length(arg) in [2..3] then
        Error("usage: CheckMatCode( <H> [, <name> ], <F> )");
    fi;
    F := arg[ Length( arg ) ];
    if IsInt( F ) then
        F := GF( F );
    fi;
    if (Length(arg[1]) <> 0) and (IsCodeword( arg[1][1] ) ) then
        arg[1] := VectorCodeword( arg[1], F );
    fi;
    if (Length(arg[1]) = 0) or (Length(arg[1][1]) = 0 ) then
        Error("use WholeSpaceCode to generate a code with redundancy 0");
    fi;
    H := VectorCodeword(arg[1], F);
    H2 := BaseMat(H);
    if Length(H2) < Length(H) then
        H := H2;
    fi;
    C := rec(
             isDomain := true,
             isCode := true,
             isLinearCode := true,
             operations := LinCodeOps,
             checkMat := H,
             wordLength := Length(H[1]),
             baseField := F );
    if Length( arg ) = 3 then
        C.name := arg[2];
    else
        C.name := "code defined by check matrix";
    fi;
    return C;
end;

#############################################################################
##
#F  RandomLinearCode( <n>, <k> [, <F>] )  . . . . . . . .  random linear code
##
RandomLinearCode := function(arg)
    local n, k, F;
    if Length(arg) = 2 then
        n := arg[1];
        k := arg[2];
        F := GF(2);
    elif Length(arg) = 3 then
        n := arg[1];
        k := arg[2];
        if IsInt(arg[3]) then
            F := GF(arg[3]);
        else
            F := arg[3];
        fi;
    else
        Error("usage: RandomLinearCode(<n>, <k> [, <F>])");
    fi;
    if k = 0 then
        return NullCode( n, F );
    else
        return GeneratorMatCode(PermutedCols(List(IdentityMat(k,F), i ->
                   Concatenation(i,List([k+1..n],j->Random(F)))),
                   Random(SymmetricGroup(n))), "random linear code", F);
    fi;
end;

#############################################################################
##
#F  HammingCode( <r> [, <F>] )  . . . . . . . . . . . . . . . .  Hamming code
##
HammingCode := function(arg)
    local r, q, H, F, H2, C, i, j, n, TupAllow, wd;

    TupAllow := function(W)
        local l;
        l := 1;
        while (W[l] = F.zero) and l < Length(W) do
            l := l + 1;
        od;
        return (W[l] = F.one);
    end;

    if Length(arg) = 1 then
        r := arg[1];
        F := GF(2);
    elif Length(arg) = 2 then
        r := arg[1];
        if IsInt(arg[2]) then
            F := GF(arg[2]);
        else
            F := arg[2];
        fi;
    else
        Error("usage: HammingCode(<r> [, <F>])");
    fi;
    q := Size(F);
    if not IsPrimePowerInt(q) then
        Error("q must be prime power");
    fi;
    H := Tuples(Elements(F), r);
    H2 := [];
    j := 1;
    for i in [1..Length(H)] do
        if TupAllow(H[i]) then
            H2[j] := H[i];
            j := j + 1;
        fi;
    od;
    n := (q^r-1)/(q-1);
    C := CheckMatCode(TransposedMat(H2), Concatenation("Hamming (", String(r),
                 ",", String(q), ") code"), F);
    C.lowerBoundMinimumDistance := 3;
    C.upperBoundMinimumDistance := 3;
    C.boundsCoveringRadius := [ 1 ];
    C.isPerfectCode := true;
    C.isSelfDualCode := false;
    C.specialDecoder := HammingDecoder;
    if q = 2 then
        C.isNormalCode := true;
    fi;
    if q = 2 then
        wd := [1, 0];
        for i in [2..n] do
            Add(wd, 1/i * (Binomial(n, i-1) - wd[i] - (n-i+2)*wd[i-1]));
        od;
        C.weightDistribution := wd;
    fi;
    return C;
end;

#############################################################################
##
#F  SimplexCode( <r>, <F> ) .  The SimplexCode is the Dual of the HammingCode
##
SimplexCode := function( arg )
    local C, F;
    if not Length(arg) in [1..2] then
        Error("usage: SimplexCode(<r> [, <F>])");
    fi;
    if Length( arg ) = 2 then
        F := arg[ 2 ];
        if IsInt( F ) then
            F := GF( F );
        fi;
    else
        F := GF(2);
    fi;
    C := DualCode( ApplyFunc( HammingCode, arg ) );
    C.name := "simplex code";
    if F = GF(2) then
        C.isNormalCode := true;
        C.boundsCoveringRadius := [ 2^( arg[ 1 ] - 1 ) - 1 ];
    fi;
    return C;
end;

#############################################################################
##
#F  ReedMullerCode( <r>, <k> )  . . . . . . . . . . . . . .  Reed-Muller code
##
##  ReedMullerCode(r, k) creates a binary Reed-Muller code of dimension k,
##  order r; 0 <= r <= k
##
ReedMullerCode := function (arg)
    local mat,c,src,dest,t,index,num,dim,C,wd, k, r,h,t,A, bcr;

    if Length(arg) = 2 then
        k := arg[2];
        r := arg[1];
    else
        Error("usage: ReedMullerCode(<r>, <k>)");
    fi;
    if r > k then 
        return ReedMullerCode(k, r); #for compatibility with older versions
                                     #of guava, It used to be
                                     #ReedMullerCode(k, r);
    fi;
    mat := [ [] ];
    num := 2^k;
    dim := Sum(List([0..r], x->Binomial(k,x)));
    for t in [1..num] do
        mat[1][t] := Z(2)^0;
    od;
    if r > 0 then
        Append(mat, TransposedMat(Tuples ([0*Z(2), Z(2)^0], k)));
        for t in [2..r] do
            for index in Combinations([1..k], t) do
                dest := List([1..2^k], i->Product(index, j->mat[j+1][i]));
                Append(mat, [dest]);
            od;
        od;
    fi;
    C := GeneratorMatCode( mat, Concatenation("Reed-Muller (", String(r), ",", 
            String(k), ") code"), GF(2) );
    C.lowerBoundMinimumDistance := 2^(k-r);
    C.upperBoundMinimumDistance := 2^(k-r);    
    C.isPerfectCode := false;
    C.isSelfDualCode := (2*r = k-1);
    if r = 0 then
        wd := List([1..num + 1], x -> 0);
        wd[1] := 1;
        wd[num+1] := 1;
        C.weightDistribution := wd;
    elif r = 1 then
        wd := List([1..num + 1], x -> 0);
        wd[1] := 1;
        wd[num + 1] := 1;
        wd[num / 2 + 1] := Size(C) - 2;
        C.weightDistribution := wd;
    elif r = 2 then
        wd := List([1..num + 1], x -> 0);
        wd[1] := 1;
        wd[num + 1] := 1;
        for h in [1..QuoInt(k,2)] do
            A := 2^(h*(h+1));
            for t in [0..2*h-1] do
                A := A*(2^(k-t)-1);
            od;
            for t in [1..h] do
                A := A/(2^(2*t)-1);
            od;
            wd[2^(k-1)+2^(k-1-h)+1] := A;
            wd[2^(k-1)-2^(k-1-h)+1] := A;
        od;
        wd[2^(k-1)+1] := Size(C)-Sum(wd);
        C.weightDistribution := wd;
    fi;
    
    bcr := BoundsCoveringRadius( C );
    
    if 0 <= r and r <= k-3 then
        if IsEvenInt( r ) then
            C.boundsCoveringRadius :=
              [ Maximum( 2^(k-r-3)*(r+4), bcr[1] )
                .. Maximum( bcr ) ];
        else
            C.boundsCoveringRadius := 
              [ Maximum( 2^(k-r-3)*(r+5), bcr[ 1 ] )
                .. Maximum( bcr ) ];
        fi;
    fi;

    if r = k then
        C.boundsCoveringRadius := [ 0 ];
    elif r = k - 1 then
        C.boundsCoveringRadius := [ 1 ];
    elif r = k - 2 then
        C.boundsCoveringRadius := [ 2 ];
    elif r = k - 3 then
        if IsEvenInt( k ) then
            C.boundsCoveringRadius := [ k + 2 ];
        else
            C.boundsCoveringRadius := [ k + 1 ];
        fi;
    elif r = 0 then
        C.boundsCoveringRadius := [ 2^(k-1) ];
    elif r = 1 then
        if IsEvenInt( k ) then
            C.boundsCoveringRadius := [ 2^(k-1) - 2^(k/2-1) ];
        elif k = 5 then
            C.boundsCoveringRadius := [ 12 ];
        elif k = 7 then 
            C.boundsCoveringRadius := [ 56 ];
        elif k >= 15 then
            C.boundsCoveringRadius := [ 2^(k-1) - 2^((k+1)/2)*27/64 
                                        .. 2^(k-1) - 2^( QuoInt(k,2)-1 ) ];
        else
            C.boundsCoveringRadius := [ 2^(k-1) - 2^((k+1)/2)/2
                                        .. 2^(k-1) - 2^( QuoInt(k,2)-1 ) ];
        fi;
    elif r = 2 then
        if k = 6 then
            C.boundsCoveringRadius := [ 18 ];
        elif k = 7 then 
            C.boundsCoveringRadius := [ 40 .. 44 ];
        elif k = 8 then
            C.boundsCoveringRadius := [ 84
              .. Maximum( bcr ) ];
        fi;
    elif r = 3 then
        if k = 7 then
            C.boundsCoveringRadius := [ 20 .. 23 ];
        fi;
    elif r = 4 then
        if k = 8 then
            C.boundsCoveringRadius := [ 22
              .. Maximum( bcr ) ];
        fi;
    fi;

    if r = 1 and
       ( IsEvenInt( k ) or k = 3 or k = 5 or k = 7 ) then
        C.isNormalCode := true;
    fi;
    
    return C;
end;

#############################################################################
##
#F  LexiCode( <M | n>, <d>, <F> )  . . . . .  Greedy code with standard basis
##
LexiCode := function(arg)
    local n, k, d, F, base, elms, i, dist, Sz, vec, word, C, one, zero, pos,
          carry; 

#    AdvanceLexicographical := function(V)
#        # AdvanceLexicographical(V) changes V into its lexicographical
#        # successor. A boolean is returned which is false after the last
#        # vector.
#        local carry, pos;
#        carry := true; pos := Length(V);
#        while (carry) and (pos > 0) do
#            carry := false;
#            if V[pos] = zero then
#                V[pos] := one;
#            else
#                V[pos] := F.root^(LogFFE(V[pos],F.root)+1);
#                if V[pos] = one then 
#                    V[pos] := zero; 
#                    carry := true;
#                fi;
#            fi;
#            pos := pos - 1;
#        od;
#        return not carry;
#    end;

    if Length(arg) <> 3 then
        Error("usage: LexiCode(<matrix | n>, <d>, <F>)");
    fi;
    if IsInt(arg[3]) then
        F := GF( arg[3] );
    else
        F := arg[3];
    fi;
    if IsList(arg[1]) then 
        base := VectorCodeword(arg[1], F);
        n := Length(base[1]);
        k := Length(base);
    else
        n := arg[1]; k := n;
        base := IdentityMat(n, F);
    fi;
    d := arg[2];
    one := F.one;
    zero := F.zero;
    elms:=[ ];
    Sz := 0;
    vec := NullVector(k,F);
    repeat
        word := vec*base;
        i := 1; 
        dist := d; 
        while (dist >= d) and (i <= Sz) do
            dist := DistanceVecFFE(word, elms[i]);
            i := i + 1; 
        od; 
        if dist >= d then 
            Add(elms,Copy(word));
            Sz := Sz + 1;
        fi;
        # generate the (lexicographical) next word in F^k
        carry := true;
        pos := k;
        while carry and (pos > 0) do
            if vec[pos] = zero then
                carry := false;
                vec[pos] := one;
            else
                vec[pos] := F.root^(LogFFE(vec[pos],F.root)+1);
                if vec[pos] = one then 
                    vec[pos] := zero;
                else
                    carry := false;
                fi;
            fi;
            pos := pos - 1;
        od;
    until carry;
    if Size(F) = 2 then  # or even (2^(2^LogInt(LogInt(q,2),2)) = q) ?
        C := GeneratorMatCode(elms, "lexicode", F);
    else
        C := ElementsCode(elms, "lexicode", F);
    fi;
    C.lowerBoundMinimumDistance := d;
    return C;
end;

#############################################################################
##
#F  GreedyCode( <M>, <d> [, <F>] )  . . . . Greedy code from list of elements
##
GreedyCode := function(arg) 
    local n, d, F, space, elms, i, j, dist, Sz, word, C; 
    if not Length(arg) in [2..3] then
    	Error("usage: GreedyCode(<matrix>, <d> [, <F>])");
    elif Length(arg) = 3 then
        F := arg[3];
        space := VectorCodeword(arg[1], F);
    else
        space := VectorCodeword(arg[1]);
        F := DefaultField(space[1]);
    fi;
    if IsInt(F) then
        F := GF(F);
    fi;
    n := Length(space[1]);
    d := arg[2];
    elms := [ space[1] ];
    Sz := 1;
    for word in space do
        i := 1;
        repeat
            dist := DistanceVecFFE(word, elms[i]);
            i := i + 1; 
        until dist < d or i > Sz;
        if dist >= d then 
            Add(elms,word);
            Sz := Sz + 1;
        fi; 
    od;
    C := ElementsCode( elms, "Greedy code, user defined basis", F );
    C.lowerBoundMinimumDistance := d;
    return C;
end;

#############################################################################
##
#F  AlternantCode( <r>, <Y> [, <alpha>], <F> )  . . . . . . .  Alternant code
##
AlternantCode := function(arg)
    local C, n, Y, F, q, i, temp, els, r;
    if Length(arg) < 3 or Length(arg) > 4 then
        Error("usage: AlternantCode(r, Y [, <alpha>], F)");
    fi;
    r := arg[1];
    Y := arg[2];
    n := Length(Y);
    if Length(arg) = 4 then
        F := arg[4];
        els := Set(VectorCodeword(arg[3],F));
        Y := VectorCodeword( Y, F );
    else
            F:=arg[3];
            els := Elements(F){[2..n+1]};
            Y := VectorCodeword( Y, F );
    fi;
    if IsInt(F) then
        F := GF(F);
    fi;
    if ForAny(Y, i-> i = F.zero) then
        Error("Y contains zero");
    elif Length(els) <> Length(Y) then
        Error("<Y> and <alpha> have inequal length or <alpha> is not distinct");
    fi;
    q := F.char;
    temp := NullMat(n, n, F);
    for i in [1..n] do
        temp[i][i] := Y[i];
    od;
    Y := temp;
    C := CheckMatCode( BaseMat(VerticalConversionFieldMat( List([0..r-1], 
                 i -> List([1..n], j-> els[j]^i)) * Y)), "alternant code", F );
    C.lowerBoundMinimumDistance := r + 1;
    return C;
end;

#############################################################################
##
#F  GoppaCode( <G>, <L | n> ) . . . . . . . . . . . . . . . . . .  Goppa code
##
GoppaCode := function(arg)
    local C, GP, F, L, n, q, m, r, zero, temp;

    if Length(arg) <> 2  then
        Error("usage: GoppaCode(<P>, <L | n>  )");
    fi;
    GP := PolyCodeword(arg[1]);
    F := GP.baseRing;
    q := F.char;
    m := F.degree;
    F := GF(q);
    zero := F.zero;
    r := Degree(GP);

    # find m
    if IsInt(arg[2]) then
        n := arg[2];
        m := Maximum(m, LogInt(n,q));
        repeat
            L := Filtered(Elements(GF(q^m)),i -> Value(GP,i) <> zero);
            m := m + 1;
        until Length(L) >= n;
        m := m - 1;
        L := L{[1..n]};
    else
        L := arg[2];
        n := Length(L);
        m := Maximum(m, DefaultField(L).degree);
    fi;
    C := CheckMatCode( BaseMat(VerticalConversionFieldMat( List([0..r-1],
                 i-> List(L, j-> (j)^i / Value(GP, j) )) )), "Goppa code", F);

    # Make the code
    temp := Factors(GP);
    if (q = 2) and (Length(temp) = Length(Set(temp))) then
        C.lowerBoundMinimumDistance := Minimum(n, 2*r + 1);
    else
        C.lowerBoundMinimumDistance := Minimum(n, r + 1);
    fi;
    return C;
end;

#############################################################################
##
#F  CordaroWagnerCode( <n> )  . . . . . . . . . . . . . . Cordaro-Wagner code
##
CordaroWagnerCode := function(arg)
    local r, C, zero, one, F, n, d;
    if Length(arg) = 1 then
        n := arg[1];
        if n < 2 then
            Error("n must be 2 or more");
        fi;
    else 
        Error("usage: CordaroWagnerCode(<n>)");
    fi;
    r := Int((n+1)/3);
    d := (2 * r - Int( (n mod 3) / 2) );
    F := GF(2);
    zero := F.zero;
    one := F.one;
    C := GeneratorMatCode( [Concatenation(List([1..r],i -> zero),
                     List([r+1..n],i -> one)),
                     Concatenation(List([r+1..n],i -> one), List([1..r],
                             i -> zero))], "Cordaro-Wagner code", F );
    C.lowerBoundMinimumDistance := d;
    C.upperBoundMinimumDistance := d;
    C.weightDistribution := List([1..n+1], i-> 0);
    C.weightDistribution[1] := 1;
    C.weightDistribution[2*r+1] := 1;
    C.weightDistribution[n-r+1] := C.weightDistribution[n-r+1] + 2;
    return C;
end;

#############################################################################
##
#F  GeneralizedSrivastavaCode( <a>, <w>, <z> [, <t>] [, <F>] )  . . . . . .  
##
GeneralizedSrivastavaCode := function(arg)
    local C, a, w, z, t, F, n, s, i, H;
    if Length(arg) < 3 or Length(arg) > 5 then
        Error("usage: GeneralizedSrivastavaCode(a, w, z [, t ][, F ])");
    fi;
    a := arg[1];
    w := arg[2];
    z := arg[3];
    if Length(arg) > 3 and IsInt(arg[4]) then
        t := arg[4];
        i := 1;
    else
        t := 1;
        i := 0;
    fi;
    if Length(arg) = 5 or ( Length(arg) = 4 and IsField(arg[3]) ) then
        F := arg[4 + i];
        a := VectorCodeword(a, F);
        w := VectorCodeword(w, F);
        z := VectorCodeword(z, F);
    else
        a := VectorCodeword(a);
        w := VectorCodeword(w);
        z := VectorCodeword(z);
        F := DefaultField(Concatenation(a, w, z));
    fi;
    n := Length(a);
    s := Length(w);
    if IsInt(F) then
        F := GF(F);
    fi;
    if Length(Set(Concatenation(a,w))) <> n + s then
        Error("<alpha> and w are not distinct");
    fi;
    if ForAny(z,i -> i = F.zero) then
        Error("<z> must be nonzero");
    fi;

    H := [];
    for i in List([1..s], index -> List([1..t], vert -> List([1..n],
            hor -> z[hor]/(a[hor] - w[index])^vert))) do
        Append(H, i);
    od;
    C := CheckMatCode( BaseMat(VerticalConversionFieldMat(H)),
                 "generalized Srivastava code", GF(F.char) );
    C.lowerBoundMinimumDistance := s + 1;
    return C;
end;

#############################################################################
##
#F  SrivastavaCode( <a>, <w> [, <mu>] [, <F>] ) . . . . . . . Srivastava code
##
SrivastavaCode := function(arg)
    local C, a, w, F, mu, n, s, i, zero, TheMat;
    if Length(arg) < 2 or Length(arg) > 4 then
        Error("usage: SrivastavaCode(a, w [, <mu>][, F])");
    fi;
    a := arg[1];
    w := arg[2];
    if Length(arg) = 4 then
        F := arg[4];
        a := VectorCodeword(a, F);
        w := VectorCodeword(w, F);
        mu := arg[3];
    elif Length(arg) = 3 then
        if IsInt(arg[3]) then
            mu := arg[3];
            a := VectorCodeword(a);
            w := VectorCodeword(w);
            F := DefaultField(Concatenation(a,w));
        else
            F := arg[3];
            a := VectorCodeword(a, F);
            w := VectorCodeword(w, F);
            mu := 1;
        fi;
    else
        a := VectorCodeword(a);
        w := VectorCodeword(w);
        F := DefaultField(Concatenation(a,w));
        mu := 1;
    fi;
    if IsInt(F) then
        F := GF(F);
    fi;
    n := Length(a);
    s := Length(w);
    if Length(Set(Concatenation(a,w))) <> n + s then
        Error("the elements of <alpha> and w are not distinct");
    fi;
    zero := F.zero;
    for i in [1.. n] do
        if a[i]^mu = zero then
            Error("z[",i,"] = ",a[i],"^",mu," = ",zero);
        fi;
    od;
    TheMat := List([1..s], j -> List([1..n], i -> a[i]^mu/(a[i] - w[j]) ));
    C := CheckMatCode( BaseMat(VerticalConversionFieldMat(TheMat)),
                 "Srivastava code", GF(F.char) );
    C.lowerBoundMinimumDistance := s + 1;
    return C;
end;

#############################################################################
##
#F  ExtendedBinaryGolayCode( )  . . . . . . . . .  extended binary Golay code
##
ExtendedBinaryGolayCode := function()
    local C;
    C := ExtendedCode(BinaryGolayCode());
    C.name := "extended binary Golay code";
    Unbind( C.history );
    C.isCyclicCode := false;
    C.isPerfectCode := false;
    C.isSelfDualCode := true;
    C.boundsCoveringRadius := [ 4 ];
    C.isNormalCode := true;
    C.weightDistribution :=
      [1,0,0,0,0,0,0,0,759,0,0,0,2576,0,0,0,759,0,0,0,0,0,0,0,1];
    #C.automorphismGroup := M24;
    C.lowerBoundMinimumDistance := 8;
    C.upperBoundMinimumDistance := 8;
    return C;
end;

#############################################################################
##
#F  ExtendedTernaryGolayCode( ) . . . . . . . . . extended ternary Golay code
##
ExtendedTernaryGolayCode := function()
    local C;
    C := ExtendedCode(TernaryGolayCode());
    C.isCyclicCode := false;
    C.isPerfectCode := false;
    C.isSelfDualCode := true;
    C.boundsCoveringRadius := [ 3 ];
    C.isNormalCode := true;
    C.name := "extended ternary Golay code";
    Unbind( C.history );
    C.weightDistribution := [1,0,0,0,0,0,264,0,0,440,0,0,24];
    #C.automorphismGroup := M12;
    C.lowerBoundMinimumDistance := 6;
    C.upperBoundMinimumDistance := 6;
    return C;
end;

#############################################################################
##
#F  BestKnownLinearCode( <n>, <k> [, <F>] ) .  returns best known linear code
#F  BestKnownLinearCode( <rec> )
##
##  L describs how to create a code. L is a list with two elements:
##  L[1] is a function and L[2] is a list of arguments for L[1].
##  One or more of the argumenst of L[2] may again be such descriptions and
##  L[2] can be an empty list.
##  The field .construction contains such a list or false if the code is not
##  yet in the apropiatelibrary file (/tbl/codeq.g)
##
BestKnownLinearCode := function(arg)
    local MakeCode, bds, C;

    # L describs how to create a code. L is a list with two elements:
    # L[1] is a function and L[2] is a list of arguments for L[1].
    # One or more of the argumenst of L[2] may again be such descriptions and
    # L[2] can be an empty list.
    MakeCode := function(L)
        #beware: this is the most beautiful function in GUAVA (according to J)
        if IsList(L) and IsBound(L[1]) and IsFunc(L[1]) then
            return ApplyFunc( L[1], List( L[2], i -> MakeCode(i) ) );
        else
            return L;
        fi;
    end;
    
    if not Length(arg) in [1..3] then
        Error("usage: BestKnownLinearCode( < n, k [, F] > | < boundsrec > )");
    fi;

    if IsRec(arg[1]) then
        if IsBound(arg[1].construction) then
            bds := arg[1];
        else
            bds := BoundsMinimumDistance(arg[1].n, arg[1].k,
                            arg[1].q);
        fi;
    else
        bds := ApplyFunc(BoundsMinimumDistance, arg);
    fi;
    if bds.construction = false then 
        Error("code not yet in library");
    else
        C := MakeCode(bds.construction);
        if LowerBoundMinimumDistance(C) > bds.lowerBound then
            Print("New table entry found!\n");
        fi;
        C.lowerBoundMinimumDistance := Maximum(bds.lowerBound,
                                               LowerBoundMinimumDistance(C));
        C.upperBoundMinimumDistance := Minimum(bds.upperBound,
                                               UpperBoundMinimumDistance(C));
        return C;
    fi;
end;

#############################################################################
##
#F  GeneratorPolCode( <G>, <n> [, <name> ], <F> ) .  code from generator poly
##
GeneratorPolCode := function (arg)
    local F, G, n, R;
    if not Length( arg ) in [3..4] then
        Error("usage: GeneratorPolCode(<G>, <n> [, <name> ], <F>)");
    fi;
    n := arg[2];
    F := arg[ Length(arg) ];
    if IsInt(F) then
        F := GF(F);
    fi;
    G := PolyCodeword( arg[1], F );
    G := Gcd(G,F.one*(X(F)^n-1));
    R := rec(
             generatorPol := G,
             baseField := F,
             wordLength := n,
             isDomain := true,
             isCode := true,
             isCyclicCode := true,
             operations := CycCodeOps
             );
    if Length( arg ) = 4 then
        R.name := arg[3];
    else
        R.name := "code defined by generator polynomial";
    fi;
    return R;
end;

#############################################################################
##
#F  CheckPolCode( <H>, <n> [, <name> ], <F> ) . .  code from check polynomial
##
CheckPolCode := function (arg)
    local F, H, n, R;
    if not Length( arg ) in [3..4] then
        Error("usage: CheckPolCode(<H>, <n> [, <name> ], <F>)");
    fi;
    n := arg[2];
    F := arg[ Length(arg) ];
    if IsInt(F) then
        F := GF(F);
    fi;
    H := PolyCodeword( arg[1], F );
    H := Gcd(H, F.one*(X(F)^n-1));
    R := rec(
             isDomain := true,
             isCode := true,
             isCyclicCode := true,
             operations := CycCodeOps,
             checkPol := H,
             baseField := F,
             wordLength := n
             );
    if Length( arg ) = 4 then
        R.name := arg[3];
    else
        R.name := "code defined by check polynomial";
    fi;
    return R;
end;

#############################################################################
##
#F  RepetitionCode( <n> [, <F>] ) . . . . . . . repetition code of length <n>
##
RepetitionCode := function(arg)
    local n, F, C, q;
    if Length(arg)=1 then
        n:=arg[1];
        F:=GF(2);
    elif Length(arg) = 2 then
        n:=arg[1];
        if IsInt(arg[2]) then
            F := GF(arg[2]);
        else
            F:=arg[2];
        fi;
    else
        Error("usage: RepetitionCode(<n> [, <F>])");
    fi;
    q :=Size(F);
    C := GeneratorPolCode(Polynomial(F, List([1..n], t-> F.one)), n,
                 "repetition code", F );
    C.lowerBoundMinimumDistance := n;
    C.upperBoundMinimumDistance := n;
    if n = 2 and q = 2 then
        C.isSelfDualCode := true;
    else
        C.isSelfDualCode := false;
    fi;
    C.weightDistribution := NullVector(n+1);
    C.weightDistribution[1] := 1;
    C.weightDistribution[n+1] := q-1;
    if n < 260 then
        C.automorphismGroup := SymmetricGroup(n);
    fi;
    C.boundsCoveringRadius := [ Minimum(n-1,QuoInt((q-1)*n,q)) ];
    if F = GF(2) then
        C.isNormalCode := true;
    fi;
    if (n mod 2 = 0) or (F <> GF(2)) then
        C.isPerfectCode:=false;
    else
        C.boundsCoveringRadius:= [ QuoInt(n,2) ];
        C.isPerfectCode:=true;
    fi;
    return C;
end;

#############################################################################
##
#F  WholeSpaceCode( <n> [, <F>] ) . . . . . . . . . . returns <F>^<n> as code
##
WholeSpaceCode := function (arg) 
    local C, index, n, F, q;
    if Length(arg) = 1 then
        F := GF(2);
    elif Length(arg) = 2 then
        if IsInt(arg[2]) then
            F := GF(arg[2]);
        else
            F := arg[2];
        fi;
    else
        Error("usage: WholeSpaceCode( <n> [, <F>] )");
    fi;
    n := arg[1];
    C := GeneratorPolCode( F.one*X(F)^0, n, "whole space code", F);
    C.lowerBoundMinimumDistance := 1;
    C.upperBoundMinimumDistance := 1;
    C.automorphismGroup := SymmetricGroup(n);
    C.boundsCoveringRadius:=[ 0 ];
    if F = GF(2) then
        C.isNormalCode := true;
    fi;
    C.isPerfectCode:=true;
    C.isSelfDualCode := false;
    q := Size(F) - 1;
    C.weightDistribution := List([0..n], i-> q^i*Binomial(n, i));
    return C;
end;

#############################################################################
##
#F  CyclicCodes( <n> )  . .  returns a list of all cyclic codes of length <n>
##
CyclicCodes := function(arg)
    local n, F, f, Pl, r, codes;
    if not Length(arg) in [2..3] then
        Error("usage: CyclicCodes( <n> [, <k> ], <F> )");
    fi;
    n := arg[1];
    F := arg[Length(arg)];
    if IsInt(F) then
        F := GF(F);
    fi;
    f := Factors(F.one*(X(F)^n-1));
    if Length(arg) = 2 then
        Pl := List(Combinations(f), c->Product(c)*X(F)^0);
    else
        Pl := [];
        r := n - arg[2];

        codes := function(f, g)
           local i, tempf;
           if Degree(g) < r then
               i := 1;
               while i <= Length(f) and Degree(g)+Degree(f[i][1]) <= r do
                   if f[i][2] = 1 then
                       tempf := f{[ i+1 .. Length(f) ]};
                   else
                       tempf := Copy( f );
                       tempf[i][2] := f[i][2] - 1;
                   fi;
                   codes( tempf, g * f[i][1] );
                   i := i + 1;
               od;
           elif Degree(g) = r then
               Add( Pl, g );
           fi;
        end;

        codes( Collected( f ), X(F)^0 );
    fi;
    return List(Pl, p->GeneratorPolCode(p,n,"enumerated code",F));
end;

#############################################################################
##
#F  NrCyclicCodes( <n>, <F>)  . . .  number of cyclic codes of length <n>
##
NrCyclicCodes := function(arg)
    local n, F;
    if Length(arg) = 2 then
        n := arg[1];
        if IsInt(arg[2]) then
            F := GF(arg[2]);
        else
            F := arg[2];
        fi;
    else
        Error("usage: NrCyclicCodes( <n>, <F> )");
    fi;
    return NrCombinations(Factors(F.one*(X(F)^n-1)));
end;

#############################################################################
##
#F  BCHCode( <n> [, <b>], <delta> [, <F>] ) . . . . . . . . . . . .  BCH code
##
##  BCHCode (n [, b ], delta [, F]) returns the BCH code over F with
##  wordlength n, designedDistance delta, constructed from powers
##  x^b, x^(b+1), ..., x^(b+delta-2), where x is a primitive n'th power root
##  of unity; b = 1 by default; the function returns a narrow sense BCH code
##  Gcd(n,q) = 1 and 2<=delta<=n-b+1
BCHCode := function (arg)
    local n, q, start, stop, m, b, C, test, Cyclo, PowerSet, t,
          zero, desdist, G, superfl, i, BCHTable;

    BCHTable := [ [31,11,11], [63,36,11], [63,30,13], [127,92,11],
                  [127,85,13], [255,223,9], [255,215,11], [255,207,13],
                  [255,187,19], [255,171,23], [255,155,27], [255,99,47],
                  [255,79,55], [255,29,95], [255,21,111] ];

    if Length(arg) = 2 then
        n := arg[1];
        q := 2;
        start := 1;
        stop := arg[2] - 1;
    elif Length(arg) = 3 then
        n := arg[1];
        if IsInt(arg[3]) then
            start := arg[2];
            stop := start + arg[3] - 2;
            q := 2;
        else
            start := 1;
            stop := arg[2] - 1;
            if IsInt(arg[3]) then
                q := arg[3];
            else
                q := Size(arg[3]);
            fi;
        fi;
    elif Length(arg) = 4 then
        n := arg[1];
        start := arg[2];
        stop := start + arg[3] - 2;
        if IsInt(arg[4]) then
            q := arg[4];
        else
            q := Size(arg[4]);
        fi;
    else
        Error ("usage: BCHCode (<n> [, <start>], <delta> [, <F>])");
    fi;
    if Gcd(n,q) <> 1 then
        Error ("n and q must be relative primes");
    fi;
    zero := GF(q).zero;
    m := OrderMod(q,n);
    b := PrimitiveUnityRoot(q,n);
    PowerSet := [start..stop];
    G := X(GF(q))^0;
    while Length(PowerSet) > 0 do
        test := PowerSet[1];
        G := G * Polynomial(GF(q), MinPol(b^test));
        t := (q*test) mod n;
        while t <> test do
            RemoveSet(PowerSet, t);
            t := (q*t) mod n;
        od;
        RemoveSet(PowerSet, test);
    od;
    C := GeneratorPolCode(G, n, GF(q));
    C.specialDecoder := BCHDecoder;
    # Calculate Bose distance:
    Cyclo := CyclotomicCosets(q,n);
    PowerSet := [];
    for t in [start..stop] do
        for test in [1..Length(Cyclo)] do
            if t in Cyclo[test] then
                AddSet(PowerSet, Cyclo[test]);
            fi;
        od;
    od;
    PowerSet := Flat(PowerSet);
    while stop + 1 in PowerSet do
    	stop := stop + 1;
    od;
    while start - 1 in PowerSet do
    	start := start - 1;
    od;
    desdist := stop - start + 2;
    if desdist > n then
        Error("invalid designed distance");
    fi;
    # In some cases the true minimumdistance is known:
    C.designedDistance := desdist;
    C.lowerBoundMinimumDistance := desdist;
    if (q=2) and (n mod desdist = 0) and (start = 1) then
        C.upperBoundMinimumDistance := desdist;
    elif q=2 and desdist mod 2 = 0 and (n=2^m - 1) and (start=1) and
      (Sum(List([0..QuoInt(desdist-1, 2) + 1], i -> Binomial(n, i))) >
       (n + 1) ^ QuoInt(desdist-1, 2)) then
        C.upperBoundMinimumDistance := desdist;
    elif (n = q^m - 1) and (desdist = q^OrderMod(q,desdist) - 1)
      and (start=1) then
        C.upperBoundMinimumDistance := desdist;
    fi;
    # Look up this code in the table
    if start=1 then
        for i in BCHTable do
            if i[1] = n and i[2] = Dimension(C) then
                C.lowerBoundMinimumDistance := i[3];
                C.upperBoundMinimumDistance := i[3];
            fi;
        od;
    fi;
    # Calculate minimum of q*desdist - 1 for primitive n.s. BCH code
    if q^m - 1 = n and start = 1 then
        PowerSet := [start..stop];
        superfl := true;
        i := PowerSet[Length(PowerSet)] * q mod n;
        while superfl do
            while i <> PowerSet[Length(PowerSet)] and not i in PowerSet do
                i := i * q mod n;
            od;
            if i = PowerSet[Length(PowerSet)] then
                superfl := false;
            else
                PowerSet := PowerSet{[1..Length(PowerSet)-1]};
                i := PowerSet[Length(PowerSet)] * q mod n;
            fi;
        od;
        C.upperBoundMinimumDistance := Minimum(UpperBoundMinimumDistance(C),
                                               q * (Length(PowerSet) + 1) - 1);
    fi;
    C.name := Concatenation("BCH code, delta=",
                      String(desdist), ", b=", String(start));
    return C;
end;

#############################################################################
##
#F  ReedSolomonCode( <n>, <d> ) . . . . . . . . . . . . . . Reed-Solomon code
##
##  ReedSolomonCode (n, d) returns a primitive narrow sense BCH code with
##  wordlength n, over alphabet q = n+1, designed distance d
ReedSolomonCode := function (arg)
    local C,b,q,wd,w, n, d;
    if Length(arg) = 2 then
    	n := arg[1];
        d := arg[2];
    else
        Error("usage: ReedSolomonCode(<n>, <d>)");
    fi;
    q := n+1;
    if not IsPrimePowerInt(q) then
        Error("q = n+1 must be a prime power");
    fi;
    b := Z(q);
    
    C := GeneratorPolCode(Product([1..d-1], i-> (X(GF(q))-b^i)), n,
                 "Reed-Solomon code", GF(q) );
    C.roots := List([1..d-1], i->b^i);
    C.lowerBoundMinimumDistance := d;
    C.upperBoundMinimumDistance := d;
    C.designedDistance := d;
    C.specialDecoder := BCHDecoder;
    IsMDSCode(C);	    # Calculate weightDistribution field
    return C;
end;

#############################################################################
##
#F  RootsCode( <n>, <list> )  . . . code constructed from roots of polynomial
##
##  RootsCode (n, rootlist) or RootsCode (n, <powerlist>, F) returns the
##  code with generator polynomial equal to the least common multiplier of
##  the minimal polynomials of the n'th roots of unity in the list.
##  The code has wordlength n
##
RootsCode := function(arg)
    local n, L, F, G, C, num, power, q, z, i, rootslist, powerlist, max, rs;
    if Length(arg) = 2 then
        n := arg[1];
        L := Set(arg[2]);
        q := Field(L).char;
        z := Field(L).root;
        F := GF(q);
        if List(L, i->i^n) <> NullVector(Length(L), F) + z^0 then
            Error("powers must all be n'th roots of unity");
        fi;
    elif Length(arg) = 3 then
        n := arg[1];
        if IsInt(arg[3]) then
            q := arg[3];
        else
            q := Size(arg[3]);
        fi;
        F := GF(q);
        z := PrimitiveUnityRoot(q, n);
        L := Set(List(arg[2], i->z^i));
    else
        Error("usage: RootsCode( <n>, <rootlist> )");
    fi;
    rs := Copy(L);
    G := X(F)^0;
    rootslist := [];
    powerlist := [];
    max := LogFFE(PrimitiveUnityRoot(q,n));
    while rs <> [] do
        power := rs[1];
        i := power;
        repeat
            AddSet(rootslist, i);
            i := i^q;
            RemoveSet(rs, i);
            AddSet(powerlist, LogFFE(i)/max);
        until i = power;
        G := G * Polynomial(F, MinPol(power));
    od;
    C := GeneratorPolCode( G, n, "code defined by roots", F );
    # Find the largest number of successive powers for BCH bound
    max := 1;
    i := 1;
    num := Length(powerlist);
    for z in [2..num] do
        if powerlist[z] <> powerlist[i] + z-i then
            max := Maximum(max, z - i);
            i := z;
        fi;
    od;
    C.lowerBoundMinimumDistance := Maximum(max, num+1 - i) + 1;
    C.roots := rootslist;
    return C;
end;

#############################################################################
##
#F  QRCode( <n> [, <F>] ) . . . . . . . . . . . . . .  quadratic residue code
##
QRCode := function (arg)
    local n, q, m, b, Q, N, t, g, lower, upper, C, F;
    if Length(arg) = 1 then
        n := arg[1];
        q := 2;
    elif Length(arg) = 2 then
        n := arg[1];
        if IsInt(arg[2]) then
            q := arg[2];
        else
            q := Size(arg[2]);
        fi;
    else
        Error("usage: QRCode( n [, F] )");
    fi;
    if Jacobi(q,n) <> 1 then
        Error("q must be a quadratic residue modulo n");
    elif not IsPrimeInt(n) then
        Error("n must be a prime");
    elif not IsPrimeInt(q) then
        Error("q must be a prime");
    fi;
    m := OrderMod(q,n);
    F := GF(q^m);
    b := PrimitiveUnityRoot(q,n);
    Q := [];
    N := [1..n];
    for t in [1..n-1] do
        AddSet(Q, t^2 mod n);
    od;
    for t in Q do
        RemoveSet(N, t);
    od;
    g := Product(Q, i -> Polynomial(F, [-b^i, b^0]));
    C := GeneratorPolCode( Polynomial(GF(q), g.coefficients), n,
                 "quadratic residue code", GF(q) );
    if RootInt(n)^2 = n then
        lower := RootInt(n);
    else
        lower := RootInt(n)+1;
    fi;
    if n mod 4 = 3 then
        while lower^2-lower+1 < n do
            lower := lower + 1;
        od;
    fi;
    if (n mod 8 = 7) and (q = 2) then
        while lower mod 4 <> 3 do
            lower := lower + 1;
        od;
    fi;
    upper := WeightCodeword(Codeword(g.coefficients));
    C.lowerBoundMinimumDistance := lower;
    C.upperBoundMinimumDistance := upper;
    return C;
end;

#############################################################################
##
#F  NullCode( <n> [, <F>] ) . . . . . . . . . . .  code consiting only of <0>
##
NullCode := function(arg)
    local n, F, C;
    if Length(arg) = 1 then
        n := arg[1];
        F := GF(2);
    elif Length(arg) = 2 then
        n := arg[1];
        if IsInt(arg[2]) then
            F := GF(arg[2]);
        else
            F := arg[2];
        fi;
    else 
        Error("usage: NullCode( <n> [, <F> ] )");
    fi;
    C := GeneratorPolCode(0*X(F), n, "nullcode", F);
    C.lowerBoundMinimumDistance := n;
    C.upperBoundMinimumDistance := n;
    C.weightDistribution := Concatenation([1], NullVector(n));
    C.automorphismGroup := SymmetricGroup(n);
    C.boundsCoveringRadius := [ n ];
    return C;
end;

#############################################################################
##
#F  FireCode( <G>, <b> )  . . . . . . . . . . . . . . . . . . . . . Fire code
##
##  FireCode (G, b) constructs the Fire code that is capable of correcting any
##  single error burst of length b or less.
##  G is a primitive polynomial of degree m
##
FireCode := function(arg)
    local G, b, m, GenPol, F, C, n;
    if Length(arg) <> 2 then
        Error("usage: FireCode( <Pol>, <Int> )");
    else
        G := PolyCodeword(arg[1]);
        b := arg[2];
    fi;
    if G.baseRing <> GF(2) then
        Error("polynomial must be over GF(2)");
    fi;
    if Length(Factors(G)) <> 1 then
        Error("polynomial G must be primitive");
    fi;
    m := Degree(G);
    n := Lcm(2^m-1,2*b-1);
    C := GeneratorPolCode( G*(X(GF(2))^(2*b-1) + 1), n,
                 Concatenation(String(b), " burst error correcting fire code"),
                 GF(2) );
    return C;
end;

#############################################################################
##
#F  BinaryGolayCode( )  . . . . . . . . . . . . . . . . . . binary Golay code
##
BinaryGolayCode := function()
    return rec(
             generatorPol := Polynomial(GF(2), 
                     Z(2)^0*[1,0,1,0,1,1,1,0,0,0,1,1]),
             baseField := GF(2),
             wordLength := 23,
             dimension := 12,
             redundancy := 11,
             size := 2^12,
             name := "binary Golay code",
             lowerBoundMinimumDistance := 7,
             upperBoundMinimumDistance := 7,
             weightDistribution :=
             [1,0,0,0,0,0,0,253,506,0,0,1288,1288,0,0,506,253,0,0,0,0,0,0,1],
             boundsCoveringRadius := [ 3 ],
             isNormalCode := true,  
             isDomain := true,
             isCode := true,
             isLinearCode := true,
             isCyclicCode := true,
             isPerfectCode := true,
             operations := CycCodeOps
             );
end;

#############################################################################
##
#F  TernaryGolayCode( ) . . . . . . . . . . . . . . . . .  ternary Golay code
##
TernaryGolayCode := function()
    return rec(				
                 generatorPol := Polynomial(GF(3),
                         Z(3)^0*[2,0,1,2,1,1]),
                 baseField := GF(3),
                 wordLength := 11,
                 dimension := 6,
                 redundancy := 5,
                 size := 3^6,
                 name := "ternary Golay code",
                 lowerBoundMinimumDistance := 5,
                 upperBoundMinimumDistance := 5,
                 weightDistribution := [1,0,0,0,0,132,132,0,330,110,0,24],
                 boundsCoveringRadius := [ 2 ],
                 isNormalCode := true,  
                 isDomain := true,
                 isCode := true,
                 isLinearCode := true,
                 isCyclicCode := true,
                 isPerfectCode := true,
                 operations := CycCodeOps
                 );
end;

