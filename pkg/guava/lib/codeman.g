#############################################################################
##
#A  codeman.g               GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains functions for manipulating codes
##
#H  $Log: codeman.g,v $
#H  Revision 1.2  1997/01/20 15:05:30  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.55  1996/05/15  09:29:55  eminkes
#H  Removed a bug I introduced in UUVCode.
#H
#H  Revision 1.54  1996/04/26  09:47:42  eminkes
#H  Made changes that were necessary for the new release of GUAVA.
#H
#H  Revision 1.53  1995/05/18  02:12:37  mschoene
#H  changed location of executables, fixed 'TmpName' call
#H
#H  Revision 1.52  1995/04/13  12:23:17  jcramwin
#H  fixed a bug in ShortenedCode
#H
#H  Revision 1.51  1995/02/14  12:37:30  jcramwin
#H  changed some C.wordLength in WordLength(C)
#H
#H  Revision 1.50  1994/11/15  10:06:17  rbaart
#H  LinearDirectSumCode: fixed C.wordLength
#H
#H  Revision 1.49  1994/11/15  08:38:36  rbaart
#H  DirectSumCode: if both codes are self orthogonal, the result is also.
#H
#H  Revision 1.48  1994/11/11  16:50:30  rbaart
#H  Changed name of ConversionFieldCode
#H
#H  Revision 1.47  1994/11/10  16:02:29  rbaart
#H  Changed names
#H
#H  Revision 1.46  1994/11/10  12:21:41  rbaart
#H  Changed some names of manipulations
#H
#H  Revision 1.45  1994/11/10  09:25:08  jcramwin
#H  fixed a bug with PuncteredCode(WholeSpacecode(4,2)) for upperboundMinDist
#H
#H  Revision 1.44  1994/11/09  14:06:49  jcramwin
#H  changed the .name of shortened
#H
#H  Revision 1.43  1994/11/09  10:37:58  rbaart
#H  EvenWeightSubcode: fixed bug with lowerBoundMinimumDistance
#H
#H  Revision 1.42  1994/11/09  10:27:04  jcramwin
#H  changed something yesterday
#H
#H  Revision 1.41  1994/11/07  15:50:07  rbaart
#H  Fixed bug: changed Cnew := Concat... into Cnew.name := Concat...
#H
#H  Revision 1.40  1994/11/07  11:24:27  jcramwin
#H  changed "u | u + v Construction" in "u|u+v construction"
#H
#H  Revision 1.39  1994/11/07  10:18:46  jcramwin
#H  changed "U | U + V" in "u | u + v"
#H
#H  Revision 1.38  1994/11/03  14:43:55  rbaart
#H  PuncturedCode: fixed bug in C.name
#H
#H  Revision 1.37  1994/11/03  10:58:17  jcramwin
#H  changed NullVector
#H
#H  Revision 1.36  1994/11/02  15:59:37  jcramwin
#H  sharpened the upperbounds for puncturing cyclic codes
#H
#H  Revision 1.35  1994/10/31  15:28:47  jcramwin
#H  Changed DescriptionForCode changed in HistoryOfCode
#H
#H  Revision 1.34  1994/10/27  14:05:50  rbaart
#H  CycCodeOps.AugmentEdCode
#H
#H  Revision 1.33  1994/10/27  08:19:01  jcramwin
#H  changed the dispatchers for the functions which manipulate two codes
#H  (according to the decisions in the meeting)
#H
#H  Revision 1.32  1994/10/24  16:00:53  rbaart
#H  Fixed small bugs in displaying
#H
#H  Revision 1.31  1994/10/24  09:58:52  rbaart
#H  PuncturedCode: fixed bug in DescriptionForCode
#H
#H  Revision 1.30  1994/10/24  09:19:19  rbaart
#H  PermutedCode now copies fields that stay the same in stead of omitting
#H  fields that change
#H
#H  Revision 1.29  1994/10/21  15:12:45  rbaart
#H  Changed .name in DescriptionForCode
#H
#H  Revision 1.28  1994/10/21  14:28:10  jcramwin
#H  did some DescriptionForCode
#H
#H  Revision 1.27  1994/10/21  11:57:07  jcramwin
#H  still working on the DescriptionForCode
#H
#H  Revision 1.26  1994/10/21  09:34:26  rbaart
#H  PermutedCode unbinds automorphism group
#H
#H  Revision 1.25  1994/10/20  18:01:25  jcramwin
#H  changed some functions with DescriptionForCode
#H
#H  Revision 1.24  1994/10/20  12:31:52  rbaart
#H  LengthenedCode: fixed bug IsCyclicCode(C)
#H
#H  Revision 1.23  1994/10/20  11:14:40  jcramwin
#H  ExtendecCode and LengthenedCode can now be called with 2 arguments
#H
#H  Revision 1.22  1994/10/20  09:12:01  jcramwin
#H  ExtendedCode and LengthenedCode are now ready to get an extra argument
#H
#H  Revision 1.21  1994/10/19  15:14:22  rbaart
#H  ExtendedCode: fixed bug with .wordLength
#H
#H  Revision 1.20  1994/10/13  15:08:49  rbaart
#H  Changed codeword functions
#H
#H  Revision 1.19  1994/10/12  12:40:54  rbaart
#H  ConstantWeightSubcode: check for wt=0
#H
#H  Revision 1.18  1994/10/11  13:39:19  rbaart
#H  Added ConstructionBCode
#H
#H  Revision 1.17  1994/10/10  08:37:12  rbaart
#H  EvenWeightSubcode: fixed problem with upperbound mindist
#H
#H  Revision 1.16  1994/10/06  14:52:57  rbaart
#H  Very small changes
#H
#H  Revision 1.15  1994/10/06  13:58:42  rbaart
#H  DualCode no longer uses KrawtchoukMat
#H
#H  Revision 1.14  1994/10/06  10:08:43  rbaart
#H  DualCode: fixed bug with weightDistribution
#H
#H  Revision 1.13  1994/10/05  09:38:48  rbaart
#H  ResidueCode: fixed minor bugs
#H
#H  Revision 1.12  1994/10/04  14:46:33  rbaart
#H  ConverionFieldCode: fixed bug w.r.t. Field
#H
#H  Revision 1.11  1994/10/04  14:14:53  rbaart
#H  Changed ResidueCode
#H
#H  Revision 1.10  1994/10/04  13:35:56  rbaart
#H  Added ResidueCode
#H
#H  Revision 1.9  1994/10/04  08:54:53  rbaart
#H  Added some used names in local
#H
#H  Revision 1.8  1994/10/04  08:32:35  jcramwin
#H  removed some unused names in local
#H
#H  Revision 1.7  1994/09/30  10:36:59  rbaart
#H  Fixed bug in DualCode (name of cyclic code)
#H
#H  Revision 1.6  1994/09/29  11:31:45  rbaart
#H  Added ConcatenationCode
#H
#H  Revision 1.5  1994/09/29  10:32:16  rbaart
#H  Fixed Ops.operations. ... bug
#H
#H  Revision 1.4  1994/09/28  16:46:14  jcramwin
#H  made DirectSumCode a dispatcher-function
#H
#H  Revision 1.3  1994/09/28  15:54:35  jcramwin
#H  made dischpachter functions of UnionCode, IntersectionCode and
#H  DirectProductCode
#H
#H  Revision 1.2  1994/09/28  11:47:45  rbaart
#H  Split ConversionFieldCode & CosetCode
#H
#H  Revision 1.1  1994/09/28  10:30:51  rbaart
#H  Initial revision
#H
##

#############################################################################
##
#F  DualCode( <C> ) . . . . . . . . . . . . . . . . . . . .  dual code of <C>
##
DualCode := function(C)
    if IsBound(C.isSelfDualCode) and C.isSelfDualCode then
        return ShallowCopy(C);
    else
        IsCyclicCode(C);
        return C.operations.DualCode(C);
    fi;
end;

CodeOps.DualCode := function(C)
    local Cnew;
    Cnew := CheckMatCode(BaseMat(VectorCodeword(Elements(C))), "dual code",
                    Field(C));
    Cnew.history := History(C);
    return Cnew;
end;

LinCodeOps.DualCode := function(C)
    local C1, Pr, n, newwd, wd, oldrow, newrow, i, j, q;
    if IsBound(C.generatorMat) then
        C1 := CheckMatCode(C.generatorMat, "dual code", Field(C));
    elif IsBound(C.checkMat) then
        C1 := GeneratorMatCode(C.checkMat, "dual code", Field(C));
    fi;
    if IsBound(C.weightDistribution) then
        n := WordLength(C);
        wd := C.weightDistribution;
        q := Size(Field(C)) - 1;
        newwd := [Sum(wd)];
        oldrow := List([1..n+1], i->1);
        newrow := [];
        for i in [2..n+1] do
            newrow[1] := Binomial(n, i-1) * q^(i-1);
            for j in [2..n+1] do
                newrow[j] := newrow[j-1] - q * oldrow[j] - oldrow[j-1];
            od;
            newwd[i] := newrow * wd;
            oldrow := Copy(newrow);
        od;
        C1.weightDistribution := newwd / ((q+1) ^ Dimension(C));
        Pr := PositionProperty(C1.weightDistribution{[2..n+1]}, i-> i <> 0);
        if Pr = false then
            Pr := n;
        fi;
        C1.lowerBoundMinimumDistance := Pr;
        C1.upperBoundMinimumDistance := Pr;
    fi;
    C1.history := History(C);
    return C1;
end;

CycCodeOps.DualCode := function(C)
    local C1, r, n, Pr, wd, q, newwd, oldrow, newrow, i, j;
    if IsBound(C.generatorPol) then
        r := ReciprocalPolynomial(C.generatorPol,Redundancy(C));
        r := r/LeadingCoefficient(r);
        C1 := CheckPolCode(r, WordLength(C), "dual code", Field(C));
    fi;
    if IsBound(C.checkPol) then
        r := ReciprocalPolynomial(C.checkPol,Dimension(C));
        r := r/LeadingCoefficient(r);
        C1 := GeneratorPolCode(r, WordLength(C), "dual code", Field(C));
    fi;
    if IsBound(C.weightDistribution) then
        n := WordLength(C);
        wd := C.weightDistribution;
        q := Size(Field(C)) - 1;
        newwd := [Sum(wd)];
        oldrow := List([1..n+1], i->1);
        newrow := [];
        for i in [2..n+1] do
            newrow[1] := Binomial(n, i-1) * q^(i-1);
            for j in [2..n+1] do
                newrow[j] := newrow[j-1] - q * oldrow[j] - oldrow[j-1];
            od;
            newwd[i] := newrow * wd;
            oldrow := Copy(newrow);
        od;
        C1.weightDistribution := newwd / ((q+1) ^ Dimension(C));
        Pr := PositionProperty(C1.weightDistribution{[2..n+1]}, i-> i <> 0);
        if Pr = false then
            Pr := n;
        fi;
        C1.lowerBoundMinimumDistance := Pr;
        C1.upperBoundMinimumDistance := Pr;
    fi;
    C1.history := History(C);
    return C1;
end;

#############################################################################
##
#F  AugmentedCode( <C> [, <L>] )  . . .  add words to generator matrix of <C>
##
AugmentedCode := function(arg)
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: AugmentedCode(<C> [, <list>])");
    fi;
    if not IsLinearCode(arg[1]) then
        Error("argument must be a linear code");
    fi;
    return ApplyFunc( arg[1].operations.AugmentedCode, arg );
end;

LinCodeOps.AugmentedCode := function (arg)
    local Cold, C, L;
    if Length(arg) = 1 then
        Cold := arg[1];
        L := NullMat(1,WordLength(Cold), Field(Cold)) + Field(Cold).one;
    else
        Cold := arg[1];
        L := VectorCodeword(arg[2], Cold );
        if not IsList(L[1]) then
            L := [L];
        else
            L := Set(L);
        fi;
    fi;
    C := GeneratorMatCode(BaseMat(Concatenation(GeneratorMat(Cold),L)),
                 Concatenation("code, augmented with ", String(Length(L)),
                         " word(s)"), Field(Cold));
    if Length(C.generatorMat) > Dimension(Cold) then
        C.upperBoundMinimumDistance := Minimum(
                 UpperBoundMinimumDistance(Cold),
                 Minimum(List(L, l-> WeightCodeword(Codeword(l)))));
        C.history := History(Cold);
        return C;
    else
        return ShallowCopy(Cold);
    fi;
end;

CycCodeOps.AugmentedCode := LinCodeOps.AugmentedCode;

#############################################################################
##
#F  EvenWeightSubcode( <C> )  . . .  code of all even-weight codewords of <C>
##
EvenWeightSubcode := function(C)
    IsCyclicCode(C);
    return C.operations.EvenWeightSubcode(C);
end;

CodeOps.EvenWeightSubcode := function(Cold)
    local C, P, n, G, Els, E, d, i, s,q, Gold;
    q := Size(Field(Cold));
    n := WordLength(Cold);
    Els := Cold.operations.Elements(Cold);
    E := []; s := 0;
    for i in [1..Size(Cold)] do
        if IsEvenInt(WeightCodeword(Els[i])) then
            Append(E, [Els[i]]);
            s := s + 1;
        fi;
    od;
    if s <> Size(Cold) then
        C := ElementsCode( E, "even weight subcode", GF(q) );
        d := [LowerBoundMinimumDistance(Cold),
              UpperBoundMinimumDistance(Cold)];
        for i in [1..2] do 
            if q=2 and IsOddInt(d[i] mod 2) then
                d[i] := Minimum(d[i]+1,n);
            fi;
        od;
        C.lowerBoundMinimumDistance := d[1];
        C.upperBoundMinimumDistance := d[2];
        if IsBound(Cold.weightDistribution) then
            C.weightDistribution := Copy(Cold.weightDistribution);
            for i in [1..QuoInt(n+1,2)] do 
                C.weightDistribution[2*i] := 0;
            od;
        fi;
        C.history := History(Cold);
        return C;
    else
        return ShallowCopy(Cold);
    fi;
end;

LinCodeOps.EvenWeightSubcode := function(Cold)
    local C, P, edited, n, G, Els, E, i, s,q, Gold;
    q := Size(Field(Cold));
    n := WordLength(Cold);
    edited := false;
    if q = 2 then
        # Why is the next line needed?
        P := NullVector(n, GF(2));
        G := [];
        Gold := GeneratorMat(Cold);
        for i in [1..Dimension(Cold)] do
            if WeightCodeword(Codeword(Gold[i])) mod 2 <> 0 then
                if not edited then
                    P := Gold[i];
                    edited := true;
                else
                    Append(G, [Gold[i]+P]);
                fi;
            else
                Append(G, [Gold[i]]);
            fi;
        od;
        if edited then
            C := GeneratorMatCode(BaseMat(G),"even weight subcode",GF(q));
        fi;
    else
        Els := Cold.operations.Elements(Cold);
        E := []; s := 0;
        for i in [1..Size(Cold)] do
            if IsEvenInt(WeightCodeword(Els[i])) then
                Append(E, [Els[i]]);
                s := s + 1;
            fi;
        od;
        edited := (s <> Size(Cold));
        if edited then
            C := ElementsCode(E, "even weight subcode", GF(q) );
        fi;
    fi;

    if edited then
        C.lowerBoundMinimumDistance := Minimum(n, 
                                             LowerBoundMinimumDistance(Cold));
        if q = 2 and IsOddInt(C.lowerBoundMinimumDistance) then
            C.lowerBoundMinimumDistance := C.lowerBoundMinimumDistance + 1;
        fi;
        if IsBound(Cold.weightDistribution) then
            C.weightDistribution := Copy(Cold.weightDistribution);
            for i in [1..QuoInt(n+1,2)] do 
                C.weightDistribution[2*i] := 0;
            od;
        fi;
        C.history := History(Cold);
        return C;
    else
        return ShallowCopy(Cold);
    fi;
end;

CycCodeOps.EvenWeightSubcode := function(Cold)
    local C, P, edited, n, Els, E, d, i, q;
    q := Size(Field(Cold));
    n := WordLength(Cold);
    edited := false;
    if (q =2) then
        P := X(GF(2))-GF(2).one;
        if Gcd(P, GeneratorPol(Cold)) <> P then
            C := GeneratorPolCode( GeneratorPol(Cold)*P, n,
                         "even weight subcode", Field(Cold));
            if IsBound(C.roots) then
                AddSet(C.roots, Z(2)^0);
            fi;
            edited := true;
        fi;
    else
        Els := Cold.operations.Elements(Cold);
        E := [];
        for i in [1..Size(Cold)] do
            if IsEvenInt(WeightCodeword(Els[i])) then
                Append(E, [Els[i]]);
            else
                edited := true;
            fi;
        od;
        if edited then
            C := ElementsCode(E, "even weight subcode", Field(Cold));
        fi;
    fi;

    if edited then
        C.lowerBoundMinimumDistance := Minimum(n,
                                            LowerBoundMinimumDistance(Cold));
        if q = 2 and IsOddInt(C.lowerBoundMinimumDistance) then
            C.lowerBoundMinimumDistance := C.lowerBoundMinimumDistance + 1;
        fi;
        if IsBound(Cold.weightDistribution) then
            C.weightDistribution := Copy(Cold.weightDistribution);
            for i in [1..QuoInt(n+1,2)] do 
                C.weightDistribution[2*i] := 0;
            od;
        fi;
        C.history := History(Cold);
        return C;
    else
        return ShallowCopy(Cold);
    fi;
end;

#############################################################################
##
#F  ConstantWeightSubcode( <C> [, <w>] )  .  all words of <C> with weight <w>
##
ConstantWeightSubcode := function(arg)
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: ConstantWeightSubcode(<C> [, <weight>]");
    fi;
    IsCyclicCode(arg[1]);
    if Length(arg) = 1 then
        return arg[1].operations.ConstantWeightSubcode(arg[1]);
    else
        return arg[1].operations.ConstantWeightSubcode(arg[1], arg[2]);
    fi;
end;

CodeOps.ConstantWeightSubcode := function (arg)
    local C, D, wt, Els;
    if Length(arg) = 1 then
        C := arg[1];
        wt := PositionProperty(WeightDistribution(C){[2..WordLength(C)+1]},
                      i->i > 0);
        if wt = false then
            wt := WordLength(C);
        fi;
    elif Length(arg) = 2 then
        C := arg[1];
        wt := arg[2];
    fi;
    Els := Filtered(Elements(C), c->WeightCodeword(c) = wt);
    if Els <> [] then
        D := ElementsCode(Els,Concatenation( "code with codewords of weight ",
                             String(wt)), Field(C) );
        D.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
        D.history := History(C);
        return D;
    else
        Error("no words of weight ",wt);
    fi;
end;

LinCodeOps.ConstantWeightSubcode := function (arg)
    local C, D, wt, incode, inV, infile, i, Els;
    if Length(arg) = 1 then
        C := arg[1];
        wt := MinimumDistance(C);
    else
        C := arg[1];
        wt := arg[2];
    fi;
    if wt = 0 then
        return NullCode(WordLength(C), Field(C));
    fi;
    if Dimension(C) = 0 then
        Error("no constant weight subcode of a null code is defined");
    fi;
    incode := TmpName(); PrintTo( incode, "\n" );
    inV := TmpName(); PrintTo( inV, "\n" );
    infile := TmpName(); PrintTo( infile, "\n" );
    GuavaToLeon(C, incode);
    ExecPkg("leon","bin/wtdist",
            Concatenation("-q ",incode,"::code ",
            String(wt), " ", inV,"::code"),".");
    ExecPkg("guava","bin/leonconv",
            Concatenation("-c ",inV," ",infile),".");
    Read(infile);
    RemoveFiles(incode,inV,infile);
    Els := [];
    for i in Elements(Field(C)){[2..Size(Field(C))]} do
        Append(Els, i * GUAVA_TEMP_VAR);
    od;
    if Els <> [] then
        D := ElementsCode(Els, Concatenation( "code with codewords of weight ",
                String(wt)), Field(C) );
        D.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
        D.history := History(C);
        return D;
    else
        Error("no words of weight ",wt);
    fi;
end;

CycCodeOps.ConstantWeightSubcode := LinCodeOps.ConstantWeightSubcode;

#############################################################################
##
#F  ExtendedCode( <C> [, <i>] ) . . . . . code with added parity check symbol
##
ExtendedCode := function( arg )
    if not Length(arg) in [1..2] then
        Error("usage: ExtendedCode(<C> [, <i>])");
    fi;
    IsCyclicCode( arg[1] );
    if Length(arg) = 1 then
        return arg[1].operations.ExtendedCode(arg[1], 1);
    elif arg[2] < 1 then
        return ShallowCopy(arg[1]);
    else
        return arg[1].operations.ExtendedCode(arg[1], arg[2]);
    fi;
end;

CodeOps.ExtendedCode := function (Cold, nrcolumns ) 
    local C, word, zeros, elements, vec, i, n, q;
    n := WordLength(Cold)+1;
    q := Size(Field(Cold));
    zeros:= List( [ 2 .. nrcolumns ], i-> Field(Cold).zero );
    elements := [];
    for word in Elements(Cold) do
        vec := word.vector;
        Add(elements, Codeword(Concatenation(vec, [-Sum(vec)], zeros)));
        
    od;
    C := ElementsCode( elements, "extended code", q );
    C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(Cold);
    C.upperBoundMinimumDistance := UpperBoundMinimumDistance(Cold);
    if q = 2 then
        if LowerBoundMinimumDistance(C) mod 2 = 1 then
            C.lowerBoundMinimumDistance := C.lowerBoundMinimumDistance + 1;
        fi;
        if UpperBoundMinimumDistance(C) mod 2 = 1 then
            C.upperBoundMinimumDistance := C.upperBoundMinimumDistance + 1;
        fi;
        if IsBound(Cold.innerDistribution) then
            C.innerDistribution := NullVector(n+1);
            C.innerDistribution[1] := Cold.innerDistribution[1];
            for i in [1 .. QuoInt( WordLength(Cold), 2 ) ] do
                C.innerDistribution[i*2+1]:= Cold.innerDistribution[i*2+1]+
                                             Cold.innerDistribution[i*2];
            od;
            if IsOddInt(WordLength(Cold)) then
                C.innerDistribution[WordLength(Cold) + 2] := 
                  Cold.innerDistribution[WordLength(Cold) + 1];
            fi;
        fi;
        if IsBound(Cold.weightDistribution) then
            C.weightDistribution := NullVector( n + 1);
            C.weightDistribution[1] := Cold.weightDistribution[1];
            for i in [1 .. QuoInt( WordLength(Cold), 2 ) ] do
                C.weightDistribution[i*2+1]:=Cold.weightDistribution[i*2+1]+
                                             Cold.weightDistribution[i*2];
            od;
            if IsOddInt(WordLength(Cold)) then
                C.weightDistribution[WordLength(Cold) + 2] := 
                  Cold.weightDistribution[WordLength(Cold) + 1];
            fi;
        fi;
        if IsBound(Cold.boundsCoveringRadius)
           and Length( Cold.boundsCoveringRadius ) = 1 then
            C.boundsCoveringRadius := 
              [ Cold.boundsCoveringRadius[ 1 ] + nrcolumns ];
        fi;
    else
        C.upperBoundMinimumDistance := UpperBoundMinimumDistance(Cold) + 1;
    fi;
    C.history := History(Cold);
    return C;
end;

LinCodeOps.ExtendedCode := function(Cold, nrcolumns) 
    local C, G, word, zeros, i, n, q;

    n := WordLength(Cold) + nrcolumns;
    zeros := List( [2 .. nrcolumns], i-> Field(Cold).zero );
    q := Size(Field(Cold));
    G := Copy(GeneratorMat(Cold));
    for word in G do
        Add(word, -Sum(word));
        Append(word, zeros);
    od;
    C := GeneratorMatCode(G, "extended code", q);
    C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(Cold);
    C.upperBoundMinimumDistance := UpperBoundMinimumDistance(Cold);
    if q = 2 then
        if IsOddInt( LowerBoundMinimumDistance(C) ) then
            C.lowerBoundMinimumDistance := C.lowerBoundMinimumDistance + 1;
        fi;
        if IsOddInt( UpperBoundMinimumDistance(C) ) then
            C.upperBoundMinimumDistance := C.upperBoundMinimumDistance + 1;
        fi;
        if IsBound(Cold.weightDistribution) then
            C.weightDistribution := NullVector(n + 1);
            C.weightDistribution[1] := 1;
            for i in [ 1 .. QuoInt( WordLength(Cold), 2 ) ] do
                C.weightDistribution[i*2+1]:=Cold.weightDistribution[i*2+1]+
                                             Cold.weightDistribution[i*2];
            od;
            if IsOddInt(WordLength(Cold)) then
                C.weightDistribution[WordLength(Cold) + 2] := 
                  Cold.weightDistribution[WordLength(Cold) + 1];
            fi;
        fi;
        if IsBound(Cold.boundsCoveringRadius)
           and Length( Cold.boundsCoveringRadius ) = 1 then
            C.boundsCoveringRadius := 
              [ Cold.boundsCoveringRadius[ 1 ] + nrcolumns ];
        fi;
    else
        C.upperBoundMinimumDistance := UpperBoundMinimumDistance(Cold) + 1;
    fi;
    C.history := History(Cold);
    return C;
end;

CycCodeOps.ExtendedCode := LinCodeOps.ExtendedCode;

#############################################################################
##
#F  ShortenedCode( <C> [, <L>] )  . . . . . . . . . . . . . .  shortened code
##
ShortenedCode := function(arg)
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: ShortenedCode(C [, L])");
    fi;
    IsCyclicCode(arg[1]);
    return ApplyFunc( arg[1].operations.ShortenedCode, arg );
end;

CodeOps.ShortenedCode := function( arg )
    local C, Cnew, L, i, e, zero, q, baseels, element, max, number, 
          temp, els, n;
    C := arg[1];
    if Length(arg) = 1 then
        L := [1];
    elif Length(arg) = 2 then
        if IsList(arg[2]) then
            L := Reversed(Set(arg[2]));
        else
            L := [arg[2]];
        fi;
    else
        Error("usage: ShortenedCode(C [, L])");
    fi;
    zero := Field(C).zero;
    baseels := Elements(Field(C));
    q := Size(Field(C));
    els := VectorCodeword(Elements(C));
    for i in L do
        temp := List(els, x -> x[i]);
        max := 0;
        for e in baseels do
            number := Length(Filtered(temp, x -> x=e));
            if number > max then
                max := number;
                element := e;
            fi;
        od;
        temp := [];
        n := Length(els[1]);
        for e in els do
            if e[i] = element then
                Add(temp,Concatenation(e{[1..i-1]},e{[i+1..n]}));
            fi;
            els := temp;
        od;
    od;
    Cnew := ElementsCode(temp, "shortened code", Field(C));
    Cnew.history := History(C);
    Cnew.lowerBoundMinimumDistance := Minimum(LowerBoundMinimumDistance(C),
                                              WordLength(Cnew));
    return Cnew;
end;

LinCodeOps.ShortenedCode := function( arg )
    local C, Cnew, L, G, i, e, zero, q, baseels, temp, n;
    C := arg[1];
    if Length(arg) = 1 then
        L := [1];
    elif Length(arg) = 2 then
        if IsList(arg[2]) then
            L := Reversed(Set(arg[2]));
        else
            L := [arg[2]];
        fi;
    fi;
    zero := Field(C).zero;
    baseels := Elements(Field(C));
    q := Size(Field(C));
    G := Copy(GeneratorMat(C));
    for i in L do
        e := 0;
        repeat
            e := e + 1;
        until (e > Length(G)) or (G[e][i] <> zero);
        if G <> [] then 
            n := Length(G[1]);
        else
            n := WordLength(C);
        fi;
        if e <= Length(G) then
            temp := G[e];
            G := Concatenation(G{[1..e-1]},G{[e+1..Length(G)]});
            G := List(G, x-> x - temp * (x[i] / temp[i]));
        fi;
        G := List(G, x->Concatenation(x{[1..i-1]},x{[i+1..n]}));
        if G = [] then 
            return NullCode(WordLength(C)-Length(L), Field(C));
        fi;
    od;
    Cnew := GeneratorMatCode( BaseMat(G), "shortened code", Field(C));
    Cnew.history := History(C);
    Cnew.lowerBoundMinimumDistance := Minimum(LowerBoundMinimumDistance(C),
                                              WordLength(Cnew));
    return Cnew;
end;

CycCodeOps.ShortenedCode := function( arg )
    local C, Cnew, L, G, i, e, zero, n, q, baseels, element, max, number,
          temp, els;
    C := arg[1];
    if Length(arg) = 1 then
        L := [1];
    elif Length(arg) = 2 then
        if IsList(arg[2]) then
            L := Reversed(Set(arg[2]));
        else
            L := [arg[2]];
        fi;
    fi;
    zero := Field(C).zero;
    baseels := Elements(Field(C));
    q := Size(Field(C));
    G := Copy(GeneratorMat(C));
    for i in L do
        e := 0;
        repeat
            e := e + 1;
        until (e > Length(G)) or (G[e][i] <> zero);
        if G <> [] then
            n := Length(G[1]);
        else
            n := WordLength(C);
        fi;
        if e <= Length(G) then
            temp := G[e];
            G := Concatenation(G{[1..e-1]},G{[e+1..Length(G)]});
            G := List(G, x-> x - temp * (x[i] / temp[i]));
        fi;
        G := List(G, x->Concatenation(x{[1..i-1]},x{[i+1..n]}));
        if G = [] then 
            return NullCode(WordLength(C)-Length(L), Field(C));
        fi;
    od;
    Cnew := GeneratorMatCode(BaseMat(G), "shortened code", Field(C));
    Cnew.history := History(C);
    Cnew.lowerBoundMinimumDistance := Minimum(LowerBoundMinimumDistance(C),
                                              WordLength(Cnew));
    Cnew.upperBoundMinimumDistance := Minimum(UpperBoundMinimumDistance(C),
                                              WordLength(Cnew));
    return Cnew;
end;

#############################################################################
##
#F  PuncturedCode( <C> [, <list>] ) . . . . . . . . . . . . .  punctured code
##
##  PuncturedCode(C [, remlist]) punctures a code by leaving out the
##  coordinates given in list remlist. If remlist is omitted, then
##  the last coordinate will be removed.
##
PuncturedCode := function(arg)
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: PuncturedCode( <C> [, <list>] )");
    fi;
    IsCyclicCode(arg[1]);
    return ApplyFunc( arg[1].operations.PuncturedCode, arg );
end;

CodeOps.PuncturedCode := function (arg)
    local Cold, C, remlist, keeplist, n;
    if Length(arg) = 1 then
        Cold := arg[1];
        n := WordLength(Cold);
        remlist := [n];
        keeplist := [1..n-1];
    elif Length(arg) = 2 then
        Cold := arg[1];
        n := WordLength(Cold);
        remlist := arg[2];
        if not IsList(remlist) then
            remlist := [remlist];
        fi;
        remlist := Set(remlist);
        keeplist := [1..n];
        SubtractSet(keeplist, remlist);
    fi;
    C := ElementsCode(
                 VectorCodeword(Elements(Cold)){[1 .. Size(Cold)]}{keeplist},
                 "punctured code", Field(Cold));
    C.history := History(Cold);
    C.lowerBoundMinimumDistance := Maximum(LowerBoundMinimumDistance(Cold) -
                                           Length(remlist), 1);
    C.upperBoundMinimumDistance := Minimum( Maximum( 1,
                                           UpperBoundMinimumDistance(Cold) ), 
                                           n-Length(remlist));
    return C;
end;

LinCodeOps.PuncturedCode := function (arg)
    local Cold, C, remlist, keeplist, n;
    if Length(arg) = 1 then
        Cold := arg[1];
        n := WordLength(Cold);
        remlist := [n];
        keeplist := [1..n-1];
    elif Length(arg) = 2 then
        Cold := arg[1];
        n := WordLength(Cold);
        remlist := arg[2];
        if not IsList(remlist) then
            remlist := [remlist];
        fi;
        remlist := Set(remlist);
        keeplist := [1..n];
        SubtractSet(keeplist, remlist);
    fi;
    C := GeneratorMatCode(
                 GeneratorMat(Cold){[1..Dimension(Cold)]}{keeplist},
                 "punctured code", Field(Cold) );
    C.history := History(Cold);
    C.lowerBoundMinimumDistance := Maximum(LowerBoundMinimumDistance(Cold) -
                                           Length(remlist), 1);
    C.upperBoundMinimumDistance := Minimum( Maximum( 1, 
                                           UpperBoundMinimumDistance(Cold) ), 
                                           n - Length(remlist));
    return C;
end;

CycCodeOps.PuncturedCode := function (arg)
    local Cold, C, remlist, keeplist, n;
    if Length(arg) = 1 then
        Cold := arg[1];
        n := WordLength(Cold);
        remlist := [n];
        keeplist := [1..n-1];
    elif Length(arg) = 2 then
        Cold := arg[1];
        n := WordLength(Cold);
        remlist := arg[2];
        if not IsList(remlist) then
            remlist := [remlist];
        fi;
        remlist := Set(remlist);
        keeplist := [1..n];
        SubtractSet(keeplist, remlist);
    fi;
    C := GeneratorMatCode(
                 GeneratorMat(Cold){[1..Dimension(Cold)]}{keeplist},
                 "punctured code", Field(Cold));
    C.history := History(Cold);
    C.lowerBoundMinimumDistance := Maximum(LowerBoundMinimumDistance(Cold) -
                                           Length(remlist), 1);

# The next line is the only different one from the linear function
# Cyclic codes allways have at least one codeword with minimal weight in the
# i-th column (for every i), so if you remove a column, that word has a
# weight of one less -> the minimumdistance is reduced by one

    C.upperBoundMinimumDistance := Minimum( Maximum( 1, 
                                           UpperBoundMinimumDistance(Cold)-1),
                                           n - Length(remlist) );
    return C;
end;

#############################################################################
##
#F  ExpurgatedCode( <C>, <L> )  . . . . .  removes codewords in <L> from code
##
##  The usual way of expurgating a code is removing all words of odd weight.
##
ExpurgatedCode := function(arg)
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: ExpurgatedCode(<C> , <list>)");
    fi;
    return ApplyFunc( arg[1].operations.ExpurgatedCode, arg );
end;

CodeOps.ExpurgatedCode := function(arg)
    if not IsLinearCode(arg[1]) then
        Error("can't expurgate a non-linear code; ",
              "consider using RemovedElementsCode");
    else
        return ApplyFunc( ExpurgatedCode, arg );
    fi;
end;

LinCodeOps.ExpurgatedCode := function(arg)
    local Cold, C, L, num, H, F;
    if Length(arg) = 1 then
        C := EvenWeightSubcode(arg[1]);
        if Dimension(C) = Dimension(arg[1]) then
            # No words of even weight, so expurgate by removing first row of
            # generator matrix
            return ExpurgatedCode(arg[1], GeneratorMat(arg[1])[1]);
        else
            return C;
        fi;
    else
        Cold := arg[1];
        L := VectorCodeword( arg[2], Field(Cold) );
        if not IsList(L[1]) then 
            L := [L];
        else
            L := Set(L);
        fi;
    fi;
    F := Field(Cold);
    L := Filtered(L, l-> l in Cold);
    H := List(L, function(l)
        local V,p;
        V := NullVector(WordLength(Cold), F);
        p := PositionProperty(l, i-> not (i = F.zero));
        if not (p = false) then
            V[p] := F.one;
        fi;
        return V;
    end);
    H := BaseMat(Concatenation(CheckMat(Cold), H));
    num := Length(H) - Redundancy(Cold);
    if num > 0 then
        C := CheckMatCode( H, Concatenation("code, expurgated with ",
                     String(num), " word(s)"), F);
        C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(Cold);
        C.history := History(Cold);
        return C;
    else
        return ShallowCopy(Cold);
    fi;
end;

CycCodeOps.ExpurgatedCode := LinCodeOps.ExpurgatedCode;

#############################################################################
##
#F  AddedElementsCode( <C>, <L> ) . . . . . . . . . .  adds words in list <L>
##
AddedElementsCode := function(arg)
    if Length(arg) <> 2 then
        Error("usage: AddedElementsCode(<C> [, <list>])");
    fi;
    return arg[1].operations.AddedElementsCode(arg[1], arg[2]);
end;

CodeOps.AddedElementsCode := function(arg)
    local C, Cnew, L, e, w, E, CalcWD, wd, num;
    C := arg[1];
    L := VectorCodeword( arg[2], C );
    if not IsList(L[1]) then
        L := [ L ];
    else
        L := Set(L);
    fi;
    E := VectorCodeword(C.operations.Elements(C));
    if IsBound(C.weightDistribution) then
        wd := Copy(C.weightDistribution);
        CalcWD := true;
    else
        CalcWD := false;
    fi;
    num := 0;
    for e in L do
        if not (e in C) then
            Add(E, e);
            num := num + 1;
            if CalcWD then
                w := WeightCodeword(Codeword(e)) + 1;
                wd[w] := wd[w] + 1;
            fi;
        fi;
    od;
    if num > 0 then
        Cnew := ElementsCode(E, Concatenation( "code with ", String(num),
                " word(s) added"), Field(C));
        if CalcWD then
            Cnew.weightDistribution := wd; 
        fi;
        Cnew.history := History(C);
        return Cnew;
    else
        return ShallowCopy(C);
    fi;
end;

LinCodeOps.AddedElementsCode := CodeOps.AddedElementsCode;

CycCodeOps.AddedElementsCode := LinCodeOps.AddedElementsCode;

#############################################################################
##
#F  RemovedElementsCode( <C>, <L> ) . . . . . . . . removes words in list <L>
##
RemovedElementsCode := function(arg)
    if Length(arg) <> 2 then
        Error("usage: RemovedElementsCode(<C>, <list>)");
    fi;
    return arg[1].operations.RemovedElementsCode(arg[1], arg[2]);
end;

CodeOps.RemovedElementsCode := function(arg)
    local C, L, E, E2, e, num, s, CalcWD, wd, w, Cnew;
    C := arg[1];
    L := VectorCodeword( arg[2], C );
    if not IsList(L[1]) then
        L := [L];
    else
        L := Set(L);
    fi;
    E := Set(VectorCodeword(C.operations.Elements(C)));
    E2 := [];
    if IsBound(C.weightDistribution) then
        wd := Copy(C.weightDistribution);
        CalcWD := true;
    else
        CalcWD := false;
    fi;
    for e in E do
        if not (e in L) then
            Add(E2, e);
        elif CalcWD then
            w := WeightCodeword(Codeword(e)) + 1;
            wd[w] := wd[w] - 1;
        fi;
    od;
    num := Size(E) - Size(E2);
    if num > 0 then
        Cnew := ElementsCode(E2, Concatenation( "code with ", String(num),
                " word(s) removed"), Field(C) );
        if CalcWD then
            Cnew.weightDistribution := wd;
        fi;
        Cnew.history := History(C);
        Cnew.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
        return Cnew;
    else
        return ShallowCopy(C);
    fi;
end;

LinCodeOps.RemovedElementsCode := CodeOps.RemovedElementsCode;

CycCodeOps.RemovedElementsCode := LinCodeOps.RemovedElementsCode;

#############################################################################
##
#F  LengthenedCode( <C> [, <i>] ) . . . . . . . . . . . . . .  lengthens code
##
LengthenedCode := function( arg )
    if not Length(arg) in [1..2] then
        Error("usage: LengtenedCode(<C> [, <i>])");
    fi;
    IsCyclicCode(arg[1]);
    if Length(arg) = 1 then
        return arg[1].operations.LengthenedCode(arg[1], 1);
    else
        return arg[1].operations.LengthenedCode(arg[1], arg[2]);
    fi;
end;

CodeOps.LengthenedCode := function(C, nrcolumns)
    local Cnew;
    Cnew := ExtendedCode(AugmentedCode(C), nrcolumns);
    Cnew.history := History(C);
    Cnew.name := Concatenation("code, lengtened with ",String(nrcolumns),
            " column(s)");
    return Cnew;
end;

LinCodeOps.LengthenedCode := CodeOps.LengthenedCode;

CycCodeOps.LengthenedCode := LinCodeOps.LengthenedCode;

#############################################################################
##
#F  ResidueCode( <C> [, <w>] )  . .  takes residue of <C> with respect to <w>
##
##  If w is omitted, a word from C of minimal weight is used
##
ResidueCode := function(arg)
    local C, w, i, d;
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: ResidueCode( <C> [, <w>] )");
    fi;
    C := arg[1];
    if Length(arg) = 2 then
        w := Codeword(arg[2], C);
    else
#        w := CodewordNr(ConstantWeightSubcode(C), 1);
        d := MinimumDistance(C);
        i := 2;
        while WeightCodeword(CodewordNr(C, i)) > d do
            i := i + 1;
        od;
        w := CodewordNr(C, i);
    fi;
    if WeightCodeword(w) = 0 then
        Error("word of weight 0 is not allowed");
    elif WeightCodeword(w) = WordLength(C) then
        Error("all-one word is not allowed");
    fi;
    return C.operations.ResidueCode(C, w);
end;

CodeOps.ResidueCode := function(C, w)
    if not IsLinearCode(C) then
        Error("argument must be a linear code");
    else
        return ResidueCode(C, w);
    fi;
end;

LinCodeOps.ResidueCode := function(C, w)
    local Cnew, q, d;
    Cnew := PuncturedCode(ExpurgatedCode(C, w), Support(w));
    q := Size(Field(C));
    d := MinimumDistance(C) - WeightCodeword(w) * ( (q-1) / q );
    if not IsInt(d) then
        d := Int(d) + 1;
    fi;
    Cnew.lowerBoundMinimumDistance := d;
    Cnew.history := History(C);
    Cnew.name := "residue code";
    return Cnew;
end;

CycCodeOps.ResidueCode := LinCodeOps.ResidueCode;

#############################################################################
##
#F  ConstructionBCode( <C> )  . . . . . . . . . . .  code from construction B
##
##  Construction B (See M&S, Ch. 18, P. 9) assumes that the check matrix has
##  a first row of weight d' (the dual distance of C). The new code has a
##  check matrix equal to this matrix, but with columns removed where the
##  first row is 1.
##
ConstructionBCode := function(arg)
    if Length(arg) <> 1 then
        Error("usage: ConstructionBCode( <C> )");
    elif Field(arg[1]) <> GF(2) then
        Error("only valid for binary codes");
    fi;
    return arg[1].operations.ConstructionBCode(arg[1]);
end;

CodeOps.ConstructionBCode := function(C)
    if not IsLinearCode(C) then
        Error("only valid for linear codes");
    else
        return ConstructionBCode( C );
    fi;
end;

LinCodeOps.ConstructionBCode := function(C)
    local i, H, dd, M, mww, Cnew, DC, keeplist;
    DC := DualCode(C);
    H := Copy(GeneratorMat(DC));        # is check matrix of C
    M := Size(DC);
    dd := MinimumDistance(DC);          # dual distance of C
    i := 2;
    repeat
        mww := CodewordNr(DC, i);
        i := i + 1;
    until WeightCodeword(mww) = dd;
    i := i - 2;
    keeplist := Set([1..WordLength(C)]);
    SubtractSet(keeplist, Support(mww));
    mww := VectorCodeword(mww);
    # make sure no row dependencies arise;
    H[Redundancy(C)-LogInt(i, Size(Field(C)))] := mww;
    H := List(H, h -> h{keeplist});
    Cnew := CheckMatCode(H, Concatenation("Construction B (",String(dd),
            " coordinates)"), Field(C));
    Cnew.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
    Cnew.history := History(DC);
    return Cnew;
end;

CycCodeOps.ConstructionBCode := LinCodeOps.ConstructionBCode;

#############################################################################
##
#F  PermutedCode( <C>, <P> )  . . . . . . . permutes coordinates of codewords
##
PermutedCode := function(C, P)
    if not IsPerm(P) then
        Error("second argument must be a permutation");
    fi;
    IsLinearCode(C);
    return C.operations.PermutedCode(C, P);
end;

CodeOps.PermutedCode := function(Cold, P)
    local C, field;
    C := ElementsCode(PermutedCols(VectorCodeword(Elements(Cold)), P),
                 "permuted code", Field(Cold));
    # Copy the fields that stay the same:
    for field in ["weightDistribution", "innerDistribution",
            "lowerBoundMinimumDistance", "upperBoundMinimumDistance",
            "boundsCoveringRadius", "isPerfectCode", "isSelfDualCode"] do 
        if IsBound(Cold.(field)) then
            C.(field) := Cold.(field);
        fi;
    od;
    C.history := History(Cold);
    return C;
end;

LinCodeOps.PermutedCode := function(Cold, P)
    local C, field;
    C := GeneratorMatCode(
                 PermutedCols(GeneratorMat(Cold), P),
                 "permuted code", Field(Cold));
    # Copy the fields that stay the same:
    for field in ["weightDistribution", "boundsCoveringRadius", 
            "isPerfectCode",
            "isSelfDualCode", "lowerBoundMinimumDistance",
            "upperBoundMinimumDistance", "minimumWeightOfGenerators"]  do 
        if IsBound(Cold.(field)) then
            C.(field) := Cold.(field);
        fi;
    od;
    C.history := History(Cold);
    return C;
end;

CycCodeOps.PermutedCode := function(Cold, P)
    local C, field;
    C := GeneratorMatCode(
                 PermutedCols(GeneratorMat(Cold), P),
                 "permuted code", Field(Cold) );
    # Copy the fields that stay the same:
    for field in ["weightDistribution", "boundsCoveringRadius", 
            "isPerfectCode",
            "isSelfDualCode", "lowerBoundMinimumDistance",
            "upperBoundMinimumDistance", "minimumWeightOfGenerators",
            "upperBoundOptimalMinimumDistance"]  do 
        if IsBound(Cold.(field)) then
            C.(field) := Cold.(field);
        fi;
    od;
    C.history := History(Cold);
    return C;
end;

#############################################################################
##
#F  StandardFormCode( <C> ) . . . . . . . . . . . . standard form of code <C>
##
StandardFormCode := function(C)
    IsCyclicCode(C);
    return C.operations.StandardFormCode(C);
end;

CodeOps.StandardFormCode := function(C)
    local Cnew;
    Cnew := ShallowCopy(C);
    Cnew.elements := Set(Elements(C));
    Unbind(C.specialDecoder);
    Cnew.name := "standard form";
    Cnew.history := History(C);
    return Cnew;
end;

LinCodeOps.StandardFormCode := function(C)
    local P, str, Cnew;
    # Make a CodeCopy?
    Cnew := Copy(C);
    str := "standard form";
    GeneratorMat(Cnew);
    P := PutStandardForm(Cnew.generatorMat, true, Field(C));
    Unbind(Cnew.checkMat);
    if P <> () then
        Append(str, Concatenation(", permuted with ", String(P)));
        Unbind(Cnew.syndromeTable);
        Unbind(Cnew.outerDistribution);
        Unbind(Cnew.isCyclicCode);
        Unbind(Cnew.automorphismGroup);
        Unbind(Cnew.elements);
        Unbind(Cnew.specialDecoder);
    fi;
    Cnew.name := str;
    Cnew.history := History(C);
    return Cnew;
end;

CycCodeOps.StandardFormCode := function(C)
    local P, str, Cnew;
    # Make a CodeCopy?
    Cnew := Copy(C);
    str := "standard form";
    GeneratorMat(Cnew);
    P := PutStandardForm(Cnew.generatorMat, true, Field(C));
    Unbind(Cnew.checkMat);
    Unbind(Cnew.generatorPol);
    Unbind(Cnew.checkPol);
    Unbind(Cnew.isCyclicCode);
    Cnew.operations := LinCodeOps;
    if P <> () then
        Append(str, ", permuted with ", String(P), ",");
        Unbind(Cnew.syndromeTable);
        Unbind(Cnew.outerDistribution);
        Unbind(Cnew.isCyclicCode);
        Unbind(Cnew.automorphismGroup);
        Unbind(Cnew.elements);
        Unbind(Cnew.specialDecoder);
    fi;
    Cnew.name := str;
    Cnew.history := History(C);
    return Cnew;
end;

#############################################################################
##
#F  ConversionFieldCode( <C> )  . . . . . converts code from GF(q^m) to GF(q)
##
ConversionFieldCode := function(C)
    IsCyclicCode(C);
    return C.operations.ConversionFieldCode(C);
end;

CodeOps.ConversionFieldCode := function(C)
    local F,x,q,m,i, ConvertElms, Cnew;

    ConvertElms := function (M)
        local res,n,k,vec,coord, ConvTable, Nul, g, zero;
        res := [];
        n := Length(M[1]);
        k := Length(M);
        g := Polynomial(GF(q), MinPol(F.root));
        zero := F.zero;
        Nul := List([1..m], i -> zero);
        ConvTable := [];
        x := Indeterminate(GF(q));
        for i in [1..Size(F) - 1] do
            ConvTable[i] := VectorCodeword(x^(i-1) mod g, m);
        od;
        for vec in [1..k] do
            res[vec] := [];
            for coord in [1..n] do
                if M[vec][coord] <> zero then
                    Append(res[vec], ConvTable[LogFFE(M[vec][coord]) + 1]);
                else
                    Append(res[vec], Nul);
                fi;
            od;
        od;
        return res;
    end;

    F := Field(C);
    q := F.char;
    m := F.degree;
    Cnew := ElementsCode(
                    ConvertElms(VectorCodeword(Elements(C))),
                    Concatenation("code, converted to basefield GF(",
                            String(q),")"),
                    F );
    Cnew.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
    Cnew.upperBoundMinimumDistance := Minimum(WordLength(C), 
                                              m*UpperBoundMinimumDistance(C));
    Cnew.history := History(C);
    return Cnew;
end;

LinCodeOps.ConversionFieldCode := function(C)
    local F, Cnew;
    F := Field(C);
    Cnew := GeneratorMatCode(
                    HorizontalConversionFieldMat( GeneratorMat(C), F),
                    Concatenation("code, converted to basefield GF(",
                            String(F.char), ")"),
                    GF(F.char));
    Cnew.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
    Cnew.upperBoundMinimumDistance := Minimum(WordLength(C), 
                                    F.degree * UpperBoundMinimumDistance(C));
    Cnew.history := History(C);
    return Cnew;
end;

CycCodeOps.ConversionFieldCode := LinCodeOps.ConversionFieldCode;

#############################################################################
##
#F  CosetCode( <C>, <f> ) . . . . . . . . . . . . . . . . . . .  coset of <C>
##
CosetCode := function(arg)
    IsCyclicCode(arg[1]);
    if Length(arg) <> 2 then
        Error("usage: CosetCode( <C>, <f> )");
    fi;
    return arg[1].operations.CosetCode(arg[1], arg[2]);
end;

CodeOps.CosetCode := function(C, f)
    local i, els, Cnew;
    f := Codeword(f, Field(C));
    els := [];
    for i in [1..Size(C)] do
        Add(els, C.elements[i] + f);
    od;
    Cnew := ElementsCode(els, "coset code", Field(C) );
    if IsBound(C.weightDistribution) and IsLinearCode(C) then
        Cnew.innerDistribution := Copy(C.weightDistribution);
    fi;
    Cnew.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
    Cnew.upperBoundMinimumDistance := UpperBoundMinimumDistance(C);
    Cnew.history := History(C);
    return Cnew;
end;

LinCodeOps.CosetCode := function(C, f)
    local Cnew, i, els, Cels;
    f := Codeword(f, Field(C));
    if f in C then
        return ShallowCopy(C);
    fi;
    els := [];
    Cels := C.operations.Elements(C);
    for i in [1..Size(C)] do
        Add(els, Cels[i] + f);
    od;
    Cnew := ElementsCode(els, "coset code", Field(C) );
    if IsBound(C.weightDistribution) and IsLinearCode(C) then
        Cnew.innerDistribution := Copy(C.weightDistribution);
    fi;
    Cnew.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C);
    Cnew.upperBoundMinimumDistance := UpperBoundMinimumDistance(C);
    Cnew.history := History(C);
    return Cnew;
end;

CycCodeOps.CosetCode := LinCodeOps.CosetCode;

#############################################################################
##
#F  DirectSumCode( <C1>, <C2> ) . . . . . . . . . . . . . . . . .  direct sum
##
##  DirectSumCode(C1, C2) creates a (n1 + n2 , M1 M2 , min{d1 , d2} ) code
##  by adding each codeword of the second code to all the codewords of the
##  first code.
##
DirectSumCode := function (C1, C2)
    local LinearDirectSumCode, UnrestrictedDirectSumCode;
     
UnrestrictedDirectSumCode := function (C1, C2) 
    local i, j, C, els, G, n, sumcr;
    els := [];
    for i in VectorCodeword(C1.operations.Elements(C1)) do
        Append(els,List(VectorCodeword(C2.operations.Elements(C2)),
                x-> Concatenation(i, x ) ) );
    od;
    C := ElementsCode( els, "direct sum code", Field(C1) );
    n := WordLength(C1) + WordLength(C2);
    if Size(C) <= 1 then
        C.lowerBoundMinimumDistance := n;
        C.upperBoundMinimumDistance := n;
    else
        C.lowerBoundMinimumDistance := Minimum(LowerBoundMinimumDistance(C1),
                                               LowerBoundMinimumDistance(C2));
        C.upperBoundMinimumDistance := Minimum(UpperBoundMinimumDistance(C1),
                                               UpperBoundMinimumDistance(C2));
    fi;
    if IsBound(C1.weightDistribution) and IsBound(C2.weightDistribution) then
        C.weightDistribution := NullVector(WordLength(C)+1);
        for i in [1..WordLength(C1)+1] do
            for j in [1..WordLength(C2)+1] do
                C.weightDistribution[i+j-1] := C.weightDistribution[i+j-1]+
                                             WeightDistribution(C1)[i] *
                                             WeightDistribution(C2)[j];
            od;
        od;
    fi;
    if IsBound(C1.innerDistribution) and IsBound(C2.innerDistribution) then
        C.innerDistribution := NullVector(WordLength(C) + 1);
        for i in [1..WordLength(C1) + 1 ] do
            for j in [1..WordLength(C2) + 1 ] do
                C.innerDistribution[i+j-1] := C.innerDistribution[i+j-1]+
                                            InnerDistribution(C1)[i] *
                                            InnerDistribution(C2)[j];

            od;
        od;
    fi;
    if IsBound(C1.boundsCoveringRadius) 
       and IsBound(C2.boundsCoveringRadius) then
        sumcr := List( C1.boundsCoveringRadius,
          x -> x + C2.boundsCoveringRadius );
        sumcr := Set( Flat( sumcr ) );
        IsRange( sumcr );
        C.boundsCoveringRadius := sumcr;
    fi;
    C.history := MergeHistories(History(C1), History(C2));
    return C;
end;

LinearDirectSumCode := function (C1, C2) 
    local i, j, C, zeros1, zeros2, G, n, sumcr;
    zeros1 := NullVector(WordLength(C1), Field(C1));
    zeros2 := NullVector(WordLength(C2), Field(C1));
    G := List(GeneratorMat(C1),x -> Concatenation(x,zeros2));
    Append(G,List(GeneratorMat(C2),x -> Concatenation(zeros1,x)));
    C := GeneratorMatCode( G, "direct sum code", Field(C1) );
    n := WordLength(C1) + WordLength(C2);
    if Size(C) <= 1 then
        C.lowerBoundMinimumDistance := n;
        C.upperBoundMinimumDistance := n;
    else
        C.lowerBoundMinimumDistance := Minimum(LowerBoundMinimumDistance(C1),
                                               LowerBoundMinimumDistance(C2));
        C.upperBoundMinimumDistance := Minimum(UpperBoundMinimumDistance(C1),
                                               UpperBoundMinimumDistance(C2));
    fi;
    if IsBound(C1.weightDistribution) and IsBound(C2.weightDistribution) then
        C.weightDistribution := NullVector(WordLength(C)+1);
        for i in [1..WordLength(C1)+1] do
            for j in [1..WordLength(C2)+1] do
                C.weightDistribution[i+j-1] := C.weightDistribution[i+j-1]+
                                             WeightDistribution(C1)[i] *
                                             WeightDistribution(C2)[j];
            od;
        od;
    fi;
    if IsBound(C1.boundsCoveringRadius) 
       and IsBound(C2.boundsCoveringRadius) then
        sumcr := List( C1.boundsCoveringRadius,
          x -> x + C2.boundsCoveringRadius );
        sumcr := Set( Flat( sumcr ) );
        IsRange( sumcr );
        C.boundsCoveringRadius := sumcr;
    fi;
    if IsBound( C1.isNormalCode ) and IsBound( C2.isNormalCode ) and
       IsNormalCode( C1 ) and IsNormalCode( C2 ) then
        C.isNormalCode := true;
    fi;
    if IsBound(C1.isSelfOrthogonalCode) and IsBound(C2.isSelfOrthogonalCode) 
       and IsSelfOrthogonalCode(C1) and IsSelfOrthogonalCode(C2) then
        C.isSelfOrthogonalCode := true;
    fi;
    C.history := MergeHistories(History(C1), History(C2));
    return C;
end;

    if Field(C1) <> Field(C2) then
        Error("Codes are not in the same basefield");
    fi;
    if IsLinearCode(C1) and IsLinearCode(C2) then
        return LinearDirectSumCode(C1, C2);
    else
        return UnrestrictedDirectSumCode(C1, C2);
    fi;
end;

#############################################################################
##
#F  ConcatenationCode( <C1>, <C2> ) . . . . .  concatenation of <C1> and <C2>
##
ConcatenationCode := function(C1, C2)
    local LinearConcatenationCode, UnrestrictedConcatenationCode;

    UnrestrictedConcatenationCode := function(C1, C2)
        local E, e, C;
        E := [];
        for e in [1..Size(C1)] do
            Add(E,Concatenation(VectorCodeword(Elements(C1)[e]),
                    VectorCodeword(Elements(C2)[e])));
        od;
        C := ElementsCode( E, "concatenation code", Field(C1) );
        C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C1) +
                                       LowerBoundMinimumDistance(C2);
        C.upperBoundMinimumDistance := UpperBoundMinimumDistance(C1) +
                                       UpperBoundMinimumDistance(C2);
        # CoveringRadius?
        C.history := MergeHistories(History(C1), History(C2));
        return C;
    end;
    
    LinearConcatenationCode := function(C1, C2)
        # When this function is called, both arguments are linear codes.
        local E, e, C, G, G1, G2;
        G := [];
        G1 := GeneratorMat(C1);
        G2 := GeneratorMat(C2);
        for e in [1..Dimension(C1)] do
            Add(G, Concatenation(G1[e], G2[e]));
        od;
        C := GeneratorMatCode(G, "concatenation code", Field(C1));
        C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C1) +
                                       LowerBoundMinimumDistance(C2);
        C.upperBoundMinimumDistance := UpperBoundMinimumDistance(C1) +
                                       UpperBoundMinimumDistance(C2);
        # CoveringRadius?
        C.history := MergeHistories(History(C1), History(C2));
        return C;
    end;
    
    if Size(C1) <> Size(C2) then
        Error("both codes must have equal size");
    elif Field(C1) <> Field(C2) then
        Error("both codes must be over the same field");
    fi;
    if IsLinearCode(C2) and IsLinearCode(C1) then
        return LinearConcatenationCode(C1, C2);
    else
        return UnrestrictedConcatenationCode(C1, C2);
    fi;
end;

#############################################################################
##
#F  DirectProductCode( <C1>, <C2> ) . . . . . . . . . . . . .  direct product
##
##  DirectProductCode constructs a new code from the direct product of two
##  codes by taking the Kronecker product of the two generator matrices
##
DirectProductCode := function (C1,C2)
    local LinearDirectProductCode;

    LinearDirectProductCode := function (C1, C2)
        local C;
        if Field(C1) <> Field(C2) then
            Error("both codes must have the same basefield");
        fi;
        C := GeneratorMatCode(
                     KroneckerProduct(GeneratorMat(C1), GeneratorMat(C2)),
                     "direct product code",
                     Field(C1));
        if Dimension(C) = 0 then
            C.lowerBoundMinimumDistance := WordLength(C);
            C.upperBoundMinimumDistance := WordLength(C);
        else
            C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C1) *
                                           LowerBoundMinimumDistance(C2);
            C.upperBoundMinimumDistance := UpperBoundMinimumDistance(C1) *
                                           UpperBoundMinimumDistance(C2);
        fi;
        C.history := MergeHistories(History(C1), History(C2));
        if IsBound( C1.boundsCoveringRadius ) and
           IsBound( C2.boundsCoveringRadius ) then
            C.boundsCoveringRadius := [
              Maximum( WordLength( C1 ) * C2.boundsCoveringRadius[ 1 ],
                       WordLength( C2 ) * C1.boundsCoveringRadius[ 1 ] )
              .. GeneralUpperBoundCoveringRadius( C ) ];
        fi;
        return C;
    end;
    
    if not IsLinearCode(C1) or not IsLinearCode(C2) then
        Error("both codes must be linear");
    else
        return LinearDirectProductCode(C1, C2);
    fi;
end;

#############################################################################
##
#F  UUVCode( <C1>, <C2> ) . . . . . . . . . . . . . . .  u | u+v construction
##
##  Uuvcode(C1, C2) # creates a ( 2n , M1 M2 , d = min{2 d1 , d2} ) code
##  with codewords  (u | u + v) for all u in C1 and v in C2
##
UUVCode := function (C1, C2) 
    local LinearUUVCode, UnrestrictedUUVCode;

    LinearUUVCode := function(C1, C2)
        local C, F, diff, zeros, zeros2, G, n; 
        F := Field(C1); 
        n := WordLength(C1)+Maximum(WordLength(C1),WordLength(C2));
        diff := WordLength(C1)-WordLength(C2);
        zeros := NullVector(WordLength(C1), F);
        zeros2 := NullVector(AbsInt(diff), F);
        if diff < 0 then
            G := List(GeneratorMat(C1),u -> Concatenation(u,u,zeros2));
            Append(G,List(GeneratorMat(C2),v -> Concatenation(zeros,v)));
        else
            G:=List(GeneratorMat(C1),u -> Concatenation(u,u));
            Append(G,List(GeneratorMat(C2),v->Concatenation(zeros,v,zeros2)));
        fi;
        C := GeneratorMatCode( G, "U|U+V construction code", F );
        if Dimension(C1) = 0 then
            if Dimension(C2) = 0 then
                C.lowerBoundMinimumDistance := n;
                C.upperBoundMinimumDistance := n;
            else
                C.lowerBoundMinimumDistance := LowerBoundMinimumDistance(C2);
                C.upperBoundMinimumDistance := UpperBoundMinimumDistance(C2);
            fi;
        elif Dimension(C2) = 0 then
            C.lowerBoundMinimumDistance := 2*LowerBoundMinimumDistance(C1);
            C.upperBoundMinimumDistance := 2*UpperBoundMinimumDistance(C1);
        else
            C.lowerBoundMinimumDistance := Minimum(
                2*LowerBoundMinimumDistance(C1),LowerBoundMinimumDistance(C2));
            C.upperBoundMinimumDistance := Minimum(
                2*UpperBoundMinimumDistance(C1),UpperBoundMinimumDistance(C2));
        fi;
        
        C.history := MergeHistories(History(C1), History(C2));
        return C;
    end;

    UnrestrictedUUVCode := function(C1, C2)
        local C, F, i, M1, diff, zeros, extended, Els1, Els2, els, n; 
        F := Field(C1); 
        n := WordLength(C1)+Maximum(WordLength(C1),WordLength(C2));
        diff := WordLength(C1)-WordLength(C2);
        M1 := Size(C1);
        zeros := NullVector(AbsInt(diff), F);
        els := [];
        Els1 := Elements(C1);
        Els2 := Elements(C2);
        if diff>0 then 
            extended := List(Els2,x->Concatenation(x.vector,zeros)); 
            for i in [1..M1] do 
                Append(els, List(extended,x-> Concatenation(Els1[i].vector, 
                        Els1[i].vector + x ) ) );
            od; 
        elif diff<0 then 
            for i in [1..M1] do 
                extended := Concatenation(Els1[i].vector,zeros);
                Append(els,List(Els2,x-> Concatenation(Els1[i].vector,
                        extended + x.vector ) ) ); 
            od; 
        else 
            for i in [1..M1] do
                Append(els,List(Els2,x-> Concatenation(Els1[i].vector,
                        Els1[i].vector + x.vector ) ) );
            od; 
        fi;
        C := ElementsCode(els, "U|U+V construction code", F);
        C.lowerBoundMinimumDistance := Minimum(2*LowerBoundMinimumDistance(C1),
                                               LowerBoundMinimumDistance(C2));
        C.upperBoundMinimumDistance := Minimum(2*UpperBoundMinimumDistance(C1),
                                               UpperBoundMinimumDistance(C2));

        C.history := MergeHistories(History(C1), History(C2));
        return C;
    end;

    if Field(C1)<>Field(C2) then 
        Error("Codes are not in the same basefield"); 
    elif IsLinearCode(C1) and IsLinearCode(C2) then
        return LinearUUVCode(C1, C2);
    else
        return UnrestrictedUUVCode(C1, C2);
    fi;
end;

#############################################################################
##
#F  UnionCode( <C1>, <C2> ) . . . . . . . . . . . . .  union of <C1> and <C2>
##
UnionCode := function (C1,C2)
    local LinearUnionCode;

    LinearUnionCode := function (C1,C2)
        local C, Els, e;
        if Field(C1) <> Field(C2) then
            Error("codes are not in the same basefield");
        elif WordLength(C1) <> WordLength(C2) then
            Error("wordlength must be the same");
        fi;
        C := AugmentedCode(C1, GeneratorMat(C2));
        C.upperBoundMinimumDistance := Minimum(UpperBoundMinimumDistance(C1),
                                               UpperBoundMinimumDistance(C2));
        C.history := MergeHistories(History(C1), History(C2));
        C.name := "union code";
        return C;
    end;

# should there be a special function for cyclic codes? Or in other words:
# If C1 and C2 are cyclic, does that mean that UnionCode(C1, C2) is cyclic?

    if IsLinearCode(C2) and IsLinearCode(C1) then
        return LinearUnionCode(C1, C2);
    else
        Error ("use AddedElementsCode for non-linear codes");
    fi;
end;

#############################################################################
##
#F  IntersectionCode( <C1>, <C2> )  . . . . . . intersection of <C1> and <C2>
##
IntersectionCode := function (C1,C2)
    local UnrestrictedIntersectionCode, LinearIntersectionCode,
          CyclicIntersectionCode;
    
    UnrestrictedIntersectionCode := function(C1, C2)
        local C, Els, e;
        Els := [];
        for e in Elements(C1) do
            if e in C2 then Add(Els, e); fi;
        od;
        if Els = [] then
            return false; # or an Error?
        else
            C := ElementsCode(Els, "intersection code", Field(C1));
        fi;
        C.lowerBoundMinimumDistance := Maximum(LowerBoundMinimumDistance(C1),
                                               LowerBoundMinimumDistance(C2));
        C.history := MergeHistories(History(C1), History(C2));
        return C;
    end;
    
    LinearIntersectionCode := function (C1,C2)
        local C;
        C := DualCode(AugmentedCode(DualCode(C1), CheckMat(C2)));
        C.lowerBoundMinimumDistance := Maximum(LowerBoundMinimumDistance(C1),
                                               LowerBoundMinimumDistance(C2));
        C.history := MergeHistories(History(C1), History(C2));
        C.name := "intersection code";
        return C;
    end;
    
    CyclicIntersectionCode := function (C1,C2)
        local C;
        C := GeneratorPolCode(
                     Lcm(GeneratorPol(C1), GeneratorPol(C2)), WordLength(C1),
                     "intersection code", Field(C1) );
        if IsBound(C1.roots) and IsBound(C2.roots) then
            C.roots := UnionSet(C1.roots, C2.roots);
        fi;
        C.lowerBoundMinimumDistance := Maximum(LowerBoundMinimumDistance(C1),
                                               LowerBoundMinimumDistance(C2));
        C.history := MergeHistories(History(C1), History(C2));
        return C;
    end;
    
    if Field(C1) <> Field(C2) then
        Error("basefield must be the same.");
    elif WordLength(C1) <> WordLength(C2) then
        Error("wordlength must be the same");
    fi;
    # The 'lowest' function will be used.
    if IsCyclicCode(C1) and IsCyclicCode(C2) then
        return CyclicIntersectionCode(C1, C2);
    elif IsLinearCode(C1) and IsLinearCode(C2) then
        return LinearIntersectionCode(C1, C2);
    else
        return UnrestrictedIntersectionCode(C2, C1);
    fi;
end;
