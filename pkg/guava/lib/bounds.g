#############################################################################
##
#A  bounds.g                GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains functions for calculating with bounds
##
#H  $Log: bounds.g,v $
#H  Revision 1.3  1997/01/21 13:45:45  gap
#H  vfelsch eliminated all calls of LOADED_PACKAGES.guava
#H
#H  Revision 1.47  1994/11/10  11:10:10  jcramwin
#H  removed CreateBoundsTable
#H
#H  Revision 1.46  1994/11/10  10:55:03  jcramwin
#H  put CreateBoundsTable back in
#H
#H  Revision 1.45  1994/11/09  17:53:17  rbaart
#H  BoundsMinimumDistance: changed An into an
#H  Moved ComplexTable to file tblgener.g (CreateBoundsTable)
#H
#H  Revision 1.44  1994/11/09  16:43:44  rbaart
#H  Added OpsOps to BoundsOps record
#H
#H  Revision 1.43  1994/11/07  15:08:39  jcramwin
#H  Oh, what did I do?
#H
#H  Revision 1.42  1994/11/07  11:26:20  jcramwin
#H  - changed the way the tables are stored and some things with printing
#H
#H  Revision 1.41  1994/11/04  16:45:05  jcramwin
#H  fixed a bug in BoundsMinimumDistance, obj.cons was not declared
#H
#H  Revision 1.40  1994/11/04  13:39:32  jcramwin
#H  added some comments
#H
#H  Revision 1.39  1994/10/27  08:30:38  rbaart
#H  BoundsMinimumDistance: puncturing -> extending in InfoLine
#H
#H  Revision 1.38  1994/10/24  14:51:40  jcramwin
#H  I spelled something wrong, sorry...
#H
#H  Revision 1.37  1994/10/21  14:51:37  jcramwin
#H  Bounds now returns a (bad) value for to big n, k or q
#H
#H  Revision 1.36  1994/10/20  17:09:25  rbaart
#H  BoundsMinimumDistance: fixed buggy in obj.lastman
#H
#H  Revision 1.35  1994/10/20  16:10:33  jcramwin
#H  changed something with the printing of bounds
#H
#H  Revision 1.34  1994/10/20  11:26:18  jcramwin
#H  BoundsMinimumDistance doesn't do ExtendedCode of ExtendedCode of ExtendedCode.
#H
#H  Revision 1.33  1994/10/19  11:05:58  rbaart
#H  BoundsOps.Print changed
#H
#H  Revision 1.32  1994/10/19  10:40:42  jcramwin
#H  changed the way .construction stores the code
#H
#H  Revision 1.31  1994/10/18  16:06:00  rbaart
#H  Added BoundsOps
#H
#H  Revision 1.30  1994/10/18  15:34:29  jcramwin
#H  Played around with ConversionFieldCode
#H
#H  Revision 1.29  1994/10/18  08:55:08  rbaart
#H  Small changes in comments
#H
#H  Revision 1.28  1994/10/17  15:54:27  jcramwin
#H  BoundsMinimumDistance now creates a field .construction
#H
#H  Revision 1.27  1994/10/14  16:19:49  jcramwin
#H  changed the name of the file for the lowerbound in codesq.g
#H
#H  Revision 1.26  1994/10/14  15:23:34  rbaart
#H  BoundsMinimumDistance: fixed bug in u|u+v construction of lower bounds
#H
#H  Revision 1.25  1994/10/14  15:10:33  rbaart
#H  BoundsMinimumDistance: fixed bug with .upperBoundMinimumDistance
#H
#H  Revision 1.24  1994/10/14  14:53:47  jcramwin
#H  the list with upperbounds are moved from /shared/points/ to /guava/tbl/
#H
#H  Revision 1.23  1994/10/14  11:23:36  rbaart
#H  Added BoundsMinimumDistance and changed it a lot
#H
#H  Revision 1.22  1994/10/14  08:41:47  jcramwin
#H  too much to mention
#H
#H  Revision 1.21  1994/10/11  14:11:31  jcramwin
#H  Made Construction B for upperbound sharper
#H
#H  Revision 1.20  1994/10/11  13:14:44  jcramwin
#H  speeded up the unupspeedable: Griesmer in InitialTable in ComplexTable
#H
#H  Revision 1.19  1994/10/11  10:02:33  jcramwin
#H  some cosmetic changed inComplexTable
#H
#H  Revision 1.18  1994/10/11  09:22:40  rbaart
#H  Fixed bugs in ComplexTable
#H
#H  Revision 1.17  1994/10/10  18:45:47  jcramwin
#H  renamed SimpleTable to ComplexTable
#H
#H  Revision 1.16  1994/10/07  15:05:15  rbaart
#H  Worked on SimpleTable
#H
#H  Revision 1.15  1994/10/06  16:14:39  rbaart
#H  Added residue to upperbounds ruleslist
#H
#H  Revision 1.14  1994/10/06  14:16:53  jcramwin
#H  de upperbound in SimpleTable geintegreerd
#H
#H  Revision 1.13  1994/10/06  11:27:37  jcramwin
#H  some small changes in SimpleTable
#H
#H  Revision 1.12  1994/10/05  15:36:26  jcramwin
#H  if a n, k came from extending an other code, the table has no entry
#H
#H  Revision 1.11  1994/10/05  14:57:51  jcramwin
#H  did some more testing
#H
#H  Revision 1.10  1994/10/05  13:24:53  jcramwin
#H  tried some more rules but non of them did any good.....
#H
#H  Revision 1.9  1994/09/30  16:20:00  jcramwin
#H  changed a path name
#H
#H  Revision 1.8  1994/09/30  15:35:47  jcramwin
#H  SimpleTable rewritten
#H
#H  Revision 1.7  1994/09/29  18:50:52  jcramwin
#H  NewTable new works! ( SimpleTable(i) = NewTable(i) )
#H
#H  Revision 1.6  1994/09/29  15:36:31  jcramwin
#H  added a new function which calculates the table in an other way
#H
#H  Revision 1.5  1994/09/29  12:40:07  rbaart
#H  Check size of table on updating with points
#H
#H  Revision 1.4  1994/09/29  11:31:39  jcramwin
#H  UUV added (or did I do that before?) to SimpleTable
#H  concatenation added to SimpleTable
#H  Changed the setup of ST to RulesList
#H
#H  Revision 1.3  1994/09/29  10:03:56  rbaart
#H  Adjusted path of pointlst.txt
#H
#H  Revision 1.2  1994/09/29  10:01:21  rbaart
#H  Simple table added
#H
#H  Revision 1.1  1994/09/28  09:52:53  jcramwin
#H  Initial revision
#H
##

#############################################################################
##
#F  UpperBoundHamming( <n>, <d>, <q> )  . . . . . . . . . . . . Hamming bound
##
UpperBoundHamming := function (n,d,q)
    return Int((q^n)/(SphereContent(n,QuoInt(d-1,2),q)));
end;

#############################################################################
##
#F  UpperBoundSingleton( <n>, <d>, <q> )  . . . . . . . . . . Singleton bound
##
UpperBoundSingleton := function (n,d,q)
    return q^(n - d + 1);
end;

#############################################################################
##
#F  UpperBoundPlotkin( <n>, <d>, <q> )  . . . . . . . . . . . . Plotkin bound
##
UpperBoundPlotkin := function (n,d,q)
    local t, fact;
    t := 1 - 1/q;
    if (q=2) and (n = 2*d) and (d mod 2 = 0) then
        return 4*d;
    elif (q=2) and (n = 2*d + 1) and (d mod 2 = 1) then
        return 4*d + 4;
    elif d >= t*n + 1 then
        return Int(d/( d - t*n));
    elif d < t*n + 1 then
        fact := (d-1) / t;
        if not IsInt(fact) then
            fact := Int(fact) + 1;
        fi;
        return Int(d/( d - t * fact)) * q^(n - fact);
    fi;
end;

#############################################################################
##
#F  UpperBoundGriesmer( <n>, <d>, <q> ) . . . . . . . . . . .  Griesmer bound
##
UpperBoundGriesmer := function (n,d,q)
    local s, den, k, add;
    den := 1;
    s := 0;
    k := 0;
    add := 0;
    while s <= n do
        if add <> 1 then
            add := QuoInt(d, den) + SignInt(d mod den);
        fi;
        s := s + add;
        den := den * q;
        k := k + 1;
    od;
    return q^(k-1);
end;

#############################################################################

##
#F  UpperBoundElias( <n>, <d>, <q> )  . . . . . . . . . . . . . . Elias bound
##
UpperBoundElias := function (n,d,q)
    local t, r, tn, tnd, tnr, num, den, curmin, res, space, totspace;
    tn := Int( (q-1)*n / q);
    tnd := Int( (q-1)*n*d / q);
    num := tnd * q^n;
    space := 1;
    totspace := 1;
    curmin := q^n;
    for r in [0..tn] do
        tnr := Int ( 2*(q-1)*n*r / q);
        if tnd - tnr + r^2 > 0 then
            den := (tnd - tnr + r^2) * totspace;
            res := Int(num / den);
            space := (space * (n-r)) / (r+1);
            totspace := totspace + space;
            if res < curmin then
                curmin := res;
            fi;
        else
            r := tn + 1;
        fi;
    od;
    return curmin;
end;

#############################################################################
##
#F  UpperBoundJohnson( <n>, <d> ) . . . . . . . . . . Johnson bound for <q>=2
##
UpperBoundJohnson := function (n,d)
    local UBConsWgt, e, num, den;
    UBConsWgt := function (n1,d1)
        local fact, e, res, t;
        e := Int((d1-1) / 2);
        res := 1;
        for t in [0..e] do
            res := Int(res * (n1 - (e-t)) / ( d1 - (e-t)));
        od;
        return res;
    end;
    e := Int((d-1) / 2);
    num := Binomial(n,e+1) - Binomial(d,e)*UBConsWgt(n,d);
    den := Int(num / Int(n / (e+1)));
    return Int(2^n / (den + SphereContent(n,e,2)));
end;

#############################################################################
##
#F  UpperBound( <n>, <d> [, <F>] )  . . . .  upper bound for minimum distance
##
##  calculates upperbound for a code C of word length n, minimum distance at
##  least d over an alphabet Q of size q, using the minimum of the Hamming,
##  Plotkin and Singleton bound.
##
UpperBound := function (arg)
    local n,d,q,MinBound, l;

    MinBound := function (n1,d1,q1)
        local mn1;
        mn1 := Minimum (UpperBoundPlotkin(n1,d1,q1),
                       UpperBoundSingleton(n1,d1,q1),
                       UpperBoundElias(n1,d1,q1));
        if q1 = 2 then
            return Minimum(mn1, UpperBoundJohnson(n1,d1));
        else
            return Minimum(mn1, UpperBoundHamming(n1,d1,q1));
        fi;
    end;

    l := Length(arg);
    if l>3 or l<2 then
        Error("usage: UpperBound (<n>, <d> [, <F> ])");
    fi;
    n := arg[1];
    d := arg[2];
    if l > 2 and IsInt(arg[3]) then
        q := arg[3];
    elif l > 2 and IsField(arg[3]) then
        q := Size(arg[3]);
    else
        q := 2;
    fi;   
    if n < d then
        return 0;
    elif n = d then
        return q;
    elif d = 1 then
        return q^n;
    fi;
    if (q=2) then
        if d mod 2 = 0 then
            return Minimum(MinBound(n,d,q), MinBound(n-1,d-1,q));
        else
            return Minimum(MinBound(n,d,q), MinBound(n+1,d+1,q));
        fi;
    else
        return MinBound(n,d,q);
    fi;
end;

#############################################################################
##
#F  IsPerfectCode( <C> )  . . . . . .  determines whether C is a perfect code
##
IsPerfectCode := function(C)
    local n, q, dist, d, t, IsTrivialPerfect, ArePerfectParameters;

    IsTrivialPerfect := function(C)
        # Checks if C has only one or zero codewords, or is the whole
        # space, or is a repetition code of odd length over GF(2).
        # These are 'trivial' perfect codes.
        return ((Size(C) <= 1) or 
                (Size(C) = Size(Field(C))^WordLength(C)) or
                ((Field(C) = GF(2)) and (Size(C) = 2) and
                 (Mod(WordLength(C),2) <> 0) and (IsCyclicCode(C))));
    end;

    ArePerfectParameters := function(q, n, M, dvec) 
        local k, r;
        # Can the parameters be of a perfect code? If they don't belong
        # to a trivial perfect code, they should be the same as a Golay
        # or Hamming code.
        k := LogInt(M, q);
        if M <> q^k then 
            return false; #nothing wrong here
        elif (q = 2) and (n = 23) then
            return (k = 12) and (7 in [dvec[1]..dvec[2]]);
        elif (q = 3) and (n = 11) then
            return (k = 6) and (5 in [dvec[1]..dvec[2]]);
        else
            r := n-k;
            return (n = ((q^r-1)/(q-1))) and 
                   (3 in [dvec[1]..dvec[2]]);
        fi;
    end;

    if IsBound(C.isPerfectCode) then
        return C.isPerfectCode;
    fi;
    n := WordLength(C);
    q := Size(Field(C));
    dist := [LowerBoundMinimumDistance(C), UpperBoundMinimumDistance(C)];
    if IsTrivialPerfect(C) then
        C.isPerfectCode := true;
        if (C.size > 1) then
            C.coveringRadius := Int(MinimumDistance(C)/2);
        else
            C.coveringRadius := n;
        fi;
    elif not ArePerfectParameters(q, n, Size(C), dist) then
        C.isPerfectCode := false;
    else
        t := List(dist, d->QuoInt(d-1, 2));
        if t[1] = t[2] then 
            d := t[1]*2 +1;
        else
            d := MinimumDistance(C);
        fi;
        C.isPerfectCode := (d mod 2 = 1) and 
                           ArePerfectParameters(q, n, Size(C), [d,d]);
        if C.isPerfectCode then
            C.lowerBoundMinimumDistance := d;
            C.upperBoundMinimumDistance := d;
            C.coveringRadius := Int(d/2);
        fi;
    fi;
    return C.isPerfectCode;
end;

#############################################################################
##
#F  IsMDSCode( <C> )  . . .  checks if C is a Maximum Distance Separable Code
##
IsMDSCode := function(C)
    local wd, w, n, d, q;
    q:= Size(Field(C));
    n:= WordLength(C);
    d:= MinimumDistance(C);
    if d = n - LogInt(Size(C),q) + 1 then
        if not IsBound(C.weightDistribution) then
            wd := List([0..n], i -> 0);
            wd[1] := 1;
            for w in [d..n] do
                # The weight distribution of MDS codes is exactly known
                wd[w+1] := Binomial(n,w)*Sum(List([0..w-d],j -> 
                                   (-1)^j * Binomial(w,j) *(q^(w-d+1-j)-1)));
            od;
            C.weightDistribution := wd;
        fi;
        return true;
    else
        return false; #this is great!
    fi;
end;

#############################################################################
##
#F  OptimalityCode( <C> ) . . . . . . . . . .  estimate for optimality of <C>
##
##  OptimalityCode(C) returns the difference between the smallest known upper-
##  bound and the actual size of the code. Note that the value of the
##  function UpperBound is not allways equal to the actual upperbound A(n,d)
##  thus the result may not be equal to 0 for all optimal codes!
##
OptimalityCode := function(C)
    return UpperBound(WordLength(C), MinimumDistance(C), Size(Field(C))) - 
Size(C);
end;

#############################################################################
##
#F  OptimalityLinearCode( <C> ) .  estimate for optimality of linear code <C>
##
##  OptimalityLinearCode(C) returns the difference between the smallest known
##  upperbound on the size of a linear code and the actual size.
##
OptimalityLinearCode := function(C)
    local q, ub;
    q := Size(Field(C));
    ub := Minimum(UpperBound(WordLength(C), MinimumDistance(C), q),
                  UpperBoundGriesmer(WordLength(C), MinimumDistance(C), q));
    return q^LogInt(ub,q) - Size(C);
end;

#############################################################################
##
#F  BoundsOps . . . . . . .  operations containig functions for bounds record
##
BoundsOps := rec( name := "BoundsOps", operations := OpsOps );

BoundsOps.String := function(R)
    local line;
    line := Concatenation("an optimal linear [", String(R.n), ",",
                    String(R.k), ",d] code over GF(", String(R.q), ") has d");
    if R.upperBound <> R.lowerBound then
        Append(line,Concatenation(" in [", String(R.lowerBound),"..",
                String(R.upperBound),"]"));
    else
        Append(line,Concatenation("=",String(R.lowerBound)));
    fi;
    return line;
end;

BoundsOps.Display := function(R, D)
    local i, ref;
    Print(String(R),"\n");
    if IsBound(R.lowerBoundExplanation) then
        for i in [1..SizeScreen()[1]-2] do Print( "-" ); od; Print( "\n" );
        for i in R.lowerBoundExplanation do
            Print(i, "\n");
        od;
    fi;
    if IsBound(R.upperBoundExplanation) then
        for i in [1..SizeScreen()[1]-2] do Print( "-" ); od; Print( "\n" );
        for i in R.upperBoundExplanation do
            Print(i, "\n");
        od;
    fi;
    if IsBound(R.references) and Length(RecFields(R.references)) > 0 then
        for i in [1..SizeScreen()[1]-2] do  Print( "-" ); od; Print( "\n" );
        for i in RecFields(R.references) do
            Print("Reference ", i, ":\n");
            for ref in R.references.(i) do
                Print(ref, "\n");
            od;
        od;
    fi;
end;

BoundsOps.Print := function(R)
    Print(String(R));
end;

#############################################################################
##
#F  BoundsMinimumDistance( <n>, <k>, <F> )  . .  gets data from bounds tables
##
##  LowerBoundMinimumDistance uses (n, k, q, true)
##  LowerBoundMinimumDistance uses (n, k, q, false)
BoundsMinimumDistance := function(arg)
    local n, k, q, RecurseBound, res, InfoLine, GLOBAL_ALERT,
          DoTheTrick, kind, obj;

    InfoLine := function(n, k, d, S, spaces, prefix)
        local K, res;
        if kind = 1 then
            K := "L";
        else
            K := "U";
        fi;
        return String(Flat([ String(prefix, spaces), K, "b(", String(n), ","
                       , String(k), ")=", String(d), ", ", S ]));
    end;

    DoTheTrick := function( obj, man, str)
        if IsBound(obj.lastman) and obj.lastman = man then
            obj.expl[Length(obj.expl)] := str;
        else
            Add(obj.expl, str);
        fi;
        obj.lastman := man;
        return obj;
    end;
#F              RecurseBound
    RecurseBound := function(n, k, spaces, prefix)
        local i, obj, obj2;
        if k = 0 then  # Nullcode
            return rec(d := n, expl := [InfoLine(n, k, n,
                           "null code", spaces, prefix) ], cons :=
                       [NullCode, [n, q]]);
        elif k = 1 then  # RepetitionCode
            return rec(d := n, expl := [InfoLine(n, k, n,
                           "repetition code", spaces, prefix) ],
                       cons := [RepetitionCode, [n, q]]);
        elif k = 2 and q = 2 then  # Cordaro-Wagner code
            obj := rec( d :=2*Int((n+1)/3) - Int(n mod 3 / 2) );
            obj.expl :=  [InfoLine(n, k, obj.d, "Cordaro-Wagner code", spaces,
                                 prefix)];
            obj.cons := [CordaroWagnerCode,[n]];
            return obj;
        elif k = n-1 then  # Dual of the repetition code
            return rec( d :=2, 
                        expl := [InfoLine(n, k, 2,
                           "dual of the repetition code", spaces, prefix) ],
                        cons := [DualCode,[[RepetitionCode, [n, q]]]]);
        elif k = n then  # Whole space code
            return rec( d :=1,
                        expl := [InfoLine(n, k, 1, Concatenation(
                           "entire space GF(", String(q), ")^",
                                   String(n)), spaces, prefix) ],
                        cons := [WholeSpaceCode, [n, q]]);
        elif not IsBound(GUAVA_BOUNDS_TABLE[kind][q][n][k]) then
            if kind = 1 then  # trivial for lower bounds
                obj := rec(d :=2, 
                           expl := [InfoLine(n, k, 2,
                                   "expurgated dual of repetition code",
                                   spaces, prefix) ],
                           cons := [DualCode,[[RepetitionCode, [n, q]]]]);
                for i in [ k .. n - 2 ] do
                    obj.cons := [ExpurgatedCode,[obj.cons]];
                od;
                return obj;
            else  # Griesmer for upper bounds
                obj := rec( d := 2);
                while Sum([0..k-1], i -> 
                        QuoInt(obj.d, q^i) + SignInt(obj.d mod q^i)) <= n do
                    obj.d := obj.d + 1;
                od;
                obj.d := obj.d - 1;
                obj.expl := [InfoLine(n, k, obj.d, "Griesmer bound", spaces,
                                    prefix)];
                return obj;
            fi;
                       #Look up construction in table
        elif IsInt(GUAVA_BOUNDS_TABLE[kind][q][n][k]) then
            i := GUAVA_BOUNDS_TABLE[kind][q][n][k];
            if i = 1 then  # Shortening
                obj := RecurseBound(n+1, k+1, spaces, "");
                if IsBound(obj.lastman) and obj.lastman = 1 then
                    Add(obj.cons[2][2], Length(obj.cons[2][2])+1);
                else
                    obj.cons := [ShortenedCode, [ obj.cons, [1] ]];
                fi;
                return DoTheTrick( obj, 1, InfoLine(n, k, obj.d,
                               "by shortening of:", spaces, prefix) );
            elif i = 2 then  # Puncturing
                obj := RecurseBound(n+1, k, spaces, "");
                obj.d := obj.d - 1;
                if IsBound(obj.lastman) and obj.lastman = 2 then
                    Add(obj.cons[2][2], Length(obj.cons[2][2])+1);
                else
                    obj.cons := [ PuncturedCode, [ obj.cons, [1] ]];
                fi;
                return DoTheTrick( obj, 2, InfoLine(n, k, obj.d,
                               "by puncturing of:", spaces, prefix) );
            elif i = 3 then  # Extending
                obj := RecurseBound(n-1, k, spaces, "");
                if q = 2 and IsOddInt(obj.d) then
                    obj.d := obj.d + 1;
                fi;
                if IsBound(obj.lastman) and obj.lastman = 3 then
                    obj.cons[2][2] := obj.cons[2][2] + 1;
                else
                    obj.cons := [ ExtendedCode, [ obj.cons, 1 ]];
                fi;
                return DoTheTrick( obj, 3, InfoLine(n, k, obj.d,
                               "by extending:", spaces, prefix) );
            # Methods for upper bounds:
            elif i = 11 then  # Shortening
                obj := RecurseBound(n-1, k-1, spaces, "");
                return DoTheTrick( obj, 11, InfoLine(n, k, obj.d,
                               "otherwise shortening would contradict:",
                               spaces, prefix) );
            elif i = 12 then  # Puncturing
                obj := RecurseBound(n-1, k, spaces, "");
                obj.d := obj.d + 1;
                return DoTheTrick( obj, 12, InfoLine(n, k, obj.d,
                               "otherwise puncturing would contradict:",
                               spaces, prefix) );
            elif i = 13 then  #Extending
                obj := RecurseBound(n+1, k, spaces, "");
                if q=2 and IsOddInt(obj.d) then
                    obj.d := obj.d - 1;
                fi;
                return DoTheTrick( obj, 13, InfoLine(n, k, obj.d,
                               "otherwise extending would contradict:",
                               spaces, prefix) );
            else
                Error("invalid table entry; table is corrupted");
            fi;
        else
            i := GUAVA_BOUNDS_TABLE[kind][q][n][k];
            if i[1] = 0 then  # Code from library
                if IsBound( GUAVA_REF_LIST.(i[3]) ) then
                    res.references.(i[3]) := GUAVA_REF_LIST.(i[3]);
                else
                    res.references.(i[3]) := GUAVA_REF_LIST.ask;
                fi;
                obj := rec( d := i[2], expl := [InfoLine(n, k, i[2],
                               Concatenation("reference: ", i[3]),
                               spaces, prefix)], cons := false );
                if kind = 1 and not GLOBAL_ALERT then
                    GUAVA_TEMP_VAR := [n, k];
                    ReadPkg( "guava", "tbl", 
                            Concatenation( "codes", String(q) ) );
                    if GUAVA_TEMP_VAR = false then
                        GLOBAL_ALERT := true;
                    fi;
                    obj.cons := GUAVA_TEMP_VAR;
                fi;
                return obj;
            elif i[1] = 4 then  # Construction B
                obj := RecurseBound(n+i[2],k+i[2]-1, spaces, "");
                obj.cons := [ConstructionBCode, [obj.cons]];
                Add(obj.expl, InfoLine(n, k, obj.d, Concatenation(
                                    "by contruction B (deleting ",String(i[2]),
                                    " coordinates of a word in the dual)"),
                                    spaces, prefix) );
                Unbind(obj.lastman);
                return obj;
            elif i[1] = 5 then  # u | u+v construction
                obj := RecurseBound(n/2,   i[2], spaces + 4, "C1: ");
                obj2 :=RecurseBound(n/2, k-i[2], spaces + 4, "C2: ");
                obj.d := Minimum( 2 * obj.d, obj2.d );
                obj.cons := [UUVCode,[obj.cons, obj2.cons]];
                obj.expl := Concatenation(obj2.expl, obj.expl);
                Add(obj.expl, InfoLine(n, k, obj.d, 
                        "u|u+v construction of C1 and C2:", spaces, prefix));
                Unbind(obj.lastman);
                return obj;
            elif i[1] = 6 then  # Concatenation
                obj  := RecurseBound(n-i[2], k, spaces + 4, "C1: ");
                obj2 := RecurseBound(i[2],   k, spaces + 4, "C2: ");
                obj.cons := [ConcatenationCode,[obj.cons, obj2.cons]];
                obj.d := obj.d + obj2.d;
                obj.expl := Concatenation(obj2.expl, obj.expl);
                Add(obj.expl, InfoLine(n, k, obj.d, 
                        "concatenation of C1 and C2:", spaces, prefix));
                Unbind(obj.lastman);
                return obj;
            elif i[1] = 7 then  # ResidueCode
                obj := RecurseBound(i[2], k+1, spaces, "");
                obj.d := QuoInt(obj.d, q) + SignInt(obj.d mod q);
                Add(obj.expl, InfoLine(n, k, obj.d, "residue code of:", spaces,
                                    prefix) );
                obj.cons := [ResidueCode, [obj.cons]];
                Unbind(obj.lastman);
                return obj;
            elif i[1] = 14 then  # Construction B
                obj := RecurseBound(n-i[2], k-i[2]+1, spaces, "");
                Add(obj.expl, InfoLine(n, k, obj.d,
                        "otherwise construction B would contradict:", spaces,
                        prefix) );
                Unbind(obj.lastman);
                return obj;
            else
                Error("invalid table entry; table is corrupted");
            fi;
        fi;
    end;
#F              Function body
    if Length(arg) < 2 or Length(arg) > 4 then
        Error("usage: OptimalLinearCode( <n>, <k> [, <F>] )");
    fi;
    n := arg[1];
    k := arg[2];
    q := 2;
    if Length(arg) > 2 then
        if IsInt(arg[3]) then
            q := arg[3];
        else
            q := Size(arg[3]);
        fi;
    fi;
    if k > n then
        Error("k must be less than or equal to n");
    fi;
    # Check that right tables are present
    if not IsBound(GUAVA_REF_LIST) or Length(RecFields(GUAVA_REF_LIST))=0 then
        ReadPkg( "guava", "tbl", "refs" );
    fi;
    res := rec(n := n,
               k := k,
               q := q,
               references := rec(),
               construction := false,
               operations := BoundsOps);
    if not ( IsBound(GUAVA_BOUNDS_TABLE[1][q]) and
             IsBound(GUAVA_BOUNDS_TABLE[2][q]) ) and
       q > 4 then
        # Left the following lines out and replaced them with the previous,
        # and the else-part of this if, because using READ in
        # this way does not work in GAP 3.5.
        # (the behaviour of LOADED_PACKAGES has changed)
#       not READ(Concatenation(LOADED_PACKAGES.guava, "tbl/bdtable",
#               String(q),".g")) then       
        res.lowerBound := 1;
        res.upperBound := n - k + 1;
        return res;
#        Error("boundstable for q = ", q, " is not implemented.");
    else
        ReadPkg( "guava", "tbl", Concatenation( "bdtable", String(q) ) );
    fi;
    if n > Length(GUAVA_BOUNDS_TABLE[1][q]) then
        # no error should be returned here, otherwise Upper and
        # LowerBoundMinimumDistance would not work. The upper bound
        # could easely be sharpened by the Griesmer bound and
        # if n - k > Sz then
        #     upperbound >= n - k + 1;
        # else
        #     upperbound <= Ub[ Sz ][ k - n + Sz ]
        # fi;
        # lowerbound >= Lb[ Sz ][Minimum(Sz, k)]

        res.lowerBound := 1;
        res.upperBound := n - k + 1;
        return res;
#        Error("no data for n > ", Length(GUAVA_BOUNDS_TABLE[1][q]));
    fi;
    if Length(arg) < 4 or arg[4] then
        kind := 1;
        GLOBAL_ALERT := (Length(arg) = 4);
        obj := RecurseBound( n, k, 0, "");
        if not GLOBAL_ALERT then
            res.construction := obj.cons;
        fi;
        res.lowerBound := obj.d;
        res.lowerBoundExplanation := Reversed( obj.expl );
    fi;
    if Length(arg) < 4 or not arg[4] then
        kind := 2;
        obj := RecurseBound( n, k, 0, "");
        res.upperBound := obj.d;
        res.upperBoundExplanation := Reversed( obj.expl );
    fi;
    return res;
end;
#############################################################################
##
