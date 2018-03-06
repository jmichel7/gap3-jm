########################################################################
##
#F  BoundsCoveringRadius( <C> )
##
##  Find a lower and an upper bound for the covering radius of C.
##

BoundsCoveringRadius := function ( C )
    
    if not IsBound( C.boundsCoveringRadius ) then
        C.boundsCoveringRadius := C.operations.BoundsCoveringRadius( C );
    fi;
    
    return C.boundsCoveringRadius;
    
end;

CodeOps.BoundsCoveringRadius := function ( C )
    
    if IsLinearCode( C ) then
        return BoundsCoveringRadius( C );
    fi;
    
    return [ 1 .. WordLength( C ) ];
    
end;

LinCodeOps.BoundsCoveringRadius := function ( C )
    
    return [ GeneralLowerBoundCoveringRadius( C ) 
             .. GeneralUpperBoundCoveringRadius( C ) ];
    
end;

CycCodeOps.BoundsCoveringRadius := LinCodeOps.BoundsCoveringRadius;

########################################################################
##
#F  CoveringRadiusSearch( <C> )
##
##  Try to compute the covering radius. Don't compute all coset
##  leaders, but increment the lower bound as soon as a coset leader
##  is found.
##

CoveringRadiusSearch := function ( C )

    local k, n, i, j, lastone, zerofound, IsCosetLeader,
          lb, we, wd, vc, i, continue, codewords,
          leaderfound, allexamined, supp, elmsC, elms, len, one, zero,
          boundscr;
    
    IsCosetLeader := function( codewords, len, word, wt, one )
        
        local i, check, cw, wcw, j;
        
        check := true;
        i := 1;
        while i <= len and check do
            cw := codewords[ i ] + word;
            wcw := 0;
            for j in [ 1 .. Length( cw ) ] do
                if cw[ j ] = one then
                    wcw := wcw + 1;
                fi;
            od;
            if wcw < wt then
                check := false;
            fi;
            i := i + 1;
            
        od;
        
        return check;
        
    end;

    if not IsCode( C ) then
        Error( "CoveringRadiusSearch: <C> must be a code" );
    fi;
    if not IsLinearCode( C ) then
        Error( "CoveringRadiusSearch: <C> must be a linear code" );
    fi;
    if Size( Field( C ) ) <> 2 then
        Error( "CoveringRadiusSearch: <C> must be a binary code" );
    fi;
    
#    if IsBound( C.coveringRadius ) then
#        Return( C.coveringRadius );
#    fi;
            
    boundscr := BoundsCoveringRadius( C );
    lb := boundscr[ 1 ];
    n := WordLength( C );
    wd := WeightDistribution( C );
    #    elmsC := Elements( C );
    elms := [];
    for i in [ 0 .. n ] do
        if wd[ i + 1 ] > 0 then
            elms[ i + 1 ] := Elements(
                                     ConstantWeightSubcode( C, i ) );
#                                     elmsC, x -> WeightCodeword(x) = i );
        fi;
    od;
    
    for i in [ 1 .. n+1 ] do
        if IsBound(elms[ i ]) then
            for j in [ 1 .. Length( elms[ i ] ) ] do
                elms[ i ][ j ] := VectorCodeword( elms[ i ][ j ] );
            od;
        fi;
    od;
    
        

    # try to find a coset leader with weight > lb
    # if found, increase lb
    one := GF(2).one;
    zero := GF(2).zero;
    
    
    continue := true;
    while continue do
        k := C.boundsCoveringRadius[ 1 ] + 1;
        Print( "Trying ", k, " ...\n" );
        
        codewords := [ NullVector(n, GF(2) ) ];
        for i in [ 1 .. Minimum( n, 2 * k - 1) ] do
            if wd[ i + 1 ] <> 0 then
                Append( codewords, elms[ i + 1 ] );
            fi;
        od;
        
        len := Length( codewords );

#----------------
        
        vc := NullVector( n, GF(2) );
        for i in [ 1 .. k ] do
            vc[ i ] := one;
        od;
        lastone := k;
        allexamined := false;
        leaderfound := false;

        while not leaderfound and not allexamined do
            if not IsCosetLeader( codewords, len, vc, k, one ) then
                if lastone = n then
                    zerofound := false;
                    i := lastone - 1;
                    while i > n - k and vc[ i ] = one do
                        i := i - 1;
                    od;
                    if i = n - k then 
                        allexamined := true;
                    else
                        if i = n - k + 1 then
                            Print( Codeword(vc), "\n" );
                        fi;
                        j := i;
                        i := i + 1;
                        while vc[ j ] = zero do
                            j := j - 1;
                        od;
                        vc[ j ] := zero;
                        vc[ j + 1 ] := one;
                        j := j + 2;
                        if i <> j then
                            while i <= lastone do
                                vc[ j ] := one;
                                vc[ i ] := zero;
                                i := i + 1;
                                j := j + 1;
                            od;
                            lastone := j - 1;
                        else
                            lastone := n;
                        fi;
                    fi;
                else
                    vc[ lastone ] := zero;
                    lastone := lastone + 1;
                    vc[ lastone ] := one;
                fi;
            else
                leaderfound := true;
            fi;
        od;

        if leaderfound then
            C.boundsCoveringRadius := 
              Filtered( C.boundsCoveringRadius, x -> x >= k );
            continue := false;
        else
            C.boundsCoveringRadius :=
              [ C.boundsCoveringRadius[ 1 ] ];
            continue := false;
        fi;
        
#---------
    od;
    
    IsRange( C.boundsCoveringRadius );
    return( C.boundsCoveringRadius );
end;
        


#############################################################################
##
#F  CoveringRadius( <C> ) . . . . . . . . . . . .  the covering radius of <C>
##
##  Not useful for large codes.
##
##  It's now a bit more useful - Eric Minkes.
##
CoveringRadius := function ( C )
    
    if IsBound( C.operations.SpecialCoveringRadius ) then
        C.boundsCoveringRadius := C.operations.SpecialCoveringRadius( C );
    fi;
    
    if Length( BoundsCoveringRadius( C ) ) = 1 then
        return C.boundsCoveringRadius[ 1 ];
    else
        if not IsLinearCode( C ) then
            C.boundsCoveringRadius := [ C.operations.CoveringRadius( C ) ];
        else
            if Redundancy( C ) < 20 then
                C.boundsCoveringRadius := 
                  [ C.operations.CoveringRadius( C ) ];
            elif Redundancy( C ) < 30 and Dimension( C ) < 20 and 
              IsLinearCode( C ) and Size( Field( C ) ) = 2 then
                C.boundsCoveringRadius :=
                  [ CoveringRadiusSearch( C ) ];
            else
                Print("CoveringRadius: warning, the covering radius of \n",
                      "this code cannot be computed straightforward. \n",
                      "Try to use FindCoveringRadius( <C> ).\n",
                      "(see the manual for more details).\n" );
                Print("The covering radius of <C> lies in the interval:\n");
                return BoundsCoveringRadius( C );
            fi;
        fi;
    fi;
    return C.boundsCoveringRadius[ 1 ];
end;

CodeOps.CoveringRadius := function (C)
    local Code, vector, d, curmax, n, q, one, count, size, t, i, j,
    LocalDistance, zero, large, AddOneToVector;

    if IsLinearCode(C) then
        return CoveringRadius(C);    
    fi;
    q := Size(Field(C));
    n := WordLength(C);
    size := Size(C);
    one := Field(C).one;
    if q = 2 then
        Code := List(VectorCodeword(Elements(C)), i -> BlistList(i, [one]));
        vector := List([1..n], i-> false);
        zero := false;
        large := true;
        AddOneToVector := function(t) vector[t] := true; end;
        LocalDistance := DistanceBlist;
    else
        zero := Field(C).zero;
        if IsPrimeInt(q) then
            Code := IntVecFFE(VectorCodeword(Elements(C)));
        else
            Code := List(VectorCodeword(Elements(C)), i -> List(i, function(j) 
                if j = zero then
                    return 0;
                else return LogFFE(j,Z(q))+1; fi;
            end));
        fi;
        vector := List([1..n], i-> 0);
        zero := 0;
        large := q-1;
        AddOneToVector := function(t) vector[t] := vector[t] + 1; end;
        LocalDistance := function(a,b)
            j := 0;
            for i in [1..n] do
                if a[i] <> b[i] then j := j + 1; fi;
            od;
            return j;
        end;
    fi;
    curmax := n;
    for t in Code do
        d := LocalDistance(t, vector);
        if d < curmax then
            curmax := d;
        fi;
    od;
    for count in [2..q^n] do
        t := n;
        while vector[t] = large do
            vector[t] := zero;
            t := t - 1;
        od;
        AddOneToVector(t);
        t := 1;
        repeat
            d := LocalDistance(Code[t], vector);
            t := t + 1;
        until d <= curmax or t > size;
        if d > curmax then
            curmax := n;
            for t in Code do
                d := LocalDistance(t, vector);
                if d < curmax then
                    curmax := d;
                fi;
            od;
        fi;
    od;
    return curmax;
end;

LinCodeOps.CoveringRadius := function ( C )
    if Redundancy(C) = 0 then
        return 0;
    else
        return Maximum(List(SyndromeTable(C), i->WeightCodeword(i[1])));
    fi;
end;

CycCodeOps.CoveringRadius := LinCodeOps.CoveringRadius;

