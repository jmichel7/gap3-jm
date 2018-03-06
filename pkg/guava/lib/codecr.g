########################################################################
##
#F  CoveringRadius( <code> )
##
##  Return the covering radius of <code>
##  In case a special algorithm for this code exist, call
##  it first.
##
##  Not useful for large codes.
##
##  That's why I changed it, see the manual for more details
##  -- eric minkes.
##

CoveringRadius := function ( code )
    
    # call the special algorithm for this code, if it exists
    if IsBound( code.operations.SpecialCoveringRadius ) then
        code.boundsCoveringRadius :=
          code.operations.SpecialCoveringRadius( code );
    fi;
    
    if Length( BoundsCoveringRadius( code ) ) = 1 then
        return code.boundsCoveringRadius[ 1 ];
    else
        if not IsLinearCode( code ) then
            # for small codes in large spaces, this may take
            # a very long time
            code.boundsCoveringRadius :=
              [ code.operations.CoveringRadius( code ) ];
        else
            if Redundancy( code ) < 20 then
                code.boundsCoveringRadius :=
                  [ code.operations.CoveringRadius( code ) ];
            else
                InfoCoveringRadius(
                  "CoveringRadius: warning, the covering radius of \n",
                  "this code cannot be computed straightforward. \n",
                  "Try to use IncreaseCoveringRadiusLowerBound( <code> ).\n",
                  "(see the manual for more details).\n",
                  "The covering radius of <code> lies in the interval:\n" );
                return BoundsCoveringRadius( code );
            fi;
        fi;
    fi;
    return code.boundsCoveringRadius[ 1 ];
end;

# the original function of GUAVA, for unrestricted codes
CodeOps.CoveringRadius := function (code)
    local Code, vector, d, curmax, n, q, one, count, size, t, i, j,
    LocalDistance, zero, large, AddOneToVector;

    if IsLinearCode(code) then
        return CoveringRadius(code);
    fi;
    q := Size(Field(code));
    n := WordLength(code);
    size := Size(code);
    one := Field(code).one;
    if q = 2 then
        Code := List(VectorCodeword(Elements(code)), i -> BlistList(i, [one]));
        vector := List([1..n], i-> false);
        zero := false;
        large := true;
        AddOneToVector := function(t) vector[t] := true; end;
        LocalDistance := DistanceBlist;
    else
        zero := Field(code).zero;
        if IsPrimeInt(q) then
            Code := IntVecFFE(VectorCodeword(Elements(code)));
        else
            Code := List(VectorCodeword(Elements(code)), i -> List(i, function(j)
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

# the original function for linear codes
LinCodeOps.CoveringRadius := function ( code )
    if Redundancy(code) = 0 then
        return 0;
    else
        return Maximum( List( SyndromeTable( code ),
          i -> WeightCodeword( i[ 1 ] ) ) );
    fi;
end;

CycCodeOps.CoveringRadius := LinCodeOps.CoveringRadius;


########################################################################
##
#F  BoundsCoveringRadius( <code> )
##
##  Find a lower and an upper bound for the covering radius of code.
##

BoundsCoveringRadius := function ( code )
    
    if not IsBound( code.boundsCoveringRadius ) then
        code.boundsCoveringRadius :=
          [ GeneralLowerBoundCoveringRadius( code )
            .. GeneralUpperBoundCoveringRadius( code ) ];
        if Length(code.boundsCoveringRadius) = 0 then
            code.boundsCoveringRadius:=[0,WordLength(code)];
        fi;
    fi;
    
    return code.boundsCoveringRadius;
end;


########################################################################
##
#F  SetCoveringRadius( <code>, <cr> )
##  SetCoveringRadius( <code>, <interval> )
##
##  Enable the user to set the covering radius him/herself.
##

SetCoveringRadius := function ( arg )

    local code, cr;
    if Length( arg ) <> 2 then
        Error( "usage: SetCoveringRadius( <code>, <cr | interval> )" );
    fi;
    code := arg[ 1 ];
    if not IsCode( code ) then
        Error( "SetCoveringRadius: <code> must be a code" );
    fi;
    cr := arg[ 2 ];
    if IsInt( cr ) then
        if cr < 0 then
            Error( "SetCoveringRadius: <cr> must be non-negative" );
        fi;
    elif not IsVector( cr ) then
        Error( "SetCoveringRadius: <cr> must be an integer or a ",
               "list of integers" );
    fi;
    if IsInt( cr ) then
        code.boundsCoveringRadius := [ cr ];
    else
        code.boundsCoveringRadius := IntersectionSet(
            BoundsCoveringRadius( code ), cr );
        if Length( code.boundsCoveringRadius ) = 0 then
            code.boundsCoveringRadius := cr;
        fi;
        IsRange( code.boundsCoveringRadius );
    fi;
end;


########################################################################
##
#F  IncreaseCoveringRadiusLowerBound(
##      <code> [, <stopdistance> ] [, <startword> ] )
##

IncreaseCoveringRadiusLowerBound := function ( arg )
    
    local code,              # the code
          n,                 # length of the code
          k,                 # dimension of the code
          genmat,            # generator matrix of the code
          field,             # field of the code
          q,                 # the size of field
          fieldels,          # elements of field
          fieldzero,         # zero element of field
          nonzeroels,        # non-zero elements of field
          stopdistance,      # if the current word is a coset leader
                             # of this weight, then stop the algorithm
          boundscr,          # current bounds on the covering radius
          lb,                # current achieved lower bound
          current,           # the element of field^n that we are
                             # currently checking
          cwcurrent,         # current in Codeword form
          distcurrenttocode, # the distance of current to the code
          currentchanged,    # did we change current during the loop ?
          counterchanged,    # number of changes we have made so far
          countertotal,      # number of iterations we have made so far
          h, i, j,           # indexes to make the first slice
          slice,             # the number of the current slice
          numberofslices,    # the elements of the code are handled in 
                             # slices of 2^10 elements
          slicedim,          # the dimension of a slice
          slicesize,         # the size of a slice ( q^slicedim )
          index,             # to enumerate the codewords, the
                             # correct generator must be added
                             # index contains the generator number
          satisfied,         # are we satisfied with the results ?
          words,             # a list of codewords in the current slice
          wordslist0,        # a list of codewords with distance
                             # distcurrentocode + 0 from current
          wordslist1,        # a list of codewords with distance
                             # distcurrentocode + 1 from current
          word,              # a index for wordslist*
          coord,             # the coordinate where current will
                             # be changed
          staychance,        # chance that we stay at the same distance
          downchance,        # chance that we move closer to the code
          bestdist,          # best distance reached so far
          bestword,          # the corresponding word
          newelement,        # the new element for current[ coord ]
          newdist;           # the new distance of the changed current
                             # to the code

    # check the arguments
    if Length( arg ) < 1 or Length( arg ) > 3 then
        Error( "Usage: IncreaseCoveringRadiusLowerBound( ",
               "<code> [, <stopdistance> ] [, <startword> ] )" );
    fi;
    
    code := arg[ 1 ];
    if not IsCode( code ) then
        Error( "IncreaseCRLB: <code> must be a code" );
    fi;
    
    # extract some code parameters
    n := WordLength( code );
    if IsLinearCode( code ) then
        k := Dimension( code );
    fi;
    
    field := Field( code );
    q := Size( field );
    
    # do we want to stop at a certain weight ?
    if Length( arg ) >= 2 then
        if IsInt( arg[ 2 ] ) then
            stopdistance := arg[ 2 ];
        else
            stopdistance := -1;
        fi;
    else
        stopdistance := -1;     # default: never stop, unless
                                # the lower bound meets the upper bound
    fi;
    
    # if we cannot compute the minimum distance of the code,
    # this algorithm will not be of much help either.
    # <lb> is a safe place to start searching
    lb := IntFloor( MinimumDistance( code ) / 2 );
    
    boundscr := BoundsCoveringRadius( code );
    if lb > boundscr[ 1 ] then
        code.boundsCoveringRadius := Filtered( boundscr, n -> lb <= n );
        IsRange( code.boundsCoveringRadius );
        boundscr := code.boundsCoveringRadius;
    fi;
    
    if Length( boundscr ) = 1 then
        # if there is nothing to compute, 
        # then just return the covering radius
        return boundscr[ 1 ];     
    fi;
    
    # make a vector of distance lb or more away
    if ( stopdistance = -1 and Length( arg ) = 2 ) or
         Length( arg ) = 3 then
        current := VectorCodeword( Codeword( arg[ Length( arg ) ] ) );
        distcurrenttocode := MinimumDistance( code, current );
    else
        current := RandomVector( n, RandomList( [0..n] ), field );
        distcurrenttocode := MinimumDistance( code, current );
        while distcurrenttocode < lb do
            current := RandomVector( n, RandomList( [0..n] ), field );
            distcurrenttocode := MinimumDistance( code, current );
        od;
    fi;
    bestdist := distcurrenttocode;
    bestword := Copy( current );
    
    # initialise some parameters and useful variables
    fieldels := Elements( field );
    fieldzero := field.zero;
    nonzeroels := Difference( fieldels, [ fieldzero ] );
    
    if IsLinearCode( code ) then
        genmat := GeneratorMat( code );
        # try to make the size of a group codewords to be about
        # CRMemSize
        slicedim := LogInt( CRMemSize, q );
        numberofslices := Maximum( 1, q^( k-slicedim ) );
        slicesize := q^slicedim;
    fi;
    
    satisfied := false;
    currentchanged := true;
    counterchanged := 0;
    countertotal := 0;
    staychance := 10; # maybe make arguments of these
    downchance := 1;
    
    # the words array contains all codewords generated by
    # genmat[ 1 ] ... genmat[ slicedim ] ( or genmat[ k ], if
    # slicedim > k )
    
    if IsLinearCode( code ) then
        words := [ NullVector( n, field ) ];
        for i in [ 1 .. Minimum( slicedim, k ) ] do
            for j in [ 1 .. q^( i-1 ) ] do
                for h in [ 1 .. q-1 ] do
                    Add( words, words[ j ] + genmat[ i ] * nonzeroels[ h ] );
                od;
            od;
        od;
    else
        # if the code is non-linear, the
        # elements field is obligatory,
        # so it fits in memory.
        # no need for all the hassle as in the linear case,
        # but the disadvantage is that we can't handle
        # large non-linear codes.
        # but then again, GUAVA can't handle these anyway.
        words := VectorCodeword( Elements( code ) );
    fi;

    # start algorithm, stop when we are satisfied with the results
    # (either we have a coset leader of weigth <stopweigth>,
    # or lb = ub)
    while not satisfied do
        countertotal := countertotal + 1;
        if countertotal mod 1000 = 0 then
            InfoCoveringRadius( "Number of runs: ",
                    countertotal, 
                    "  best distance so far: ", bestdist, "\n" );
        fi;
        if currentchanged then
            # current word has changed, generate three lists of
            # codewords that have distance distcurrenttocode,
            # distcurrenttocode + 1 and distcurrenttocode + 2
            
            cwcurrent := Codeword( current );
            
            if distcurrenttocode > bestdist then
                bestdist := distcurrenttocode;
                bestword := Copy( current );
                InfoCoveringRadius( 
                        "New best distance: ", bestdist, "\n" );
            fi;
                
            wordslist0 := [];
            wordslist1 := [];

            if IsLinearCode( code ) then
                for slice in [ 0 .. numberofslices-1 ] do
                    if slice > 0 then
                        i := k - slicedim - 1;
                        while EuclideanRemainder( slice, q^i ) <> 0 do
                            i := i - 1;
                        od;
                        index := slicedim + i + 1;
                        words := List( words, x -> x + genmat[ index ] );
                    fi;
                    
                    Append( wordslist0, 
                            Filtered( words, x -> 
                                    DistanceVecFFE( x, current ) 
                                    = distcurrenttocode ) );
                    Append( wordslist1,
                            Filtered( words, x ->
                                    DistanceVecFFE( x, current ) 
                                    = distcurrenttocode + 1 ) );
                    
                od;
            else
                wordslist0 := Filtered( words, x->
                                      DistanceVecFFE( x, current )
                                      = distcurrenttocode );
                wordslist1 := Filtered( words, x->
                                      DistanceVecFFE( x, current )
                                      = distcurrenttocode + 1 );
            fi;
            
            currentchanged := false;
            counterchanged := counterchanged + 1;
            if EuclideanRemainder( counterchanged, 100 ) = 0 then
                InfoCoveringRadius( "Number of changes: ", 
                        counterchanged, "\n" );
            fi;
        fi;

        # pick a coordinate
        # the algorithm will look what happens if we change this coordinate
        coord := RandomList( [ 1 .. n ] );
        
        # the possible new element for this coordinate
        # is picked at random from the field elements
        newelement := RandomList( Difference( fieldels,
                              [ current[ coord ] ] ) );
        
        # check the new word against the codewords that are
        # at distance distcurrenttocode + 0 from current
        # the result can be: 
        #   1) the distance to all words in wordslist0 is 
        #      one more than distcurrenttocode 
        #      (this is the situation that we hope for)
        #   2) there is at least one word in wordslist0 that
        #      has distance distcurrenttocode - 1 to the
        #      new current
        #   3) if the field is not GF(2):
        #      there is a word in wordslist0 that stays
        #      at the same distance
        newdist := distcurrenttocode + 1;
        for word in wordslist0 do
            if word[ coord ] <> current[ coord ] then
                if word[ coord ] = newelement then
                    newdist := distcurrenttocode - 1;
                else
                    newdist := distcurrenttocode;
                fi;
            fi;
        od;
        
        # only check against other words if the previous tests
        # did not fail
        if newdist > distcurrenttocode then
            
            # check the new word against the codewords that are at 
            # distance distcurrenttocode + 1 from current
            # again, two results are possible
            #   1) all words in wordslist1 are at distance 
            #      distcurrenttocode + 1 or + 2 (this is good)
            #   2) there is a word in wordslist1 that now has
            #      distance distcurrenttocode to current
            #      this means we did not find an improvement
            
            for word in wordslist1 do
                if word[ coord ] = newelement then
                    newdist := distcurrenttocode;
                fi;
            od;
        fi;
        
        if newdist > distcurrenttocode then
            # we found a new coset leader with larger weight
            
            # now change current 
            current[ coord ] := newelement;
            currentchanged := true;
            distcurrenttocode := newdist;
            
            # also check whether the covering radius lower bound
            # can be increased, this is what the whole
            # algorithm is about !
            
            if distcurrenttocode > boundscr[ 1 ] then
                
                # write directly to the code to make the change
                # permanent, even if the user interrupted us
                code.boundsCoveringRadius :=
                  Filtered( boundscr, x -> x >= distcurrenttocode );
                # make it a range together if possible
                IsRange( code.boundsCoveringRadius );
                boundscr := code.boundsCoveringRadius;

                # maybe we have reached the upper bound
                # then we can stop altogether !
                if Length( boundscr ) = 1 then
                    satisfied := true;
                fi;
            fi;
        elif newdist = distcurrenttocode then
            # the change to the word did not change the
            # distance to the code
            if RandomList( [ 1 .. 100 ] ) <= staychance then
                current[ coord ] := newelement;
                currentchanged := true;
            fi;
        else
            # make it a 1 in 100 chance to get closer anyway
            # because we do not want to get stuck in a
            # suboptimal coset
            if RandomList( [ 1 .. 100 ] ) <= downchance then
                current[ coord ] := newelement;
                currentchanged := true;
                distcurrenttocode := newdist;
            fi;
        fi;
        
        # maybe the distance of current to the code 
        # is high enough for the user
        # then we should stop
        if distcurrenttocode = stopdistance then
            satisfied := true;
        fi;
    od;
    
    # return the new covering radius bounds, and a coset leader
    # that has weight equal to the lower bound
    return rec( boundsCoveringRadius := code.boundsCoveringRadius,
                cosetLeader := Codeword( current ) );
end;


########################################################################
##
#F  ExhaustiveSearchCoveringRadius( <code> )
##
##  Try to compute the covering radius. Don't compute all coset
##  leaders, but increment the lower bound as soon as a coset leader
##  is found.
##

ExhaustiveSearchCoveringRadius := function ( arg )

    local k, n, i, j, lastone, zerofound, IsCosetLeader,
          lb, we, wd, vc, i, continue, codewords,
          leaderfound, allexamined, supp, elmsC, elms, len, one, zero,
          boundscr, code, stopsoon;
    
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

    if Length( arg ) < 1 or Length( arg ) > 2 then
        Error(
          "usage: ExhaustiveSearchCoveringRadius( <code> [, <quick> ] )" );
    fi;
    code := arg[ 1 ];
    if not IsCode( code ) then
        Error( "CoveringRadiusSearch: <code> must be a code" );
    fi;
    if not IsLinearCode( code ) then
        Error( "CoveringRadiusSearch: <code> must be a linear code" );
    fi;
    if Size( Field( code ) ) <> 2 then
        Error( "CoveringRadiusSearch: <code> must be a binary code" );
    fi;
    if Length( arg ) = 2 then
        stopsoon := arg[ 2 ];
        if stopsoon <> false then
            stopsoon := true;
        fi;
    else
        stopsoon := true;
    fi;

    boundscr := BoundsCoveringRadius( code );
    if Length( boundscr ) = 1 then
        return boundscr[ 1 ];
    fi;
    
    lb := boundscr[ 1 ];
    n := WordLength( code );
    wd := WeightDistribution( code );
    elms := [];
    for i in [ 0 .. n ] do
        if wd[ i + 1 ] > 0 then
            elms[ i + 1 ] := Elements(
                                     ConstantWeightSubcode( code, i ) );
        fi;
    od;
    
    for i in [ 1 .. n+1 ] do
        if IsBound( elms[ i ] ) then
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
        k := code.boundsCoveringRadius[ 1 ] + 1;
        InfoCoveringRadius( "Trying ", k, " ...\n" );
        codewords := [ NullVector(n, GF(2) ) ];
        for i in [ 1 .. Minimum( n, 2 * k - 1) ] do
            if wd[ i + 1 ] <> 0 then
                Append( codewords, elms[ i + 1 ] );
            fi;
        od;
        len := Length( codewords );

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
            code.boundsCoveringRadius :=
              Filtered( code.boundsCoveringRadius, x -> x >= k );
            if stopsoon then
                continue := false;
            fi;
        else
            code.boundsCoveringRadius :=
              [ code.boundsCoveringRadius[ 1 ] ];
            continue := false;
        fi;
    od;
    
    IsRange( code.boundsCoveringRadius );
    return( code.boundsCoveringRadius );
end;


########################################################################
##
#F  CoveringRadiusLowerBoundTable
##

CoveringRadiusLowerBoundTable := [
    [ 3, 2, , , , , , , , ,                     # n = 13
       ,  , , , , , , , , ,
       ,  , , , , , , , , ,
       ,  , , , , , , , , ,
       ,  , , , , , , , ],
    [ 3, 3, , , , , , , , ,                     # n = 14
       ,  , , , , , , , , ,
       ,  , , , , , , , , ,
       ,  , , , , , , , , ,
       ,  , , , , , , , ],
    [ 4, 3, 3, , , , , , , ,                    # n = 15
       ,  ,  , , , , , , , ,
       ,  ,  , , , , , , , ,
       ,  ,  , , , , , , , ,
       ,  ,  , , , , , , ],
    [ 4, 4, 3, 3, , , , , , ,                   # n = 16
       ,  ,  ,  , , , , , , ,
       ,  ,  ,  , , , , , , ,
       ,  ,  ,  , , , , , , ,
       ,  ,  ,  , , , , , ],
    [ 4, 4, 3, 3, 3, , , , , ,                  # n = 17
       ,  ,  ,  ,  , , , , , ,
       ,  ,  ,  ,  , , , , , ,
       ,  ,  ,  ,  , , , , , ,
       ,  ,  ,  ,  , , , , ],
    [ 5, 4, 4, 3, 3, 3, , , , ,                 # n = 18
       ,  ,  ,  ,  ,  , , , , ,
       ,  ,  ,  ,  ,  , , , , ,
       ,  ,  ,  ,  ,  , , , , ,
       ,  ,  ,  ,  ,  , , , ],
    [ 5, 4, 4, 4, 3, 3, 2, , , ,                # n = 19
       ,  ,  ,  ,  ,  ,  , , , ,
       ,  ,  ,  ,  ,  ,  , , , ,
       ,  ,  ,  ,  ,  ,  , , , ,
       ,  ,  ,  ,  ,  ,  , , ],
    [ 6, 5, 4, 4, 4, 3, 3, 2, , ,               # n = 20
       ,  ,  ,  ,  ,  ,  ,  , , ,
       ,  ,  ,  ,  ,  ,  ,  , , ,
       ,  ,  ,  ,  ,  ,  ,  , , ,
       ,  ,  ,  ,  ,  ,  ,  , ],
    [ 6, 5, 5, 4, 4, 3, 3, 3, , ,               # n = 21
       ,  ,  ,  ,  ,  ,  ,  , , ,
       ,  ,  ,  ,  ,  ,  ,  , , ,
       ,  ,  ,  ,  ,  ,  ,  , , ,
       ,  ,  ,  ,  ,  ,  ,  , ],
    [ 6, 6, 5, 5, 4, 4, 3, 3, 3, ,              # n = 22
       ,  ,  ,  ,  ,  ,  ,  ,  , ,
       ,  ,  ,  ,  ,  ,  ,  ,  , ,
       ,  ,  ,  ,  ,  ,  ,  ,  , ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 7, 6, 5, 5, 4, 4, 3, 3, 3, 3,             # n = 23
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 7, 6, 6, 5, 5, 4, 4, 3, 3, 3,             # n = 24
      3,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 7, 7, 6, 6, 5, 5, 4, 4, 3, 3,             # n = 25
      3, 2,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 8, 7, 7, 6, 6, 5, 5, 4, 4, 3,             # n = 26
      3, 3, 2,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 8, 8, 7, 6, 6, 5, 5, 4, 4, 4,             # n = 27
      3, 3, 3, 2,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 9, 8, 7, 7, 6, 6, 5, 5, 4, 4,             # n = 28
      4, 3, 3, 3, 2,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 9, 8, 8, 7, 7, 6, 6, 5, 5, 4,             # n = 29
      4, 4, 3, 3, 3, 2,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 9, 9, 8, 7, 7, 6, 6, 5, 5, 5,             # n = 30
      4, 4, 4, 3, 3, 3,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
       ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 10, 9, 8, 8, 7, 7, 6, 6, 5, 5,            # n = 31
       5, 4, 4, 3, 3, 3, 3,  ,  ,  ,
        ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
        ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
        ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 10, 9, 9, 8, 8, 7, 7, 6, 6, 5,            # n = 32
       5, 4, 4, 4, 3, 3, 3, 3,  ,  ,
        ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
        ,  ,  ,  ,  ,  ,  ,  ,  ,  ,
        ,  ,  ,  ,  ,  ,  ,  ,  ],
    [ 11, 10, 9, 9, 8, 7, 7, 6, 6, 6,           # n = 33
       5,  5, 4, 4, 4, 3, 3, 3, 3,  ,
        ,   ,  ,  ,  ,  ,  ,  ,  ,  ,
        ,   ,  ,  ,  ,  ,  ,  ,  ,  ,
        ,   ,  ,  ,  ,  ,  ,  ,  ],
    [ 11, 10, 10, 9, 8, 8, 7, 7, 6, 6,          # n = 34
       5,  5,  5, 4, 4, 4, 3, 3, 3, 2,
        ,   ,   ,  ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,  ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,  ,  ,  ,  ,  ,  ],
    [ 11, 11, 10, 9, 9, 8, 8, 7, 7, 6,          # n = 35
       6,  5,  5, 5, 4, 4, 4, 3, 3, 3,
       2,   ,   ,  ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,  ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,  ,  ,  ,  ,  ,  ],
    [ 12, 11, 10, 10, 9, 9, 8, 8, 7, 7,         # n = 36
       6,  6,  5,  5, 5, 4, 4, 4, 3, 3,
       3,  2,   ,   ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,   ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,   ,  ,  ,  ,  ,  ],
    [ 12, 11, 11, 10, 9, 9, 8, 8, 7, 7,         # n = 37
       6,  6,  6,  5, 5, 4, 4, 4, 4, 3,
       3,  3,  2,   ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,   ,  ,  ,  ,  ,  ,  ,
        ,   ,   ,   ,  ,  ,  ,  ,  ],
    [ 13, 12, 11, 10, 10, 9, 9, 8, 8, 7,        # n = 38
       7,  6,  6,  6,  5, 5, 4, 4, 4, 3,
       3,  3,  3,  2,   ,  ,  ,  ,  ,  ,
        ,   ,   ,   ,   ,  ,  ,  ,  ,  ,
        ,   ,   ,   ,   ,  ,  ,  ,  ],
    [ 13, 12, 12, 11, 10, 10, 9, 9, 8, 8,       # n = 39
       7,  7,  6,  6,  5,  5, 5, 4, 4, 4,
       3,  3,  3,  3,  2,   ,  ,  ,  ,  ,
        ,   ,   ,   ,   ,   ,  ,  ,  ,  ,
        ,   ,   ,   ,   ,   ,  ,  ,  ],
    [ 14, 13, 12, 11, 11, 10, 9, 9, 8, 8,       # n = 40
       7,  7,  7,  6,  6,  5, 5, 5, 4, 4,
       4,  3,  3,  3,  3,  2,  ,  ,  ,  ,
        ,   ,   ,   ,   ,   ,  ,  ,  ,  ,
        ,   ,   ,   ,   ,   ,  ,  ,  ],
    [ 14, 13, 12, 12, 11, 10, 10, 9, 9, 8,      # n = 41
       8,  7,  7,  6,  6,  6,  5, 5, 5, 4,
       4,  4,  3,  3,  3,  3,  2,  ,  ,  ,
        ,   ,   ,   ,   ,   ,   ,  ,  ,  ,
        ,   ,   ,   ,   ,   ,   ,  ,  ],
    [ 14, 14, 13, 12, 11, 11, 10, 10, 9, 9,     # n = 42
       8,  8,  7,  7,  6,  6,  6,  5, 5, 5,
       4,  4,  4,  3,  3,  3,  3,  2,  ,  ,
        ,   ,   ,   ,   ,   ,   ,   ,  ,  ,
        ,   ,   ,   ,   ,   ,   ,   ,  ],
    [ 15, 14, 13, 12, 12, 11, 10, 10, 9, 9,     # n = 43
       8,  8,  7,  7,  7,  6,  6,  6, 5, 5,
       5,  4,  4,  4,  3,  3,  3,  3, 2,  ,
        ,   ,   ,   ,   ,   ,   ,   ,  ,  ,
        ,   ,   ,   ,   ,   ,   ,   ,  ],
    [ 15, 14, 14, 13, 12, 11, 11, 10, 10, 9,    # n = 44
       9,  8,  8,  7,  7,  7,  6,  6,  5, 5,
       5,  4,  4,  4,  4,  3,  3,  3,  3,  ,
        ,   ,   ,   ,   ,   ,   ,   ,   ,  ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 16, 15, 14, 13, 12, 12, 11, 11, 10, 10,   # n = 45
       9,  9,  8,  8,  7,  7,  7,  6,  6,  5,
       5,  5,  4,  4,  4,  4,  3,  3,  3,  3,
        ,   ,   ,   ,   ,   ,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 16, 15, 14, 14, 13, 12, 12, 11, 11, 10,   # n = 46
       9,  9,  8,  8,  8,  7,  7,  6,  6,  6,
       5,  5,  5,  4,  4,  4,  4,  3,  3,  3,
       3,   ,   ,   ,   ,   ,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 16, 16, 15, 14, 13, 13, 12, 11, 11, 10,   # n = 47
      10,  9,  9,  8,  8,  8,  7,  7,  6,  6,
       6,  5,  5,  5,  4,  4,  4,  3,  3,  3,
       3,  2,   ,   ,   ,   ,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 17, 16, 15, 14, 14, 13, 12, 12, 11, 11,   # n = 48
      10, 10,  9,  9,  8,  8,  7,  7,  7,  6,
       6,  6,  5,  5,  5,  4,  4,  4,  3,  3,
       3,  3,  2,   ,   ,   ,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 17, 16, 16, 15, 14, 13, 13, 12, 12, 11,   # n = 49
      11, 10,  9,  9,  9,  8,  8,  7,  7,  7,
       6,  6,  6,  5,  5,  5,  4,  4,  4,  3,
       3,  3,  3,  2,   ,   ,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 18, 17, 16, 15, 14, 14, 13, 13, 12, 11,   # n = 50
      11, 10, 10,  9,  9,  8,  8,  8,  7,  7,
       7,  6,  6,  6,  5,  5,  5,  4,  4,  4,
       3,  3,  3,  3,  2,   ,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 18, 17, 16, 16, 15, 14, 13, 13, 12, 12,   # n = 51
      11, 11, 10, 10,  9,  9,  8,  8,  8,  7,
       7,  6,  6,  6,  5,  5,  5,  5,  4,  4,
       4,  3,  3,  3,  3,  2,   ,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 19, 18, 17, 16, 15, 15, 14, 13, 13, 12,   # n = 52
      12, 11, 11, 10, 10,  9,  9,  8,  8,  8,
       7,  7,  6,  6,  6,  5,  5,  5,  5,  4,
       4,  4,  3,  3,  3,  3,  2,   ,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 19, 18, 17, 16, 16, 15, 14, 14, 13, 12,   # n = 53
      12, 11, 11, 10, 10,  9,  9,  9,  8,  8,
       7,  7,  7,  6,  6,  6,  5,  5,  5,  4,
       4,  4,  4,  3,  3,  3,  3,  2,   ,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 19, 18, 18, 17, 16, 15, 15, 14, 13, 13,   # n = 54
      12, 12, 11, 11, 10, 10,  9,  9,  9,  8,
       8,  7,  7,  7,  6,  6,  6,  5,  5,  5,
       4,  4,  4,  4,  3,  3,  3,  3,  2,   ,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 20, 19, 18, 17, 16, 16, 15, 14, 14, 13,   # n = 55
      13, 12, 12, 11, 11, 10, 10,  9,  9,  8,
       8,  8,  7,  7,  7,  6,  6,  6,  5,  5,
       5,  4,  4,  4,  4,  3,  3,  3,  3,  2,
        ,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 20, 19, 18, 18, 17, 16, 15, 15, 14, 14,   # n = 56
      13, 12, 12, 11, 11, 10, 10, 10,  9,  9,
       8,  8,  8,  7,  7,  7,  6,  6,  6,  5,
       5,  5,  4,  4,  4,  4,  3,  3,  3,  3,
       2,   ,   ,   ,   ,   ,   ,   ,   ],
    [ 21, 20, 19, 18, 17, 16, 16, 15, 14, 14,   # n = 57
      13, 13, 12, 12, 11, 11, 10, 10,  9,  9,
       9,  8,  8,  7,  7,  7,  6,  6,  6,  5,
       5,  5,  5,  4,  4,  4,  4,  3,  3,  3,
       3,  2,   ,   ,   ,   ,   ,   ,   ],
    [ 21, 20, 19, 18, 18, 17, 16, 15, 15, 14,   # n = 58
      14, 13, 13, 12, 12, 11, 11, 10, 10,  9,
       9,  9,  8,  8,  7,  7,  7,  6,  6,  6,
       5,  5,  5,  5,  4,  4,  4,  4,  3,  3,
       3,  3,  2,   ,   ,   ,   ,   ,   ],
    [ 22, 21, 20, 19, 18, 17, 17, 16, 15, 15,   # n = 59
      14, 14, 13, 12, 12, 11, 11, 11, 10, 10,
       9,  9,  8,  8,  8,  7,  7,  7,  6,  6,
       6,  5,  5,  5,  5,  4,  4,  4,  4,  3,
       3,  3,  3,  2,   ,   ,   ,   ,   ],
    [ 22, 21, 20, 19, 18, 18, 17, 16, 16, 15,   # n = 60
      14, 14, 13, 13, 12, 12, 11, 11, 10, 10,
      10,  9,  9,  8,  8,  8,  7,  7,  7,  6,
       6,  6,  5,  5,  5,  5,  4,  4,  4,  3,
       3,  3,  3,  3,  2,   ,   ,   ,   ],
    [ 23, 21, 20, 20, 19, 18, 17, 17, 16, 15,   # n = 61
      15, 14, 14, 13, 13, 12, 12, 11, 11, 10,
      10, 10,  9,  9,  8,  8,  8,  7,  7,  7,
       6,  6,  6,  5,  5,  5,  5,  4,  4,  4,
       3,  3,  3,  3,  3,  2,   ,   ,   ],
    [ 23, 22, 21, 20, 19, 18, 18, 17, 16, 16,   # n = 62
      15, 15, 14, 14, 13, 12, 12, 12, 11, 11,
      10, 10,  9,  9,  9,  8,  8,  8,  7,  7,
       7,  6,  6,  6,  5,  5,  5,  4,  4,  4,
       4,  3,  3,  3,  3,  3,  2,   ,   ],
    [ 23, 22, 21, 20, 20, 19, 18, 17, 17, 16,   # n = 63
      16, 15, 14, 14, 13, 13, 12, 12, 11, 11,
      11, 10, 10,  9,  9,  9,  8,  8,  7,  7,
       7,  6,  6,  6,  6,  5,  5,  5,  4,  4,
       4,  4,  3,  3,  3,  3,  3,  2,   ],
    [ 24, 23, 22, 21, 20, 19, 18, 18, 17, 16,   # n = 64
      16, 15, 15, 14, 14, 13, 13, 12, 12, 11,
      11, 10, 10, 10,  9,  9,  8,  8,  8,  7,
       7,  7,  6,  6,  6,  6,  5,  5,  5,  4,
       4,  4,  4,  3,  3,  3,  3,  3,  2 ]
];

########################################################################
##
#F  GeneralLowerBoundCoveringRadius( <n>, <size> [, <F> ] )
##  GeneralLowerBoundCoveringRadius( <code> )
##

GeneralLowerBoundCoveringRadius := function ( arg )

    local code, n, size, field, l, glbcr, GLBCRforCodes,
          GLBCRforParameters;

    GLBCRforCodes := function( code )
        local n, size, field, listofbounds, k;
        n := WordLength( code );
        size := Size( code );
        field := Field( code );
        if not IsLinearCode( code ) then
            listofbounds := [
              LowerBoundCoveringRadiusSphereCovering( n, size, field, false ),
            ];
            if field = GF(2) then
                Append( listofbounds, [
                  LowerBoundCoveringRadiusVanWee1( n, size, field, false ),
                  LowerBoundCoveringRadiusVanWee2( n, size, false ),
                  LowerBoundCoveringRadiusCountingExcess( n, size, false )
                ] );
                if n <= 200 then
                    Append( listofbounds, [
                      LowerBoundCoveringRadiusEmbedded1( n, size, field, false ),
                      LowerBoundCoveringRadiusEmbedded2( n, size, field, false )
                    ] );
                fi;
            fi;
            return Maximum( listofbounds );
        elif field = GF(2) then
            k := Dimension( code );
            # small codimensions (n-k)
            if k = n then
                return 0;
            elif k = n - 1 then
                return 1;
            elif k = n - 2 then
                return 1;
            elif k = n - 3 and n >= 7 then
                return 1;
            elif k = n - 4 and n >= 5 then
                if n >= 15 then
                    return 1;
                else
                    return 2;
                fi;
            elif k = n - 5 and n >= 9 then
                if n >= 31 then
                    return 1;
                else
                    return 2;
                fi;
            elif k = n - 6 then
                if n >= 63 then
                    return 1;
                elif n >= 14 then
                    return 2;
                elif n = 12 then
                    return 3;
                fi;
            elif k = n - 7 and n >= 21 and n <= 64 then
                return 2;
            elif k = n - 8 and n >= 30 and n <= 64 then
                return 2;
            elif k = n - 9 and n >= 44 and n <= 64 then
                return 2;
            elif k = 1 then
                return IntFloor( n/2 );
            elif k = 2 then
                return IntFloor( (n-1)/2 );
            elif k = 3 then
                return IntFloor( (n-2)/2 );
            elif k = 4 then
                return IntFloor( (n-4)/2 );
            elif k = 5 then
                return IntFloor( (n-5)/2 );
            fi;

            # use the table of Cohen, Litsyn, Lobstein and Mattson
            if n >= 13 and n <= 64 and k>=6 and k<= 64 then
                if IsBound( 
                     CoveringRadiusLowerBoundTable[ n - 12 ][ k - 5 ] ) then
                    return CoveringRadiusLowerBoundTable[ n - 12 ][ k - 5 ];
                fi;
            fi;
            
            # not in the table. use the bounds
            listofbounds := [
              LowerBoundCoveringRadiusSphereCovering( n, size, field, false ),
              LowerBoundCoveringRadiusVanWee1( n, size, field, false ),
              LowerBoundCoveringRadiusVanWee2( n, size, false ),
              LowerBoundCoveringRadiusCountingExcess( n, size, false )
            ];
            if n <= 200 then
                Append( listofbounds, [
                  LowerBoundCoveringRadiusEmbedded1( n, size, field, false ),
                  LowerBoundCoveringRadiusEmbedded2( n, size, field, false ),
                ] );        
            fi;
            return Maximum( listofbounds );
        else  # field is not GF(2)
            listofbounds := [
              LowerBoundCoveringRadiusSphereCovering( n, size, field, false ),
            ];
            return Maximum( listofbounds );
        fi;
    end;

    GLBCRforParameters := function( n, size, field )
        local listofbounds;
        listofbounds := [
          LowerBoundCoveringRadiusSphereCovering( n, size, field, false ),
          LowerBoundCoveringRadiusVanWee1( n, size, field, false ),
          LowerBoundCoveringRadiusEmbedded1( n, size, field, false ),
          LowerBoundCoveringRadiusEmbedded2( n, size, field, false )
        ];
        if field = GF(2) then
            Append( listofbounds, [
              LowerBoundCoveringRadiusVanWee2( n, size, false ),
              LowerBoundCoveringRadiusCountingExcess( n, size, false )
            ] );
        fi;
        return Maximum( listofbounds );
    end;
    
    l := Length( arg );
    if l < 1 or l > 3 then
        Error( "usage: GeneralLBCR( <n>, <size> [, <F> ] | <code> )" );
    fi;
    if l = 1 then
        if IsCode( arg[ 1 ] ) then
            code := arg[ 1 ];
            glbcr := GLBCRforCodes( arg[ 1 ] );
            if IsBound( code.boundsCoveringRadius ) then
                glbcr := Maximum( [ glbcr,
                  code.boundsCoveringRadius[ 1 ] ] );
            fi;
            return glbcr;
        else
            Error( "GeneralLBCR: <code> must be a code" );
        fi;
    else
        n := arg[ 1 ];
        size := arg[ 2 ];
        if Length( arg ) = 3 then
            field := arg[ 3 ];
        else
            field := GF( 2 );
        fi;
        if not IsInt( n ) or
           not IsInt( size ) or
           not IsField( field ) then
            Error( "GeneralLBCR: <n>, <size> must be integers, ",
                   "<field> a field" );
        fi;
        if n < 1 or size < 1 then
            Error( "GeneralLBCR: <n> and <size> must be positive" );
        fi;
        glbcr := GLBCRforParameters( n, size, field );
        if IsBound( code.boundsCoveringRadius ) then
            glbcr := Maximum( [ glbcr,
              code.boundsCoveringRadius[ 1 ] ] );
        fi;
        return glbcr;
    fi;
end;


########################################################################
##
#F  LowerBoundCoveringRadiusSphereCovering( <n>, <r> [, <F> ] [, true ] )
##

LowerBoundCoveringRadiusSphereCovering := function ( arg )
    
    local n, m, q, lb, ub, tmp, tmpcr, tmpsize, t;
    
    if arg[ Length( arg ) ] = false then
        # last argument is false, try to find a lower bound for
        # the covering radius, given length and size of the code, and
        # optional the field. default is GF(2)
        if Length( arg ) = 4 then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCRSphereCovering: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;
        
        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCRSphereCovering: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error( "LBCRSphereCovering: <n> must be > 0" );
        fi;
        
        m := arg[ 2 ];
        if not IsInt( m ) then
            Error( "LBCRSphereCovering: <m> must be an integer" );
        fi;
        if m <= 0 or m > q^n then
            Error( "LBCRSphereCovering: <m> must be > 0 and <= q^n" );
        fi;
        
        # everything is set up. now compute the bound
        
        lb := 0;
        ub := n;
        while lb <> ub do
            tmpcr := IntFloor( ( lb + ub ) / 2 );
            tmpsize := IntCeiling( q^n / SphereContent( n, tmpcr, q ) );
            if tmpsize > m then
                lb := tmpcr + 1;
            else
                ub := tmpcr;
            fi;
        od;
        return ub;
    else
        # the last argument is not false
        # now it is assumed that the first argument is the length
        # of the code and the second argument is the covering radius
        # of the code. a lower bound for the minimal size of the
        # code is returned
        if Length( arg ) = 4 or
           ( Length( arg ) = 3 and arg[ 3 ] <> true ) then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCRSphereCovering: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;
        
        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCRSphereCovering: <n> must be an integer" );
        fi;
        if n<=0 then
            Error( "LBCRSphereCovering: <n> must be positive" );
        fi;
        
        t := arg[ 2 ];
        if not IsInt( t ) then
            Error( "LBCRSphereCovering: <t> must be an integer" );
        fi;
        if t < 0 or t > n then
            Error( "LBCRSphereCovering: <t> must be >= 0 and <= <n>" );
        fi;
            
        return IntCeiling( q^n / SphereContent( n, t, q ) );
    fi;
end;
        

########################################################################
##
#F  LowerBoundCoveringRadiusVanWee1( ... )
##

LowerBoundCoveringRadiusVanWee1 := function ( arg )
    
    local n, q, m, lb, ub, tmpcr, tmpsize, t, tmp;
    
    if arg[ Length( arg ) ] = false then
        # last argument is false, try to find a lower bound for
        # the covering radius, given length and size of the code, and
        # optional the field. default is GF(2)

        if Length( arg ) = 4 then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCRVanWee1: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;

        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCRVanWee1: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error( "LBCRVanWee1: <n> must be > 0" );
        fi;
        
        m := arg[ 2 ];
        if not IsInt( m ) then
            Error( "LBCRVanWee1: <m> must be an integer" );
        fi;
        if m <= 0 or m > q^n then
            Error( "LBCRVanWee1: <m> must be > 0 and <= q^n" );
        fi;
        
        # everything is set up. now compute the bound
        
        lb := 0;
        ub := n;
        while lb <> ub do
            tmpcr := IntFloor( (ub+lb) / 2 );
            tmpsize := LowerBoundCoveringRadiusVanWee1( n, tmpcr, q );
            if tmpsize > m then
                lb := tmpcr + 1;
            else
                ub := tmpcr;
            fi;
        od;
        return ub;
    else
        # the last argument is not false
        # now it is assumed that the first argument is the length
        # of the code and the second argument is the covering radius
        # of the code. a lower bound for the minimal size of the
        # code is returned
        if Length( arg ) = 4 or ( Length( arg ) = 3 and arg[ 3 ] <> true ) then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCRVanWee1: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;
        
        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCRVanWee1: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error(" LBCRVanWee1: <n> must be > 0" );
        fi;
        
        t := arg[ 2 ];
        if not IsInt( t ) then
            Error( "LBCRVanWee1: <t> must be an integer" );
        fi;
        if t < 0 or t > n then
            Error( "LBCRVanWee1: <t> must be >= 0 and <= <n>" );
        fi;
        if t = n then
            return 1;
        fi;
        
        tmp := (Binomial(n,t))/(IntCeiling((n-t)/(t+1)));
        tmp := tmp * (IntCeiling((t+1)/(t+1)) - (t+1)/(t+1));
        tmp := SphereContent( n, t, q ) - tmp;
        return IntCeiling( q^n / tmp );
    fi;
end;


#############################################################################
##
#F  LowerBoundCoveringRadiusVanWee2( <n>, <r> ) Counting Excess bound
##

LowerBoundCoveringRadiusVanWee2 := function ( arg )
    
    local n, m, lb, ub, tmpcr, eps, tmpsize, t, tmpb1, tmpb2;
    
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "Usage: LBCRVanWee2( <n>, <t|size> [, true|false ] )" );
    fi;
    
    n := arg[ 1 ];
    if not IsInt( n ) then
        Error( "LBCRVanWee2: <n> must be an integer" );
    fi;
    if n<=0 then
        Error( "LBCRVanWee2: <n> must be positive" );
    fi;
    
    if Length( arg ) = 3 and arg[ 3 ] = false then
        m := arg[ 2 ];
        if not IsInt( m ) then
            Error( "LBCRVanWee2: <m> must be an integer" );
        fi;
        if m <= 0 or m > 2^n then
            Error( "LBCRVanWee2: <m> must be > 0 and <= 2^n" );
        fi;
        lb := 0;
        ub := IntFloor( n/2 );
        while lb <> ub do
            tmpcr := IntFloor( ( ub + lb ) / 2 );
            tmpb1 := (tmpcr+2)*(tmpcr+1)/2; # Binomial(tmpcr+2,2)
            tmpb2 := (n-tmpcr+1)*(n-tmpcr)/2; # Binomial(n-tmpcr+1,2)
            eps := tmpb1 * IntCeiling( tmpb2/tmpb1 ) - tmpb2;
            tmpsize := 2^n * ( SphereContent( n, 2, 2 ) + eps
                               - 1/2*( tmpcr + 2 )*( tmpcr - 1 ) );
            tmpsize := tmpsize / ( SphereContent( n, tmpcr, 2 ) * 
                               ( SphereContent( n, 2, 2 ) 
                                 - 1/2*( tmpcr + 2 )*( tmpcr - 1 ) )
                               + eps * SphereContent( n, tmpcr - 2 ) );
            if tmpsize > m then
                lb := tmpcr + 1;
            else
                ub := tmpcr;
            fi;
        od;
        return ub;
    else
        if Length( arg ) = 3 and arg[ 3 ] <> true then
            Error("Usage: LBCRVanWee2( <n>, <t> [, <true|false> ] )");
        fi;
        t := arg[ 2 ];
        if not IsInt( t ) then
            Error( "LBCRVanWee2: <t> must be an integer" );
        fi;
        if t < 0 or t > n then
            Error( "LBCRVanWee2: <t> must be >= 0 and <= <n>" );
        fi;
        if 2 * t > n then
            return 0;
        fi;
        tmpb1 := (t+2)*(t+1)/2;
        tmpb2 := (n-t+1)*(n-t)/2;
        eps := tmpb1 * IntCeiling( tmpb2/tmpb1 ) - tmpb2;
        tmpsize := 2^n * ( SphereContent( n, 2, 2 ) + eps
                           - 1/2*( t + 2 )*( t - 1 ) );
        tmpsize := tmpsize / ( SphereContent( n, t, 2 ) 
                           * ( SphereContent( n, 2, 2 ) 
                               - 1/2*( t + 2 )*( t - 1 ) ) 
                           + eps * SphereContent( n, t-2, 2 ) );
        return IntCeiling(tmpsize);
    fi;
end;


#############################################################################
##
#F  LowerBoundCoveringRadiusCountingExcess( <n>, <r> )
##
LowerBoundCoveringRadiusCountingExcess := function ( arg )
    
    local n, m, lb, ub, tmpcr, tmpsize, t, rho, eps;
    
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "Usage: LBCRCountingExcess( <n>, <t|m> [, true|false ] )" );
    fi;
    
    n := arg[ 1 ];
    if not IsInt( n ) then
        Error( "LBCRCountingExcess: <n> must be an integer" );
    fi;
    if n<=0 then
        Error( "LBCRCountingExcess: <n> must be positive" );
    fi;
    
    if Length( arg ) = 3 and arg[ 3 ] = false then
        m := arg[ 2 ];
        if not IsInt( m ) then
            Error( "LBCRCountingExcess: <m> must be an integer" );
        fi;
        if m<=0 or m > 2^n then
            Error( "LBCRCountingExcess: <m> must be > 0 and <= 2^n" );
        fi;
        
        lb := 0;
        ub := IntFloor( ( n-1 ) / 2 );
        while lb <> ub do
            tmpcr := IntFloor( ( ub + lb ) / 2 );
            tmpsize := LowerBoundCoveringRadiusCountingExcess(
                               n, tmpcr, true );
            if tmpsize > m then
                lb := tmpcr + 1;
            else
                ub := tmpcr;
            fi;
        od;
        return ub;
    else
        if Length( arg ) = 3 and arg[ 3 ] <> true then
            Error( "Usage: LBCRCountingExcess( <n>, <t|m> [, true|false ] )");
        fi;
        t := arg[ 2 ];
        if not IsInt( t ) then
            Error( "LBCRCountingExcess: <t> must be an integer" );
        fi;
        if t < 0 or t > n then
            Error( "LBCRCountingExcess: <t> must be >=0 and <= <n>" );
        fi;
        
        if t < 2 or 2 * t + 1 > n then
            return 0;
        fi;
        if t = 2 then
            rho := n - 3 + 2/n;
        else
            rho := n - t - 1;
        fi;
        eps := ( t+1 ) * IntCeiling( ( n+1 ) / ( t+1 ) ) - ( n+1 );
        if eps > t then 
            return 0;
        fi;
        tmpsize := 2^n * ( rho + eps );
        tmpsize := tmpsize / ( rho * SphereContent( n, t ) 
                           + eps * SphereContent( n, t-1 ) );
        
        return IntCeiling( tmpsize );
    fi;
end;


########################################################################
##
#F  LowerBoundCoveringRadiusEmbedded1( <n>, <r> [, <givesize> ] )
##

LowerBoundCoveringRadiusEmbedded1 := function ( arg )
    
    local q, n, m, lb, ub, tmpcr, tmpsize, t, upperb;
    
    if arg[ Length( arg ) ] = false then
        # last argument is false, try to find a lower bound for
        # the covering radius, given length and size of the code, and
        # optional the field. default is GF(2)

        if Length( arg ) = 4 then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCREmbedded1: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;
        
        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCREmbedded1: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error( "LBCREmbedded1: <n> must be > 0" );
        fi;
        
        m := arg[ 2 ];
        if not IsInt( m ) then
            Error( "LBCREmbedded1: <m> must be an integer" );
        fi;
        if m <= 0 or m > q^n then
            Error( "LBCREmbedded1: <m> must be > 0 and <= q^n" );
        fi;
        
        # everything is set up. now compute the bound
        
        if m = q^n then
            return 0;
        fi;
        if n = 1 then
            return 0;
        fi;
        
        lb := 1;   
        ub := n;
        while lb <> ub do
            tmpcr := IntFloor( (ub+lb) / 2 );
            tmpsize := SphereContent( n, tmpcr, q )
                       - Binomial( 2*tmpcr, tmpcr ); 
            if tmpsize = 0 then
                lb := lb + 1;
            elif tmpsize < 0 then
                ub := tmpcr;
            else
                if 2*tmpcr+1 > n then
                    upperb := 1;
                else
                    upperb := UpperBound( n, 2*tmpcr+1, q );
                fi;
                tmpsize := IntCeiling( ( q^n - upperb
                                   * Binomial( 2 * tmpcr, tmpcr ) )
                                   / tmpsize );
                if tmpsize > m then
                    lb := tmpcr + 1;
                else
                    ub := tmpcr;
                fi;
            fi;
        od;
        return ub;
    else
        # the last argument is not false
        # now it is assumed that the first argument is the length
        # of the code and the second argument is the covering radius
        # of the code. a lower bound for the minimal size of the
        # code is returned
        if Length( arg ) = 4 or
           ( Length( arg ) = 3 and arg[ 3 ] <> true ) then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCREmbedded1: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;

        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCREmbedded1: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error(" LBCREmbedded1: <n> must be > 0" );
        fi;
        
        t := arg[ 2 ];
        if not IsInt( t ) then
            Error( "LBCREmbedded1: <t> must be an integer" );
        fi;
        if t < 0 or t > n then
            Error( "LBCREmbedded1: <t> must be >= 0 and <= <n>" );
        fi;
        
        tmpsize := SphereContent( n, t, q ) - Binomial( 2*t, t );
        if tmpsize <= 0 then
            return 0;
        else
            if 2 * t + 1 > n then
                upperb := 1;
            else
                upperb := UpperBound( n, 2*t+1, q );
            fi;
            return IntCeiling( ( q^n - upperb
              * Binomial( 2*t, t ) ) / tmpsize );
        fi;
    fi;
end;


########################################################################
##
#F  LowerBoundCoveringRadiusEmbedded2( <n>, <r> [, <givesize> ] )
##

LowerBoundCoveringRadiusEmbedded2 := function ( arg )
    
    local q, n, m, lb, ub, tmpcr, tmpsize, t, upperb;
    
    if arg[ Length( arg ) ] = false then
        # last argument is false, try to find a lower bound for
        # the covering radius, given length and size of the code, and
        # optional the field. default is GF(2)
        if Length( arg ) = 4 then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCREmbedded2: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;
        
        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCREmbedded2: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error( "LBCREmbedded2: <n> must be > 0" );
        fi;
        
        m := arg[ 2 ];
        if not IsInt( m ) then
            Error( "LBCREmbedded2: <m> must be an integer" );
        fi;
        if m <= 0 or m > q^n then
            Error( "LBCREmbedded2: <m> must be > 0 and <= q^n" );
        fi;
        
        # everything is set up. now compute the bound
        
        if m = q^n then
            return 0;
        fi;
        if n = 1 or n = 2 then
            return 0;
        fi;
        
        lb := 1;   
        ub := n;
        while lb <> ub do
            
            tmpcr := IntFloor( (ub+lb) / 2 );
            tmpsize := SphereContent( n, tmpcr, q )
                       - 3/2 * Binomial( 2*tmpcr, tmpcr ); 
            if tmpsize = -1/2 then
                lb := lb + 1;
            elif tmpsize <= 0 then
                ub := tmpcr;
            else
                if 2*tmpcr+1 > n then
                    upperb := 1;
                else
                    upperb := UpperBound( n, 2*tmpcr+1, q );
                fi;
                tmpsize := IntCeiling( ( q^n - 2 * upperb
                                   * Binomial( 2*tmpcr, tmpcr ) )
                                   / tmpsize );
                if tmpsize > m then
                    lb := tmpcr + 1;
                else
                    ub := tmpcr;
                fi;
            fi;
        od;
        return ub;
    else
        # the last argument is not false
        # now it is assumed that the first argument is the length
        # of the code and the second argument is the covering radius
        # of the code. a lower bound for the minimal size of the
        # code is returned
        if Length( arg ) = 4 or
           ( Length( arg ) = 3 and arg[ 3 ] <> true ) then
            if IsInt( arg[ 3 ] ) then
                q := arg[ 3 ];
            else
                if not IsField( arg[ 3 ] ) or not IsFinite( arg[ 3 ] ) then
                    Error( "LBCREmbedded2: ", arg[3],
                           " is not a finite field" );
                else
                    q := Size( arg[ 3 ] );
                fi;
            fi;
        else
            q := 2;
        fi;
        
        n := arg[ 1 ];
        if not IsInt( n ) then
            Error( "LBCREmbedded2: <n> must be an integer" );
        fi;
        if n <= 0 then
            Error(" LBCREmbedded2: <n> must be > 0" );
        fi;
        
        t := arg[ 2 ];
        if not IsInt( t ) then
            Error( "LBCREmbedded2: <t> must be an integer" );
        fi;
        if t < 0 or t > n then
            Error( "LBCREmbedded2: <t> must be >= 0 and <= <n>" );
        fi;
        
        tmpsize := SphereContent( n, t, q ) - 3/2*Binomial( 2*t, t );
        if tmpsize <= 0 then
            return 0;
        else
            if 2*t+1 > n then
                upperb := 1;
            else
                upperb := UpperBound( n, 2*t+1, q );
            fi;
            
            return IntCeiling( ( q^n - 2*upperb
                           * Binomial( 2*t, t ) ) / tmpsize );
        fi;
        
    fi;
    
end;


#############################################################################
##
#F  LowerBoundCoveringRadiusInduction( <n>, <r> ) Induction bound
##

LowerBoundCoveringRadiusInduction := function( arg )
    
    local n, t;
    
    if Length( arg ) <> 2 then
        Error( "Usage: LBCRInduction( <n>, <t> )" );
    fi;
    
    n := arg[ 1 ];
    if not IsInt( n ) then
        Error( "LBCRInduction: <n> must be an integer" );
    fi;
    if n <= 0 then 
        Error( "LBCRInduction: <n> must be positive" );
    fi;

    t := arg[ 2 ];
    if not IsInt( t ) then
        Error( "LBCRInduction: <t> must be an integer" );
    fi;
    if t < 0 or t > n then
        Error( "LBCRInduction: <t> must be >= 0 and <= <n>" );
    fi;
    
    if n = 2 * t + 2 and t >= 1 then
        return 4;
    elif n = 2 * t + 3 and t >= 1 then
        return 7;
    elif n = 2 * t + 4 and t >= 4 then
        return 8;
    else
        return 0;
    fi;
    
end;

        
########################################################################
##
#F  GeneralUpperBoundCoveringRadius( <code> )
##

GeneralUpperBoundCoveringRadius := function ( code )
    
    local listofbounds, ubdelsarte;

    if not IsCode( code ) then
        Error( "GeneralUpperBoundCoveringRadius: <code> must ",
               "be a code" );
    fi;

    listofbounds := [ UpperBoundCoveringRadiusStrength( code ) ];
    if WordLength( code ) <= 100 then
        Append( listofbounds, [
          UpperBoundCoveringRadiusDelsarte( code )
                ] );
    fi;
    if IsLinearCode( code ) then
        Append( listofbounds, [
          UpperBoundCoveringRadiusRedundancy( code ),
        ] );
        if Field( code ) = GF(2) then
            if ( IsBound( code.lowerBoundMinimumDistance ) and
                 IsBound( code.upperBoundMinimumDistance ) ) and
               code.lowerBoundMinimumDistance = 
                 code.upperBoundMinimumDistance 
               then
                Append( listofbounds, [
                        UpperBoundCoveringRadiusGriesmerLike( code )
                        ] );
            fi;
        fi;
    fi;
    if IsCyclicCode( code ) and Field( code ) = GF(2) then
        Append( listofbounds, [
          UpperBoundCoveringRadiusCyclicCode( code )
        ] );
    fi;
    if IsBound( code.boundsCoveringRadius ) then
        Add( listofbounds, Maximum( code.boundsCoveringRadius ) );
    fi;
    return Minimum( listofbounds );
end;


########################################################################
##
#F  UpperBoundCoveringRadiusRedundancy( <code> )
##
##  Return the redundancy of the code as an upper bound for
##  the covering radius.
##
##  Only for linear codes.
##

UpperBoundCoveringRadiusRedundancy := function ( code )

    if not IsCode( code ) or not IsLinearCode( code ) then
        Error( "UBCRRedundancy: <code> must be a linear code" );
    fi;

    return Redundancy( code );

end;


########################################################################
##
#F  UpperBoundCoveringRadiusDelsarte( <code> )
##

UpperBoundCoveringRadiusDelsarte := function ( code )
    
    local dual, wddual, p;
    
    if IsLinearCode( code ) then
        if not IsBound( code.boundsCoveringRadius ) then
            # avoid recursion
            code.boundsCoveringRadius := [ 0 .. WordLength( code ) ];
        fi;
        dual := DualCode( code );
        wddual := WeightDistribution( dual );
        return WeightVector( wddual ) - 1;
    else
        p := CodeMacWilliamsTransform( code );
        p := p.coefficients;
        return WeightVector( p ) - 1;
    fi;
end;


########################################################################
##
#F  UpperBoundCoveringRadiusStrength( <code> )
##
##  Return (q-1)n/q as an upper bound for <code>, if it
##  has strength 1 (i.e. every coordinate contains each element
##  of the field the same number of times).
##

UpperBoundCoveringRadiusStrength := function ( code )

    local q, n, UnrestrictedStrength1, LinearStrength1;

    if not IsCode( code ) then
        Error( "UBCRStrength: <code> must be a code" );
    fi;

    UnrestrictedStrength1 := function( code )
        local stris1, i, j, n, q, els, fieldels, coordlist, number, zerocols;
        stris1 := true;
        i := 1;
        n := WordLength( code );
        els := VectorCodeword( Elements( code ) );
        q := Size( Field( code ) );
        if not ( Length( els ) in List( [ 1 .. n ], x -> q^x ) ) then
            stris1 := false;
        fi;
        fieldels := Elements( Field( code ) );
        zerocols := 0;
        while stris1 and i <= n do
            coordlist := List( els, x -> x[ i ] );
            for j in fieldels do
                number := Length( Filtered( coordlist, x -> x = j ) );
                if number = Length( els ) then
                    zerocols := zerocols + 1;
                elif number <> Length( els ) / q then
                    stris1 := false;
                fi;
            od;
            i := i + 1;
        od;
        if stris1 then
            return IntFloor((q-1)*(n-zerocols)/q)+zerocols;
        else
            return n;
        fi;
    end;

    LinearStrength1 := function( code )
        local zerocols, i, j, genmat, n, k, zero, onlyzeroes;

        genmat := GeneratorMat( code );
        n := WordLength( code );
        k := Dimension( code );
        zerocols := 0;
        zero := Field(code).zero;
        for i in [ 1 .. n ] do
            onlyzeroes := true;
            j := 1;
            while onlyzeroes and j <= k do
                onlyzeroes := ( genmat[ j ][ i ] = zero );
                j := j + 1;
            od;
            if onlyzeroes then
                zerocols := zerocols + 1;
            fi;
        od;

        return IntFloor((q-1)*(n-zerocols)/q) + zerocols;
    end;

    q := Size( Field( code ) );
    n := WordLength( code );
    if not IsLinearCode( code ) then
        return UnrestrictedStrength1( code );
    else
        if Dimension( code ) > 0 then
            return LinearStrength1( code );
        else
            return n;
        fi;
    fi;
end;



########################################################################
##
#F  UpperBoundCoveringRadiusGriesmerLike( <code> )
##

UpperBoundCoveringRadiusGriesmerLike := function ( code )

    local q;

    if not IsCode( code ) or not IsLinearCode( code ) then
        Error( "UBCRGriesmerLike: <code> is not a linear code" );
    fi;

    q := Size( Field( code ) );
    return WordLength( code ) - Sum( [ 1 .. Dimension( code ) ],
                   x -> IntCeiling( MinimumDistance( code ) / q^x ) );
end;


########################################################################
##
#F  UpperBoundCoveringRadiusCyclicCode( <code> )
##

UpperBoundCoveringRadiusCyclicCode := function ( code )
    
    if not IsCode( code ) or not IsCyclicCode( code ) then
        Error( "UBCRCyclicCode: <code> is not a cyclic code" );
    fi;
    
    return WordLength( code ) - Dimension( code ) + 1 -
           IntCeiling( WeightCodeword(
                Codeword( GeneratorPol( code ) ) ) / 2 );
end;


