##  mytz.g,  version 27/ 8/96

##############################################################################
##
#F  TzCommutatorPair( <Tietze record>, <relator> ) . . . . . . . . . . . . . . 
##                                test whether a relator is a commutator [x,y]
##
TzCommutatorPair := function( T, rel )

    local  tietze, numgens, invs, numinvs, 
           x, ax, ix, px, y, ay, iy, py, pair;

    # check the given argument to be a Tietze record.
    if not ( IsBound( T.isTietze ) and T.isTietze ) then
        Error( "argument must be a Tietze record" );
    fi;
    tietze := T.tietze;
    invs := tietze[TZ_INVERSES];
    numgens := tietze[TZ_NUMGENS];
    numinvs := 1 + 2 * numgens;
    if not ( IsList( rel ) and ( Length( rel ) = 4 ) ) then
        Print( "Second parameter must be a relator of length 4.\n" );
        return [ 0, 0 ];
    fi;
    x := rel[1];
    y := rel[2];
    ax := AbsInt( x );            
    ay := AbsInt( y );
    px := Position( invs, x );
    ix := invs[numinvs + 1 - px];
    py := Position( invs, y );
    iy := invs[numinvs + 1 - py];
    if ( ( rel[3] = ix ) and ( rel[4] = iy ) ) then
        pair := Set( [ ax, ay ] );
    else
        pair := [ 0, 0 ];
    fi;
    return pair;    
end;

##############################################################################
##
#F  TzPartition( <Tietze record> ) partition generators into commuting subsets
##
TzPartition := function( T )

    local  tietze, numgens, gens, partition, count, bipartite,
           invs, allgens, others, rels, numrels, lengths, numpart, numold,
           numinvs, i, j, k, x, ax, y, ay, z, r, lr, s, ok,
           commpairs, c, d, pair, commutators, comm1, comm2,
           left, right, rest, powerof, powers, freq, fnum, L, S;

    # check the given argument to be a Tietze record.
    if not ( IsBound( T.isTietze ) and T.isTietze ) then
        Error( "argument must be a Tietze record" );
    fi;
    # if IsBound( T.partition ) then
    #     return T.partition;
    # fi;

    tietze := T.tietze;
    invs := tietze[TZ_INVERSES];
    rels := tietze[TZ_RELATORS];
    lengths := tietze[TZ_LENGTHS];
    numrels := tietze[TZ_NUMRELS];
    numgens := tietze[TZ_NUMGENS];
    numinvs := 1 + 2 * numgens;
    allgens := [1..numgens];
    left := allgens;
    partition := [ Set( allgens ) ];

    # phase1: find explicit commutators
    commpairs := [ ];
    for j in [1..numgens] do
        Add( commpairs, [ ] );
    od;
    commutators := 0 * [1..numrels];
    comm1 := 0;
    comm2 := 0;
    for j in [1..numrels] do
        if ( lengths[j] <> 4 ) then
            commutators[j] := false;
        else
            r := rels[j];
            pair := TzCommutatorPair( T, r );
            if ( pair <> [ 0, 0 ] ) then
                commutators[j] := true;
                comm2 := comm2 + 1;
                ax := pair[1];
                ay := pair[2];
                if ( ( ax in allgens ) and ( ay in allgens ) ) then
                    Add( commpairs[ax], ay );
                    Add( commpairs[ay], ax );
                fi;
            else
                commutators[j] := false;
            fi;
        fi;
    od;
    commpairs := List( commpairs, L -> Set( L ) );
    if ( ( T.printLevel > 1 ) and ( comm2 > 0 ) ) then
        Print( "There were ", comm2, " commutators found in phase 1\n" );
    fi;

    # phase2: find implicit commutators
    while ( comm1 < comm2 ) do
        comm1 := comm2;
        for j in [1..numrels] do 
            r := rels[j];
            lr := Length( r );
            if ( lr > 1 ) then
                freq := Collected( List( r, x -> AbsInt( x ) ) );
                fnum := Length( freq );
                for x in [1..fnum] do
                    L := freq[x];
                    if ( L[2] = 1 ) then
                        S := allgens;
                        for y in [1..fnum] do
                            if ( x <> y ) then
                                z := freq[y][1];
                                S := Intersection( S, commpairs[z] );
                            fi;
                        od;
                        y := L[1];
                        commpairs[y] := Union( commpairs[y], S );
                        for z in S do
                            commpairs[z] := Union( commpairs[z], [y] );
                        od;
                    fi;
                od;
            fi;
        od;
    comm2 := Sum( List( commpairs, L -> Length( L ) ) )/2;
    if ( ( T.printLevel > 1 ) and ( comm2 > comm1 ) ) then
        Print( "There were ", comm2 - comm1 );
        Print( " commutators found in phase 2\n" );
    fi;
    od;

    # phase3: separate into parts
    numpart := 1;
    numold := 0;
    while ( ( numold <> numpart ) and ( Length( partition[1] ) > 1 ) ) do
        if ( numpart > 1 ) then
            commpairs := List( commpairs, 
                                  L -> Difference( L, partition[2] ) );
        fi;
        gens := partition[1];
        others := Difference( allgens, gens );
        rest := Sublist( partition, [ 2..Length( partition ) ] );
        bipartite := false;
        left := commpairs[ gens[1] ];
        if ( left <> [ ] ) then
            right := Difference( gens, left );
            bipartite := true;
            for c in commpairs do
                d := Difference( c, left );
                if not ( c = [ ] ) then
                    if not ( ( c = left ) or ( d = right ) ) then
                        bipartite := false;
                    fi;
                fi;
            od;
        fi;
        if bipartite then
            # reorder the letters in each relator (but not the commutators)
            partition := Concatenation( [ left ], [ right ], rest );
            for j in [1..numrels] do
                if not commutators[j] then
                    r := rels[j];
                    s := [ ];
                    for x in r do
                        ax := AbsInt( x );
                        if ( ax in left ) then
                            Add( s, x );
                        fi;
                    od;
                    for y in r do
                        ay := AbsInt( y );
                        if ( ay in right ) then
                            Add( s, y );
                        fi;
                    od;
                    for x in r do
                        ax := AbsInt( x );
                        if ( ax in others ) then
                            Add( s, x );
                        fi;
                    od;
                    if ( r <> s ) then
                        rels[j] := s;
                        tietze[TZ_MODIFIED] := true;
                    fi;
                fi;
            od;
            if ( T.printLevel > 1 ) then
                Print( "#I factoring ", gens, " into " );
                Print( left, ", ", right, "\n" );
            fi;
        else
            left := [ ];
        fi;
    numold := numpart;
    numpart := Length( partition );
    od;
    T.partition := partition;
    return partition;
end;

##############################################################################
##
#F  FactorsPresentation( <T> [,<print level>] ) . . . . . . . . . . . . . . . 
##                            factor T as commuting presentations [T1,T2,...]
##                               with one factor for each part of T.partition
##
FactorsPresentation := function( arg )

    local  T, printlevel, tietze, gens, rels, numrels, numgens, invs, numinvs,
           lengths, flags, partition, part, numpart, i, j, k, rel, diff, F,
           commutator, fx, fy, rel, pair, factor, fac, len, chosen, posn,
           subrels, subnumi, subtot, sublen, subflags, subtriv, subT;

    # check that the first argument is a tietze record
    T := arg[1];
    if not ( IsBound( T.isTietze ) and T.isTietze ) then
        Error( "argument must be a presentation record" );
    fi;
    # check that the second argument is an integer
    printlevel := 1;
    if ( Length( arg ) = 2 ) then
        printlevel := arg[2];
    fi;
    if not IsInt( printlevel ) then
        Error( "second argument must be an integer" );
    fi;

    partition := TzPartition( T );
    if ( Length( partition ) = 1 ) then
        return [ T ];
    fi;
    tietze := T.tietze;
    invs := tietze[TZ_INVERSES];
    numgens := tietze[TZ_NUMGENS];
    numinvs := 1 + 2 * numgens;
    gens := T.generators;
    rels := tietze[TZ_RELATORS];
    lengths := tietze[TZ_LENGTHS];
    flags := tietze[TZ_FLAGS];
    numrels := tietze[TZ_NUMRELS];
    numpart := Length( partition );
    factor := 0 * [1..numgens ];
    for i in [1..numpart] do
        part := partition[i];
        for j in [1..Length(part)] do
            factor[ part[j] ] := i;
        od;
    od;
    chosen := List( [1..numpart], i -> [ ] );
    for j in [1..numrels] do
        rel := rels[j];
        len := lengths[j];
        commutator := false;
        if ( len = 4 ) then
            pair := TzCommutatorPair( T, rel );
            if ( pair <> [ 0, 0 ] ) then
                fx := factor[ pair[1] ];
                fy := factor[ pair[2] ];
                commutator := ( fx <> fy );
            fi;
        fi;
        if not commutator then
            fac := List( rel, i -> factor[ AbsInt( i ) ] );
            # check that fac is increasing
            for i in [2..len] do
                if ( fac[i] < fac[i-1] ) then
                    Print( "relator = ", rel, "\n" );
                    Error( "factors mixed in relator" );
                fi;
            od;
            for i in Set( fac ) do
                Add( chosen[i], j );
            od;           
        fi;
    od;
    # construct Tietze records for each factor
    subflags := 0 * [1..numpart];
    F := List( [1..numpart], i -> Copy( T ) );
    for i in [1..numpart] do
        subrels := [ ];
        for j in chosen[i] do
            rel := rels[j];
            len := lengths[j];
            fac := List( rel, k -> factor[ AbsInt( k ) ] );
            posn := Filtered( [1..len], k -> fac[k]=i );
            Add( subrels, Sublist( rel, posn ) );
        od;
        sublen := List( subrels, L -> Length( L ) );
        subflags := List( chosen[i], k -> flags[k] );
        F[i].printLevel := printlevel;
        subT := F[i].tietze;
        diff := Difference( [1..numgens], partition[i] );
        subtriv := List( diff, x -> [x] );
        subrels := Concatenation( subtriv, subrels );
        sublen := Concatenation( List( subtriv, x -> 1 ), sublen );
        subflags := Concatenation( List( subtriv, x -> 0 ), subflags );
        subnumi := Length( subrels );
        subtot := 0;
        for j in [1..subnumi] do
            subtot := subtot + sublen[j];
        od;
        subT[TZ_TOTAL] := subtot;
        subT[TZ_RELATORS] := subrels;
        subT[TZ_NUMRELS] := subnumi;
        subT[TZ_LENGTHS] := sublen;
        subT[TZ_FLAGS] := subflags;
        subT[TZ_STATUS] := [ subT[1], subT[2], subT[3] ];
    od;
    return F;
end;

##############################################################################
##
#F  TzRenumber( <T>, <L> ) . . . . . . renumber generators of T according to L
##

TzRenumber := function( T, L )

    local  R, S, P, Q, M, ngen, i, j, k, invT, invR, relT, relR, numrel,
           r, rel, lens, len, partT, partR, part, imT, imR, im;

    if not ( IsBound( T.isTietze ) and T.isTietze ) then
        Error( "first argument must be a tietze record" );
    fi;
    if not IsList( L ) then
        Error( "second argument must be a list" );
    fi;
    ngen := Length( T.generators );
    S := Set( L );
    if not ( ( Length( L ) = ngen ) and ( S = [1..ngen] ) ) then
        Error( "<L> must determine a permutation of the components" );
    fi;
    P := PermList( L );
    Q := P^-1;
    M := ListPerm( Q );
    R := Copy( T );
    invT := T.tietze[ TZ_INVERSES ];
    invR := R.tietze[ TZ_INVERSES ];
    for i in [1..ngen] do
        j := ngen + 1 + i;
        R.generators[i] := T.generators[ L[i] ];
        k := invT[j];
        if ( k = -i ) then
            invR[ ngen + 1 + M[i] ] := - M[i];
        elif ( k = i ) then
            invR[ ngen + 1 + M[i] ] := M[i];
        else
            Error( "cannot determine inverses" );
        fi;
    od;
    numrel := T.tietze[ TZ_NUMRELS ];
    lens := T.tietze[ TZ_LENGTHS ];
    relT := T.tietze[ TZ_RELATORS ];
    relR := R.tietze[ TZ_RELATORS ];
    for r in [1..numrel] do
        len := lens[r];
        rel := relT[r];
        relR[r] := 0 * [1..len];
        for i in [1..len] do
            j := rel[i];
            if ( j > 0 ) then
                relR[r][i] := M[j];
            else
                relR[r][i] := -M[-j];
            fi;
        od;
    od;
    if IsBound( T.partition ) then
        partT := T.partition;
        partR := R.partition;
        for i in [ 1..Length( partT ) ] do
            part := partT[i];
            len := Length( part );
            partR[i] := 0 * [1..len];
            for j in [1..len] do
                partR[i][j] := L[ part[j] ];
            od;
        od;
    fi;
    if IsBound( T.imagesOldGens ) then
        imT := T.imagesOldGens;
        imR := R.imagesOldGens;
        for i in [ 1..Length( imT ) ] do
            im := imT[i];
            len := Length( im );
            imR[i] := 0 * [1..len];
            for j in [1..len] do
                imR[i][j] := M[ im[j] ];
            od;
        od;
        imT := T.preImagesNewGens;
        imR := R.preImagesNewGens;
        for i in [ 1..Length( imT ) ] do
            imR[i] := imT[ L[i] ];
        od;
    fi;
    return R;
end;
