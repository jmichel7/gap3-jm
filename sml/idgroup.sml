if not IsBound( IdGroupTree ) then
    IdGroupTree := rec( fp := [ 1 .. 1000 ], next := [ ] );
fi;
if not IsBound(StandardPresentation) then StandardPresentation := 0; fi;
if not IsBound( InfoIdGroup1 ) then InfoIdGroup1 := Ignore; fi;

#############################################################################
##
#F  AgGroupCode( code ) . . . . . . . .construct ag-group from numerical code
##
AgGroupCode := function( code )
    local n1, size, indices, f, l, mi, n, t1, indices, gens, rels, g, i, 
          uc, ll, rr, t, j, z, z2, result;

    # set up
    g    := ["a", "b", "c", "d", "e", "f", "g", "h"];
    gens := [];
    rels := [];
    rr   := [];

    # get the size
    size:= code mod 1000;
    if size = 0 then 
        size := 1000;
        code := code - 1000;
    fi;
    n := QuoInt( code, 1000 );

    # single out 1-group
    if size = 1 then
        return AgGroupFpGroup( rec( generators := [ ] ,relations := [ ] ) );
    fi;

    # get relative orders
    f := Factors(size);
    l := Length(f);
    mi:= Maximum(f)-1;
    if Length(Set(f)) > 1 then
        if mi^l < 2^28 then
            indices:=CoefficientsInt(List([1..l], x->mi),n mod (mi^l))+2;
        else
            indices := [ ];
            n1 := n mod (mi^l);
            for i in Reversed( [1..l] ) do
                indices[ i ] := ( n1 mod mi ) + 2;
                n1 := QuoInt( n1, mi );
            od;
        fi;
        n:=QuoInt(n,mi^l);
    else
        indices := f;
    fi;

    # construct generators and set up for power relations
    for i in [1..l] do
        gens[i] := AbstractGenerator(g[i]);
        rels[i] := gens[i]^indices[i];
    od;

    # set up non-trivial relators
    ll:=l*(l+1)/2-1;
    if ll < 28 then
        uc:=Reversed(CoefficientsInt(List([1..ll],x->2),n mod (2^ll)));
    else
        uc := [];
        n1 := n mod (2^ll);
        for i in [1..ll] do
            uc[i] := n1 mod 2;
            n1 := QuoInt( n1, 2 );
        od;
    fi;
    n:=QuoInt(n,2^ll);

    # construct non-trivial relators - get tails
    for i in [1..Sum(uc)] do
        t:=CoefficientsInt(indices,n mod size);
        g:=gens[1]^0;
        for j in [1..l] do
            if t[j] > 0 then 
                g:=g*gens[j]^t[j];
            fi;
        od;
        Add(rr,g);
        n:=QuoInt(n,size);
    od;

    # compute non-trivial power relators
    z:=1;
    for i in [1..l-1] do
        if uc[i] = 1 then
            rels[i] := rels[i]/rr[z];
            z:=z+1;
        fi;
    od;
    z2:=l-1;

    # compute non-trivial commutator relators
    for i in [1..l] do
        for j in [i+1..l] do
            z2:=z2+1;
            if uc[z2] = 1 then
                Add(rels,Comm(gens[j],gens[i])/rr[z]);
                z:=z+1;
            fi;
        od;
    od;

    # return
    return AgGroupFpGroup( rec(generators:=gens,relations:=rels) );
end;

#############################################################################
##
#F  CodeAgGroup( <G> ) . . . . . . . . . . .convert small ag group to integer
##
CodeAgGroup := function( G )
    local code, indices, l, mi, i, base, nt, r, j, e;

    if not IsAgGroup( G ) then
        Error("CodeAgGroup: <G> must be AgGroup\n");
    elif Size(G) > 1000 or Size(G) in [512, 768] then
        Error("CodeAgGroup: <G> has wrong size\n");
    fi;
 
    # basic structures
    l := Length( G.generators );
    indices := List( G.generators, RelativeOrderAgWord );
    mi := Maximum( indices ) - 1;
    code := 0;
    base := 1;

    # code indices of ag-series for non-p-groups
    if Length( Set( indices ) ) > 1 then
        for i in Reversed( [ 1 .. l ] ) do
            code := code + base * ( indices[ i ] - 2 );
            base := base * mi;
        od;
    fi;

    #  code which powers are not trivial and collect values into nt
    nt := [];
    for i in [ 1 .. l - 1 ] do
        r := G.generators[ i ] ^ indices[ i ];
        if r <> G.identity then
            Add( nt, r );
            code := code + base;
        fi;
        base := base * 2;
    od;

    # ... and commutators
    for i in [ 1 .. l - 1 ] do
        for j in [ i + 1 .. l ] do
            r := Comm( G.generators[ j ], G.generators[ i ] );
            if r <> G.identity then
                Add( nt, r );
                code := code + base;
            fi;
            base := base * 2;
        od;
    od;

    # code now non-trivial words
    e := Enumeration( G );
    # if gap-3.5 is used
    if not IsBound( e.number ) then
        e.number := e.NumberElement;
    fi;

    for i in nt do
        code := code + base * (e.number( e, i ) - 1 );
        base := base * Size( G );
    od;

    # and note down the size
    code := code * 1000 + Size( G );

    return code;
end;

#############################################################################
##
#F InitRandomIsomorphismChecking( <G> ) . . . . store some useful information
##
InitRandomIsomorphismChecking := function( G )
    local S, genlen, i, j, cands, orders, ppow;

    if IsPermGroup( G ) then
        S := SpecialAgGroup( AgGroup( G ) );
    else
        S := SpecialAgGroup( G );
    fi;

    # Leedham-Green subgroups
    S.lgSubgroups := [ ];
    S.lgSubElements := [ ];
    # elements of Leedham-Green subgroups of minimal order
    S.lgSubMinElem := [ ];
    genlen := Length( S.generators );

    for j in Reversed( [ 1 .. Length( S.first ) ] ) do
        i := S.first[ j ];
        S.lgSubgroups[ i ] := AgSubgroup( S, S.generators{[ i .. genlen ]}, 
                              true );
        if i <= Length( G.generators ) then
            cands := Difference( Elements( S.lgSubgroups[ i ] ),
                             Elements( S.lgSubgroups[ S.first[ j + 1 ] ] ) );
            orders := List( cands, x -> Order( S, x ) );
            ppow := Intersection( List( [ 1 .. genlen ],
                                    x -> S.weights[ i ][ 3 ] ^ x ), orders );
            S.lgSubElements[ i ] := cands{ Filtered( [ 1 .. Length(cands) ],
                                                x -> orders[ x ] in ppow ) };
            S.lgSubMinElem[ i ] := cands{ Filtered( [1 .. Length( cands ) ],
                                            x -> orders[ x ] = ppow[ 1 ] ) };
        fi;
    od;
    return S;
end;

#############################################################################
##
#F  RandomSpecialPres( <S> ) . . . . . . . . . .random special ag pres of <S>
##
##  <S> has to be the output of InitRandomIsomorphismChecking.
##
RandomSpecialPres := function( S )
    local Word, gens, rels, genlen, compser, newgens, elems, i, j, k, m;

    # construct word in abstract generators coresponding w
    Word := function( w )
        local word, d;

        word := gens[ 1 ] ^ 0;
        for d in [ 1 .. genlen ] do
            while not w in compser[ d + 1 ] do
                w := newgens[ d ] ^ -1 * w;
                word := word * gens[ d ];
            od;
        od;
        return word;
    end;
    
    # set up
    genlen  := Length( S.generators );
    gens    := [ "a", "b", "c", "d", "e", "f", "g", "h"];
    gens    := List( gens{[ 1..genlen ]}, AbstractGenerator );
    compser := [ ];
    newgens := [ ];

    compser[ genlen + 1 ] := S.lgSubgroups[ genlen + 1 ];
    # 
    for m in Reversed( [1 .. Length( S.first ) - 1 ] ) do
        i := S.first[ m ];
        j := S.first[ m + 1 ];
        elems := S.lgSubMinElem[ i ];

        # if the Leedham-Green subgroups do not have index p
        for k in Reversed( [ i + 1 .. j - 1 ] ) do
            # find the next new genenerator
            newgens[ k ] := Random( elems );
            # get the next composition subgroup
            compser[ k ] := MergedIgs( S, compser[k+1], [newgens[k]], false);
            # remaining elements for follwing generators
            elems := Difference( elems, Elements( compser[ k ] ) );
            if Length( elems ) = 0 then
                # no element of minimal order is left
                elems := Difference( S.lgSubElements[ i ],
                                     Elements( compser[ k ] ) );
            fi;
        od;

        newgens[ i ] := Random( elems );
        compser[ i ] := S.lgSubgroups[ i ];
    od;

    # find the relations corresponding to the pc-generators "newgens"
    # power relators
    rels := List( [ 1 .. genlen ], x -> gens[ x ] ^ S.weights[ x ][ 3 ] /
                                Word( newgens[ x ] ^ S.weights[ x ][ 3 ] ) );

    # commutator relations
    for i in [ 1 .. genlen - 1 ] do
        for j in [ i + 1 .. genlen ] do
            if Comm( newgens[ j ], newgens[ i ] ) <> S.identity then
                Add( rels, Comm( gens[ j ], gens[ i ] ) /
                                Word( Comm( newgens[ j ], newgens[ i ] ) ) );
            fi;
        od;
    od;
    return CodeAgGroup( SpecialAgGroup(
          AgGroupFpGroup( rec( generators := gens, relations := rels ) ) ) );
end;

#############################################################################
##
#F  RandomIsomorphismChecking( <G>, <H>, <l> ) . . . . . . . try to find isom
##
##  Try to find isomorphism between <G> and <H>. Note that <G> and <H> must
##  be given as AgGroups of size in the catalogue. If the result is "true"
##  then the groups are isomorphic. If the result is "unknown" then the
##  algorithm could not find an isomorphism. Note that in none of the cases
##  an isomorphism will be given explicitly. The integer <l> is used to
##  stop the algorithm at a certain point. The larger <l> is, the higher
##  is the probability that the result "unknown" means in fact that the
##  groups are not isomorphic. For $l = 10$ the probablity is already
##  quite high.
##
RandomIsomorphismChecking := function( G, H, l )
    local s1, s2, i, l1, l2, d1, d2, r;

    s1 := InitRandomIsomorphismChecking( G );
    s2 := InitRandomIsomorphismChecking( H );
    i  := 0;
    l1 := [ ];
    l2 := [ ];
    d1 := 0;
    d2 := 0;
    repeat

        # the first group
        r := RandomSpecialPres( s1 );
        if r in l2 then return true; fi;
        if r in l1 then d1 := d1 + 1; 
        else Add( l1, r );
        fi;

        # the second group
        r := RandomSpecialPres( s2 );
        if r in l1 then return true; fi;
        if r in l2 then d2 := d2 + 1; 
        else Add( l2, r );
        fi;

        # check the stopping criteria
        i := i + 1;
        InfoIdGroup1( d1, "    ", d2, "    nach    ", i, "\n" );
    until Minimum( d1, d2 ) >= l;
    return "unknown";
end;

#############################################################################
##
#F  IdGroupRandomTest( <G>, <list> ) . . . . . . . . . . . . . . . . . .local
##
##  <G> is the group in question and <list> is the list packed presentations
##  of the posible candidates for <G>
##
IdGroupRandomTest := function( g, c )
    local str, str1, str2, i, l, l1, l2, r;

    # prepare groups for guessing presentations
    str  := InitRandomIsomorphismChecking( AgGroup( g ) );
    str1 := InitRandomIsomorphismChecking( AgGroupCode( c[ 1 ] ) );
    str2 := InitRandomIsomorphismChecking( AgGroupCode( c[ 2 ] ) );

    # init lists of (numerical coded) presenatations
    l  := [ ];
    l1 := [ ];
    l2 := [ ];

    # repeat until two identical presentations are found
    repeat
        r := RandomSpecialPres( str );
        if not r in l then Add( l, r ); fi;

        r := RandomSpecialPres( str1 );
        if r in l then return c[ 1 ]; fi;
        if not r in l1 then Add( l1, r ); fi;

        r := RandomSpecialPres( str2 );
        if r in l then return c[ 2 ]; fi;
        if not r in l2 then Add( l2, r ); fi;

        InfoIdGroup1( Length(l), " ", Length(l1), " ", Length(l2), "\n" );
    until false;
end;

#############################################################################
##
#F  IdGroupSpecialFp( <G>, <integer> ). . . . . . . . . . . . . . . . . local
##
##  Compute a special fingerprint.
##
IdGroupSpecialFp := function( g, i )
    local p, S, classbound, stanpres;

    p := [ 2, 3, 5, 2 ];
    S := SylowSubgroup( g, p[ i ] );

    if i in [ 1 .. 3 ] then
        # get the standard presenatation of the relevant sylow subgroup
        if StandardPresentation = 0 then
            RequirePackage( "anupq" );
        fi;
        classbound := [ 8, 6, 4 ];
        InfoIdGroup1( "#I  IdGroupSpecialFp: anu-pq ist started\n" );
        stanpres := StandardPresentation( FpGroup( S ), p[ i ],
                                         "ClassBound", classbound[ i ] );
        return QuoInt( CodeAgGroup( stanpres ), 1000 ) mod 1000;

    else # i = 4
        # investigate the operation of g on its 2-sylow subgroup
        # this fingerprint differs 2 frattinifree groups of type 2^6:7
        return Sum( List( Orbits( g, Elements( S ) ),
                          x -> Size( Subgroup( g, x ) ) ) );
    fi;
end;


#############################################################################
##
#F  EvalFpCoc( g, coc, desc ) . . . . . . . . . . . . . . . . . . . . . local
##
EvalFpCoc := function( g, coc, desc )
    local powers, exp, targets, result, i, j, g1, g2, fcd4, pos;

    if desc[ 1 ] = 1 then
        # test, if g^i in cl(g)
        return List( coc[ desc[ 2 ] ],
                     function( x )
                     if x[ 1 ] ^ desc[ 3 ] in x then return 1; fi; return 0;
                     end );

    elif desc[ 1 ] = 2 then
        # test, if cl(g) is root of cl(h)
        exp := QuoInt( Order( g, coc[ desc[ 2 ] ][ 1 ][ 1 ] ),
                       Order( g, coc[ desc[ 3 ] ][ 1 ][ 1 ] ) );
        powers := Flat( coc[ desc[ 3 ] ] );
        return List( coc[ desc[ 2 ] ],
                     function(x)
                     if x[ 1 ] ^ exp in powers then return 1; fi; return 0;
                     end );

    elif desc[ 1 ] = 3 then
        # test, if cl(g) is power of cl(h)
        exp := QuoInt( Order( g, coc[ desc[ 3 ] ][ 1 ][ 1 ] ),
                       Order( g, coc[ desc[ 2 ] ][ 1 ][ 1 ] ) );
        # just one representative for each class of power-candidates
        powers := List( coc[ desc[ 2 ] ], x -> x[ 1 ] );
        result := List( powers, x -> 0 );
        for i in List( Flat( coc[ desc[ 3 ] ] ), x -> x ^ exp ) do
            for j in [ 1 .. Length( powers ) ] do
                if i = powers[ j ] then
                    result[ j ] := result[ j ] + 1;
                fi;
            od;
        od;
        return result;

    else 
        # test how often the word [ a, b ] * a^2 is hit
        targets := List( coc[ desc[ 2 ] ], x -> x[ 1 ] );
        result := List( targets, x -> 0 );
        fcd4 := Flat( coc[ desc[ 4 ] ] );
        for g1 in Flat( coc[ desc[ 3 ] ] ) do
            for g2 in fcd4 do
                if desc[ 1 ] = 4 then 
                    pos := Position( targets, Comm( g1, g2 ) * g1 ^ 2 );
                else 
                # desc[ 1 ] = 5
                    pos := Position( targets, Comm( g1, g2 ) * g1 ^ 3 );
                fi;
                if pos <> false then
                    result[ pos ] := result[ pos ] + 1;
                fi;
            od;
        od;
        return result;
    fi;
end;

#############################################################################
##
#F  IdSmallGroup( <G> ) . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| <= 1000 and |G| not in {256, 
##  512, 768} and |G| consists of more than 3 primes.
##
IdSmallGroup := function( G )
    local level, branch, indices, fp, elementOrders, setElOrders, l, L, i, j,
          collElOrders, coc, desc, pos, filename, ldesc, Pack, sfp, newcls,
          classes, classtyps, sclasstyps;

    # packs information from a list - with posible loos of information
    Pack := function( list )
        local r, i;

        if Length( list ) = 0 then 
            return 0;
        fi;
        list := Flat( list );
        r := list[ 1 ] mod 99661;
        for i in list{[ 2 .. Length( list ) ]} do
            r := (r * 10 + i ) mod 99661;
        od;
        return r;
    end;

    # set up
    level := 1;
    branch := IdGroupTree;
    indices := [ ];
    l := "abcdefghijklmnopqrstuvwxyz";
	# all letters in file names have to be small case
	# case-sensitivity is not needed for the used names
    # L := "ABCDEFGHIJK";
    L := "abcdefghijkl";


    # main loop
    while not IsInt( branch ) do
        
        if level = 1 then
            fp := Size( G );

        elif level = 2 then 
            fp := Pack( List( DerivedSeries( G )
                      {[ 2 .. Length( DerivedSeries( G ) ) ]}, Size ) );

        elif level = 3 then 
            if IsAbelian( G ) then 
                fp := Pack( AbelianInvariants( G ) );
            else 
                elementOrders := List( Elements( G ), x -> Order( G, x ) );
                setElOrders := Set( elementOrders );
                fp := Pack( setElOrders{[ 2 .. Length( setElOrders ) ]} );
            fi;

        elif level = 4 then 
            collElOrders := Collected( elementOrders );
            fp := Pack( List( collElOrders{[ 2 .. Length( collElOrders ) ]},
                              x -> x[ 2 ] ) );

        elif level = 5 then 
            # on level 5 the tests on conjugacy classes start
            classes := Orbits( G, Elements( G ) );
            classtyps := List( classes,
                               x -> [ Order( G, x[ 1 ] ), Length( x ) ] );
            sclasstyps := Set( classtyps );
            # coc is   Clusters Of Conjugacy   classes
            coc := List( sclasstyps, x-> [ ] );
            for i in [ 1 .. Length( sclasstyps ) ] do
                for j in [ 1 .. Length( classes ) ] do
                    if sclasstyps[ i ] = classtyps[ j ] then
                        Add( coc[ i ], classes[ j ] );
                    fi;
                od;
            od;
            
            fp := Pack( List( coc{[ 2 .. Length( coc ) ]},
                              x -> [ Length( x[ 1 ] ), Length( x ) ] ) );

        elif not IsList( branch.desc ) then
            if branch.desc = 0 then
                # this special case could apear only on level >= 6
                fp := IdGroupRandomTest( G, branch.fp );
            
            else
                # use a special fingerprint, apears just for level >= 6
                fp := IdGroupSpecialFp( G, branch.desc );
            fi;

        else
            # usuall case for level >= 6
            for desc in branch.desc do
                # reconstruct orignial description list of the test
                ldesc := [ desc mod 1000 ];
                desc := QuoInt( desc, 1000 );
                while desc > 0 do
                    Add( ldesc, desc mod 100 );
                    desc := QuoInt( desc, 100 );
                od;
                desc := Reversed( ldesc );
                
                # evaluate the test
                fp := EvalFpCoc( G, coc, desc );

                # split up clusters of classes acording to the result of test
                sfp := Set( fp );
                newcls := List( sfp, x-> [ ] );
                for i in [ 1 .. Length( sfp ) ] do
                    for j in [ 1 .. Length( fp ) ] do
                        if sfp[ i ] = fp[ j ] then
                            Add( newcls[ i ], coc[ desc[ 2 ] ][ j ] );
                        fi;
                    od;
                od;
                coc := Concatenation( coc{[ 1 .. desc[ 2 ] -1 ]}, newcls,
                                   coc{[ desc[ 2 ] + 1 .. Length( coc ) ]} );
            od;

            # make fingerprint independ from the rowing of conj-classes
            fp := Pack( Collected( fp ) );
        fi;

        pos := Position( branch.fp, fp );
        if pos = false then
            Error( "IDGROUP: fatal Error. Please mail group to\n",
                   "Hans-Ulrich.Besche@math.rwth-aachen.de" );
        fi;
        Add( indices, pos );

        # load required branch of 'IdGroupTree' if it is not in memory
        if not IsBound( branch.next[ pos ] ) then
            filename := "";
            if Size( G ) < 10 then 
                Append( filename, "00" );
            elif Size( G ) < 100 then 
                Add( filename, '0' );
            fi;
            Append( filename, String( Size( G ) ) );
            for i in indices{[ 2 .. Length( indices ) ]} do
                if i > 26 then 
                    Add( filename, L[ QuoInt( i - 1, 26 ) ] );
                fi;
                Add( filename, l[ ( i - 1 ) mod 26 + 1 ] );
            od;
            filename := ConcatenationString( "idlib/id", 
                filename{[ 1 .. Minimum( Length( filename ) - 1, 6 ) ]}, ".",
                filename{[ Minimum( Length( filename ), 7 ) ..
                Length( filename ) ]} );
            ReadSml( filename );
            InfoIdGroup1( "#I IdSmallGroup reads ", filename,
                     " SIZE( IdGroupTree ) = ", SIZE( IdGroupTree ) , "\n" );
        fi;

        branch := branch.next[ pos ];
        level := level + 1;
    od;

    # branch is now a integer
    return branch;
end;

#############################################################################
##
#F  IdP1Q1R1Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p*q*r.
##
IdP1Q1R1Group := function( G )
    local  n, p, q, r, s1, s2, s3, typ, PO, A, C, AC, x, s, y, B, BC, id;

    # get primes
    n := Size(G);
    p := Factors(n);
    r := p[3];
    q := p[2];
    p := p[1];

    # compute the sylow subgroups
    s1  := SylowSubgroup( G, r );
    s2  := SylowSubgroup( G, q );
    s3  := SylowSubgroup( G, p );

    if IsAbelian( G ) then
        typ := "pqr";

    elif IsNormal( G, s3 ) then
        typ := "Dqr x p";

    elif not IsAbelian( Closure( s1, s2 ) ) then
        typ := "Hpqr";

    elif IsAbelian( Closure( s1, s3 ) ) then
        typ := "Dpq x r";

    elif IsAbelian( Closure( s2, s3 ) ) then
        typ := "Dpr x q";

    else 

        # find <A> and <C>
        C  := SylowSubgroup(G,p).generators[1];
        A  := SylowSubgroup(G,r).generators[1];
        AC := A^C;
        PO := 0;
        s  := 1;
        while AC <> PO  do
            s := s+1;
            if s mod r <> 1 and s^p mod r = 1  then
                PO := A^s;
            fi;
        od;

        # correct <C>
        x := First( [2..r-1], t -> t mod r <> 1 and t^p mod r = 1 );
        s := LogMod( x, s, r );
        C := C^s;

        # now find <B>
        B  := SylowSubgroup(G,q).generators[1];
        BC := B^C;
        PO := 0;
        s  := 1;
        while BC <> PO  do
            s := s+1;
            if s mod r <> 1 and s^p mod r = 1  then
                PO := B^s;
            fi;
        od;

        # and find <s>
        y := First( [2..q-1], t -> t mod q <> 1 and t^p mod q = 1 );
        s := LogMod( s, y, q );

        typ := s; # the typ is Gpqr( s )
    fi;

    # find the types existing for Size( g ) and count up id
    id := 1;

    if typ = "Hpqr" then
        return id;
    fi;
    if r mod (p*q) = 1 then id := id + 1; fi;

    if typ = "Dqr x p" then
        return id;
    fi;
    if r mod q = 1 then id := id + 1; fi;

    if typ = "Dpq x r" then
        return id;
    fi;
    if q mod p = 1 then id := id + 1; fi;

    if typ = "Dpr x q" then
        return id;
    fi;
    if r mod p = 1 then id := id + 1; fi;

    if IsInt( typ ) then 
        # g is Gpqr( typ )
        return  id - 1 + typ;
    fi; 
    if ( r mod p = 1 ) and ( q mod p = 1 ) then id := id + p - 1; fi;

    # remaining typ is pqr
    return id;
end;

#############################################################################
##
#F  IdP2Q1Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p^2*q.
##
IdP2Q1Group := function( G )
    local  n, p, q, id, s1, s2, typ;

    # get primes
    n := Size(G);
    p := Factors(n);
    q := p[3];
    p := p[1];

    # compute the sylow subgroups
    s1  := SylowSubgroup( G, q );
    s2  := SylowSubgroup( G, p );

    if IsAbelian( G ) then
        if IsCyclic( s2 ) then
            typ := "p2 x q";
        else
            typ := "p x pq";
        fi;

    elif IsElementaryAbelian( s2 ) then
        if n = 12 and IsNormal( G, s2 ) then
            typ := "a4";
        else
            typ := "Dpq x p";
        fi;

    elif not IsTrivial( Centralizer( s2, s1 ) ) then
        typ := "Gp2q";
    else 
        typ := "Hp2q";
    fi;

    id := 1;

    if typ = "Gp2q" then
        return id; fi;
    if q mod p = 1 then id := id + 1; fi;

    if typ = "p2 x q" then
        return id; fi;
    id := id + 1;

    if typ = "Hp2q" then
        return id; fi;
    if q mod (p*p) = 1 then id := id + 1; fi;

    if typ = "a4" then
        return 3; fi;
    if n = 12 then id := id + 1; fi;

    if typ = "Dpq x p" then
        return id; fi;
    if q mod p = 1 then id := id + 1; fi;

    return id;
end;

#############################################################################
##
#F  IdP1Q2Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p*q^2.
##
IdP1Q2Group := function( G )

    local  n, p, q, s1, s2, typ, lat, nor, non, x, s, A, B, C, AC, BC, PO, id;

    # get primes
    n := Size(G);
    p := Factors(n);
    q := p[3];
    p := p[1];

    # compute the sylow subgroups
    s1  := SylowSubgroup( G, q );
    s2  := SylowSubgroup( G, p );

    if IsAbelian( G ) then
        if IsCyclic( s1 ) then
            typ := "p x q2";
        else
            typ := "pq x q";
        fi;

    elif ( p <> 2 ) and ( (q+1) mod p = 0 ) then
        typ := "Npq2";

    elif IsCyclic( s1 ) then
        typ := "Mpq2";

    elif not IsTrivial( Centralizer( s1, s2 ) ) then
        typ := "Dpq x q";

    else    
        lat := List( ConjugacyClassesSubgroups(s1), Representative );
        lat := Filtered( lat, t -> Size(t) = q );
        nor := [];
        non := [];
        x   := 1;
        while x <= Length(lat) and 0 = Length(non) and Length(nor) < 3 do
            if IsNormal( G, lat[x] )  then
                Add( nor, lat[x] );
            else
                Add( non, lat[x] );
            fi;
            x := x + 1;
        od;
        if 0 = Length(non) and 2 < Length(nor)  then
            typ := "Kpq2";
        else
            while x <= Length(lat) and Length(nor) < 2  do
                if IsNormal( G, lat[x] )  then
                    Add( nor, lat[x] );
                fi;
                x := x + 1;
            od;
            A  := nor[1].generators[1];
            B  := nor[2].generators[1];
            C  := s2.generators[1];
            AC := A^C;
            x  := 1;
            PO := 0;
            while PO <> AC  do
                x := x+1;
                if x^p mod q = 1  then
                    PO := A^x;
                fi;
            od;
            BC := B^C;
            s  := 1;
            PO := 0;
            while s < q and PO <> BC  do
                s := s+1;
                if s mod p <> 0 and s mod p <> 1  then
                    PO := B^((x^s) mod q);
                fi;
            od;
            s := s mod p;
            if ((1/s) mod p) < s  then s := (1/s) mod p;  fi;
            typ := s;
        fi;
    fi;

    id := 1;

    if typ = "Mpq2" then
        return id;
    fi;
    if q mod p = 1 then id := id + 1; fi;

    if typ = "p x q2" then
        return id;
    fi;
    id := id + 1;

    if typ = "Dpq x q" then
        return id;
    fi;
    if q mod p = 1 then id := id + 1; fi;

    if typ = "Npq2" then
        return id;
    fi;
    if (p<>2) and ( (q+1) mod p = 0 ) then id := id + 1; fi;

    if typ = "Kpq2" then
        return id;
    fi;
    if q mod p = 1 then id := id + 1; fi;

    if IsInt( typ ) then
        # g is Lpq2(typ)
        return id - 1 +
            Position( Filtered( [ 2 .. p-1 ], x -> x <= 1/x mod p), typ );
    fi;
    if ( p <> 2 ) and ( q mod p = 1 ) then id := id + ( p - 1 ) / 2; fi;

    return id;
end;

#############################################################################
##
#F  IdP1Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p.
##
IdP1Group := function( G )

    return 1;
end;

#############################################################################
##
#F  IdP2Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p^2.
##
IdP2Group := function( G )

    if IsCyclic( G ) then
        return 1;
    fi;
    return 2;
end;

#############################################################################
##
#F  IdP3Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p^3.
##
IdP3Group := function( G )

    if IsAbelian( G ) then
        if IsCyclic( G ) then
            return 1;
        elif IsElementaryAbelian( G ) then
            return 5;
        fi;
        return 2;

    else 
        if Size( G ) = 8 then
            if Length(Filtered(Elements(G),x->Order(G,x)=2)) = 1 then
                return 4;
            fi;
            return 3; 
        fi;

        if Maximum( List( Generators(G), x -> Order(G,x) ) ) =
                    Factors( Size( G ) )[ 1 ] then
            return 3;
        fi;    
        return 4;
    fi;
end;



#############################################################################
##
#F  IdP1Q1Group( <G> ). . . . . . . . . . . . . . . . . . . . . . . . . local
##
##  Compute identification for <G> where |G| = p*q.
##
IdP1Q1Group := function( G )
    local typ, p, q, id;

    p := Factors( Size( G ) );
    q := p[ 2 ];
    p := p[ 1 ];

    if IsAbelian( G ) then
        typ := "pq";
    else 
        typ := "Dpq";
    fi;

    id := 1;

    if typ = "Dpq" then 
        return id;
    fi;
    if q mod p = 1 then id := id + 1; fi;

    # typ is pq
    return id;
end;

#############################################################################
##
#F  IdGroup( <G> ) . . . . . . . . . . . . . . . .identify group, if possible
## 
##  It will be possible, if |G| <= 1000 and |G| <> 256, 512, 768 or if
##  |G| is a product of at most 3 primes. 
##  G should be an AgGroup or a PermGroup
##
IdGroup := function( G )
    local size, primes, sprimes, result;

    # set up group
    if not IsAgGroup( G ) and not IsPermGroup( G ) then
        Error("IdGroup: <G> must be AgGroup or PermGroup\n");
    fi;

    # set up
    size := Size( G );
    primes := Factors( size );
    sprimes := Set( primes );
    
    # catch the trivial case
    if size = 1 then
        result := 1;

    # p-groups with size <= p ^ 3
    elif ( Length( sprimes ) = 1 ) and ( Length( primes ) <= 3 ) then
        if Length( primes ) = 1 then
            result := IdP1Group( G );
        elif Length( primes ) = 2 then
            result := IdP2Group( G );
        else 
            result := IdP3Group( G );
        fi;

    # pq-groups of typ pq, ppq and pqq
    elif ( Length( sprimes ) = 2 ) and ( Length( primes ) <= 3 ) then
        if Length( primes ) = 2 then 
            result := IdP1Q1Group( G );
        else 
            if primes[ 1 ] = primes[ 2 ] then
                result := IdP2Q1Group( G );
            else
                result := IdP1Q2Group( G );
            fi;
        fi;

    # pqr-groups of typ pqr
    elif ( Length( sprimes ) = 3 ) and ( Length( primes ) = 3 ) then
        result := IdP1Q1R1Group( G );

    # if Size( G ) is not as above usuall test restrictions on size
    elif ( size > 1000 ) or ( size in [ 256, 512, 768 ] ) then
        Error( "IdGroup: Size(G) restricted to 1000, except 256, 512, 768" );

    # final case
    else
        result := IdSmallGroup( G );
    fi;

    # return
    return [ Size( G ), result ];
end;
