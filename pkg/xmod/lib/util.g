###############################################################################
##
##  util.g                       for GAP 3.4                   version 10/ 1/97
##
###############################################################################
##
#A  util.g                       GAP library                      Chris Wensley
## 
#Y  Copyright
##
#H  $Log: util.g,v $
#H  Revision 1.1  1997/03/27 13:35:53  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

##############################################################################
##
#F  GroupOpsZeroMorphism                  maps G to the identity subgroup of H
##

GroupOpsZeroMorphism := function( G, H )

    local  genG, imH, z;

    if not ( IsGroup( G ) and IsGroup( H ) and IsBound( H.identity ) ) then
        Error( "argument must be a pair of groups with identity" );
    fi;
    genG := G.generators;
    imH := List( genG, x -> H.identity );
    z := GroupHomomorphismByImages( G, H, genG, imH );
    return z;
end;

##############################################################################
##
#F  GroupOpsInclusionMorphism                 includes subgroup H in a group G
##

GroupOpsInclusionMorphism := function( H, G )

    local  genH, hom;

    genH := H.generators;
    if ( G = H ) then
        hom := IdentityMapping( G );
        hom.generators := genH;
        hom.genimages := genH;
    elif ( Length( genH ) > 0 ) then
        hom := GroupHomomorphismByImages( H, G, genH, genH );
    else
        hom := GroupHomomorphismByImages( H, G, [H.identity], [G.identity] );
    fi;
    return hom;
end;

###############################################################################
##
#F  GroupOpsIsAutomorphismGroup           test whether a group of automorphisms
##

GroupOpsIsAutomorphismGroup := function( G )

    local  ok, id, g, src, ops;

    id := G.identity;
    ok := ( IsRec( id ) and IsBound( id.operations ) );
    if ok then
        ops := id.operations;
    else
        return false;
    fi;
    ok := ( ops = IdentityGroupHomomorphismOps );
    if not ok then
        Print( "G.identity.operations = ", ops, "\n" );
        return false;
    fi;
    if IsBound( G.morphismDomain ) then
        src := G.morphismDomain;
    else
        src := id.source;
    fi;
    for g in G.generators do
        if not ( ( g.source = src ) and ( g.range = src )
                                  and IsAutomorphism( g ) ) then
            ok := false;
        fi;
    od;
    if ( ok and not IsBound( G.morphismDomain ) ) then
        G.morphismDomain := src;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsAutomorphismPair     test whether an automorphism group / perm group pair
##
##  Problem (5/7/96):  Testing that a2p, p2a are isomorphisms takes too long
##                     and so has been removed from this function.

IsAutomorphismPair := function( R )

    if not IsRec( R ) then
        return false;
    fi;
    if IsBound( R.isAutomorphismPair ) then
        return R.isAutomorphismPair;
    fi;
    if (   ( not (      IsBound( R.auto )
                    and IsAutomorphismGroup( R.auto ) ) ) 
        or ( not (      IsBound( R.perm )
                    and IsPermGroup( R.perm ) ) )
        or ( not (      IsBound( R.p2a )
                    # and IsIsomorphism( R.p2a )
                    and ( R.p2a.source = R.perm ) 
                    and ( R.p2a.range = R.auto ) ) )
        or ( not (      IsBound( R.a2p )
                    # and IsIsomorphism( R.a2p )
                    and ( R.a2p.source = R.auto ) 
                    and ( R.a2p.range = R.perm ) ) ) ) then
        return false;
    fi;
    R.IsAutomorphismPair := true;
    return true;
end;

##############################################################################
##
#F  AutomorphismPair                       set up AutoGroup <-> PermGroup pair
##

AutomorphismPair := function( A )

    local  G, A, a, genA, genG, eG1, eG2, new, num1, num2, eG1, eG2,
           ima, P, genP, a2p, p2a, imag, pair;

    if not IsAutomorphismGroup( A ) then
        Error( "<A> must be a group of automorphisms" );
    fi;
    G := A.morphismDomain;
    a := A.identity;
    genA := A.generators;
    genG := G.generators;
    # determine the closure of genG under A
    eG2 := genG;
    num1 := 0;
    num2 := Length( genG );
    new := genG;
    while ( num1 <> num2 ) do
        num1 := num2;
        eG1 := eG2;
        for a in genA do
            ima := List( new, g -> Image( a,g ) );
            eG2 := Union( eG2, ima );
        od;
        new := Difference( eG2, eG1 );
        num2 := Length( eG2 );
    od;
    # construct perm group P using action of A on this closure
    P := Operation( A, eG2 );
    genP := P.generators;
    if IsBound( G.name ) then
        if not IsBound( A.name ) then
            if ( IsBound( G.automorphismGroup ) 
                 and ( A = G.automorphismGroup ) ) then
                A.name := Concatenation( "Aut(", G.name, ")" );
            else
                A.name := Concatenation( "SubAut(", G.name, ")" );
            fi;
        fi;
        P.name := Concatenation( "Perm", A.name );
    fi;
    a2p := OperationHomomorphism( A, P );
    imag := List( A.generators, a -> Image( a2p, a ) );
    a2p.genimages := imag;
    p2a := GroupHomomorphismByImages( P, A, imag, genA );
    pair := rec( );
    pair.auto := A;
    pair.perm := P;
    pair.a2p := a2p;
    pair.p2a := p2a;
    pair.isAutomorphismPair := true;
    A.automorphismPair := pair;
    return pair;
end;

##############################################################################
##
#F  GroupOpsAutomorphismPermGroup             AutomorphismPermGroup for groups
##

GroupOpsAutomorphismPermGroup := function( G )

    local  A, pair;

    A := AutomorphismGroup( G );
    pair := AutomorphismPair( A );
    return pair.perm;
end;

##############################################################################
##
#F  GroupOpsInnerAutomorphismGroup        group of inner automorphism of N < G
##  

GroupOpsInnerAutomorphismGroup := function( arg )

    local  nargs, N, G, id, genG, genN, genA, A, g, conj, name;

    nargs := Length( arg );
    if ( ( nargs < 1 ) or ( nargs > 2 ) ) then
        Error( "\nUsage:  InnerAutomorphismGroup( G [, N] );\n" );
    fi;
    G := arg[1];
    if IsBound( G.innerAutomorphismGroup ) then
        return G.innerAutomorphismGroup;
    fi;
    if ( nargs = 2 ) then
        N := arg[2];
        if not IsNormal( G, N ) then
            Error( "Second parameter must be normal subgroup of the first\n" );
        fi;
    else
        N := G;
    fi;
    genG := G.generators;
    genN := N.generators;
    id := InclusionMorphism( N, N );
    genA := [ ];
    for g in genG do
        conj := InnerAutomorphism( N, g );
        conj.generators := genN;
        conj.genimages := List( genN, n -> Image( conj, n ) );
        Add( genA, conj );
    od;
    A := Group( genA, id );
    A.morphismDomain := G;
    A.isAutomorphismGroup := true;
    return A;
end;

###############################################################################
##
#F  EndomorphismClasses                         finds all homomorphisms  G -> G
##

EndomorphismClasses := function( arg )

    local  nargs, valid, case, switch, oldcase, G, normG, rcosG, N, Nnum,
           quotG, genG, Qnum, lat, reps, Rnum, Q, R, iso, L, im, phi,
           nat, proj, comp, normal, cosets, conj, Cnum, g, gim, zero,
           i, j, k, l, Lnum, aut, idgp, E, Enum, classes, class, name, ok;

    nargs := Length( arg );
    valid := ( nargs < 3 );
    G := arg[1];
    if not ( IsBound( G.isPermGroup ) and G.isPermGroup ) then
        Error( "parameter must be a permutation group" );
    fi;
    if ( nargs = 2 ) then
        case := arg[2];
        valid := ( valid and ( IsInt( case ) and ( case in [0..3] ) ) );
        if not valid then
            Print( "\nUsage:  EndomorphismClasses( G [, case] );\n" );
            Print( " choose  case = 1  to include automorphisms & zero,\n" );
            Print( "default  case = 2  to exclude automorphisms & zero,\n" );
            Print( "         case = 3  when  N meet H  trivial,\n" );
            return false;
        fi;
    else
        case := 2;
    fi;
    switch := false;
    if IsBound( G.endomorphismClasses ) then
        E := G.endomorphismClasses;
        if E.areNonTrivial then
            oldcase := 2;
        else
            oldcase := 1;
        fi;
        # case may be 0 if called from EndomorphismImages
        if ( case = 0 ) then
            case := oldcase;
        fi;
        switch := not ( case = oldcase );
        if not switch  then
            return E;
        elif ( case = 2 ) then
            Enum := Length( E.classes );
            E.classes := Sublist( E.classes, [2..Enum-1] );
            E.areNonTrivial := true;
            return E;
        fi;
    fi;
    genG := G.generators;
    idgp := Subgroup( G, [ ] );
    if IsBound( G.name ) then
        name := G.name;
    else
        name := "G";
    fi;

    # determine non-monomorphic endomorphisms by kernel
    normG := NormalSubgroups( G );
    Nnum := Length( normG );
    for i in [2..Nnum-1] do
        normG[i].name := Concatenation( name, ".N", String( i ) );
    od;
    normG[1].name := Concatenation( name, ".id" );
    Nnum := Length( normG );
    if ( XModPrintLevel > 1 ) then
        Print( name, " has ", Nnum, " normal subgroups.\n" );
    fi;
    lat := Lattice( G );
    reps := List( lat.classes, class -> class.representative );
    Rnum := Length( reps );
    if ( XModPrintLevel > 1 ) then
        Print( name, " has ", Rnum, " subgroup classes.\n" );
    fi;    
    for i in [2..Rnum-1] do
        reps[i].name := Concatenation( name, ".H", String( i ) );
    od;
    rcosG := List( normG, N -> RightCosets( G, N ) );
    quotG := List( rcosG, cosets -> Operation( G, cosets, OnRight ) );
    Qnum := Length( quotG );
    for i in [2..Nnum-1] do
        quotG[i].name := Concatenation( name, ".Q", String( i ) );
    od;
    normal := [ true ];
    for j in [2..(Rnum-1)] do
        Add( normal, ( lat.classes[j].normalizer = G ) );
    od;
    Add( normal, true );
    # the zero class is needed here if switching
    Q := Group( () );
    zero := rec( );
    zero.quotient := Q;
    zero.projection := ZeroMorphism( G, Q );
    zero.autoGroup := Group( InclusionMorphism( Q, Q ) );
    zero.rangeNumber := 1;
    zero.isomorphism := ZeroMorphism( Q, idgp );
    zero.conj := [ () ];

    classes := [ ];
    if ( case = 1 ) then
        # all automorphisms are endomorphisms
        aut := AutomorphismGroup( G );
        class := rec( );
        class.quotient := G;
        class.projection := InclusionMorphism( G, G );
        class.autoGroup := aut;
        class.rangeNumber := Rnum;
        class.isomorphism := InclusionMorphism( G, G );
        class.conj := [ () ];
        Add( classes, class );
        if switch then
            E.areNonTrivial := false;
            E.classes := Concatenation( [ class ], E.classes, [ zero ] );
            return E;
        fi; 
    fi;
    for i in [2..(Qnum-1)] do
        N := normG[i];
        Q := quotG[i];
        proj := OperationHomomorphism( G, Q );
        aut := AutomorphismGroup( Q );
        for j in [2..(Rnum-1)] do
            R := reps[j];
            if ( case < 3 ) then
                ok := true;
            else
                ok := ( Intersection( N, R ) = idgp );
            fi;
            if ok then
                iso := IsomorphismGroups( Q, R );
            fi;
            if ( ok and not ( iso = false ) ) then
                # (unnecessary?) check that this gives a homomorphism
                comp := proj * iso;
                im := List( genG, x -> Image( comp, x ) );
                phi := GroupHomomorphismByImages( G, R, genG, im );
                if not IsHomomorphism( phi ) then
                    Error( "phi not a homomorphism" );
                fi;
                class := rec( );
                class.quotient := Q;
                class.projection := proj;
                class.autoGroup := aut;
                class.rangeNumber := j;
                class.isomorphism := iso;
                if not normal[j] then
                    cosets := Cosets( G, lat.classes[j].normalizer );
                    conj := List( cosets, C -> C.smallest );
                    class.conj := conj;
                else
                    class.conj := [ () ];
                fi;
            Add( classes, class );
            fi;
        od;
    od;
    if ( case = 1 ) then
        # finally, add the zero map
        Add( classes, zero );
    fi;

    E := rec( );
    E.isDomain := true;
    E.isEndomorphismClasses := true;
    E.areNonTrivial := ( case > 1 );
    E.intersectionFree := (case = 3 );
    E.classes := classes;
    E.group := G;
    E.latticeLength := Rnum;
    E.latticeReps := reps;
    G.endomorphismClasses := E;
    return E;
end;

###############################################################################
##
#F  EndomorphismImages           finds all homomorphisms  G -> G  by converting
##                                 EndomorphismClasses into a list of genimages
##

EndomorphismImages := function( arg )

    local  nargs, valid, case, G, genG, Q, R, k, iso, proj, comp, im,
           aut, autos, Anum, rho, psi, phi, conj, Cnum, g, gim, l, L,
           LR, Lnum, latlen, Rnum, reps, E, Enum, c, class, classes;

    nargs := Length( arg );
    valid := ( nargs < 3 );
    G := arg[1];
    if not ( IsBound( G.isPermGroup ) and G.isPermGroup ) then
        Error( "parameter must be a permutation group" );
    fi;
    if ( nargs = 2 ) then
        case := arg[2];
        valid := ( valid and ( IsInt( case ) and ( case in [1..3] ) ) );
        if not valid then
            Print( "\nUsage:  EndomorphismImages( G [, case] );\n" );
            Print( " choose  case = 1  to include automorphisms & zero,\n" );
            Print( "default  case = 2  to exclude automorphisms & zero,\n" );
            Print( "         case = 3  when  N meet H  trivial,\n" );
            return false;
        fi;
    else
        case := 0;
    fi;
    genG := G.generators;
    E := EndomorphismClasses( G, case );
    classes := E.classes;
    latlen := E.latticeLength;
    reps := E.latticeReps;
    Enum := Length( classes );
    L := [ ];
    for c in [1..Enum] do
        class := classes[c];
        Rnum := class.rangeNumber;
        R := reps[Rnum];
        conj := class.conj;
        Cnum := Length( conj );
        Q := class.quotient;
        proj := class.projection;
        iso := class.isomorphism;
        comp := proj * iso;
        aut := class.autoGroup;
        autos := Elements( aut );
        Anum := Length( autos );
        LR := [ ];
        for k in [1..Anum] do
            rho := autos[k];
            psi := proj * rho * iso;
            im := List( genG, x -> Image( psi, x ) );
            Add( LR, im );
        od;
        Lnum := Length( LR );
        for k in [2..Cnum] do
            g := conj[k];
            for l in [1..Lnum] do
                im := LR[l];
                gim := List( im, x -> x^g );
                Add( LR, gim );
            od;
        od;
        L := Concatenation( L, LR ); 
    od;
    return L;
end;

###############################################################################
##
#F  IdempotentImages           finds all idempotent homs  G -> G  by converting
##                                 EndomorphismClasses into a list of genimages
##

IdempotentImages := function( arg )

    local  valid, nargs, case, G, genG, Q, R, k, iso, proj, comp, im,
           aut, autos, Anum, rho, psi, phi, conj, Cnum, g, gim, l, 
           L, LI, Lnum, E, Enum, Rnum, reps, latlen, c, class, classes, ok;

    nargs := Length( arg );
    valid := ( nargs < 3 );
    G := arg[1];
    if not ( IsBound( G.isPermGroup ) and G.isPermGroup ) then
        Error( "first parameter must be a permutation group" );
    fi;
    if ( nargs = 2 ) then
        case := arg[2];
        valid := ( valid and ( IsInt( case ) and ( case in [1..4] ) ) );
    fi;
    if not valid then
        Print( "\nUsage: IdempotentImages( G [, case] );\n" );
        Print( "where  case = 1  for ALL idempotent images,\n" );
        Print( "       case = 2  for all non-trivial images,\n" );
        Print( "       case = 3  for case 2 & one group per conj class,\n" );
        Print( "       case = 4  for case 3 and sorted into images.\n" );
        return false;
    fi;
    genG := G.generators;
    if ( case = 1 ) then
        E := EndomorphismClasses( G, 1 );
    else
        E := EndomorphismClasses( G, 2 );
    fi;    
    classes := E.classes;
    Enum := Length( classes );
    latlen := E.latticeLength;
    reps := E.latticeReps;
    if ( case = 4 ) then
        L := List( [1..latlen], r -> [ ] );
    else
        L := [ ];
    fi;
    for c in [1..Enum] do
        class := classes[c];
        Rnum := class.rangeNumber;
        R := reps[Rnum];
        conj := class.conj;
        Cnum := Length( conj );
        Q := class.quotient;
        proj := class.projection;
        iso := class.isomorphism;
        comp := proj * iso;
        aut := class.autoGroup;
        autos := Elements( aut );
        Anum := Length( autos );
        LI := [ ];
        for k in [1..Anum] do
            rho := autos[k];
            psi := proj * rho * iso;
            im := List( genG, x -> Image( psi, x ) );
            if ( psi * psi = psi ) then
                Add( LI, im );
            fi;
        od;
        Lnum := Length( LI );
        if ( case < 3 ) then
            for k in [2..Cnum] do
                g := conj[k];
                for l in [1..Lnum] do
                    im := LI[l];
                    psi := GroupHomomorphismByImages( G, G, genG, im );
                    gim := List( im, x -> Image( psi, x^(g^-1) )^g );
                    phi := GroupHomomorphismByImages( G, G, genG, gim );
                    Add( LI, gim );
                od;
            od;
        fi;
        if ( case = 4 ) then
            L[Rnum] := Concatenation( L[Rnum], LI );
        else
            L := Concatenation( L, LI ); 
        fi;
    od;
    return L;
end;

###############################################################################
##
#F  IsFpPair                         test whether an fp group / perm group pair
##

IsFpPair := function( R )

    if not IsRec( R ) then
        return false;
    fi;
    if IsBound( R.isFpPair ) then
        return R.isFpPair;
    fi;

    if (   ( not (      IsBound( R.fp )
                    and IsFpGroup( R.fp ) ) ) 
        or ( not (      IsBound( R.perm )
                    and IsPermGroup( R.perm ) ) )
        or ( not (      IsBound( R.f2p )
                    and IsIsomorphism( R.f2p )
                    and ( R.f2p.source = R.fp ) 
                    and ( R.f2p.range = R.perm ) ) ) ) then
        Print( "!! Error: invalid  FpPair !! \n" );         
        return false;
    fi; 

    R.isFpPair := true;
    return true;
    end;

###############################################################################
##
#F  IsSemidirectPair        test whether a semidirect product / perm group pair
##

IsSemidirectPair := function( R )

    if not IsRec( R ) then
        return false;
    fi;
    if IsBound( R.isSemidirectPair ) then
        return R.isSemidirectPair;
    fi;
    if (   ( not (      IsBound( R.sdp )
                    and IsSemidirectProduct( R.sdp ) ) ) 
        or ( not (      IsBound( R.perm )
                    and IsPermGroup( R.perm ) ) )
        or ( not (      IsBound( R.p2s )
                    # and IsIsomorphism( R.p2s )
                    and ( R.p2s.source = R.perm ) 
                    and ( R.p2s.range = R.sdp ) ) )
        or ( not (      IsBound( R.s2p )
                    # and IsIsomorphism( R.s2p )
                    and ( R.s2p.source = R.sdp ) 
                    and ( R.s2p.range = R.perm ) ) ) ) then
        return false;
    fi;
    R.isSemidirectPair := true;
    return true;
end;

###############################################################################
##
#F  PairIsomorphism                 uses an isomorphism to set up a Pair record
##
##  Problem (5/7/96):  IsHomomorphism(iso)  takes a long time
##                     so isomorphism testing has been removed.

PairIsomorphism := function( iso )

    local  ok, G, P, R, imag, gens, inv;

    ### isomorphism testing disabled: ###
    # if not IsIsomorphism( iso ) then
    #     Error( "input parameter must be an isomorphism" );
    # fi;

    G := iso.source;
    if IsFpGroup( G ) then
        if IsBound( G.fpPair ) then
            return G.fpPair;
        fi;
    elif IsSemidirectProduct( G ) then
        if IsBound( G.semidirectPair ) then
            return G.semidirectPair;
        fi;
    elif IsAutomorphismGroup( G ) then
        if IsBound( G.autoPair ) then
            return G.autoPair;
        fi; 
    else
        Error( "Source must be an FpGroup, AutoGroup or SemidirectProduct" );
    fi;
    P := iso.range;
    if not IsPermGroup( P ) then
        Error( "Range group must be a permutation group" );
    fi;

    # set up the inverse isomorphism using generators and their images
    gens := G.generators;
    # imag := List( gens, x -> Image( iso, x ) );
    imag := iso.genimages;
    inv := GroupHomomorphismByImages( P, G, imag, gens );
    inv.generators := imag;
    inv.genimages := gens;
    R := rec( );
    R.perm := P;
    if IsFpGroup( G ) then
        R.fp := G;
        R.f2p := iso;
        R.p2f := inv;
        ok := IsFpPair( R );
    elif IsSemidirectProduct( G ) then
        R.sdp := G;
        R.s2p := iso;
        R.p2s := inv;
        ok := IsSemidirectPair( R );
    elif IsAutomorphismGroup( G ) then
        R.auto := G;
        R.a2p := iso;
        R.p2a := inv;
        ok := IsAutomorphismPair( R );
    else
        Error( "wrong type of source group" );
    fi;
    return R;
end;

###############################################################################
##
#F  RegularFpPair                         construct the Regular Representation
##                                                FpPair of an FpGroup  G

RegularFpPair := function( G )

    local id, reg, P, phi;

    if not IsFpGroup( G ) then
        Error( "Input parameter should be an FpGroup" );
    fi;

    if IsBound( G.regularPair ) then
        return G.regularPair;
    else
        id := Subgroup( G, [ ] );
        reg := OperationCosetsFpGroup( G, id );
        phi := OperationHomomorphism( G, reg );
      # phi := GroupHomomorphismByImages(G, reg, G.generators, reg.generators);
        IsHomomorphism( phi );
        if IsBound( G.name ) then
            reg.name := Concatenation( "Reg(", G.name, ")" );
        fi;
        P := PairIsomorphism( phi );
        P.isRegularPair := true;
        P.degree := Size( reg );
        P.position := 1;
        G.regularPair := P;
        return P;
    fi;
end;

###############################################################################
##
#F  MinTransitiveFpPair            find a minimal core-free subgroup of K
#F                             and hence a faithful permutation representation
#F                                             isoKRK --> RK        

MinTransitiveFpPair := function( K )

    local  genK, ind, posn, regK, pairK, ngK, genregK, lat, llat, oK, i, j,
           S, cosS, elcosS, J, oJ, kj, jj, genJ, isoKJ, res, rep;

    if not IsFpGroup( K ) then
        Error( "parameter should be an FpGroup" );
    fi;
    genK := K.generators;
    pairK := RegularFpPair( K );
    regK := pairK.perm;
    genregK := pairK.p2f.generators;
    ngK := Length( genregK );
    if ( XModPrintLevel > 1 ) then
        Print( "\nGenerators of Regular Representation of K :-\n" );
        for i in [1..ngK] do
            Print( genregK[i], "\n" );
        od;
    fi;
    lat := Lattice( regK );
    oK := Size( regK );
    ind := oK;
    llat := Length( lat.classes );
    for i in Reversed( [1..llat] ) do
        S := lat.classes[i].representative;
        ind := Index( regK, S );
        posn := i;
        cosS := Cosets( regK, S );
        elcosS := 0*[1..ind];
        for j in [1..ind] do
            elcosS[j] := Elements( cosS[j] );
        od;
        genJ := 0*[1..ngK];
        for j in [1..ngK] do
            kj := genregK[j];
            jj := Permutation( kj, elcosS, OnRightCosets );
            genJ[j] := jj; 
        od;
        J := Group( genJ, () );
        oJ := Size( J );
        if ( oJ = oK ) then
            if ( XModPrintLevel > 1 ) then
                Print( "\nMinimal degree representation found at class ", i );
                Print( " with index ", ind, "\n" );
                Print( "\nRepresentative subgroup has generators:\n" );
                Print( J.generators, "\n" );
            fi;
            isoKJ := GroupHomomorphismByImages( K, J, genK, genJ );
            rep := PairIsomorphism( isoKJ );
            rep.isMinTransitivePair := true;
            rep.generators := genJ;
            rep.degree := ind;
            rep.position := posn;
            K.fpPair := rep;
            return rep;
        fi;
    od;
end;

###############################################################################
##
#F  MinBitransitiveFpPair            find a minimal degree representation for K
#F                                     using the union of two transitive K-sets

MinBitransitiveFpPair := function( K )

    local  H1, H2, lat, ll, ll1, i, i1, i2, i12, L, J, j1, j2,
           x, x1, x2, d1, d2, P, oP, degP, isoKP, genP, genK, degK, posn,
           corelist, icore, ind12, oK, k1, k2, gK, ngK, regK, pairK,
           S1, S2, cosS1, cosS2, elS12, rep;

    if not IsFpGroup( K ) then
        Error( "parameter should be an FpGroup" );
    fi;
    genK := K.generators;
    if IsBound( K.fpPair ) then
        pairK := K.fpPair;
    else
        pairK := RegularFpPair( K );
    fi;
    degK := pairK.degree;
    posn := pairK.position;
    regK := pairK.perm;
    k1 := regK.generators[1];
    k2 := regK.generators[2];
    # Print( "\nIn  MinBitransitiveFpPair  with posn, k1, k2 = " );
    # Print( posn, "  ", k1, "  ", k2, "\n" );
    oK := Size( regK );
    ngK := Length( genK );
    lat := Lattice( regK );
    ll := Length( lat.classes );
    ll1 := ll-1;
    if ( XModPrintLevel > 1 ) then
        Print( "\nSubgroup Lattice has ", ll, " classes.\n\n" );
    fi;
    H1 := lat.classes[ posn ];
    H2 := lat.classes[ ll ];
    degP := degK;
    corelist := 0 * [ 1..ll1 ];
       icore := 0 * [ 1..ll1 ]; 
    for i in [1..ll1] do
        L := lat.classes[i].representative;
        J := Core( regK, L );
        corelist[i] := J;
        icore[i] := Index( regK, J );
    od;
    j1 := 0;
    j2 := 0;
    for i1 in [2..ll1] do
        if (icore[i1] < oK) then
            for i2 in [(i1+1)..ll1] do
                if ( icore[i2] < oK ) then
                    J := Intersection( corelist[i1], corelist[i2] );
                    if ( Size(J) = 1 ) then
                        ind12 := icore[i1] + icore[i2];
                        if ( ind12 < degP ) then
                            degP := ind12;
                            j1 := i1;
                            j2 := i2;
                        fi;
                    fi;
                fi;
            od;
        fi;
    od;

    if ( degP < degK ) then
        Print( "\nLower degree representation found using pairs: " );
        Print( " degP := ", degP, "\n" );
        S1 := lat.classes[j1].representative;
        S2 := lat.classes[j2].representative;
        cosS1 := Cosets( regK, S1 );
        cosS2 := Cosets( regK, S2 );
        i1 := icore[j1];
        i2 := icore[j2];
        i12 := i1+i2;
        elS12 := 0*[1..i12];
        for i in [1..i1] do
            elS12[i] := Elements( cosS1[i] );
        od;
        for i in [1..i2] do
            elS12[i1+i] := Elements( cosS2[i] );
        od;
        x1 := Permutation( k2, elS12, OnRightCosets );
        x2 := Permutation( k1, elS12, OnRightCosets );
        P := Group( x1, x2 );
        genP := P.generators;
        Print( "\nP = ", P, "\n" );
        isoKP := GroupHomomorphismByImages( K, P, genK, genP );
        oP := Size( P );
        if ( oP <> oK ) then
            Print( "\n\n !!! oK <> oP !!! \n\n" );
        fi;
        Print( "\nGenerators for the representation :- \n" );
        for gK in [1..ngK] do
            x := genK[gK];
            Print( "x,isoKP(x) = ", x, "  ", Image( isoKP, x ), "\n" );
        od;
        rep := PairIsomorphism( isoKP );
        rep.isMinBitransitivePair := true;
        rep.degree := degP;
        K.fpPair := rep;
    else
        rep := pairK;
    fi;
    return rep;
end;

###############################################################################
##
#F  MinTransitiveCyclicFpPair    find a minimal core-free cyclic subgroup of K
#F                               & hence a faithful permutation representation
#F                                               isoKRK --> RK

MinTransitiveCyclicFpPair := function( K )

    local  genK, posn, regK, pairK, ngK, genregK, ccK, nccK, oK,
           i, j, Id, S, cosS, elcosS, g, og, core, deg, mindeg,
           J, oJ, kj, jj, genJ, isoKJ, rep;

    Print( "\nSeeking transitive representation with a cyclic subgroup.\n" );
    genK := K.generators;
    pairK := RegularFpPair( K );
    regK := pairK.perm;
    genregK := regK.generators;
    ngK := Length(genregK);
    ccK := ConjugacyClasses( regK );
    oK := Size( regK );
    mindeg := oK;
    posn := 1;
    Id := Subgroup( regK, [ ] );
    nccK := Length( ccK );
    if ( XModPrintLevel > 1 ) then
        Print( "Trivial core found at classes :- \n" );
    fi;
    for i in Reversed( [2..nccK] ) do
        g := ccK[i].representative;
        S := Subgroup( regK, [g] );
        core := Core( regK, S );  
        if ( core = Id ) then
            og := Order( regK, g );
            deg := oK/og;
            if ( XModPrintLevel > 1 ) then
                Print( i, " " );
            fi;
            if ( deg < mindeg ) then
                mindeg := deg;
                posn := i;
            fi;
        fi;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "\nmindeg = ", mindeg, "  at class ", posn, "\n" );
    fi;
    if ( posn > 1 ) then
        g := ccK[posn].representative;
        og := Order( regK, g );
        deg := oK/og;
        S := Subgroup( regK, [g] );
        cosS := Cosets( regK, S );
        elcosS := 0 * [1..deg];
        for j in [1..deg] do
            elcosS[j] := Elements( cosS[j] );
        od;
        genJ := 0 * [1..ngK];
        for j in [1..ngK] do
            kj := genregK[j];
            jj := Permutation( kj, elcosS, OnRightCosets );
            genJ[j] := jj;
        od;
        Print( "\nList of generators found :- \n", genJ, "\n" );
        J := Group( genJ, () );
        oJ := Size( J );
        if ( oJ = oK ) then
            if ( XModPrintLevel > 1 ) then
                Print( "\nMinimal degree representation found at class ", i );
                Print( " with degree ",mindeg,"\n" );
                Print( "\nRepresentative subgroup has generators:\n" );
                Print( J.generators, "\n" );
            fi;
            isoKJ := GroupHomomorphismByImages( K, J, genK, genJ );
            rep := PairIsomorphism( isoKJ );
            rep.isMinTransitiveCyclicPair := true;
            rep.degree := mindeg;
        fi;
    else
        rep := pairK;
    fi;
    return rep;
end;

##############################################################################
##
#F  SmallFpPair    call one or more of the functions:
#F                 MinTransitiveRep, MinBitransitiveRep,MinTransitiveCyclicRep
#F                                   to find a perm rep of  K  of small degree

SmallFpPair := function( K )

    local   rep, rep1, RK1, genRK1, isoKRK1, degRK1, minpos,
            rep2, RK2, genRK2, isoKRK2, degRK2, pairs,
            genK, pairK, regK, oK,  RK,  genRK,  isoKRK,  degRK ;

    if IsBound( K.fpPair ) then
        return K.fpPair;
    fi;
    if not ( IsBound( K.isFpGroup ) and K.isFpGroup ) then
        Error( "SmallFpPair requires an FpGroup as input parameter" );
    fi;

    genK := K.generators;
    pairK := RegularFpPair( K );
    oK := Size( K );
    if (oK<64) then
        rep1 := MinTransitiveFpPair( K );
        RK1 := rep1.perm;
        genRK1 := rep1.generators;
        isoKRK1 := rep1.f2p;
        degRK1 := rep1.degree;
        minpos := rep1.position;
    else
        rep1 := MinTransitiveCyclicFpPair( K );
        RK1 := rep1.perm;
        isoKRK1 := rep1.f2p;
        degRK1 := rep1.degree;
    fi;
    pairs := ( oK < 31 );
    if ( IsAbelian( RK1 ) ) then
        Print( "Group is abelian \n");
        pairs := false;
    fi;
    if pairs then 
        rep2 := MinBitransitiveFpPair( K );
        if ( XModPrintLevel > 2 ) then
            Print( "In SmallFpPair, at pairs, where rep2 = \n", rep2, "\n" );
        fi;
        degRK2 := rep2.degree;
        if ( degRK2> 0 ) then
            RK2 := rep2.perm;
            isoKRK2 := rep2.f2p;
        fi;
    fi;
    if ( (pairs = false) or ( degRK2 = 0 ) or ( degRK2 > degRK1 ) ) then
        RK := RK1;
        degRK := degRK1;
        isoKRK := isoKRK1;
    else 
        RK := RK2;
        degRK := degRK2;
        isoKRK := isoKRK2;
    fi;
    rep := PairIsomorphism( isoKRK );
    rep.degree := degRK;
    K.fpPair := rep;
    return rep;
end;

##############################################################################
##
#F  FpPairPermGroup                 PermGroup --> pair [ FpGroup -> PermGroup]
##

FpPairPermGroup := function( P )

    local pres, phi, degP, G, R;

    if IsBound( P.fpPair ) then
        return P.fpPair;
    fi;
    if not IsPermGroup( P ) then
        Error( "Input parameter should be a Permutation Group" );
    fi;

    degP := Length( PermGroupOps.MovedPoints( P ) );
    pres := PresentationViaCosetTable( P );
    G := FpGroupPresentation( pres );

    phi := GroupHomomorphismByImages( G, P, G.generators, P.generators );
    if not IsHomomorphism( phi ) then
        Error( "phi fails to be a homomorphism within FpPairPermGroup" );
    fi;

    R := PairIsomorphism( phi );
    R.degree := degP;
    R.presentation := pres;
    if IsBound( P.name ) then
        G.name := Concatenation( P.name, "Fp" );
        R.name := Concatenation( P.name, "Pair" );
    fi;
    P.fpPair := R;
    return R;
end;

##############################################################################
##
#F  FpPair               PermGroup or FpGroup --> pair [ FpGroup -> PermGroup]
##

FpPair := function( J )

    local R;

    if IsBound( J.fpPair ) then
        return J.fpPair;
    fi;

    if ( not IsPermGroup( J ) and not IsFpGroup( J ) ) then
        Error( "Input parameter should be a PermGroup or an FpGroup" );
    fi;

    if IsPermGroup( J ) then
        R := FpPairPermGroup( J );
    else
        R := SmallFpPair( J );
    fi;
    return R;
end;

##############################################################################
##
#F  SemidirectPair                 set up SemidirectProduct <-> PermGroup pair
##
#A  Altered version by Derek Holt, 6/9/96, attempts to achieve smaller degree

SemidirectPair := function( G )

    local  ok, g, oG, genG, img, P, oP, genP, L1, L2, L3, perm, phiGP, PI,
           num1, num2, new, R, r, S, s, ok, map, ordS, genno, LS, LR;

    if not IsSemidirectProduct( G ) then
        Error( "Input parameter must be a semidirect product group." );
    fi;
    if IsBound( G.semidirectPair ) then
        return G.semidirectPair;
    fi;

    genG := G.generators;
    oG := Size( G );
    # if R>=S then produce rep of degree Size(R)
    R := G.groups[1];
    S := G.groups[2];

    # Beginning of modified code by dfh
    # Take the action of R on S  by conjugation.
    LS := Elements( S );
    map := G.mapping;
    genP := [ ];
    for r in R.generators do
        L2 := List( LS, x -> Image(Image(map, r), x ));
        L3 := List( L2, x -> Position( LS, x ) );
        perm := PermList( L3 );
        Add( genP, perm );
    od;
    # Check order of group at this stage, to see if the action is faithful.
    # If not, adjoin the action of R as a permutation group (if it is one)
    # or by right multiplication on itself if not.
    P := Group( genP, () );
    if Size(P) <> Size(R) then
        ordS := Size(S);
        LR := Elements( R );
        genno := 0;
        for r in R.generators do
            if IsPermGroup(R) then
                L2 := ListPerm( r );
                L3 := List( L2, x -> x + ordS );
            else
                L2 := List( LR, x -> x*r );
                L3 := List( L2, x -> Position( LR, x ) + ordS );
            fi;
            L3 := Concatenation( [1..ordS], L3 );
            perm := PermList( L3 );
            genno := genno + 1;
            genP[genno] := genP[genno] * perm;
        od;
    fi;
    # Take the action of S on S by right multiplication.
    for s in S.generators do
        L2 := List( LS, x -> x*s );
        L3 := List( L2, x -> Position( LS, x ) );
        perm := PermList( L3 );
        Add( genP, perm );
    od;
    P := Group( genP, () );
    phiGP := GroupHomomorphismByImages( G, P, genG, genP );
    # phiGP.genimages := List( G.generators, g -> Image( phiGP, g ) );
    # Why not following, which is much quicker?
    phiGP.genimages := genP;
    # End of modified code by dfh

    PI := PairIsomorphism( phiGP );
    P.semidirect := G;
    if IsBound( G.name ) then
        P.name := Concatenation( "Perm(", G.name, ")" );
    fi;
    G.semidirectPair := PI;
    return( PI );
end;

##############################################################################
##
#F  OldSemidirectPair              set up SemidirectProduct <-> PermGroup pair
##

OldSemidirectPair := function( G )

    local  ok, g, oG, genG, img, P, oP, genP, L1, L2, L3, perm,
           phiGP, PI, num1, num2, new, R, r, S, s, ok;

    if not IsSemidirectProduct( G ) then
        Error( "Input parameter must be a semidirect product group." );
    fi;
    if IsBound( G.semidirectPair ) then
        return G.semidirectPair;
    fi;

    genG := G.generators;
    oG := Size( G );
    # if R>=S then produce rep of degree Size(R)
    R := G.groups[1];
    S := G.groups[2];
    if IsSubgroup( R, S ) then
        L1 := Elements( R );
        genP := [ ];
        for r in R.generators do
            L2 := List( L1, x -> x^r );
            L3 := List( L2, x -> Position( L1, x ) );
            perm := PermList( L3 );
            Add( genP, perm );
        od;
        for s in S.generators do
            L2 := List( L1, x -> x*s );
            L3 := List( L2, x -> Position( L1, x ) );
            perm := PermList( L3 );
            Add( genP, perm );
        od;
        P := Group( genP, () );
        phiGP := GroupHomomorphismByImages( G, P, genG, genP );
    else
        # determine the closure of genG under conjugation
        L2 := genG;
        num1 := 0;
        num2 := Length( genG );
        new := genG;
        while ( num1 <> num2 ) do
            num1 := num2;
            L1 := L2;
            for g in genG do
                img := List( new, x -> x^g );
                L2 := Union( L2, img );
            od;
            new := Difference( L2, L1 );
            num2 := Length( L2 );
        od;
        # construct perm group P using conjugation
        P := Operation( G, L2 );
        genP := P.generators;
        phiGP := OperationHomomorphism( G, P );
    fi;
    oP := Size( P );
    ok := true;
    # too expensive to use:
    # ok := IsHomomorphism( phiGP );
    if not ( ok and ( oG = oP ) )then
        # construct the regular representation
        L2 := Elements( G );
        P := Operation( G, L2, OnRight );
        genP := P.generators;
        phiGP := OperationHomomorphism( G, P );
    fi;
    phiGP.genimages := List( G.generators, g -> Image( phiGP, g ) );
    PI := PairIsomorphism( phiGP );
    P.semidirect := G;
    if IsBound( G.name ) then
        P.name := Concatenation( "Perm(", G.name, ")" );
    fi;
    G.semidirectPair := PI;
    return( PI );
end;

##############################################################################
##
#F  PrintList        
##

PrintList := function( L )

    local  ok, len, i;

    ok := IsList( L );
    if ok then
        len := Length( L );
        for i in [1..len] do
            Print( L[i], "\n" );
        od;
    else
        Print( "Function  PrintList  can only print a list.\n" );
    fi;
end;

##############################################################################
##
#F  DistinctRepresentatives                 Hall's algorithm on a list of sets
##

DistinctRepresentatives := function( L )

    local  n, rep, U, len, i, j, k, used, found, S, T, M, P, x, y, z;

    if not  ( IsList( L ) and
              ( ForAll( L, IsList ) or ForAll( L, IsSet ) ) ) then
        Error( "argument should be a list of sets" );
    fi;
    n := Length( L );
    U := [1..n];
    len := 0 * U;
    for i in U do
        S := L[i];
        if IsList( S ) then
            S := Set( S );
        fi;
        len[i] := Length( S );
        if ( len[i] = 0 ) then
            Error( "subsets must be non-empty" );
        fi;
        if not ForAll( S, j -> ( j in U ) ) then
            Error( "each set must be a subset of [1..n]" );
        fi;
    od;
    rep := 0 * U;
    used := 0 * U;
    rep[1] := L[1][1];
    used[ rep[1] ] := 1;
    for i in [2..n] do
        found := false;
        S := L[i];
        j := 0;
        while ( ( j < len[i] ) and not found ) do 
            j := j+1;
            x := S[j];
            if ( used[x] = 0 ) then
                rep[i] := x;
                used[x] := i;
                found := true;
            fi;
        od;
        # construct the graph component
        T := Copy( S );
        M := List( T );
        P := 0 * U;
        for x in M do
            P[x] := i;
        od;
        j := 0;
        while not found do
            j := j+1;
            x := M[j];
            k := used[x];
            if ( k = 0 ) then
                # reassign representatives
                y := P[x];
                while ( y <> i ) do
                    z := rep[y];
                    rep[y] := x;
                    used[x] := y;
                    x := z;
                    y := P[x];
                od;
                rep[i] := x;
                used[x] := i;
                found := true;
            else
                for y in L[k] do
                    if not ( y in T ) then
                        Add( M, y );
                        P[y] := k;
                        T := Union( T, [y] );
                    fi;
                od;
            fi;
            if ( ( not found ) and ( j = Length( M ) ) ) then
                Print( "Hall condition not satisfied!\n" );
                return false;
            fi;
        od;
    od;
    return rep;
end;

##############################################################################
##
#F  CommonRepresentatives                    for two lists of sets/lists  J, K
##                  returns [ rep, L ] where rep is a list of representatives,
##                            and L determines the appropriate reordering of K

CommonRepresentatives := function( J, K )

    local  U, i, j, k, m, n, lenJ, lenK, S, L, I, rep, perm, common;

    if not  ( IsList( J ) and 
              ( ForAll( J, IsList ) or ForAll( J, IsSet ) ) ) then
        Error( "first argument should be a list of sets" );
    fi;
    m := Length( J );
    if not  ( IsList( K ) and
              ( ForAll( K, IsList ) or ForAll( K, IsSet ) ) ) then
        Error( "second argument should be a list of sets" );
    fi;
    n := Length( K );
    if not ( m = n ) then
        Error( "lists <J> and <K> have unequal length" );
    fi;
    U := [1..n];
    lenJ := 0 * U;
    lenK := 0 * U;
    for i in U do
        S := J[i];
        if IsList( S ) then
            S := Set( S );
        fi;
        lenJ[i] := Length( S );
        if ( lenJ[i] = 0 ) then
            Error( "sets must be non-empty" );
        fi;
        S := K[i];
        if IsList( S ) then
            S := Set( S );
        fi;
        lenK[i] := Length( S );
        if ( lenK[i] = 0 ) then
            Error( "sets must be non-empty" );
        fi;
    od;
    L := List( U, x -> [ ] );
    for i in U do
        S := J[i];
        for j in U do
            I := Intersection( S, K[j] );
            if ( Length( I ) > 0 ) then
                Add( L[i], j );
            fi;
        od;
    od;
    rep := DistinctRepresentatives( L );
    perm := PermList( rep );
    K := Permuted( K, perm^-1 );
    common := 0 * U;
    for i in U do
        I := Intersection( J[i], K[i] );
        common[i] := I[1];
    od;
    return [ common, rep ];
end;

##############################################################################
##
#F  IsCommonTransversal         test if a common left/right cosets transversal
##

IsCommonTransversal := function( arg )

    local  nargs, usage, G, H, T,
           eG, eH, T, oG, oH, g, h, t, pos, ind, found;

    usage := "\nUsage:  IsCommonTransversal( G, H, T )\n";
    nargs := Length( arg );
    if ( nargs <> 3 ) then
        Print( usage );
        return false;
    fi;
    G := arg[1];
    H := arg[2];
    if not ( IsPermGroup( G ) and IsPermGroup( H ) 
                              and IsSubgroup( G, H ) ) then
        Print( "second perm group parameter must be subgroup of first" );
        return false;
    fi;
    T := arg[3];
    if not IsList( T ) then
        Print( "third parameter must be a list of representatives" );
        return false;
    fi;
    oG := Size( G );
    oH := Size( H );
    eG := Elements( G );
    eH := Elements( H );
    ind := oG/oH;
    found := 0 * [1..oG];
    for t in T do
        if not ( t in eG ) then
            Print( "element of T not in G" );
            return false;
        fi;
        for h in eH do
            g := t*h;
            pos := Position( eG, g );
            found[pos] := found[pos] +1;
            g := h*t;
            pos := Position( eG, g );
            found[pos] := found[pos] + 1;
        od;
    od;  
    for t in [1..oG] do
        if not ( found[t] = 2 ) then
            Print( eG[t], " found ", found[t], " times\n" );
            return false;
        fi;
    od;
    return true;
end;

##############################################################################
##
#F  CommonTransversal               find a common left/right coset transversal
##                                    using system of common representatives

CommonTransversal := function( G, H )

    local  L, EL, R, ER, T;

    if not IsSubgroup( G, H ) then
        Error( "<H> must be a subgroup of <G>" );
    fi;
    L := LeftCosets( G, H );
    R := RightCosets( G, H );
    EL := List( L, x -> Elements( x ) );
    ER := List( R, x -> Elements( x ) );
    T := CommonRepresentatives( EL, ER );
    return T[1];
end;

##  end of file  util.g
