###############################################################################
##
##  deriv.g                     for GAP 3.4                    version 10/ 1/97
##
###############################################################################
##
#A  deriv.g                     GAP library                       Chris Wensley
#A                                                                    Murat Alp
#Y  Copyright
##
##  This file contains functions compute derivation of crossed modules
##
#H  $Log: deriv.g,v $
#H  Revision 1.1  1997/03/27 13:35:46  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

##############################################################################
##
#F  Cat1SectionByImagesOps
##

Cat1SectionByImagesOps :=
   OperationsRecord( "Cat1SectionByImagesOps", GroupHomomorphismByImagesOps );

Cat1SectionByImagesOps.Print := function( xi )
    Print( "Cat1SectionByImages( ", xi.source, ", ",
            xi.range, ", ", xi.generators, ", ", xi.genimages, " )" );
end;

Cat1SectionByImagesOps.IsRegular := function( xi )

    local  C, R, genR, genrng, h, im, him, rho, ok;

    C := xi.cat1;
    R := C.range;
    genR := xi.generators;
    genrng := [ 1..Length( genR ) ];
    h := C.head;
    im := xi.genimages;
    him := List( im, g -> Image( h, g ) );
    rho := GroupHomomorphismByImages( R, R, genR, him );
    ok := IsIsomorphism( rho );
    return ok;
end;

##############################################################################
##
#F  Cat1SectionsOps
##

Cat1SectionsOps := OperationsRecord( "Cat1SectionsOps", DomainOps );

Cat1SectionsOps.Print := function( sect )

    local  C, L, size, reg, isreg, isall, str;

    if not AreSections( sect ) then
        Print( "Not a sections record!\n" );
        return false;
    fi;
    C := sect.cat1;
    L := sect.genimageList;
    isreg := ( IsBound( sect.isReg ) and sect.isReg );
    isall := ( IsBound( sect.isAll ) and sect.isAll );
    size := Length( L );
    if IsBound( sect.regular ) then
        reg := sect.regular;
    else
        reg := 0;
    fi;
    if isreg then
        str := "Regular";
    elif isall then
        str := "All";
    else
        str := "";
    fi;
    Print( str, "Sections record for cat1-group " );
    Print( Cat1Name( C ), ",\n" );
    if isreg then
        Print( ": ", reg, " regular sections, others not found." );
    elif ( reg > 0 ) then
        Print( ": ", reg, " regular sections,  " );
        Print( size - reg, " irregular ones found." );
    else
        Print( ": ", size, " sections found but unsorted." );
    fi;
end;

##############################################################################
##
#F  XModDerivationByImagesOps
##

XModDerivationByImagesOps :=
    OperationsRecord( "XModDerivationByImagesOps", MappingOps );

XModDerivationByImagesOps.Print := function( chi )
    Print( "XModDerivationByImages( ", chi.source, ", ",
            chi.range, ", ", chi.generators, ", ", chi.genimages, " )" );
end;

XModDerivationByImagesOps.IsRegular := function( chi )

    local  X, S, bdy, R, genR, genrng, im, imrho, rho, ok;

    X := chi.xmod;
    S := X.source;
    bdy := X.boundary;
    R := X.range;
    genR := R.generators;
    im := chi.genimages;
    genrng := [ 1..Length( genR ) ];
    imrho := List( genrng, i -> genR[i] * Image( bdy, im[i] ) );
    rho := GroupHomomorphismByImages( R, R , genR, imrho);
    ok := IsIsomorphism( rho );
    return ok;
end;

##############################################################################
##
#F  XModDerivationsOps
##

XModDerivationsOps := OperationsRecord( "XModDerivationsOps", DomainOps );

XModDerivationsOps.Print := function( der )

    local  X, L, size, reg, isreg, isall, str;

    if not AreDerivations( der ) then
        Print( "Not a derivations record!\n" );
        return false;
    fi;
    X := der.xmod;
    L := der.genimageList;
    isreg := ( IsBound( der.isReg ) and der.isReg );
    isall := ( IsBound( der.isAll ) and der.isAll );
    size := Length( L );
    if IsBound( der.regular ) then
        reg := der.regular;
    else
        reg := 0;
    fi;
    if isreg then
        str := "Regular";
    elif isall then
        str := "All";
    else
        str := "";
    fi;
    Print( str, "Derivations record for crossed module " );
    Print( XModName( X ), ",\n" );
    if isreg then
        Print( ": ", reg, " regular derivations, others not found." );
    elif ( reg > 0 ) then
        Print( ": ", reg, " regular derivations,  " );
        Print( size - reg, " irregular ones found." );
    else
        Print( ": ", size, " derivations found but unsorted." );
    fi;
end;

##############################################################################
##
#F  AreSections                          test to see whether a sections record
##

AreSections := function( S )

    local isect;

    if not IsRec( S ) then
        return false;
    fi;
    if not ( IsBound( S.generators ) and IsBound( S.cat1 ) ) then
        return false;
    fi;
    if not ( S.generators = S.cat1.range.generators ) then
        Print( "Warning! <C.range.generators> has changed\n" );
        if ( XModPrintLevel > 1 ) then
            Print( "           S.generators = ", S.generators, "\n" );
            Print( "S.cat1.range.generators = ",
                    S.cat1.range.generators, "\n" );
        fi;
        return false;
    fi;
    if IsBound( S.areSections ) then
        return S.areSections;
    fi;
    isect := ( IsBound( S.genimageList ) and IsBound( S.cat1 )
              and ( IsBound( S.isAll ) or IsBound( S.isReg ) ) );
    S.areSections := isect;
    return isect;
end;

##############################################################################
##
#F  AreDerivations                   tests to see whether a derivations record
##

AreDerivations := function( D )

    local  isder;

    if not IsRec( D ) then
        return false;
    fi;
    if not ( IsBound( D.generators ) and IsBound( D.xmod ) ) then
        return false;
    fi;
    if not ( D.generators = D.xmod.range.generators ) then
        Print( "Warning! <D.range.generators> has changed\n" );
        if ( XModPrintLevel > 1 ) then
            Print( "           D.generators = ", D.generators, "\n" );
            Print( "D.xmod.range.generators = ",
                    D.xmod.range.generators, "\n" );
        fi;
        return false;
    fi;
    if IsBound( D.areDerivations ) then
        return D.areDerivations;
    fi;
    isder := ( IsBound( D.genimageList ) and IsBound( D.xmod )
              and ( IsBound( D.isAll ) or IsBound( D.isReg ) ) );
    D.areDerivations := isder;
    return isder;
end;

##############################################################################
##
#F  IsSection                        tests the section axioms for a cat1-group
##

IsSection := function( arg )

    local  nargs, usage, xi, C, im, Crng, h, t, xih, xit, idrng, ok, reg;

    usage := "\nUsage: IsSection( xi ) or IsSection( C, im );\n";
    nargs := Length( arg );
    if not IsRec( arg[1] ) then
        return false;
    fi;
    if ( nargs = 1 ) then
        xi := arg[1];
        if not ( IsBound( xi.genimages ) and IsBound( xi.cat1 ) ) then
            return false;
        fi;
        C := xi.cat1;
        im := xi.genimages;
        if ( IsBound( C.sections ) and
             ( C.sections.generators <> C.range.generators ) ) then
            Print( "<C.range.generators> has changed\n" );
            return false;
        fi;
        if IsBound( xi.isSection ) then
            return xi.isSection;
        fi;
    elif ( nargs = 2 ) then
        C := arg[1];
        if not ( IsCat1( C ) and IsList( arg[2] ) ) then
            Print( usage );
            return false;
        fi;
        im := arg[2];
        xi := Cat1SectionByImages( C, im );
    else
        Print( usage );
        return false;
    fi;
    if not ( IsHomomorphism( xi ) and
             ( xi.source = C.range ) and ( xi.range = C.source ) ) then
        Print( "<xi> not a homomorphism: C.range -> C.source\n" );
        return false;
    fi;
    Crng := C.range;
    idrng := InclusionMorphism( Crng, Crng );   
    h := C.head;
    t := C.tail;
    xit := xi * t;
    ok := ( xit = idrng );
    xih := xi * h;
    reg := IsMonomorphism( xih );
    xi.isSection := ok;
    xi.isRegular := reg;
    return ok;
end;

##############################################################################
##
#F  Cat1SectionByImages                               sets up GroupHomByImages
##

Cat1SectionByImages := function( arg )

    local  xi, R, G, genR, ngR, nargs, usage, C, im, isect;

    usage := "\nUsage: Cat1SectionByImages( C, im, [, true|false ] );\n";
    nargs := Length( arg );
    if ( (nargs < 2 ) or ( nargs > 4 ) ) then
        Print( usage );
        return false;
    fi;
    isect := ( ( nargs = 3 ) and IsBool( arg[3] ) and arg[3] );
    C := arg[1];
    im := arg[2];
    if not IsCat1( C ) then
        Error( "<C> must be a cat1-group" );
    fi;
    G := C.source;
    R := C.range;
    genR := R.generators;
    ngR := Length( genR );
    if not ( IsList( im ) and ( Length( im ) = ngR )
                   and ForAll( im, x -> ( x in G ) ) ) then
        Error( "<im> must be a list of |genR| elements in <G>" );
    fi;
    xi := GroupHomomorphismByImages( R, G, genR, im );
    if isect then
        xi.isSection := true;
    fi;
    xi.cat1 := C;
    xi.operations := Cat1SectionByImagesOps;
    return xi;
end;

##############################################################################
##
#F  IsDerivation                    checks bdy*chi is a hom & action preserved
##

IsDerivation := function( arg )

    local  nargs, usage, chi, im, inv, X, bdy, R, r, genR, invR, ord, S,
           genrng, rho, imrho, ok, g, i, j, r, s, t, u, v, w, aut, act,
           pair, fp, pres, T, numrels, lenrels, rels, rel, len, triv;

    if not IsRec( arg[1] ) then
        return false;
    fi;
    usage := " \nUsage: IsDerivation( chi ) or IsDerivation( X, im )\n" ;
    nargs := Length( arg );
    if ( nargs = 1 ) then
        chi := arg[1];
        if not ( IsBound( chi.genimages ) and IsBound( chi.xmod ) ) then
            return false;
        fi;
        X := chi.xmod;
        im := chi.genimages;
        if ( IsBound( X.derivations ) and
             ( X.derivations.generators <> X.range.generators ) ) then
            Print( "<X.range.generators> has changed\n" );
            return false;
        fi;
        if IsBound( chi.isDerivation ) then
            return chi.isDerivation;
        fi;
    elif ( nargs = 2 ) then
        X := arg[1];
        if not ( IsXMod( X ) and IsList( arg[2] ) ) then
            Print( usage );
            return false;
        fi;
        im := arg[2];
        chi := XModDerivationByImages( X, im );
    else
        Print( usage );
        return false;
    fi;
    S := X.source;
    R := X.range;
    genR := R.generators;
    triv := List( genR, r -> () );
    if ( chi.genimages = triv ) then
        return true;
    fi;
    invR := List( genR, r -> r^(-1) );
    genrng := [ 1..Length( genR ) ];
    bdy := X.boundary;
    act := X.action;

    # calculate  chi(r^-1)  for each generator  r
    inv := 0 * genrng;
    for j in genrng do
        r := genR[j];
        s := im[j];
        aut := Image( act, r^(-1) );
        inv[j] := Image( aut, s )^(-1);
    od;
    if ( XModPrintLevel > 3 ) then
        Print( "  Images = ", im, "inverses = ", inv, "\n" );
    fi;

    pair := FpPair( R );
    fp := pair.fp;
    pres := pair.presentation;
    T := pres.tietze;
    numrels := T[TZ_NUMRELS];
    rels := T[TZ_RELATORS];
    lenrels := T[TZ_LENGTHS];
    w := ();
    v := ();
    for i in [1..numrels] do
        rel := rels[i];
        len := lenrels[i];
        for j in Reversed( [1..len] ) do
            g := rel[j];
            if ( g > 0 ) then
                r := genR[g];
                s := im[g];
            else
                r := invR[-g];
                s := inv[-g];
            fi;
            aut := Image( act, v );
            u := Image( aut, s );
            w := u * w;
            v := r * v;
        od;
        if ( w <> () ) then
            if ( XModPrintLevel > 2 ) then
                Print( "chi(rel) <> ()  when rel = ", rel, "\n" );
            fi;
            return false;
        fi;
    od;
    chi.isDerivation := true;
    return true;
end;

##############################################################################
##
#F  XModDerivationByImages                                 sets up the mapping
##

XModDerivationByImages := function( arg )

    local  nargs, usage, X, im, isder, chi, R, S, genR, ngR;

    usage := "\nUsage: XModDerivationByImages( X, im, [, true|false ] );\n";
    nargs := Length( arg );
    if ( ( nargs < 2 ) or ( nargs > 4 ) ) then
        Print( usage );
        return false;
    fi;
    isder := ( ( nargs = 3 ) and IsBool( arg[3] ) and arg[3] );
    X := arg[1];
    im := arg[2];
    if not IsXMod( X ) then
        Error( "<X> must be a crossed module" );
    fi;
    S := X.source;
    R := X.range;
    genR := R.generators;
    ngR := Length( genR );
    if not ( IsList( im ) and ( Length( im ) = ngR )
                   and ForAll( im, x -> ( x in S ) ) ) then
        Error( "<im> must be a list of |genR| elements in S" );
    fi;

    chi := rec( );
    chi.isGeneralMapping := true;
    chi.domain := Mappings;
    chi.source := R;
    chi.range := S;
    chi.generators := genR;
    chi.genimages := im;
    chi.xmod := X;
    if isder then
        chi.isDerivation := true;
    fi;
    chi.operations := XModDerivationByImagesOps;
    return chi;
end;

###############################################################################
##
#F  DerivationSection    construct an XMod derivation from a cat1-group section
##

DerivationSection := function( xi )

    local  C, imxi, eR, eK, R, genR, ngR, S, X, imchi, r, er, s, i, chi;

    if not IsSection( xi ) then
        Error( "Parameter must be a section of a cat1-group" );
    fi;
    C := xi.cat1;
    imxi := xi.genimages;
    X := XModCat1( C );
    S := C.kernel;
    R := C.range;
    genR := xi.generators;
    ngR := Length( genR );
    eR := C.embedRange;
    eK := C.embedKernel;
    imchi := 0 * [ 1..ngR ];
    for i in [ 1..Length( genR ) ] do
        r := genR[i];
        er := Image( eR, r );
        s := er^(-1) * imxi[i];
        imchi[i] := PreImagesRepresentative( eK, s );
        if ( XModPrintLevel > 2 ) then
            Print( "In xi->chi :- ", [ i, r, er, s], "\n" );
        fi;
    od;
    chi := XModDerivationByImages( X, imchi, true );
    chi.isDerivation := IsDerivation( chi );
    return chi;
end;

##############################################################################
##
#F  SectionDerivation        the cat1-group section determined by a derivation
##

SectionDerivation := function( chi )

    local  R, genR, ngR, xi, imchi, imxi, i, r, er, eR, eK, eKchi, g, C;

    if not IsDerivation( chi ) then
        Error( "Parameter must be a derivation of a crossed module" );
    fi;
    X := chi.xmod;
    R := X.range;
    genR := R.generators;
    ngR := Length( genR );
    C := Cat1XMod( X );
    eR := C.embedRange;
    eK := C.embedKernel;
    imchi := chi.genimages;
    eKchi := List( imchi, s -> Image( eK, s ) );
    imxi := 0 * [ 1..ngR ];
    for i in [ 1..Length( genR ) ] do
        r := genR[i];
        er := Image( eR, r );
        g := er * eKchi[i];
        imxi[i] := g;
    od;
    xi := Cat1SectionByImages( C, imxi, true );
    return xi;
end;

##############################################################################
##
#F  CompositeDerivation                 Whitehead composite of two derivations
##

CompositeDerivation := function( chi, chj )

    local  X, imi, imj, R, genR, numrng, rng, r, k,
           s, si, sj, bdy, bsj, imcomp, comp;

    if not ( IsDerivation( chi ) and IsDerivation( chj ) ) then
        Error( "function requires two derivations as input" );
    fi;
    X := chi.xmod;
    if not ( X = chj.xmod ) then
        Error( "<chi>,<chj> must be derivations for the SAME xmod" );
    fi;
    R := X.range;
    genR := chi.generators;
    if not ( chj.generators = genR ) then
        Error( "<chi>,<chj> have different sets of generators" );
    fi;
    numrng := Length( genR );
    rng := [ 1..numrng ];
    imi := chi.genimages;
    imj := chj.genimages;
    bdy := X.boundary;
    imcomp := 0 * rng;
    for k in rng do
        r := genR[k];
        sj := imj[k];
        si := imi[k];
        bsj := Image( bdy, sj );
        s := DerivationImage( chi, bsj );
        imcomp[k] := si * sj * s;
    od;
    comp := XModDerivationByImages( X, imcomp, true );
    return comp;
end;

##############################################################################
##
#F  CompositeSection                       Whitehead composite of two sections
##

CompositeSection := function( xi, xj )

    local  C, R, genR, h, e, hxj, ehxj, xihxj, r, im1, im2, im3,
           numrng, k, imcomp, comp;

    if not ( IsSection( xi ) and IsSection( xj ) ) then
        Error( "function requires two sections as input" );
    fi;
    C := xi.cat1;
    if not ( C = xj.cat1 ) then
        Error( "<xi>,<xj> must be sections of the SAME cat1-group" );
    fi;
    R := C.range;
    genR := xi.generators;
    if not ( genR = xj.generators ) then
        Error( "<xi>,<xj> have different generating sets" );
    fi;
    h := C.head;
    e := C.embedRange;
    hxj := xj * h;
    ehxj := hxj * e;
    xihxj := hxj * xi;
    numrng := Length( genR );
    imcomp := 0 * [1..numrng];
    for k in [1..numrng] do
        r := genR[k];
        im1 := Image( xj, r );
        im2 := Image( ehxj, r )^(-1);
        im3 := Image( xihxj, r );
        imcomp[k] := im1 * im2 * im3;
    od;
    comp := Cat1SectionByImages( C, imcomp, true );
    return comp;
end;

##############################################################################
##
#F  GenerationOrder         elements of G generated as words in the generators
##

GenerationOrder := function( G )

    local  G, oG, eG, ord, ngG, genG, g, i, j, k, P, pos, n, x, y;

    if not IsPermGroup( G ) then
        Error( "<G> must be a permutation group" );
    fi;
    if IsBound( G.generationOrder ) then
        return G.generationOrder;
    fi;
    oG := Size( G );
    eG := Elements( G );
    genG := G.generators;
    ngG := Length( genG );
    ord := 0 * [1..oG];
    ord[1] := 1;
    P := 0 * [1..oG];
    P[1] := [ 1, 0 ];
    for n in [1..ngG] do
        g := genG[n];
        pos := Position( eG, g );
        P[pos] := [ 1, n ];
        ord[ n+1 ] := pos;
    od;
    n := ngG + 1;
    k := 2;
    while ( k < oG ) do
        i := ord[k];
        x := eG[i];
        for j in [1..ngG] do
            g := genG[j];
            y := x * g;
            pos := Position( eG, y );
            if ( P[pos] = 0 ) then
                P[pos] := [ i, j ];
                n := n+1;
                ord[n] := pos;
            fi;
        od;
        k := k+1;
    od;
    G.generationOrder := ord;
    G.generationPairs := P;
    return ord;
end;

##############################################################################
##
#F  CheckGenerationPairs             G.generationPairs, G.generationOrder ok ?
##

CheckGenerationPairs := function( G )

    local  eG, oG, genG, P, i, g, x;

    if not ( IsBound( G.elements ) and IsBound( G.generationOrder )
                                   and IsBound( G.generationPairs ) ) then
        return false;
    fi;
    eG := Elements( G );
    oG := G.size;
    genG := G.generators;
    P := G.generationPairs;
    if ( P[1] <> [1,0] ) then
        return false;
    fi;
    for i in [2..oG] do
        x := eG[ P[i][1] ];
        g := genG[ P[i][2] ];
        if ( x*g <> eG[i] ) then
            if ( XModPrintLevel > 1 ) then
                Print( x, " * ", g, " <> ", eG[i], "\n" );
            fi;
            return false;
        fi;
    od;
    return true;
end;

##############################################################################
##
#F  DerivationTable                  returns list of lists of DerivationImages
##

DerivationTable := function( arg )

    local  nargs, X, str, T, i, chi, L, D, size;

    nargs := Length( arg );
    if not ( IsXMod( arg[1] ) or AreDerivations( arg[1] ) ) then
        Error( "Parameter must be a crossed module or derivations record" );
    fi;
    if AreDerivations( arg[1] ) then
        D := arg[1];
        X := D.xmod;
    else
        X := arg[1];
    fi;
    if IsBound( X.derivations ) then
        D := X.derivations;
        if ( IsBound( D.isReg ) and D.isReg ) then
            str := "reg";
        elif ( IsBound( D.isAll ) and D.isAll ) then
            str := "all";
        fi;
    else
        str := "reg";
    fi;
    if ( ( nargs = 2 ) and IsString( arg[2] ) ) then
        str := arg[2];
        if ( ( str <> "reg" ) and ( str <> "all" ) ) then
            Error( "<str> must be \"reg\" or \"all\" " );
        fi;
    fi;
    if ( str = "reg" ) then
        D := RegularDerivations( X );
    else
        D := AllDerivations( X );
    fi;
    if IsBound( D.table ) then
        return D.table;
    fi;
    L := D.genimageList;
    if ( str = "reg" ) then
        size := D.regular;
    else
        size := Length( L );
    fi;
    T := 0 * [1..size];
    for i in [1..size] do
        chi := XModDerivationByImages( X, L[i], true );
        T[i] := DerivationImages( chi );
    od;
    D.table := T;
    return T;
end;

##############################################################################
##
#F  DerivationImages                  returns list of positions of chi(r) in S
##

DerivationImages := function( arg )

    local  nargs, chi, X, R, S, elR, genR, ngR, ord, P, oR, elS,
           i, j, k, x, y, L, ok, imchi, act, agen, a;

    nargs := Length( arg );
    if ( nargs = 1 ) then
        chi := arg[1];
    elif ( ( nargs = 2 ) and IsXMod( arg[1] ) and IsList( arg[2] ) ) then
        chi := XModDerivationByImages( arg[1], arg[2] );
    else
        return false;
    fi;
    if not IsDerivation( chi ) then
        Error( "<chi> must be a derivation of a crossed module" );
    fi;
    if IsBound( chi.imagePositions ) then
        return chi.imagePositions;
    fi;
    if ( nargs = 1 ) then
        X := chi.xmod;
    fi;
    R := chi.source;
    genR := R.generators;
    ngR := Length( genR );
    if ( genR <> chi.generators ) then
        Error( "genR <> chi.generators" );
    fi;
    imchi := chi.genimages;
    act := X.action;
    agen := List( genR, r -> Image( act, r ) );
    S := chi.range;
    elR := Elements( R );
    elS := Elements( S );
    ord := GenerationOrder( R );
    P := R.generationPairs;
    if ( XModPrintLevel > 1 ) then
        ok := CheckGenerationPairs( R );
        Print( "checking generation pairs: ", ok, "\n" );
    fi;
    oR := Size( R );
    L := 0 * [1..oR];
    L[1] := Position( elS, () );
    for k in [2..oR] do
        i := ord[k];
        j := P[i][2];
        x := elS[ L[ P[i][1] ] ];
        a := agen[j];
        y := Image( a, x ) * imchi[j]; 
        L[i] := Position( elS, y );
    od;
    chi.imagePositions := L;
    return L;
end;

##############################################################################
##
#F  DerivationImage                  image of  r \in R  by the derivation  chi
##

DerivationImage := function( chi, r )

    local  X, S, R, genR, imchi, elR, elS, ngR, genrng, j, g, s, u, v,
           rpos, spos, ord, P, a, act;

    if not IsDerivation( chi ) then
        Error( "first parameter must be a crossed module derivation" );
    fi;
    if ( r = () ) then
        return ();
    fi;
    X := chi.xmod;
    S := X.source;
    R := X.range;
    elR := Elements( R );
    elS := Elements( S );
    genR := chi.generators;
    imchi := chi.genimages;
    if not ( r in R ) then
        Error( "second parameter must be an element of chi.source" );
    fi;
    rpos := Position( elR, r );
    if IsBound( chi.imagePositions ) then
        spos := chi.imagePositions[ rpos ];
        return elS[ spos ];
    fi;
    ord := GenerationOrder( R );
    P := R.generationPairs;
    j := P[rpos][1];
    g := P[rpos][2];
    if ( j = 1 ) then   # r is a generator
        return imchi[g];
    fi;
    act := X.action;
    u := imchi[g];
    v := genR[g];
    while ( j > 1 ) do
        g := P[j][2];
        j := P[j][1];
        s := imchi[g];
        a := Image( act, v );
        u := Image( a, s ) * u;
        v := genR[g] * v;
    od;
    return u;
end;

###############################################################################
##
#F  BacktrackSectionsJ          recursion used by RegularSections & AllSections
##

BacktrackSectionsJ := function( sectdata, k, im, j )

    local  partim, ok, num, g, phi;
    
    if k = 0 then
        im := ShallowCopy( im );
        partim := [ im ];
    else
        partim := [ ];
        for g in sectdata.preimages[ j ] do
            im[j] := sectdata.embedgen[j] * g;
            phi := GroupHomomorphismByImages
                  ( sectdata.subgps[j], sectdata.range, 
                    sectdata.subgen[j], Sublist( im, [1..j] ) ); 
            ok := IsHomomorphism( phi );
            if ok then
               Append( partim, BacktrackSectionsJ( sectdata, k-1, im, j+1 ) );
            fi;
        od;
    fi;    
    return partim;
end;

##############################################################################
##
#F  SectionsByEndoClasses    finds subgroups of G isomorphic to quotients of G
##

SectionsByEndoClasses := function( C, str, exists )

    local  R, genR, ngR, H, genH, endR, rho, rhoR, dispR, G, bdy, classes,
           imbdy, elimbdy, sect, images, sectdata, genim, reps, ok, ok1,
           count, d, xi, LS, LG, sectrec, regnum, k, iso, proj, im, aut,
           numaut, autos, a, rho, psi, numconj, conj, g, gim, l, b,
           auto, endo, c0, c, class;

    G := C.source;
    R := C.range;
    genR := R.generators;
    ngR := Length( genR );
    bdy := C.boundary;
    if IsSurjective( bdy ) then
        imbdy := R;
    else
        imbdy := Subgroup( R, bdy.genimages );
    fi;
    elimbdy := Elements( imbdy );
    if ( Length( elimbdy ) = 1 ) then
        C.isRModule := true;
    fi;
    if exists then
        c0 := 2;
    else
        c0 := 1;
    fi;
    if ( c0 = 2 ) then    # use existing regular images
        images := C.sections.genimageList;
        regnum := Length( images );
    else
        images := [ ];
    fi;
    sectdata := rec( );
    sectdata.range := G;
    sectdata.subgen := List( [1..ngR], i -> Sublist( genR, [1..i] ) ); 
    sectdata.subgps :=
             List( [1..ngR], i -> Subgroup( R, sectdata.subgen[i] ) );
    sectdata.embedgen := List( genR, r -> Image( C.embedRange, r ) );

    if ( str = "reg" ) then
        auto := rec( );
        auto.quotient := R;
        auto.projection := InclusionMorphism( R, R );
        auto.autoGroup := AutomorphismGroup( R );
        auto.isomorphism := auto.projection;
        auto.conj := [ () ];
        classes := [ auto ];
    else
        endo := EndomorphismClasses( R, 1 );
        classes := endo.classes;
        reps := endo.latticeReps;
    fi;
    count := 0;
    for c in [ c0..Length( classes ) ] do
        class := classes[c];
        if ( str = "reg" ) then
            H := R;
        else
            H := reps[ class.rangeNumber ];
        fi;
        conj := class.conj;
        proj := class.projection;
        iso := class.isomorphism;
        aut := class.autoGroup;
        autos := Elements( aut );
        numaut := Length( autos );
        numconj := Length( conj );
        endR := 0 * [ 1..(numaut*numconj) ];
        for a in [1..numaut] do
            rho := autos[a];
            psi := proj * rho * iso;
            im := List( genR, x -> Image( psi, x ) );
            endR[a] := im;
        od;
        for l in [ 2..numconj ] do
            g := conj[l];
            b := (l-1) * numaut;
            for a in [ 1..numaut ] do
                im := endR[a];
                gim := List( im, x -> x^g );
                endR[ b+a ] := gim;
            od;
        od;
        for rhoR in endR do
            rho := GroupHomomorphismByImages( R, R, genR, rhoR );
            dispR := List( genR, r -> r^(-1) * Image( rho, r ) );
            ok := true;
            for d in dispR do
                if ( Position( elimbdy, d ) = false ) then
                    ok := false;
                fi;
            od;
            if ok then
                LS := List( dispR, d -> PreImages( bdy, d ) );
                LG := List( LS, s -> Image( C.embedKernel, s ) );
                sectdata.preimages := List( LG, c -> Elements( c ) );
                sect := BacktrackSectionsJ( sectdata, ngR, [ ], 1 );
                for k in [ 1..Length( sect ) ] do
                    genim := sect[k];
                    xi := Cat1SectionByImages( C, genim );
                    ok1 := IsSection( xi );
                    if ok1 then
                        count := count + 1;
                    fi;
                od;
                Append( images, sect );
            fi;
        od;
        if ( c = 1 ) then
            regnum := count;
        fi;
        if ( XModPrintLevel > 1 ) then
            Print( "after class ", class, "\n" );
            Print( " there were ", count, " sections found.\n" );
        fi;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "\nThere were ", count, " sections found\n" );
    fi;
    sectrec := rec( );
    sectrec.areSections := true;
    sectrec.isReg := ( str = "reg" );
    sectrec.isAll := ( str = "all" );
    sectrec.regular := regnum;
    sectrec.genimageList := images;
    sectrec.generators := Copy( genR );
    sectrec.cat1 := C;
    sectrec.operations := Cat1SectionsOps;
    C.sections := sectrec;
    return sectrec;
end;

##############################################################################
##
#F  Sections                                  intermediate dispatcher function
##

Sections := function( C, str, how )

    local  exists, old, X, D, i, num, imder, imsec, chi, xi, sectrec;

    if not IsCat1( C ) then
        Error( "Parameter must be a cat1-group" );
    fi;
    if not ( IsString( str ) and str in [ "all", "reg" ] ) then
        Error( "Invalid all/reg string" );
    fi;
    if not ( IsString( how ) and ( how in [ "endo", "xmod", "dflt" ] ) ) then
        Error( "Invalid method string" );
    fi;
    if IsBound( C.sections ) then
        old := C.sections;
        if ( IsBound( old.isReg ) and old.isReg ) then
            exists := "reg";
        elif ( IsBound( old.isAll ) and old.isAll ) then
            exists := "all";
        else
            Error( "Invalid field  C.sections" );
        fi;
    else
        exists := "none";
    fi;
    if ( exists = str ) then
        return old;
    fi;
    if ( ( exists = "all" ) and ( str = "reg" ) ) then
        if ( XModPrintLevel > 1 ) then
            Print( "returning existing ALL sections\n" );
        fi;
        return old;
    fi;
    if ( how = "xmod" ) then
        X := XModCat1( C );
        D := BacktrackDerivations( X, str );
    fi;
    if ( ( exists = "none" ) and ( how = "dflt" ) ) then
        how := "endo";
    fi;
    if ( ( exists = "none" ) or ( how = "xmod" ) ) then
        sectrec := rec( );
        sectrec.areSections := true;
        # see if xmod.derivations exists
        if ( IsBound( C.xmod ) and IsBound( C.xmod.derivations ) ) then
            X := C.xmod;
            D := X.derivations;
            if ( AreDerivations( D ) and DerivationsSorted( D ) ) then
                imder := D.genimageList;
                if ( str = "reg" ) then
                    num := D.regular;
                else
                    num := Length( imder );
                fi;
                imsec := 0 * [1..num];
                for i in [1..num] do
                    chi := XModDerivationByImages( X, imder[i], true );
                    xi := SectionDerivation( chi );
                    imsec[i] := xi.genimages;
                od;
                sectrec.regular := D.regular;
                sectrec.isReg := ( str = "reg" );
                sectrec.isAll := ( str = "all" );
                sectrec.genimageList := imsec;
                sectrec.generators := Copy( C.range.generators );
                sectrec.cat1 := C;
                sectrec.operations := Cat1SectionsOps;
                C.sections := sectrec;
                return sectrec;
            fi;
        fi;
    fi;
    exists := ( exists = "reg" );
    sectrec := SectionsByEndoClasses( C, str, exists );
    return sectrec;    
end;

##############################################################################
##
#F  RegularSections              find all invertible sections for a cat1-group
##

RegularSections := function( arg )

    local  nargs, C, how, ok;

    nargs := Length( arg );
    ok := true;
    C := arg[1];
    if not ( ( nargs < 3 ) and IsCat1( C ) ) then
        ok := false;
    elif ( nargs = 2 ) then
        how := arg[2];
        if not ( IsString( how ) and ( how in [ "endo", "xmod" ] ) ) then
            ok := false;
        fi;
    else
        how := "dflt";
    fi;
    if not ok then
        Error( "Usage: RegularSections( C, [, \"endo\" | \"xmod\" ] );" );
    fi;
    return Sections( C, "reg", how );
end;

##############################################################################
##
#F  AllSections                             find all sections for a cat1-group
##

AllSections := function( arg )

    local  nargs, C, how, ok;

    nargs := Length( arg );
    ok := true;
    C := arg[1];
    if not ( ( nargs < 3 ) and IsCat1( C ) ) then
        ok := false;
    elif ( nargs = 2 ) then
        how := arg[2];
        if not ( IsString( how ) and ( how in [ "endo", "xmod" ] ) ) then
            ok := false;
        fi;
    else
        how := "dflt";
    fi;
    if not ok then
        Error( "Usage: AllSections( C, [, \"endo\" | \"xmod\" ] );" );
    fi;
    return Sections( C, "all", how );
end;

##############################################################################
##
#F  BacktrackDerivationsJ          recursive function for BacktrackDerivations
##

BacktrackDerivationsJ := function( X, subs, imrho, imchi, j, str )
    
    local  Xsrc, Xrng, derivgen, s, ok, rho, el, bdy, k, J, genJ,
           ord, r, aut, w, t, i, chi;

    Xsrc := X.source;
    Xrng := X.range;
    el := Elements(Xsrc);
    k := Length( Xrng.generators );
    bdy := X.boundary;
    if ( k < j ) then
        # complete list of images found: if a derivation, add to genimagesList
        imchi := ShallowCopy( imchi );
        derivgen := [ imchi ];
        chi := XModDerivationByImages( X, imchi );
        ok := IsDerivation( chi );
        if ( ok and ( str = "reg" ) ) then
            ok := IsRegular( chi );
        fi;
        if not ok then
            derivgen := [ ];
        fi;
    else
        J := subs[j];
        genJ := J.generators;
        derivgen := [ ];
        for s in el do
            imchi[ j ] := s;
            imrho[ j ] := genJ[j] * Image(bdy,s);
            rho := GroupHomomorphismByImages(J,Xrng,genJ,imrho);
            ok := IsHomomorphism(rho);
            if ok then
                r := genJ[j];
                ord := Order( Xrng, r );
                w := s;
                t := s;
                aut := Image( X.action, r );
                for i in [1..ord-1] do
                    t := Image( aut, t );
                    w := t * w;
                od;
                if ( w <> () ) then
                    ok := false;
                    if ( XModPrintLevel > 1 ) then
                        Print( "test fails at j,r,s,w = ", [j,r,s,w], "\n" );
                    fi;
                fi;
            fi;
            if ok then
                Append( derivgen,
                    BacktrackDerivationsJ( X, subs, imrho, imchi, j+1, str ) );
            fi;
            imrho := Sublist( imrho, [1..j] );
        od;
    fi;
    return derivgen;
end;

##############################################################################
##
#F  BacktrackDerivations            recursive construction for all derivations
##

BacktrackDerivations := function( X, str )

    local  R, len, genR, subs, images, sorted, derivrec;

    if not IsXMod( X ) then
        Error( "parameter must be a crossed module" );
    fi;
    if not ( IsString( str ) and str in [ "reg", "all" ] ) then
        Error( "Invalid reg|all string" );
    fi;
    R := X.range;
    genR := R.generators;
    len := Length( genR );
    subs := List( [1..len], i -> Subgroup( R, Sublist( genR, [1..i] ) ) ); 
    images := BacktrackDerivationsJ( X, subs, [], [], 1, str );
    derivrec := rec( );
    derivrec.areDerivations := true;
    derivrec.isReg := ( str = "reg" );
    derivrec.isAll := ( str = "all" );
    derivrec.genimageList := images;
    derivrec.operations := XModDerivationsOps;
    derivrec.xmod := X;
    derivrec.generators := Copy( genR );
    X.derivations := derivrec;
    if ( str = "reg" ) then
        sorted := DerivationsSorted( X );
    fi;
    return derivrec;
end;

##############################################################################
##
#F  DerivationsSorted          place regular derivations first in genimageList
##

DerivationsSorted := function( arg )

    local  D, reg, L, M, i, size, chi, pos;

    if not ( ( Length( arg ) = 1 ) and
             ( IsXMod( arg[1] ) or AreDerivations( arg[1] ) ) ) then
        Print( "Input must be a crossed module or derivations record\n" );
        return false;
    fi;
    if IsXMod( arg[1] ) then
        X := arg[1];
        if not IsBound( X.derivations ) then
            Print( "No derivations record exists for X\n" );
            return false;
        fi;
        D := X.derivations;
    else
        D := arg[1];
        X := D.xmod;
    fi;
    if ( IsBound( D.regular ) and ( D.regular > 0 ) ) then
        return true;
    fi;
    L := D.genimageList;
    size := Length( L );
    M := 0 * [1..size];
    for i in [1..size] do
        chi := XModDerivationByImages( X, L[i], true );
        M[i] := IsRegular( chi );
    od;
    SortParallel( M, L );
    pos := Position( M, false );
    if ( pos = false ) then
        D.regular := size;
    else
        D.regular := pos - 1;
    fi;
    return true;
end;

##############################################################################
##
#F  Derivations                               intermediate dispatcher function
##

Derivations := function( X, str, how )

    local  exists, old, C, S, i, num, imder, imsec, chi, xi, derivrec, sorted;

    if not IsXMod( X ) then
        Error( "First parameter must be a crossed module" );
    fi;
    if not ( IsString( str ) and str in [ "all", "reg" ] ) then
        Error( "Invalid all/reg string" );
    fi;
    if not ( IsString( how ) and ( how in [ "back", "cat1", "dflt" ] ) ) then
        Error( "Invalid method string" );
    fi;
    if IsBound( X.derivations ) then
        old := X.derivations;
        if ( IsBound( old.isReg ) and old.isReg ) then
            exists := "reg";
        elif ( IsBound( old.isAll ) and old.isAll ) then
            exists := "all";
        else
            Error( "Invalid field  X.derivations" );
        fi;
    else
        exists := "none";
    fi;
    if ( exists = str ) then
        return old;
    fi;
    if ( ( exists = "all" ) and ( str = "reg" ) ) then
        sorted := DerivationsSorted( X );
        if ( XModPrintLevel > 1 ) then
            Print( "returning existing ALL derivations\n" );
        fi;
        return old;
    fi;
    if ( how = "cat1" ) then
        C := Cat1XMod( X );
        if not IsBound( C.sections ) then
            S := SectionsByEndoClasses( C, str, false );
        fi;
    fi;
    if ( ( exists = "none" ) and ( how = "dflt" ) ) then
        how := "back";
    fi;
    if ( ( exists = "none" ) or ( how = "cat1" ) ) then
        derivrec := rec( );
        derivrec.areSections := true;
        # see if cat1.sections exists
        if ( IsBound( X.cat1 ) and IsBound( X.cat1.sections ) ) then
            C := X.cat1;
            S := C.sections;
            if AreSections( S ) then
                imsec := S.genimageList;
                if ( str = "reg" ) then
                    num := S.regular;
                else
                    num := Length( imsec );
                fi;
                imder := 0 * [1..num];
                for i in [1..num] do
                    xi := Cat1SectionByImages( C, imsec[i], true );
                    chi := DerivationSection( xi );
                    imder[i] := chi.genimages;
                od;
                derivrec.regular := S.regular;
                derivrec.isReg := ( str = "reg" );
                derivrec.isAll := ( str = "all" );
                derivrec.genimageList := imder;
                derivrec.generators := Copy( X.range.generators );
                derivrec.xmod := X;
                derivrec.operations := XModDerivationsOps;
                X.derivations := derivrec;
                return derivrec;
            fi;
        fi;
    fi;
    exists := ( exists = "reg" );
    derivrec := BacktrackDerivations( X, str );
    return derivrec;    
end;

##############################################################################
##
#F  RegularDerivations    find all invertible derivations for a crossed module
##

RegularDerivations := function( arg )

    local  nargs, X, how, ok;

    nargs := Length( arg );
    ok := true;
    X := arg[1];
    if not ( ( nargs < 3 ) and IsXMod( X ) ) then
        ok := false;
    elif ( nargs = 2 ) then
        how := arg[2];
        if not ( IsString( how ) and ( how in [ "back", "cat1" ] ) ) then
            ok := false;
        fi;
    else
        how := "dflt";
    fi;
    if not ok then
        Error( "Usage: RegularDerivations( C, [, \"back\" | \"cat1\" ] );" );
    fi;
    return Derivations( X, "reg", how );
end;

##############################################################################
##
#F  AllDerivations                   find all derivations for a crossed module
##

AllDerivations := function( arg )

    local  nargs, X, how, ok;

    nargs := Length( arg );
    ok := true;
    X := arg[1];
    if not ( ( nargs < 3 ) and IsXMod( X ) ) then
        ok := false;
    elif ( nargs = 2 ) then
        how := arg[2];
        if not ( IsString( how ) and ( how in [ "back", "cat1" ] ) ) then
            ok := false;
        fi;
    else
        how := "dflt";
    fi;
    if not ok then
        Error( "Usage: AllDerivations( C, [, \"back\" | \"cat1\" ] );" );
    fi;
    return Derivations( X, "all", how );
end;

##############################################################################
##
#F  WhiteheadMonoidTable( X )                 Table of products of derivations
##

WhiteheadMonoidTable := function( X )

    local  C, D, L, i, j, chi, images, J, M, size;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
    fi;
    D := AllDerivations( X );
    if IsBound( D.monoidTable ) then
        return D.monoidTable;
    fi;
    L := D.genimageList;
    size := Length( L );
    M := 0 * [1..size];
    C := 0 * [1..size];
    for i in [1..size] do
        C[i] := XModDerivationByImages( X, L[i], true );
        images := DerivationImages( C[i] );
    od;
    for i in [1..size] do
        J := 0 * [1..size];
        for j in [1..size] do
            chi := CompositeDerivation( C[i], C[j] );
            J[j] := Position( L, chi.genimages );
        od;
        M[i] := J;
    od;
    D.monoidTable := M;
    return M;
end;

###############################################################################
##
#F  WhiteheadGroupTable( X )           Table of products of regular derivations
##

WhiteheadGroupTable := function( X )

    local  C, D, L, i, j, chi, images, J, M, reg;

    if not IsXMod ( X ) then
        Error( " Input parameter must be a crossed module" );
    fi;
    D := RegularDerivations( X );
    if IsBound( D.groupTable ) then
        return D.groupTable;
    fi;
    reg := D.regular;
    L := D.genimageList;
    M := 0 * [1..reg];
    C := 0 * [1..reg];
    for i in [1..reg] do
        C[i] := XModDerivationByImages( X, L[i], true );
        images := DerivationImages( C[i] );
    od;
    for i in [1..reg] do
        J := 0 * [1..reg];
        for j in [1..reg] do
            chi := CompositeDerivation( C[i], C[j] );
            J[j] := Position( L, chi.genimages );
        od;
        M[i] := J;
    od;
    D.groupTable := M;
    return M;
end;

##############################################################################
##
#F  TestWhiteheadIsomorphism                  checks compatibility of products
##

TestWhiteheadIsomorphism := function( X )

    local  W, oW, EW, D, L, reg, i, j, p, r, ok, chi, chj, comp;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
    fi;

    ok := true;
    D := RegularDerivations( X );
    L := D.genimageList;
    reg := D.regular;
    W := WhiteheadPermGroup( X );
    oW := Size( W );
    EW := Elements( W );
    for i in [1..oW] do
        chi := XModDerivationByImages( X, L[i], true );
        for j in [1..oW] do
            chj := XModDerivationByImages( X, L[j], true );
            comp := CompositeDerivation( chi, chj );
            p := Position( L, comp.genimages );
            r := Position( EW, EW[i] * EW[j] );
            if ( p <> r ) then
                Print( "failure at ", i, ", ", j, "\n" );
                ok := false;
                return ok;
            fi;
        od;
    od;
    return ok;
end;

#############################################################################
##
#F  InverseDerivations      Finds all semigroup inverses XJ for derivation Xi
##                                                i.e.  XiXjXi=Xi & XjXiXj=Xj

InverseDerivations := function( chi )

    local  X, x, D, L, size, j, chj, chk, chl, chh, chg, inv;

    if not IsDerivation( chi ) then
        Error( "Parameter must be a derivation of a crossed module" );
    fi;
    X := chi.xmod;
    D := AllDerivations( X );
    L := D.genimageList;
    size := Length( L );
    inv := [ ];
    for j in [1..size] do
        chj := XModDerivationByImages( X, L[j], true );
        chk := CompositeDerivation( chi, chj );
        chl := CompositeDerivation( chk, chi );
        chh := CompositeDerivation( chj, chi );
        chg := CompositeDerivation( chh, chj );
        if ( ( chi.genimages = chl.genimages )
            and ( chg.genimages = chj.genimages ) ) then
            Add( inv, Position( L, chj.genimages ) );
        fi;
    od;
    return inv;
end;

##############################################################################
##
#F  ListInverseDerivations               List all inverses for each derivation
##

ListInverseDerivations := function( X )

    local  L, i, M, size, D, chi, inv;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
    fi;
    D := AllDerivations( X );
    L := D.genimageList;
    size := Length( L );
    M := 0 * [1..size];
    if IsBound( D.inverses ) then
        return D.inverses;
    fi;
    for i in [1..size] do
        chi := XModDerivationByImages( X, L[i], true );
        inv := InverseDerivations( chi );
        M[i] := inv;
    od;
    D.inverses := M;
    return M;
    end;

##############################################################################
##
#F  InnerDerivation               derivations of the form \eta_s(r) = s^r s^-1
##

InnerDerivation := function( X, s )
    
    local  Xact, Xsrc, Xrng, genrng, q, r, imeta, eta, nrng, a, j;

    if not ( IsXMod ( X ) ) then
        Error( "Input parameter must be a crossed module" );
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    Xact := X.action;
    if not ( s in Xsrc ) then
        Error( "Second parameter must be an element of the source of X" );
    fi;
    genrng := Xrng.generators;
    nrng := Length( genrng );
    imeta := 0 * [1..nrng];
    for j in [1..nrng] do
        r := genrng[j];
        a := Image( Xact, r );
        q := Image( a, s );
        imeta[j] := q * s^(-1);
    od;
    eta := XModDerivationByImages( X, imeta, true );
    return eta;
end;

##############################################################################
##
#F  ListInnerDerivations                            list all inner derivations
##

ListInnerDerivations := function( X )
    
    local  L, D, Xsrc, numsrc, elsrc, s, eta, i;

    if not ( IsXMod ( X ) ) then
        Error( "Input parameter must be a crossed module" );
    fi;
    if IsBound( X.derivations ) then
        D := X.derivations;
        if IsBound( X.derivations.innerImageList ) then
            return X.derivations.innerImageList;
        fi;
    fi;
    Xsrc := X.source;
    numsrc := Size( Xsrc );
    elsrc := Elements( Xsrc );
    L := 0 * [1..numsrc];
    for i in [1..numsrc] do
        s := elsrc[i];
        eta := InnerDerivation( X, s );
        L[i] := eta.genimages;
    od;
    if IsBound( X.derivations ) then
        D.innerImageList := L;
    fi;
    return L;
end;

##############################################################################
##
#F  InnerDerivationGroup                       Subgroup of the Whitehead Group
##

InnerDerivationGroup := function( X )

    local  D, L, W, eW, S, genS, im, pos, inn, inner, phi, ok;

    if not ( IsXMod ( X ) ) then
        Error( "Input parameter must be a crossed module" );
    fi;
    D := RegularDerivations( X );
    L := D.genimageList;
    W := WhiteheadPermGroup( X );
    eW := Elements( W );
    if IsBound( W.inner ) then
        return W.inner;
    fi;
    S := X.source;
    genS := S.generators;
    inn := List( genS, s -> InnerDerivation( X, s ) );
    pos := List( inn, chi -> Position( L, chi.genimages ) );
    im := List( pos, i -> eW[i] );
    inner := Subgroup( W, im );
    phi := GroupHomomorphismByImages( S, inner, genS, im );
    ok := IsHomomorphism( phi );
    if not ok then
        Error( "homomorphism failure,  S -> inner <= W" );
    fi;
    W.inner := inner;
    W.innerHomomorphism := phi;
    return inner;
end;

##############################################################################
##
#F  SourceEndomorphismDerivation               produces a homomorphism  \sigma
##                                                    from a derivation \chi

SourceEndomorphismDerivation := function( arg )

    local  nargs, usage, ok, X, im, chi, Xsrc, Xrng, Xbdy, elsrc, elrng,
           images, srcgen, srcrng, imbdy, pos, imchi, imsigma, sigma; 

    nargs := Length( arg );
    usage := "\nUsage: SourceEndomorphismDerivation( chi | X, im );";
    ok := true;
    if ( nargs = 1 ) then
        chi := arg[1];
    elif ( nargs = 2 ) then
        X := arg[1];
        im := arg[2];
        chi := XModDerivationByImages( X, im );
    else
        ok := false;
    fi;
    ok := ( ok and IsDerivation( chi ) );
    if not ok then
        Print( usage );
        return false;
    fi;
    if ( nargs = 1 ) then
        X := chi.xmod;
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    Xbdy := X.boundary;
    elsrc := Elements( Xsrc );
    elrng := Elements( Xrng );
    images := DerivationImages( chi );
    srcgen := Xsrc.generators;
    srcrng := [ 1..Length( srcgen ) ];
    imbdy := List( srcgen, s -> Image( Xbdy, s ) );
    pos := List( imbdy, r -> Position( elrng, r ) );
    imchi := List( pos, j -> images[j] );
    imsigma := List( srcrng, j -> srcgen[j] * elsrc[ imchi[j] ] );
    sigma := GroupHomomorphismByImages( Xsrc, Xsrc, srcgen, imsigma );
    return sigma;
end;

##############################################################################
##
#F  RangeEndomorphismDerivation                produces a homomorphism  \sigma
##                                                   from a derivation  \chi

RangeEndomorphismDerivation := function( arg )

    local  nargs, usage, ok, chi, X, im, Xrng, Xbdy,
           rnggen, rngrng, imbdy, imchi, imrho, rho; 

    nargs := Length( arg );
    usage := "\nUsage: RangeEndomorphismDerivation( chi | X, im );";
    ok := true;
    if ( nargs = 1 ) then
        chi := arg[1];
    elif ( nargs = 2 ) then
        X := arg[1];
        im := arg[2];
        chi := XModDerivationByImages( X, im );
    else
        ok := false;
    fi;
    ok := ( ok and IsDerivation( chi ) );
    if not ok then
        Print( usage );
        return false;
    fi;
    if ( nargs = 1 ) then
        X := chi.xmod;
    fi;
    Xrng := X.range;
    Xbdy := X.boundary;
    rnggen := chi.generators;
    rngrng := [ 1..Length( rnggen ) ];
    imchi := chi.genimages;
    imbdy := List( imchi, s -> Image ( Xbdy, s ) );
    imrho := List( rngrng, i -> rnggen[i] * imbdy[i] );
    rho := GroupHomomorphismByImages( Xrng, Xrng, rnggen, imrho );
    return rho;
end;

##############################################################################
##
#F  XModEndomorphismDerivation              chi -> XModMorphism( [sigma,rho] )
##

XModEndomorphismDerivation := function( arg )

    local  nargs, usage, hom, sigma, rho, ok, X, im, chi;

    nargs := Length( arg );
    usage := "\nUsage: XModEndomorphismDerivation( chi | X, im );";
    ok := true;
    if ( nargs = 1 ) then
        chi := arg[1];
    elif ( nargs = 2 ) then
        X := arg[1];
        im := arg[2];
        chi := XModDerivationByImages( X, im );
    else
        ok := false;
    fi;
    ok := ( ok and IsDerivation( chi ) );
    if not ok then
        Print( usage );
        return false;
    fi;
    if ( nargs = 1 ) then
        X := chi.xmod;
    fi;
    sigma := SourceEndomorphismDerivation( chi );
    rho := RangeEndomorphismDerivation( chi );
    hom := XModMorphism( X, X, [ sigma, rho ] );
    ok := IsXModMorphism( hom );
    return hom;
end;

##############################################################################
##
#F  ImageAutomorphismDerivation       image of an automorphism on a derivation
##

ImageAutomorphismDerivation := function( mor, chi )

    local  X, R, genR, imj, rho, imrho, sigma, invrho,
           rngrng, k, r, rr, crr, chj;

    if not ( IsXModMorphism( mor ) and IsDerivation( chi ) ) then 
        Error( "Parameters must be an xmod morphism and a derivation" );
    fi;
    X := mor.source;
    if not ( ( X = mor.range ) and ( X = chi.xmod ) ) then
        Error( "either mor.range or chi.xmod <> X" );
    fi;
    sigma := mor.sourceHom;
    rho := mor.rangeHom;
    R := X.range;
    genR := chi.generators;
    rngrng := [ 1..Length( genR ) ];
    imrho := rho.genimages;
    invrho := GroupHomomorphismByImages( R, R, imrho, genR );
    imj := 0 * rngrng;
    for k in rngrng do
        r := genR[k];
        rr := Image( invrho, r );
        crr := DerivationImage( chi, rr );
        imj[k] := Image( sigma, crr );
    od;
    chj := XModDerivationByImages( X, imj, true );
    return chj;
end;

##############################################################################
##
#F  SourceEndomorphismSection                     section xi => gamma : G -> G
##

SourceEndomorphismSection := function( arg )

    local  nargs, usage, C, im, xi, G, genG, R, t, h, e, eh, xit, ehxit, xih,
           gamma, imgamma, ok;

    nargs := Length( arg );
    usage := "\nUsage: SourceEndomorphismSection( xi | C, im );\n";
    ok := true;
    if ( nargs = 1 ) then
        xi := arg[1];
        C := xi.cat1;
    elif ( nargs = 2 ) then
        C := arg[1];
        im := arg[2];
        xi := Cat1SectionByImages( C, im );
    else
        ok := false;
    fi;
    if not IsCat1( C ) then
        Error( "first parameter must be a cat1-group" );
    fi;
    ok := ( ok and IsSection( xi ) );
    if not ok then
        Print( usage );
        return false;
    fi; 
    G := C.source;
    R := C.range;
    t := C.tail;
    h := C.head;
    e := C.embedRange;
    eh := CompositionMapping( e, h );
    xit := CompositionMapping( xi, t );
    xih := CompositionMapping( xi, h );
    ehxit := CompositionMapping( eh, xit );
    genG := G.generators;
    imgamma := List( genG, g -> Image( ehxit, g ) * Image( xit, g^-1 ) * g
                                    * Image( eh, g^-1 ) * Image( xih, g ) );
    gamma := GroupHomomorphismByImages( G, G, genG, imgamma );
    if not IsEndomorphism( gamma ) then
        Error( "gamma fails to be an endomorphism" );
    fi;
    return gamma;
end;

##############################################################################
##
#F  RangeEndomorphismSection                        section xi => rho : R -> R
##

RangeEndomorphismSection := function( arg )

    local  nargs, usage, xi, C, im, ok, chi, rho;

    nargs := Length( arg );
    usage := "\nUsage: SourceEndomorphismSection( xi | C, im );\n";
    ok := true;
    if ( nargs = 1 ) then
        xi := arg[1];
        C := xi.cat1;
    elif ( nargs = 2 ) then
        C := arg[1];
        im := arg[2];
        xi := Cat1SectionByImages( C, im );
    else
        ok := false;
    fi;
    if not IsCat1( C ) then
        Error( "first parameter must be a cat1-group" );
    fi;
    ok := ( ok and IsSection( xi ) );
    if not ok then
        Print( usage );
        return false;
    fi; 
    chi := DerivationSection( xi );
    rho := RangeEndomorphismDerivation( chi );
    return rho;
end;

##############################################################################
##
#F  Cat1EndomorphismSection                  xi -> Cat1Morphism( [gamma,rho] )
##

Cat1EndomorphismSection := function( arg )

    local  nargs, usage, hom, gamma, rho, ok, C, im, xi;

    nargs := Length( arg );
    usage := "\nUsage: Cat1EndomorphismSection( xi | C, im );";
    ok := true;
    if ( nargs = 1 ) then
        xi := arg[1];
    elif ( nargs = 2 ) then
        C := arg[1];
        im := arg[2];
        xi := Cat1SectionByImages( C, im );
    else
        ok := false;
    fi;
    ok := ( ok and IsSection( xi ) );
    if not ok then
        Print( usage );
        return false;
    fi;
    if ( nargs = 1 ) then
        C := xi.cat1;
    fi;
    gamma := SourceEndomorphismSection( xi );
    rho := RangeEndomorphismSection( xi );
    hom := Cat1Morphism( C, C, [ gamma, rho ] );
    ok := IsCat1Morphism( hom );
    return hom;
end;

##############################################################################
##
#F  TableSourceEndomorphismDerivations  produces table of homomorphisms \sigma
##                                         from a table of derivations  \chi

TableSourceEndomorphismDerivations := function( X )

    local  D, L, len, M, i, chi, sigma;

    if not ( IsXMod ( X ) ) then
        Error( "First input parameter must be a crossed module" );
    fi;
    if not IsBound( X.derivations ) then
        Print( "No X.derivations record exists.\n" );
        return [ ];
    fi;
    D := X.derivations;
    L := D.genimageList;
    len := Length( L );
    M := 0 * [1..len];
    for i in [1..len] do
        chi := XModDerivationByImages( X, L[i], true );
        sigma := SourceEndomorphismDerivation( chi );
        M[i] := sigma.genimages;
    od;
    return M;
end;
 
##############################################################################
##
#F  TableRangeEndomorphismDerivations     produces table of homomorphisms \rho
##                                           from a table of derivations  \chi

TableRangeEndomorphismDerivations := function( X )
  
    local  D, L, len, M, i, chi, rho;

    if not ( IsXMod ( X ) ) then
        Error( "First input parameter must be a crossed module" );
    fi;
    if not IsBound( X.derivations ) then
        Print( "No X.derivations record exists.\n" );
        return [ ];
    fi;
    D := X.derivations;
    L := D.genimageList;
    len := Length( L );
    M := 0 * [1..len];
    for i in [1..len] do
        chi := XModDerivationByImages( X, L[i], true );
        rho := RangeEndomorphismDerivation( chi );
        M[i] := rho.genimages;
    od;
    return M;
end;

##  end of file  deriv.g
