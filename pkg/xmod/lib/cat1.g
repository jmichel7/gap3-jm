##############################################################################
##
##  cat1.g                   for GAP 3.4                      version 10/ 1/97
##
##############################################################################
##
#A  cat1.g                   GAP library                         Chris Wensley
#A                                                                   Murat Alp
#Y  Copyright
##
##  This file contains functions that manipulate cat-1 groups
##  and associated constructions.
##
#H  $Log: cat1.g,v $
#H  Revision 1.1  1997/03/27 13:35:42  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

##############################################################################
##
#F  IsCat1                                  checks conditions for a cat1-group
##

IsCat1 := function( C )

    local  Csrc, Crng, gensrc, genrng, x, e, f, t, h,
           eh, he, et, te, teh, het, kert, kerh, kerth;

    if not IsRec( C ) then
        return false;
    fi;
    if IsBound ( C.isCat1 ) then
        return C.isCat1;
    fi;
    if (   not ( IsBound( C.source ) and IsBound( C.source.generators ) )
        or not ( IsBound( C.range ) and IsBound( C.range.generators ) ) ) then
        Error( "Source and/or Range fields not defined." );
    fi; 
    Csrc := C.source;
    Crng := C.range;
    gensrc := Csrc.generators;
    genrng := Crng.generators;
    if ( not IsBound( C.tail )  or  not IsBound( C.head ) ) then
         Print( "head and tail homomorphisms not defined\n" );
         return false;
    fi;
    h := C.head;
    t := C.tail;
    if not IsHomomorphism( t ) then
        Print( "C.tail is not a homomorphism\n" );
        return false;
    fi;
    if not IsHomomorphism ( h ) then
        Print( "C.head is not a homomorphism\n" );
        return false;
    fi;
    if not ( IsBound( C.embedRange ) and IsBound( C.embedKernel ) ) then
         Print( "embedding homomorphisms not defined\n" );
         return false;
    fi;
    e := C.embedRange;
    f := C.embedKernel;

    # checking the first condition of cat-1 group
    eh := CompositionMapping( e, h );
    he := CompositionMapping( h, e );
    et := CompositionMapping( e, t );
    te := CompositionMapping( t, e );
    teh := CompositionMapping( t, eh );
    if not ( teh = h ) then
        Print( "condition  teh = h  is not satisfied \n" );
        return false;
    fi;
    het := CompositionMapping( he, t );
    if not ( het = t ) then
        Print( "condition  het = t, is not satisfied \n" );
        return false;
    fi;

    # checking the second condition of cat-1 group
    kert := Kernel( t );
    kerh := Kernel( h );
    kerth := CommutatorSubgroup( kert, kerh );
    if not ( Size( kerth ) = 1 ) then
        Print("condition  [kert,kerh] = 1  is not satisfied \n");
        return false;
    fi;
    if not ( ( f.source = C.kernel ) and ( f.range = C.source ) ) then
        Print( "Warning: C.embedKernel incorrectly defined?\n" );
    fi; 
    if not ( IsBound( C.boundary ) and IsBound( C.kernel) ) then
        Print( "Warning: C.boundary and/or C.kernel not yet defined!\n" );
    fi;

    C.isCat1 := true;
    return true;
end;

##############################################################################
##
#F  Cat1Ops                                       cat1-group record operations
##

Cat1Ops := OperationsRecord( "Cat1Ops", DomainOps );

##############################################################################
##
#F  Cat1Ops.\=                                  define equality of cat1-groups
##

Cat1Ops.\= := function( C, D )

    local  isEql;

    if not ( IsCat1( C ) and IsCat1( D ) ) then
        Error( "C,D not cat1-groups" );
    fi;

    isEql := (    ( C.source = D.source )
              and ( C.range = D.range )
              and ( C.tail = D.tail )
              and ( C.head = D.head )   );
    return isEql;
end;

##############################################################################
##
#F  Cat1Ops.Size                         [ Size( C.source ), Size( C.range ) ]
##

Cat1Ops.Size := function( C )

    local  L;
 
    if not ( IsCat1( C ) ) then
        Error( "Group must be a cat1 group" );
    fi;
    if IsBound( C.size ) then
        return C.size;
    fi;
    L := [ Size( C.source ), Size( C.range ) ];
    C.size := L;
    return L;
end;

##############################################################################
##
#F  Cat1Ops.Elements             [ Elements( C.source ), Elements( C.range ) ]
##

Cat1Ops.Elements := function( C )

    local  eC;
 
    if not ( IsCat1( C ) ) then
        Error( "Group must be a cat1 group" );
    fi;
    if IsBound( C.elements ) then
        return C.elements;
    fi;
    eC := [ Elements( C.source ), Elements( C.range ) ];
    C.elements := eC;
    return eC;
end;

##############################################################################
##
#F  Cat1Ops.Print                       special print commands for cat1-groups
##

Cat1Ops.Print := function( C )

    local name;

    if IsBound( C.name ) then
        name := C.name;
    else
        name := "[?==>?]";
    fi;
    Print( "cat1-group ", name, " " );
end;

##############################################################################
##
#F  Cat1Print                                    print details of a cat1-group
##

Cat1Print := function( C )

    local  name, Csrc, Crng, head, tail, gensrc, genrng, genker, imbdy,
           imrng, imker, imhead, imtail, Arec;

    if not ( IsBound( C.isCat1 ) and C.isCat1 ) then
        Error( "Cat1Print can only print a cat1-group" );
    fi;
    if not IsBound( C.name ) then
        name := "[?==>?]";
    else
        name := C.name;
    fi;
    Csrc := C.source;
    Crng := C.range;
    gensrc := Csrc.generators;
    genrng := Crng.generators;
    genker := C.kernel.generators;
    imtail := List( gensrc, x -> Image( C.tail, x ) );
    imhead := List( gensrc, x -> Image( C.head, x ) );
    imrng := List( genrng, x -> Image( C.embedRange, x ) );
    imbdy := List( genker, x -> Image( C.boundary, x ) );
    imker := List( genker, x -> Image( C.embedKernel, x ) );

    Print( "\ncat1-group ", name, " :- \n" );
    Print( ": source group has generators:\n" );
    Print( "  ", Csrc.generators, "\n" );
    Print( ":  range group has generators:\n" );
    Print( "  ", Crng.generators, "\n" );
    Print( ": tail homomorphism maps source generators to:\n" );
    Print( "  ", imtail, "\n" );
    Print( ": head homomorphism maps source generators to:\n" );
    Print( "  ", imhead, "\n" );
    Print( ": range embedding maps range generators to:\n" );
    Print( "  ", imrng, "\n" );
    Print( ": kernel has generators:\n" );
    Print( "  ", genker, "\n" );
    Print( ": boundary homomorphism maps generators of kernel to:\n" );
    Print( "  ", imbdy, "\n" );
    Print( ": kernel embedding maps generators of kernel to:\n" );
    Print( "  ", imker, "\n" );
    if IsBound( C.xmod ) then
        Print( ": associated crossed module is ", C.xmod, "\n" );
    fi;
    if IsBound( C.sections ) then
        Print( ": ", C.sections, "\n" );
    fi;
    if IsBound( C.actorSquare ) then
        Arec := C.actorSquare;
        Print( ": ActorSquare: [" );
        if IsBound( Arec.WhiteheadPermGroup ) then
                  Print( "WhiteheadPermGroup, " ); fi;
        if IsBound( Arec.automorphismPermGroup ) then
                  Print( "AutomorphismPermGroup, " ); fi;
        if IsBound( Arec.Whitehead ) then Print( "Whitehead, " ); fi;
        if IsBound( Arec.Lue ) then Print( "Lue, " ); fi;
        if IsBound( Arec.Norrie ) then Print( "Norrie, " ); fi;
        if IsBound( Arec.actor ) then Print( "Actor, " ); fi;
        if IsBound( Arec.innerActor ) then Print( "InnerActor, " ); fi;
        if IsBound( Arec.innerMorphism ) then Print( "InnerMorphism, " ); fi;
        if IsBound( Arec.centre ) then Print( "Centre, " ); fi;
        Print( "]\n" );
    fi;
    Print( "\n" );
end;

##############################################################################
##
#F  Cat1Name              construct a name for a cat1-group or a morphism
##

Cat1Name := function( C )

    local  nsrc, nrng, name, mor;

    if not IsCat1( C ) then
        Error( "input parameter must be a cat-1 group\n" );
    fi;
    if IsBound( C.source.name ) then
        nsrc := C.source.name;
    else
        nsrc := "?";
    fi;
    if IsBound( C.range.name ) then
        nrng := C.range.name;
    else
        nrng := "?";
    fi;
    name := Concatenation( "[", nsrc, " ==> ", nrng, "]" );
    C.name := name;
    return name;
end;

##############################################################################
##
#F  Cat1                Construct a Cat1Group from tail and head morphisms
##

Cat1 := function( arg )

    local  nargs, usage, G, t, h, e, genG, imh, imt, hres, tres, R,
           kert, kerh, kerth, kertgen, bdy, imbdy, incR, incK, C, ok;
    
    nargs := Length( arg );
    usage := "\nUsage:  Cat1( G, tail, head [, embed] );\n";
    if ( ( nargs < 3 ) or ( nargs > 4 ) ) then
        Print( usage );
        return false;
    fi;
    G := arg[1];
    if not IsPermGroup( G ) then
        Error( "first parameter must be a permutation group" );
    fi;
    genG := G.generators;
    # Check that  t,h  are homomorphisms
    t := arg[2];
    h := arg[3];
    if ( not IsHomomorphism( t ) or not IsHomomorphism( h )
           or ( t.source <> G ) or ( h.source <> G )
           or ( t.range <> h.range ) )  then
        Print( usage );
        return false;
    fi;
    imh := List( genG, x -> Image( h, x ) );
    imt := List( genG, x -> Image( t, x ) );
    if IsSurjective( t ) then
        R := t.range;
    else
        R := Subgroup( G, imt );
    fi;
    if ( IsBound( G.name ) and not IsBound( R.name ) ) then
        R.name := "R";
    fi;
    kert := Kernel ( t );
    kerh := Kernel ( h );
    if ( nargs = 4 ) then
        incR := arg[4];
        if ( not IsHomomorphism( incR ) or ( incR.source <> R )
               or ( incR.range <> G ) )  then
            Print( usage );
            return false;
        fi;
    else
        incR := InclusionMorphism( R, G );
    fi;
    incK := InclusionMorphism( kert, G );
    hres := GroupHomomorphismByImages( G, R, genG, h.genimages );
    tres := GroupHomomorphismByImages( G, R, genG, t.genimages );
    kertgen := kert.generators;
    imbdy := List( kertgen, x -> Image( h, x) );
    bdy := GroupHomomorphismByImages( kert, R, kertgen, imbdy );
    
    C := rec ( );
    C.source := G;
    C.range := R;
    C.tail := tres;
    C.head := hres;
    C.embedRange := incR;
    C.kernel := kert;
    C.boundary := bdy;
    C.embedKernel := incK;
    ok := IsCat1( C );
    if ok then
        C.isDomain := true;
        C.operations := Cat1Ops;
        C.name := Cat1Name( C );
    fi;
    if ( G = R ) then
        C.isIdentity := true;
    fi;
    return C;
end;

##############################################################################
##
#F  Cat1XMod                  construct the cat1-group associated to a crossed
##                                   module with a permutation group as source
##

Cat1XMod := function( X )

    local  Xsrc, Xrng, Xact, Xbdy, S, genS, gensrc, genrng, imbdy,
           P, G, genG, t, h, e, f, imt, imh, ime, imf, C;
  
    if IsBound(X.cat1) then
        return X.cat1;
    fi;
    if not IsXMod( X ) then
        Error( "X is not a crossed module" );
    fi;
    
    Xsrc := X.source;
    Xrng := X.range;
    Xact := X.action;
    Xbdy := X.boundary;
    if XModOps.IsTrivialAction( X ) then
        G := DirectProduct( Xrng, Xsrc );
        if ( IsBound( Xsrc.name ) and IsBound( Xrng.name ) ) then
            G.name := Concatenation( Xrng.name, Xsrc.name );
        fi;
        genG := G.generators;
        gensrc := Xsrc.generators;
        genrng := Xrng.generators;
        imbdy := List( gensrc, s -> Image( X.boundary, s ) );
        imt := Concatenation( genrng, List( gensrc, s -> () ) );
        imh := Concatenation( genrng, imbdy );
        t := GroupHomomorphismByImages( G, Xrng, genG, imt );
        h := GroupHomomorphismByImages( G, Xrng, genG, imh );
        e := GroupHomomorphismByImages( Xrng, G, genrng, genrng );
        imf := List( gensrc, s -> s^G.perms[2] );
        f := GroupHomomorphismByImages( Xsrc, G, gensrc, imf );
    else
        S := SemidirectProduct( Xrng, Xact, Xsrc );
        if ( IsBound( Xsrc.name ) and IsBound( Xrng.name ) ) then
             S.name := Concatenation ( Xrng.name, " |X ", Xsrc.name );
        else
             S.name := "SDP";
        fi;
        genS := S.generators;
        imt := List( genS, x -> x.element[1] );
        imh := List( genS, x -> x.element[1] * Image( Xbdy, x.element[2]) );
        P := SemidirectPair( S );
        G := P.perm;
        genG := G.generators;
        t := GroupHomomorphismByImages( G, Xrng, genG, imt );
        h := GroupHomomorphismByImages( G, Xrng, genG, imh );
        e := Embedding( Xrng, S, 1 ) * P.s2p;
        f := Embedding( Xsrc, S, 2 ) * P.s2p;
        ime := List( e.source.generators, r -> Image( e, r ) );
        imf := List( f.source.generators, s -> Image( f, s ) );
        e := GroupHomomorphismByImages( Xrng, G, Xrng.generators, ime );
        f := GroupHomomorphismByImages( Xsrc, G, Xsrc.generators, imf );
        G.pair := P;
    fi;

    C := rec ( );
    C.source := G;
    C.range := Xrng;
    C.tail := t;
    C.head := h;
    C.embedRange := e;
    C.kernel := Xsrc;
    C.boundary := Xbdy;
    C.embedKernel := f;
    C.isDomain := true;
    C.operations := Cat1Ops;
    Cat1Name(C);
    X.cat1 := C;
    C.xmod := X;
    return C;
end; 

##############################################################################
##
#F  SemidirectCat1XMod        construct the cat1-group associated to a crossed
##                            module with a semidirect product group as source

SemidirectCat1XMod := function( X )

    local  Xsrc, Xrng, Xact, Xbdy, G, genG, h, t, e, f, C;
  
    if IsBound(X.cat1) then
        return X.cat1;
    fi;
    if not IsXMod( X ) then
        Error( "X is not a crossed module" );
    fi;
    
    Xsrc := X.source;
    Xrng := X.range;
    Xact := X.action;
    Xbdy := X.boundary;
    G := SemidirectProduct( Xrng, Xact, Xsrc );
    if ( IsBound( Xsrc.name ) and IsBound( Xrng.name ) ) then
         G.name := Concatenation ( Xrng.name, " |X ", Xsrc.name );
    else
         G.name := "G";
    fi;
    genG := G.generators;
    t := GroupHomomorphismByImages( G, Xrng, genG,
                                           List( genG, x -> x.element[1] ) );
    h := GroupHomomorphismByImages( G, Xrng, genG,
              List( genG, x -> x.element[1] * Image( Xbdy, x.element[2]) ) );
    e := Embedding( Xrng, G, 1 );
    f := Embedding( Xsrc, G, 2 );

    C := rec ( );
    C.source := G;
    C.range := Xrng;
    C.tail := t;
    C.head := h;
    C.embedRange := e;
    C.kernel := Xsrc;
    C.boundary := Xbdy;
    C.embedKernel := f;
    C.isDomain := true;
    C.operations := Cat1Ops;
    Cat1Name(C);
    X.cat1 := C;
    C.xmod := X;
    return C;
end; 

##############################################################################
##
#F  XModCat1           construct the crossed module associated to a cat1-group
##
 
XModCat1 := function( C )
 
    local  Csrc, Crng, gensrc, genrng, genker, kert, innaut, autgen,
           imautgen, idkert, a, aut, act, phi, j, r, X, Cek, Cer;
     
    if not IsCat1( C ) then
        Error( "C is not a cat1-group." );
    fi; 
    if IsBound( C.xmod ) then
        return C.xmod;
    fi;
    Csrc := C.source;
    Crng := C.range;
    Cer := C.embedRange;
    Cek := C.embedKernel;
    kert := C.kernel;
    if ( ( not IsBound( kert.name ) ) and IsBound( C.name ) ) then
        kert.name := Concatenation( "ker(", C.name, ")" );
    fi;
    gensrc := Csrc.generators;
    genrng := Crng.generators;
    genker := kert.generators;

    if ( IsBound( C.isIdentity ) and C.isIdentity ) then
        # X has trivial source and action
        aut := Group( InclusionMorphism( kert, kert ) );
        aut.name := "triv_aut";
        aut.identity.genimages := genker;
        act := ZeroMorphism( Crng, aut );
        act.name := "zero";
    else
        autgen := [ ];
        for r in genrng do
            imautgen := List( genker, s -> Image( Cek, s ) );
            imautgen := List( imautgen, g -> g^( Image( Cer, r ) ) );
            imautgen := List( imautgen,
                              g -> PreImagesRepresentative( Cek, g ) );
            a := GroupHomomorphismByImages( kert, kert, genker, imautgen );
            Add( autgen, a );
        od;
        idkert := InclusionMorphism( kert, kert );
        aut := Group( autgen, idkert );
        act := GroupHomomorphismByImages( Crng, aut, genrng, autgen );
        if IsBound( kert.name ) then
            aut.name := Concatenation( "innaut(", kert.name, ")" );
        else
            aut.name := "aut";
        fi;
        if not IsHomomorphism( act ) then
            Error( "act is not a homomorphism" );
        fi;
        if IsBound( Crng.name ) then
            act.name := Concatenation( "act(", Crng.name, ")" );
        else
            act.name := "act";
        fi;
    fi;

    X := XMod( C.boundary, act );
    C.xmod := X;
    X.cat1 := C;
    return X;
end;

###############################################################################
##
#F  ConjugationCat1                          construct cat1-group  G |X N ==> G
##

ConjugationCat1 := function( arg )

    local  nargs, S, R, X, C, Csrc;

    nargs := Length( arg );
    S := arg[1];
    if ( nargs = 1 ) then
        R := S;
    else
        R := arg[2];
    fi;
    if not ( IsPermGroup( S ) and IsPermGroup( R ) ) then
        Error( "the 2 arguments must both be PermGroups" );
    fi;
    if not IsSubgroup( S, R ) then
        Error( "Second parameter must be a normal subgroup of the first" );
    fi;

    X := ConjugationXMod( S, R );
    C := Cat1XMod( X );
    Csrc := C.source;
    if ( IsBound( R.name ) and IsBound( S.name ) ) then
        Csrc.name := Concatenation( S.name, " |X ", R.name );
    fi;
    return C;
end;

##############################################################################
##
#F  DirectProductCat1                    t&h1 x t*h2 : (G1 x G2) --> (R1 x R2)
##

DirectProductCat1 := function( C, D )

    local  C, D, B, G, R, genG, genR, t, h, imt, imh, e, ime, name,
           im1, im2, perm1, perm2, em1, em2;

    if (    not ( IsBound( C.isCat1 ) and C.isCat1 )
         or not ( IsBound( D.isCat1 ) and D.isCat1 ) ) then
        Error( "Both parameters must be cat1-groups" );
    fi;
    G := DirectProduct( C.source, D.source );
    G.name := Concatenation( C.source.name, "x", D.source.name );
    R := DirectProduct( C.range, D.range );
    R.name := Concatenation( C.range.name, "x", D.range.name );
    genG := G.generators;
    genR := R.generators;
    name := Concatenation( C.name, "x", D.name );

    perm1 := R.perms[1];
    perm2 := R.perms[2];
    im1 := C.tail.genimages;
    im2 := D.tail.genimages;
    em1 := List( im1, r -> r^perm1 );
    em2 := List( im2, r -> r^perm2 );
    imt := Concatenation( em1, em2 );
    t := GroupHomomorphismByImages( G, R, genG, imt );

    im1 := C.head.genimages;
    im2 := D.head.genimages;
    em1 := List( im1, r -> r^perm1 );
    em2 := List( im2, r -> r^perm2 );
    imh := Concatenation( em1, em2 );
    h := GroupHomomorphismByImages( G, R, genG, imh );

    perm1 := G.perms[1];
    perm2 := G.perms[2];
    im1 := C.embedRange.genimages;
    im2 := D.embedRange.genimages;
    em1 := List( im1, g -> g^perm1 );
    em2 := List( im2, g -> g^perm2 );
    ime := Concatenation( em1, em2 );
    e := GroupHomomorphismByImages( R, G, genR, ime );

    B := Cat1( G, t, h, e );
    B.isDirectProduct := true;
    return B;
end;

##############################################################################
##
#F  AllCat1sJ                            recursive function called by AllCat1s
##

AllCat1sJ := function( len, k, tup, genim, j, G )

    local  ok, nextgen, num, g, phi, genH, imH, H, genG, elG;
    
    if k = 0 then
        tup := ShallowCopy( tup );
        genim := ShallowCopy( genim );
        nextgen := [ tup ];
    else
        elG := Elements( G );
        genG := G.generators;
        nextgen := [ ];
        for num in len do
            tup[ j ] := num;
            genim[ j ] := elG[num];
            genH := Sublist( genG, [1..j] );
            imH := Sublist( genim, [1..j] );
            H := Subgroup( G, genH );
            phi := GroupHomomorphismByImages( H, G, genH, imH );
            ok := IsHomomorphism( phi );
            if ok then
                Append( nextgen, 
                        AllCat1sJ( len, k-1, tup, genim, j+1, G ) );
            fi;
        od;
    fi;
    return nextgen;
end;

#############################################################################
##
#F  BacktrackEndomorphismImages                      ??? what is this for ???
##

BacktrackEndomorphismImages := function( G )

    local  eG, oG, genG, ngG, set, im, num;

    genG := G.generators;
    eG := Elements( G );
    ngG := Length( genG );
    oG := Length( eG );
    set := [ 1..oG ];
    im := AllCat1sJ( set, ngG, [ ], [ ], 1, G );
    num := Length( im );
    return im;
end;

##############################################################################
##
#F  AllCat1s                        all cat1-structures on a permutation group
##

AllCat1s := function( G )

    local  el, genG, len, ngG, set, tup, lent, L, image, gens, phi, im, ok, 
           D, C, lenl, K, imim, images, t, h, prod, cond, kert, kerh, kerth,
           CC, g, i, j, k, AICG, imagelist, pos, pos1, Llist, lenim,
           homnum, homcount, catnum, A, imageclass, classnum, class,
           found, K, genK, alpha, subclass, gen, K2, catcount, classlen,
           bdy, kerbdy, okaut, okend, name, Endo, Enum, classes, LI, c, Q, R,
           conj, proj, Cnum, aut, Anum, autos, iso, comp, im, rho, psi, phi,
           Cnum, g, gim, Inum, Rnum, latlen, reps, idem, ilen;

    if not ( IsBound( G.isPermGroup ) and G.isPermGroup ) then
        Error( "parameter must be a permutation group" );
    fi;
    if IsBound( G.name ) then
        name := G.name;
    else
        name := "G";
    fi;
    genG := G.generators;
    ngG := Length( genG );
    okaut := IsBound( G.automorphismGroup );
    if not okaut then
        Print( "Calculating the automorphism group of ", name, ".\n" );
    fi;
    A := AutomorphismGroup( G );
    okend := IsBound( G.endomorphismClasses );
    if not okend then
        Print( "Calculating non-trivial endomorphism classes,\n" );
        Print( "which satisfy the condition  N meet H = id.\n" );
    fi;
    Endo := EndomorphismClasses( G, 3 );
    classes := Endo.classes;
    latlen := Endo.latticeLength;
    reps := Endo.latticeReps;
    Enum := Length( classes );
    Print( "There are ", Enum, " endomorphism classes.\n" );
    Print( "Calculating idempotent endomorphisms.\n" );
    idem := IdempotentImages( G, 4 );
    ilen := List( idem, i -> Length( i ) );
    ilen[1] := 1;
    ilen[latlen] := 1;
    Print( "# idempotents mapping to lattice class representatives:\n" );
    Print( ilen, "\n" );
    C := [ ];
    Cnum := 0;
    for i in [1..latlen] do
        R := reps[i];
        LI := idem[i];
        Inum := Length( LI );
        for j in [1..Inum] do
            h := GroupHomomorphismByImages( G, R, genG, LI[j] );
            kerh := Kernel( h );
            # (e;t,h:G->R) isomorphic to (e;h,t:G->R), so take k>=j
            for k in [j..Inum] do
                t := GroupHomomorphismByImages( G, R, genG, LI[k] );
                kert := Kernel( t );
                kerth := CommutatorSubgroup( kert, kerh );
                if ( Size( kerth ) = 1 ) then
                    if ( XModPrintLevel > 2 ) then
                        Print( "trivial kerth at i,j,k = ", [i,j,k], "\n" );
                    fi;
                    CC := Cat1( G, t, h );
                    g := Cnum; 
                    AICG := false;
                    while ( not AICG and ( g > 0 ) ) do
                        # this is expensive?:
                        AICG := AreIsomorphicCat1s( C[g], CC );
                        g := g-1;
                    od;
                    if ( AICG = false ) then
                        Add( C, CC );
                        Cnum := Cnum + 1;
                    fi; 
                fi;
            od;
        od;
    od; 

    # when G is abelian, the zero map gives a cat1-structure
    if IsAbelian( G ) then
        Q := reps[1];
        t := ZeroMorphism( G, Q );
        CC := Cat1( G, t, t );
        C := Concatenation( [ CC ], C );
        Cnum := Cnum + 1;
    fi;
    # the identity is the only surjective cat1-structure up to isomorphism
    t := GroupHomomorphismByImages( G, G, genG, genG );
    CC := Cat1( G, t, t );    
    Add( C, CC );
    Cnum := Cnum + 1;
    for i in [1..Cnum] do
        CC := C[i];
        Print( "Isomorphism class ", i, "\n" );
        if ( Size( CC.kernel ) < 101 ) then
            Print( ": kernel of tail = ", GroupId( CC.kernel ).names, "\n" );
        else
            Print( ": kernel has size = ", Size( CC.kernel ), "\n" );
        fi;
        if ( Size( CC.range ) < 101 ) then
            Print( ":    range group = ", GroupId( CC.range ).names, "\n" );
        else
            Print( ":  range has size = ", Size( CC.range ), "\n" );
        fi;
        bdy := CC.boundary;
        kerbdy := Kernel( bdy );
        if ( kerbdy <> CC.kernel ) then
            if ( Size( kerbdy ) < 101 ) then
                Print( ":  kernel of bdy = ", GroupId( kerbdy ).names, "\n" );
            else
                Print( ": kernel of bdy has size = ", Size( kerbdy ), "\n" );
            fi;
        fi;
    Print( "\n" );
    od;
    return C;
end;

##############################################################################
##
#F  BacktrackAllCat1sall               alternative method for all cat1-strures
##

BacktrackAllCat1s := function( G )

    local  el, genG, len, lenel, set, tup, lent, L, image, gens, phi, im, ok, 
           D, C, lenl, K, imim, images, t, h, prod, cond, kert, kerh, kerth,
           CC, g, i, j, k, lenc, AICG, imagelist, pos, pos1, Llist, lenim,
           homnum, homcount, catnum, A, elA, imageclass, classnum, class,
           found, K, genK, alpha, subclass, gen, K2, catcount, classlen,
           bdy, kerbdy;

    el := Elements( G );
    genG := G.generators;
    len := Length( genG );
    lenel := Length( el );
    A := AutomorphismGroup( G );
    elA := Elements( A );
    set := [ 1..lenel ];
    tup := AllCat1sJ( set, len, [ ], [ ], 1, G );
    lent := Length( tup );  
    Print( "     Number of homomorphisms, phi : G -> G = ", lent, "\n" );
    L := [ ];
    image := [ ];
    imagelist := [ ];
    imageclass := [ [ G ] ];
    Llist := [ [ ] ];
    classnum := 1;
    homnum := 0;
    homcount := 0; 
    for i in [1..lent] do
        gens := List( tup[i], x -> el[x] );
        phi := GroupHomomorphismByImages( G, G, genG, gens );
        ok := IsHomomorphism( phi );
        if ok then
            if ( phi * phi = phi ) then
                homnum := homnum + 1;
                im := Image( phi, G );
                class := 0; 
                found := false;
                while ( not ( found ) and ( class < classnum ) ) do
                    class := class + 1;
                    imagelist := imageclass[ class ];
                    pos := Position( imagelist, im );
                    found := not( pos = false );
                od;
                if ( pos = 1 ) then
                    Add( Llist[class], i ); 
                    homcount := homcount + 1;
                elif ( pos = false ) then
                    Add( imageclass, [ im ] );  
                    classnum := Length( imageclass );
                    Add( Llist, [ i ] );
                    K := im;
                    homcount := homcount + 1;
                    genK := K.generators;
                    subclass := [ K ];
                    for alpha in elA do
                        gen := List( genK, g -> g^alpha );
                        K2 := Subgroup( G, gen );
                        pos := Position( subclass, K2 );
                        if ( pos = false ) then
                            Add( subclass, K2 );
                        fi;
                    od;
                    imageclass[ classnum ] := subclass;
                fi; 
            fi;
        fi;
    od;   
    Print( "Number of homomorphisms with (phi^2 = phi) = ", homnum, "\n" ); 
    Print( "         Number of homomorphisms processed = ",homcount,"\n" );
    if ( XModPrintLevel > 2 ) then
        Print( "imageclass = \n", imageclass, "\n" );
        Print( "Llist =\n", Llist, "\n" );
    fi;
    G.imageclass := imageclass;
    catcount := 0 * [1..classnum];
    classlen := List( imageclass, L -> Length( L ) );
    C := [ ];
    lenc := Length( C );
    for k in [1..classnum] do
        if ( XModPrintLevel > 2 ) then
            Print( "Class ", k, " :-\n" );
        fi;
        K := imageclass[ k ][ 1 ];  
        if ( K = G ) then 
            t := GroupHomomorphismByImages( G, G, genG, genG );
            CC := Cat1( G, t, t );    
            Add( C, CC );
            catcount[k] := catcount[k] + 1;
            lenc := lenc + 1;
        else
            L := Llist[ k ]; 
            lenl := Length( L );
            for i in [1..lenl] do
                images := List( tup[L[i]], x -> el[x] ); 
                h := GroupHomomorphismByImages( G, G, genG, images );
                kerh := Kernel( h );
                for j in [i..lenl] do
                    if ( XModPrintLevel > 2 ) then
                        Print( "(", i, ",", j, ") " );
                    fi;
                    imim := List( tup[L[j]],x -> el[x] );
                    t := GroupHomomorphismByImages( G, G, genG, imim );
                    kert := Kernel( t );
                    kerth := CommutatorSubgroup( kert, kerh );
                    if ( Size( kerth ) = 1 ) then
                        CC := Cat1( G, t, h );
                        catcount[k] := catcount[k] + 1;
                        lenc := Length( C );
                        if ( lenc = 0 ) then
                            Add( C, CC );
                            lenc := lenc + 1;
                        else
                            g := lenc; 
                            AICG := false;
                            while ( not AICG and ( g > 0 ) ) do
                                AICG := AreIsomorphicCat1s( C[g], CC );
                                g := g - 1;
                            od;
                            if ( AICG = false ) then
                                Add( C, CC );
                                lenc := lenc + 1;
                                if ( t <> h ) then
                                     Print( t, "\n", h, "\n\n" );
                                fi;
                            fi;
                        fi;  
                    fi;  
                od;
                if ( XModPrintLevel > 2 ) then
                    Print( "\n" );
                fi;
            od; 
        fi;
    od;  
    catnum := 0;
    for i in [1..classnum] do
        catnum := catnum + catcount[i] * classlen[i];
    od;
    Print( "                       Number of cat-1-groups = ", catnum, "\n");
    Print( "Number of isomorphism classes of cat-1-groups = ", lenc, "\n\n" );
    for i in [1..lenc] do
        CC := C[i];
        Print( "Isomorphism class ", i, "\n" );
        Print( ":    range group = ", GroupId( CC.range ).names, "\n" );
        Print( ": kernel of tail = ", GroupId( CC.kernel ).names, "\n" );
        bdy := CC.boundary;
        kerbdy := Kernel( bdy );
        if ( kerbdy <> CC.kernel ) then
            Print( ":  kernel of bdy = ", GroupId( kerbdy ).names, "\n" );
        fi;
    od;
    return C;
end;

##############################################################################
##
#F  SubCat1                sub-cat1-group determined by subgroup of the source
##

SubCat1 := function( C, Dsrc )

    local  Csrc, Crng, tail, head, Drng, Drng1, embed, eDrng, genDsrc, 
           imh, imeh, Dhead, imt, imet, Dtail, D, i, ok;
          
    if not ( IsCat1( C ) ) then
        Error( "Argument is not a cat1-group" );
    fi;   
    Csrc := C.source;
    if not ( IsSubgroup( Csrc, Dsrc ) ) then
        Print("Dsrc is not a subgroup of Csrc\n");
        return false;
    fi;   
    Crng := C.range;
    tail := C.tail;
    head := C.head;
    embed := C.embedRange;
    Drng := Image( tail, Dsrc );
    Drng1 := Image( head, Dsrc );
    if not ( Drng = Drng1 ) then
        Print( "t(Dsrc) <> h(Dsrc)\n");
        return false;
    fi;
    eDrng := Image( embed, Drng );
    if not ( IsSubgroup( Dsrc, eDrng ) ) then
        if ( XModPrintLevel > 2 ) then
            Print( "t(Dsrc) does not embed in Dsrc\n" );
        fi;
        return false;
    fi;   
    genDsrc := Dsrc.generators;
    imh := List( genDsrc, i -> Image( head, i ) );
    imeh := List( imh, i -> Image( embed, i ) );
    Dhead := GroupHomomorphismByImages( Dsrc, Dsrc, genDsrc, imeh );
    imt := List( genDsrc, i -> Image( tail, i ) );
    imet := List( imt, i -> Image( embed, i ) );
    Dtail := GroupHomomorphismByImages( Dsrc, Dsrc, genDsrc, imet );
    D := Cat1( Dsrc, Dhead, Dtail );
    ok := IsCat1( D );
    if not ok then
        return false;
    else
        D.parent := C;
        D.name := Concatenation( "[Sub", C.name, "]" );
        return D;
    fi;
end;   

##############################################################################
##
#F  IdentitySubCat1                    identity sub-cat1-group of a cat1-group
##

IdentitySubCat1 := function( C )

    local  Csrc, idelsrc, Isub;

    if not ( IsCat1 ( C ) ) then
	Error( "Input parameter must be a cat1-group" );
    fi;
    Csrc := C.source;
    idelsrc := Elements( Csrc )[1];
    Isub := SubCat1( C, Subgroup( Csrc, [idelsrc] ) );
    Isub.name := Concatenation( "[Id", C.name, "]" );
    return Isub;
end; 

##############################################################################
##
#F  RepSubCat1s       representative SubCat1 for a conjugacy class of C.source
##

RepSubCat1s := function( C )

    local  Csrc, L, class, len, D, i, rep, sub;
    
    if not ( IsCat1( C ) ) then
        Error("Argument is not a cat1-group");
    fi;
    Csrc := C.source;
    L := Lattice( Csrc );
    class := L.classes;
    len := Length( class );
    D := [ ];
    for i in [1..len] do 
        rep := class[i].representative;
        sub := SubCat1( C, rep );
        if IsRec( sub ) then
            Add( D, sub );
        fi;
    od;
    return D;
end;  

##############################################################################
##
#F  NormalSubCat1s                                  all normal sub-cat1-groups
##

NormalSubCat1s := function( C )

    local  Csrc, normal, normsub, norm, sub;
    
    if not ( IsCat1( C ) ) then
        Error( "Argument is not a cat1-group" );
    fi;
    Csrc := C.source;
    normal := NormalSubgroups( Csrc );
    normsub := [ ];
    for norm in normal do
        sub := SubCat1( C, norm );
        if IsRec( sub ) then
            Add( normsub, sub );
        fi;
    od;
    return normsub;
end;

##############################################################################
##
#F  Cat1Select                           get xmod or cat1-group from data file
##

Cat1Select := function( arg )

    local  type, size, gpnum, num, norm, nargs, usage, maxsize,
           start, iso, count, pos, pos2, names,
           i, j, ncat1, G, genG, M, L, genR, R, t, kert, h, C, X;

    maxsize:= Cat1ListMaxSize;
    nargs := Length( arg );
    usage := 
       "\nUsage:  Cat1Select( size, gpnum, num )\n\n";
    if ( ( nargs = 0 ) or not 
            ( ( arg[1] > 0 ) and ( arg[1] <= maxsize ) ) ) then
        Print( usage );
        return 0;
    fi;

    # find starting positions of iso classes of groups of size <= maxsize
    iso := NumbersOfIsomorphismClasses;
    count := 1;
    start := [ 1 ];
    for j in iso do
        count := count + j;
        Add( start, count );
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "  iso = ", iso, "\n  start = ", start, "\n" );
    fi;

    size := arg[1];
    if ( ( size < 1 ) or ( size > maxsize ) ) then
        Error( "only groups of order up to ", maxsize, " in Cat1List");
        return false;
    fi;
    pos := start[ size ];
    if ( size < maxsize ) then
        pos2 := start[ size + 1 ] - 1;
    else
        pos2 := Length( Cat1List );
    fi;
    names := List( [ pos..pos2], n -> Cat1List[n][4] );
    if ( ( nargs = 1 ) or not ( arg[2] > 0 ) ) then
        Print( usage );
        return names;
    fi;

    gpnum := arg[2];
    if ( gpnum > iso[size] ) then
        Print( "# classes of size ", size, " is ", iso[size], "\n" );
        Print( "second parameter ", gpnum ); 
        Print( " > # classes of groups with this size" );
        Print( usage );
        return names;
    fi;

    j := pos + gpnum - 1;
    M := Cat1List[j];
    G  := Group(M[3], ( ));
    G.name := M[4];
    ncat1 := Length( M[5] ) + 1;
    genG := G.generators;
    if ( ( nargs = 2 ) or not
         ( ( arg[3] >= 1 ) and ( arg[3] <= ncat1 ) ) ) then
        Print( "\nThere are ", ncat1, " cat1-structures for the group ");
        Print( G, ".\n" ); 
        Print( "[ [range gens], source & range names," ); 
        Print( " [tail.genimages], [head.genimages] ]" );
        Print( "  :- \n" );
        Print( "[ ", genG, ",  tail = head = identity mapping ]\n" ); 
        for i in [2..ncat1] do
            Print( M[5][i-1], "\n" );
        od;
        Print( usage, "Group has generators " );
        return( genG );
    fi;

    num := arg[3];
    if ( num = 1 ) then
        L := [ genG, G.name, G.name, genG, genG ];
    else
        L := M[5][num-1];
    fi;
    genR := L[1];
    R := Subgroup( G, genR );
    R.name := L[3];
    t := GroupHomomorphismByImages( G, R, genG, L[4] );
    h := GroupHomomorphismByImages( G, R, genG, L[5] );
    kert := Kernel( t );
    kert.name := L[2];
    C := Cat1( G, t, h );
    
    if ( nargs = 3 ) then
        return C;
    fi;
    return names;
end;

##############################################################################
##
#F  ReverseCat1                                  C -> C~,  with  t,h  reversed
##

ReverseCat1 := function( C )

    local  G, R, t, h, e, kert, kertgen, imbdy, bdy, incK, D, ok;

    if not IsCat1( C ) then
        Error( "parameter must be a cat1-group" );
    fi;
    G := C.source;
    R := C.range;
    t := C.head;
    h := C.tail;
    e := C.embedRange;
    kert := Kernel( t );
    kertgen := kert.generators;
    imbdy := List( kertgen, x -> Image( h, x ) );
    bdy := GroupHomomorphismByImages( kert, R, kertgen, imbdy );
    incK := InclusionMorphism( kert, G );

    D := rec ( );
    D.source := G;
    D.range := R;
    D.tail := t;
    D.head := h;
    D.embedRange := e;
    D.kernel := kert;
    D.boundary := bdy;
    D.embedKernel := incK;
    ok := IsCat1( D );
    if ok then
        D.isDomain := true;
        D.operations := Cat1Ops;
        D.name := Cat1Name( D );
    fi;
    if ( G = R ) then
        D.isIdentity := true;
    fi;
    return D;
end;

##############################################################################
##
#F  ReverseIsomorphismCat1                      hom : C -> C~, g -> g~, r -> r
##

ReverseIsomorphismCat1 := function( C )

    local  G, genG, R, idR, t, h, e, et, eh, im, phi, D, mor;

    if not IsCat1( C ) then
        Error( "parameter must be a cat1-group" );
    fi;
    G := C.source;
    t := C.tail;
    h := C.head;
    e := C.embedRange;
    et := CompositionMapping( e, t );
    eh := CompositionMapping( e, h );
    genG := G.generators;
    im := List( genG, g -> Image( eh, g ) * g^-1 * Image( et, g ) );
    phi := GroupHomomorphismByImages( G, G, genG, im );
    if not IsAutomorphism( phi ) then
        Error( "phi fails to be an automorphism" );
    fi;
    R := C.range;
    idR := InclusionMorphism( R, R );
    idR.genimages := R.generators;
    D := ReverseCat1( C );
    mor := Cat1Morphism( C, D, [ phi, idR ] );
    return mor;
end;

##############################################################################
##
#F  Cat1Ops.InclusionMorphism                    identity morphism for an Cat1
##

Cat1Ops.InclusionMorphism := function( C, D )

    local  idX, idsrc, idrng, ok;

    idsrc := InclusionMorphism( C.source, D.source );
    idrng := InclusionMorphism( C.range, D.range );
    idX := Cat1Morphism( C, D, [idsrc, idrng] );
    ok := IsCat1Morphism( idX );
    return idX;
end;

##############################################################################
##
#F  Cat1MorphismOps.InverseMorphism                inverse of an Cat1 morphism
##

Cat1MorphismOps.InverseMorphism := function( mor )

    local  invmor, invsrc, invrng, ok;

    if not IsCat1Morphism( mor ) then
        Error( "parameter must be a morphism of cat1-group" );
    fi;
    if not IsMonomorphism( mor ) then
        Error( "parameter must be a monomorphism of cat1-groups" );
    fi;
    if not (     IsIsomorphism( mor.sourceHom )
             and IsIsomorphism( mor.rangeHom ) ) then
    Error( "source and range homomorphisms not isomorphisms" );
    fi;
    invsrc := InverseMapping( mor.sourceHom );
    invrng := InverseMapping( mor.rangeHom );
    invmor := Cat1Morphism( mor.range, mor.source, [invsrc,invrng] );
    ok := IsCat1Morphism( invmor );
    return invmor;
end;

##############################################################################
##
#F  Cat1Ops.WhiteheadPermGroup( X )                  group of regular sections
##

Cat1Ops.WhiteheadPermGroup := function( C )

    local  X, S, D, W;

    X := XModCat1( C );
    S := RegularSections( C );
    D := RegularDerivations( X );
    W := WhiteheadPermGroup( X );
    return W;
end;

##############################################################################
##
#F  Cat1Ops.Actor                         constructs actor of C via xmod actor
##

Cat1Ops.Actor := function( C )

    local  X, Y, name, D;

    X := XModCat1( C );
    Y := Actor( X );
    D := Cat1XMod( Y );
    name := Cat1Name( C );
    D.name := Concatenation( "Actor", name );
    return D;
end;

##############################################################################
##
#F  Cat1Ops.InnerActor                  group of inner actions of a cat1-group
##

Cat1Ops.InnerActor := function( C )

    local X, Y, D;

    X := XModCat1( C );
    Y := InnerActor( X );
    D := Cat1XMod( Y );
    return Y;
end;

##  end of file  cat1.g
