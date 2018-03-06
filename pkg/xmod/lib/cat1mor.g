##############################################################################
##
##  cat1mor.g                for GAP 3.4                      version  9/ 1/97
##
##############################################################################
##
#A  cat1mor.g                GAP library                         Chris Wensley
#A                                                                   Murat Alp
#Y  Copyright
##
##  This file contains functions that manipulate cat-1 groups
##  and associated constructions.
##
#H  $Log: cat1mor.g,v $
#H  Revision 1.1  1997/03/27 13:35:45  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

##############################################################################
##
#F  IsCat1Morphism                     check whether a morphism of cat1-groups
##

IsCat1Morphism := function( phi )

    local  phisrc, phirng, C, D, Csrc, Crng, Dsrc, Drng,
           c2, d1, b1, gensrc, genrng, phisrc, phirng;

    if not IsRec( phi ) then
        return false;
    fi;
    if not ( IsBound( phi.domain ) and ( phi.domain = Mappings ) ) then
        return false;
    fi;
    if IsBound( phi.isCat1Morphism ) then
        return phi.isCat1Morphism;
    fi;

    phi.isCat1Morphism := false;

    if (   not ( IsBound( phi.source ) and IsCat1( phi.source ) )
        or not ( IsBound( phi.range ) and IsCat1( phi.range ) ) ) then
      Print( "Error: phi.source & phi.range must exist as cat1-groups! \n" );
      Print( phi.source, " -> ", phi.range, "\n" );
      return false;
    fi;

    C := phi.source;
    D := phi.range;
    Csrc := C.source;
    Crng := C.range;
    Dsrc := D.source;
    Drng := D.range;
    if not (     IsBound( phi.rangeHom )
             and IsHomomorphism( phi.rangeHom )
             and ( phi.rangeHom.source = Crng )
             and ( phi.rangeHom.range = Drng )  ) then
        if ( XModPrintLevel > 2 ) then
            Print( "phi.rangeHom not a hom or incorrect source/range\n" );
        fi;
        return false;
    fi;

    # now check that the homomorphisms commute
    phisrc := phi.sourceHom;
    phirng := phi.rangeHom;
    gensrc := Csrc.generators;
    genrng := Crng.generators;
    if ( XModPrintLevel > 2 ) then
        Print( "checking that the diagram commutes : -\n" );
        Print( "D.tail(phisrc(x)) = phirng(C.tail(x)) \n" );
        Print( "D.head(phisrc(x)) = phirng(C.head(x)) \n" );
    fi;
    for c2 in gensrc do
        d1 := ( c2 ^ phisrc ) ^ D.tail;
        b1 := ( c2 ^ C.tail ) ^ phirng;
        if not (d1 = b1) then
            if ( XModPrintLevel > 2 ) then
                Print( "tails do not commute with phi! \n" );
                Print( "when c2= ", c2, " D.tail(phisrc(c2))= ", d1, "\n" );
                Print( "              and phirng(C.tail(c2))= ", b1, "\n" );
            fi;
            return false;
        fi;
    od;
    for c2 in gensrc do
        d1 := ( c2 ^ phisrc ) ^ D.head;
        b1 := ( c2 ^ C.head ) ^ phirng;
        if not (d1 = b1) then
            if ( XModPrintLevel > 2 ) then
                Print( "tails do not commute with phi! \n" );
                Print( "when c2= ", c2, " D.head(phisrc(c2))= ", d1, "\n" );
                Print( "              and phirng(C.head(c2))= ", b1, "\n" );
            fi;
            return false;
        fi;
    od;

    phi.isCat1Morphism := true;
    return true;
end;

##############################################################################
##
#F  Cat1MorphismName                concatenates the names of source and range
##

Cat1MorphismName := function( mor )

    local nsrc, nrng, name;

    if not IsCat1Morphism( mor ) then
        Error( "Input parameter must be a Cat1 morphism\n" );
    fi;
    if IsCat1Morphism( mor ) then
        if IsBound( mor.source.name ) then
            nsrc := mor.source.name;
        else
            nsrc := "[?]";
        fi;
        if IsBound(mor.range.name) then
            nrng := mor.range.name;
        else
            nrng := "[?]";
        fi;
        name := Concatenation("<", nsrc, "-->", nrng, ">");
        mor.name := name;
    fi;
    return name;
end; 

##############################################################################
##
#V  Cat1MorphismOps                      cat1-group morphism record operations
##

Cat1MorphismOps := OperationsRecord( "Cat1MorphismOps", MappingOps );


##############################################################################
##
#F  Cat1MorphismOps.\=             define equality of morphisms of cat1-groups
##

Cat1MorphismOps.\= := function( mor1, mor2 )

    local  isEql;

    isEql :=     ( mor1.source = mor2.source )
             and ( mor1.range = mor2.range )
             and ( mor1.sourceHom = mor2.sourceHom )
             and ( mor1.rangeHom = mor2.rangeHom );
    return isEql;
end;

##############################################################################
##
#F  Cat1MorphismOps.IsMonomorphism
##

Cat1MorphismOps.IsMonomorphism := function( mor )

    local  ok;

    if not IsCat1Morphism( mor ) then
        return false;
    fi;
    if IsBound( mor.isMonomorphism ) then
        return mor.isMonomorphism;
    fi;
    ok := (IsMonomorphism( mor.sourceHom ) and IsMonomorphism( mor.rangeHom ));
    mor.isMonomorphism := ok;
    return ok;
end;

##############################################################################
##
#F  Cat1MorphismOps.IsEpimorphism
##

Cat1MorphismOps.IsEpimorphism := function( mor )

    local  ok;

    if not IsCat1Morphism( mor ) then
        return false;
    fi;
    if IsBound( mor.isEpimorphism ) then
        return mor.isEpimorphism;
    fi;
    ok := ( IsEpimorphism( mor.sourceHom ) and IsEpimorphism( mor.rangeHom ) );
    mor.isEpimorphism := ok;
    return ok;
end;

##############################################################################
##
#F  Cat1MorphismOps.IsIsomorphism
##

Cat1MorphismOps.IsIsomorphism := function( mor )

    local  ok;

    if not IsCat1Morphism( mor ) then
        return false;
    fi;
    if IsBound( mor.isIsomorphism ) then
        return mor.isIsomorphism;
    fi;
    ok := ( IsMonomorphism( mor ) and IsEpimorphism( mor ) );
    mor.isIsomorphism := ok;
    return ok;
end;

##############################################################################
##
#F  Cat1MorphismOps.IsEndomorphism
##

Cat1MorphismOps.IsEndomorphism := function( mor )

    local  ok;

    if not IsCat1Morphism( mor ) then
        return false;
    fi;
    if IsBound( mor.isEndomorphism ) then
        return mor.isEndomorphism;
    fi;
    ok := (IsEndomorphism( mor.sourceHom ) and IsEndomorphism( mor.rangeHom ));
    mor.isEndomorphism := ok;
    return ok;
end;

##############################################################################
##
#F  Cat1MorphismOps.IsAutomorphism
##

Cat1MorphismOps.IsAutomorphism := function( mor )

    local  ok;

    if not IsCat1Morphism( mor ) then
        return false;
    fi;
    if IsBound( mor.isAutomorphism ) then
        return mor.isAutomorphism;
    fi;
    ok := ( IsEndomorphism( mor ) and IsIsomorphism( mor ) );
    mor.isAutomorphism := ok;
    return ok;
end;

##############################################################################
##
#F  Cat1MorphismOps.Print      special print commands for car1-group morphisms
##

Cat1MorphismOps.Print := function( mor )

    local name;

    if IsBound( mor.name ) then
        name := mor.name;
    else
        name := "<[?]-->[?]>";
    fi;
    Print( "Morphism of cat1-groups ", name );
end;

##############################################################################
##
#F  Cat1Morphism                    cat1 morphism from a pair of homomorphisms
##

Cat1Morphism := function( arg )

    local  C, D, L, Csrc, Dsrc, Crng, Drng, Csrcgen, f, g, h, t, tt, hh, e,
           nargs, usage, ee, fe, ttf, gt, gtim, hhf, gh, ghim, mor, ismor;

    nargs := Length( arg );
    usage := "\nUsage: Cat1Morphism( C, D, [gamma,rho], [true|false] );\n";
    if not ( ( nargs > 2 ) and ( nargs < 5 ) ) then
        Print( usage );
        return false;
    fi;
    C := arg[1];
    D := arg[2];
    if (    not ( IsBound( C.isCat1 ) and C.isCat1 )
         or not ( IsBound( D.isCat1 ) and D.isCat1 ) ) then
        Print( usage );
        Error( "first two parameters must be cat1-groups" );
    fi;
    Csrc := C.source;
    Dsrc := D.source;
    Crng := C.range;
    Drng := D.range;
    Csrcgen := Csrc.generators;
    L := arg[3];
    if ( not ( IsList( L ) and ( Length( L ) = 2 ) ) ) then
        Error( "the third parameter must be a pair of homomorphisms" );
    fi;
    f := L[1];
    g := L[2];
    t := C.tail;
    h := C.head;
    tt := D.tail;
    hh := D.head;
    e := C.embedRange;
    ee := D.embedRange;
    fe := CompositionMapping( f, e );
    ttf := CompositionMapping( tt, f );
    gtim := List( t.genimages, x -> Image( g, x ) );
    gt := GroupHomomorphismByImages( Csrc, Drng, Csrcgen, gtim );
    if ( ( XModPrintLevel > 1 ) and ( ttf = gt ) ) then
        Print( "condition  Dtail.f = Dhead.f.CembedRange.Ctail  satisfied.\n");
    fi;
    hhf := CompositionMapping( hh, f );
    ghim := List( h.genimages, x -> Image( g, x ) );
    gh := GroupHomomorphismByImages( Csrc, Drng, Csrcgen, ghim );
    if ( ( XModPrintLevel > 1 ) and ( hhf = gh ) ) then
        Print( "condition  Dhead.f = Dhead.f.CembedRange.Chead  satisfied.\n");
    fi;
    ismor := ( ( nargs = 4 ) and IsBool( arg[4] ) and arg[4] );
    mor := rec ( );
    mor.source := C;
    mor.range := D;
    mor.sourceHom := f;
    mor.rangeHom := g;
    mor.domain := Mappings;
    if not ismor then
        ismor := IsCat1Morphism( mor );
    fi;
    if ismor then
        mor.operations := Cat1MorphismOps;
        mor.name := Cat1MorphismName( mor );
    fi;
    return mor;
end;

##############################################################################
##
#F  Cat1MorphismPrint               print details of a morphism of cat1-groups
##

Cat1MorphismPrint := function( mor )

    local  C, Csrc, Crng, D, Dsrc, Drng, ok, name, morsrc, morrng;

    if not IsBound( mor.isCat1Morphism ) then
        ok := IsCat1Morphism( mor );
    fi;
    if not mor.isCat1Morphism then
        Print("Cat1MorphismPrint can only print a morphism of cat1-groups\n");
        return false;
    fi;
    if not IsBound( mor.name ) then
        Cat1MorphismName( mor );
    fi;
    name := mor.name;

    C := mor.source;
    Csrc := C.source;
    Crng := C.range;
    D := mor.range;
    morsrc := mor.sourceHom;
    morrng := mor.rangeHom;
    if not ( IsBound( morsrc.genimages ) and IsBound( morrng.genimages ) ) then
        Print( "Missing genimages in mor.source(range)Hom\n" );
    fi;
    Print( "Morphism of cat1-groups := \n" );
    Print( ": Source = ", C, "\n" );
    Print( ":  Range = ", D, "\n" );
    Print( ": Source homomorphism maps source generators to:\n" );
    Print( "  ", morsrc.genimages, "\n" );
    Print( ": Range homomorphism maps range generators to:\n" );
    Print( "  ", morrng.genimages, "\n" );
end;

##############################################################################
##
#F  Cat1MorphismOps.CompositeMorphism       composite of 2 Cat1 morphisms
##

Cat1MorphismOps.CompositeMorphism := function( phi, psi )

    local C, D, Csrc, Crng, Dsrc, Drng, genCsrc, genCrng,
          phisrc, phirng, psisrc, psirng, homgensrc, homsrc,
          homgenrng, homrng, hom, x;

    C := phi.source;
    D := psi.range;
    Csrc := C.source;
    Crng := C.range;
    Dsrc := D.source;
    Drng := D.range;
    genCsrc := Csrc.generators;
    genCrng := Crng.generators;
    phisrc := phi.sourceHom;
    phirng := phi.rangeHom;
    psisrc := psi.sourceHom;
    psirng := psi.rangeHom;
    homgensrc := List(genCsrc, x -> Image(psisrc, Image(phisrc,x)));
    homsrc := GroupHomomorphismByImages(Csrc,Dsrc,genCsrc,homgensrc);
    homgenrng := List(genCrng, x -> Image(psirng, Image(phirng,x)));
    homrng := GroupHomomorphismByImages(Crng,Drng,genCrng,homgenrng);
    hom := Cat1Morphism(C,D,[homsrc,homrng]);
    return hom;
end;

##############################################################################
##
#F  Cat1MorphismSourceHomomorphism    morphism of cat1-groups using f : G1->G2
##

Cat1MorphismSourceHomomorphism := function( C, D, f )

    local  Csrc, Dsrc, fe, res, mor;

    if (    not ( IsBound( C.isCat1 ) and C.isCat1 )
         or not ( IsBound( D.isCat1 ) and D.isCat1 ) ) then
        Error( "first two parameters must be cat1-groups" );
    fi;
    Csrc := C.source;
    Dsrc := D.source;
    if (    not IsHomomorphism( f )
         or not ( f.source = Csrc )
         or not ( f.range = Dsrc ) ) then
        Error( "the third parameter must be a hom from Csrc to Dsrc" );
    fi;

    fe := CompositionMapping( f, C.embedRange );
    res := CompositionMapping( D.head, fe );
    res.genimages := List( C.range.generators, x -> Image( res, x ) );
    mor := Cat1Morphism( C, D, [ f, res ] );
    return mor;
end;

#############################################################################
##
#F  Cat1MorphismRangeHomomorphism
##

Cat1MorphismRangeHomomorphism := function( C, D, phi )

    local  Csrc, Crng, Dsrc, Drng, psi, Csrcgen, Dsrcim,
           Cmap, Dmap, g, g1, gm, g2, mor;

    if not ( IsCat1( C ) and IsCat1( D ) ) then
        Error( "first 2 parameters must be cat1-groups" );
    fi;
    if not (     IsHomomorphism( phi ) 
             and ( phi.source = C.range )
             and ( phi.range = D.range ) ) then
        Error( "invalid data" );
    fi;

    Csrc := C.source;
    Crng := C.range;
    Dsrc := D.source;
    Drng := D.range;
    Cmap := Csrc.mapping;
    Dmap := Dsrc.mapping;
    Csrcgen := Csrc.generators;
    Dsrcim := [ ];
    for g in Csrcgen do
        g1 := Image( phi, g.element[1] );
        gm := Image( Dmap, g1 );
        g2 := Image( phi, g.element[2] );
        Add( Dsrcim, SemidirectProductElement( g1, gm, g2 ) );
    od;
    psi := GroupHomomorphismByImages( Csrc, Dsrc, Csrcgen, Dsrcim );
    if not IsHomomorphism( psi ) then
        Error( "psi is not a homomorphism!" );
    fi;
    mor := Cat1Morphism( C, D, [ psi, phi ] );
    return mor;
end;

###############################################################################
##
#F  InclusionMorphismSubCat1
## 

InclusionMorphismSubCat1 := function( C, D )

    local Csrc, Crng, Dsrc, Drng, genCsrc, genCrng, 
           genDsrc, genDrng, morsrc, morrng, phi, ok;

    if ( not IsCat1( C ) ) or ( not IsCat1( D ) ) then
        Error( "Both arguments must be cat1 groups " );
    fi;
    Csrc := C.source;
    Crng := C.range;
    Dsrc := D.source;
    Drng := D.range;
    genCsrc := Csrc.generators;
    genCrng := Crng.generators;
    genDsrc := Dsrc.generators;
    genDrng := Drng.generators;
    morsrc := GroupHomomorphismByImages(Csrc,Dsrc,genCsrc,genCsrc);
    morrng := GroupHomomorphismByImages(Crng,Drng,genCrng,genCrng);
    phi := rec ( );
    phi.source := C;
    phi.range := D;
    phi.sourceHom := morsrc;
    phi.rangeHom := morrng;
    phi.domain := Mappings;
    ok := IsCat1Morphism( phi );
    if ok then
        phi.operations := Cat1MorphismOps;
        phi.name := Cat1MorphismName(phi);
    fi;
    return phi;
end;

############################################################################
##
#F  XModMorphismCat1Morphism
##

XModMorphismCat1Morphism := function(  mor )

    local  C, D, X, Y, Xsrc, Xrng, Ysrc, Yrng, genXsrc, genXrng,
           x, morimage, imgen, im, images, phis, phir, phi;

    if not ( IsCat1Morphism( mor ) ) then
        Error( "Argument  must be cat1-group morphism" );
        return false;
    fi;
    if IsBound( mor.xmodMorphism ) then
         return mor.xmodMorphism;
    fi;
    C := mor.source;
    D := mor.range;
    X := XModCat1(C);
    Y := XModCat1(D);
    if not ( IsXMod( X ) and IsXMod( Y ) ) then
        Error( "Arguments must be crossed module" );
        return false;
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    Ysrc := Y.source;
    Yrng := Y.range;
    genXsrc := Xsrc.generators;
    genXrng := Xrng.generators;
    morimage := mor.rangeHom.genimages;
    imgen := List( genXsrc, x -> SemidirectProductElement( (), (), x ) );
    im := List( imgen, x -> Image( mor.sourceHom, x ) );
    images := List( im, x -> x.element[2] );
    phis := GroupHomomorphismByImages( Xsrc, Ysrc, genXsrc, images);
    phir := GroupHomomorphismByImages( Xrng, Yrng, genXrng, morimage);
    phi := XModMorphism( X, Y, [phis,phir] );
    mor.xmodMorphism := phi;
    phi.cat1Morphism := mor;
    return phi;
end;

###########################################################################
##
#F Cat1MorphismXModMorphism
##

Cat1MorphismXModMorphism := function( mor )

    local  X, Y, C, D, Yact, Csrc, Crng, Dsrc, Drng, Ds2p, 
           genCsrc, genCrng, genXsrc, genXrng, imager, imgen, m, images,
           phis, phir, phi, morsrc, morrng, e, ek, g, ig, eg, pg, esrc, erng;

    if not ( IsXModMorphism( mor ) ) then
        Error( "Argument  must be crossed module morphism" );
    fi;
    if IsBound( mor.cat1Morphism ) then
        return mor.cat1Morphism;
    fi;
    X := mor.source;
    Y := mor.range;
    C := Cat1XMod(X);
    D := Cat1XMod(Y);
    if not ( IsCat1( C ) and IsCat1( D ) ) then
        Error( "Arguments must be cat1-groups" );
    fi;
    morsrc := mor.sourceHom;
    morrng := mor.rangeHom;
    e := D.embedRange;
    ek := D.embedKernel;
    Yact := Y.action;
    Csrc := C.source;
    Crng := C.range;
    Dsrc := D.source;
    Drng := D.range;
    if not ( IsPermGroup( Dsrc ) and IsBound( Dsrc.semidirect ) ) then
        Error( "<Dsrc> must be a premutation rep of a semidirect product" );
    fi;
    genCsrc := Csrc.generators;
    genCrng := Crng.generators;
    genXsrc := X.source.generators;
    genXrng := X.range.generators;
    imager := mor.rangeHom.genimages;
    phir := GroupHomomorphismByImages( Crng, Drng, genCrng, imager ); 
    if not IsHomomorphism( phir ) then
        Error( "<phir> not a homomorphism" );
    fi;
    images := [ ];
    erng := Dsrc.semidirect.embeddings[1];
    esrc := Dsrc.semidirect.embeddings[2];
    Ds2p := Dsrc.semidirect.semidirectPair.s2p;
    for g in genXrng do
        ig := Image( morrng, g );
        eg := Image( erng, ig );
        pg := Image( Ds2p, eg );
        Add( images, pg );
    od;
    for g in genXsrc do
        ig := Image( morsrc, g );
        eg := Image( esrc, ig );
        pg := Image( Ds2p, eg );
        Add( images, pg );
    od;
    phis := GroupHomomorphismByImages( Csrc, Dsrc, genCsrc, images );
    if not IsHomomorphism( phis ) then
        Error( "<phis> nit a homomorphism" );
    fi;
    phi := Cat1Morphism( C, D, [phis,phir] );
    mor.cat1Morphism := phi;
    phi.xmodMorphism := mor;
    return phi;
end;

##############################################################################
##
#F  AreIsomorphicCat1s
##

AreIsomorphicCat1s := function( C1, C2 )

    local  t1, h1, t2, h2, G1, G2, R1, R2, e1, A, elA, iso, alpha, ok, ko,
           phi, comp1, comp2, comp3, comp4, same;

    if not ( IsCat1( C1 ) and IsCat1( C2 ) ) then
        Print( "Arguments must be cat1-groups\n" );
        return false;
    fi;         
    G1 := C1.source;
    G2 := C2.source;
    same := ( G1 = G2 );
    if not same then
        phi := IsomorphismGroups( G1, G2 );
    else
        phi := InclusionMorphism( G1, G1 );
    fi;
    if ( phi = false ) then
        return false;
    fi;
    R1 := C1.range;
    R2 := C2.range;
    if not ( Size( R1 ) = Size( R2 ) ) then
        return false;
    fi;
    t1 := C1.tail;
    h1 := C1.head;
    t2 := C2.tail;
    h2 := C2.head;
    e1 := C1.embedRange;
    ko := true;
    A := AutomorphismGroup( G1 );
    elA := Elements( A );
    for iso in elA do
        if not same then
            iso := iso * phi;
        fi;
        comp1 := iso * h2;
        comp3 := iso * t2;
        alpha := e1 * comp3;
        ok := IsIsomorphism( alpha );
        if ok then 
            comp2 := h1 * alpha;
            comp4 := t1 * alpha; 
            if ( ( comp1 = comp2 ) and ( comp3 = comp4 ) ) then
                return true;
            fi;
        fi;
    od;
    return false;
end;          

##############################################################################
##
#F  ReverseIsomorphismCat1                 hom : C -> C~, g -> g~, r -> r
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
#F  Cat1Ops.InclusionMorphism               identity morphism for a cat1-group
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

##  end of file  cat1mor.g
