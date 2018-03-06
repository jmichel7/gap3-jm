###############################################################################
##
##  xmodmor.g                   for GAP 3.4                    version 10/ 1/97
##
###############################################################################
##
#A  xmodmor.g                   GAP library                       Chris Wensley
#A                                                                    Murat Alp
#Y  Copyright
##
##  This file contains functions that manipulate crossed modules,
##  pre-crossed modules, group graphs and associated constructions.
##
#H  $Log: xmodmor.g,v $
#H  Revision 1.1  1997/03/27 13:35:56  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

###############################################################################
##
#F  IsXModMorphism                         check if crossed module homomorphism
##

IsXModMorphism := function( phi )

    local phisrc, phirng, X, Y, Xrng, Xsrc, Yrng, Ysrc,
          x2, x1, y2, z2, y1, z1, gensrc, genrng;

    if not IsRec( phi ) then
        return false;
    fi;
    if not ( IsBound( phi.domain ) and ( phi.domain = Mappings ) ) then
        return false;
    fi;
    if IsBound ( phi.isXModMorphism ) then
        return phi.isXModMorphism;
    fi;

    if (   ( not ( IsBound( phi.source ) or IsXMod( phi.source ) ) ) 
        or ( not ( IsBound( phi.range )  or IsXMod( phi.range ) ) ) ) then
      Print( "Error: phi.source & phi.range must exist as XModPerms! \n" );
      Print( phi.source, " -> ", phi.range, "\n" );
      return false;
    fi;

    X := phi.source;
    Y := phi.range;
    if not ( IsBound( X.source ) and IsBound( X.range ) and
             IsBound( Y.source ) and IsBound( Y.range ) ) then
        return false;
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    Ysrc := Y.source;
    Yrng := Y.range;

    if not (     IsBound( phi.rangeHom ) 
             and IsHomomorphism( phi.rangeHom )
             and ( phi.rangeHom.source = Xrng )
             and ( phi.rangeHom.range = Yrng )   ) then
        if ( XModPrintLevel > 2 ) then
            Print( "phi.rangeHom not a hom or incorrect source/range\n" );
            fi;
        return false;
    fi;
    if not (     IsBound( phi.sourceHom ) 
             and IsHomomorphism( phi.sourceHom )
             and ( phi.sourceHom.source = Xsrc )
             and ( phi.sourceHom.range = Ysrc )   ) then
        if ( XModPrintLevel > 2 ) then
            Print( "phi.sourceHom not a hom or incorrect source/range\n" );
            fi;
        return false;
    fi;

    # now check that the homomorphisms commute
    phisrc := phi.sourceHom;
    phirng := phi.rangeHom;
    gensrc := Xsrc.generators;
    genrng := Xrng.generators;
    if ( XModPrintLevel > 2 ) then
        Print( "Checking that the diagram commutes :- \n" );
        Print( "         Y.boundary(phisrc(x)) = phirng(X.boundary(x)) \n" );
    fi;
    for x2 in gensrc do
        y1 := ( x2 ^ phisrc ) ^ Y.boundary;
        z1 := ( x2 ^ X.boundary ) ^ phirng;
        if not (y1 = z1) then
            if ( XModPrintLevel > 2 ) then
                Print( "Square does not commute! \n" );
                Print( "when x2= ", x2, " Ybdy(phisrc(x2))= ", y1, "\n" );
                Print( "              and phirng(Xbdy(x2))= ", z1, "\n" );
            fi;
            return false;
        fi;
    od;

    # now check that the actions commute:
    if ( XModPrintLevel > 2 ) then
        Print( "Checking:  phisrc(x2^x1) = phisrc(x2)^(phirng(x1))\n" );
    fi;
    for x2 in gensrc do
        for x1 in genrng do
            y2 := (x2 ^ (x1 ^ X.action)) ^ phisrc;
            z2 := (x2 ^ phisrc) ^ ((x1 ^ phirng) ^ Y.action);
            if not (y2 = z2) then
                if ( XModPrintLevel > 1 ) then
                    Print( "Actions do not commute! \n" );
                    Print( "When  x2 = ", x2, "  and  x1 = ", x1, "\n" );
                    Print( "           phisrc(x2^x1) = ", y2, "\n" );
                    Print( "and phisrc(x2)^(phirng(x1) = ", z2, "\n" );
                fi;
                return false;
            fi;
        od;
    od;

    phi.isXModMorphism := true;
    return true;
end;

###############################################################################
##
#F  XModMorphismName                Select a name for a crossed module morphism
##

XModMorphismName := function( phi )

    local  X, Xname, Y, Yname, name;

    # any existing name is NOT automatically returned

    if not ( IsBound( phi.source ) and IsBound( phi.range ) ) then
        Error( " Input parameter must have source and range fields" );
        return "false";
    fi;
    X := phi.source;
    Y := phi.range;
    if not ( IsXMod( X ) and IsXMod( Y ) ) then
        Error( "the source and range of phi must be crossed modules" );
        return "false";
    fi;

    # phi.isXModMorphism need NOT be bound, nor be true

    Xname := XModName( X );
    Yname := XModName( Y );
    name := Concatenation( "<", Xname, " >-> ", Yname, ">" );
    phi.name := name;
    return name;
end;

#############################################################################
##
#V  XModMorphismOps                 crossed module morphism record operations
##

XModMorphismOps := OperationsRecord( "XModMorphismOps", MappingOps );

#############################################################################
##
#F  XModMorphismOps.\=            define equality of crossed module morphisms
##

XModMorphismOps.\= := function( mor1, mor2 )

    local isEql;

    isEql :=     ( mor1.source = mor2.source )
             and ( mor1.range = mor2.range )
             and ( mor1.sourceHom = mor2.sourceHom )
             and (mor1.rangeHom = mor2.rangeHom );
    return isEql;
end;

##############################################################################
##
#F  XModMorphismOps.IsMonomorphism
##

XModMorphismOps.IsMonomorphism := function( mor )

    local  ok;

    if not IsXModMorphism( mor ) then
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
#F  XModMorphismOps.IsEpimorphism
##

XModMorphismOps.IsEpimorphism := function( mor )

    local  ok;

    if not IsXModMorphism( mor ) then
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
#F  XModMorphismOps.IsIsomorphism
##

XModMorphismOps.IsIsomorphism := function( mor )

    local  ok;

    if not IsXModMorphism( mor ) then
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
#F  XModMorphismOps.IsEndomorphism
##

XModMorphismOps.IsEndomorphism := function( mor )

    local  ok;

    if not IsXModMorphism( mor ) then
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
#F  XModMorphismOps.IsAutomorphism
##

XModMorphismOps.IsAutomorphism := function( mor )

    local  ok;

    if not IsXModMorphism( mor ) then
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
#F  XModMorphismOps.Order         order of an automorphism of a crossed module
##

XModMorphismOps.Order := function( phi )

    local  X, phisrc, phirng, Xsrc, Xrng, Xsrcauto, Xrngauto, ord;

    if not ( phi.isXModMorphism ) then
        Error( "Parameter is not a morphism of crossed modules" );
        return false;
    fi;
    if not IsAutomorphism( phi ) then
       Error( "Parameter is not an automorphism" );
       return false;
    fi;
    if IsBound( phi.order ) then
        return phi.order;
    fi;

    X := phi.source;
    Xsrc := X.source;
    Xrng := X.range;
    phisrc := phi.sourceHom;
    phirng := phi.rangeHom;
    if not IsBound( Xsrc.autoGroup ) then
        Xsrcauto := AutomorphismGroup( Xsrc );
        Xsrcauto.isAutomorphismGroup := true;
    else
        Xsrcauto := Xsrc.autoGroup;
    fi;
    if not IsBound( Xrng.autoGroup ) then
        Xrngauto := AutomorphismGroup( Xrng );
        Xrngauto.isAutomorphismGroup := true;
    else
        Xrngauto := Xrng.autoGroup;
    fi;    
    ord := Lcm( Integers, Order( Xsrcauto,phisrc ), Order( Xrngauto,phirng ) );
    phi.order := ord;
    return ord;
end;

##############################################################################
##
#F  XModMorphismOps.Print            special print commands for XMod morphisms
##

XModMorphismOps.Print := function( mor )

    local name;

    if IsBound( mor.name ) then
        name := mor.name;
    else
        name := "<[?]>->[?]>";
    fi;
    Print( "Morphism of crossed modules ", name );
end;

#############################################################################
##
#F  XModMorphismPrint              print details of a crossed module morphism
##

XModMorphismPrint := function( mor )

    local  morsrc, morrng, X, Y, name, ok;

    if IsBound( mor.isXModMorphism ) then
        ok := mor.isXModMorphism;
    else
        ok := IsXModMorphism( mor );
    fi;
    if not IsBound( mor.name ) then
        XModMorphismName( mor );
    fi;
    name := mor.name;

    X := mor.source;
    Y := mor.range;
    morsrc := mor.sourceHom;
    morrng := mor.rangeHom;
    if not ok then
        Print( "WARNING: not a morphism of crossed modules! :-\n" );
    else
        Print( "Morphism of crossed modules :- \n" );
    fi;
    Print( ": Source = ", X, " with generating sets:\n  " );
    Print( morsrc.source.generators, "\n  ", morrng.source.generators, "\n" );
    if ( Y = X ) then
        Print( ": Range = Source\n" );
    else
        Print( ":  Range = ", Y, " with generating sets:\n  " );
        Print( morsrc.range.generators, "\n  ",
               morrng.range.generators, "\n" );
    fi;
    Print( ": Source Homomorphism maps source generators to:\n" );
    Print( "  ", morsrc.genimages, "\n" );
    Print( ": Range Homomorphism maps range generators to:\n" );
    Print( "  ", morrng.genimages, "\n" );
    Print( ": isXModMorphism? ", ok, "\n" );
end; 

###############################################################################
##
#F  XModMorphism                          construct crossed module homomorphism
##

XModMorphism := function( arg )

    local  X, Y, L, ismor, mor, morsrc, morrng, Xsrc, Ysrc, Xrng, Yrng, 
           nargs, usage, x2, x1, y2, z2, y1, z1, gensrc, genrng;

    nargs := Length( arg );
    usage := "\nUsage: XModMorphism( X, Y, [sigma,rho], [true|false] );\n";
    if not ( ( nargs > 2 ) and ( nargs < 5 ) ) then
        Print( usage );
        return false;
    fi;
    X := arg[1];
    Y := arg[2];
    if ( not IsXMod( X ) ) or ( not IsXMod( Y ) ) then
        Print( usage );
        Error( "First two arguments must be crossed modules" );
        return false;
    fi;
    Xsrc := X.source;
    Ysrc := Y.source;
    Xrng := X.range;
    Yrng := Y.range;
    L := arg[3];
    if ( not IsList( L ) ) or ( not ( Length( L ) = 2 ) ) then
        Print( usage );
        Error( "Third argument must be a list with 2 morphisms" );
        return false;
    fi;
    morsrc := L[1];
    morrng := L[2];
    if    ( not IsHomomorphism( morsrc ) )
       or ( not ( morsrc.source = Xsrc ) )
       or ( not ( morsrc.range = Ysrc ) ) then
        Print( usage );
        Error( "L[1] must be a homomorphism from Xsrc to Ysrc" );
        return false;
    fi;
    if    ( not IsHomomorphism( morrng ) )
       or ( not ( morrng.source = Xrng ) )
       or ( not ( morrng.range = Yrng ) ) then
        Print( usage );
        Error( "L[2] must be a homomorphism from Xrng to Yrng" );
        return false;
    fi;
    ismor := ( ( nargs = 4 ) and IsBool( arg[4] ) and arg[4] );
    mor := rec( );
    mor.sourceHom := morsrc;
    mor.rangeHom := morrng;
    mor.source := X;
    mor.range := Y;
    mor.name := XModMorphismName( mor );
    mor.domain := mor.sourceHom.domain;   # = Mappings ?
    if not ismor then
        ismor := IsXModMorphism( mor );
    fi;
    if ismor then
        mor.isXModMorphism := true;
        mor.operations := XModMorphismOps;
    fi;
    return mor;
end;

##############################################################################
##
#F  XModMorphismOps.CompositeMorphism             composite of 2 XModMorphisms
##

XModMorphismOps.CompositeMorphism := function( phi, psi )

    local  XX, ZZ, Xsrc, Xrng, Zsrc, Zrng, genXsrc, genXrng,
           phisrc, phirng, psisrc, psirng, homgensrc,  
           homgenrng, homsrc, homrng, hom, ok;

    if not ( IsXModMorphism( phi ) and IsXModMorphism( psi ) ) then
        Error( "parameters must both be morphisms of crossed modules" );
        return false;
    fi; 
    if (not (phi.range = psi.source) ) then
        Error( "Range of first morphism <> source of second morphism" );
        return false;
    fi;

    XX := phi.source;
    ZZ := psi.range;
    Xsrc := XX.source;
    Xrng := XX.range;
    Zsrc := ZZ.source;
    Zrng := ZZ.range;

    genXsrc := Xsrc.generators;
    genXrng := Xrng.generators;
    phisrc := phi.sourceHom;
    phirng := phi.rangeHom;
    psisrc := psi.sourceHom;
    psirng := psi.rangeHom;
    homgensrc := List( genXsrc, g -> Image( psisrc, Image( phisrc, g ) ) );
    homgenrng := List( genXrng, g -> Image( psirng, Image( phirng, g ) ) );
    homsrc := GroupHomomorphismByImages( Xsrc, Zsrc, genXsrc, homgensrc );
    homrng := GroupHomomorphismByImages( Xrng, Zrng, genXrng, homgenrng );
    hom := XModMorphism( XX, ZZ, [homsrc,homrng] );
    ok := IsXModMorphism( hom );   
    return hom;
end;

##############################################################################
##
#F  XModOps.IdentityMorphism                     identity morphism for an XMod
##

XModOps.IdentityMorphism := function( X )

    local  idX, idsrc, idrng, ok;

    if not IsXMod( X ) then
        Error( "parameter must be a crossed module" );
    fi;
    idsrc := InclusionMorphism( X.source, X.source );
    idrng := InclusionMorphism( X.range, X.range );
    idX := XModMorphism( X, X, [idsrc, idrng] );
    ok := IsXModMorphism( idX );
    return idX;
end;

##############################################################################
##
#F  XModOps.InnerAutomorphism                    action of  r in Xrange  on  X
##

XModOps.InnerAutomorphism := function( X, r )

    local  Xsrc, Xrng, mor, morsrc, gensrc, imsrc, autr, genrng, morrng, ok;

    if not IsXMod( X ) then
        Error( "first parameter must be a crossed module" );
    elif not ( r in X.range ) then
        Error( "second parameter must be an element of the range of X" );
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    autr := Image( X.action, r );
    gensrc := Xsrc.generators;
    imsrc := List( gensrc, s -> Image( autr, s ) );
    morsrc := GroupHomomorphismByImages( Xsrc, Xsrc, gensrc, imsrc );
    genrng := Xrng.generators;
    morrng := InnerAutomorphism( Xrng, r );
    morrng.generators := genrng;
    morrng.genimages := List( genrng, x -> Image( morrng, x ) );
    mor := XModMorphism( X, X, [morsrc, morrng] );
    ok := IsXModMorphism( mor );
    return mor;
end;

##############################################################################
##
#F  XModMorphismOps.InverseMorphism                inverse of an XMod morphism
##

XModMorphismOps.InverseMorphism := function( mor )

    local  invmor, invsrc, invrng, ok;

    if not IsXModMorphism( mor ) then
        Error( "parameter must be a morphism of crossed modules" );
    fi;
    if not IsMonomorphism( mor ) then
        Error( "parameter must be a monomorphism of crossed modules" );
    fi;
    if not (     IsIsomorphism( mor.sourceHom )
             and IsIsomorphism( mor.rangeHom ) ) then
    Error( "source and range homomorphisms not isomorphisms" );
    fi;
    invsrc := InverseMapping( mor.sourceHom );
    invrng := InverseMapping( mor.rangeHom );
    invmor := XModMorphism( mor.range, mor.source, [invsrc,invrng] );
    ok := IsXModMorphism( invmor );
    return invmor;
end;

##############################################################################
##
#F  XModMorphismOps.Kernel                find Kernel XMod of an XMod morphism
##

XModMorphismOps.Kernel := function( mor )

    local  XX, ZZ, J, K, morrng, morsrc ;

    if ( not IsXModMorphism( mor ) ) then
        Error( "requires a crossed module morphism as input" );
        return false;
    fi;
    if IsBound( mor.kernel ) then
        return mor.kernel;
    fi;
    if not IsBound( mor.name ) then
        XModMorphismName( mor );
    fi;

    XX := mor.source;
    morsrc := mor.sourceHom;
    morrng := mor.rangeHom;
    J := Kernel( morsrc );
    J.parent := Parent( XX.source );
    K := Kernel( morrng );
    K.parent := Parent( XX.range );
    ZZ := SubXMod( XX, J, K );
    ZZ.name := Concatenation( "Ker", mor.name );
    mor.kernel := ZZ;
    return ZZ;
end;

###############################################################################
##
#F  ImageXModMorphism            find Image of a SubXMod under an XMod morphism
##

ImageXModMorphism := function( mor, ZZ )

    local  XX, YY, Xsrc, Xrng, Zsrc, Zrng, W, P, Q, J, K, morsrc, morrng ;

    if ( not IsXModMorphism( mor ) ) then
        Error( "first parameter must be a crossed module morphism" );
        return false;
    fi;
    if not IsBound( mor.name ) then
        XModMorphismName( mor );
    fi;
    XX := mor.source;
    YY := mor.range;
    if not IsXMod( ZZ ) then
        Error( "second parameter must be a crossed module" );
        return false;
    fi;
    if ( ( XX = ZZ ) and ( Size( XX ) = Size( YY ) ) ) then
        return YY;
    fi;
    Xsrc := XX.source;
    Xrng := XX.range;
    Zsrc := ZZ.source;
    Zrng := ZZ.range;
    morsrc := mor.sourceHom;
    morrng := mor.rangeHom;
    J := Image( morsrc, Zsrc );
    K := Image( morrng, Zrng );
    if ( ( Size( J ) = 1 ) and ( Size( K ) = 1 ) ) then
        W := IdentitySubXMod( YY );
    else
        W := SubXMod( YY, J, K );
    fi;
    W.name := Concatenation( "[Im(", ZZ.name, ") by ", mor.name, "]" );
    mor.image := W;
    return W;
end;

###############################################################################
##
#F  SourceXModXPModMorphism       top XMod from a morphism of crossed P-modules
##

SourceXModXPModMorphism := function( mor )
        
    local  X, Y, Xsrc, Xrng, Ysrc, Yrng, Ybdy, y, z, S, Sbdy, Saut, 
           genXsrc, genXrng, genYsrc, ngXrng, ngYsrc, inn, innaut,
           act, images, idXsrc, idXrng;

    if ( not IsXModMorphism( mor ) ) then
       Error( "Input parameter must be a crossed module morphism" );
       return false;
    fi;

    X := mor.source;
    Y := mor.range;
    Xsrc := X.source;
    Ysrc := Y.source;
    Xrng := X.range;
    Yrng := Y.range;
    idXsrc := InclusionMorphism( Xsrc, Xsrc );
    idXrng := InclusionMorphism( Xrng, Xrng );

    if ( mor.rangeHom = InclusionMorphism( Xrng, Xrng ) ) then
        mor.rangeHom := idXrng;
    fi;
    if   not ( X.range = Y.range)
      or not ( mor.rangeHom = idXrng ) then
        Error("Not a morphism of crossed modules having a common range");
    fi;
    genXsrc := Xsrc.generators;
    genXrng := Xrng.generators;
    genYsrc := Ysrc.generators;
    ngXrng := Length( genXrng);
    ngYsrc := Length( genYsrc );
    Ybdy := Y.boundary;
    innaut := [ ];
    for y in genYsrc do
        z := Image( Ybdy, y );
        images := List( genXsrc, x -> x^z );
        inn := GroupHomomorphismByImages( Xsrc, Xsrc, genXsrc, images );
        Add( innaut, inn );
    od;
    Saut := Group( innaut, idXsrc );
    act := GroupHomomorphismByImages( Ysrc, Saut, genYsrc, innaut );
    S := XMod( mor.sourceHom, act );
    return S; 
end;

##############################################################################
##
#F  XModMorphismAutoPerm                         ( p \in P  -->  mor: X -> X )
##

XModMorphismAutoPerm := function( X, p )

    local  Xsrc, Xrng, gensrc, genrng, P, p2arng, p2asrc,
           sigma0, imsigma, sigma, rho0, imrho, rho, mor;

    if not IsXMod( X ) then
        Error( "First parameter must be a crossed module" );
    fi;
    if not ( IsBound( X.actorSquare ) and
             IsBound( X.actorSquare.automorphismPermGroup ) ) then
        Error( "<X> has no automorphismPermGroup so far" );
    fi;
    P := AutomorphismPermGroup( X );
    if not ( p in P ) then
        Error( "Parameter 2 must be an element of <P>" );
    fi;
    X := P.xmod;
    Xsrc := X.source;
    Xrng := X.range;
    gensrc := Xsrc.generators;
    genrng := Xrng.generators;
    p2arng := Xrng.automorphismGroup.automorphismPair.p2a;
    p2asrc := Xsrc.automorphismGroup.automorphismPair.p2a;
    sigma0 := Image( p2asrc, Image( P.projsrc, p ) );
    if ( sigma0 = IdentityMapping( Xsrc ) ) then
        sigma0 := InclusionMorphism( Xsrc, Xsrc );
    fi;
    imsigma := List( gensrc, s -> Image( sigma0, s ) );
    sigma := GroupHomomorphismByImages( Xsrc, Xsrc, gensrc, imsigma );
    rho0 := Image( p2arng, Image( P.projrng, p ) );
    if ( rho0 = IdentityMapping( Xrng ) ) then
        rho0 := InclusionMorphism( Xrng, Xrng );
    fi;
    imrho := List( genrng, r -> Image( rho0, r ) );
    rho := GroupHomomorphismByImages( Xrng, Xrng, genrng, imrho );
    mor := XModMorphism( X, X, [ sigma, rho ], true );
    return mor;
end;

##############################################################################
##
#F  IsomorphicXMod                        xmod X + isos [sigma,rho] ==> xmod Y
##

IsomorphicXMod := function( arg )

    local  nargs, usage, X, sigma, rho, XS, XR, YS, YR, Xbdy, genXS, genYS,
           genXR, genYR, Xact, genXact, actlen, imXbdy, imYbdy, Ybdy,
           genYaut, Yaut, Yact, idYS, r, a, ima, i, r, s, Y;

    nargs := Length( arg );
    usage := "\nUsage: IsomorphicXMod( X, [sigma,rho] );\n";
    if not ( (nargs = 2) and IsXMod( arg[1] ) and IsList( arg[2] )
                         and ( Length( arg[2] ) = 2 ) ) then
        Print( usage );
        return false;
    fi;
    X := arg[1];
    sigma := arg[2][1];
    rho := arg[2][2];
    if not ( IsIsomorphism( sigma ) and IsIsomorphism( rho )
           and ( sigma.source = X.source ) and ( rho.source = X.range ) ) then
        Print( usage );
        return false;
    fi;
    Xbdy := X.boundary;
    genXR := rho.generators;
    genYR := rho.genimages;
    genXS := sigma.generators;
    genYS := sigma.genimages;
    YS := sigma.range;
    YR := rho.range;
    imXbdy := List( genXS, s -> Image( Xbdy, s ) );
    imYbdy := List( imXbdy, r -> Image( rho, r ) );
    Ybdy := GroupHomomorphismByImages( YS, YR, genYS, imYbdy );
    Xact := X.action;
    genXact := Xact.generators;
    actlen := Length( genXact );
    genYaut := 0 * [1..actlen];
    for i in [1..actlen] do
        r := genXR[i];
        a := Image( Xact, r );
        ima := List( genXS, s -> Image( sigma, Image( a, s ) ) );
        genYaut[i] := GroupHomomorphismByImages( YS, YS, genYS, ima );
    od;
    idYS := InclusionMorphism( YS, YS );
    Yaut := Group( genYaut, idYS );
    Yact := GroupHomomorphismByImages( YR, Yaut, genYR, genYaut );
    Y := XMod( Ybdy, Yact );
    return Y;
end;

# end of file xmodmor.g
