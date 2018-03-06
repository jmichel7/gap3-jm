###############################################################################
##
##  xmod.g                      for GAP 3.4                    version 10/ 1/97
##
###############################################################################
##
#A  xmod.g                      GAP library                       Chris Wensley
#A                                                                    Murat Alp
#Y  Copyright
##
##  This file contains functions that manipulate crossed modules,
##  pre-crossed modules, group graphs and associated constructions.
##
#H  $Log: xmod.g,v $
#H  Revision 1.1  1997/03/27 13:35:54  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

###############################################################################
##
#F  IsXMod                            check that the crossed module axioms hold
##

IsXMod := function( X )
    
    local Xsrc, Xrng, hom, a, aut, act, gensrc, ngsrc, genrng, ngrng, 
          ssrc, x1, y1, z1, x2, y2, z2, w2;

    if not IsRec( X ) then
       return false;
    fi;
    if IsBound( X.isXMod ) then
        return X.isXMod;
    fi;
    if ( not ( IsBound( X.range )  ) or not ( IsBound( X.source ) ) ) then
        return false;
    fi;
    Xrng := X.range;
    genrng := X.range.generators;
    ngrng := Length( genrng );
    Xsrc := X.source;
    gensrc := X.source.generators;
    ngsrc := Length( gensrc );
    ssrc := Size( Xsrc );
    
    # Check  X.boundary: X.source -> X.range
    if not (     IsBound( X.boundary ) 
             and ( X.boundary.source = X.source )
             and ( X.boundary.range = X.range ) ) then
        if ( XModPrintLevel > 1 ) then
            Print( "Error: require  X.boundary : X.source -> X.range \n" );
        fi;
        return false;
    fi;
    # checking  IsHomomorphism(hom) gives cokernel error when Xsrc = [ ]
    if ( ssrc > 1 ) then
        if not IsHomomorphism( X.boundary ) then
            if ( XModPrintLevel > 1 ) then
                Print( "Error:  the boundary map is NOT a homomorphism!\n" );
            fi;
            return false;
        fi;
    fi;
    hom := X.boundary;

    # Check  X.aut  is a group of automorphisms  X.source -> X.source
    if not ( IsBound( X.aut ) and IsGroup( X.aut ) ) then
        if ( XModPrintLevel > 1 ) then
            Print( "Error: group of actions on X.source does not exist! \n" );
        fi;
        return false;
    fi;
    aut := X.aut;
    if ( aut = Group( InclusionMorphism( Xsrc, Xsrc ) ) ) then
       X.isTrivialAction := true;
    else
        a := aut.generators[1];
        if not ( ( a.source = X.source )
                 and ( a.range = X.source ) 
                 and IsIsomorphism( a )  ) then
            if ( XModPrintLevel > 1 ) then
                Print( "Require automorphism  X.aut.1  on  X.source\n" );
            fi;
            return false;
        fi;
    fi;

    # Check  X.action: X.range -> X.aut
    if not (     IsBound( X.action ) 
             and ( X.action.source = X.range )
             and ( X.action.range = aut )    ) then
        if ( XModPrintLevel > 1 ) then
            Print( "Error: require  X.action : X.range -> X.aut \n" );
        fi;
        return false;
    fi;
    if ( ( XModPrintLevel > 1 ) and ( Size( aut ) = 1 ) ) then
        Print( "X.action trivial => not checking a homomorphism!\n" );
    else
        if not IsHomomorphism( X.action ) then
            if ( XModPrintLevel > 1 ) then
                Print( " X.action is not a homomorphism|\n" );
            fi;
            return false;
        fi;
    fi;
    act := X.action;

    if ( XModPrintLevel > 2 ) then
        Print( "Checking  CM1) hom(x2^x1) = (hom(x2))^x1 \n" );
    fi;
    for x1 in genrng do
        for x2 in gensrc do
            # Print( "x1,x2 = ", x1, ",  ", x2, "\n" );
            y1 := ( x2 ^ ( x1^act ) ) ^ hom;
            z1 := ( x2 ^ hom ) ^ x1;
            if ( y1 <> z1 ) then
                if ( XModPrintLevel > 1 ) then
                    Print( "CM1) fails at  x1 = ", x1, ",  x2 = ", x2, "\n" );
                    Print( "  hom(x2^x1) = ", y1, "\n" );
                    Print( "(hom(x2))^x1 = ", z1, "\n" );
                fi;
                return false;
            fi;
        od;
    od;

    if ( XModPrintLevel > 2 ) then
        Print( "Checking  CM2) x2^(hom(y2)) = x2^y2 \n" );
    fi;
    for x2 in gensrc do
        for y2 in gensrc do
            # Print( "x2,y2 = ", x2, ",  ", y2, "\n" );
            z2 := x2 ^ ((y2 ^ hom) ^ act);
            w2 := x2 ^ y2;
            if ( z2 <> w2 ) then
                if ( XModPrintLevel > 1 ) then
                    Print( "CM2) fails at  x2 = ", x2, ",  y2 = ", y2, "\n" );
                    Print( "x2^(hom(y2)) = ", z2, "\n" );
                    Print( "       x2^y2 = ", w2, "\n" );
                fi;
                return false;
            fi;
        od;
    od;
    X.isXMod := true;
    return true;
end;

#############################################################################
##
#F  XModPrint                               print details of a crossed module
##

XModPrint := function( X )

    local  name, bdy, act, aut, len, i, ispar, Xsrc, Xrng, gensrc, genrng,
           mor, triv, a, Arec;

    if not ( IsBound( X.isXMod ) and X.isXMod ) then
        Print( "XModPrint can only print a crossed module\n" );
        return false;
    fi;
    if not IsBound( X.name ) then
        X.name := "[?]";
    fi;
    name := X.name;
    Xsrc := X.source;
    Xrng := X.range;
    gensrc := Xsrc.generators;
    genrng := Xrng.generators;

    Print( "\nCrossed module ", name, " :- \n" );
    ispar := not IsBound( Xsrc.parent );
    if ( ispar and IsBound( Xsrc.name ) ) then
        Print( ": Source group ", Xsrc );
    elif ( IsBound( Xsrc.parent ) and IsBound( Xsrc.parent.name ) ) then
        Print( ": Source group has parent ( ", Xsrc.parent, " ) and" );
    else
        Print( ": Source group" );
    fi;
    Print( " has generators:\n" );
    Print( "  ", gensrc, "\n" );
    ispar := not IsBound( Xrng.parent );
    if ( ispar and IsBound( Xrng.name ) ) then
        Print( ": Range group = ", Xrng );
    elif ( IsBound( Xrng.parent ) and IsBound( Xrng.parent.name ) ) then
        Print( ": Range group has parent ( ", Xrng.parent, " ) and" );
    else
        Print( ": Range group" );
    fi;
    Print( " has generators:\n" );
    Print( "  ", genrng, "\n" );
    Print( ": Boundary homomorphism maps source generators to:\n" );
    bdy := X.boundary;
    if not ( bdy.generators = gensrc ) then
        Print( "! Warning: Xsrc.generators <> Xbdy.generators!\n" );
    fi;
    Print( "  ", bdy.genimages, "\n" );
    act := X.action;
    if not ( act.generators = genrng ) then
        Print( "! Warning: Xact.generators <> Xrng.generators!\n" );
    fi;
    aut := X.aut;
    triv := ( aut = Group( InclusionMorphism( Xsrc, Xsrc ) ) );
    if IsBound( act.isEmbedding ) then
        Print( ": Action homomorphism is the SemidirectProduct Embedding\n" );
    else
        len := Length( genrng );
        if ( len = 0 ) then
            triv := true;
        else
            for i in [1..len] do
                a := act.genimages[i];
                if ( a = IdentityMapping( Xsrc ) ) then
                    if not IsBound( a.generators ) then
                        a.generators := gensrc;
                    fi;
                    if not IsBound( a.genimages ) then
                        a.genimages := a.generators;
                    fi;
                fi;
                if not ( a.generators = gensrc ) then
                    Print( "! Warning: Xautgen[",i,"] <> Xsrc.generators!\n" );
                fi;
            od;
        fi;
        if not triv then
            Print( ": Action homomorphism maps" );
            Print( " range generators to automorphisms:\n" );
            for i in [1..len] do
                Print( "  ", genrng[i], " --> { source gens --> " );
                Print( act.genimages[i].genimages, " }\n" );
            od;
        fi;
        if triv then
            Print( "  The automorphism group is trivial\n" );
        else
            if ( len = 1 ) then
                Print( "  This automorphism generates" );
            else
                Print( "  These ", len, " automorphisms generate" );
            fi;
            Print( " the group of automorphisms.\n" );
        fi;
    fi;
    if ( IsBound( X.kernel ) and IsBound( X.kernel.generators ) ) then
        Print( ": Kernel of the crossed module has generators:\n" );
        Print( "  ", X.kernel.generators, "\n" );
    fi;
    if IsBound( X.cat1 ) then
        Print( ": Associated cat1-group = ", X.cat1, "\n" );
    fi;
    if IsBound( X.isInducedXMod ) then
        Print( ": Induced XMod from ", X.xmod, "  with source morphism:\n" );
        mor := X.morphism;
        Print( "  ", mor.sourceHom.generators, "\n   --> " );
        Print( mor.sourceHom.genimages, "\n" );
        # Print( "  ", mor.rangeHom.generators, "\n   --> " );
        # Print( mor.rangeHom.genimages, "\n" );
        if IsBound( X.sourceId ) then
            Print( "  Source has GroupId = " );
            Print( X.sourceId.catalogue, ", ", X.sourceId.names, "\n" );
        fi;
    fi;
    if IsBound( X.derivations ) then
        Print( ": ", X.derivations, "\n" );
    fi;
    if IsBound( X.actorSquare ) then
        Arec := X.actorSquare;
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

###############################################################################
##
#F  XModName                             Select a name for the crossed module X
##

XModName := function( X )

    local  sn, rn, name;

    # any existing name if NOT automatically returned
    # unless X is one of the actor square crossed modules
    if not IsXMod ( X ) then
        Error( " Input parameter must be a crossed module" );
        return "false";
    fi;
    if ( IsBound( X.name ) and
         ( IsBound( X.isActor ) or IsBound( X.isInnerActor ) or
           IsBound( X.isLue ) or IsBound( X.isWhitehead ) or
           IsBound( X.isNorrie ) or IsBound( X.isCentre ) ) ) then
        return X.name;
    fi;
    if IsBound( X.source.name ) then
        sn := X.source.name;
    else
        sn := "?";
    fi;
    if IsBound( X.range.name ) then
        rn := X.range.name;
    else
        rn := "?";
    fi;
    name := Concatenation( "[", sn, "->", rn, "]" );
    X.name := name;
    return name;
end;

###############################################################################
##
#F  XMod              crossed module from given source, range, boundary, action
##

XMod := function( arg )

    local  nargs, usage, R, S, genR, genS, bdy, aut, act0, imact0, act, imact,
           new, i, a0, ima, a, X, name, ok;

    nargs := Length( arg );
    usage := "Usage: XMod( boundary, action )";
    if ( nargs <> 2 ) then
        Error( usage );
    fi;
    bdy := arg[1];
    if not IsHomomorphism( bdy ) then
        if ( XModPrintLevel > 1 ) then
            Print( "Boundary is not a homomorphism:\n", bdy, "\n" );
        fi;
        Error( usage );
    fi;
    S := bdy.source;
    genS := S.generators;
    R := bdy.range;
    genR := R.generators;
    act0 := arg[2];
    if not ( IsHomomorphism( act0 ) and ( act0.source = R ) ) then
        if ( XModPrintLevel > 1 ) then
            if not ( act0.source = R ) then
                Print( "The range group is not the source of the action.\n" );
            else
                Print( "<act> is not a homomorphism:\n", act0, "\n" );
            fi;
        fi;
        Error( usage );
    fi;
    aut := act0.range;
    if not (     IsAutomorphismGroup( aut )
             and ( aut.identity = InclusionMorphism( S, S ) ) ) then
        if ( XModPrintLevel > 1 ) then
            if not IsAutomorphismGroup( aut ) then
                Print( "<aut> is not an automorphism group:\n", aut, "\n" );
            else
                Print( "aut.identity <> InclusionMorphism( S, S )\n" );
            fi;
        fi;
        Error( usage );
    fi;
    imact0 := act0.genimages;
    if not ( act0.generators = genR ) then
        imact := List( genR, r -> Image( act0, r ) );
        act0 := GroupHomomorphismByImages( R, aut, genR, imact );
    fi; 
    new := false;
    imact0 := act0.genimages;
    imact := Copy( imact0 );
    for i in [ 1..Length( imact0 ) ] do
        a0 := imact0[i];
        if ( a0.generators <> genS ) then
            new := true;
            ima := List( genS, s -> Image( a0, s ) );
            a := GroupHomomorphismByImages( S, S, genS, ima );
            imact[i] := a;
        fi;
    od;
    if new then
        act := GroupHomomorphismByImages( R, aut, genR, imact );
    else
        act := act0;
    fi;

    X := rec( );
    X.isDomain := true;
    X.isParent := true;
    X.source := S;
    X.range := R;
    X.boundary := bdy;
    X.action := act;
    X.aut := aut;
    ok := IsXMod( X );
    if not ok then
        Error( "data fails to define a crossed module" );
    fi;
    X.operations := XModOps;
    name := XModName( X );
    return X;
end;

##############################################################################
##
#F  ConjugationXMod                 create a crossed module from normal R in S
##

ConjugationXMod := function( arg )

    local  nargs, X, bdy, act, aut, genR, ngR, genS, ngS, idS, name, a,
           autgen, imautgen, phi, j, r, s, R, S, genN, f2pN, imgenN;

    # check the number of arguments and their type
    nargs := Length( arg );
    if ( ( nargs < 1 ) or ( nargs > 3 ) ) then
        Error( "\nUsage: ConjugationXMod( S [, R ] )\n" );
    fi;
    S := arg[1];
    if ( nargs = 1 ) then
        R := S;
    else
        S := arg[2];
        R := arg[1];
        if not ( IsPermGroup( S ) and IsPermGroup( R ) ) then
            Error( "the 2 arguments must both be PermGroups" );
        fi;
        # silly fix for trivial subgroups ??
        if ( ( Size( R ) = 1 ) and ( Size( S ) = 1 ) ) then
            R := S;
        fi;
        if not ( IsSubgroup( R, S ) and IsNormal( R, S ) ) then
            Error( "Second parameter must be a normal subgroup of the first" );
        fi;
    fi;
    genR := R.generators;
    genS := S.generators;
    ngR := Length( genR );
    ngS := Length( genS );
    idS := InclusionMorphism( S, S );
    bdy := InclusionMorphism( S, R );
    bdy.name := "bdy";
    autgen := [ ];
    for r in genR do
        imautgen := List( genS, s -> s^r );
        a := GroupHomomorphismByImages( S, S, genS, imautgen );
        Add( autgen, a );
    od;
    aut := Group( autgen, idS );
    if ( XModPrintLevel > 1 ) then
        Print( "Group of conjugations has size ", Size(aut), "\n" );
    fi;
    act := GroupHomomorphismByImages( R, aut, genR, autgen );
    if ( not IsHomomorphism( act ) ) then
        Error( "action is not a homomorphism" );
    fi;

    X := XMod( bdy, act );
    X.isConjugationXMod := true;
    if ( Length( aut.generators ) = 0 ) then
        X.isTrivialAction := true;
    fi;
    return X;
end;

###############################################################################
##
#F  TrivialActionXMod      crossed module from  f : abelian G -> Z <= centre(R)
## 

TrivialActionXMod := function( f )

    local  R, ZR, S, X, aut, act, name;

    if not IsHomomorphism( f ) then
        Error( "Input parameter must be a homomorphism" );
    fi;
    S := f.source;
    if not IsAbelian( S ) then
        Error( "the source of  f  must be abelian" );
    fi;
    R := f.range;
    ZR := Centre( R );
    if not IsSubgroup( ZR, Image( f, S ) ) then
        Error( "image of source must lie in the centre of range" );
    fi;
    aut := Group( InclusionMorphism( S, S ) );
    act := ZeroMorphism( R, aut );

    X := XMod( f, act );
    X.isTrivialAction := true;
    return X;
end;

###############################################################################
##
#F  CentralExtensionXMod   create a crossed module from central surjection S->R
##

CentralExtensionXMod := function( f )

    local  R, S, ZS, ker, X, bdy, act, aut, genR, ngR, genS,
           name, autgen, imautgen, id, phi, j, s;

    if ( not IsHomomorphism( f ) or not IsSurjective( f ) ) then
        Error( "Input parameter must be a surjective homomorphism" );
    fi;
    S := f.source;
    R := f.range;
    ZS := Centre( S );
    ker := Kernel( f );
    if not IsSubgroup( ZS, ker ) then
        Error( "Kernel of surjection is not central" );
    fi;
    genS := S.generators;
    genR := f.genimages;
    ngR := Length( genR );
    autgen := [ ];
    for s in genS do
        Add( autgen, List( genS, t -> t^s ) );
    od;
    imautgen := [ ];
    for j in [1..ngR] do 
        phi := GroupHomomorphismByImages( S, S, genS, autgen[j] );
        Add( imautgen, phi );
    od;
    id := InclusionMorphism( S, S );
    aut := Group( imautgen, id );
    if ( XModPrintLevel > 1 ) then
        Print( "Group of conjugations has size ", Size(aut), "\n" );
    fi;
    act := GroupHomomorphismByImages( R, aut, genR, imautgen );
    if ( not IsHomomorphism( act ) ) then
        Error( "action is not a homomorphism" );
    fi;

    X := XMod( f, act );
    X.isCentralExtensionXMod := true;
    if ( Length( aut.generators ) = 0 ) then
        X.isTrivialAction := true;
    fi;
    return X;
end;

###############################################################################
##
#F  AutomorphismXMod       create a crossed module  G --> A \in [Inn(G),Aut(G)]
##  

AutomorphismXMod := function( arg )

    local  nargs, G, IA, A, pairA, a, X, P, genG, inn, img, inngen,
           g, idG, abelian;

    # check the number of arguments and their type
    nargs := Length( arg );
    if ( ( nargs < 1 ) or ( nargs > 2 ) ) then
        Error( "\nUsage:  AutomorphismXMod( G [, N] );\n" );
    fi;
    G := arg[1];
    if IsBound( G.automorphismXMod ) then
        return G.automorphismXMod;
    fi;
    if not IsPermGroup( G ) then
        Error( "First parameter must be a permutation group" );
        return false;
    fi;
    if ( nargs = 1 ) then
        A := AutomorphismGroup( G );
        A.isAutomorphismGroup := true;
        pairA := AutomorphismPair( A );
        P := pairA.perm;
    else
        A := arg[2];
        IA := InnerAutomorphismGroup( G );
        for inn in IA.generators do
            if not ( inn in A ) then
                Error( "Second parameter must include all inner autos.\n" );
            fi;
        od;
        pairA := AutomorphismPair( A );
        P := pairA.perm;
    fi;
    abelian := IsAbelian( G );
    genG := G.generators;
    if abelian then
        inngen := List( genG, g -> () );
    else
        inngen := [ ];
        for g in genG do
            a := InnerAutomorphism( G, g );
            img := Image( pairA.a2p, a );
            Add( inngen, img );
        od;
    fi;
    inn := GroupHomomorphismByImages( G, P, genG, inngen );

    X := XMod( inn, pairA.p2a );
    X.isAutomorphismXMod := true;
    if ( Length( A.generators ) = 0 ) then
        X.isTrivialAction := true;
    fi;
    if ( nargs = 1 ) then
        G.automorphismXMod := X;
    fi;
    return X;
end;

###############################################################################
##
#F  InnerAutomorphismXMod                 create a crossed module  G --> Inn(G)
##  

InnerAutomorphismXMod := function( G )

    local  IA, X;

    if not IsPermGroup( G ) then
        Error( "Parameter must be a permutation group" );
        return false;
    fi;
    IA := InnerAutomorphismGroup( G );
    if not IsBound( IA.morphismDomain ) then
        IA.morphismDomain := G;
    fi;
    X := AutomorphismXMod( G, IA );
    XModName( X );
    return X;
end;

##############################################################################
##
#F  IsRModuleRecord    tests for automorphism group acting on an abelian group
##

IsRModuleRecord := function( Rmod )

    local  M, A, R, ok, genA, message, pairA;

    if not IsRec( Rmod ) then
        return false;
    fi;
    if IsBound( Rmod.isRModule ) then
        return Rmod.isRModule;
    fi;
    message:="\nAn R-module is a record Rmod with fields .module .auto\n";
    if not ( IsBound( Rmod.module ) and IsGroup( Rmod.module )
                                    and IsAbelian( Rmod.module ) ) then
        Print( message );
        return false;
    fi;
    M := Rmod.module;
    if not IsBound( Rmod.auto ) then
        Print( message );
        return false;
    fi;
    A := Rmod.auto;
    ok := IsAutomorphismGroup( A );
    if not ok then
        Print( "Rmod.auto is not a group of automorphisms\n" );
        return false;
    fi;
    if Size( A ) = 1 then
        Print( "The group of automorphisms is trivial.\n" );
    else
        genA := A.generators[1];
        if not ( ( genA.source = M ) and ( genA.range = M ) ) then
            Print( "Rmod.auto is not an automorphism group" );
            Print( " of Rmod.module\n" );
            return false;
        fi;
    fi;
    pairA := AutomorphismPair( A );
    Rmod.perm := pairA.perm;
    Rmod.isRModule := true;
    return true;
end;

##############################################################################
##
#F  RModuleXMod                       0 : M --> G  where  G  acts on abelian M
##

RModuleXMod := function( R )

    local  M, G, X, bdy, ok, P, PG, gensrc, imsrc;

    if not IsRModule( R ) then
        Error( "Input parameter must be an R-module." );
    fi;
    M := R.module;
    G := R.auto;
    ok := IsAutomorphismGroup( G );
    P := R.perm;
    bdy := ZeroMorphism( M, P );
    if ( IsBound( M.name ) and not IsBound( P.name ) ) then
        P.name := Concatenation( "act(", M.name, ")" );
    fi;

    X := rec( );
    X.isDomain := true;
    X.isParent := true;
    X.source := M;
    X.boundary := bdy;
    X.range := P;
    X.aut := G;
    X.action := G.automorphismPair.p2a;
    X.isXMod := IsXMod( X );
    X.operations := XModOps;
    X.isRModuleXMod := true;
    if ( Length( G.generators ) = 0 ) then
        X.isTrivialAction := true;
    fi;
    X.name := XModName( X );
    return X;
end;

##############################################################################
##
#F  XModSelect                       construct standard type of crossed module
##

XModSelect := function( arg )

    local  size, gpnum, norm, type, nargs, usage, maxsize, names,
           pos, pos2, j, M, G, X, NG, N, iso, start, count, aut, R;
    
    maxsize := Cat1ListMaxSize;
    usage := "\nUsage:  XModSelect( size, gpnum, [, type, norm] )\n";
    nargs := Length( arg );
    if ( ( nargs = 0 ) or not
            ( ( arg[1] > 0 ) and ( arg[1] < maxsize ) ) ) then
        Print( usage );
        return false;
    fi;
    # find starting position of iso classes of groups of order size
    iso := NumbersOfIsomorphismClasses;
    count := 1;
    start := [ 1 ];
    for j in iso do
        count := count + j;
        Add( start, count );
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "  iso = ", iso, "\nstart = ", start, "\n" );
    fi;
    size := arg[1];
    pos := start[ size ];
    if ( size < maxsize ) then
        pos2 := start[ size + 1 ] - 1;
    else
        pos2 := Length( Cat1List );
    fi;
    names := List( [pos..pos2], n -> Cat1List[n][4] );

    if ( ( nargs = 1 ) or not ( arg[2] > 0 ) ) then
        Print( usage );
        return names;
    fi;
    gpnum := arg[2];
    j := pos + gpnum - 1;
    M := Cat1List[j];
    G  := Group(M[3], ( ));
    G.name := M[4];

    if ( nargs = 2 ) then
        type := "aut";
        norm := 0;
    fi;
    if ( nargs > 2 ) then
        type := arg[3];
        if not ( ( type = "conj" ) or ( type = "aut" ) 
                                 or ( type = "rmod" ) ) then
            Print( usage );
            Print( "third parameter must be a string:  conj, aut, or rmod\n" );
            return names;
        fi;
    fi;
    if ( type = "conj" ) then
        if ( nargs = 4 ) then
            norm := arg[4];
        else
            norm := 0;
        fi;
    fi;

    if ( gpnum > iso[size] ) then
        Print( "# classes of size ", size, " is ", iso[size], "\n" );
        Error( "gpnum > # classes of groups with this size" );
        return false;
    fi;
    if ( type = "aut" ) then 
        X := AutomorphismXMod( G );
    elif ( type = "conj" ) then
        NG := NormalSubgroups(G);
        if ( norm > Length(NG) ) then
            Error(" norm is bigger than length of NormalSubgroup" );
            return false;
        fi;
        if ( norm = 0 ) then
            norm := Length( NG );
        fi;
        N := NG[norm];
        if ( ( norm <> 1 ) and ( norm <> Length( NG ) ) ) then
            N.name := Concatenation( "N", String( norm ) );
        fi;
        X := ConjugationXMod( G, N );
    elif ( type = "rmod" ) then
        if not IsAbelian( G ) then
            Print( "G must be abelian to be an R-module\n" );
            return false;
        fi;
        aut := AutomorphismGroup( G );
        R := rec( );
        R.module := G;
        R.auto := aut;
        X := RModuleXMod( R );        
    else 
        Print( usage );
        return false;
    fi;
    return X;
end;

##############################################################################
##
#F  IsSubXMod
##

IsSubXMod := function( X, S )

    local  ok, Ssrc, Srng, gensrc, genrng, s, r, r1, r2, im1, im2;

    if not ( IsRec( X ) and IsRec( S ) ) then
        return false;
    fi;
    if not ( IsXMod( X ) and IsXMod( S ) ) then
        return false;
    fi;
    if ( IsBound( S.parent ) and ( S.parent = X ) ) then
        return true;
    fi;
    Ssrc := S.source;
    Srng := S.range;
    if not (     IsSubgroup( X.source, Ssrc )
             and IsSubgroup( X.range, Srng ) ) then
        if ( XModPrintLevel > 2 ) then
            Print( "IsSubgroup failure in IsSubXMod\n" );
        fi;
        return false;
    fi;
    ok := true;
    gensrc := Ssrc.generators;
    genrng := Srng.generators;
    for s in gensrc do
        if ( Image( X.boundary, s ) <> Image( S.boundary, s ) ) then
            ok := false;
        fi;
    od;
    if not ok then
        if ( XModPrintLevel > 2 ) then
            Print( "boundary maps have different images\n" );
        fi;
        return false;
    fi;
    for r in genrng do
        r1 := Image( X.action, r );
        r2 := Image( S.action, r );
        for s in gensrc do
            im1 := Image( r1, s );
            im2 := Image( r2, s );
            if ( im1 <> im2 ) then
                ok := false;
                if ( XModPrintLevel > 2 ) then
                    Print( "s,im1,im2 = ", [s,im1,im2], "\n" );
                fi;
            fi;
        od;
    od;
    if not ok then
        if ( XModPrintLevel > 2 ) then
            Print( "actions have different images\n" );
        fi;
        return false;
    fi;
    if IsParent( X ) then
        S.parent := X;
    else
        S.parent := X.parent;
    fi;
    return true;
end;

##############################################################################
##
#F  SubXMod                       creates SubXMod from Ssrc<=Xsrc & Srng<=Xrng
##

SubXMod := function( X, Ssrc, Srng )

    local  Xsrc, Xrng, Xbdy, Xact, Xaut, genSsrc, genSrng, Xname, Sname,
           S, Sbdy, Saut, Sact, r, innaut, genXrng, genXsrc, ssrc,
           idSsrc, imact, imgen, imbdy, imSsrc, imalpha, alpha, ok;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
        return false;
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    Xbdy := X.boundary;
    Xaut := X.aut;
    Xact := X.action;
    ok := true;
    if not ( IsSubgroup( Xsrc, Ssrc ) ) then
        Print( "Ssrc is not a subgroup of Xsrc\n" );
        ok := false;
    fi;
    if not ( IsSubgroup( Xrng, Srng ) ) then
        Print( "Srng is not a subgroup of Xrng\n" );
        ok := false;
    fi;

    ssrc := Size( Ssrc );
    genXsrc := X.source.generators;
    genXrng := X.range.generators;
    genSsrc := Ssrc.generators;
    genSrng := Srng.generators;
    idSsrc := InclusionMorphism( Ssrc, Ssrc );

    imgen := List( genSsrc, x -> Image( Xbdy, x ) );
    imSsrc := Subgroup( Xrng, imgen );
    if not IsSubgroup( Srng, imSsrc ) then
        if ( XModPrintLevel > 1 ) then
            Print( "Xbdy(Ssrc) is not a subgroup of Srng\n" );
        fi;
        ok := false;
    fi;

    if ok then
        Sbdy:= GroupHomomorphismByImages( Ssrc, Srng, genSsrc, imgen );
        innaut := [ ];
        for r in genSrng do
            alpha := Image( Xact, r );
            imgen := List( genSsrc, x -> Image( alpha, x ) );
            imalpha := Subgroup( Ssrc, imgen );
            if not ( IsSubgroup( Ssrc, imalpha ) ) then
                if ( XModPrintLevel > 1 ) then
                    Print( "Srng does not act correctly on Ssrc\n" );
                fi;
                ok := false;
            fi;
            if ok then
                alpha:=GroupHomomorphismByImages( Ssrc, Ssrc, genSsrc, imgen );
            fi;
            Add( innaut, alpha );
        od;
    fi;

    if ok then
        if ( ssrc = 1 ) then
            Saut := Group( InclusionMorphism( Ssrc, Ssrc ) );
            innaut := List( genSrng, s -> InclusionMorphism( Ssrc, Ssrc ) );
        else
            Saut := Group( innaut, idSsrc );
        fi;
        Sact := GroupHomomorphismByImages( Srng, Saut, genSrng, innaut );
        
        if ( not IsHomomorphism( Sact ) ) then
            Print( "Sact is not a homomorphism\n" );
            ok := false;
        fi;
    fi;

    S := rec();
    S.isDomain := true;
    if IsParent( X ) then
        S.parent := X;
    else
        S.parent := X.parent;
    fi;
    S.source := Ssrc;
    S.range := Srng;
    if not ok then
        S.isXMod := false;
    else
        S.boundary := Sbdy;
        S.action := Sact;
        S.aut := Saut;
        ok := IsXMod( S );
    fi;

    if ok then
        S.operations := XModOps;
        if ( ( Size( Srng ) = 1 ) and not IsBound( Srng.name ) ) then
            Srng.name := "I";
        fi;
        Xname := XModName( X );
        if (ssrc = 1 ) then
            Ssrc.name := "I";
            Sname := XModName( S );
        elif ( ( Ssrc = Xsrc ) and ( Srng = Xrng ) ) then
            Sname := Xname;
        elif ( IsBound( Ssrc.name ) or IsBound( Srng.name ) ) then
            Sname := XModName( S );
        elif IsBound( X.name ) then
            Sname := Concatenation( "[Sub", X.name, "]" );    
        else
            Sname := "[Sub[?->?]]";
        fi;
        S.name := Sname;
    fi;
    return S;
end;

##############################################################################
##
#F  IdentitySubXMod                           creates sub-crossed module []->[]
##

IdentitySubXMod := function( X )

    local  Xsrc, Xrng, idsrc, idrng, I;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    idsrc := Xsrc.identity;
    idrng := Xrng.identity;
    I := SubXMod( X, Subgroup( Xsrc, [idsrc] ), Subgroup( Xrng, [idrng] ) );
    I.name := Concatenation( "[Id", X.name, "]" );
    return I;
end;

##############################################################################
##
#F  IsNormalSubXMod                  tests to see whether a sub-xmod is normal
##

IsNormalSubXMod := function( X, S )

    local  xr, a, ss, im, xs, sr;

    if not ( IsRec( X ) and IsRec( S ) ) then
        return false;
    fi;
    if not IsSubXMod( X, S ) then
        return false;
    fi;
    for xr in X.range.generators do
        a := Image( X.action, xr );
        for ss in S.source.generators do
            im := Image( a, ss );
            if not ( im in S.source ) then
                if ( XModPrintLevel > 2 ) then
                    Print( "ss,xr,ss^xr = ", [ss,xr,im], "\n" );
                fi;
                return false;
            fi;
        od;
    od;
    for sr in S.range.generators do
        a := Image( X.action, sr );
        for xs in X.source.generators do
            im := xs^(-1) * Image( a, xs );
            if not ( im in S.source ) then
                if ( XModPrintLevel > 2 ) then
                    Print( "sr,xs,sr^(-1)*xs^sr = ", [sr,xs,im], "\n" );
                fi;
                return false;
            fi;
        od;
    od;
    return true;
end;

##############################################################################
##
#F  NormalSubXMods                                 Ysrc <| Xsrc & Yrng <| Xrng
##

NormalSubXMods := function( X )

    local  Xsrc, Xrng, Y, i, j, slen, rlen, norm, normsrc, normrng;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
        return false;
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    norm := [ ];
    normsrc := NormalSubgroups( Xsrc );
    normrng := NormalSubgroups( Xrng );
    slen := Length( normsrc );
    rlen := Length( normrng );
    for i in [ 1..slen ] do
        for j in [ 1..rlen ] do
            if ( ( i = 1 ) and ( j = 1 ) ) then
                Y := IdentitySubXMod( X );
            elif ( ( i = slen ) and ( j = rlen ) ) then
                Y := X;
            else
                Y := SubXMod( X, normsrc[i], normrng[j] );
            fi;
            if ( IsXMod( Y ) and IsNormalSubXMod( X, Y ) ) then   
                Add( norm, Y );
            fi;
        od;
    od;
    return norm;
end;

##############################################################################
##
#F  FactorXMod         constructs quotient of X by a normal sub-crossed module
##

FactorXMod := function( X, N )

    local  Xsrc, Xrng, Nsrc, Nrng, Qsrc, Qrng, Psrc, Prng, cosetsrc, cosetrng,
           repsrc, reprng, homsrc, homrng, genQsrc, genQrng, genPsrc, genPrng,
           q2psrc, q2prng, Xbdy, genNsrc, imsrc, im1, im2, imrng, Xact, Xaut,
           r, autgen, alpha, repim, id, ima, a, Pbdy, Paut, Pact, P,
           genXsrc, genXrng, natsrc, natrng, projsrc, projrng, nat, ok;

    if not ( IsXMod( X ) and IsXMod( N ) ) then
        Error( "functions requires two crossed modules as parameters" );
    fi;
    if not IsNormalSubXMod( X, N ) then
        Error( "2nd parameter not a normal sub-crossed module of the 1st" );
    fi;
    Xsrc := X.source;
    Xrng := X.range;
    Nsrc := N.source;
    Nrng := N.range;
    Qsrc := Xsrc/Nsrc;
    Qrng := Xrng/Nrng;
    cosetsrc := RightCosets( Xsrc, Nsrc );
    Psrc := Operation( Xsrc, cosetsrc, OnRight );
    homsrc := OperationHomomorphism( Xsrc, Psrc );
    cosetrng := RightCosets( Xrng, Nrng );
    Prng := Operation( Xrng, cosetrng, OnRight );
    homrng := OperationHomomorphism( Xrng, Prng );
    if ( IsBound( Xsrc.name ) and IsBound( Nsrc.name ) ) then
        Qsrc.name := Concatenation( "(", Xsrc.name, " / ", Nsrc.name, ")" );
        Psrc.name := Concatenation( "Perm(", Qsrc.name, ")" );
    fi;
    if ( IsBound( Xrng.name ) and IsBound( Nrng.name ) ) then
        Qrng.name := Concatenation( "(", Xrng.name, " / ", Nrng.name, ")" );
        Prng.name := Concatenation( "Perm(", Qrng.name, ")" );
    fi;
    genXsrc := Xsrc.generators;
    genXrng := Xrng.generators;
    genQsrc := Qsrc.generators;
    genQrng := Qrng.generators;
    genPsrc := Psrc.generators;
    genPrng := Prng.generators;
    q2psrc := GroupHomomorphismByImages( Qsrc, Psrc, genQsrc, genPsrc );
    q2prng := GroupHomomorphismByImages( Qrng, Prng, genQrng, genPrng );
    if not ( IsIsomorphism( q2psrc ) and IsIsomorphism( q2prng ) ) then
        Error( "isomorphism failure between quotients and reps" );
    fi;
    genNsrc := Nsrc.generators;
    Xbdy := X.boundary;
    repsrc := List( genQsrc, q -> q.element.representative );
    im2 := List ( repsrc, s -> FactorGroupElement( Nsrc, s ) );
    imsrc := List( im2, s -> Image( q2psrc, s ) );
    reprng := List( genQrng, q -> q.element.representative );
    im1 := List( repsrc, s -> Image( Xbdy, s ) );
    im2 := List( im1, r -> FactorGroupElement( Nrng, r ) );
    imrng := List( im2, c -> Image( q2prng, c ) ); 
    Pbdy := GroupHomomorphismByImages( Psrc, Prng, imsrc, imrng );
    Xaut := X.aut;
    Xact := X.action;
    autgen := [ ];
    for r in reprng do
        alpha := Image( Xact, r );
        repim := List( repsrc, s -> Image( alpha, s ) );
        im2 := List( repim, s -> FactorGroupElement( Nsrc, s ) );
        ima := List( im2, c -> Image( q2psrc, c ) );
        a := GroupHomomorphismByImages( Psrc, Psrc, imsrc, ima );
        Add( autgen, a );
    od;
    id := InclusionMorphism( Psrc, Psrc );
    Paut := Group( autgen, id );
    Paut .isAutomorphismGroup := true;
    Pact := GroupHomomorphismByImages( Prng, Paut, genPrng, autgen );
    P := XMod( Pbdy, Pact );
    natsrc := NaturalHomomorphism( Xsrc, Qsrc );
    im1 := List( genXsrc, s -> Image( natsrc, s ) );
    im2 := List( im1, c -> Image( q2psrc, c ) );
    projsrc := GroupHomomorphismByImages( Xsrc, Psrc, genXsrc, im2 );
    natrng := NaturalHomomorphism( Xrng, Qrng );
    im1 := List( genXrng, r -> Image( natrng, r ) );
    im2 := List( im1, c -> Image( q2prng, c ) );
    projrng := GroupHomomorphismByImages( Xrng, Prng, genXrng, im2 );    
    nat := XModMorphism( X, P, [ projsrc, projrng ] );
    ok := IsXModMorphism( nat );
    if not ok then
        Error( "natural morphism failure" );
    fi;
    P.isFactorXMod := true;
    P.numerator := X;
    P.denominator := N;
    P.naturalMorphism := nat;
    return P;
end;

#############################################################################
##
#V  XModOps                                  crossed module record operations
##

XModOps := OperationsRecord( "XModOps", DomainOps );

#############################################################################
##
#F  XModOps.\=                             define equality of crossed modules
##

XModOps.\= := function( X, Y )

    local isEql;

    if not ( IsXMod( X ) and IsXMod( Y ) ) then
        Error( "X,Y not XMods" );
    fi;

    isEql := (    ( X.source = Y.source )
              and ( X.range = Y.range )
              and ( X.boundary = Y.boundary )
              and ( X.aut = Y.aut )
              and ( X.action = Y.action )   );
    return isEql;
end;

##############################################################################
##
#F  XModOps.Print                  special print commands for crossed modules
##

XModOps.Print := function( X )

    local name;

    if IsBound( X.name ) then
        name := X.name;
    else
        name := "[?]";
    fi;
    if not ( name = "?" ) then
        Print( "Crossed module ", name );
    fi;
end;

##############################################################################
##
#F  XModOps.IsParent

XModOps.IsParent := function( X )
    return  not IsBound( X.parent );
end;

##############################################################################
##
#F  XModOps.Parent

XModOps.Parent := function( L )

    local  X, Y;

    # get the parent of the first xmod
    if IsBound( L[1].parent ) then
        X := L[1].parent;
    else
        X := L[1];
    fi;
    # check that all other xmods have the same parent
    for Y in L do
        if ( IsBound( Y.parent ) and ( Y.parent <> X ) ) then
            Error( "<X> and <Y> must have the same parent xmod" );
        elif ( not IsBound( Y.parent ) and ( Y <> X ) ) then
            Error( "<X> and <Y> must have the same parent xmod" );
        fi;
    od;
    return X;
end;

##############################################################################
##
#F  XModOps.Size                                 [ Size(Source), Size(Range) ]
##

XModOps.Size := function( X )
    return [ Size( X.source ), Size( X.range ) ];
end;

##############################################################################
##
#F  XModOps.Elements                     [ Elements(Source), Elements(Range) ]
##

XModOps.Elements := function( X )

    return [ Elements( X.source ), Elements( X.range ) ];
end;

##############################################################################
##
#F  XModOps.Kernel                 the kernel of X is the kernel of X.boundary
##

XModOps.Kernel := function( X )

    return Kernel( X.boundary );
end;

##############################################################################
##
#F  XModOps.IsAspherical                         tests that X.boundary is mono
##

XModOps.IsAspherical := function( X )

    return IsMonomorphism( X.boundary );
end;

##############################################################################
##
#F  XModOps.IsSimplyConnected                    tests that X.boundary is onto
##

XModOps.IsSimplyConnected := function( X )

    local  ok;

    ok := IsEpimorphism( X.boundary );
    if ok then
        X.isSimplyConnected := true;
    fi;
    return ok;
end;

##############################################################################
##
#F  XModOps.IsConjugation                      tests for normal inclusion xmod
##

XModOps.IsConjugation := function( X )

    local  Xsrc, Xrng, Xbdy, Xact, Xaut, s, sa, sb, sc, si,
           gensrc, genim, autgen, alpha, r, conj, imbdy, inv;

    Xsrc := X.source;
    Xrng := X.range;
    Xbdy := X.boundary;
    Xact := X.action;
    Xaut := X.aut;
    # the boundary map should be an inclusion
    if not IsMonomorphism( X.boundary ) then
        return false;
    fi;
    if ( Xsrc.generators = [ ] ) then
        return true;
    fi;
    gensrc := Xsrc.generators;
    genim := List( gensrc, s -> Image( Xbdy, s ) );
    imbdy := Subgroup( Xrng, genim );
    inv := GroupHomomorphismByImages( imbdy, Xsrc, genim, gensrc );

    for r in Xrng.generators do
        alpha := Image( Xact, r );
        conj := InnerAutomorphism( imbdy, r );
        for s in gensrc do
            sa := Image( alpha, s );
            sb := Image( Xbdy, s );
            sc := Image( conj, sb );
            si := Image( inv, sc );
            if ( sa <> si ) then
                return false;
            fi;
        od;
    od;
    return true;
end;

##############################################################################
##
#F  XModOps.IsTrivialActionXMod         Image(X.action) in X.aut is identity ?
##

XModOps.IsTrivialAction := function( X )

    local  act, zero;

    act:= X.action;
    zero := ZeroMorphism( X.range, X.aut );
    if ( act <> zero ) then
        return false;
    fi;
    return true;
end;

###############################################################################
##
#F  XModOps.IsCentralExtension           X.bdy surjective with central kernel ?
##

XModOps.IsCentralExtension := function( X )

    local  bdy, ZS, ker, S, genS, r, s, t, ta, tc, act, alpha, conj, cgen;

    bdy := X.boundary;
    if not IsSurjective( bdy ) then
        return false;
    fi;
    ker := Kernel( bdy );
    S := X.source;
    act := X.action;
    genS := S.generators;
    ZS := Centre( S );
    if not IsSubgroup( ZS, ker ) then
        return false;
    fi;
    # now test the action
    for s in genS do
        r := Image( bdy, s );
        alpha := Image( act, r );
        cgen := List( genS, t -> t^s );
        conj := GroupHomomorphismByImages( S, S, genS, cgen );
        for t in genS do
            ta := Image( alpha, t );
            tc := Image( conj, t );
            if ( ta <> tc ) then
                return false;
            fi;
        od;
    od;
    return true;
end;

##############################################################################
##
#F  XModOps.DirectProduct              (bdy1 x bdy2) : (S1 x S2) --> (R1 x R2)
##

XModOps.DirectProduct := function( X, Y )

    local  Xsrc, Xrng, Ysrc, Yrng, genXrng, lenXrng, genYrng, lenYrng, j,
           XY, S, R, genS, genR, lenR, aut, autgen, act, imbdy, bdy, name,
           Xbdy, Ybdy, imXbdy, imYbdy, permX, permY, emX, emY,
           genXsrc, genYsrc, a, ema, ima, b, r;

    if (    not ( IsBound( X.isXMod ) and X.isXMod )
         or not ( IsBound( Y.isXMod ) and Y.isXMod ) ) then
        Error( "Both parameters must be crossed modules" );
    fi;
    Xsrc := X.source;
    Ysrc := Y.source;
    genXsrc := Xsrc.generators;
    genYsrc := Ysrc.generators;
    S := DirectProduct( Xsrc, Ysrc );
    if ( IsBound( Xsrc.name ) and IsBound( Ysrc.name ) ) then
        S.name := Concatenation( Xsrc.name, "x", Ysrc.name );
    fi;
    Xrng := X.range;
    Yrng := Y.range;
    genXrng := Xrng.generators;
    genYrng := Yrng.generators;
    lenXrng := Length( genXrng );
    lenYrng := Length( genYrng );
    R := DirectProduct( Xrng, Yrng );
    if ( IsBound( Xrng.name ) and IsBound( Yrng.name ) ) then
        R.name := Concatenation( Xrng.name, "x", Yrng.name );
    fi;
    if not ( IsBound( R.perms ) and IsBound( S.perms ) ) then
        Error( "Direct products not permutation groups" );
    fi;
    genS := S.generators;
    genR := R.generators;
    lenR := Length( genR );
    name := Concatenation( X.name, "x", Y.name );

    Xbdy := X.boundary;
    Ybdy := Y.boundary;
    imXbdy := Xbdy.genimages;
    imYbdy := Ybdy.genimages;
    permX := R.perms[1];
    permY := R.perms[2];
    emX := List( imXbdy, r -> r^permX );
    emY := List( imYbdy, r -> r^permY );
    imbdy := Concatenation( emX, emY );
    if ( XModPrintLevel > 1 ) then
        Print( Concatenation( imXbdy, imYbdy ), " -> \n", imbdy, "\n" );
    fi;
    bdy := GroupHomomorphismByImages( S, R, genS, imbdy );
    permX := S.perms[1];
    permY := S.perms[2];
    emX := List( genXsrc, s -> s^permX );
    emY := List( genYsrc, s -> s^permY );
    autgen := 0 * [1..lenR ];
    for j in [1..lenXrng] do
        r := genXrng[j];
        a := Image( X.action, r );
        ima := List( genXsrc, s -> Image( a, s ) );
        ema := Concatenation( List( ima, s -> s^permX ), emY );
        b := GroupHomomorphismByImages( S, S, genS, ema );
        autgen[j] := b;
    od;
    for j in [ (lenXrng+1)..lenR] do
        r := genYrng[ j - lenXrng ];
        a := Image( Y.action, r );
        ima := List( genYsrc, s -> Image( a, s ) );
        ema := Concatenation( emX, List( ima, s -> s^permY ) );
        b := GroupHomomorphismByImages( S, S, genS, ema );
        autgen[j] := b;
    od;
    aut := Group( autgen, InclusionMorphism( S, S ) );
    act := GroupHomomorphismByImages( R, aut, genR, autgen );
    XY := XMod( bdy, act );
    XY.isDirectProduct := true;
    return XY;
end;

###############################################################################
##
#F  XModOps.IsAutomorphismXMod       test X of form: G -> A \in [Inn(G),Aut(G)]
##

XModOps.IsAutomorphismXMod := function( X )

    local  Xsrc, Xrng, IA, a, aut;

    Xsrc := X.source;
    Xrng := X.range;
    if ( IsBound( Xsrc.automorphismXMod ) and
                ( Xsrc.automorphismXMod = Xrng ) ) then
        return true;
    fi;
    aut := X.aut;
    if not ( IsBound( aut.automorphismPair ) and
           ( aut.automorphismPair = Xrng ) ) then
        return false;
    fi;
    IA := InnerAutomorphismGroup( Xsrc );
    for a in IA.generators do
        if not ( a in aut ) then
            return false;
        fi;
    od;
    return true;
end;

###############################################################################
##
#F  XModOps.IsZeroBoundary
##

XModOps.IsZeroBoundary := function( X )

    local  bdy, zero;

    bdy := X.boundary;
    zero := ZeroMorphism( X.source, X.range );
    if ( bdy <> zero ) then
        return false;
    fi;
    return true;
end;

###############################################################################
##
#F  XModOps.IsRModule
##

XModOps.IsRModule := function( X )

    local  ok, S;

    S := X.source;
    ok := ( IsZeroBoundary( X ) and IsAbelian( S ) );
    return ok;
end;

##############################################################################
##
#F  WhatTypeXMod           tests X to see if one of the standard constructions
##

WhatTypeXMod := function( X )

    local  ok, okt, okz, okc, oke, oka, okm, L;

    if not IsXMod( X ) then
        Error( "parameter must be a crossed module" );
    fi;
    L := [ ];
    okt := IsTrivialAction( X );
    if okt then
        Add( L, " triv, " );
    fi;
    okz := IsZeroBoundary( X );
    if okz then
        Add( L, " zero, " );
    fi;
    okc := IsConjugation( X );
    if okc then
        Add( L, " conj, " );
    fi;
    oke := IsCentralExtension( X );
    if oke then
        Add( L, " extn, " );
    fi;
    oka := IsAutomorphismXMod( X );
    if oka then
        Add( L, " auto, " );
    fi;
    okm := IsRModule( X );
    if okm then
        Add( L, " RMod, " );
    fi;
    ok := ( okt or okz or okc or oke or oka or okm );
    if not ok then
        Add( L, " undetermined " );
    fi;
    return L;
end;

###############################################################################
##
#F  XModOps.InclusionMorphism                construct crossed module inclusion
##

XModOps.InclusionMorphism := function( S, X )

    local phi, phisrc, phirng, Xsrc, Ssrc, Xrng, Srng;

    if ( not IsXMod( X ) ) or ( not IsXMod( S ) ) then
        Error( "Both arguments must be crossed modules" );
        return false;
    fi;
    if not IsSubXMod( X, S ) then
        Error( "S must be a sub-crossed module of X in IsSubXMod(S,X)" );
    fi;
    Xsrc := X.source;
    Ssrc := S.source;
    Xrng := X.range;
    Srng := S.range;

    phisrc := InclusionMorphism( Ssrc, Xsrc );
    phirng := InclusionMorphism( Srng, Xrng );

    phi := rec( );
    phi.sourceHom := phisrc;
    phi.rangeHom := phirng;
    phi.source := S;
    phi.range := X;
    phi.name := XModMorphismName( phi );
    phi.isXModMorphism := true;
    phi.operations := XModMorphismOps;
    phi.domain := Mappings;
    return phi;
end;

##############################################################################
##
#F  XModOps.WhiteheadPermGroup                    group of regular derivations
##

XModOps.WhiteheadPermGroup := function( X )

    local  L, M, N, i, j, k, D, reg, e, f, g, elts, Arec,
           gens, genpos, used, found, pos, W;

    if not ( IsXMod ( X ) ) then
        Error( " Input parameter must be a crossed module" );
    fi;
    D := RegularDerivations( X );
    if ( IsBound( X.actorSquare ) and
         IsBound( X.actorSquare.WhiteheadPermGroup ) ) then
        return X.actorSquare.WhiteheadPermGroup;
    fi;
    reg := D.regular;
    M := WhiteheadGroupTable( X );
    N := [ ];
    for k in [1..reg] do
        L := [1..reg];
        for j in [1..reg] do
            L[j] := M[j][k];
        od;
        Add( N, L );
    od;
    elts := List( N, x -> PermList( x ) );
    gens := [ ];
    genpos := [ ];
    used := 0 * [1..reg];
    used[1] := 1;
    for i in [1..reg] do
        if ( used[i] = 0 ) then
            e := elts[i];
            Add( gens, e );
            Add( genpos, i );
            used[i] := 1;
            found := Concatenation( [ ( ) ], gens );
            for e in found do
                for g in gens do
                    f := e*g;
                    if not ( f in found ) then
                        Add( found, f );
                        k := Position( elts, f );
                        used[k] := 1;
                    fi;
                od;
            od;
        fi;
    od;
    # find which derivations are generators
    if ( Length( genpos ) = 0 ) then
        Print( "D.genpos is empty!  Adding the identity!\n" );
        genpos := [ 1 ];
    fi;
    W := Group( gens, () );
    W.elements := elts;
    if IsBound( X.name ) then
        W.name := Concatenation( "WG(", X.name, ")" );
    fi;
    Arec := ActorSquareRecord( X );
    Arec.WhiteheadPermGroup := W;
    D.genpos := genpos;
    return W;
end;

##############################################################################
##
#F  XModOps.Whitehead            (InnerSourceHom : X.range -> Whitehead Group)
##

XModOps.Whitehead := function( X )

    local  Xsrc, gensrc, imeta, dersrc, D, L, psrc, W, genW, ngW, elW, eta,
           genpos, imchi, chi, sigma, j, s, autS, invgen, imact0, act0,
           imact, act, name, WX;

    Xsrc := X.source;
    D := RegularDerivations( X );
    L := D.genimageList;
    W := WhiteheadPermGroup( X );
    genW := W.generators;
    genpos := D.genpos;
    ngW := Length( genpos );
    elW := Elements( W );
    gensrc := Xsrc.generators;
    # determine the boundary map = innerSourceHom
    imeta := [ ];
    for s in gensrc do
        chi := InnerDerivation( X, s );
        j := Position( L, chi.genimages );
        Add( imeta, elW[j] );
    od;
    eta := GroupHomomorphismByImages( Xsrc, W, gensrc, imeta );
    if not IsHomomorphism( eta ) then
        Error( "Whitehead boundary fails to be a homomorphism" );
    fi;
    # now calculate the action homomorphism
    autS := AutomorphismGroup( X.source );
    imact0 := 0 * [1..ngW];
    for j in [1..ngW] do
        imchi := L[ genpos[j] ];
        sigma := SourceEndomorphismDerivation( X, imchi );
        imact0[j] := sigma; 
    od;
    invgen := List( genW, g -> g^-1 );
    act0 := GroupHomomorphismByImages( W, autS, invgen, imact0 );
    imact := List( genW, w -> Image( act0, w ) );
    act := GroupHomomorphismByImages( W, autS, genW, imact );
    WX := XMod( eta, act );
    name := XModName( X );
    WX.name := Concatenation( "Whitehead", name );
    WX.isWhitehead := true;
    return WX;
end;

##############################################################################
##
#F  XModOps.Norrie                         crossed module  ( Xrng --> Aut(X) )
##

XModOps.Norrie := function( X )

    local  Xsrc, Xrng, genrng, inn, im, PX, genP, bdy, ok, r, autr, conjr,
           psrc, prng, emsrc, emrng, a2psrc, a2prng, p2arng, act, f, NX,
           Arng, act, imrng, imact, autrng, i, p, projp, proja, ima, a, name;

    Xsrc := X.source;
    Xrng := X.range;
    genrng := Xrng.generators;
    PX := AutomorphismPermGroup( X );
    genP := PX.generators;
    autrng := [ 1..Length( genP ) ];
    a2psrc := Xsrc.automorphismGroup.automorphismPair.a2p;
    a2prng := Xrng.automorphismGroup.automorphismPair.a2p;
    p2arng := Xrng.automorphismGroup.automorphismPair.p2a;
    Arng := a2prng.source;
    # determine the boundary map
    im := [ ];
    for r in genrng do
        autr := Image( X.action, r );
        psrc := Image( a2psrc, autr );
        emsrc := Image( PX.embedSourceAuto, psrc );
        conjr := InnerAutomorphism( Xrng, r );
        prng := Image( a2prng, conjr );
        emrng := Image( PX.embedRangeAuto, prng );
        Add( im, emrng * emsrc );
    od;
    bdy := GroupHomomorphismByImages( Xrng, PX, genrng, im );
    ok := IsHomomorphism( bdy );
    if not ok then
        Error( "<bdy> is not a group homomorphism" );
    fi;

    # determine the action
    imact := 0 * autrng;
    for i in autrng do
        p := genP[i];
        projp := Image( PX.projrng, p );
        proja := Image( p2arng, projp );
        ima := List( genrng, r -> Image( proja, r ) );
        a := GroupHomomorphismByImages( Xrng, Xrng, genrng, ima );
        imact[i] := a;
    od;
    act := GroupHomomorphismByImages( PX, Arng, genP, imact );
    for f in act.genimages do
        if ( f = IdentityMapping( Xrng ) ) then
            f := InclusionMorphism( Xrng, Xrng );
        fi;
    od;

    NX := XMod( bdy, act );
    name := XModName( X );
    NX.name := Concatenation( "Norrie", name );
    NX.isNorrie := true;
    return NX;
end;

##############################################################################
##
#F  XModOps.Lue                            crossed module  ( Xrng --> Aut(X) )
##

XModOps.Lue := function( X )

    local  Xbdy, Nbdy, NX, Lbdy, LX, PX, Xsrc, gensrc, genP, act, f, Asrc,
           p2asrc, imsrc, imact, autrng, i, p, projp, proja, ima, a, name;

    NX := Norrie( X );
    Nbdy := NX.boundary;
    Xbdy := X.boundary;
    Lbdy := Xbdy * Nbdy;

    PX := AutomorphismPermGroup( X );
    genP := PX.generators;
    autrng := [ 1..Length( genP ) ];
    Xsrc := X.source;
    gensrc := Xsrc.generators;
    p2asrc := Xsrc.automorphismGroup.automorphismPair.p2a;
    Asrc := p2asrc.range;
    imact := 0 * autrng;
    for i in autrng  do
        p := genP[i];
        projp := Image( PX.projsrc, p );
        proja := Image( p2asrc, projp );
        ima := List( gensrc, s -> Image( proja, s ) );
        a := GroupHomomorphismByImages( Xsrc, Xsrc, gensrc, ima );
        imact[i] := a;
    od;
    act := GroupHomomorphismByImages( PX, Asrc, genP, imact );
    for f in act.genimages do
        if ( f = IdentityMapping( Xsrc ) ) then
            f := InclusionMorphism( Xsrc, Xsrc );
        fi;
    od;
    LX := XMod( Lbdy, act );
    name := XModName( X );
    LX.name := Concatenation( "Lue", name );
    LX.isLue := true;
    return LX;
end;

###############################################################################
##
#F  XModOps.Actor
##

XModOps.Actor := function( X )

    local  W, eW, PX, D, genpos, ngW, genW, invW, imbdy0, bdy0, imbdy, bdy,
           i, j, mor, a2psrc, a2prng, emsrc, emrng, imsrc, imrng, ActX, L,
           chi, chj, alpha, impos, imgen, imact, act, aut, mor, phi, id, name;

    D := RegularDerivations( X );
    L := D.genimageList;
    W := WhiteheadPermGroup( X );
    eW := Elements( W );
    PX := AutomorphismPermGroup( X );
    genpos := D.genpos;
    ngW := Length( genpos );
    # determine the boundary map
    genW := List( genpos, i -> eW[i] );
    invW := List( genW, g -> g^-1 );
    imbdy0 := 0 * genpos;
    a2psrc := X.source.automorphismGroup.automorphismPair.a2p;
    a2prng := X.range.automorphismGroup.automorphismPair.a2p;
    emsrc := PX.embedSourceAuto;
    emrng := PX.embedRangeAuto;
    for i in [1..ngW] do
        j := genpos[i];
        mor := XModEndomorphismDerivation( X, L[j] );
        imsrc := Image( emsrc, Image( a2psrc, mor.sourceHom ) );
        imrng := Image( emrng, Image( a2prng, mor.rangeHom ) );
        imbdy0[i] := imsrc * imrng;
    od;
    bdy0 := GroupHomomorphismByImages( W, PX, invW, imbdy0 );
    imbdy := List( genW, w -> Image( bdy0, w ) );
    bdy := GroupHomomorphismByImages( W, PX, genW, imbdy );

    # determine the action
    imact := [ ];
    for alpha in PX.autogens do
        mor := XModMorphism( X, X, alpha );
        impos := 0 * [1..ngW];
        for i in [1..ngW] do
            j := genpos[i];
            chi := XModDerivationByImages( X, L[j], true );
            chj := ImageAutomorphismDerivation( mor, chi );
            impos[i] := Position( L, chj.genimages );
        od;
        imgen := List( impos, i -> eW[i] );
        phi := GroupHomomorphismByImages( W, W, genW, imgen );
        Add( imact, phi );
    od;
    id := InclusionMorphism( W, W );
    aut := Group( imact, id );
    aut.name := "Aut(W)";
    act := GroupHomomorphismByImages( PX, aut, PX.generators, imact );
    ActX := XMod( bdy, act );
    name := XModName( X );
    ActX.name := Concatenation( "Actor", name );
    ActX.isActor := true;
    return ActX;
end;

##############################################################################
##
#F  XModOps.InnerMorphism                           maps X into its actor xmod
##

XModOps.InnerMorphism := function( X )

    local  A, WX, NX, eta, theta, inn;

    WX := Whitehead( X );
    eta := WX.boundary;
    NX := Norrie( X );
    theta := NX.boundary;
    A := Actor( X );
    inn := XModMorphism( X, A, [eta,theta] );
    if not IsXModMorphism( inn ) then
        Print( "Inner morphism fails to be a morphism!\n" );
        inn := false;
    fi;
    inn.isInnerMorphism := true;
    return inn;
end;

##############################################################################
##
#F  XModOps.Centre              the kernel of the inner morphism  X -> X.actor
##

XModOps.Centre := function( X )

    local K, name;

    K := Kernel( InnerMorphism( X ) );
    name := XModName( X );
    K.name := Concatenation( "Centre", name );
    K.isCentre := true;
    return K;
end;

##############################################################################
##
#F  XModOps.InnerActor           the image of the inner morphism  X -> X.actor
##

XModOps.InnerActor := function( X )

    local  InnX, mor, name, ActX;
      
    mor := InnerMorphism( X );
    ActX := mor.range;
    InnX := ImageXModMorphism( mor, X );
    if ( InnX = ActX ) then
        InnX := ActX;
    else
        name := XModName( X );
        InnX.name := Concatenation( "InnerActor", name );
    fi;
    InnX.isInnerActor := true;
    return InnX;
end;

##############################################################################
##
#F  AutomorphismPermGroup            subgroup of  Aut(X.source) x Aut(X.range)
##

XModOps.AutomorphismPermGroup := function( X )

    local  Xsrc, gensrc, Xrng, genrng, Xbdy, Asrc, Psrc, Arng, Prng, genA,
           D, genD, P, genP, p2, p2asrc, p2arng, Esrc, Erng, imsrc, imrng,
           phi0, imphi, phi, psi0, impsi, psi, num, es, er, e, mor, ismor,
           genPsrc, genPrng, projsrc, projrng, egensrc, egenrng,
           embedPsrc, embedPrng, newsrc, newrng, Arec, srcpair, rngpair;

    Xsrc := X.source;
    gensrc := Xsrc.generators;
    Xrng := X.range;
    genrng := Xrng.generators;
    Xbdy := X.boundary;
    Asrc := AutomorphismGroup( Xsrc );
    Asrc.isAutomorphismGroup := true;
    srcpair := AutomorphismPair( Asrc );
    Psrc := srcpair.perm;
    p2asrc := srcpair.p2a;
    Esrc := Elements( Psrc );
    Arng := AutomorphismGroup( Xrng );
    Arng.isAutomorphismGroup := true;
    rngpair := AutomorphismPair( Arng );
    Prng := rngpair.perm;
    p2arng := rngpair.p2a;
    Erng := Elements( Prng );
    D := DirectProduct( Psrc, Prng );
    genD := D.generators;
    P := Subgroup( D, [ ] );
    p2 := D.perms[2];
    newsrc := D.news[1];
    newrng := D.news[2];
    if ( IsBound( Psrc.name ) and IsBound( Prng.name ) ) then
        D.name := Concatenation( Psrc.name,"x", Prng.name );
    fi;
    num := 0;
    genA := [ ];
    for es in Esrc do
        phi0 := Image( p2asrc, es );
        imphi := List( gensrc, s -> Image( phi0, s ) );
        phi := GroupHomomorphismByImages( Xsrc, Xsrc, gensrc, imphi );
        for er in Erng do
            psi0 := Image( p2arng, er );
            impsi := List( genrng, r -> Image( psi0, r ) );
            psi := GroupHomomorphismByImages( Xrng, Xrng, genrng, impsi );
            if ( Xbdy * psi  =  phi * Xbdy ) then
                num := num + 1;
                e := es * er^p2;
                if not ( e in P ) then
                    mor := XModMorphism( X, X, [phi,psi] );
                    ismor := IsXModMorphism( mor );
                        if ismor then
                        Add( genA, [ phi, psi ] );
                        P := Closure( P, e );
                    fi;
                fi;
            fi;
        od;
    od;
    genP := P.generators;
    imsrc := List( genP, g -> PermList( List( newsrc, x->x^g) ) );
    imrng := List( genP, g -> PermList( List( newrng, x->(x^g)^(p2^-1) ) ) );
    projsrc := GroupHomomorphismByImages( P, Psrc, genP, imsrc );
    projrng := GroupHomomorphismByImages( P, Prng, genP, imrng );
    genPsrc := Psrc.generators;
    embedPsrc := GroupHomomorphismByImages( Psrc, D, genPsrc, genPsrc );
    genPrng := Prng.generators;
    egenrng := List( genPrng, p -> p^p2 );
    embedPrng := GroupHomomorphismByImages( Prng, D, genPrng, egenrng );
    P.embedSourceAuto := embedPsrc;
    P.embedRangeAuto := embedPrng;
    P.projsrc := projsrc;
    P.projrng := projrng;
    P.autogens := genA;
    P.xmod := X;
    Arec := ActorSquareRecord( X );
    Arec.automorphismPermGroup := P;
    return P;
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

## end of file xmod.g
