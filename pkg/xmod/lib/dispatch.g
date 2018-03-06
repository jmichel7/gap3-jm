###############################################################################
##
##  dispatch.g                  for GAP 3.4                    version  9/ 1/97
##
###############################################################################
##
#A  dispatch.g                   GAP library                      Chris Wensley
#A                                                                    Murat Alp
#Y  Copyright
##
##  This file contains functions that manipulate crossed modules,
##  pre-crossed modules, group graphs and associated constructions.
##
#H  $Log: dispatch.g,v $
#H  Revision 1.1  1997/03/27 13:35:47  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H

##############################################################################
##
#F  ZeroMorphism                   mapping a structure to the identity element

ZeroMorphism := function( S, D )

    if not ( IsDomain( S ) and IsDomain( D ) ) then
        Error( "both arguments must be domains" );
    elif not IsBound( D.identity ) then
        Error( "second argument must have an identity sub-object" );
    fi;
    if ( IsGroup( S ) and IsGroup( D ) ) then
        return GroupOpsZeroMorphism( S, D );
    fi; 
    if IsBound( D.operations.ZeroMorphism ) then
        return D.operations.ZeroMorphism( S, D );
    else
        Error( "<D>.operations does not have ZeroMorphism defined" );
    fi;
end;

##############################################################################
##
#F  InclusionMorphism                              inclusion of a substructure
##

InclusionMorphism := function( S, D )

    if not ( IsDomain( S ) and IsDomain( D ) ) then
        Error( "both arguments must be domains" );
    elif not ( Parent( S ) = Parent( D ) ) then
        Error( "the two arguments must have a common parent" );
    fi;
    if ( IsGroup( S ) and IsGroup( D ) ) then
        return GroupOpsInclusionMorphism( S, D );
    fi; 
    if IsBound( D.operations.InclusionMorphism ) then
        return D.operations.InclusionMorphism( S, D );
    else
        Error( "<D>.operations does not have InclusionMorphism defined" );
    fi;
end;

##############################################################################
##
#F  IsAutomorphismGroup                tests if automorphism group of a domain
##

IsAutomorphismGroup := function( D )

    local ok;

    if IsDomain( D ) then
        if IsBound( D.isAutomorphismGroup ) then
            return D.isAutomorphismGroup;
        elif IsGroup( D ) then
            return GroupOpsIsAutomorphismGroup( D );
        fi;
        ok := D.operations.IsAutomorphismGroup( D );
        D.isAutomorphismGroup := ok;
    else
        return false;
    fi;
    return ok;
end;

##############################################################################
##
#F  AutomorphismPermGroup       construct perm rep of a group of automorphisms
##

AutomorphismPermGroup := function( D )

    local  A, L, P;

    if not IsDomain( D ) then
        Error( "<D> must be a domain" );
    fi;
    if IsBound( D.automorphismPermGroup ) then
        return D.automorphismPermGroup;
    fi;
    if ( IsBound( D.actorSquare ) and
         IsBound( D.actorSquare.automorphismPermGroup ) ) then
        return D.actorSquare.automorphismPermGroup;
    fi;
    if IsGroup( D ) then
        P := GroupOpsAutomorphismPermGroup( D );
        D.automorphismPermGroup := P;
    else
        P := D.operations.AutomorphismPermGroup( D );
    fi;
    #  saved by operations record functions in  D.actorSquare
    #  D.automorphismPermGroup := P;
    if IsBound( D.name ) then
        P.name := Concatenation( "PermAut(", D.name, ")" );
    fi;
    return P;
end;

##############################################################################
##
#F  InnerAutomorphismGroup            group of inner automorphisms of a domain
##

InnerAutomorphismGroup := function( D )

    local  A, I;
    if not IsDomain( D ) then
        Error( "<D> must be a domain" );
    fi;
    if IsBound( D.innerAutomorphismGroup ) then
        return D.innerAutomorphismGroup;
    fi;
    if IsGroup( D ) then
        I := GroupOpsInnerAutomorphismGroup( D );
    else
        I := D.operations.InnerAutomorphismGroup( D );
    fi;
    if IsBound( D.automorphismGroup ) then
        A := D.automorphismGroup;
        I := AsSubgroup( A, I );
        I.parent := A;
    else
        I.isParent := true;
    fi;
    if IsBound( D.name ) then
        I.name := Concatenation( "Inn(", D.name, ")" );
    fi;
    D.innerAutomorphismGroup := I;
    return I;
end;

##############################################################################
##
#F  IsAspherical                                 tests that X.boundary is mono
##

IsAspherical := function( D )

    local ok;

    if IsDomain( D ) then
        if IsBound( D.isAspherical ) then
            return D.isAspherical;
        fi;
        ok := D.operations.IsAspherical( D );
        D.isAspherical := ok;
    else
        return false;
    fi;
    return ok;
end;

##############################################################################
##
#F  IsSimplyConnected                            tests that X.boundary is onto
##

IsSimplyConnected := function( D )

    local ok;
    if IsDomain( D ) then
        if IsBound( D.isSimplyConnected ) then
            return D.isSimplyConnected;
        fi;
        ok := D.operations.IsSimplyConnected( D );
        D.isSimplyConnected := ok;
    else
        return false;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsConjugation
##

IsConjugation := function( D )

    local ok;
    if IsDomain( D ) then
        if IsBound( D.isConjugation ) then
            return D.isConjugation;
        fi;
        ok := D.operations.IsConjugation( D );
        D.isConjugation := ok;
    else
        return false;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsTrivialAction
##

IsTrivialAction := function( D )

    local ok;
    if IsDomain( D ) then
        if IsBound( D.isTrivialAction ) then
            return D.isTrivialAction;
        fi;
        ok := D.operations.IsTrivialAction( D );
        D.isTrivialAction := ok;
    else
        return false;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsCentralExtension
##

IsCentralExtension := function( D )

    local ok;
    if IsDomain( D ) then
        if IsBound( D.isCentralExtension ) then
            return D.isCentralExtension;
        fi;
        ok := D.operations.IsCentralExtension( D );
        D.isCentralExtension := ok;
    else
        return false;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsAutomorphismXMod          test X has the form: G -> A \in [Inn(G),Aut(G)]
##

IsAutomorphismXMod := function( D )

    local ok;
    if IsDomain( D ) then
        if IsBound( D.isAutomorphismXMod ) then
            return D.isAutomorphismXMod;
        fi;
        ok := D.operations.IsAutomorphismXMod( D );
        D.isAutomorphismXMod := ok;
    else
        return false;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsZeroBoundary
##

IsZeroBoundary := function( D )

    local ok;
    if IsDomain( D ) then
        if IsBound( D.isZeroBoundary ) then
            return D.isZeroBoundary;
        fi;
        ok := D.operations.IsZeroBoundary( D );
        D.isZeroBoundary := ok;
    else
        return false;
    fi;
    return ok;
end;

###############################################################################
##
#F  IsRModule
##

IsRModule := function( D )

    local ok;

    if IsDomain( D ) then
        if IsBound( D.isRModule ) then
            return D.isRModule;
        fi;
        ok := D.operations.IsRModule( D );
        D.isRModule := ok;
    elif IsRec( D ) then
        ok := IsRModuleRecord( D );
        D.isRModule := ok;
    else
        return false;
    fi;
    return ok;
end;

##############################################################################
##  functions for derivations, sections and related constructions
##############################################################################

##############################################################################
##
#F  ActorSquareRecord                       sets up subrecord for xmod or cat1
##

ActorSquareRecord := function( D )

    local  Arec;
    if IsDomain( D ) then
        if IsBound( D.actorSquare ) then
            return D.actorSquare;
        fi;
    else
        Error( "<D> must be a domain" );
    fi;
    Arec := rec( );
    Arec.isActorSquare := true;
    if ( D.operations = XModOps ) then
        Arec.xmod := D;
    elif ( D.operations = Cat1Ops ) then
        Arec.cat1 := D;
    fi;
    D.actorSquare := Arec;
    return Arec;
end;

##############################################################################
##
#F  WhiteheadPermGroup                   group of regular derivations/sections

WhiteheadPermGroup := function( D )

    local  W, Arec;

    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.WhiteheadPermGroup ) then
            return Arec.WhiteheadPermGroup;
        fi;
        W := D.operations.WhiteheadPermGroup( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.WhiteheadPermGroup := W;
    return W;
end;

##############################################################################
##
#F  Whitehead                            construct Whitehead XMod or Cat1Group
##

Whitehead := function( D )

    local  Arec, W;
    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.Whitehead ) then
            return Arec.Whitehead;
        fi;
        W := D.operations.Whitehead( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.Whitehead := W;
    return W;
end;

##############################################################################
##
#F  Norrie                                  construct Norrie XMod or Cat1Group
##

Norrie := function( D )

    local  Arec, N;
    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.Norrie ) then
            return Arec.Norrie;
        fi;
        N := D.operations.Norrie( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.Norrie := N;
    return N;
end;

##############################################################################
##
#F  Lue                                        construct Lue XMod or Cat1Group
##

Lue := function( D )

    local  Arec, L;
    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.Lue ) then
            return Arec.Lue;
        fi;
        L := D.operations.Lue( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.Lue := L;
    return L;
end;

##############################################################################
##
#F  Actor                              construct actor of an XMod or Cat1Group
##

Actor := function( D )

    local  A, Arec;
    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.actor ) then
            return Arec.actor;
        fi;
        A := D.operations.Actor( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.actor := A;
    return A;
end;

##############################################################################
##
#F  InnerMorphism                      map an XMod or Cat1Group into its actor
##

InnerMorphism := function( D )

    local  Arec, mor;
    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.innerMorphism ) then
            return Arec.innerMorphism;
        fi;
        mor := D.operations.InnerMorphism( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.innerMorphism := mor;
    return mor;
end;

##############################################################################
##
#F  InnerActor                   the image of the inner morphism  D -> D.actor
##

InnerActor := function( D )

    local  A, Arec;
    if IsDomain( D ) then
        Arec := ActorSquareRecord( D );
        if IsBound( Arec.innerActor ) then
            return Arec.innerActor;
        fi;
        A := D.operations.InnerActor( D );
    else
        Error( "<D> must be a domain" );
    fi;
    Arec.innerActor := A;
    return A;
end;

## end of file  dispatch.g
