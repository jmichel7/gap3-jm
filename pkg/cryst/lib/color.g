#############################################################################
##
#A  color.g                 CrystGap library                     Bettina Eick
#A                                                              Franz G"ahler
#A                                                              Werner Nickel
##
#Y  Copyright 1990-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  CrystGap - the crystallographic groups package for GAP (color groups)
##  

#############################################################################
##
#F  IsColorGroup( G ) . . . . . . . . . . . . . . . . . .is it a color group?
##
IsColorGroup := function( G )
    return IsBound( G.isColorGroup ) and G.isColorGroup;
end;


#############################################################################
##
#F  ColorSubgroup( G ) . . . . . . . . . . . . . . extract the color subgroup
##
ColorSubgroup := function( G )
    if not IsColorGroup( G ) then
        Error("G must be a color group");
    fi;
    return G.colorSubgroup;
end;


#############################################################################
##
#F  ColorCosets( G ) . . . . . . . . . . . . . . . . . . . . . . color cosets
##
ColorCosets := function( G )

    local c;

    # G must be a color group
    if not IsColorGroup( G ) then
        Error("G must be a color group");
    fi;

    # if not known, compute them
    if not IsBound( G.colorCosets ) then
        G.colorCosets := RightCosets( G, G.colorSubgroup );        
    fi;

    return G.colorCosets;

end;


#############################################################################
##
#F  ColorOfElement( G, elem ) . . . . . . . . . . . . . . color of an element
##
ColorOfElement := function( G, elem )

    local cos, i;

    cos := ColorCosets( Parent( G ) );
    for i in [1..Length( cos )] do
        if elem in cos[i] then
            return i;
        fi;
    od;
    Error("elem must be an element of G");

end;


#############################################################################
##
#F  ColorPermGroup( G ) . . . . . . . . . . . . . . . . . . . color PermGroup
##
ColorPermGroup := function( G )

    local P, C, gens, permgens, H;

    # G must be a color group
    if not IsColorGroup( G ) then
        Error("G must be a color group");
    fi;

    if not IsBound( G.colorPermGroup ) then

        P := Parent( G );
        G.colorPermGroup    := Operation( G, ColorCosets( P ), OnRight );
        G.colorHomomorphism := OperationHomomorphism( G, G.colorPermGroup );

        if IsCrystGroup( G ) then
            G.colorHomomorphism.operations := CrystGroupOpHomOps;
        fi;

    fi;

    return G.colorPermGroup;

end;


#############################################################################
##
#F  ColorHomomorphism( G ) . . . . . . . . . .homomorphism to color PermGroup
##
ColorHomomorphism := function( G )
    local P;
    if not IsBound( G.colorHomomorphism ) then
        P := ColorPermGroup( G );
    fi;
    return G.colorHomomorphism;
end;


#############################################################################
##
#F  ColoredSubgroup( G, gens ) . . . . . . . . . . . . . . . colored subgroup
##
ColoredSubgroup := function( G, gens )

    local P, U, V;

    P := Parent( G );

    # first get an uncolored subgroup
    U := P.operations.UncoloredSubgroup( P, gens );

    # compute color subgroup of U; 
    # we need a copy of U which does not have a colored parent
    V := Stabilizer( Group( U ), P.colorCosets[1], OnRight );
    U.colorSubgroup := P.operations.UncoloredSubgroup( P, V.generators ); 

    # then color U
    U.operations := P.operations;
    U.isColorGroup := true;

    return U;

end;


#############################################################################
##
#F  ColorPointGroupOps . . . . . . . . . . . . . . . . . . ColorPointGroupOps
##
ColorPointGroupOps                     := ShallowCopy( PointGroupOps );
ColorPointGroupOps.name                := "ColorPointGroupOps";
ColorPointGroupOps.UncoloredSubgroup   := ColorPointGroupOps.Subgroup;
ColorPointGroupOps.Subgroup            := ColoredSubgroup;


#############################################################################
##
#F  ColoredPointGroup( S ) . . . . . . . . . . . . . . . .colored point group
##
ColoredPointGroup := function( S )

    local d, P, H, U;

    # get first the uncolored point group
    P := S.operations.UncoloredPointGroup( S );
    S.pointGroup := P;

    # if S and its ColorSubgroup are lattice-equal, 
    # the point group can be colored
    if TranslationsCrystGroup( S ) = 
                 TranslationsCrystGroup( ColorSubgroup( S ) ) then
        d := S.dimension - 1;
        P.colorSubgroup := P.operations.Subgroup( Parent( P ), 
            List( S.colorSubgroup.generators, x -> x{[1..d]}{[1..d]} ) );
        P.colorCosets := List( ColorCosets( S ), 
            x -> RightCoset( P.colorSubgroup, 
                             x.representative{[1..d]}{[1..d]} ) );
        P.isColorGroup := true;
        P.operations := ColorPointGroupOps;
    fi;
    return P;

end;


#############################################################################
##
#F  ColorCrystGroupOps . . . . . . . . . . . . . . . . . . ColorCrystGroupOps
##
ColorCrystGroupOps                     := ShallowCopy( CrystGroupOps );
ColorCrystGroupOps.name                := "ColorCrystGroupOps";
ColorCrystGroupOps.UncoloredSubgroup   := ColorCrystGroupOps.Subgroup;
ColorCrystGroupOps.Subgroup            := ColoredSubgroup;
ColorCrystGroupOps.UncoloredPointGroup := ColorCrystGroupOps.PointGroup;
ColorCrystGroupOps.PointGroup          := ColoredPointGroup;


#############################################################################
##
#F  ColorGroup( G, H ) . . . . . . . . . . . . . . . . . . make a color group
##
ColorGroup := function( G, H )

    local NewG, C;

    # H must be a subgroup of G
    if not IsSubgroup( G, H ) then
        Error("H must be a subgroup of G");
    fi;

    # we allow only parent groups; subgroups inherit coloring from parent
    if not IsParent( G ) then
        Error("G must be a parent");
    fi;

    # since G may contain uncolored information components, make a new group
    NewG := Group( G.generators, G.identity );

    # if G is a CrystGroup, NewG should be one, too
    if IsCrystGroup( G ) then
        NewG := CrystGroup( NewG );
        AddTranslationsCrystGroup( NewG, TranslationsCrystGroup( G ) );
    fi;

    # make NewG a color group
    NewG.colorSubgroup := NewG.operations.Subgroup( NewG, H.generators );
    NewG.isColorGroup  := true;

    # to fix a numbering of the colors, we fix the color cosets
    C := ColorCosets( NewG );

    # set the operations record; color CrystGroups get ColorCrystGroupPps, 
    # other kinds of color groups their own unique operations record
    # (we better don't change the operations record of, say, all PermGroups)
    if IsCrystGroup( G ) then
        Unbind( NewG.pointGroup );
        Unbind( NewG.pointHomom );
        NewG.operations := ColorCrystGroupOps;
    else
        NewG.operations := ShallowCopy( NewG.operations );
        NewG.operations.name := 
                Concatenation( "Color", NewG.operations.name );
        NewG.operations.UncoloredSubgroup := NewG.operations.Subgroup;
        NewG.operations.Subgroup := ColoredSubgroup;
    fi;

    return NewG;

end;



