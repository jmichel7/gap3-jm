#############################################################################
##
#A  cryst.g                 CrystGap library                     Bettina Eick
#A                                                              Franz G"ahler
#A                                                              Werner Nickel
##
#Y  Copyright 1990-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  CrystGap - the crystallographic groups package for GAP (part 1)
##  

#############################################################################
##
## first some common utility routines, most of which deal with
## integral matrices
##
#############################################################################


#############################################################################
##
#F  RowEchelonForm  . . . . . . . . . . row echelon form of an integer matrix
##
RowEchelonForm := function( M )
    
    local a, i, j, k, m, n, r, Cleared;
    
    if M = [] then return []; fi;
    
    M := Copy( M );
    m := Length( M ); n := Length( M[1] );
    
    i := 1; j := 1;
    while i <= m and j <= n do
        k := i; while k <= m and M[k][j] = 0 do k := k+1; od;
        if k <= m then
            r := M[i]; M[i] := M[k]; M[k] := r;
            for k in [k+1..m] do
                a := AbsInt( M[k][j] );
                if a <> 0 and a < AbsInt( M[i][j] ) then
                    r := M[i]; M[i] := M[k]; M[k] := r;
                fi;
            od;
            if M[i][j] < 0 then M[i] := -1 * M[i]; fi;
            Cleared := true;
            for k in [i+1..m] do
                a := QuoInt(M[k][j],M[i][j]);
                if a <> 0 then  M[k] := M[k] - a * M[i]; fi;
                if M[k][j] <> 0 then Cleared := false; fi;
            od;
            if Cleared then i := i+1; j := j+1; fi;
        else 
            j := j+1;
        fi;
    od;
    return M{[1..i-1]};
end;


#############################################################################
##
#F  RowEchelonFormVector  . . . . . . . . . . . .row echelon form with vector 
##
RowEchelonFormVector := function( M, b )
    
    local a, i, j, k, m, n, r, Cleared;
    
    M := Copy( M );
    m := Length( M );
    if m = 0 then return M; fi;
    n := Length( M[1] );
    
    i := 1; j := 1;
    while i <= m and j <= n do
        k := i; while k <= m and M[k][j] = 0 do k := k+1; od;
        if k <= m then
            r := M[i]; M[i] := M[k]; M[k] := r;
            r := b[i]; b[i] := b[k]; b[k] := r;
            for k in [k+1..m] do
                a := AbsInt( M[k][j] );
                if a <> 0 and a < AbsInt( M[i][j] ) then
                    r := M[i]; M[i] := M[k]; M[k] := r;
                    r := b[i]; b[i] := b[k]; b[k] := r;
                fi;
            od;
            if M[i][j] < 0 then M[i] := -1 * M[i]; b[i] := -1 * b[i]; fi;
            Cleared := true;
            for k in [i+1..m] do
                a := QuoInt(M[k][j],M[i][j]);
                if a <> 0 then
                    M[k] := M[k] - a * M[i];
                    b[k] := b[k] - a * b[i];
                fi;
                if M[k][j] <> 0 then Cleared := false; fi;
            od;
            if Cleared then i := i+1; j := j+1; fi;
        else 
            j := j+1;
        fi;
    od;
    return M{[1..i-1]};
end;


#############################################################################
##
#F  RowEchelonFormT . . . . . . . row echelon form with transformation matrix
##
RowEchelonFormT := function( M, T )
    
    local a, i, j, k, m, n, r, Cleared;
    
    M := Copy( M );
    m := Length( M ); n := Length( M[1] );
    
    i := 1; j := 1;
    while i <= m and j <= n do
        k := i; while k <= m and M[k][j] = 0 do k := k+1; od;
        if k <= m then
            r := M[i]; M[i] := M[k]; M[k] := r;
            r := T[i]; T[i] := T[k]; T[k] := r;
            for k in [k+1..m] do
                a := AbsInt( M[k][j] );
                if a <> 0 and a < AbsInt( M[i][j] ) then
                    r := M[i]; M[i] := M[k]; M[k] := r;
                    r := T[i]; T[i] := T[k]; T[k] := r;    
                fi;
            od;
            if M[i][j] < 0 then 
                M[i] := -1 * M[i]; 
                T[i] := -1 * T[i]; 
            fi;
            Cleared := true;
            for k in [i+1..m] do
                a := QuoInt(M[k][j],M[i][j]);
                if a <> 0 then 
                    M[k] := M[k] - a * M[i]; 
                    T[k] := T[k] - a * T[i];
                fi;
                if M[k][j] <> 0 then Cleared := false; fi;
            od;
            if Cleared then i := i+1; j := j+1; fi;
        else 
            j := j+1;
        fi;
    od;
    return M{[1..i-1]};
end;


#############################################################################
##
#F  HermiteNormalForm . . . . . . .  Hermite normal form of an integer matrix
##
HermiteNormalForm := function( M )
    local   k,  h,  j;
    
    M := RowEchelonForm( M );
    
    for k in Reversed([1..Length(M)]) do
        h := 1; while M[k][h] = 0 do h := h+1; od;
        for j in [1..k-1] do
            M[j] := M[j] - Int( M[j][h]/M[k][h] ) * M[k];
        od;
    od;
    
    h := 1;
    for k in [1..Length(M)] do
        while M[k][h] = 0 do h := h+1; od;
        for j in [1..k-1] do
            M[j] := M[j] - Int( M[j][h]/M[k][h] ) * M[k];
            if M[j][h] < 0 then M[j] := M[j] + M[k]; fi;
        od;
    od;
    
    return M;
end;


#############################################################################
##
#F  FractionModOne  . . . . . . . . . . . . . . . . . . a fraction modulo one
##
FractionModOne := function( q )
    q := q - Int(q);
    if q < 0 then q := q+1; fi;
    return q;
end;


#############################################################################
##
## Now the operations record for CrystGroups
##
#############################################################################


#############################################################################
##
#F  CrystGroupOps  . . . . . . . . . . . . initially is a copy of MatGroupOps
##
CrystGroupOps := OperationsRecord( "CrystGroupOps", MatGroupOps );


#############################################################################
##
#F  IsCrystGroup( <S> )  . .  everything that has CrystGroupOps as op. record
##
IsCrystGroup := function ( S )
    return IsGroup( S ) and IsBound( S.isCrystGroup ) and S.isCrystGroup;
end;


#############################################################################
##
#F  PointGroup( <S> )  . . . . . . . .  point group of crystallographic group
##
PointGroup := function ( S )

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if not IsBound( S.pointGroup ) then
        S.pointGroup := S.operations.PointGroup( S );
    fi;
    return S.pointGroup;

end;


#############################################################################
##
#F  ConjugateInternalGenerators( S ) . . .conjugate int. gens with int. basis
##
ConjugateInternalGenerators := function ( S )

    local d, B, P, gens, conj, inv, i;

    d := S.dimension - 1;
    B := S.internalBasis;
    P := PointGroup( S );

    # conjugate the internal generators if necessary
    if B <> IdentityMat( d ) then

        # conjugate internal point group generators
        gens := P.internalGenerators;
        inv  := B^-1;
        for i in [1..Length(gens)] do        
            gens[i] := B*gens[i]*inv;
        od;

        # conjugate internal space group generators
        conj := IdentityMat( d+1 );
        conj{[1..d]}{[1..d]} := B;
        inv := conj^-1;
        gens := S.internalGenerators;
        for i in [1..Length(gens)] do        
            gens[i] := conj*gens[i]*inv;
        od;

    fi;

end;


#############################################################################
##
#F  ReducedTranslationBasis .  reduce basis of transl. lattice to normal form
##
ReducedTranslationBasis := function ( trans )

    local tmp, den;

    if trans = [] then
        return trans;
    fi;

    tmp   := Flat( trans );
    if ForAll( tmp, IsInt ) then
        trans := HermiteNormalForm( trans );
    else
        den := Lcm( List( tmp, x -> Denominator( x ) ) );
        trans := HermiteNormalForm( den*trans )/den;
    fi;
    return trans;

end;


#############################################################################
##
#F  AddInternalBasis( S ) . . . . . . . . . . . . . . . . .add internal basis
##
AddInternalBasis := function( S )

    local d, T, comp, i, j, k;

    d := S.dimension - 1;
    T := S.translations;
    if Length( T ) = d then
        S.internalBasis := T;
    elif Length( T ) = 0 then
        S.internalBasis := IdentityMat( d );
    else
        comp := NullMat( d - Length(T), d );
        i:=1; j:=1; k:=1;
        while i <= Length( T ) do
            while T[i][j] = 0 do
                comp[k][j] := 1;
                k := k+1; j:=j+1;
            od;
            i := i+1; j := j+1;
        od;
        while j <= d do
            comp[k][j] := 1;
            k := k+1; j:=j+1;
        od;            
        S.internalBasis := Concatenation( T, comp );
    fi;

end;


#############################################################################
##
#F  AddTranslationsCrystGroup( S, basis ) . .add basis of translation lattice
##
AddTranslationsCrystGroup := function ( S, basis )

    local d, P, H, I, i;

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;

    # if translations are already known, do nothing
    if not IsBound( S.translations ) then

        d := S.dimension - 1;
        S.translations := ReducedTranslationBasis( basis );
        AddInternalBasis( S );
        S.isStandard   := S.internalBasis = IdentityMat( d );
        S.isSpaceGroup := Length( S.translations ) = d;

        # conjugate internal generators if necessary
        P := PointGroup( S );    # binds internal generators if necessary
        ConjugateInternalGenerators( S );

    fi;

end;


#############################################################################
##
#F  FindTranslationsCrystGroup( S ) . . . . . . find basis of transl. lattice
##
FindTranslationsCrystGroup := function ( S )

    local d, P, I, Sgens, Pgens, trans, g, m, i, j, F, r, mat;

    d := S.dimension - 1;
    P := PointGroup( S );
    I := IdentityMat( d );
    Sgens := [];
    Pgens := [];
    trans := [];

    # first the obvious translations
    for g in S.generators do
        m := g{[1..d]}{[1..d]};
        if m = I then
            Add( trans, g[d+1]{[1..d]} );
        else
            Add( Sgens, g );
            Add( Pgens, m );
        fi;
    od;

    # the translations from multiple point group generators
    if Length( Set( Pgens ) ) < Length( Pgens ) then
        SortParallel( Pgens, Sgens );
        i := 1;
        while i <= Length( Pgens ) do
            j := i+1;
            while j <= Length( Pgens ) and Pgens[i] = Pgens[j] do
                Add( trans, Sgens[i][d+1]{[1..d]} - Sgens[j][d+1]{[1..d]} );
                j := j+1;
            od;
            i := j;
        od;
    fi;

    # then the translations from the point group relators
    F := FpGroup( P );
    for r in F.relators do
        mat := MappedWord( r, F.generators, S.internalGenerators );
        Add( trans, mat[d+1]{[1..d]} );
    od;

    # make translations invariant under point group
    trans := Set( Union( Orbits( P, trans ) ) );
    return ReducedTranslationBasis( trans );

end;


#############################################################################
##
#F  TranslationsCrystGroup( <S> )  . . basis of transl. latt. of cryst. group
##
TranslationsCrystGroup := function ( S )

    local T;

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;

    if not IsBound( S.translations ) then
        T := FindTranslationsCrystGroup( S );
        AddTranslationsCrystGroup( S, T );
    fi;

    return S.translations;

end;


#############################################################################
##
#F  CheckTranslations( <S> ) . . . . . . . . . . . . . . . .CheckTranslations
##
CheckTranslations := function ( S )

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;

    if IsBound( S.translations ) then 
        if S.translations <> FindTranslationsCrystGroup( S ) then   
            Print("#W  Warning: Translations are INCORRECT - you better\n",
                  "#W           start again with a fresh group!\n"); 
        fi;
    fi;

end;


#############################################################################
##
#F  IsSpaceGroup( <S> )  . . . . . . . . . . . . . . . .  is S a space group?
##
IsSpaceGroup := function ( S )
    local T;
    if not IsCrystGroup( S ) then
        return false;
    elif not IsBound( S.isSpaceGroup ) then
        T := TranslationsCrystGroup( S );
        S.isSpaceGroup := Length( T ) = S.dimension - 1;
    fi;
    return S.isSpaceGroup;
end;


#############################################################################
##
#F  IsStandardCrystGroup( <S> )  . . . . . . . . . . . is S in standard form?
##
IsStandardCrystGroup := function ( S )
    local T;
    if not IsCrystGroup( S ) then
        return false;
    elif not IsBound( S.isStandard ) then
        T := TranslationsCrystGroup( S );
        S.isStandard := S.internalBasis = IdentityMat( S.dimension - 1 );
    fi;
    return S.isStandard;
end;


#############################################################################
##
#F  IsStandardSpaceGroup( <S> )  . . . . . . . .space group in standard form?
##
IsStandardSpaceGroup := function( S )
    return IsStandardCrystGroup( S ) and IsSpaceGroup( S );
end; 


#############################################################################
#F  IsSymmorphicSpaceGroup( S ) . . . . . . . . . . . . . . .is S symmorphic? 
##
IsSymmorphicSpaceGroup := function( S )
    if not IsSpaceGroup( S ) then
        return false;
    fi;
    if not IsBound( S.isSymmorphic ) then
        S.isSymmorphic := S.operations.IsSymmorphic( S );
    fi;
    return S.isSymmorphic;
end;


#############################################################################
##
#F  CrystGroup  . .  construct cryst group with appropriate operations record
##
## CrystGroup accepts as argument either a matrix group or a list of 
## generating matrices, or an argument identifying a group from the 
## crystallographic groups library.
## 
CrystGroup := function( arg )

    local S, d, P, lastcolumn, ok, m;

    # first check the arguments
    if IsMatGroup( arg[1] ) then
        S := Copy( arg[1] );
    elif IsMat(arg[1]) then
        S := Group( arg, IdentityMat( Length( arg[1] ) ) );
    elif Length( arg ) = 2 and IsList( arg[1] ) and IsMat( arg[2] ) then
        S := Group( arg[1], arg[2] );
    elif IsInt( arg[1] ) then
        if Length( arg ) = 2 then
            S := TransposedSpaceGroup( arg[1], arg[2] );
        elif Length( arg ) = 5 then
            S := TransposedSpaceGroup( arg[1], arg[2], arg[3],
                                               arg[4], arg[5] );
        else
            Error("wrong number of arguments");
        fi;
    elif IsString( arg[1] ) then
        S := TransposedSpaceGroup( arg[1] );
    else
        Error("wrong argument type");
    fi;

    d := S.dimension - 1;

    # make S a crystallograhic group
    S.operations := CrystGroupOps;
    S.isCrystGroup := true;

    # if S is from the library, we know a few things
    if IsBound( S.crTransposedSpaceGroupTyp ) then

        AddTranslationsCrystGroup( S, IdentityMat( d ) );
        S.name := Concatenation( "Cryst", S.name{[16..Length(S.name)]} );
        Unbind( S.fpGroup );

    elif IsBound( S.crSpaceGroupTyp ) then

        Error("use TransposedSpaceGroup for space groups from library");

    else

        # do a few basic checks
        ok := true;
        lastcolumn := List( [1..d], x -> 0 );
        Add( lastcolumn, 1 );
        for m in S.generators do
            ok := ok and m{[1..d+1]}[d+1] = lastcolumn;
        od;
        if not ok then
            Error("this is not a crystallographic group");
        fi;

    fi;

    return S;

end;


#############################################################################
##
#F  WyckoffPositionOps . . . . . . . . . . . . . . . . . . WyckoffPositionOps
##
WyckoffPositionOps := OperationsRecord( "WyckoffPositionOps" );


#############################################################################
##
#F  WyckoffPositionOps.Print . . . . . . . . . . . . WyckoffPositionOps.Print
##
WyckoffPositionOps.Print := function( w ) 
    Print( "[ class: ", w.class, ", translation: ",w.translation, 
           ", basis: ", w.basis, " ]" );
end;


#############################################################################
##
#F  IsWyckoffPosition . . . . . . . . . . . . . . . . . . . IsWyckoffPosition
##
IsWyckoffPosition := function( w )
    if not IsRec( w ) or not IsBound( w.isWyckoffPosition) then
        return false; 
    else
        return w.isWyckoffPosition;
    fi;
end;


#############################################################################
##
#F  WyckoffSpaceGroup . . . . . . . . . . . . . spaceGroup of WyckoffPosition
##
WyckoffSpaceGroup := function( w ) 
    if not IsWyckoffPosition( w ) then
        Error("w must be a Wyckoff position");
    else
        return w.spaceGroup; 
    fi;
end;


#############################################################################
##
#F  WyckoffTranslation . . . . . . . . . . . . translation of WyckoffPosition
##
WyckoffTranslation := function( w ) 
    if not IsWyckoffPosition( w ) then
        Error("w must be a Wyckoff position");
    else
        return w.translation; 
    fi;
end;


#############################################################################
##
#F  WyckoffBasis . . . . . . . . . . . . . . . . . . Basis of WyckoffPosition
##
WyckoffBasis := function( w ) 
    if not IsWyckoffPosition( w ) then
        Error("w must be a Wyckoff position");
    else
        return w.basis; 
    fi;
end;


#############################################################################
##
#F  WyckoffStabilizer . . . . . . . . . . . . . Stabilizer of WyckoffPosition
##
WyckoffStabilizer := function( w ) 
    if not IsWyckoffPosition( w ) then
        Error("w must be a Wyckoff position");
    else
        return w.stabilizer; 
    fi;
end;


#############################################################################
##
#F  WyckoffPosClass . . . . . . . . . . . . . class number of WyckoffPosition
##
WyckoffPosClass := function( w ) 
    if not IsWyckoffPosition( w ) then
        Error("w must be a Wyckoff position");
    else
        return w.class; 
    fi;
end;


#############################################################################
##
#F  ConjugatedCrystGroup( <S>, <conj> ) . change basis of translation lattice
##
ConjugatedCrystGroup := function ( S, conj )

    local d, c, ci, C, Ci, gens, i, gen, im, I, g, m, G, P, R, w;

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if not IsParent( S ) then
        Error("only parent groups should be conjugated");
    fi;

    # get the conjugators;
    d := S.dimension - 1;
    if Length( conj ) = d then
        c  :=conj;
        ci :=conj^-1;
        C  :=IdentityMat( d+1 ); C {[1..d]}{[1..d]} := c;
        Ci :=IdentityMat( d+1 ); Ci{[1..d]}{[1..d]} := ci;
    else
        C  := conj;
        Ci := conj^-1;
        c  := C {[1..d]}{[1..d]};
        ci := Ci{[1..d]}{[1..d]};
    fi;

    # conjugate generators of S
    gens := ShallowCopy( S.generators );
    for i in [1..Length(gens)] do
        gens[i] := C*gens[i]*Ci;
    od;
    R := CrystGroup( gens, S.identity );

    # bind the translations
    if IsBound( S.translations ) then
        AddTranslationsCrystGroup( R, S.translations*ci );
    fi;

    # bind Wyckoff positions
    if IsBound( S.wyckoffPositions ) then
        R.wyckoffPositions := List( S.wyckoffPositions, w -> rec(
            basis             := w.basis*ci,
            translation       := w.translation*ci,
            class             := w.class,
            stabilizer        := R.operations.Subgroup( Parent( R ),
                              List( w.stabilizer.generators, x -> C*x*Ci ) ), 
            operations        := WyckoffPositionOps,
            spaceGroup        := R,
            isWyckoffPosition := true ) );
    fi;
    return R;

end;
    

#############################################################################
##
#F  CrystGroupOps.Size( S ) . . . . . . . . .size for crystallographic groups
##
CrystGroupOps.Size := function( S )

   if Length( TranslationsCrystGroup( S ) ) > 0 then
       return "infinity";
   else
       return MatGroupOps.Size( PointGroup( S ) );
   fi;

end;


#############################################################################
##
#F  CrystGroupOps.Elements( S ) . . . . . . . . . . . . . . . . . . .Elements
##
CrystGroupOps.Elements := function( S )

   if not IsFinite( S ) then
       Error("S is infinite");
   else
       return MatGroupOps.Elements( S );
   fi;

end;


#############################################################################
##
#F  CrystGroupOps.\in( <e>, <S> ) .check membership in crystallographic group
##
CrystGroupOps.\in := function( e, S )

   local d, g, s, t, T, sol, x;

   if not IsMat( e ) then
       return false;
   fi;

   d := S.dimension - 1;
   g := e{[1..d]}{[1..d]};
   if not g in PointGroup( S ) then 
       return false; 
   fi;

   s := PreImagesRepresentative( S.pointHomom, g );
   t := s[d+1]{[1..d]} - e[d+1]{[1..d]};   
   T := TranslationsCrystGroup( S );

   if Length( T ) = 0 then
       return ForAll( t, x -> x=0 );
   else
       sol := SolutionMat( T, t );
       if IsBool( sol ) then
           return false;
       else
           return ForAll( sol, IsInt );
       fi;
   fi;

end;


#############################################################################
##
#F  CrystGroupOps.Subgroup( <parent>, <subgens> ) . .Subgroup of a CrystGroup
##
CrystGroupOps.Subgroup := function( parent, subgens )

    local S;

    S := GroupOps.Subgroup( parent, subgens );
    S.operations   := CrystGroupOps;
    S.dimension    := parent.dimension;
    S.field        := parent.field;
    S.isCrystGroup := true;

    return S;

end;


#############################################################################
##
#F  PointGroupOps . . . . . . . . . . . . . . . . . . . . . . . PointGroupOps
##
PointGroupOps := ShallowCopy( MatGroupOps );


#############################################################################
##
#F  PointGroupOps.FpGroup( P ) . . . . . . . . . . . FpGroup of a point group
##
PointGroupOps.FpGroup := function ( P )

    local S, d, U, A, F, pres, T, gens, inv, i, conj;

    if not IsBound( P.isPointGroup ) and P.isPointGroup then
        Error("P must be the point group of a CrystGroup");
    fi;

    # get the CrystGroup of P
    S := P.crystGroup;
    d := S.dimension - 1;

    # compute an isomorphic permutation group
    U := PermGroup( P );

    # distinguish between solvable and non-solvable case
    if IsSolvable( U ) then

        # switch to an AgGroup
        if not IsBound( U.agGroup ) then
            U.agGroup := AgGroup( U );
        fi;
        A := U.agGroup;

        # store new internalGenerators of P and add info
        P.internalGenerators := List( A.generators, x -> 
                        Image( U.bijection, Image( A.bijection, x ) ) );

        # get the new non-translational space group generators
        S.internalGenerators := List( P.internalGenerators, 
                      x -> PreImagesRepresentative( S.pointHomom, x ) ); 

        # conjugate the internal generators if necessary
        if IsBound( S.translations ) then
            ConjugateInternalGenerators( S );
        fi;

        P.agGroup    := A;
        P.isSolvable := true;
        S.isSolvable := true;

        return FpGroup( A );

    else

        # we don't want to have any generators eliminated
        F    := U.operations.FpGroup( U );
        pres := PresentationFpGroup( F, 0 );
        pres.generatorsLimit := Length( U.generators );
        SimplifyPresentation( pres );
        F := FpGroupPresentation( pres );

        # add info to P
        P.isSolvable := false;
        S.isSolvable := false;

        return F;

    fi;

end;


#############################################################################
##
#F  CrystGroupOps.PointGroup( S ) . . . . . . . . point group of a CrystGroup
##
CrystGroupOps.PointGroup := function ( S )

    local d, im, I, Pgens, Sgens, i, P; 

    d  := S.dimension-1;
    im := List( S.generators, m -> m{[1..d]}{[1..d]} );
    I  := IdentityMat( d );

    Pgens := [];
    Sgens := [];
    for i in [1..Length( im )] do
        if im[i] <> I and not im[i] in Pgens then
            Add( Pgens, im[i] );
            Add( Sgens, S.generators[i] );
        fi;
    od;

    P := Group( Pgens, IdentityMat( d ) );
    P.operations   := PointGroupOps;
    P.crystGroup   := S;
    P.isPointGroup := true;
    P.internalGenerators := Pgens;
    S.internalGenerators := Sgens;

    # set the point homomorphism
    S.pointHomom := GroupHomomorphismByImages( S, P, S.generators, im );
    S.pointHomom.isMapping := true;
    S.pointHomom.isGroupHomomorphism := true;
    S.pointHomom.isSurjetive := true;

    return P;

end;    


#############################################################################
##
#F  CrystGroupOps.FpGroup( S ) . . . . . . . . . . . . .FpGroup of CrystGroup
##
CrystGroupOps.FpGroup := function( S )

    local d, P, F, T, gensP, relsP, matsP, gensT, matsT, 
          gens, rels, rel, tail, vec, word, i, j, R;

    d := S.dimension-1;
    P := PointGroup( S );
    F := FpGroup( P );
    T := TranslationsCrystGroup( S );

    # if S is finite, and thus isomorphic to its PointGroup
    if Length( T ) = 0 then
        return F;
    fi;

    gensP := F.generators;
    relsP := F.relators;
    matsP := S.internalGenerators;
    
    gensT := AbstractGenerators( "t", Length( T ) );
    matsT := List( gensT, x -> IdentityMat( d+1 ) );
    for i in [1..Length( matsT )] do
        matsT[i][d+1][i] := 1;
    od;

    gens  := Concatenation( gensP, gensT );
    rels  := [];

    # compute tails
    for rel in relsP do
        tail := MappedWord( rel, gensP, matsP );
        vec  := - tail[d+1]{[1..d]};
        word := Product( List( [1..d], x -> gensT[x]^vec[x] ) );   
        Add( rels, rel * word );
    od;

    # compute operation
    for i in [1..Length( gensP )] do
        for j in [1..Length( gensT )] do
            rel  := Comm( gensT[j], gensP[i] );
            tail := Comm( matsT[j], matsP[i] );
            vec  := - tail[d+1]{[1..d]};
            word := Product( List( [1..d], x -> gensT[x]^vec[x] ) );   
            Add( rels, rel * word );
        od;
    od;

    # compute presentation of T
    for i in [1..Length(T)-1] do
        for j in [i+1..Length(T)] do
            Add( rels, Comm( gensT[j], gensT[i] ) );
        od;
    od;
    
    R := Group( gens, IdWord );
    R.relators := rels;
    return R;

end;



#############################################################################
##
## Routines for the determination of space groups for a given a point group
##
#############################################################################


#############################################################################
##
#F  FlattenMatMat( < MatMat > ). . . . . . . . . . Flatten matrix of matrices
##
FlattenMatMat := function( mat )
   # flatten a matrix whose entries are matrices to a normal matrix
   local m;
   m := mat;
   m := List( [1..Length(m[1])], 
              j -> Concatenation( List([1..Length(m)], i -> m[i][j] ) ) );
   m := TransposedMat( Concatenation( List( [1..Length(m)], 
                                      i -> TransposedMat(m[i]) ) ) );
   return m;
end;


#############################################################################
##
#F  NullMatMat( <d>, <d1>, <d2> ). . . . . . . d1xd2-matrix of d-NullMatrices
##
NullMatMat := function( d, d1, d2 )
   # return d1 x d2 matrix, whose entries are d x d NullMatrices
   return List( [1..d1], i -> List( [1..d2], j -> NullMat(d,d) ) );
end;


#############################################################################
##
#F  AugmentedMatrix( <matrix>, <trans> ). . . . .  construct augmented matrix
##
AugmentedMatrix := function( m, b )
   local g, t, x;
   g := Copy(m); for x in g do Add(x,0); od;
   t := Copy(b); Add(t,1);
   Add(g,t);
   return g;
end;


#############################################################################
##
#F  MakeSpaceGroup( <matgrp>, <trans> ). . . . . . . .  construct space group
##
MakeSpaceGroup := function( P, t )
   # construct space group from point group and translation vector
   local d, Pgens, Sgens, i, m, S;

   d     := P.dimension;
   Pgens := P.generators;

   # first the non-translational generators
   Sgens := List( [1..Length(Pgens)], 
                  i -> AugmentedMatrix( Pgens[i], t{[(i-1)*d+1..i*d]} ) );

   # the pure translation generators
   for i in [1..d] do
      m := IdentityMat( d+1 );
      m[d+1][i] := 1;
      Add( Sgens, m );
   od;

   # make the space group and return it
   S := CrystGroup( Sgens, IdentityMat(d+1) );
   AddTranslationsCrystGroup( S, IdentityMat( d ) );

   return S;

end;


#############################################################################
##
#F  GroupExtEquations( <d>, <gens>, <rels> ) . equations for group extensions
##
GroupExtEquations := function( d, gens, rels )
   # construct equations which determine the non-primitive translations
   local mat, i, j, r, prod;

   mat := NullMatMat( d, Length(gens), Length(rels) );
   for i in [1..Length(rels)] do
      prod := IdentityMat(d);
      r    := rels[i];
      for j in Reversed([1..Length(r)]) do
         if r[j]>0 then
            mat[ r[j] ][i] := mat[ r[j] ][i]+prod;
            prod := gens[ r[j] ]*prod;
         else
            prod := gens[-r[j] ]^-1*prod;
            mat[-r[j] ][i] := mat[-r[j] ][i]-prod;
         fi;
      od;
   od;
   return FlattenMatMat( mat );
end;


#############################################################################
##
#F  StandardTranslation( <trans>, <nullspace> ) . .reduce to std. translation
##
StandardTranslation := function( L, NN )
   # reduce non-primitive translations to "standard" form
   local N, j, k;

   # first apply "continuous" translations
   for N in NN[1] do
      j := PositionProperty( N, x -> x=1 );
      L := L-L[j]*N;
   od;
   L := List( L, FractionModOne );

   # and then "discrete" translations
   for N in NN[2] do
      j := PositionProperty( N, x -> x<>0 );
      k := Int(L[j]/N[j]);
      if k>0 then
         L := List( L-k*N, FractionModOne );
      fi;
   od;

   return L;

end;


#############################################################################
##
#F  SolveHomEquationsModZ( <mat> ) . . . . . . . . . . .  solve x*mat=0 mod Z
##
SolveHomEquationsModZ := function( M )

    local Q, L, N, N2, mul;

    Q := IdentityMat( Length(M) );
    
    # first diagonalize M
    M := TransposedMat(M);
    M := RowEchelonForm( M );
    while not IsDiagonalMat(M) do
        M := TransposedMat(M);
        M := RowEchelonFormT(M,Q);
        if not IsDiagonalMat(M) then
            M := TransposedMat(M);
            M := RowEchelonForm(M);
        fi;
    od;

    # and then determine the solutions of x*M=0 mod Z
    if Length(M)>0 then
        L := List( [1..Length(M)], i -> [ 0 .. M[i][i]-1 ] / M[i][i] );
        L := List( Cartesian( L ), l -> l * Q{[1..Length(M)]} );
    else
        L := NullMat( 1, Length(Q) );
    fi;

    # we later need the space in which one can freely shift
    # non-primitive translations; first the translations which 
    # can be applied with rational coefficients

    if Length(M)<Length(Q) then
        N := Q{[Length(M)+1..Length(Q)]};
        TriangulizeMat( N );
    else
        N := [];
    fi; 

    # and now those which allow only integral coefficients
    if N<>[] then
       N2  := List( N, n -> List( n, FractionModOne ) );
       mul := Lcm( Integers, List( Flat(N2), Denominator ) );
       N2  := HermiteNormalForm( N2*mul ) / mul;
       N2  := List( N2, n -> List( n, FractionModOne ) );
       N2  := Filtered( N2, n -> n<>0*N[1] );
    else
       N2  := [];
    fi;

    # reduce non-primitive translations to standard form
    L := Set( List( L, x -> StandardTranslation( x, [ N, N2 ] ) ) );

    return [ L, [ N, N2 ] ];

end;


#############################################################################
##
#F  EliminateEquivExtensions( <trans>, <nullspace>, <norm>, <grp> ) . . . . .
#F  . .  . . remove extensions equivalent by conjugation with elems from norm
##
EliminateEquivExtensions := function( ll, nn, norm, grp )

   # check for conjugacy with generators of the normalizer of grp in GL(n,Z)

   local cent, d, gens, sgens, res, orb, x, y, c, n, i, j, sg, h, m;

   cent := Filtered( norm, x -> ForAll( grp.generators, g -> x*g=g*x ) );
   SubtractSet( norm, cent );

   d     := grp.dimension;
   gens  := grp.generators;
   sgens := List( gens, g -> AugmentedMatrix( g, List( [1..d], x -> 0 ) ) );

   res := [ ];
   while ll<>[] do
      Add( res, ll[1] );
      orb := [ ll[1] ]; 
      for x in orb do

         # first the generators which are in the centralizer
         for c in cent do
            y := List([1..Length(gens)], i -> x{ [(i-1)*d+1..i*d] }*c );
            y := StandardTranslation( Concatenation(y), nn );
            if not y in orb then 
               Add( orb, y ); 
            fi;
         od;

         # then the remaining ones; this is more complicated
         for n in norm do
            for i in [1..Length(gens)] do
               for j in [1..d] do
                  sgens[i][d+1][j]:=x[(i-1)*d+j];
               od;
            od;
            sg := Group( sgens, IdentityMat( d+1 ) );
            sg.size := "infinity";
            h :=GroupHomomorphismByImages( sg, grp, sgens, gens );
            y :=[];
            for i in [1..Length(gens)] do
               m := PreImagesRepresentative( h, n*gens[i]*n^-1 );
               Append( y, m[d+1]{[1..d]}*n );
            od;
            y := StandardTranslation( y, nn );
            if not y in orb then
               Add( orb, y ); 
            fi;
         od;

      od;
      SubtractSet( ll, orb );
   od;

   return res;

end;


#############################################################################
##
#F  SpaceGroupsPointGroup( <grp> [, <norm>] ) . . compute group extensions,
#F     . . . . . . . . . . . inequivalent by conjugation with elems from norm
##
SpaceGroupsPointGroup := function( arg )

# construct group extensions of Z^d by grp

   local d, grp, norm, F, pres, rels, gens, mat, ext;

   grp := arg[1];
   if Length(arg)>1 then 
      norm := arg[2]; 
   else 
      norm := []; 
   fi;
   d := grp.dimension;

   # catch the trivial case
   if IsTrivial(grp) then
#      return [ [ List([1..d], i -> 0 ) ], [ IdentityMat(d), [] ] ];
      return [ MakeSpaceGroup( grp, [] ) ];
   fi;

   # first get group relators for grp
   # we don't want to have any generators eliminated
   F    := FpGroup( PermGroup( grp ) );
   pres := PresentationFpGroup( F, 0 );
   pres.generatorsLimit := Length( grp.generators );
   SimplifyPresentation( pres );

   # pres := PresentationViaCosetTable( PermGroup(grp) );
   rels := pres.tietze[6];
   gens := grp.generators;

   # construct equations which determine the non-primitive translations
   mat := GroupExtEquations( d, gens, rels );

   # now solve them modulo integers
   ext := SolveHomEquationsModZ( mat );
   
   # eliminate group extensions which are equivalent as space groups
   if Length( ext[1] )>2 then
      norm := Set( Filtered( norm, x -> not x in Elements(grp) ) );
      if Length(norm)>0 then
         ext[1] := EliminateEquivExtensions( ext[1], ext[2], norm, grp );
      fi;
   fi;
   
#   return ext;
   return List( ext[1], x -> MakeSpaceGroup( grp, x ) );

end;


#############################################################################
##
## Routines for the determination of maximal subgroups of space groups
##
#############################################################################


#############################################################################
##
#F  UseSmash() . . . . . . . . . . . .switch to using the Smash share package
##
UseSmash := function()
    RequirePackage("smash");
    MatGroupOps.BasesMaximalSubmodules 
          := MatGroupOps.BasesMaximalSubmodules1;
end;


#############################################################################
##
#F  NoSmash() . . . . don't use the Smash share package - this is the default
##
NoSmash := function()
    MatGroupOps.BasesMaximalSubmodules 
          := MatGroupOps.BasesMaximalSubmodules2;
end;


# if Smash is not loaded, get around smash definitions
if not IsBound( GModule ) then GModule := false; fi;
if not IsBound( DualGMod ) then DualGMod := false; fi;
if not IsBound( ChopGMod ) then ChopGMod := false; fi;
if not IsBound( MinSubGMods ) then MinSubGMods := false; fi;


#############################################################################
##
#F  CoefficientsMod( base, v ) . . . . . . . coefficients of v in factorspace
##
CoefficientsMod := function( base, v )

    local head, i, zero, coeff, w, j, h;

    if not IsBound( base.fullbase ) then
        base.fullbase := Concatenation( base.subspace, base.factorspace );
    fi;
    if not IsBound( base.depth ) then
        head := [];
        for i in [1..Length( base.fullbase )] do
            head[i] := DepthVector( base.fullbase[i] );
        od;
        base.depth := head;
    fi;

    zero  := base.fullbase[1] * base.field.zero;
    coeff := Copy( zero );
    w     := v;
    while w <> zero do
        j := DepthVector( w );
        h := Position( base.depth, j );
        coeff[h] := coeff[h] + w[j];
        w := w - w[j] * base.fullbase[h];
    od;
    return coeff{[Length(base.subspace)+1..Length(base.fullbase)]};
end;


#############################################################################
##
#F  TriangularizeMatVector  . . . . . . . . . . compute triangularized matrix
##
## This function computes the upper triangular form of the integer
## matrix M via elementary row operations and performs the same
## operations on the (column) vector b.
##
## The function works in place.
##
TriangularizeMatVector := function( M, b )
    local   zero,  c,  r,  i,  t;
    
    zero := M[1][1] * 0;
    
    c := 1; r := 1;
    while c <= Length(M[1]) and r <= Length(M) do
        i := r; while i <= Length(M) and M[i][c] = zero do i := i+1; od;

        if i <= Length(M) then
            t := b[r]; b[r] := b[i]; b[i] := t;
            t := M[r]; M[r] := M[i]; M[i] := t;

            b[r] := b[r] / M[r][c];
            M[r] := M[r] / M[r][c]; 

            for i in [1..r-1] do 
                b[i] := b[i] - b[r] * M[i][c];
                AddCoeffs( M[i], M[r], -M[i][c] );
            od;
            for i in [r+1..Length(M)] do
                b[i] := b[i] - b[r] * M[i][c];
                AddCoeffs( M[i], M[r], -M[i][c] );
            od;
            r := r+1;    
        fi;
        c := c+1;
    od;
    
    for i in [r..Length(M)] do Unbind(M[i]); od;
end;


#############################################################################
##
#F  SolutionInhomEquations  . . .  solve an inhomogeneous system of equations
##
## This function computes the set of solutions of the equation
##
##                           X * M = b.
##
SolutionInhomEquations := function( M, b )
    local   zero,  one,  i,  c,  d,  r,  heads,  S,  v;
    
    zero := M[1][1] * 0;
    one  := M[1][1] ^ 0;
    
    M := TransposedMat( M );
    d := Length(M[1]);
    b := Copy( b );
    
    TriangularizeMatVector( M, b );
    for i in [Length(M)+1..Length(b)] do
        if b[i] <> zero then return false; fi;
    od;

    # determine the null space
    c := 1; r := 1; heads := []; S := [];
    while r <= Length(M) and c <= d do
        
        while c <= d and M[r][c] = zero do 
            v := zero * [1..d];
            v{heads} := M{[1..r-1]}[c];  v[c] := -one;
            Add( S, v );
            c := c+1;
        od;
        
        if c <= d then
            Add( heads, c ); c := c+1; r := r+1;
        fi;
    od;
    while c <= d do
        v := zero * [1..d];
        v{heads} := M{[1..r-1]}[c];  v[c] := -one;
        Add( S, v );
        c := c+1;
    od;
    
    # one particular solution
    v := zero * [1..d];
    v{heads} := b{[1..Length(M)]};
    TriangulizeMat( S );
    return rec( basis := S, translation := v );
end;


#############################################################################
##
#F  SolutionHomEquations  . . . . . . solve a homogeneous system of equations
##
## This function computes the set of solutions of the equation
##
##                           X * M = 0.
##
SolutionHomEquations := function( M )
    
    return SolutionInhomEquations( M, 
                   0 * Field(M[1][1]).zero * [1..Length(M)] );
end;


#############################################################################
##
#F  MatJacobianMatrix( G, mats )
##
MatJacobianMatrix := function( G, mats )
    local   J,  D,      # the result,  the one `column' of J
            imats,      # list of inverted matrices
            d,          # dimension of matrices
            j,  k,  l,  # loop variables
            r,          # run through the relators
            h, p;       # generator in r,  its position in G.generators
    
    
    if Length(G.generators) = 0 then return []; fi;
    
    d := Length( mats[1] );
    imats := List( mats, m->m^-1 );
    
    J := List( [1..d*Length(G.generators)], k->[] );
    for j in [1..Length(G.relators)] do
        r := G.relators[j];
        D := NullMat( d*Length(G.generators), d );
    
        for k in [1..LengthWord(r)] do
            h := Subword( r,k,k );
            p := Position( G.generators, h );
            
            if p <> false then
                D := D * mats[p];
                for l in [1..d] do
                    D[(p-1)*d+l][l] := D[(p-1)*d+l][l] + 1;
                od;
            else
                p := Position( G.generators, h^-1 );
                for l in [1..d] do
                    D[(p-1)*d+l][l] := D[(p-1)*d+l][l] - 1;
                od;
                D := D * imats[p];
            fi;

            J{[1..Length(G.generators)*d]}{[(j-1)*d+1..j*d]} := D;
        od;
    od;
    return J;
end;


#############################################################################
##
#F  OneCoboundariesSG( <G>, <mats> )  . . . . . . . . . . . . . . B^1( G, M )
##
OneCoboundariesSG := function( G, mats )
    local   d,  I,  S,  i;
    
    d := Length( mats[1] );
    I := IdentityMat( d );
    
    if Length(mats) <> Length(G.generators) then
        return Error( "As many matrices as generators expected" );
    fi;
    
    S := List( [1..d], i->[] );
    for i in [1..Length(mats)] do
        S{[1..d]}{[(i-1)*d+1..i*d]} := mats[i] - I;
    od;

    TriangulizeMat( S );

    return S;
end;


#############################################################################
##
#F  OneCocyclesVector( <G>, <mats>, <b> ) . . . . . . . . . . . . . . . . . . 
##
OneCocyclesVector := function( G, mats, b )
    local   J,  L;
    
    J := MatJacobianMatrix( G, mats );
    
    b := Concatenation( b );
    L := SolutionInhomEquations( J, b );

    return L;
end;


#############################################################################
##
#F  OneCocyclesSG( <G>, <mats> )  . . . . . . . . . . . . . . . . Z^1( G, M )
##
OneCocyclesSG := function( G, mats )
    local   J,  L;
    
    J := MatJacobianMatrix( G, mats );
    L := SolutionHomEquations( J ).basis;
        
    return L;
end;


#############################################################################
##
#F  OneCohomology( <G>, <mats> )  . . . . . . . . . . . . . . . . H^1( G, M )
##
OneCohomologySG := function( G, mats )
    
    return rec( cocycles     := OneCocyclesSG( G, mats ),
                coboundaries := OneCoboundariesSG( G, mats ) );
    
end;


#############################################################################
##
#F  ListOneCohomology( <H> ) . . . . . . . . . . . . . . . . . . . . list H^1
##
##    Run through the triangularized basis of Z and find those head
##    entries which do not occur in the basis of B.  For each such
##    vector in the basis of Z we need to run through the coefficients
##    0..p-1. 
##
ListOneCohomology := function( H )
    local   Z,  B,  C,  zero,  coeffs,  h,  j,  i;
    
    Z := H.cocycles;
    B := H.coboundaries;
    if Length(Z) = 0 then return B; fi;
    
    C := Elements( Field( Z[1][1] ) );
    zero := Z[1][1]*0;
    
    coeffs := [];
    h := 1; j := 1;
    for i in [1..Length(Z)] do
        while Z[i][h] = zero do h := h+1; od;
        if j > Length(B) or B[j][h] = zero then
            coeffs[i] := C;
        else
            coeffs[i] := [ zero ];  j := j+1;
        fi;
    od;

    return List( Cartesian( coeffs ), t->t*Z );
end;


#############################################################################
##
#F  ComplementsSG . . . . . . . . . . . . compute complements up to conjugacy
##
ComplementsSG := function( G, mats, b )
    local   Z,  B,  C,  d,  n;
    
    if Length(G.generators) = 0 then return [ [] ]; fi;
    
    Z := OneCocyclesVector( G, mats, b );
    if Z = false then return []; fi;
    
    B := OneCoboundariesSG( G, mats );
    
    C := ListOneCohomology( 
                 rec( cocycles := Z.basis,
                      coboundaries := B ) );
    
    C := List( C, c->c + Z.translation );
    
    d := Length( mats[1] );
    n := Length(G.generators);
    
    return List( C, c->List( [0..n-1], x->c{ [1..d] + x*d } ) );
end;


#############################################################################
##
#F  BasesMaximalSubmodules1( G ) . . . . . .maximal $G$-submodules with smash
##
MatGroupOps.BasesMaximalSubmodules1 := function( G )
    local matlist, field, dim, dual, mini, module, chop, subm;

    # set up
    matlist := G.generators;
    field   := G.field;
    dim     := G.dimension;

    # the trivial cases
    if dim = 1 then
        return [[]];
    elif Length( matlist ) = 0 then
        return List( NormedVectors( RowSpace( field,
               IdentityMat( dim, field ) ) ), x -> [x] );
    fi;

    module := GModule( matlist, field );
    dual := DualGMod( module );
    mini := [];
    chop := ChopGMod( dual );
    for subm in chop do
        Append( mini, MinSubGMods( subm[1], dual ) );
    od;
    if mini = [[]] then return [[]]; fi;
    return List( mini, x -> NullspaceMat( TransposedMat( x ) ) );
end;


#############################################################################
##
#F  BasesMaximalSubmodules2( G ) . . . . . . maximal $G$-submodules - generic 
##
MatGroupOps.BasesMaximalSubmodules2 := function( G )

    local mats, Mdescr,          # dual representation
          modules, module,       # list of bases of modules
          comp,                  # composition factors of the dual
          a,                     # peakword 
          kern,                  # kernel of peakwords
          dim,                   # dim of kernel
          spin,                  # spinning of a vector
          j;

    if Length( G.generators ) = 0 then
        modules := List( NormedVectors( RowSpace( G.field, G.dimension ) ),
                   x -> [x] ); 
    else

        # go over to dual represenation
        mats := List( G.generators, x -> TransposedMat( x ) );
        Mdescr := [mats, G.field, G.dimension];
        
        # compute composition factors
        comp  := MTXOps.CompositionFactors( Mdescr );

        # if the dual is irreducible
        if Length( comp.factors ) = 1  then
            return [[]];
        fi;

        # otherwise calculate all minimal submodules of Mdescr
        modules := [];
        a := MTXOps.Peakwords( Mdescr, comp );
        for j  in [ 1 .. Length(a.peaks) ]  do
            if not IsBool( a.peaks[j] )  then
                kern := NullspaceMat( a.peaks[j] );
                dim  := a.dims[j];
                spin := MTXOps.SpinedBase( Mdescr, kern, dim );
                Append( modules, Set( spin ) );
            else
                Print("#W  could not find all maximal submodules \n");
            fi;
        od;
    fi;

    # convert them to maximal ones 
    for j  in [1..Length(modules)]  do
        modules[j] := NullspaceMat( TransposedMat( modules[j] ) );
    od;
    return modules;
end;


#############################################################################
##
#F  MatGroupOps.BasesMaximalSubmodules( G ) . . . .either of the two versions
##
MatGroupOps.BasesMaximalSubmodules := MatGroupOps.BasesMaximalSubmodules2;


#############################################################################
##
#F  MaximalSubgroupsRepsTG( < S > ) . .translationengleiche maximal subgroups
##
## This function computes conjugacy class representatives of the maximal 
## subgroups of $S$ which contain $T$. Note that this function may be
## slow, if $S$ is not solvable.
##
MaximalSubgroupsRepsTG := function( S )

    local  d, P, U, T, F, A, sp, max, ind, l, i, gens, g, exp, h, j, 
           trans, conj, invconj, sub; 

    d := S.dimension - 1;
    P := PointGroup( S );
    U := PermGroup( P );

    # catch a trivial case
    if Size( U ) = 1 then
        return [];
    fi;

    # compute the translation generators
    T     := TranslationsCrystGroup( S );
    trans := List( T, x -> IdentityMat( d+1 ) );
    for i in [1..Length(T)] do
        trans[i][d+1]{[1..d]} := T[i];
    od;

    # first the solvable case
    if IsSolvable(U) then

        F := FpGroup( P );
        A := P.agGroup;

        # compute maximal subgroups in ag group
        sp  := SpecialAgGroup( A );
        max := Concatenation( sp.operations.RepsMaximalSubgroups( sp ) );
        max := List( max, x -> Image( sp.bijection, x ) ); 
    
        # construct conjugators if necessary
        if not S.isStandard then
            conj                 := IdentityMat( d+1 );
            conj{[1..d]}{[1..d]} := S.internalBasis;
            invconj              := conj^-1;
        fi;

        # compute preimages in space group and construct subgroups
        for i in [1..Length(max)] do

            gens := [];
            for g in max[i].generators do
                exp := Exponents( A, g );
                h := Copy( S.identity );
                for j in [1..Length(exp)] do
                    h := h * S.internalGenerators[j]^exp[j];
                od;
                Add( gens, h );
            od;

            if not S.isStandard then
                for j in [1..Length(gens)] do
                    gens[j] := invconj*gens[j]*conj;
                od;
            fi;
            Append( gens, trans );
            sub := S.operations.Subgroup( Parent( S ), gens );
            AddTranslationsCrystGroup( sub, T );

            if IsParent( S ) then
                sub.index := Index( A, max[i]);
            elif IsBound( S.index ) then
                sub.index := S.index*Index( A, max[i]);
            fi;
            max[i] := sub;

        od;
        return max;

    fi;

    # now the non-solvable case
    max := U.operations.ConjugacyClassesMaximalSubgroups( U );
    max := List( max, x -> Representative( x ) );
    ind := List( max, x -> Index( U, x ) );
    max := List( max, x -> Image( U.bijection, x ) );

    # go back to S, and construct the subgroups
    for i in [1..Length(max)] do 

        gens := List( max[i].generators, 
                      x -> PreImagesRepresentative( S.pointHomom, x ) );
        Append( gens, trans );
        sub  := S.operations.Subgroup( Parent( S ), gens ); 
        AddTranslationsCrystGroup( sub, T );

        if IsParent( S ) then
            sub.index := ind[i];
        elif IsBound( S.index ) then
            sub.index := S.index*ind[i];
        fi;
        max[i] := sub;

    od;
    return max;

end; 


#############################################################################
##
#F  CocycleInfo( <G> ) . . . . . compute info about the extension of P with T
##
CocycleInfo := function( G )

    local F, d, coc, rel, new;

    F := FpGroup( PointGroup( G ) );
    d := G.dimension - 1;

    coc := [];
    for rel in F.relators do
        new := MappedWord( rel, F.generators, G.internalGenerators );
        Add( coc, new[d+1]{[1..d]} );
    od;
    return coc;
end;


#############################################################################
##
#F  SimpleGenerators( d, gens ) . . . . . . . . . . . simplify the generators 
##
SimpleGenerators := function( d, gens )

    local I, new, g, trans, t, m;

    I     := IdentityMat( d );
    new   := [];
    trans := [];
    for g in gens do
        if g{[1..d]}{[1..d]} = I then
            Add( trans, g[d+1]{[1..d]} );
        else
            Add( new, g );
        fi;
    od;
    trans := ReducedTranslationBasis( trans );

    # add the new translation generators
    for t in trans do
        m := IdentityMat( d+1 );
        m[d+1]{[1..d]} := t;
        Add( new, m );
    od;
    return [ new, trans ];

end; 


#############################################################################
##
#F  MaximalSubgroupsRepsKG( < G >, <ps> ) . .klassengleiche maximal subgroups
##
## This function computes represenatives of the conjugacy classes of maximal
## subgroups of $G$ which have $p$-power index for some $p$ in the list $ps$
## and do not contain $T$. 
## In the case that $G$ is solvable it is more efficient to use the function
## 'MaximalSubgroupSG' and filter the corresponding maximal subgroups.
##
MaximalSubgroupsRepsKG := function( G, primes )

    local pres, coc, rep, d, n, maximals, p, field, V, repp, cocp, 
          matgrp, mods, module, U, F, hom, cocin, repin, comp, c,
          modgens, powers, vec, elm, basis, cocpre, gens, i, j, h,
          base, primeslist, d, Ggens, conj, invconj, T, trans;

    # check argument
    if IsInt( primes ) then
        primeslist := [primes];
    else
        primeslist := primes;
    fi;

    T := TranslationsCrystGroup( G );
    n := Length( T );

    # extract the point group
    pres  := FpGroup( PointGroup( G ) );
    coc   := List( CocycleInfo( G ), x -> x{[1..n]} );
    rep   := List( G.pointGroup.internalGenerators, x -> x{[1..n]}{[1..n]} );
    d     := G.dimension - 1;
    Ggens := G.internalGenerators;
    trans := List( T, x -> IdentityMat( d+1 ) );
    for i in [1..n] do
        trans[i][d+1][i] := 1;
    od; 

    # determine the conjugators if necessary
    if not G.isStandard then
        conj := IdentityMat( d+1 );
        conj{[1..d]}{[1..d]} := G.internalBasis;
        invconj := conj^-1;
    fi;

    # view them as matrices over GF(p)
    maximals := [];
    for p in primeslist do
        field := GF(p);
        V     := RowSpace( field, n );
        repp  := List( rep, x -> x * field.one );
        cocp  := List( coc, x -> x * field.one );
        matgrp:= Group( repp, IdentityMat( n, field ) );
        mods  := matgrp.operations.BasesMaximalSubmodules( matgrp );
        powers:= List( trans, x -> x^p );

        # compute induced operation on T/maxmod and induced cocycle
        for module in mods do

            # compute group of translations of maximal subgroup
            modgens := [];
            for vec in module do
                elm := G.identity;
                for j in [1..Length( vec )] do
                    elm := elm * trans[j]^IntVecFFE(vec, j);
                od;
                Add( modgens, elm );
            od;
            Append( modgens, powers );

            # compute subspace
            if Length( module ) = 0 then
                U := Subspace( V, [Zero(V)] );
            else 
                U := Subspace( V, module );
            fi;
           
            # compute quotient space
            F := V / U;
            basis := Basis( F );
            base := BaseSteinitz( V, U );
            cocin := List( cocp, x -> CoefficientsMod( base, x ) );
            repin := InducedActionSpaceMats( basis, repp );

            # use complement routine
            comp := ComplementsSG( pres, repin, cocin );

            # compute generators of G corresponding to complements
            for i in [1..Length( comp )] do
                cocpre := List( comp[i], x -> x * base.factorspace );
                gens := [];
                for j in [1..Length( cocpre )] do
                    elm := Ggens[j];
                    for h in [1..Length( cocpre[j] )] do
                        elm := elm*trans[h]^IntVecFFE(cocpre[j], h);
                    od;
                    Add( gens, elm );     
                od;
        
                # append generators of group of translations 
                Append( gens, modgens );

                # conjugate generators if necessary
                if not G.isStandard then
                    for j in [1..Length(gens)] do
                        gens[j] := invconj*gens[j]*conj;
                    od;
                fi;

                # construct subgroup and append index
                gens := SimpleGenerators( d, gens );
                comp[i] := G.operations.Subgroup( Parent( G ), gens[1] );
                AddTranslationsCrystGroup( comp[i], gens[2] );
                if IsParent( G ) then
                    comp[i].index := p^(n - Length( module ));
                elif IsBound( G.index ) then
                    comp[i].index := G.index*p^(n - Length( module ));
                fi;
            od;
            Append( maximals, comp );
        od;
    od;
    return maximals;
end;


#############################################################################
##
#F  MaximalSubgroupsRepsSG( <G>, <p> ) . . .maximal subgroups of solvable <G>
##
## This function computes representatives of the conjugacy classes of the 
## maximal subgroups of $p$-power index in $G$ in the case that $G$ is 
## solvable.
##
MaximalSubgroupsRepsSG := function( G, p )

    local P, U, T, Ggens, F, n, d, t, A, S, kernel, imgs, max, i, w, M, m, 
          gens, exp, h, g, conj, invconj, j;

    P := PointGroup( G );
    U := PermGroup( P );
    if not IsSolvable(U) then
        Error("G must be solvable \n");
    fi;

    F := Copy( FpGroup( G ) );
    T := TranslationsCrystGroup( G );

    Ggens := G.internalGenerators;
    n := Length( Ggens );
    d := G.dimension - 1;
    t := Length( T );
    Append( F.relators, List( [n+1..n+t], x -> F.generators[x]^p ) );
    A := AgGroupFpGroup( F );
    S := SpecialAgGroup( A );

    # compute generators of kernel G -> A and preimages
    kernel := List( [1..t], x -> IdentityMat( d+1 ) );
    for i in [1..t] do
        kernel[i][d+1][i] := 1;
    od;
    imgs   := Concatenation( Ggens, Copy( kernel ) );
    for i in [1..t] do
        kernel[i][d+1][i] := p;
    od;
    
    # compute maximal subgroups of S
    max := [];
    for i in [1..Length(S.first)-1] do
        w := S.weights[S.first[i]];
        if w[2] = 1 and w[3] = p then
            M := S.operations.MatGroupSagGroup( S, i );
            m := M.operations.BasesMaximalSubmodules( M );
            Append( max, S.operations.SubgroupsModules( S, i, m ) );
        fi;
    od; 

    # determine the conjugators if necessary
    if not G.isStandard then
        conj := IdentityMat( d+1 );
        conj{[1..d]}{[1..d]} := G.internalBasis;
        invconj := conj^-1;
    fi;

    # compute corresponding subgroups in G
    for i in [1..Length(max)] do

        M := Image( S.bijection, max[i] );
        gens := [];
        for g in Cgs(M) do
            exp := Exponents( A, g );
            h := Product( List( [1..Length(exp)], x -> imgs[x]^exp[x] ) );
            Add( gens, h );
        od;
        Append( gens, kernel );

        if not G.isStandard then
            for j in [1..Length(gens)] do
                gens[j] := invconj*gens[j]*conj;
            od;
        fi;
        gens := SimpleGenerators( d, gens );
        M := G.operations.Subgroup( Parent( G ), gens[1] );
        AddTranslationsCrystGroup( M, gens[2] );

        if IsParent( G ) then
            M.index := Index( S, max[i] );
        elif IsBound( G.index ) then
            M.index := G.index*Index( S, max[i] );
        fi;
        max[i] := M;

    od;
    return max;
end;


#############################################################################
##
#F  MaximalSubgroupsRepresentatives( S [, flag] [, index] )  conj. class reps
##
MaximalSubgroupsRepresentatives := function( arg )

    local S, i, flag, index, str, IsTG, sub, root;

    S := arg[1];

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;

    # get flag and index
    for i in [2..Length( arg )] do
        if IsString( arg[i] ) then 
            flag  := arg[i];
        elif IsInt( arg[i] ) then
            index := [ arg[i] ];
        elif IsList( arg[i] ) and IsInt( arg[i][1] ) then
            index := arg[i];
        else
            Error("wrong argument types");
        fi;
    od;

    # check arguments
    str := Concatenation("if flag 'translationEqual' is not present, ",
                         "an index must be specified" );
    if IsBound( flag ) then
        if flag <> "translationEqual" then
            if not IsBound( index ) then
                Error( str );
            fi;
            if flag <> "classEqual"  then
                Error( 
                "flag must be either 'translationEqual' or 'classEqual'");
            fi;
        fi;
    else
        if not IsBound( index ) then
            Error( str );
        fi;            
    fi;

    # check if subgroup is TG with S
    IsTG := function( x )
        return TranslationsCrystGroup( x ) = TranslationsCrystGroup( S );
    end;

    # if i = prime^k, return prime, else false
    root := function( i )
        local lst;
        lst := FactorsInt( i );
        if Length( lst ) = 1 then
            return lst[1];
        else
            return false;
        fi;
    end;

    # we need only TG subgroups
    if IsBound( flag ) and flag = "translationEqual" then
        sub := MaximalSubgroupsRepsTG( S );
        if IsBound( index ) then
            sub := Filtered( sub, x -> root( x.index ) in index );
        fi;

    # S is solvable, and we need also KG subgroups
    elif IsSolvable( PermGroup( PointGroup( S ) ) ) then
        sub := Concatenation( List( index, 
                    p -> MaximalSubgroupsRepsSG( S, p ) ) );
        if IsBound( flag ) then
            sub := Filtered( sub, 
                       x -> IsTG( x ) = (flag = "translationEqual") ); 
        fi;

    # S is not solvable, and we need all subgroups
    elif not IsBound( flag ) then 
        sub := MaximalSubgroupsRepsTG( S );
        sub := Filtered( sub, x -> root( x.index ) in index );
        sub := Concatenation( sub, MaximalSubgroupsRepsKG( S, index ) );

    # S is not solvable, and we need only KG subgroups
    else
        sub := MaximalSubgroupsRepsKG( S, index );
    fi;

    return sub;

end;


#############################################################################
##
## Routines for the determination of Wyckoff positions
##
#############################################################################


#############################################################################
##
#F  NormalizeAffineSubspace . . . . . . . . normal form of an affine subspace
##
NormalizeAffineSubspace := function( U )
    local  n, N, N2, mul, v, t, d, T, j, k;
    
    v := List( U.translation, FractionModOne );
    U.basis := HermiteNormalForm( U.basis );
    N := ShallowCopy( U.basis );
   
    #  Reduce the translation modulo the basis
    if N<>[] then

        # set components of v to zero where possible
        TriangulizeMat( N );
        for n in N do
            j := PositionProperty( n, x -> x=1 );
            v := v-v[j]*n;
        od;
        v := List( v, FractionModOne );

        # actually, we wanted to set those components only to an integer
        N   := List( N, n -> List( n, FractionModOne ) );
        mul := Lcm( Integers, List( Flat(N), Denominator ) );
        N   := HermiteNormalForm( N*mul ) / mul;
        N   := List( N, n -> List( n, FractionModOne ) );
        N   := Filtered( N, n -> n<>0*v );
        for n in N do
            j := PositionProperty( n, x -> x<>0 );
            k := Int(v[j]/n[j]);
            if k<>0 then
                v := List( v-k*n, FractionModOne );
            fi;
        od;

    fi;

    t := v - U.translation;
    U.translation := v;

    # adjust the stabilizer    
    if IsBound( U.stabilizer ) then
        d := Length(t);
        T := IdentityMat( d+1 ); T[d+1]{[1..d]} := t;
        U.stabilizer.generators
          := List( U.stabilizer.generators, g->T*g*T^-1 );
    fi;
    
end;


#############################################################################
##
#F  ImageEquivAffSpace  . . . . . . . . . . . . . image of an affine subspace
##                                                under a space group element
##
ImageEquivAffSpace := function( S, g )
    local   U,  d,  t,  l;
    
    U := rec();
    
    d := Length(g)-1;
    l := g{[1..d]}{[1..d]};
    t := g[d+1]{[1..d]};
    
    U.basis := List( S.basis, v->v*l );
    U.translation  := S.translation * l + t;
    
    NormalizeAffineSubspace( U );
    return U;
end;


#############################################################################
##
#F  InAffineSpaces  . . . . . .  test membership in a set of affine subspaces
##
InAffineSpaces := function( lst, sp )
    local   l;
    
    for l in lst do
        if l.translation = sp.translation and
           l.basis = sp.basis then
            return true;
        fi;
    od;
    return false;
end;


#############################################################################
##
#F  OrbitAffineSpaces . . . . . .  orbit of a space group on affine subspaces
##
OrbitAffineSpaces := function( S, sp )
    local   orb, gen, g, img;
    
    orb := [ sp ];
    gen := S.internalGenerators;
    for sp in orb  do
        for g in gen  do
            img := ImageEquivAffSpace( sp, g );
            if not InAffineSpaces( orb, img ) then
                Add( orb, img );
            fi;
        od;
    od;
    
    return orb;
end;


#############################################################################
##
#F  SolveOneInhomEquationModZ . . . . . . . .  solve one inhom equation mod Z
##
##  Solve the inhomogeneous equation
##
##            a x = b (mod Z).
##
##  The set of solutions is
##                    {0, 1/a, ..., (a-1)/a} + b/a.
##  Note that 0 < b <  1, so 0 < b/a and (a-1)/a + b/a < 1.
##
SolveOneInhomEquationModZ := function( a, b )
    
    return [0..a-1] / a + b/a;
end;

        
#############################################################################
##
#F  SolveInhomEquationsModZ . . . . .solve an inhom system of equations mod Z
##
##  This function computes the set of solutions of the equation
##
##                           x * M = b  (mod Z).
##
##  RowEchelonFormT() returns a matrix Q such that Q * M is in row echelon
##  form.  This means that (modulo column operations) we have the equation
##         x * Q^-1 * D = b       with D a diagonal matrix.
##  Solving y * D = b we get x = y * Q.
##
SolveInhomEquationsModZ := function( M, b )
    local   Q,  j,  L,  space,  i,  v;
    
    b := ShallowCopy(b);
    Q := IdentityMat( Length(M) );
    
    M := TransposedMat(M);
    M := RowEchelonFormVector( M,b );
    while not IsDiagonalMat(M) do
        M := TransposedMat(M);
        M := RowEchelonFormT(M,Q);
        if not IsDiagonalMat(M) then
            M := TransposedMat(M);
            M := RowEchelonFormVector(M,b);
        fi;
    od;

    ##  Now we have D * y = b with  y =  Q * x

    ##  Check if we have any solutions modulo Z.
    for j in [Length(M)+1..Length(b)] do
        if not IsInt( b[j] ) then
            return [ [], [] ];
        fi;
    od;

    ##  Solve each line in D * y = b separately.
    L := List( [1..Length(M)], i->SolveOneInhomEquationModZ( M[i][i],b[i] ) );
    
    L := Cartesian( L );
    L := List( L, l->Concatenation( l,  0 * [Length(M)+1..Length(Q)] ) );
    L := List( L, l-> l * Q );

    L := List( L, l->List( l, q->FractionModOne(q) ) );
    
    return [ L, Q{[Length(M)+1..Length(Q)]} ];
end;


#############################################################################
##
#F  FixedPointsModZ  . . . . . . fixed points up to translational equivalence
##
##  This function takes a space group and computes the fixpoint spaces of
##  this group modulo the translation subgroup.  It is assumed that the
##  translation subgroup has full rank.
##
FixedPointsModZ := function( gens, d )
    local   I,  M,  b,  i,  g,  f,  F;
    
    #  Solve x * M + t = x modulo Z for all pairs (M,t) in the generators.
    #  This leads to the system
    #        x * [ M_1 M_2 ... ] = [ b_1 b_2 ... ]  (mod Z)

    M := List( [1..d], i->[] ); b := []; i := 0;
    I := IdentityMat(d+1);
    for g in gens do
        g := g - I;
        M{[1..d]}{[1..d]+i*d} := g{[1..d]}{[1..d]};
        Append( b, -g[d+1]{[1..d]} );
        i := i+1;
    od;

    # Catch trivial case
    if Length(M[1]) = 0 then M := List( [1..d], x->[0] ); b := [1..d]*0; fi;
    
    ##  Compute the spaces of points fixed modulo translations.
    F := SolveInhomEquationsModZ( M, b );
    b := HermiteNormalForm( F[2] );
    F := List( F[1], f -> rec( translation := f, basis := b ) );
    for f in F do
        NormalizeAffineSubspace(f);
    od;
    return F;

end;
    

#############################################################################
#F  CrystGroupOps.IsSymmorphic( S ) . . . . . . . . . . . . .is S symmorphic? 
##
CrystGroupOps.IsSymmorphic := function( S )

    local T, G, gen, d;

    # check if we have indeed a space group
    T := TranslationsCrystGroup( S );
    if not S.isSpaceGroup then
        Error("S must be a space group");
    fi;

    gen := S.internalGenerators;
    d   := S.dimension - 1;

    # symmorphic if point group trivial or 
    # there exist Wyckoff positions with full symmetry
    if Length( gen ) = 0 then
        return true;
    else
        return Length( FixedPointsModZ( gen, d ) ) > 0;
    fi;

end;          


#############################################################################
##
#F  AdjustStabilizer . . . . . . . .adjust translations of Wyckoff stabilizer
##
AdjustStabilizer := function( S, w, gens, size )

    ##  Adjust the translations of the subgroup generators 
    ##  such that w is fixed

    local d, newgens, g, t;
    
    d := S.dimension - 1;
    newgens := [];
    for g in gens do
        t := Copy(w.translation);  Add(t,1);
        t := t * g - t;
        g := Copy(g); g[d+1] := g[d+1] - t;
        Add( newgens, g );
    od;
    w.stabilizer := S.operations.Subgroup( Parent( S ), newgens );
    w.stabilizer.size := size;

end;


#############################################################################
##
#F  IdentifyWyckoffOrbits( w, norm ) . .  identify Wyckoff pos. in same orbit
##
IdentifyWyckoffOrbits := function( w, norm )

    local orbs, orb, n, f;  

    # identify positions in the same orbit
    orbs := [];
    while w <> [] do
        orb := [];
        for n in norm do
            f := ImageEquivAffSpace( w[1], n );
            if not InAffineSpaces( orb, f ) then
                Add( orb, f );
            fi;
        od;
        Add( orbs, w[1] );
        SubtractSet( w, orb );
    od;
    return orbs;

end;


#############################################################################
##
#F  CheckWyckoffOrbits( w, cosrep ) . . . . . .check length of Wyckoff orbits
##
CheckWyckoffOrbits := function( w, cosrep )

    # eliminate orbits which are too short

    local d, b, x, y, z, orbs, orb, short, i;  

    d := Length( w[1].translation );
    b := w[1].basis;

    # filter group elements which keep basis pointwise invariant
    if b<>[] then
        cosrep := Filtered( cosrep, x -> b*x{[1..d]}{[1..d]}=b );
    fi;

    # check orbit length
    orbs := [];
    for x in w do
        y := Copy( x.translation ); Add( y, 1 );
        short := false;
        orb := []; i := 1;
        while not short and i <= Length( cosrep ) do
            z := List( y*cosrep[i], FractionModOne ); 
            if z in orb then 
                short := true;
            else
                Add( orb, z );
            fi;
            i := i+1;
        od;
        if not short then Add( orbs, x ); fi;
    od;
    return orbs;

end;


#############################################################################
##
#F  ConjugateWyckoffPositions( W, trans ) . . . conj. Wyckoff pos. with basis
##
ConjugateWyckoffPositions := function( W, trans )

    local d, conj, invconj, w, gens, i;

    # prepare conjugators
    d := Length( trans );
    conj := IdentityMat( d+1 );
    conj{[1..d]}{[1..d]} := trans;
    invconj := conj^-1;

    # do the conjugation
    for w in W do
        w.basis       := w.basis*trans;
        w.translation := w.translation*trans;
        gens := w.stabilizer.generators;
        for i in [1..Length(gens)] do
            gens[i] := invconj*gens[i]*conj;
        od;
    od;

end;


#############################################################################
##
#F  WyckoffOrbit( W )  . . . . . . .  space group orbit of a Wyckoff position
##
WyckoffOrbit := function( W )
    local  S, d, T, gen, orb, u, g, img;
    
    if not IsWyckoffPosition( W ) then
        Error("W must be a Wyckoff position");
    fi;

    S   := W.spaceGroup;
    d   := S.dimension-1;
    T   := TranslationsCrystGroup( S );
    gen := S.internalGenerators;
    gen := List( gen, g -> [ g, g{[1..d]}{[1..d]}, g[1+d]{[1..d]} ] );
    orb := [ W ];

    # conjugate the Wyckoff position if necessary
    if not S.isStandard then
        ConjugateWyckoffPositions( orb, T^-1 );
    fi;

    for u in orb do
        for g in gen  do
            img := rec();
            img.translation := List( u.translation*g[2]+g[3], 
                                     FractionModOne );
            img.basis := u.basis*g[2];
            if not InAffineSpaces( orb, img ) then
                img.class      := u.class;
                img.stabilizer := S.operations.Subgroup( Parent( S ),
                    List( u.stabilizer.generators, x -> g[1]^-1*x*g[1] ) ); 
                img.operations := WyckoffPositionOps;
                img.spaceGroup := S;
                img.isWyckoffPosition := true;
                Add( orb, img );
            fi;
        od;
    od;

    # undo the conjugation
    if not S.isStandard then
        ConjugateWyckoffPositions( orb, T );
    fi;

    return orb;

end;


#############################################################################
##
#F  WyPos( S, A, r, h, pos ) . . . . . . . . . . . . . . . .Wyckoff positions
##
WyPos := function( S, A, r, h, pos )

    local d, subsizes, subgens, W, i, j, n, norm, cs, w;

    d := S.dimension - 1;

    # get subgroup generators in S
    subgens  := List( r, g -> List( g.generators, 
                              x -> PreImagesRepresentative( h, x ) ) );

    # the special positions
    W := Concatenation( [pos], List( [2..Length(r)], 
               i -> Set( FixedPointsModZ( subgens[i], d ) ) ) );

    # eliminate multiple copies of the same space; keep the last,
    # which is the one with the biggest stabilizer
    for i in [2..Length(r)-1] do
        for j in [i+1..Length(r)] do
            if Length( W[i] ) > 0 and Length( W[j] ) > 0 and
                Length( W[i][1].basis ) = Length( W[j][1].basis ) then
                W[i] := Filtered( W[i], x -> not InAffineSpaces( W[j], x ) );
            fi;
        od;
    od;

    # get one representative per orbit; loop over conjugacy classes
    for i in [2..Length(r)] do

        if Length( W[i] ) > 0 then

            # normalizer of subgroup in A
            n    := Normalizer( A, r[i] );
            norm := List( RightCosets( n, r[i] ), x ->
                     PreImagesRepresentative( h, Representative( x ) ) );

            # identify positions in the same orbit
            if Length(norm)>1 then
                W[i] := IdentifyWyckoffOrbits( W[i], norm );
                W[i] := CheckWyckoffOrbits( W[i], norm );
            fi;

            # filter positions with the right orbit length
            if Length(W[i])>0 then
                cs := List( RightCosets( A, n ), x ->
                       PreImagesRepresentative( h, Representative( x ) ) );
                W[i]  := CheckWyckoffOrbits( W[i], cs );
            fi;

            # adjust the stabilizers
            for w in W[i] do
                AdjustStabilizer( S, w, subgens[i], Size( r[i] ) );
                w.class := i;
                Add( pos, w );
            od;

        fi;

    od;

end; 


#############################################################################
##
#F  WyckoffPositions( S ) . . . . . . . . . . . . . . . . . Wyckoff positions
##
WyckoffPositions := function( S )

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if not IsBound( S.wyckoffPositions ) then
        S.wyckoffPositions := S.operations.WyckoffPositions( S );
    fi;
    return S.wyckoffPositions;

end;


#############################################################################
##
#F  CrystGroupOps.WyckoffPositions( S ) . . . . . . . . . . Wyckoff positions
##
CrystGroupOps.WyckoffPositions := function( S )

    local T, G, P, d, pos, F, A, Si, h, r, w;

    # check if we have indeed a space group
    T := TranslationsCrystGroup( S );
    if not S.isSpaceGroup then
        Error("S must be a space group");
    fi;

    # get point group G, and its permutation representation P
    G := PointGroup( S );
    P := PermGroup( G );
    d := S.dimension - 1;

    # first the general position
    pos := [ rec( basis       := IdentityMat( d ), 
                  translation := Flat( NullMat( 1, d ) ),
                  stabilizer  := S.operations.Subgroup( Parent( S ), [] ),
                  class       := 1,
                  spaceGroup  := S,
                  operations  := WyckoffPositionOps,
                  isWyckoffPosition := true ) ];

    # catch the trivial case
    if Size( P ) = 1 then 
        return pos;
    fi;

    # set up homomorphism from std representation to PermGroup or AgGroup
    if IsSolvable( P ) then     # switch to an AgGroup
        if not IsBound( P.agGroup ) then 
            F := FpGroup( G );
        fi;
        A := P.agGroup;
    else
        A := P;
    fi;
    Si := Group( S.internalGenerators, S.identity );
    h  := GroupHomomorphismByImages( Si, A, Si.generators, A.generators ); 
    r := List( ConjugacyClassesSubgroups(A), Representative );

    # now get the remaining Wyckoff positions
    WyPos( S, A, r, h, pos );

    # conjugate Wyckoff postions if necessary
    if not S.isStandard then
        ConjugateWyckoffPositions( pos, T );
    fi;

    # make the elements in pos WyckoffPositions
    for w in pos do
        w.operations        := WyckoffPositionOps;
        w.spaceGroup        := S;
        w.isWyckoffPosition := true;
    od;

    return pos;

end;          


#############################################################################
##
#F  WyckoffPositionsByStabilizer( S, stabs ) . . Wyckoff pos. for given stabs 
##
WyckoffPositionsByStabilizer := function( arg )

    local S, stabs, T, G, P, d, pos, F, A, Si, h, w;

    # check the arguments
    S := arg[1];
    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if IsGroup( arg[2] ) then
        stabs := [ arg[2] ];
    else
        stabs := arg[2];
    fi;
    
    # check if we have indeed a space group
    T := TranslationsCrystGroup( S );
    if not S.isSpaceGroup then
        Error("S must be a space group");
    fi;

    # get point group G, and its permutation representation P
    G := PointGroup( S );
    P := PermGroup( G );
    d := S.dimension - 1;

    # set up homomorphism from std representation to PermGroup or AgGroup
    stabs := List( stabs, x -> PreImage( P.bijection, x ) );
    if IsSolvable( P ) then     # switch to an AgGroup
        if not IsBound( P.agGroup ) then 
            F := FpGroup( G );
        fi;
        A := P.agGroup;
        stabs := List( stabs, x -> PreImage( A.bijection, x ) );
    else
        A := P;
    fi;

    # we have to catch the trivial stabilizer
    Sort( stabs, function(a,b) return Size(a) < Size(b); end );
    if Size( stabs[1] ) = 1 then
        pos := [ rec( 
                   basis       := IdentityMat( d ), 
                   translation := Flat( NullMat( 1, d ) ),
                   stabilizer  := S.operations.Subgroup( Parent( S ), [] ),
                   class       := 1 ) 
               ];
        stabs := Filtered( stabs, x -> Size( x ) > 1 );
    else
        pos := [];
    fi;
    stabs := Concatenation( [ TrivialSubgroup( A ) ], stabs );

    Si := Group( S.internalGenerators, S.identity );
    h  := GroupHomomorphismByImages( Si, A, Si.generators, A.generators ); 

    # now get the remaining Wyckoff positions
    WyPos( S, A, stabs, h, pos );

    # conjugate Wyckoff postions if necessary
    if not S.isStandard then
        ConjugateWyckoffPositions( pos, T );
    fi;

    # make the elements in pos WyckoffPositions
    for w in pos do
        w.spaceGroup        := S;
        w.operations        := WyckoffPositionOps;
        w.isWyckoffPosition := true;
    od;

    # return the result without storing it
    return pos;

end;          


#############################################################################
##
#F  WyckoffPositionsQClass( S, S1 ) . . . . . .  Wyckoff positions for QClass
##
WyckoffPositionsQClass := function( S, S1 )

    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if not IsBound( S.wyckoffPositions ) then
        S.wyckoffPositions := S.operations.WyckoffPositionsQClass( S, S1 );
    fi;
    return S.wyckoffPositions;

end;


#############################################################################
##
#F  CrystGroupOps.WyckoffPositionsQClass( S, S1 ) . . Wyckoff pos. for QClass
##
CrystGroupOps.WyckoffPositionsQClass := function( S, S1 )

    local T, P, d, pos, A, P1, gen, map, Si, h, r, w;

    # check if we have indeed a space group
    T := TranslationsCrystGroup( S );
    if not S.isSpaceGroup then
        Error("S must be a space group");
    fi;

    # get point group and its permutation representation P
    P := PermGroup( PointGroup( S ) );
    d := S.dimension - 1;

    # first the general position
    pos := [ rec( basis       := IdentityMat( d ), 
                  translation := Flat( NullMat( 1, d ) ),
                  stabilizer  := S.operations.Subgroup( Parent( S ), [] ),
                  class       := 1,
                  spaceGroup  := S,
                  operations  := WyckoffPositionOps,
                  isWyckoffPosition := true ) ];

    # catch the trivial case
    if Size( P ) = 1 then 
        return pos;
    fi;

    # set up homomorphism from std representation to PermGroup or AgGroup;
    # we use the PermGroup or AgGroup determined for S1
    if IsSolvable( P ) then
        A   := S1.pointGroup.agGroup;
        gen := List( A.generators, x -> Image( A.bijection, x ) );
    else
        A := S1.pointGroup.permGroup;
        gen := A.generators;
    fi;

    P1  := PermGroup( PointGroup( S1 ) );
    map := GroupHomomorphismByImages( P1, P,
                P1.generators, P.generators);
    gen := List( gen, x -> Image( map, x ) );
    gen := List( gen, x -> Image( P.bijection, x ) );
    gen := List( gen, x -> PreImagesRepresentative( S.pointHomom, x ) );
    if not S.isStandard then
       gen := List( gen, x -> x^T );
    fi;

    Si  := Group( gen, S.identity );
    h   := GroupHomomorphismByImages( Si, A, gen, A.generators ); 
    r   := List( ConjugacyClassesSubgroups(A), Representative );

    # now get the remaining Wyckoff positions
    WyPos( S, A, r, h, pos );

    # conjugate Wyckoff postions if necessary
    if not S.isStandard then
        ConjugateWyckoffPositions( pos, T );
    fi;

    # make the elements in pos WyckoffPositions
    for w in pos do
        w.operations        := WyckoffPositionOps;
        w.isWyckoffPosition := true;
    od;

    return pos;

end;          



