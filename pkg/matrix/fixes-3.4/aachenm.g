#############################################################################
##
#W  aachenm.g		   	Matrix Packages                  Frank Celler
##
#H  @(#)$Id: aachenm.g,v 1.1 1997/03/10 13:51:38 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the standard meataxe interface for the aachen meataxe.
##
RevisionMatrix.aachenm_g :=
    "@(#)$Id: aachenm.g,v 1.1 1997/03/10 13:51:38 gap Exp $";


if not IsBound(fail)  then
    fail := false;
fi;


#############################################################################
##

#F  MeatAxeModuleByGroup( <grp>, <field> )
##
MeatAxeModuleByGroup := function( arg )
    local   grp,  field,  basis,  gens,  dim;

    grp   := arg[1];
    field := arg[2];
    if Length(arg) = 3  then
        basis := arg[3];
    fi;

    # get the generators
    if IsGroup(grp)  then
        gens  := grp.generators;
        dim   := Length(grp.identity);
        if not IsBound(basis)  then
            basis := grp.identity;
        fi;
    elif IsGModule(grp)  then
        gens := GeneratorsFlag(grp);
        dim  := DimensionFlag(grp);
        if not IsBound(basis)  then
            if 0 = Length(gens)  then
                Error( "basis is unknown" );
            fi;
            basis := gens[1]^0;
        fi;
    elif IsList(grp)  then
        gens := grp;
        if 0 = Length(gens)  then
            Error( "dimension is unknown" );
        fi;
        dim := Length(gens[1]);
        if not IsBound(basis)  then
            basis := gens[1]^0;
        fi;
    fi;

    # set up the meataxe module record
    return rec( field := field, generators := gens, dimension := dim,
                basis := basis );

end;


#############################################################################
##

#V  AachenMeatAxe
##
AachenMeatAxe := rec();


#############################################################################
##
#F  AachenMeatAxeInfoFlag( <module> )
##
AachenMeatAxeInfoFlag := function( module )
    local   subs,  basis,  aachen,  mat,  m,  b;

    # if the info is already there return it
    if IsBound(module.aachen_meataxe)  then
        return module.aachen_meataxe;
    fi;

    # construct the action
    subs  := RowSpace( module.field, module.basis );
    basis := Basis( subs, module.basis );
    if basis.vectors <> module.basis  then
        Error( "basis mismatch" );
    fi;
        
    # set up the matrices
    aachen := rec();
    aachen.generators := [];
    for mat  in module.generators  do
        m := [];
        for b  in module.basis  do
            Add( m, Coefficients( basis, b * mat ) );
        od;
        DisplayMat(m);
        Print("\n");
        Add( aachen.generators, MeatAxeMat( m, module.field ) );
    od;

    # create an algebra
    aachen.algebra := UnitalAlgebra( module.field, aachen.generators );

    # and the natural module
    aachen.module := NaturalModule( aachen.algebra );

    # that's it
    module.aachen_meataxe := aachen;
    return module.aachen_meataxe;
    
end;


#############################################################################
##

#F  AachenMeatAxe.IsIrreducible( <module> )
##
AachenMeatAxe.IsIrreducible := function( module )
    local   aachen;

    # check if we know if the module is irreducible
    if IsBound(module.IsIrreducible)  then
        return module.IsIrreducible;
    fi;

    # get the info record
    aachen := AachenMeatAxeInfoFlag(module);

    # check for irreducibility
    module.IsIrreducible := IsIrreducible(aachen.module);

    # and return the result
    return module.IsIrreducible;
    
end;


#############################################################################
##
#F  AachenMeatAxe.IsAbsolutelyIrreducible( <module> )
##
AachenMeatAxe.IsAbsolutelyIrreducible := function( module )
    local   aachen;

    # check if we know if the module is abs irreducible
    if IsBound(module.IsAbsolutelyIrreducible)  then
        return module.IsAbsolutelyIrreducible;
    fi;

    # get the info record
    aachen := AachenMeatAxeInfoFlag(module);

    # check for abs irreducibility
    module.IsAbsolutelyIrreducible := IsAbsolutelyIrreducible(aachen.module);

    # and return the result
    return module.IsAbsolutelyIrreducible;
    
end;


#############################################################################
##
#F  AachenMeatAxe.DegreeSplittingField( <module> )
##
AachenMeatAxe.DegreeSplittingField := function( module )
    local   aachen,  field;

    # check if we know the degree
    if IsBound(module.DegreeSplittingField)  then
        return module.DegreeSplittingField;
    fi;

    # get the info record
    aachen := AachenMeatAxeInfoFlag(module);

    # compute the splitting field
    field := SplittingField(aachen.module);
    module.DegreeSplittingField := LogInt( Size(field),
                                           Characteristic(field) );

    # and return the result
    return module.DegreeSplittingField;
    
end;


#############################################################################
##
#F  AachenMeatAxe.ProperSubmodule( <module> )
##
AachenMeatAxe.ProperSubmodule := function( module )
    local   aachen,  modu,  subs,  basis;

    # check if we know a proper submodule
    if IsBound(module.ProperSubmodule)  then
        return module.ProperSubmodule;
    fi;

    # get the info record
    aachen := AachenMeatAxeInfoFlag(module);

    # compute an invariant subspace
    modu := [ List( aachen.generators, GapObject ), module.field,
              module.dimension ];
    subs := MTXOps.InvariantSubspace(modu);

    # if we have found a module construct a new description
    if IsInt(subs)  then
        subs := fail;
    else
        subs := rec( field := module.field, generators := module.generators,
                     basis := Basis(subs).vectors*module.basis );
        subs.dimension := Length(subs.basis);
    fi;
    module.ProperSubmodule := subs;

    # and return the result
    return module.ProperSubmodule;
    
end;


#############################################################################
##
#F  AachenMeatAxe.MinimalSubmodule( <module> )
##
AachenMeatAxe.MinimalSubmodule := function( module )
    local   aachen,  modu,  subs,  basis,  m;

    # check if we know a minimal submodule
    if IsBound(module.MinimalSubmodule)  then
        return module.MinimalSubmodule;
    fi;

    # get the info record
    aachen := AachenMeatAxeInfoFlag(module);

    # compute a proper submodule
    modu := [ module.generators, module.field, module.dimension ];
    subs := MTXOps.InvariantSubspace(modu);
    if IsInt(subs)  then
        module.MinimalSubmodule := fail;
        return module.MinimalSubmodule;
    else
        m := rec( field := module.field, generators := module.generators,
                  basis := Basis(subs).vectors );
        m.dimension := Length(m.basis);
    fi;

    # compute a minimal submodule
    repeat
        Print( "dimension: ", m.dimension, "\n" );
        modu := [ m.generators, m.field, m.dimension ];
        subs := MTXOps.InvariantSubspace(modu);
        if IsInt(subs)  then
            module.MinimalSubmodule := m;
            return module.MinimalSubmodule;
        else
            m := rec( field := m.field, generators := m.generators,
                      basis := Basis(subs).vectors );
            m.dimension := Length(m.basis);
        fi;
    until m.dimension = 1;

    # the submodule is one-dimensional, so it must be minimal
    module.MinimalSubmodule := m;
    return module.MinimalSubmodule;
    
end;


#############################################################################
##
#F  AachenMeatAxe.InducedActionSubmodule( <module>, <basis> )
##
AachenMeatAxe.InducedActionSubmodule := function( module, basis )

    # if the module has no generators there is nothing to do
    if 0 = Length(module.generators)  then
        return rec( field := module.field, dimension := module.dimension,
                    basis := IdentityMat( module.dimension, module.field ),
                    generators := [] );
    fi;

    # construct the action
    vs := VectorSpace( basis
end;


#############################################################################
##

#E  aachenm.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
