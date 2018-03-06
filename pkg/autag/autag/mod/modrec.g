###########################################################################
##
#A  modrec.g                 autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing decompositions of
##  modules into indecomposable summands, as well as computing generating
##  sets for module automorphisms.
##
##  It is an integral part of the soluble group automorphism package.
##
###########################################################################


# module records
#
# This file contains functions for storing and retrieving information
# about finite dimensional modules.
#
# It is part of a package of GAP files written for the computation of
# automorphism groups of group given by special AG-presentations.


ReadLocalDir(LOCALNAMEmod, "matrices");
ReadLocalDir(LOCALNAMEmod, "modaut");
ReadLocalDir(LOCALNAMEmod, "modhom");
ReadLocalDir(LOCALNAMEmod, "modiso");


Gmodule := function ( arg )
    local matrices, F, r;

    if Number(arg) = 1  then
        if IsGroup(arg[1]) then 
            matrices := arg[1].generators;
        else 
            matrices := arg[1];
        fi;
        F:=Field(matrices[1][1][1]);
        # NB: this means that for prime power fields, must pass field
    elif Number(arg) = 2  then
        if IsGroup(arg[1]) then 
            matrices := arg[1].generators;
        else 
            matrices := arg[1];
        fi;
        F := arg[2];
    else
        return Error("usage: Gmodule ( <matrices> [, <F>] )");
    fi;

    r := rec( field := F,
              dim := Length(matrices[1]),
              dimension := Length(matrices[1]),
              generators := matrices,
              matrices := matrices, # for smash
              isGModule := true,
              operations := GModRecOps
              );
    return r;
end;


GModOps.IsGModule := function ( module )
    if IsRec(module)=false or
       IsBound(module.isGModule)=false then
        return false;
    fi;
    return module.isGModule;
end;


# ---------------------------------------------------------------------------

GModOps.Generators := function ( module )
    return module.generators;
end;

# ---------------------------------------------------------------------------

GModOps.Field := function ( module )
    if IsBound(module.field)=false then return "unknown"; fi;
    return module.field;
end;

# ---------------------------------------------------------------------------

GModOps.DimFlag := function ( module )
    if IsBound(module.dim)=false then return "unknown"; fi;
    return module.dim;
end;

# ---------------------------------------------------------------------------

GModOps.SetBasis := function ( module,basis )
    module.basis := basis;
end;

GModOps.Basis := function (M)
    #
    if not IsBound(M.basis) then
        return IdentityMat(GModOps.DimFlag(M), GModOps.Field(M));
    fi;
    return M.basis;
end;

GModOps.SetEndAlgBasisFlag := function ( module, base )
    module.endomorphisms := base;
end;

GModOps.EndAlgBasisFlag := function ( module )
    if IsBound(module.endomorphisms)=false then return "unknown"; fi;
    return module.endomorphisms;
end;

GModOps.SetIndecomposableFlag := function ( module, flag )
    module.isIndecomposable := flag;
end;

GModOps.IndecomposableFlag := function ( module )
    if IsBound(module.isIndecomposable)=false then return "unknown"; fi;
    return module.isIndecomposable;
end;

GModOps.SetIndecomposablesFlag := function (M, comps)
    M.indecomposables := comps;
end;

GModOps.IndecomposablesFlag := function (M)
    if IsBound(M.indecomposables) then
        return M.indecomposables;
    fi;
    return "unknown";
end;

GModOps.SetMultiplicitiesFlag := function (M, mults)
    M.multiplicities := mults;
end;

GModOps.MultiplicitiesFlag := function (M)
    if IsBound(M.multiplicities) then
        return M.multiplicities;
    fi;
    return "unknown";
end;


GModOps.SetEndAlgResidueFlag := function (M, root, ord)
    M.rootResidue := root;
    M.rootResidueOrder := ord;
end;

GModOps.SetEndAlgRadicalFlag := function (M, base)
    M.endoRadical := base;
end;

GModOps.EndAlgResidueFlag := function (M)
    if IsBound(M.rootResidue) then
        return M.rootResidue;
    fi;
    return "unknown";
end;

GModOps.EndAlgRadicalFlag := function (M)
    if IsBound(M.endoRadical) then
        return M.endoRadical;
    fi;
    return "unknown";
end;

GModOps.SetEchelonisationFlag := function (M, ech)
    #
    M.echelonIndices := ech;
end;


# ---------------------------------------------------------------------------
#
# Gmodule operations

GModRecOps.Print := function ( r )
    local i, d;
    Print("Gmodule(");
    Print("\n  field := ", r.field);
    Print("\n  dim := ", r.dim);
    if IsBound(r.generators) then
        Print("\n  generators : ",Length(r.generators)," generator",
              Plural(Length(r.generators)));
    fi;
    if IsBound(r.endomorphisms) then
        Print("\n  endomorphisms : ", Length(r.endomorphisms), 
              " basis element", Plural(Length(r.endomorphisms)));
    fi;
    if IsBound(r.isIndecomposable) then
        Print("\n  indecomposable := ", r.isIndecomposable);
    fi;
    if IsBound(r.indecomposables) then
        Print("\n    indecomposables : (dim)  ");
        for i in [1..Length(r.indecomposables)] do
            d := GModOps.DimFlag(r.indecomposables[i]);
            if d < 10 then Print(" "); fi;
            Print(d, "  ");
        od;
    fi;
    if IsBound(r.multiplicities) then
        Print("\n    multiplicities :         ");
        for i in [1..Length(r.multiplicities)] do
            d := r.multiplicities[i];
            if d < 10 then Print(" "); fi;
            Print(d, "  ");
        od;
    fi;
    if IsBound(r.irreducible) then
        Print("\n  irreducible := ", r.irreducible);
    fi;
    if IsBound(r.basis) then
        Print("\n  basis : defined");
    fi;
    if IsBound(r.groupgens) then
        Print("\n  gens : group generators are defined");
    fi;
    if IsBound(r.rep) then
        Print("\n  rep : representation");
    fi;
    if IsBound(r.rootResidue) then
        Print("\n  rootResidue: primitive root of residue field");
        Print(" (size ", r.rootResidueOrder + 1, ")");
    fi;
    if IsBound(r.endoRadical) then
        Print("\n  endoRadical: ",
              " basis for radical of endomorphism algebra, dim ",
              Length(r.endoRadical));
    fi;
    if IsBound(r.auts) then
        Print("\n  auts : ", Length(r.auts), " generators for ",
              "module automorphism group");
    fi;
    if IsBound(r.autorder) then
        Print("\n  autorder : module automorphism group has order ", 
              r.autorder);
    fi;
    Print(")\n");
end;


############################################################################
##
#F GModRecOps.\^ ( l, r )
##
GModRecOps.\^ := function ( l, r )

    # change basis on module via matrix r

    local rinv, M;
    
    rinv := r^-1;
    
    if (not GModOps.IsGModule(l)) then
        Error("LHS must be Gmodule");
    elif not IsMat(r) then
        Error("RHS must be matrix");
    elif Length(r) <> GModOps.DimFlag(l) then
        Error("dimensions must match up");
    elif GModOps.Field(l) <> Field(r[1][1]) then
        Error("must be over same field");
    elif DeterminantMat(r) = GModOps.Field(l).zero then
        Error("singular matrix used for basis change");
    fi;

    M := Gmodule(List(GModOps.Generators(l), x -> x^rinv));

    M.generators := List(GModOps.Generators(l), x -> x^rinv);
    
    GModOps.SetBasis(M, GModOps.Basis(l)*r);

    if GModOps.EndAlgBasisFlag(l) <> "unknown" then
        GModOps.SetEndAlgBasisFlag(M, 
                List(GModOps.EndAlgBasisFlag(l), x -> x^rinv));
    fi;

    if GModOps.EndAlgRadicalFlag(l) <> "unknown" then
        GModOps.SetEndAlgRadicalFlag(M, 
                List(GModOps.EndAlgRadicalFlag(l), x -> x^rinv));
    fi;
    
    if GModOps.EndAlgResidueFlag(l) <> "unknown" then
        GModOps.SetEndAlgResidueFlag(M, GModOps.EndAlgResidueFlag(l)^rinv);
    fi;    

    if GModOps.IndecomposablesFlag(l) <> "unknown" then
        GModOps.SetIndecomposablesFlag(M, GModOps.IndecomposablesFlag(l));
    fi;
    
    if GModOps.MultiplicitiesFlag(l) <> "unknown" then
        GModOps.SetMultiplicitiesFlag(M, GModOps.MultiplicitiesFlag(l));
    fi;
        
    return M;
end;


############################################################################
##
#F GModRecOps.\= (M, N)
##
GModRecOps.\= := function (M, N)
    return not GModOps.IsomModules(M,N) = false;
end;



## Local Variables:
## mode:gap
## mode:outline-minor
## End:
