###########################################################################
##
#A  modiso.g                 autag package                 Michael J Smith
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


# module isomorphism and decomposition routines
#
# This file contains functions for computing with modules, including:
#
#   (1) computing a direct sum decomposition of a module into
#   indecomposable summands.
#
#   (2) deciding module isomorphism using the decomposition.
#
# The algorithm for deciding indecomposability is based on the algorithm
# described by G. Schneider in the Journal of Symbolic Computation, 
# Volume 9, Numbers 5 & 6, 1990


GModOps.ModuleHomBasis := function (M1, M2)
    #
    # Compute a basis for the vector space of module homomorphisms from
    # M1 to M2
    #
    # If the dimension of M1 is not too small, call a new spinning based
    # algorithm, otherwise just use the simple intertwining equations one.
    #
    # Returns an echelonised basis.
    
    local basis;

    if GModOps.DimFlag(M1) > 5 then
        basis := SpinHom(M1, M2);
    else
        basis := GModOps.Intertwine(M1, M2);
    fi;
    return basis;
end;





GModOps.ChopModule := function (M)
    #
    # Take a module and break it into two pieces if possible.
    #
    # The function searches for a decomposition of the module M while
    # attempting to prove indecomposability at the same time.  Of course,
    # only one of these will succeed.
    #
    local proveIndecomposability, addnilpotent, n, F, echelon, basis, 
          enddim, nildim, p, maxorder, maxa, nilbase, nilech, null, i, 
          remain, used, coeffs, a, order, pos, newa, lastdim;
    
    proveIndecomposability := function ()
        #
        # Check whether we have found the indecomposability proof. That is,
        # see whether our regular element generates a subalgebra which
        # complements the current nilpotent ideal (the approximation to
        # radical)
        #
        local maxaord;
        # NB: <maxa> is not local

        if enddim - nildim = LogInt(maxorder + 1,p) then
            #
            # Yes, found the residue field root and proved indecomposability!
            maxaord := OrderMat(maxa);
            while maxaord > maxorder do
                maxa := maxa^p;
                maxaord := maxaord / p;
            od;
            GModOps.SetEndAlgResidueFlag(M, maxa, maxaord);
            GModOps.SetEndAlgRadicalFlag(M, nilbase);
            GModOps.SetIndecomposableFlag(M, true);

            return true;
        fi;
        return false;
    end;

    addnilpotent := function (a)
        #
        # take a new nilpotent element and sift against current nilpotent
        # ideal basis. If it does not lie in the space spanned so far,
        # add it to nilbasis and.

        local i, r, c, k, done, l;
        # NB: <remain> and <nildim> are not local

        for i in [1..nildim] do
            r := echelon[nilech[i]][1]; c := echelon[nilech[i]][2];
            if a[r][c] <> F.zero then
                a := a - a[r][c] * nilbase[i] / nilbase[i][r][c];
            fi;
        od;
        
        # find which echelon index to remove due to this new element
        k := 1; done := false;
        while not done and k <= Length(remain) do
            l := remain[k];
            r := echelon[l][1]; c := echelon[l][2];
            if a[r][c] <> F.zero then
                done := true;
            else
                k := k + 1;
            fi;
        od;

        if k > Length(remain) then
            # in nilpotent ideal already, return
            return false;
        fi;
        
        # We now know this nilpotent element is a new one
        #
        Add(nilbase, a);

        # the k-th basis element was used to make the new element a. So
        # remove it from future random element calculations
        #
        Add(nilech, remain[k]);
        remain := Difference(remain, [remain[k]]);
        nildim := nildim + 1;
        i := 1;
        return true;
    end;

    n := GModOps.DimFlag(M);
    F := GModOps.Field(M);
    InfoDecompose( "\nchopModule called with module of dimension ", n, "\n");

    if n = 1 then
        #
        # A 1-dimensional module is always indecomposable
        #
        GModOps.SetIndecomposableFlag(M, true);
        GModOps.SetEndAlgResidueFlag(M, [[ F.root ]], F.size - 1);
        GModOps.SetEndAlgRadicalFlag(M, []);
        return [ M ];
    fi;

    if Length(GModOps.EndAlgBasis(M)) = 1 then
        #
        # if endomorphism algebra has dimension 1 then indecomposable
        #
        GModOps.SetIndecomposableFlag(M, true);
        GModOps.SetEndAlgResidueFlag(M, F.root * GModOps.EndAlgBasisFlag(M)[1],
                F.size - 1);
        GModOps.SetEndAlgRadicalFlag(M, []);
        return [ M ];
    fi;

    echelon := GModOps.Echelonisation(M); # echelon indices for endalg basis
    basis := GModOps.EndAlgBasis(M);
    enddim := Length(basis);            # dim of endo algebra
    nildim := 0;                        # dim of current approx to radical
    p := Size(F);
    maxorder := 1;                      # order of largest order regular elmt
                                        #   found so far
    maxa := IdentityMat(n,F);           # the regular elmt with order maxorder
    nilbase := [];                      # basis for approx to radical
    nilech := [];
    
    null := NullMat(n,n,F);

    i := 1;

    # We will "quotient" out the nilpotent subspace as we go. The elements
    # of remain tell us which (echelonised) basis elements of the
    # endomorphism algebra we will take use in our random linear
    # combination.
    #
    remain := [1..enddim];
    used := [];

    # we will loop until too many passes without an improvement in knowledge
    repeat
        #
        # randomly sample endomorphism algebra
        # 
        repeat
            coeffs := List([1..enddim], x -> Random(F));
        until coeffs{remain} <> 0*coeffs{remain};

        a := GModOps.DotProduct(coeffs, basis);

        if DeterminantMat(a) <> F.zero then
            #
            # a regular element, check to see whether its order is
            # larger than previously known, and if so whether it
            # generates the residue field modulo current nilpotent ideal
            #
            order := OrderMat(a);

            while (order mod p = 0) do
                order := order / p;
            od;
            if order > maxorder then
                maxorder := order;
                maxa := a;
                if proveIndecomposability() then
                    return [ M ];
                fi;
                i := 1;
            else
                i := i + 1;
            fi;

        elif GModOps.IsFittingMat(a) <> false then

            return GModOps.SplitModule(M,a);

        elif addnilpotent(a) then
            #
            # new nilpotent element, added to nilbasis. Now close nilbasis to
            # basis for an ideal.

            # keep a pointer to the first new element added to nilbase
            #
            pos := nildim; # a was just added

            # first add powers of a
            #
            newa := a^2;
            repeat
                lastdim := nildim;
                addnilpotent(newa);
                newa := newa * a;
            until lastdim = nildim or newa = null;
            
            # now close nilbase to make ideal basis
            #
            repeat
                for i in [1..enddim] do
                    a := nilbase[pos] * basis[i];
                    if GModOps.IsFittingMat(a) <> false then
                        return GModOps.SplitModule(M,a);
                    fi;
                    addnilpotent(a);
                od;
                pos := pos + 1;
            until pos = nildim + 1;

        fi;

        if proveIndecomposability() then
            return [ M ];
        fi;

    until (i >= GModOps.ChopLimit);
        
    Error(Cat(
            "Unable to ascertain module decomposition within time limits.\n",
            "Increase the value of GModOps.ChopLimit ",
            "and try again.\n"));
    
end;

GModOps.ChopLimit := 50;


GModOps.SplitModule := function (M, a)
    #
    # Take a Fitting element and use it to split M into a direct sum
    # of submodules. Return the submodules.
    #
    local runtime, t, n, F, p, r, rt, B, Bgens, Bendo, M1, M2;

    runtime := Runtime();

    t := GModOps.IsFittingMat(a);
    if not IsMatrix(a) then
        Error("Cannot split module with non-Fitting element\n");
    fi;

    n := GModOps.DimFlag(M);
    F := GModOps.Field(M);
    p := Size(GModOps.Field(M));

    r := RankMat(t);
    if (r >= n) then
        Error("Oops! Obviously the element found wasn't a Fitting element.\n");
    fi;

    InfoDecompose( "Decomposition found into two submodules, dimensions ");
    InfoDecompose( r,", ",n-r,".\n");

    rt := Runtime();
    B := Cat(BaseMat(t), NullspaceMat(t));
    InfoDecompose("   B...",Seconds(Runtime()-rt),"\n");
    rt := Runtime();
    Bgens := B * GModOps.Generators(M) * B^-1;
    if Length(GModOps.EndAlgBasisFlag(M)) > 0 then
        Bendo := B * GModOps.EndAlgBasisFlag(M) * B^-1;
    else
        Bendo := [];
    fi;
    InfoDecompose("   Bgens, Bendo...",Seconds(Runtime()-rt),"\n");

    GModOps.SetIndecomposableFlag(M, false);
    GModOps.SetBasis(M, B);

    rt := Runtime();
    M1 := Gmodule(List(Bgens, x -> GModOps.ExtractMat(x,1,1,r)),
                  GModOps.Field(M));

    M2 := Gmodule(List(Bgens, x -> GModOps.ExtractMat(x,r+1,r+1,n-r)),
                  GModOps.Field(M));
    InfoDecompose("   M1, M2...",Seconds(Runtime()-rt),"\n");


    # Extract the endomorphisms for the new submodules from the blocks
    # of M's endomorphisms. Of course, we will not have linear independence
    # any more, so call the echelonisation routine immediately to make
    # sure the set is trimmed down to a basis.
    #
    rt := Runtime();
    GModOps.SetEndAlgBasisFlag(M1, List(Bendo, x -> 
            GModOps.ExtractMat(x,1,1,r)));
    InfoDecompose("   EndAlg(M1)...",Seconds(Runtime()-rt),"\n");
    rt := Runtime();
    GModOps.EchEndAlg(M1);
    InfoDecompose("   ech(EndAlg(M1))...",Seconds(Runtime()-rt),"\n");

    rt := Runtime();
    GModOps.SetEndAlgBasisFlag(M2, List(Bendo, 
            x -> GModOps.ExtractMat(x,r+1,r+1,n-r)));
    InfoDecompose("   EndAlg(M2)...",Seconds(Runtime()-rt),"\n");
    rt := Runtime();
    GModOps.EchEndAlg(M2);
    InfoDecompose("   ech(EndAlg(M2))...",Seconds(Runtime()-rt),"\n");

    InfoDecompose("GModOps.SplitModule: runtime ",Seconds(Runtime()-runtime), 
            "\n");

    # Return the decomposition
    #
    return [ M1, M2 ];
end;




GModOps.ComputeDecomp := function (M)

    local parts;

    if GModOps.IndecomposableFlag(M) <> true then
        parts := GModOps.ComputeDecompRecurse(M);
        if GModOps.IndecomposableFlag(M) = false then
            GModOps.ParseDecomp(M, parts);
        fi;
    fi;

end;



GModOps.ComputeDecompRecurse := function (M)
    #
    # Take a module and try to find its decomposition into indecomposables.
    #
    local n, F, decomp, parts, M1, n1, M2, n2, parts1, parts2, B;

    n := GModOps.DimFlag(M);
    F := GModOps.Field(M);
    decomp := GModOps.ChopModule(M);

    if Length(decomp) = 1 then

        parts := [ decomp[1] ];

    else

        M1 := decomp[1]; n1 := GModOps.DimFlag(M1);
        M2 := decomp[2]; n2 := GModOps.DimFlag(M2);

        parts1 := GModOps.ComputeDecompRecurse(M1);
        parts2 := GModOps.ComputeDecompRecurse(M2);

        B := IdentityMat(n, GModOps.Field(M));

        GModOps.InsertMat( B, GModOps.Basis(M1), 1, 1 );
        GModOps.InsertMat( B, GModOps.Basis(M2), 1 + n1, 1 + n1);

        B := B * GModOps.Basis(M);

        GModOps.SetBasis(M, B);
        parts := Cat(parts1, parts2);

    fi;

    return parts;

end;



GModOps.ParseDecomp := function (M, parts)
    #
    # Examine the decomposition of M and determine the homogeneous
    # components of M. Will find a list of isomorphism class representatives
    # for the indecomposable summands of M and their corresponding
    # multiplicities. The basis flag of M is changed to exhibit this
    # structural information.
    # 
    # pass in previously computed list of submodules, parts.
    #
    local n, F, B, k, dims, index, i, remain, hcomps, hbases, thiscomp, 
          newremain, thisbases, j, b, B2, decomp, mults, pos, m, mats, dim;

    if (IsBound(M.multiplicities)) then
        #
        # Module decomposition already parsed
        
        return;
    fi;

    n := GModOps.DimFlag(M);
    F := GModOps.Field(M);
    B := GModOps.Basis(M);
    M.oldBasis := B;

    k := 1;
    dims := []; 
    index := [];
    for i in [1..Length(parts)] do
        index[i] := k;
        dims[i] := GModOps.DimFlag(parts[i]);
        k := k + dims[i];
    od;
    #
    # Now dims is a list of dimensions and index points to the start of
    # each block

    remain := [1..Length(parts)];    hcomps := [];
    hbases := [];
    while Length(remain) > 0 do
        thiscomp := [remain[1]];
        newremain := [];
        thisbases := [ IdentityMat(GModOps.DimFlag(parts[remain[1]]),F) ];
        for j in [2..Length(remain)] do
            b := GModOps.IsomIndecModules(
                         parts[remain[1]], parts[remain[j]]);
            if IsMat(b) then
                Add (thiscomp, remain[j]);
                Add (thisbases, b);
            else
                Add (newremain, remain[j]);
            fi;
        od;
        Add (hcomps, thiscomp);
        Add (hbases, thisbases);
        remain := newremain;
    od;

    B2 := NullMat(n,n,F);
    decomp := [];
    mults := [];
    pos := 1;
    for i in [1..Length(hcomps)] do
        m := GModOps.DimFlag(parts[hcomps[i][1]]);
        Add(decomp,  parts[hcomps[i][1]] );
        Add(mults, Length(hcomps[i]));
        for j in [1..Length(hcomps[i])] do
            GModOps.InsertMat(B2, hbases[i][j], pos, index[hcomps[i][j]]);
            pos := pos + m;
        od;
    od;

    B := B2 * B;

    GModOps.SetBasis(M,B);
    GModOps.SetIndecomposablesFlag(M, decomp);
    GModOps.SetMultiplicitiesFlag(M, mults);

    if IsBound(M.groupgens) then
        for i in [1..Length(M.indecomposables)] do
            M.indecomposables[i].groupgens := M.groupgens;
        od;
    fi;

    if IsBound(M.genimages) then
        mats := B * M.genimages * B^-1;
        for i in [1..Length(M.indecomposables)] do
            pos := GModOps.DecompIndex(M,i,1);
            dim := M.indecomposables[i].dim;
            M.indecomposables[i].genimages := 
              List(mats, x -> GModOps.ExtractMat(x, 
                      pos, pos, dim));
        od;
    fi;
    
end;


GModOps.IsomIndecModules := function (M1, M2)
    #
    # Check isomorphism of indecomposable modules.
    #
    # If they are isomorphic then the homomorphism space between them is a
    # disguised copy of the endomorphism algebra. This is a local algebra,
    # and hence all singular elements are nilpotent. Certainly it cannot
    # have a basis consisting entirely of nilpotent elements (a theorem of
    # Wedderburn), so at least one basis element for Hom(M1,M2) must be an
    # isomorphism if they are isomorphic.
    #
    local base, i;

    if not (GModOps.IsIndecomposable(M1) and GModOps.IsIndecomposable(M2)) then
        Error("GModOps.IsomIndecModules: requires indecomposable modules\n");
    fi;

    # module dimensions certainly must match
    if GModOps.DimFlag(M1) <> GModOps.DimFlag(M2) then return false; fi;

    # their endomorphism algebras must have same dimension
    if Length(GModOps.EndAlgBasis(M1)) <> 
       Length(GModOps.EndAlgBasis(M2)) then
        return false;
    fi;

    # radicals in endomorphisms algebra must have same dimension
    if Length(GModOps.EndAlgRadical(M1)) <>
       Length(GModOps.EndAlgRadical(M2)) then
        return false;
    fi;

    # the easy options have run out

    # Last case, both modules are idecomposable but not irreducible.
    # In this case, compute Hom and look for isom in the basis.

    base := GModOps.ModuleHomBasis(M1, M2);
    
    for i in [1..Length(base)] do
        if RankMat(base[i]) = GModOps.DimFlag(M1) then
            return base[i];
        fi;
    od;
    
    return false;

end;




GModOps.IsomModules := function (M1, M2)
    #
    # Test for isomorphism of modules. Will return one of:
    #
    # (1) the isomorphism as an F-matrix between M1 and M2
    # (2) false if the two modules are definitely not isomorphic
    #
    # Note that the isomorphism X is such that conjugating each generator
    # acting on M1 by X gives the corresponding action on M2. Therefore
    # X^-1 is a matrix whose rows correspond to a new basis of M1 that
    # duplicates the action of M2 on M1.
    #
    # If necessary, uses the decomposition into indecomposable summands.  A
    # homogeneous component is a direct sum of multiple copies of a single
    # indecomposable summand. The homogeneous components must match between
    # each module, with their multiplicities.

    local F, n, comps1, mults1, comps2, mults2, nc, X, remain, i, j, found, 
          b, k;

    F := GModOps.Field(M1);
    n := GModOps.DimFlag(M1);
    
    if F <> GModOps.Field(M2) then
        Error("Modules are over different fields.\n"); 

    elif Length(GModOps.Generators(M1)) <> Length(GModOps.Generators(M2)) then
        Error("Modules have different numbers of generators.\n"); 

    elif n <> GModOps.DimFlag(M2) then
        # Modules have different dimensions
        return false;

    elif Length(GModOps.EndAlgBasis(M1)) 
      <> Length(GModOps.EndAlgBasis(M2)) then
        # different endomorphism algebra dimensions
        return false;
    fi;

    comps1 := GModOps.Indecomposables(M1);   
    mults1 := GModOps.Multiplicities(M1); 

    comps2 := GModOps.Indecomposables(M2);   
    mults2 := GModOps.Multiplicities(M2);

    nc := Length(comps1);
    if nc <> Length(comps2) then
        # Modules have different number of homogeneous components
        return false;
    fi;
    
    # X shall accumulate the isomorphism as we check
    X := NullMat(n,n,F);

    remain := [1..nc];
    for i in [1..nc] do
        j := 1; found := false;
        while j <= nc and not found do
            if remain[j] <> 0 and mults1[i] = mults2[j] then
                # multiplicities match, now check isomorphism
                b := GModOps.IsomIndecModules(comps1[i], comps2[j]);
                if IsMat(b) then
                    found := true; remain[j] := 0;
                    for k in [1..mults1[i]] do
                        GModOps.InsertMat(X, b, 
                                GModOps.DecompIndex(M1,i,k),
                                GModOps.DecompIndex(M2,j,k)
                                );
                    od;
                fi;
            fi;
            j := j + 1;
        od;
        
        if not found then
            # No match in M2 for i-th component of M1
            return false;
        fi;
    od;

    return GModOps.Basis(M1)^-1 * X * GModOps.Basis(M2);
    #
    # Now gen(M1)^answer = gen(M2).

end;




GModOps.EchEndAlg := function (M)
    
    local ret;
    
    ret := GModOps.EcheloniseMats(GModOps.EndAlgBasis(M), GModOps.Field(M));
    
    GModOps.SetEndAlgBasisFlag(M, ret[1]);
    GModOps.SetEchelonisationFlag(M, ret[2]);

end;
                




GModOps.DecompIndex := function (M,comp,copy)
    #
    # Return the index of the first basis element of M.basis which
    # corresponds to the j-th copy of the indecomposable in the i-th
    # homogeneous component of M
    #
    local comps, mults, idx, s, i, d, j;
    
    comps := GModOps.Indecomposables(M);
    mults := GModOps.Multiplicities(M);

    if not IsBound(M.index) then
        # build the index to the starting positions of the submodules
        # and store it in M's record
        #
        idx := [];
        s := 1;
        for i in [1..Length(comps)] do
            idx[i] := [];
            d := GModOps.DimFlag(comps[i]);
            for j in [1..mults[i]] do
                Add(idx[i], s);
                s := s + d;
            od;
        od;
        M.index := idx;
    fi;
    
    return M.index[comp][copy];
    
end;




GModOps.DotProduct := function (arg)
    #
    # compute the linear combination of elements defined by the coefficient
    # vector coeffs. Optional permutation argument applies a permutation
    # to the elements.
    #
    local coeffs, elements, perm, n, zero, j, ans, i;
    
    coeffs := arg[1];
    elements := arg[2];
    if Length(arg) > 2 then
        perm := arg[3];
    else
        perm := [1..Length(elements)];
    fi;

    n := Length(coeffs);
    zero := coeffs[1]*0;
    j := 1;
    while coeffs[j] = zero and j <= n do
        j := j + 1;
    od;
    if j > n then
        return zero * elements[perm[1]];
    fi;
    
    ans := coeffs[j] * elements[perm[j]];
    for i in [j+1..n] do
        if coeffs[i] <> zero then
            ans := ans + coeffs[i] * elements[perm[i]];
        fi;
    od;

    return ans;
end;



#---------------------------------------------------------------------------
#
# record maintenance routines


GModOps.Indecomposables := function (M)
    if GModOps.IsIndecomposable(M) then
        return [ M ];
    else
        return GModOps.IndecomposablesFlag(M);
    fi;
end;

GModOps.Multiplicities := function (M)
    if GModOps.IsIndecomposable(M) then
        return [ 1 ];
    else
        return GModOps.MultiplicitiesFlag(M);
    fi;
end;

GModOps.EndAlgBasis := function (M)
    #
    if GModOps.EndAlgBasisFlag(M) = "unknown" then
        GModOps.ComputeEndAlgBasis(M);
    fi;
    return GModOps.EndAlgBasisFlag(M);
end;

GModOps.Decomposition := function (M)
    #
    local comps, mults;
    comps := GModOps.Indecomposables(M);
    mults := GModOps.Multiplicities(M);
    return Cat(List([1..Length(comps)],
                   i -> List([1..mults[i]], j -> comps[i])));
end;

GModOps.EndAlgRadical := function (M)
    #
    if GModOps.EndAlgRadicalFlag(M) = "unknown" then
        if GModOps.IsIndecomposable(M) <> true then
            Error("GModOps.EndAlgRadical: Module is not indecomposable!\n");
        fi;
    fi;
    return GModOps.EndAlgRadicalFlag(M);
end;

GModOps.RootResidue := function (M)
    #
    if not IsBound(M.rootResidue) then
        if GModOps.IsIndecomposable(M) <> true then
            Error("GModOps.RootResidue: Module is not indecomposable!\n");
        fi;
    fi;
    return M.rootResidue;
end;

GModOps.RootResidueOrder := function (M)
    #
    if not IsBound(M.rootResidueOrder) then
        if GModOps.IsIndecomposable(M) <> true then
            Error("GModOps.RootResidueOrder: Module is not indecomposable!\n");
        fi;
    fi;
    return M.rootResidueOrder;
end;

GModOps.Echelonisation := function (M)
    #
    if not IsBound(M.echelonIndices) then
        GModOps.EchEndAlg(M);
    fi;
    return M.echelonIndices;
end;

GModOps.IsIndecomposable := function ( module )
    if GModOps.IndecomposableFlag( module ) ="unknown" then
        GModOps.ComputeDecomp(module);
    fi;
    return GModOps.IndecomposableFlag( module );
end;


GModOps.ComputeEndAlgBasis := function (M)
    #
    # Take a module and find its endomorphism algebra.
    #
    GModOps.SetEndAlgBasisFlag(M, GModOps.ModuleHomBasis(M,M));
end;



