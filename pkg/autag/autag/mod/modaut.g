###########################################################################
##
#A  modaut.g                 autag package                 Michael J Smith
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


# Compute generators for the automorphism group of a module, or, equivalently,
# a generating set for the centraliser in the general linear group of the
# matrix group acting on the module.


GModOps.ModuleAutGens := function (M)
    if not IsBound(M.auts) then
        GModOps.ComputeModuleAuts(M);
    fi;
    return M.auts;
end;


GModOps.ModuleAutSize := function (M)
    if not IsBound(M.autorder) then
        GModOps.ComputeModuleAuts(M);
    fi;
    return M.autorder;
end;


CentraliserMatGroup := function (G)
    # 
    # Use module routines to compute centraliser of matrix group
    # in GL
    
    local M, C;
    
    M := Gmodule(G.generators);
    GModOps.ComputeModuleAuts(M);

    C := Group(GModOps.ModuleAutGens(M), G.identity);
    C.size := GModOps.ModuleAutSize(M);

    return C;
end;
CentralizerMatGroup := CentraliserMatGroup;




GModOps.ComputeModuleAuts := function (M)
    #
    # Calculate a generating set for the module automorphism group.
    # These are non-singular matrices that centralise the matrix group
    # defining M.
    
    local GL, SizeGL, circleGenerators, auts, autorder, n, F, comps, mults, 
          i, w, q, r, Fq, gl, g, X, j, k, nilbase, d, homs;

    GL := function( n, q )
        local	F, i, mat1, mat2;
        F := GF( q );
        if n = 1 then return Group([[F.root]]); fi;
        if q = 2 then return SL( n, 2 ); fi;
        if n = 2 then
            mat1 := [ [ F.root, 0*F.one ], [ 0*F.one,   F.one ] ];
            mat2 := [ [ -F.one,   F.one ], [  -F.one, 0*F.one ] ];
        else
            mat1 := IdentityMat( n, F ); mat1[ 1 ][ 1 ] := F.root;

            mat2 := NullMat( n, n, F );
            mat2[ 1 ][ 1 ] := -F.one; mat2[ 1 ][ n ] := F.one;
            for i in [2..n] do mat2[ i ][ i-1 ] := -F.one; od;
        fi;
        return Group( mat1, mat2 );
    end;

    SizeGL := function( n, q )
        return Product( [0..n-1], x->(q^n-q^x) );
    end;

    
    circleGenerators := function (matalg)
        local ret;
        if Length(matalg) = 0 then
            return [];
        fi;
        ret := GModOps.NilpotentBasis(matalg);
        # conjugate to convert back to original basis
        #
        return ret[2]^-1 * ret[1] * ret[2];
    end;


    auts := [];
    autorder := 1;
    
    n := GModOps.DimFlag(M);
    F := GModOps.Field(M);
    comps := GModOps.Indecomposables(M);
    mults := GModOps.Multiplicities(M);
    
    # start by building those automorphisms that fix the homogeneous
    # components - ie, do not involve maps from M_i to M_j unless
    # M_i is the same isomorphism type as M_j

    for i in [1..Length(comps)] do
        
        w := GModOps.RootResidue(comps[i]);
        q := GModOps.RootResidueOrder(comps[i]) + 1;
        r := mults[i];
        Fq := GF(q);
        gl := GL(r,q);
        autorder := autorder * SizeGL(r,q);
        
        # first the subgroup GL(multiplicity, residue field)
        #
        for g in gl.generators do
            
            X := IdentityMat(n, F);
            
            for j in [1..r] do
                for k in [1..r] do
                    if g[j][k] <> F.zero then
                        GModOps.InsertMat(X, w^LogFFE(g[j][k]),
                                GModOps.DecompIndex(M,i,j),
                                GModOps.DecompIndex(M,i,k));
                    elif j = k then
                        GModOps.InsertMat(X, w * 0,
                                GModOps.DecompIndex(M,i,j),
                                GModOps.DecompIndex(M,i,k));
                    fi;
                od;
            od;

            Add(auts, X);
        od;
        
        # now the subgroup { I + Y | Y in S } where S generates the radical
        # of the endomorphism algebra as a circle group
        #
        nilbase := circleGenerators(GModOps.EndAlgRadical(comps[i]));
        d := GModOps.DimFlag(comps[i]);
        
        autorder := autorder * (F.size ^ Length(nilbase)) ^ (r^2);

        for j in [1..Length(nilbase)] do

            X := IdentityMat(n,F);
            GModOps.InsertMat(X, IdentityMat(d,F) + nilbase[j],
                    GModOps.DecompIndex(M,i,1), GModOps.DecompIndex(M,i,1));

            Add(auts, X);
        od;
        
    od;


    # Now the automorphisms that act trivially when restricted to
    # each homogeneous component, but which include action between
    # homogeneous components via elements of Hom(M_i, M_j)
    #
    for i in [1..Length(comps)] do
        for j in [1..Length(comps)] do
            
            if i <> j then
                
                homs := GModOps.ModuleHomBasis(comps[i], comps[j]);

                if Length(homs) > 0 then
                    autorder := autorder * (F.size ^ Length(homs)) ^
                                (mults[i] * mults[j]);
                fi;
                
                for k in [1..Length(homs)] do
                    
                    X := IdentityMat(n,F);
                    GModOps.InsertMat(X, homs[k],
                            GModOps.DecompIndex(M,i,1), 
                            GModOps.DecompIndex(M,j,1));
                    Add(auts, X);

                od;
            fi;
        od;
    od;

    # Now transform basis and return automorphisms
    #
    X := GModOps.Basis(M);

    M.autorder := autorder;
    M.auts := X^-1 * auts * X;   # auts is never empty
    # changed above line from X * auts * X^-1, 17:11 Sun 25 Feb 1996
end;



GModOps.NilpotentBasis := function (matalg)
    #
    # this function is loosely based on MatGroupOps.CompositionFactors
    # 
    # compute a change of basis that exhibits the matrix algebra
    # defined by the basis 'matalg' in triangular form.
    
    local decompose, field, Y, mats, newbase;

    decompose := function ( m, b )

        local n, subs, vs, rep, newm;

        if Length(m) = 0 then
            #
            # all action is now zero, so append current full basis and
            # finish up
            
            Append(Y, b);
            
        else
            
            n := Length(m[1][1]);

            # find the intersection of the nullspaces
            #
            subs := Intersection(List(
                        m, x -> VectorSpace(
                                NullspaceMat(x), field)));
            
            # compute the quotient of the vector space by the intersection
            # of the nullspaces
            #
            vs := VectorSpace(IdentityMat(n,field),field) mod subs;

            # Use matrix group routine to compute action of nilpotent
            # matrices on the quotient vectorspace
            #
            rep := MTXOps.MatRepresentation( m, vs );

            # Take a copy of the non-zero matrices acting on the quotient space
            #
            newm := Filtered(rep.range, x -> not x = 0*x);
            
            Append(Y, Base(subs) * b);
            decompose( newm, rep.base * b );

        fi;

    end;

    # return empty list if empty matrix list
    #
    if Length(matalg) = 0 then return []; fi;

    field := Field(matalg[1][1][1]);

    Y   := [];

    decompose( matalg, IdentityMat(Length(matalg[1][1]), field));
    #
    # Y is the change of basis matrix

    if Length(matalg) > 0 then
        mats := Y * matalg * Y^-1;
    fi;
    #
    # mats is now a list of matrices in lower triangular form
    
    # echelonise them along lower diagonals
    #
    newbase := GModOps.EcheloniseNilpotentMatAlg(mats, field)[1];
    
    return [newbase, Y];
    
end;



GModOps.EcheloniseNilpotentMatAlg := function (matalg, F)
    #
    # Note: matalg is a basis for a nilpotent matrix algebra whose elements
    # are all in lower diagonal form (zeros on the main diagonal).
    #
    # Echelonisation indices are chosen as the earliest non-zero entries
    # running down diagonals below the main diagonal: 
    #   [2,1], [3,2], [4,3], ..., [3,1], [4,2], ..., [n-1,1], [n, 2], [n,1]

    local n, flags, base, ech, k, diff, i, j, found, l;
    
    n := Length(matalg[1][1]);
    flags := NullMat(n,n);
    
    base := matalg;
    ech := [];
    k := 1;

    while k <= Length(base) do
        diff := 1;
        i := 2; j := i - diff;
        found := false;
        while not found and diff < n do
            if (base[k][i][j] <> F.zero) and
               (flags[i][j] = 0) then
                found := true;
            else
                i := i + 1; j := i - diff;
                if (i > n) then
                    diff := diff + 1;
                    i := diff + 1; j := i - diff;
                fi;
            fi;
        od;

        if found then

            # Now basis element k will have echelonisation index [i,j]
            Add(ech, [i,j]);

            # First normalise the [i,j] position to 1
            base[k] := base[k] / base[k][i][j];

            # Now zero position [i,j] in all other basis elements
            for l in [1..Length(base)] do
                if (l <> k) and (base[l][i][j] <> F.zero) then
                    base[l] := base[l] - base[k] * base[l][i][j];
                fi;
            od;
            k := k + 1;

        else
            # no non-zero element found, delete from list
            base := base{ Cat([1..k-1], [k+1..Length(base)])};
        fi;
    od;
    return [base, ech];
    
end;
