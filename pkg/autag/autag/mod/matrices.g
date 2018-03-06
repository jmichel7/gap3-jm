###########################################################################
##
#A  matrices.g               autag package                 Michael J Smith
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


# Miscellaneous functions involving matrices


GModOps.EcheloniseMats := function (base, F)
    #
    # compute a semi-echelonised basis for a matrix algebra
    #
    # If a linearly dependent set of elements is supplied, this
    # routine will trim it down to a basis.
    #
    local n, m, flags, ech, k, i, j, found, l;
    
    if Length(base) = 0 then
        return [ [], [] ];
    fi;

    n := Length(base[1]);
    m := Length(base[1][1]);
    flags := NullMat(n,m);
    
    ech := [];
    k := 1;

    while k <= Length(base) do
        i := 1; j := 1;
        found := false;
        while not found and i <= n do
            if (base[k][i][j] <> F.zero) and
               (flags[i][j] = 0) then
                found := true;
            else
                j := j + 1;
                if (j > m) then
                    j := 1; i := i + 1;
                fi;
            fi;
        od;

        if found then

            # Now basis element k will have echelonisation index [i,j]
            Add(ech, [i,j]);

            # First normalise the [i,j] position to 1
            base[k] := base[k] / base[k][i][j];

            # Now zero position [i,j] in all other basis elements
            for l in [k+1..Length(base)] do
                if (l <> k) and (base[l][i][j] <> F.zero) then
                    base[l] := base[l] - base[k] * base[l][i][j];
                fi;
            od;
            k := k + 1;

        else
            # no non-zero element found, delete from list
            base{[k..Length(base)-1]} := base{[k+1..Length(base)]};
            Unbind(base[Length(base)]);
            # WAS: base := base{ Cat([1..k-1], [k+1..Length(base)])};
        fi;
    od;

    return [base, ech];
    
end;


GModOps.InsertMat := function (mat1, mat2, i, j)
    #
    # This function inserts mat2 at the i,j position of mat1

    local r, c;
    r := Length(mat2);
    if r = 0 then
        return;
    else
        c := Length(mat2[1]);
        mat1{[i..i+r-1]}{[j..j+c-1]} := mat2;
    fi;
end;


GModOps.ExtractMat := function (mat, i, j, n)
    #
    # This function returns the n-by-n submatrix at position (i,j) from the
    # input matrix.

    return mat{[i..i+n-1]}{[j..j+n-1]};
end;



GModOps.UnwrapMat := function (mat)
    #
    # This function unwraps a square matrix to a vector by appending the
    # rows together.

    local vec, row;
    vec := [];
    for row in mat do
        Append(vec, row);
    od;
    #IsVector(vec);
    return vec;
end;



GModOps.WrapMat := function (arg)
    #
    # this function wraps a n^2 vector into a nxn matrix.

    local n, m, mat, i;
    if Length(arg) = 1 then
        n := RootInt(Length(arg[1]));
        m := n;
    else
        n := arg[2];
        m := arg[3];
    fi;
    if Length(arg[1]) <> m*n then
        Error("Vectors cannot be a ",m,"x",n," matrix");
    fi;
    mat := [];
    for i in [1..n] do
        Add(mat, arg[1]{[(i-1)*m+1..i*m]} );
    od;
    #List(mat, IsVector);
    return mat;
end;



GModOps.IsFittingMat := function (mat)
    #
    # check to see if a matrix is a Fitting element.  Square the matrix a
    # number of times until the rank remains the same and stays > 0.
    
    local limit, r2, r1;

    if not Int(DeterminantMat(mat)) = 0 then
        return false;
    fi;
    limit := 10;
    r2 := RankMat(mat);
    repeat
        r1 := r2;
        mat := mat^2;
        r2 := RankMat(mat);
        if r1 = r2 and r2 > 0 then
            return mat;
        fi;
        limit := limit - 1;
    until limit = 0 or r2 = 0;
    return false;
end;


