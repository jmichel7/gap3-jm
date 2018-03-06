###########################################################################
##
#A  eqns.g                   autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################


# These routines are designed to accumulate a system of linear equations
#
#    M_1 X = V_1,  M_2 X = V_2 ...  M_t X = V_t
#
# Where each M_i is an m_i*n matrix, X is the unknown length n vector, and
# each V is an length m_i vector.  The equations can be added as each batch
# is calculated. Here is some pseudo-code to demonstrate:
#
#   eqns := newEqns (n, field);
#   i := 1;
#   repeat
#     <calculate M_i and V_i>
#     addEqns(M_i, V_i)
#     increment i;
#   until  i > t  or  eqns.fail;
#   if not eqns.fail then
#     S := solveEqns(eqns);
#   fi;
#
# As demonstrated by the example, an early notification of failure is
# available by checking ".fail".  All new equations are sifted with respect
# to the current set, and only added if they are independent of the current
# set. If a new equation reduces to the zero row and a nonzero vector
# entry, then there is no solution and this is immediately returned by
# setting eqns.fail to true.  The function solveEqns has an already
# triangulised system of equations, so it simply reduces above the pivots
# and returns the solution vector.


EqnOps := rec();


InfoEqns := Ignore;
InfoEqns2 := Ignore;


EqnOps.NewEqns := function (arg)

    local X, n, F, V, eqns;

    InfoEqns2( "#T newEqns: entering\n" );

    if Length(arg) <2 then
        Error("NewEqns(dim, field) or NewEqns(X, V)");
    fi;

    if IsInt(arg[1]) then
        X := false;
        n := arg[1];
        F := arg[2];
    else
        X := arg[1];
        V := arg[2];
        n := Length(X[1]);
        F := Field(X[1][1]); # Note: prime field only
    fi;

    eqns := rec();
    eqns.dim := n;              # number of variables
    eqns.field := F;            # field over which the equation hold
    eqns.mat := [];             # left-hand sides of system
    eqns.weights := [];         # echelon weights for lhs matrix
    eqns.vec := [];             # right-hand sides of system
    eqns.fail := false;         # flag to indicate inconsistent system
    eqns.index := [];           # index for row ordering
    eqns.operations := rec();
    eqns.operations.Print := function (r)
        Print("Fail: ", r.fail, "\n");
        Print("Dimension: ",r.dim, "\n");
        Print("Field: ",F, "\n");
        Print("Matrix:\n");
        Nice(r.mat);
        Print("Vector:\n");
        Nice(r.vec);
        Print("Weights:", r.weights, "\n");
        Print("Index:", r.index, "\n");
    end;
    InfoEqns2( "#T newEqns: leaving\n" );

    if IsMat(X) then
        EqnOps.AddEqns(eqns, X, V);
    fi;

    eqns.AddTime := 0; # accumulate time used to add equations

    return eqns;
end;
    



EqnOps.AddEqns := function ( eqns, newmat, newvec)
    #
    # Add a bunch of equations to the system of equations in <eqns>.  Each
    # row of <newmat> is the left-hand side of a new equation, and the
    # corresponding row of <newvec> the right-hand side. Each equation in
    # filtered against the current echelonised system stored in <eqns> and
    # then added if it is independent of the system.  As soon as a
    # left-hand side reduces to 0 with a non-zero right-hand side, the flag
    # <eqns.fail> is set.

    local t0, n, zero, weights, mat, vec, NextPositionProperty, ReduceRow, 
          k, t, newweight, newrow, newrhs, i, l;

    InfoEqns2( "#T addEqns: entering\n" );
    t0 := Runtime();

    n := eqns.dim;
    zero := eqns.field.zero;
    weights := eqns.weights;
    mat := eqns.mat;
    vec := eqns.vec;

    NextPositionProperty := function (list, func, start )
        local   i;
        for i  in [ start .. Length( list ) ]  do
            if func( list[ i ] )  then
                return i;
            fi;
        od;
        return false;
    end;

    ReduceRow := function (lhs, rhs)
        #
        # reduce the (lhs,rhs) against the semi-echelonised current matrix,
        # and return either: (1) the reduced rhs if the lhs reduces to zero,
        # or (2) a list containing the new echelon weight, the new row and
        # the new rhs for the system, and the row number that this
        # equation should placed.

        local lead, i, z;
        lead := PositionProperty(lhs, IsNonZeroElt);
        if lead = false then
            return rhs;
        fi;
        for i in [1..Length(weights)] do
            if weights[i] = lead then
                z := lhs[lead];
                lhs := lhs - z * mat[i]; rhs := rhs - z * vec[i];
                lead := NextPositionProperty(lhs, IsNonZeroElt, lead);
                if lead = false then
                    return rhs;
                fi;
            elif weights[i] > lead then
                return [lead, lhs, rhs, i];
            fi;
        od;
        return [lead, lhs, rhs, Length(weights)+1];
    end;
    
      
    for k in [1..Length(newmat)] do
        
        t := ReduceRow(newmat[k], newvec[k]);
        
        if IsList(t) then
            # new equation
            newweight := t[1];
            newrow := t[2];
            newrhs := t[3];
            i := t[4]; # position for new row

            # normalise so that leading entry is 1
            newrhs := newrhs / newrow[newweight]; 
            newrow := newrow / newrow[newweight]; # NB: in this order

            if i = Length(mat)+1 then
                # add new equation to end of list
                Add(mat, newrow);
                Add(vec, newrhs);
                Add(weights, newweight);
            else
                l := Length(mat);
                # move down other rows to make space for this new one...
                mat{[i+1..l+1]} := mat{[i..l]};  
                vec{[i+1..l+1]} := vec{[i..l]};
                # and then slot it in
                mat[i] := newrow;
                vec[i] := newrhs;
                weights{[i+1..l+1]} := weights{[i..l]};
                weights[i] := newweight;
            fi;

        else
            # no new equation, check whether inconsistent due to
            # nonzero rhs reduction
            
            if IsNonZeroElt(t) then
                InfoEqns2( "addEqns: FAIL!\n" );
                eqns.fail := true;
                return eqns; # return immediately
            fi;
        fi;
    od;
    eqns.AddTime := eqns.AddTime + Runtime() - t0;
    InfoEqns2( "#T addEqns: leaving\n" );
end;




EqnOps.KillAbovePivotsEqns := function (eqns)
    #
    # Eliminate entries above pivots. Note that the pivot entries are
    # all 1 courtesy of AddEqns.

    local m, n, zero, i, c, j, factor;

    InfoEqns2( "#T killAbovePivotsEqns: entering\n" );
    m := Length(eqns.mat);
    n := eqns.dim;
    if m > 0 then
        zero := eqns.field.zero;
        for i in [1..m] do
            c := eqns.weights[i];
            for j in [1..i-1] do
                if eqns.mat[j][c] <> zero then
                    InfoEqns2( "solveEqns: kill mat[",j,",",c,"]\n");
                    factor := eqns.mat[j][c];
                    eqns.mat[j] := eqns.mat[j] - factor*eqns.mat[i];
                    eqns.vec[j] := eqns.vec[j] - factor*eqns.vec[i];
                fi;
            od;
        od;
    fi;
    InfoEqns2( "#T killAbovePivotsEqns: leaving\n" );
end;



EqnOps.SolveEqns := function (eqns)
    #
    # Solve the equations.
    
    local t, m, n, zero, ans, i;
    
    InfoEqns2( "#T solveEqns: entering\n" );
    t := Runtime();

    m := Length(eqns.mat);
    n := eqns.dim;
    zero := eqns.field.zero;
    
    if eqns.fail then
        return false;
    fi;

    EqnOps.KillAbovePivotsEqns(eqns);
    
    ans := [1..n]*zero;
    for i in [1..m] do
        ans[eqns.weights[i]] := eqns.vec[i];
    od;
    
    InfoEqns2( "#T solveEqns: leaving\n" );
    InfoEqns( "#T solveEqns: ",Time(Runtime()-t)," (with ",
            Time(eqns.AddTime), " in AddEqns)\n" );
    return ans;
end;


EqnOps.NullspaceEqns := function (e)
    #
    # Take the matrix stored in equation record <e> and compute a basis
    # for its nullspace, ie  x  such that  mat * x = 0.  Note that the
    # vector is on the other side of the matrix from GAP's NullspaceMat.
    # This means we get to skip the Transposing that occurs at the top
    # of that function (a bonus!).
    #
    # This function is a modified version NullspaceMat in matrix.g

    local mat, m, n, zero, one, empty, i, k, nullspace, row;

    EqnOps.KillAbovePivotsEqns(e);
    mat := e.mat;

    m := Length(mat);
    n := e.dim;
    zero := e.field.zero;
    one  := e.field.one;

    # insert empty rows to bring the leading term of each row on the diagonal
    empty := zero*[1..n];
    i := 1;
    while i <= Length(mat)  do
        if i < n  and mat[i][i] = zero  then
            for k in Reversed([i..Minimum(Length(mat),n-1)])  do
                mat[k+1] := mat[k];
            od;
            mat[i] := empty;
        fi;
        i := i+1;
    od;
    for i  in [ Length(mat)+1 .. n ]  do
        mat[i] := empty;
    od;

    # The following comment from NullspaceMat:
    # 'mat' now  looks  like  [ [1,2,0,2], [0,0,0,0], [0,0,1,3], [0,0,0,0] ],
    # and the solutions can be read in those columns with a 0 on the diagonal
    # by replacing this 0 by a -1, in  this  example  [2,-1,0,0], [2,0,3,-1].
    nullspace := [];
    for k in [1..n] do
        if mat[k][k] = zero  then
            row := [];
            for i  in [1..k-1]  do row[i] := -mat[i][k];  od;
            row[k] := one;
            for i  in [k+1..n]  do row[i] := zero;  od;
            Add( nullspace, row );
        fi;
    od;

    return nullspace;
end;

