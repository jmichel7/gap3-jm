########################################################################
##
#F  CodeWeightEnumerator( <C> )
##
##  Returns a polynomial with the weight distribution
##  as coefficients.

CodeWeightEnumerator := function ( C )
    local x;
    x:=Indeterminate(Rationals);
    x.name:="x";
    if not IsCode(C) then
        Error( "<C> must be a code." );
    fi;
    return Polynomial(Rationals,WeightDistribution(C));
end;

########################################################################
##
#F  CodeDistanceEnumerator( <C>, <w> )
##
##  Returns a polynomial. The coefficient of x^i
##  is the number of codewords with distance i to w.

CodeDistanceEnumerator := function ( C, w )
    local x;
    x:=Indeterminate(Integers);
    x.name:="x";
    if not IsCode(C) then
        Error( "<C> must be a code.");
    fi;
    return Polynomial(Integers,DistancesDistribution(C, w));
end;

########################################################################
##
#F  CodeMacWilliamsTransform( <C> )
##
##  Returns a polynomial with as coefficients the weight
##  distribution of the dual code.

CodeMacWilliamsTransform := function ( C )
    local wd, tra, j, tmp, i, sc, wl;
    wd := WeightDistribution(C);
    wl := WordLength(C);
    sc := Size(C);
    tra := List([1..wl+1], x->0);
    for j in [0..wl] do
        tmp := 0;
        for i in [0..wl] do
            tmp := tmp + wd[i+1] * Krawtchouk(j, i, wl, 2);
        od;
        tra[j+1] := tmp / sc;
    od;
    return Polynomial(Rationals, tra);
end;

########################################################################
##
#F  WeightVector( <v> )
##
##  Returns the number of non-zeroes in a vector.

WeightVector := function ( v )
    local i, tel;
    tel := 0;
    for i in [1..Length(v)] do
        if v[i] <> 0 then
            tel := tel + 1;
        fi;
    od;
    return tel;
end;

########################################################################
##
#F  IsSelfComplementaryCode( <C> )
##
##  Return true if C is a complementary code, false otherwise.
##  A code is called complementary if for every v \in C
##  also 1 - v \in C (where 1 is the allone-word).

IsSelfComplementaryCode := function ( C )
    local els, t, i, alloneword;
    els := Elements(C);
    t := true;
    alloneword := List( [ 1 .. WordLength( C ) ], x -> 1 );
    for i in [ 1 .. Length( els ) ] do
        if not( ( alloneword - els[ i ] ) in els ) then
            t := false;
        fi;
    od;
    return t;
end;


########################################################################
##
#F  AllOneVector( <n> )
##
##  Return a vector with all ones.

AllOneVector := function ( n )
    return List([1..n], x->1);
end;

########################################################################
##
#F  AllOneCodeword( <n>, <F> )
##
##  Return a codeword with all ones.

AllOneCodeword := function ( n, F )
    return Codeword( AllOneVector( n ), F );
end;

#############################################################################
##
#F  IntCeiling( <r> )
##
##  Return the smallest integer greater than or equal to r.
##  3/2 => 2,  -3/2 => -1.

IntCeiling := function ( r )
    
    if IsInt(r) then
        # don't round integers
        return r;
    else
        if r > 0 then
            # round positive numbers to smallest integer 
            # greater than r (3/2 => 2)
            return Int(r)+1;
        else
            # round negative numbers to smallest integer
            # greater than r (-3/2 => -1)
            return Int(r);
        fi;
    fi;
end;

########################################################################
##
#F  IntFloor( <r> ) 
##
##  Return the greatest integer smaller than or equal to r.
##  3/2 => 1, -3/2 => -2.

IntFloor := function ( r )
    if IsInt( r ) then
        # don't round integers
        return r;
    else
        if r > 0 then
            # round positive numbers to largest integer
            # smaller than r (3/2 => 1)
            return Int(r);
        else
            # round negative numbers to largest integer
            # smaller than r (-3/2 => -2)
            return Int(r-1);
        fi;
    fi;
end;

########################################################################
##
#F  KroneckerDelta( <i>, <j> )
##
##  Return 1 if i = j,
##         0 otherwise

KroneckerDelta := function ( i, j )
    
    if i = j then
        return 1;
    else
        return 0;
    fi;
    
end;

#############################################################################
##
#F  SemiStandardForm( <mat>, <s> )
##
##  Put first s coordinates of mat in standard form. 
##  Return e as the rank of the s x s left upper
##  matrix. The coordinates s+1, ..., n are not permuted.

SemiStandardForm := function ( mat, s )
    local k, n, zero,
          stop, found,
          g, h, i, j,
          row, e, tau;
    k := Length(mat);     # number of rows: dimension
    n := Length(mat[1]);  # number of columns: wordlength
    zero := GF(2).zero;
    stop := false;
    e := 0;
    tau := ( );
    for j in [ 1..s ] do
        if not stop then
            if mat[j][j] = zero then
                # start looking for another pivot
                i := j;
                found := false;
                while ( i <= s ) and not found do
                    h := j;
                    while ( h <= k ) and not found do
                        if mat[h][i] <> zero then
                            found := true;
                        else
                            h := h + 1;
                        fi;  # if mat[h][j] <> zero
                    od;  # while ( h <= k ) and not found
                    if not found then
                        i := i + 1;
                    fi;  # if not found
                od;  # while ( i <= s ) and not found
                if not found then
                    stop := true;
                else
                    # pivot found at position (h,i)
                    # increase subrank 
                    e := e + 1;
                    # permutate the matrix so that (h,i) <-> (j,j)
                    if h <> j then
                        row := mat[h];
                        mat[h] := mat[j];
                        mat[j] := row;
                    fi;  # if h <> j
                    if i <> j then
                        tau := tau * (i,j);
                        for g in [ 1 .. k ] do
                            mat[g] := Permuted( mat[g], (i,j) );
                        od;  # for g in [ 1..k ]
                    fi;  # if i <> j
                fi;  # if not found
            else
                e := e + 1;
            fi;  # if mat[j][j] = zero
           
            if not stop then
                for i in [ 1..k ] do
                    if i <> j then
                        if mat[i][j] <> zero then
                            mat[i] := mat[i] + mat[j];
                        fi;  # if mat[i][j] <> zero
                    fi;  # if i <> j
                od;  # for i in [ 1..k ]
            fi;  # if not stop
        fi;  # if not stop
    od;  # for j in [ 1..s ] do
    return [ e, tau ];
end;


