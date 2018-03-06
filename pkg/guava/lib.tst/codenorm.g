########################################################################
##
#F  SubCoordinateCode( <C>, <i>, <e> )
##
##  Return the subcode of C, that has elements
##  with an e in coordinate position i.
##  If no elements has an e in position i, return false.

SubCoordinateCode := function ( C, i, e )
    local els;
    els := Elements( C );
    els := VectorCodeword( els );
    els := Filtered( els, x -> x[ i ] = e );
    if Length(els) = 0 then
        return false;
    else
        return ElementsCode( els, Field( C )  );
    fi;
end;

########################################################################
##
#F  SubCode( <C>, <D> ) 
##
##  Return true if C is a subcode of D, false otherwise.
##  C is a subcode of D if all the elements of C are
##  also elements of D.

SubCode := function ( C, D )
    local cels, dels, x, subset;
    cels := Elements( C );
    dels := Elements( D );
    subset := true;
    for x in cels do
        if not (x in dels) then
            subset := false;
        fi;
    od;
    return subset;
end;

########################################################################
##
#F  CoordNorm( <C>, <i> ) 
##  
##  Returns the norm of C with respect to coordinate i.

CoordNorm := function ( C, i )
    local cnorm, ci;
    if not IsBound( C.coordNorm ) then
        cnorm := List( [1..WordLength( C )], x -> -1 );
        ci := C.operations.CoordNorm( C, i );
        cnorm[i] := ci;
        C.coordNorm := cnorm;
    fi;
    if C.coordNorm[i] = -1 then
        C.coordNorm[i] := C.operations.CoordNorm( C, i );
        if C.coordNorm[i] = -1 then
            Error( "<C> must be linear" );
        fi;
    fi;
    return C.coordNorm[i];
end;

CodeOps.CoordNorm := function ( C, i )
    if IsLinearCode( C ) then
        return CoordNorm( C, i );
    else
        Error( "the code must be linear" );
    fi;
end;

LinCodeOps.CoordNorm := function ( C, i )
    local f0, f1, max, C0, C1, w, f;
    max := -1;
    C0 := SubCoordinateCode( C, i, Field(C).zero );
    C1 := SubCoordinateCode( C, i, Field(C).one );
    if (C0 = false) or (C1 = false) then
        Error( "Code <C> cannot be split into 2 subcodes on coord. ", i);
    fi;
    for w in Codeword(CosetLeadersMatFFE(CheckMat(C),Size(Field(C)))) do
        f0 := MinimumDistance(C0, w);
        f1 := MinimumDistance(C1, w);
        f := f0 + f1;
        if f > max then
            max:=f;
        fi;
    od;
    return max;
end;

CycCodeOps.CoordNorm := LinCodeOps.CoordNorm;

########################################################################
##
#F  CodeNorm( <C> ) 
##
##  Return the norm of C.
##  The norm of C is the minimum of the coordinate norms
##  of C with respect to i = 1, ..., n.

CodeNorm := function ( C )
    if not IsBound( C.codeNorm ) then
        C.codeNorm := C.operations.CodeNorm( C );
    fi;
    return C.codeNorm;
end;

CodeOps.CodeNorm := function ( C )
    if IsLinearCode( C ) then
        return CodeNorm( C );
    else
        Error( "The norm of a code only exists for linear codes." );
    fi;
end;

LinCodeOps.CodeNorm := function ( C )
    local min, i, ni;
    min := 10^100;
    for i in [1..WordLength(C)] do
        ni := CoordNorm( C, i );
        if ni < min then
            min := ni;
        fi;
    od;
    return min;
end;

CycCodeOps.CodeNorm := LinCodeOps.CodeNorm;

########################################################################
##
#F  CoordinateAcceptable( <C>, <i> )
##
##  Test whether coordinate i of code C is acceptable.
##  (a coordinate is acceptable if the norm of C with respect to 
##   that coordinate is less than or equal to one plus two times the 
##   covering radius of C).

CoordinateAcceptable := function ( C, i )
    
    if CoordNorm( C, i ) <= 2 * CoveringRadius( C ) + 1 then
        return true;
    else
        return false;
    fi;
    
end;

########################################################################
##
#F  GeneralizedCodeNorm( <C>, <C1>, <C2>, ... , <Ck> ) 
## 
##  Compute the k-norm of C with respect to the k subcode
##  C1, C2, ... , Ck.

GeneralizedCodeNorm := function ( arg )
    local k, i, cels, max, x, mi, j, mj, ma, min;
    
    if Length(arg) < 2 then
        Error( "At least one subcode must be specified." );
    fi;
    
    k := Length(arg) - 1;
    cels := Elements(arg[1]);
    for i in [1 .. k] do
        if not SubCode( arg[i+1], arg[1] ) then
            Error( "C",i," must be a subcode of C." );
        fi;
    od;
    max := -1;
    for x in cels do
        min := 10^100;
        for i in [1..k] do
            mi := MinimumDistance( arg[i+1], x );
            if mi < min then
                min := mi;
            fi;
        od;
        ma := -1;
        for j in [1..k] do
            mj := MinimumDistance( arg[j+1], x );
            if mj > ma then
                ma := mj;
            fi;
        od;
        if min + ma > max then
            max := min + ma;
        fi;
    od;
    return max;
end;

########################################################################
##
#F  IsNormalCode( <C> ) 
##
##  Returns true if C is a normal code, false otherwise.
##  A code is calles normal if its norm is smaller than or
##  equal to two times its covering radius + one.

IsNormalCode := function ( C )
    local a;
    if not IsBound( C.isNormalCode ) then
        C.isNormalCode := C.operations.IsNormalCode( C );
    fi;
    return C.isNormalCode;
end;

CodeOps.IsNormalCode := function ( C )
    if IsLinearCode( C ) then
        return IsNormalCode(C);
    else
        Error( "Only linear codes can be tested for normality." );
    fi;
end;

LinCodeOps.IsNormalCode := function ( C )
    local n, k, d, r;
    n := WordLength( C );
    k := Dimension( C );
    d := MinimumDistance( C );
    r := CoveringRadius( C );
    if Field( C ) = GF( 2 ) then
        if n <= 15 or
           k <= 5 or
           d <= 4 or
           r <= 3 or
           n-k <= 7 or
           IsPerfectCode( C ) or
           d >= 2 * r or
           (r = 1 and n <= 8) then
            return true;
        fi;
    fi;
    if CodeNorm(C) <= 2*CoveringRadius(C) + 1 then
        return true;
    else
        return false;
    fi;
end;

CycCodeOps.IsNormalCode := LinCodeOps.IsNormalCode;

