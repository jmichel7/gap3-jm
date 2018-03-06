########################################################################
##
#F  CoordinateSubCode( <code>, <i>, <element> )
##
##  Return the subcode of <code>, that has elements
##  with an <element> in coordinate position <i>.
##  If no elements have an <element> in position <i>, return false.
##

CoordinateSubCode := function ( code, i, element )

    local els;

    if not IsCode( code ) then
        Error( "CoordinateSubCode: <code> must be a code" );
    fi;
    if not IsInt( i ) or i < 1 or i > WordLength( code ) then
        Error( "CoordinateSubCode: <i> must lie in the range [ 1 .. n ]" );
    fi;
    if not ( element in Elements( Field( code ) ) ) then
        Error( "CoordinateSubCode: <element> must be an element of ",
                Field( code ) );
    fi;

    els := Elements( code );
    els := VectorCodeword( els );
    els := Filtered( els, x -> x[ i ] = element );
    if Length( els ) = 0 then
        return false;
    else
        return ElementsCode( els, "subcoordinate code", Field( code )  );
    fi;
end;


########################################################################
##
#F  CoordinateNorm( <code>, <i> )
##  
##  Returns the norm of code with respect to coordinate i.
##

CoordinateNorm := function ( code, i )

    local max, els, subcode, f, j, w, n;

    if not IsCode( code ) then
        Error( "CoordinateNorm: <code> must be code" );
    fi;
    if i < 1 or i > WordLength( code ) then
        Error( "CoordinateNorm: <i> must lie in the range [ 1 .. n ]" );
    fi;

    if not IsBound( code.coordinateNorm ) then
        code.coordinateNorm := List( [ 1 .. WordLength( code ) ], x -> -1 );
    fi;

    if code.coordinateNorm[ i ] = -1 then
        max := -1;
        els := Elements( Field( code ) );
        subcode := [ 1 .. Length( els ) ];
        f := [ 1 .. Length( els ) ];
        for j in [ 1 .. Length( els ) ] do
            subcode[ j ] := CoordinateSubCode( code, j, els[ j ] );
        od;
        for w in Codeword( CosetLeadersMatFFE( CheckMat( code ),
                           Size( Field( code ) ) ) ) do
            for j in [ 1 .. Length( els ) ] do
                if subcode[ j ] = false then
                    f[ j ] := WordLength( code );
                else
                    f[ j ] := MinimumDistance( subcode[ j ], w );
                fi;
            od;
            n := Sum( f );
            if n > max then
                max := n;
            fi;
        od;
        code.coordinateNorm[ i ] := max;
    fi;
    return code.coordinateNorm[i];
end;


########################################################################
##
#F  CodeNorm( <code> )
##
##  Return the norm of code.
##  The norm of code is the minimum of the coordinate norms
##  of code with respect to i = 1, ..., n.
##

CodeNorm := function ( code )
    if not IsBound( code.codeNorm ) then
        code.codeNorm := Minimum( List( [ 1 .. WordLength( code ) ],
                                    x -> CoordinateNorm( code, x ) ) );
    fi;
    return code.codeNorm;
end;


########################################################################
##
#F  IsCoordinateAcceptable( <code>, <i> )
##
##  Test whether coordinate i of <code> is acceptable.
##  (a coordinate is acceptable if the norm of code with respect to
##   that coordinate is less than or equal to one plus two times the 
##   covering radius of code).

IsCoordinateAcceptable := function ( code, i )
    
    local cr;
    
    cr := CoveringRadius( code );
    if IsInt( cr ) then
        if CoordinateNorm( code, i ) <= 2 * cr + 1 then
            return true;
        else
            return false;
        fi;
    else
        Error( "IsCoordinateAcceptable: Sorry, the covering radius is ",
               "not known and not easy to compute." );
    fi;
    
end;


########################################################################
##
#F  IsNormalCode( <code> )
##
##  Return true if code is a normal code, false otherwise.
##  A code is called normal if its norm is smaller than or
##  equal to two times its covering radius + one.
##

IsNormalCode := function ( code )
    if not IsCode( code ) or Field( code ) <> GF(2) then
        Error( "IsNormalCode: <code> must be a binary code" );
    fi;
    if not IsBound( code.isNormalCode ) then
        code.isNormalCode := code.operations.IsNormalCode( code );
    fi;
    return code.isNormalCode;
end;

CodeOps.IsNormalCode := function ( code )

    local n, k, d, r, i, isnormal;

    if IsLinearCode( code ) then
        return IsNormalCode( code );
    else
        n := WordLength( code );
        k := Dimension( code );
        d := MinimumDistance( code );
        r := CoveringRadius( code );
        if not IsInt( r ) then
            r := -1;
        fi;
        if d = 2 * r
           or ( d = 2 * r - 1 and EuclideanRemainder( n, r ) <> 0 )
           or ( r = 1 and n <= 9 )
           or ( r = 1 and Size( code ) <= 95 )
           then
            return true;
        else
            if r >= 0 then
                i := 1;
                isnormal := false;
                while i <= n and not isnormal do
                    isnormal := IsCoordinateAcceptable( code, i );
                    i := i + 1;
                od;
                return isnormal;
            else
                Error( "IsNormalCode: sorry, the covering radius for ",
                       "this code has not yet been computed." );
            fi;
        fi;
    fi;
end;

LinCodeOps.IsNormalCode := function ( code )

    local n, k, d, r, i, isnormal;

    n := WordLength( code );
    k := Dimension( code );
    d := MinimumDistance( code );
    r := CoveringRadius( code );
    if not IsInt( r ) then
        r := -1;
    fi;
    if d = 2 * r
       or ( d = 2 * r - 1 and EuclideanRemainder( n, r ) <> 0 )
       or ( r = 1 and n <= 9 )
       or ( r = 1 and Size( code ) <= 95 )
       # the following conditions are only valid for linear codes
       or ( n <= 15 )
       or ( k <= 5 )
       or ( n-k <= 7 )
       or ( d <= 4 )
       or ( r >= 0 and r <= 3 )
       or ( IsPerfectCode( code ) )
       then
        return true;
    else
        if r >= 0 then
            # a code is normal if one of the coordinates is acceptable
            i := 1;
            isnormal := false;
            while i <= n and not isnormal do
                isnormal := IsCoordinateAcceptable( code, i );
                i := i + 1;
            od;
            return isnormal;
        else
            Error( "IsNormalCode: sorry, the covering radius for ",
                   "this code has not yet been computed." );
        fi;
    fi;
end;

CycCodeOps.IsNormalCode := LinCodeOps.IsNormalCode;


########################################################################
##
#F  GeneralizedCodeNorm( <code>, <code1>, <code2>, ... , <codek> )
## 
##  Compute the k-norm of code with respect to the k subcode
##  code1, code2, ... , codek.
##

GeneralizedCodeNorm := function ( arg )
    local i, mindist, min, max, globalmax, x, union,word;
    
    if Length( arg ) < 2 then
        Error( "GeneralizedCodeNorm: no subcodes are specified" );
    fi;
    if not IsCode( arg[ 1 ] ) then
        Error( "GeneralizedCodeNorm: <code> must be a code" );
    fi;
    union := arg[ 2 ];
    for i in [ 1 .. Length( arg ) - 1 ] do
        if WordLength( arg[ i + 1 ] ) <> WordLength( arg[ 1 ] ) then
            Error( "GeneralizedCodeNorm: length of code ", i,
                   " is not equal to the length of <code>" );
        fi;
        if not ( arg[ i + 1 ] in arg[ 1 ] ) then
            Error( "GeneralizedCodeNorm: code ", i,
                   " is not a subcode of code." );
        fi;
        if i > 1 then
            union := AddedElementsCode( union,
                       Elements( arg[ i + 1 ] ) );
        fi;
    od;
    if arg[ 1 ] <> union then
        Error( "GeneralizedCodeNorm: <code> is not the union of the ",
               "subcodes" );
    fi;

    globalmax := -1;
    for word in Elements( arg[ 1 ] ) do
        mindist := List( [ 2 .. Length( arg ) ],
                          x -> MinimumDistance( arg[ x ], word ) );
        min := Minimum( mindist );
        max := Maximum( mindist );
        if min + max > globalmax then
            globalmax := min + max;
        fi;
    od;
    return globalmax;
end;

