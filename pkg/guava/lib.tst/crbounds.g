########################################################################
##
#F  LowerBoundCoveringRadiusVanWee( <n>, <r> [, <givesize> ] )
##
##  Van Wee Bound
##
LowerBoundCoveringRadiusVanWee := function ( arg )
    local i, num, den, tmp, n, r, sizewanted, m;
    
     n := arg[ 1 ];

    if Length( arg ) = 3 then
        sizewanted := arg[ 3 ];
    else
        sizewanted := true;
    fi;
    
    if sizewanted then
        r := arg[ 2 ];
        tmp := (Binomial(n,r))/(IntCeiling((n-r)/(r+1)));
        tmp := tmp * (IntCeiling((n+1)/(r+1)) - (n+1)/(r+1));
        num := SphereContent( n, r, 2 ) - tmp;
        den := 2^n;
        return IntCeiling(den/num);
    else
        m := arg[ 2 ];
        if m < 0 then
            Error( "LBCRVanWee: Size must be non-negative." );
        fi;
        if m = 0 or m = 1 then
            return n;
        else
            r := 0;
            tmp := 1;
            while m * tmp < 2^n do
                r := r + 1;
                tmp := (Binomial(n,r))/(IntCeiling((n-r)/(r+1)));
                tmp := tmp * (IntCeiling((n+1)/(r+1)) - (n+1)/(r+1));
                tmp := SphereContent( n, r, 2 ) - tmp;
            od;
            return r;
        fi;
    fi;
end;

########################################################################
##
#F  LowerBoundCoveringRadiusEmbedded1( <n>, <r> [, <givesize> ] )
##

LowerBoundCoveringRadiusEmbedded1 := function ( arg )
    local num, den, tmp, n, r, sizewanted, m, valid, old;
    
     n := arg[ 1 ];

    if Length( arg ) = 3 then
        sizewanted := arg[ 3 ];
    else
        sizewanted := true;
    fi;
    
    if sizewanted then
        r := arg[ 2 ];
        num := 2^n - UpperBound( n, 2 * r + 1 ) * Binomial( 2 * r, r );
        den := SphereContent( n, r, 2 ) - Binomial( 2 * r, r );
        if den <= 0 then 
            return 0;
        else
            return IntCeiling( num / den );
        fi;
    else
        m := arg[ 2 ];
        if m < 0 then
            Error( "LBCREmbedded1: Size must be non-negative." );
        fi;
        if m = 0 or m = 1 then
            return n;
        elif m = 2^n then
            return 0;
        else
            r := 0;
            tmp := 0;
            valid := false;
            while ( m < tmp or not valid ) and r < n do
                old := tmp;
                r := r + 1;
                tmp := SphereContent( n, r, 2 ) - Binomial( 2 * r, r );
                if tmp > 0 then 
                    tmp := IntCeiling( ( 2^n - UpperBound( n, r ) ) / tmp );
                    valid := ( old >= m );
                fi;
            od;
            if r = n then 
                return 0;
            else
                return r;
            fi;
        fi;
    fi;
end;

########################################################################
##
#F  LowerBoundCoveringRadiusEmbedded2( <n>, <r> [, <givesize> ] )
##
LowerBoundCoveringRadiusEmbedded2 := function ( n, r )
    local num, den;
    den := SphereContent( n, r, 2) - 3/2 * Binomial( 2*r, r );
    if den <= 0 then
        return 0;
    fi;
    num := 2^n - 2 * UpperBound( n, 2*r+1 ) * Binomial( 2*r, r );
    return IntCeiling( num / den );
end;

#############################################################################
##
#F  LowerBoundCoveringRadiusExcess1( <n>, <r> ) Counting Excess bound
##
LowerBoundCoveringRadiusExcess1 := function ( n, t )
    local eps, num, den;
    
    if n <= t then
        Error( "<n> must be greater than than <t>" );
    fi;
    
    eps := (t+1) * IntCeiling((n+1)/(t+1)) - (n+1);
    num := (n-t+eps) * 2^n;
    den := (n-t) * SphereContent(n,t,2) + eps * SphereContent(n,t-1,2);
    
    return IntCeiling(num/den);
end;

#############################################################################
##
#F  LowerBoundCoveringRadiusExcess2( <n>, <r> ) Counting Excess bound
##
LowerBoundCoveringRadiusExcess2 := function ( n, t )
    local eps, num, den;
    
    if n < 2*t then
        Error ( "<n> must be greater than or equal to 2*<t>" );
    fi;
    
    eps := ((t+2)*(t+1)/2) * IntCeiling(((n-t+1)*(n-t)/2)/((t+2)*(t+1)/2));
    eps := eps - ((n-t+1)*(n-t)/2);
    num := SphereContent(n,2) - ((t+2)*(t+1)/2) + eps;
    num := num * 2^n;
    den := SphereContent(n,2) - ((t+2)*(t+1)/2);
    den := den * SphereContent(n,t) + eps * SphereContent(n,t-2);
    
    return IntCeiling(num/den);
end;

#############################################################################
##
#F  LowerBoundCoveringRadiusExcess3( <n>, <r> ) Counting Excess bound
##
LowerBoundCoveringRadiusExcess3 := function ( n, t )
    local eps, rho, num, den;
    
    if t < 2 then
        Error( "The covering radius <t> must be at least 2" );
    fi;
    if n <= 2*t then
        Error( "<n> must be greater than 2*<t>+1" );
    fi;
    
    eps := (t+1) * IntCeiling((n+1)/(t+1)) - (n+1);
    if t=2 then
        rho := n - 3 + 2/n;
    else
        rho := n - t - 1;
    fi;
    num := (rho + eps) * 2^n;
    den := rho * SphereContent(n,t) + eps * SphereContent(n, t-1);
    
    return IntCeiling(num/den);
end;

#############################################################################
##
#F  LowerBoundCoveringRadiusInduction1( <n>, <r> ) Induction bound
##
LowerBoundCoveringRadiusInduction1 := function (t)
    if t < 1 then
        Error( "The covering radius <t> must be greater than or equal to 1");
    fi;
    
    return 4;
end;

#############################################################################
##
#F  LowerBoundCoveringRadiusInduction2( <n>, <r> ) Induction bound
##
LowerBoundCoveringRadiusInduction2 := function ( t )
    if t < 1 then
        Error( "The covering radius <t> must be greater than or equal to 1");
    fi;
    
    return 7;
end;

#############################################################################
##
#F  LowerBoundCoveringRadiusInduction3( <n>, <r> ) Induction bound
##
LowerBoundCoveringRadiusInduction3 := function ( t )
    if t < 4 then
        Error( "The covering radius <t> must be greater than or equal to 4");
    fi;
    
    return 8;
end;

########################################################################
##
#F  LowerBoundCoveringRadiusSphereCovering( <n>, <r> [, <F> ] )
##
## Compute a lower bound for the size of a code with length n
## and covering radius t over the field (with size) F.
## This function corresponds with propositions 2.1 and 2.13

LowerBoundCoveringRadiusSphereCovering := function ( arg )
    local n, t, q;
    
    # correct usage ?
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "usage: LowerBoundSphereCovering( <n>, <t> [, <F>] )" );
    fi;
    
    if Length( arg ) = 2 then   
        # no 3rd argument, assume GF(2)
        q := 2;
    else
        if IsInt( arg[ 3 ] ) then
            # 3rd argument is Int
            q := arg[ 3 ];
        else
            if IsField( arg[ 3 ] ) then
                # 3rd argument is Field, take size
                q := Size( arg[ 3 ] );
                if not IsInt( q ) then
                    Error( "the field must be finite" );
                fi;
            else
                # 3rd argument is wrong
                Error( "third argument must be a field ", 
                       "or the size of a field" );
            fi;
        fi;
    fi;
    if q < 2 then
        Error( "the size of the field must be greater than 1" );
    fi;
    
    n := arg[ 1 ];
    t := arg[ 2 ];
    
    # check the first two arguments
    if not IsInt( n ) or n <= 0 then 
        Error( "the length must be a positive integer" );
    fi;
    if not IsInt( t ) or t < 0 then
        Error( "the covering radius must be a non-negative integer" );
    fi;
    if t > n then
        Error( "the length must be greater than or equal to ",
               "the covering radius" );
    fi;
    
    # calculate the sphere-covering (aka sphere-packing 
    # or Johnson) bound
    return IntCeiling( ( q^n ) / ( SphereContent( n, t, q ) ) );
end;

########################################################################
##
#F  GeneralLowerBoundCoveringRadius( <n>, <size> )
##  GeneralLowerBoundCoveringRadius( <C> )
##

GeneralLowerBoundCoveringRadius := function ( arg )

    local n, size;

    if Length( arg ) = 1 and IsCode( arg[ 1 ] ) then
        n := WordLength( arg[ 1 ] );
        size := Size( arg[ 1 ] );
    elif Length( arg ) = 2 then
        n := arg[ 1 ];
        size := arg[ 2 ];
    else
        Error( "usage: GeneralLowerBoundCoveringRadius( <C> | <n>, <size> )" );
    fi;
    
    if Length( arg ) = 1 and IsBound( arg[ 1 ].boundsCoveringRadius ) then
        return arg[ 1 ].boundsCoveringRadius[ 1 ];
    fi;
    
    return Maximum(
                   LowerBoundCoveringRadiusVanWee        ( n, size, false ),
                   LowerBoundCoveringRadiusEmbedded1     ( n, size, false )
                   # LowerBoundCoveringRadiusEmbedded2     ( n, size, false ),
                   # LowerBoundCoveringRadiusExcess1       ( n, size, false ),
                   # LowerBoundCoveringRadiusExcess2       ( n, size, false ),
                   # LowerBoundCoveringRadiusExcess3       ( n, size, false ),
                   # LowerBoundCoveringRadiusInduction1    ( n, size, false ),
                   # LowerBoundCoveringRadiusInduction2    ( n, size, false ),
                   # LowerBoundCoveringRadiusInduction3    ( n, size, false ),
                   # LowerBoundCoveringRadiusSphereCovering( n, size, false ),
                   );
end;

########################################################################
##
#F  UpperBoundCoveringRadiusDelsarte( <C> )
##
UpperBoundCoveringRadiusDelsarte := function ( C )
    return WeightVector(VectorCodeword(CodeMacWilliamsTransform(C)));
end;

########################################################################
##
#F  GeneralUpperBoundCoveringRadius( <C> )
##

GeneralUpperBoundCoveringRadius := function ( C )
    
    if IsBound( C.boundsCoveringRadius ) then
        return C.boundsCoveringRadius[ Length( C.boundsCoveringRadius ) ];
    else
        return Redundancy( C );
    fi;
end;
        

