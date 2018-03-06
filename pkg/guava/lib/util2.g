########################################################################
##
#F  AllOneVector( <n> [, <field> ] )
##
##  Return a vector with all ones.
##

AllOneVector := function ( arg )
    
    local n, field;
    
    if not( Length( arg ) in [ 1, 2 ] ) then
        Error( "usage: AllOneVector( <n> [, <field> ] )" );
    fi;
    
    n := arg[ 1 ];
    if not IsInt( n ) then
        Error( "AllOneVector: <n> must be an integer" );
    fi;
    if n <= 0 then
        Error( "AllOneVector: <n> must be a positive integer" );
    fi;
    
    if Length( arg ) = 2 then
        field := arg[ 2 ];
        if IsInt( field ) then
            field := GF(field);
        fi;
        if not IsField( field ) then
            Error( "AllOneVector: <field> must be a field" );
        fi;
    else
        field := Rationals; 
    fi;
    
    return List( [ 1 .. n ], x -> field.one );
end;

########################################################################
##
#F  AllOneCodeword( <n>, <field> )
##
##  Return a codeword with <n> ones.
##

AllOneCodeword := function ( arg )
    
    local n, field;
    
    if not( Length( arg ) in [ 1, 2 ] ) then
        Error( "usage: AllOneCodeword( <n> [, <field> ]" );
    fi;
    
    n := arg[ 1 ];
    if not IsInt( n ) then
        Error( "AllOneCodeword: <n> must be an integer" );
    fi;
    if n <= 0 then
        Error( "AllOneCodeword: <n> must be a positive integer" );
    fi;
    
    if Length( arg ) = 2 then
        field := arg[ 2 ];
        if IsInt( field ) then
            field := GF(field);
        fi;
        if not IsField( field ) then
            Error( "AllOneCodeword: <field> must be a field" );
        fi;
    else
        # standard field for codes and codewords is GF(2)
        field := GF(2);
    fi;
    
    return Codeword( AllOneVector( n, field ), field );
end;


#############################################################################
##
#F  IntCeiling( <r> )
##
##  Return the smallest integer greater than or equal to r.
##  3/2 => 2,  -3/2 => -1.
##

IntCeiling := function ( r )
    
    if not IsRat( r ) then
        Error( "IntCeiling: <r> must be a rational number" );
    fi;
    
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
##

IntFloor := function ( r )
    
    if not IsRat( r ) then
        Error( "IntFloor: <r> must be a rational number" );
    fi;
    
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
##

KroneckerDelta := function ( i, j )
    
    if i = j then
        return 1;
    else
        return 0;
    fi;
    
end;


########################################################################
##
#F  BinaryRepresentation( <elements>, <length> )
##
##  Return a binary representation of an element
##  of GF( 2^k ), where k <= length.
##  
##  If elements is a list, then return the binary
##  representation of every element of the list.
##
##  This function is used to make to Gabidulin codes.
##  It is not intended to be a global function, but including
##  it in all five Gabidulin codes is a bit over the top
##
##  Therefore, no error checking is done.
##

BinaryRepresentation := function ( elementlist, length )
    
    local field, i, log, vector, element;
    
    if IsList( elementlist ) then
        return( List( elementlist,
                      x -> BinaryRepresentation( x, length ) ) );
    else
        
        element := elementlist;
        field := Field( element );

        vector := NullVector( length, GF(2) );
    
        if element = field.zero then
            # exception, log is not defined for zeroes
            return vector;
        else
            log := LogFFE( element ) + 1;
        
            for i in [ 1 .. length ] do
                if log >= 2^( length - i ) then
                    vector[ i ] := GF(2).one;
                    log := log - 2^( length - i );
                fi;
            od;
        
            return vector;
        fi;
    fi;
end;


########################################################################
##
#F  SortedGaloisFieldElements( <size> )
##
##  Sort the field elements of size <size> according to
##  their log.
##
##  This function is used to make to Gabidulin codes.
##  It is not intended to be a global function, but including
##  it in all five Gabidulin codes is not a good idea.
##

SortedGaloisFieldElements := function ( size )
    
    local field, els, sortlist;
    
    if IsInt( size ) then
        field := GF( size );
    else
        field := size;
        size := Size( field );
    fi;
    
    els := Elements( field );
    sortlist := NullVector( size );
    # log 1 = 0, so we add one to each log to avoid
    # conflicts with the 0 for zero.

    sortlist := List( els, function( x )
        if x = field.zero then
            return 0;
        else
            return LogFFE( x ) + 1;
        fi;
        end );

    sortlist{ [ 2 .. size ] } := List( els { [ 2 .. size ] },
                                       x -> LogFFE( x ) + 1 );
    SortParallel( sortlist, els );
    
    return els;
end;

