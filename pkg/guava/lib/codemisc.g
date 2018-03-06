########################################################################
##
#F  CodeWeightEnumerator( <code> )
##
##  Returns a polynomial over the rationals 
##  with degree not greater than the length of the code.
##  The coefficient of x^i equals
##  the number of codewords of weight i.
##

CodeWeightEnumerator := function ( code )
    
    local x;
    
    # these lines make sure that the polynomial is
    # printed like:
    #    x^2 + 3*x + 1
    # instead of:
    #    X(Rationals)^2 + 3*X(Rationals) + 1
    
    x := Indeterminate( Rationals );
    x.name := "x";
    
    if not IsCode( code) then
        Error( "<code> must be a code." );
    fi;
    
    return Polynomial( Rationals, WeightDistribution( code ) );
    
end;

########################################################################
##
#F  CodeDistanceEnumerator( <code>, <word> )
##
##  Returns a polynomial over the rationals
##  with degree not greater than the length of the code.
##  The coefficient of x^i equals 
##  the number of codewords with distance i to <word>.

CodeDistanceEnumerator := function ( code, word )
    
    local x;
    
    # these lines make sure that the polynomial is
    # printed like:
    #    x^2 + 3*x + 1
    # instead of:
    #    X(Rationals)^2 + 3*X(Rationals) + 1
    
    x := Indeterminate( Rationals );
    x.name := "x";
    
    if not IsCode( code ) then
        Error( "<code> must be a code.");
    fi;
    
    word := Codeword( word, code );
    
    return Polynomial( Rationals, DistancesDistribution( code, word ) );
    
end;

########################################################################
##
#F  CodeMacWilliamsTransform( <code> )
##
##  Returns a polynomial with the weight
##  distribution of the dual code as
##  coefficients.
##

CodeMacWilliamsTransform := function ( code )
    
    local weightdist, transform, size, n, x, i, j, tmp;
    
    if not IsCode( code ) then
        Error( "CodeMacWilliamsTransform: <code> must be a code" );
    fi;
    
    # these lines make sure that the polynomial is
    # printed like:
    #    x^2 + 3*x + 1
    # instead of:
    #    X(Rationals)^2 + 3*X(Rationals) + 1
    
    x := Indeterminate( Rationals );
    x.name := "x";
    
    n := WordLength( code );
    
    # if dimension < n/2, or if non-linear code, 
    # use weightdistribution of code,
    # else use weightdistribution of dual code
    if not IsLinearCode( code ) or Dimension( code ) < n / 2 then
        weightdist := WeightDistribution( code );
        size := Size( code );
        transform := List( [ 1 .. n+1 ], x -> 0 );
        for j in [ 0 .. n ] do
            tmp := 0;
            for i in [ 0 .. n ] do
                tmp := tmp + weightdist[ i+1 ] * Krawtchouk( j, i, n, 2 );
            od;
            transform[ j+1 ] := tmp / size;
        od;
    else
        transform := WeightDistribution( DualCode( code ) );
    fi;

    return Polynomial( Rationals, transform );
end;

########################################################################
##
#F  WeightVector( <vector> )
##
##  Returns the number of non-zeroes in a vector.

WeightVector := function ( vector )
    
    local pos, number, fieldzero;
    
    if not IsVector( vector ) then
        Error( "WeightVector: <vector> must be a vector" );
    fi;
    
    number := 0;
    fieldzero := Field( vector ).zero;
    
    for pos in [ 1 .. Length( vector ) ] do
        
        if vector[ pos ] <> fieldzero then
            number := number + 1;
        fi;
        
    od;
    
    return number;
    
end;


########################################################################
##
#F  RandomVector( <len> [, <weight> [, <field> ] ] )
##

RandomVector := function( arg )
    
    local len, field, wt, vec, coord, coordlist, elslist, i;
    
    if not( Length( arg ) in [ 1 .. 3 ] ) then
        Error( "usage: RandomVector( <length> [, <weight> [, <field> ] ] )" );
    fi;
    
    len := arg[ 1 ];
    if not IsInt( len ) or len <= 0 then 
        Error( "RandomVector: length must be a positive integer" );
    fi;
    
    if Length( arg ) > 1 then
        wt := arg[ 2 ];
        if not IsInt( wt ) or wt < 0 or wt > len then
            Error( "RandomVector: <weight> must be an integer in the range",
                   " 0 .. ", len );
        fi;
    else
        wt := -1;
    fi;
    
    if Length( arg ) > 2 then
        field := arg[ 3 ];
        if IsInt( field ) then
            field := GF( field );
        fi;
        if not IsField( field ) then
            Error( "RandomVector: <field> must be a field" );
        fi;
    else
        field := GF(2);
    fi;
    
    vec := NullVector( len, field );
    if wt > 0 then
        coordlist := [ 1 .. len ];
        elslist := Difference( Elements( field ), [ field.zero ] );
        # make wt elements of the vector non-zero,
        # choosing uniformly between the other field-elements
        for i in [ 1 .. wt ] do
            coord := RandomList( coordlist );
            SubtractSet( coordlist, [ coord ] );
            vec[ coord ] := RandomList( elslist );
        od;
    # do nothing if w = 0
    elif wt = -1 then
        # for each coordinate, choose uniformly from
        # all field elements, including zero
        elslist := Elements( field );
        for i in [ 1 .. len ] do
            vec[ i ] := RandomList( elslist );
        od;
    fi;

    return vec;
end;


########################################################################
##
#F  IsSelfComplementaryCode( <code> )
##
##  Return true if <code> is a complementary code, false otherwise.
##  A code is called complementary if for every v \in <code>
##  also 1 - v \in <code> (where 1 is the all-one word).
##

IsSelfComplementaryCode := function ( code )

    if not IsCode( code ) then
        Error( "IsSelfComplementaryCode: <code> is not a code" );
    fi;
    if Characteristic( code ) <> 2 then
        Error( "IsSelfComplementaryCode: <code> is not a binary code" );
    fi;

    if not IsBound( code.isSelfComplementaryCode ) then
        code.isSelfComplementaryCode :=
          code.operations.IsSelfComplementaryCode( code );
    fi;
    return code.isSelfComplementaryCode;
end;

CodeOps.IsSelfComplementaryCode := function ( code )
    
    local size, els, selfcompl, alloneword, newword; 
    
    if IsLinearCode( code ) then
        return IsSelfComplementaryCode( code );
    else
        els := Elements( code );
        selfcompl := true;
        alloneword := AllOneCodeword( WordLength( code ), GF(2) );
        while Length( els ) > 0 and selfcompl = true do
            newword := alloneword - els[ 1 ];
            if newword <> els[ 1 ] then
                if newword in els then
                    els := Difference( els, [ newword ] );
                else
                    selfcompl := false;
                fi;
                els := Difference( els, [ els[ 1 ] ] );
            fi;
        od;
        return selfcompl;
    fi;
end;

LinCodeOps.IsSelfComplementaryCode := function ( code )

    return( AllOneCodeword( WordLength( code ), GF(2) ) in code );

end;

CycCodeOps.IsSelfComplementaryCode :=
  LinCodeOps.IsSelfComplementaryCode;


########################################################################
##
#F  IsAffineCode( <code> )
##
##  Return true if <code> is affine, i.e. a linear code or
##  a coset of a linear code, false otherwise.
##

IsAffineCode := function ( code )
    
    if not IsCode( code ) then
        Error( "IsAffineCode: <code> must be a code" );
    fi;
    
    if not IsBound( code.isAffineCode ) then
        code.isAffineCode := code.operations.IsAffineCode( code );
    fi;
    
    return code.isAffineCode;
    
end;

CodeOps.IsAffineCode := function( code )
    
    if IsLinearCode( code ) then
        return IsAffineCode( code );
    elif NullWord( WordLength( code ) ) in code then
        # code cannot be a coset code of a linear code
        return false;
    elif not ( Size( code ) in List( [ 0 .. WordLength( code ) ],
            x -> Characteristic( code ) ^ x ) ) then
        # the code must have a "dimension"
        return false;
    else
        # subtract the first codeword from all codewords.
        # if the resulting code is linear, then the
        # original code is affine.
        return IsLinearCode( 
                       CosetCode( code, NullWord( WordLength( code ) ) 
                               - CodewordNr( code, 1 ) ) );
    fi;
    
end;

LinCodeOps.IsAffineCode := function( code )
    
    # all linear codes are affine
    return true;
    
end;

# and so are all cyclic codes
CycCodeOps.IsAffineCode := LinCodeOps.IsAffineCode;

        
########################################################################
##
#F  IsAlmostAffineCode( <code> )
##
##  Return true if <code> is almost affine, false otherwise.
##  A code is called almost affine if the size of any punctured
##  code is equal to q^r for some integer r, where q is the
##  size of the alphabet of the code.
##

IsAlmostAffineCode := function( code )
    
    local n, i, j, subcode, sizelist, coordlist, almostaffine;
    
    if not IsCode( code ) then
        Error( "IsAlmostAffineCode: <code> must be a code" );
    fi;
    
    if not IsBound( code.isAlmostAffineCode ) then
        if IsAffineCode( code ) then
            # every affine code is also almost affine
            almostaffine := true;
        else
            # not affine
            almostaffine := true;
            
            # however, any code over GF(2) or GF(3) is affine
            # if it is almost affine.
            # so non-affine codes with q=2,3 are also not almost affine
            if Size( Field( code ) ) = 2 
               or Size( Field( code ) ) = 3 then
                almostaffine := false;
            fi;
            
            n := WordLength( code );
            sizelist := List( [ 0 .. n ], x -> Characteristic( code ) ^ x );

            # first check whether the code itself is of size q^r
            if not ( Size( code ) in sizelist ) then
                almostaffine := false;
            else
                # now check for all possible puncturings
                i := 1;
                while almostaffine and i < n do
                    coordlist := List( Tuples( [ 1 .. n ], i ), 
                                       x -> Difference( [ 1 .. n ], x ) );
                    j := 1;
                    while almostaffine 
                      and j < Length( coordlist ) do
                        subcode := PuncturedCode( code, coordlist[ j ] );
                        # one fault is enough !
                        if not Size( subcode ) in sizelist then
                            almostaffine := false;
                        fi;
                        j := j + 1;
                    od;
                    i := i + 1;
                od;
            fi;
        fi;
        code.isAlmostAffineCode := almostaffine;
    fi;
    return code.isAlmostAffineCode;
end;


########################################################################
##
#F  IsGriesmerCode( <code> )
##
##  Return true if <code> is a Griesmer code, i.e. if
##  n = \sum_{i=0}^{k-1} d/(q^i), false otherwise.
##

IsGriesmerCode := function( code )
    
    local n, k, d, q;
    
    if not IsCode( code ) or not IsLinearCode( code ) then
        Error( "IsGriesmerCode: <code> must be a linear code" );
    fi;
    
    n := WordLength( code );
    k := Dimension( code );
    d := MinimumDistance( code );
    q := Size( Field( code ) );
    return n = Sum( [ 0 .. k-1 ], x -> IntCeiling( d / q^x ) );
end;


########################################################################
##
#F  CodeDensity( <code> )
##
##  Return the density of <code>, i.e. M*V_q(n,r)/(q^n).
##

CodeDensity := function ( code )
    
    local n, q, cr;

    if not IsCode( code ) then
        Error( "CodeDensity: <code> must be a code" );
    fi;

    cr := CoveringRadius( code );
    if not IsInt( cr ) then
        Error( "CodeDensity: the covering radius of <code> is unknown" );
    fi;
    
    n := WordLength( code  );
    q := Size( Field( code ) );
    return Size( code )
           * SphereContent( n, CoveringRadius( code ), q )
           / q^n;
    
end;


########################################################################
##
#F  DecreaseMinimumDistanceUpperBound( <C>, <s>, <iteration> )
##
##  Tries to compute the minimum distance of C.
##  The algorithm is Leon's, see for more
##  information his article.

DecreaseMinimumDistanceUpperBound := function ( arg )
    
    local C, s, iteration,
          
          trials,    # the number of trials so far
          n, k,      # some parameters of the code C
          genmat,    # the generator matrix of C
          d,         # the minimum distance so far
          continue,  # have we computed enough trials ?      
          N,         # the set { 1, ..., n }
          S,         # a random s-subset of N
          h, i, j,   # some counters
          sigma,     # permutation, mapping of N, mapping S on {1,...,s}
          tau,       # permutation, for eliminating first s columns of Emat
          Emat,      # genmat ^ sigma
          Dmat,      # (k-e,n-s) right lower submatrix of Emat
          e,         # rank of k * s submatrix of Emat
          nullrow,   # row of zeroes, for appending to Emat
          res,       # result from PutSemiStandardForm
          w,         # runs through all words spanned by Dmat
          t,         # weight of the current codeword
          v,         # word with current lowest weight
          Bmat,      # (e, n-s) right upper submatrix of Emat
          Bsupp,     # supports of differences of rows of Bmat
          Bweight,   # weights of rows of Bmat
          sup1, sup2,# temporary variables holding supports
          Znonempty, # true if e < s, false otherwise (indicates whether
                     # Zmat is a real matrix or not
          Zmat,      # ( s-e, e ) middle upper submatrix of Emat
          Zweight,   # weights of differences of rows of Zmat
          wsupp,     # weight of the current codeword w of D
          ij1,       # 0: i<>1 and j<>1 1: i=1 xor j=1 2: i=1 and j=1
          t,         # weight of current codeword
          nullw,     # nullword of length s, begin of w
          PutSemiStandardForm,   # local function for partial Gaussian
                                 # elimination
          sups,      # the supports of the elements of B
          found;     # becomes true if a better minimum distance is
                     # found
    
    if Length( arg ) <> 3 then
        Error( "usage: DecreaseMinimumDistanceUB( <C>, <s>, <iteration> )" );
    fi;
    
    
    C := arg[ 1 ];           # the code to compute the min. dist. for
    s := arg[ 2 ];           # the parameter to help find words with
                             # small weight
    iteration := arg [ 3 ];  # number of iterations to perform
    
    # check the arguments
    if not IsCode( C ) then
        Error( "DecreaseMinimumDistanceUB: <C> must be a code." );
    fi;
    if not IsLinearCode( C ) then
        Error( "DecreaseMinimumDistanceUB: <C> must be a linear code." );
    fi;
    if not IsInt( s ) then
        Error( "DecreaseMinimumDistanceUB: <s> must be an integer." );
    fi;
    if s < 1 or s > Dimension( C ) then
        Error( "DecreaseMinimumDistanceUB: <s> must lie between 1 and the ",
               "dimension of <C>." );
    fi;
    if not IsInt( iteration ) then
        Error( "DecreaseMinimumDistanceUB: <iteration> must be an integer." );
    fi;
    if iteration < 1 then 
        Error( "DecreaseMinimumDistanceLB: <iteration> must be at least zero." );
    fi;

    # the function PutSemiStandardForm is local
    #############################################################################
    ##
    #F  PutSemiStandardForm( <mat>, <s> )
    ##
    ##  Put first s coordinates of mat in standard form.
    ##  Return e as the rank of the s x s left upper
    ##  matrix. The coordinates s+1, ..., n are not permuted.
    ##
    ##  This function is based on PutStandardForm.
    ##
    ##  (maybe it's better to make this function local
    ##   in DecreaseMinimumDistanceUpperBound)
    ##

    PutSemiStandardForm := function ( mat, s )

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
                            fi;  # if mat[h][i] <> zero
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

    n := WordLength( C );
    k := Dimension( C );
    genmat := GeneratorMat( C );
    
    # step 1. initialisation
    trials := 0;
    d := n;
    continue := true; 
    found := false;

    while continue do
        
        # step 2. 
        trials := trials + 1;
        InfoMinimumDistance( "Trial nr. ", trials, "   distance: ", d, "\n" );
        
        # step 3.  choose a random s-elements subset of N
        N := [ 1 .. WordLength( C ) ];
        S := [ ];
        for i in [ 1 .. s ] do
            S[ i ] := Random( N );  # pick a random element from N
            RemoveSet( N, S[ i ] ); # and remove it from N
        od;
        Sort( S );                  # not really necessary, but
                                    # it doesn't hurt either
        
        # step 4.  choose a permutation sigma of N,
        #          mapping S onto { 1, ..., s }
        Append( S, N );
        sigma := PermList( S ) ^ (-1);
        
        # step 5.  Emat := genmat^sigma (genmat is the generator matrix C)
        Emat := [ ];
        for i in [ 1 .. k ] do
            Emat[ i ] := Permuted( genmat[ i ], sigma );
        od;
        
        # step 6.  apply elementary row operations to E
        #          and perhaps a permutation tau so that
        #          we get the following form:
        #          [ I | Z | B ]
        #          [ 0 | 0 | D ]
        #          where I is the e * e identity matrix,
        #          e is the rank of the k * s left submatrix of E
        #          the permutation tau leaves { s+1, ..., n } fixed
        
        InfoMinimumDistance( "Gaussian elimination of E ... \n");
        
        res := PutSemiStandardForm( Emat, s );
        e := res[ 1 ];    # rank (in most cases equal to s)
        tau := res[ 2 ];  # permutation of { 1, ..., s } 
        
        # append null-row to Emat (at front)
        nullrow := NullMat( 1, n, GF(2) );
        Append( nullrow, Emat );
        Emat := nullrow;
        
        InfoMinimumDistance( "Gaussian elimination of E ... done. \n" );
        
        # retrieve Dmat from Emat
        Dmat := [ ];
        for i in [ e + 1 .. k ] do
            Dmat[ i - e ] := List( [ s+1 .. n ], x -> Emat[ i+1 ][ x ] );
        od;
        
        # retrieve Bmat from Emat
        # we only need the support of the differences of the
        # rows of B
        Bmat := [ ];
        Bmat[ 1 ] := NullVector( n-s, GF(2) );
        for j in [ 2 .. e+1 ] do
            Bmat[ j ] := List( [ s+1 .. n ], x -> Emat[ j ][ x ] );
        od;
        
        InfoMinimumDistance( "Computing supports of B  ... \n" );
        sups := List( [ 1 .. e+1 ], x -> Support( Codeword( Bmat[ x ] ) ) );
        
        # compute supports of differences of rows of Bmat
        # and the weights of these supports
        # do this once every trial, instead of for each codeword,
        # to save time
        Bsupp := List( [ 1 .. e ], x -> [ ] );
        Bweight := List( [ 1 .. e ], x -> [ ] );
        for i in [ 1 .. e ] do
            sup1 := sups[ i ];
#            Bsupp[ i ] := List( [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ],
#                                x -> Difference( Union( sup1, sups[ x ] ),
#                                        Intersection( sup1, sups[ x ] ) ) );
            
            for j in [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ] do
                sup2 := sups[ j ];
                Bsupp[ i ][ j ] := Difference( Union( sup1, sup2 ),
                                               Intersection( sup1, sup2 ) );
                Bweight[ i ][ j ] := Length( Bsupp[ i ][ j ] );
            od;
        od;
        InfoMinimumDistance( "Computing supports of B  ... done. \n" );
        
        # retrieve Zmat from Emat
        # in this case we only need the weights of the supports of
        # the differences of the rows of Zmat
        # because we don't have to add them to codewords
        
        if e < s then
            
            InfoMinimumDistance( "Computing weights of Z   ... \n" );
            Znonempty := true;
            Zmat := List( [ 1 .. e ], x -> [ ] );
            Zmat[ 1 ] := NullVector( s-e, GF(2) );
            for i in [ 2 .. e+1 ] do
                Zmat[ i ] := List( [ e+1 .. s ], x -> Emat[ i ][ x ] );
            od;
            Zweight := List( [ 1 ..e ], x -> [ ] );
            for i in [ 1 .. e ] do
                for j in [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ] do
                    Zweight[ i ][ j ] := 
                      WeightCodeword( Codeword( Zmat[ i ] + Zmat[ j ] ) );
                od;
            od;
            InfoMinimumDistance( "Computing weights of Z   ... done. \n" );
            
        else
            Znonempty := false;
        fi;
              
        # step 7.  for each w in (n-s, k-e) code spanned by D
        for w in Elements( GeneratorMatCode( Dmat, GF(2) ) ) do
            wsupp := Support( w );
            
            # step 8.
            for i in [ 1 .. e ] do
                
                # step 9.
                for j in [ i + 1 - KroneckerDelta( i, 1 ) .. e+1 ] do
                    
                    ij1 := KroneckerDelta( i, 1 ) + KroneckerDelta( j, 1 );
                    
                    # step 10.
                    if Znonempty then
                        t := Zweight[ i ][ j ];
                    else
                        t := 0;
                    fi;
                    
                    # step 11.
                    if t <= ij1 then
                        
                        # step 12.
                        t := t  
                             + Bweight[ i ][ j ]
                             + Length( wsupp ) 
                             - 2 * Length( Intersection(
                                     Bsupp[ i ][ j ], wsupp ) );
                        t := t + ( 2 - ij1 );
                        if 0 < t and t < d then
                            
                            found := true;
                            
                            # step 13.
                            d := t;
                            C.upperBoundMinimumDistance :=
                              Minimum( UpperBoundMinimumDistance( C ), t );
                            # step 14.
                            nullw := NullVector( s, GF(2) );
                            Append( nullw, VectorCodeword( w ) );
                            v := Emat[ i ] + Emat[ j ] + nullw;
                            v := Permuted( v, tau ^ (-1) );
                            v := Permuted( v, sigma ^ (-1) );
                        fi;
                    fi;
                od;
            od;
        od;
        if iteration <= trials then
            continue := false;
        fi;
    od;
    if found then
        return rec( mindist := d, word := v );
    fi;
end;

