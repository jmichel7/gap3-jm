########################################################################
##
#F  AmalgatedDirectSumCode( <C>, <D>, [, <check> ] )
##
##  Return the amalgated direct sum code of C en D.
##  
##  This construction is derived from the direct sum construction,
##  but it saves one coordinate over the direct sum.
##  
##  The amalgated direct sum construction goes as follows:
##
##  Put the generator matrices G and H of C respectively D
##  in standard form as follows:
##     
##     G => [ G' | I ]     and    H => [ I | H' ]
##
##  The generator matrix of the new code then has the following form:
##     
##      [          1 0 ... 0   0 | 0 0 ............ 0 ]
##      [          0 1 ... 0   0 | 0 0 ............ 0 ] 
##      [          .........   . | .................. ] 
##      [   G'     0 0 ... 1   0 | 0 0 ............ 0 ]
##      [                    |---------------|--------]
##      [          0 0 ... 0 | 1 | 0 ... 0 0          ]
##      [--------|-----------|---|                    ]
##      [ 0 0 ............ 0 | 0   1 ... 0 0    H'    ]
##      [ .................. | 0   .........          ]
##      [ 0 0 ............ 0 | 0   0 ... 1 0          ]
##      [ 0 0 ............ 0 | 0   0 ... 0 1          ]
##
##  The codes resulting from [ G' | I ] and [ I | H' ] must
##  be acceptable in the last resp. the first coordinate.
##  Checking whether this is true takes a lot of time, however,
##  and is only performed when the boolean variable check is true.
##

AmalgatedDirectSumCode := function ( arg )
    
    local C, D, check,
          LinearAmalgatedDirectSum,
          UnrestrictedAmalgatedDirectSum;

    # check the arguments
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "Usage: AmalgatedDirectSumCode( <C>, <D>, [, <check> ] )" );
    fi;
    if Length( arg ) = 2 then
        check := false;
    else
        check := arg[ 3 ];
        if check <> true and check <> false then
            Error( "AmalgatedDirectSumCode: <check> must be true or false" );
        fi;
    fi;
    C := arg[ 1 ];
    D := arg[ 2 ];
    if not IsCode( C ) or not IsCode( D ) then
        Error( "AmalgatedDirectSumCode: <C> and <D> must be codes" );
    fi;
    if Field( C ) <> Field( D ) then
        Error( "AmalgatedDirectSumCode: <C> and <D> must be codes ",
                "over the same field" );
    fi;

    LinearAmalgatedDirectSum := function( C, D, check )

        local G, H, Cstandard, Dstandard, NewG, NewH, Nulmat,
              NewC, i;
        G := Copy( GeneratorMat( C ) );
        H := Copy( GeneratorMat( D ) );
        # standard form: G => [ G' | I ]
        PutStandardForm( G, false );
        # standard form: H => [ I | H' ]
        PutStandardForm( H, true );

        # check whether the construction is allowed;
        # this is (at this time) a lot of work
        # maybe it will disappear later
        if check then
            Cstandard := GeneratorMatCode( G, Field( C ) );
            Dstandard := GeneratorMatCode( H, Field( C ) );
            # is the last coordinate of the standardcode C' acceptable ?
            if not IsCoordinateAcceptable( Cstandard, WordLength( Cstandard ) ) then
                Error( "AmalgatedDirectSumCode: Standard form of <C> is not ",
                       "acceptable in the last coordinate.");
            fi;
            # is the last coordinate of the standardcode D' acceptable ?
            if not IsCoordinateAcceptable( Dstandard, 1 ) then
                Error( "AmalgatedDirectSumCode: Standard form of <D> is not ",
                       "acceptable in first coordinate.");
            fi;
        fi;

        NewG := G;
        # build upper part of the new generator matrix
        # append n_D - 1 zeroes to all rows of G'
        for i in NewG do
            Append( i, List( [ 1..WordLength( D ) - 1 ], x -> Field(C).zero ) );
        od;
        # concatenate the last row of G' and the first row of H'
        for i in [ 2 .. WordLength( D ) ] do
            NewG[ Dimension( C ) ][ WordLength( C ) + i - 1 ] :=
              H[ 1 ][ i ];
        od;
        # throw away the first row of H' (it is already appended)
        NewH := List( [ 2..Length( H ) ], x -> H[ x ] );
        # throw away the first column of H' (it is (1, 0, ..., 0),
        # the one is already present in the new generator matrix)
        NewH := List( [ 1..Length( NewH ) ], x -> List(
                     [ 2 .. WordLength( D ) ], y -> NewH[ x ][ y ] ) );
        # fill the lower left part with zeroes
        Nulmat := List( [ 1 .. Dimension( D ) - 1 ],
                        x -> NullVector( WordLength( C ) , Field(C) ) );
        # and put the rest of H' in the lower right corner
        for i in [ 1 .. Length( Nulmat ) ] do
            Append( Nulmat[ i ], NewH[ i ] );
        od;
        # paste the lower and the upper part
        Append( NewG, Nulmat );

        # construct the new code with the new generator matrix
        NewC := GeneratorMatCode( NewG, "amalgated direct sum code",
          Field( C ) );
        # write history
        NewC.history := MergeHistories( History( C ), History( D ) );
        NewC.history := List( NewC.history, x -> String( x ) );

        # the new covering radius is at most the sum of the
        # two old covering radii, if the old codes are normal
        if IsBound( C.isNormalCode )
               and IsBound( D.isNormalCode )
               and IsNormalCode( C )
               and IsNormalCode( D ) then
            if IsBound( C.boundsCoveringRadius )
                   and IsBound( D.boundsCoveringRadius ) then
                NewC.boundsCoveringRadius := [
                  GeneralLowerBoundCoveringRadius( NewC )
                  .. Maximum( C.boundsCoveringRadius )
                     + Maximum( D.boundsCoveringRadius ) ];
            fi;
        fi;
        return NewC;
    end;

    UnrestrictedAmalgatedDirectSum := function( C, D, check )

        local elsC, elsD, i, NewC,
              newels;

        # check whether the construction is allowed;
        # this is (at this time) a lot of work
        # maybe it will disappear later
        if check then
            # is the last coordinate of C acceptable ?
            if not IsCoordinateAcceptable( C, WordLength( C ) ) then
                Error( "AmalgatedDirectSum: <C> is not ",
                       "acceptable in the last coordinate.");
            fi;
            # is the first coordinate of D acceptable ?
            if not IsCoordinateAcceptable( D, 1 ) then
                Error( "AmalgatedDirectSum: <D> is not ",
                       "acceptable in first coordinate.");
            fi;
        fi;

        # find the elements of the new code
        newels := [];
        for i in Elements( Field( C ) ) do
            elsC := VectorCodeword( Elements(
              CoordinateSubCode( C, WordLength( C ), i ) ) );
            elsD := VectorCodeword( Elements(
              CoordinateSubCode( D, 1, i ) ) );
            if Length( elsC ) > 0 and
               Length( elsD ) > 0 then
                elsD := List( elsD,
                  x -> x{ [ 2 .. WordLength( D ) ] } );
                for i in elsC do
                    Append( newels, List( elsD,
                      x -> Concatenation( i, x ) ) );
                od;
            fi;
        od;

        if Length( newels ) = 0 then
            Error( "AmalgatedDirectSumCode: there are no ",
                   "codewords satisfying the conditions" );
        fi;
        NewC := ElementsCode( newels, "amalgated direct sum code",
          Field( C ) );

        # write history
        NewC.history := MergeHistories( History( C ), History( D ) );

        # the new covering radius is at most the sum of the
        # two old covering radii, if the old codes are normal
        if IsBound( C.isNormalCode )
               and IsBound( D.isNormalCode )
               and IsNormalCode( C )
               and IsNormalCode( D ) then
            if IsBound( C.boundsCoveringRadius )
                   and IsBound( D.boundsCoveringRadius ) then
                NewC.boundsCoveringRadius := [
                  GeneralLowerBoundCoveringRadius( NewC )
                  .. Maximum( C.boundsCoveringRadius )
                     + Maximum( C.boundsCoveringRadius ) ];
            fi;
        fi;
        return NewC;
    end;
    
    if IsLinearCode( C ) and IsLinearCode( D ) then
        # this is much faster, because it only uses
        # the generator matrices
        return LinearAmalgatedDirectSum( C, D, check );
    else
        return UnrestrictedAmalgatedDirectSum( C, D, check );
    fi;
end;

########################################################################
##
#F  BlockwiseDirectSumCode( <C1>, <L1>, <C2>, <L2> )
##
##  Return the blockwise direct sum of C1 and C2 with respect to 
##  the cosets defined by the codewords in L1 and L2.
##

BlockwiseDirectSumCode := function ( arg )
    
    local C1, C2, L1, L2, k, CosetCodeBlockwiseDirectSumCode,
          SubCodeBlockwiseDirectSumCode, unioncode1, unioncode2, i;

    # check the arguments
    if Length( arg ) <> 4 then
        Error( "Usage: BlockwiseDirectSumCode( <C1>, <L1>, <C2>, <L2> )" );
    fi;
    C1 := arg[ 1 ]; C2 := arg[ 3 ];
    L1 := arg[ 2 ]; L2 := arg[ 4 ];
    if not IsCode( C1 ) or not IsCode( C2 ) then
        Error( "BlockwiseDirectSumCode: <C1> and <C2> must be codes" );
    fi;
    if Field( C1 ) <> Field( C2 ) then
        Error( "BlockwiseDirectSumCode: the fields of <C1> and <C2>",
               "must be the same" );
    fi;
    if not IsList( L1 ) or not IsList( L2 ) then
        Error( "BlockwiseDirectSumCode: <L1> and <L2> must be lists ",
               "of codewords or of codes" );
    fi;
    if Length( L1 ) <> Length( L2 ) then
        Error( "BlockwiseDirectSumCode: <L1> and <L2>",
               "must have equal lengths" );
    fi;
    k := Length( L1 );
    if k = 0 then
        Error( "BlockwiseDirectSumCode: the lists are empty" );
    fi;
    
    CosetCodeBlockwiseDirectSumCode := function( C1, L1, C2, L2 )

        local newels, i, k, subcode1, subcode2, newcode, sum;

        k := Length( L1 );
        # make the elements of the new code
        newels := [];
        for i in [ 1 .. k ] do
            # build subcode 1, using the GUAVA-function CosetCode,
            # which adds the word L1[i] to all codewords of C1
            subcode1 := CosetCode( C1, L1[ i ] );
            subcode2 := CosetCode( C2, L2[ i ] );
            # now build the the direct sum of the two subcodes
            sum := DirectSumCode( subcode1, subcode2 );
            # add the elements of the direct sum-code to
            # the elements we already found
            # in the previous steps
            newels := Union( newels, Elements( sum ) );
        od;
        # finally build the new code with the computed elements
        newcode := ElementsCode( newels, "blockwise direct sum code",
                           Field( C1 ) );
        # write history
        newcode.history := MergeHistories( History( C1 ), History( C2 ) );
        return newcode;
    end;

    SubCodeBlockwiseDirectSumCode := function( C1, L1, C2, L2 )

        local unioncode, newels, newcode;

        newels := Elements( DirectSumCode( L1[ 1 ], L2[ 1 ] ) );
        for i in [ 2 .. Length( L1 ) ] do
            newels := Union( newels, Elements(
              DirectSumCode( L1[ i ], L2[ i ] ) ) );
        od;
        newcode := ElementsCode( newels, "blockwise direct sum code",
          Field( C1 ) );
        newcode.history := MergeHistories(
          History( C1 ), History( C2 ) );
        return newcode;
    end;

    if IsCode( L1[ 1 ] ) and IsCode( L2[ 1 ] ) then
        unioncode1 := L1[ 1 ];
        unioncode2 := L2[ 1 ];
        for i in [ 2 .. k ] do

            if not IsCode( L1[ i ] ) or not IsCode( L2[ i ] ) then
                Error( "BlockwiseDirectSumCode: not all elements of ",
                       "the lists are codes" );
            fi;
            unioncode1 := AddedElementsCode( unioncode1,
              Elements( L1[ i ] ) );
            unioncode2 := AddedElementsCode( unioncode2,
              Elements( L2[ i ] ) );
            if unioncode1 <> C1 then
                Error( "BlockwiseDirectSumCode: <C1> must be the ",
                       "union of the codes in <L1>" );
            fi;
            if unioncode2 <> C2 then
                Error( "BlockwiseDirectSumCode: <C2> must be the ",
                       "union of the codes in <L2>" );
            fi;
        od;
        return SubCodeBlockwiseDirectSumCode( C1, L1, C2, L2 );
    else
        L1 := Codeword( L1 );
        L2 := Codeword( L2 );
        return CosetCodeBlockwiseDirectSumCode( C1, L1, C2, L2 );
    fi;
end;

########################################################################
##
#F  ExtendedDirectSumCode( <L>, <B>, m )
##
##  The construction as described in the article of Graham and Sloane,
##  section V.
##  ("On the Covering Radius of Codes", R.L. Graham and N.J.A. Sloane,
##    IEEE Information Theory, 1985 pp 385-401)
##

ExtendedDirectSumCode := function ( L, B, m )
    
    local LinearExtendedDirectSumCode,
          UnrestrictedExtendedDirectSumCode;
    
    # check the arguments
    if not IsCode( L ) or not IsCode( B ) then
        Error( "L and B must be codes" );
    fi;
    if WordLength( L ) <> WordLength( B ) then
        Error( "L and B must have the same length" );
    fi;
    if not IsInt( m ) or m < 1 then
        Error( "m must be a positive integer" );
    fi;

    LinearExtendedDirectSumCode := function( L, B, m )

        local m0, GL, GB, NewG, i, j, G, firstzeros, lastzeros, newcode;

        # if L is a subset of B, then
        # skip the last mth copy of L
        # (it doesn't add anything new to the code )
        if L in B and m > 1 then
            m0 := m - 1;
        else
            m0 := m;
        fi;

        GL := Copy( GeneratorMat( L ) );
        GB := Copy( GeneratorMat( B ) );
        # the new generator matrix, fill with zeros first
        NewG := List( [ 1 .. Dimension( L ) * m0 + Dimension( B ) ],
                      x -> NullWord( WordLength( L ) * m ) );
        # first m0 * Dimension(L) rows,
        # form: [ GL 0  ... 0  ]
        #       [ 0  GL ... 0  ]
        #       [ ...........  ]
        #       [ 0  0      GL ] (this row is omitted if L in B)
        for i in [ 1 .. m0 ] do
            # construct rows (i-1)*Dimension(L) till i*Dimension(L)-1
            firstzeros := NullVector( ( i-1 ) * WordLength( L ), Field( L ) );
            lastzeros := NullVector( ( m-i ) * WordLength( L ), Field( L ) );
            for j in [ 1..Dimension( L ) ] do
                NewG[ j + ( i-1 ) * Dimension( L ) ] :=
                  Concatenation( firstzeros, GL[j], lastzeros );
            od;
        od;
        # last row of new generator matrix
        # [ GB GB GB GB ... GB ]
        for i in [ 1..Dimension( B ) ] do
            NewG[ i + m * Dimension( L ) ] :=
              Concatenation( List( [ 1..m ], x->GB[ i ] ) );
        od;
        newcode := GeneratorMatCode( NewG,
                           Concatenation( String( m ),
                                   "-fold extended direct sum code" ),
                           Field( L ) );
        newcode.history := MergeHistories( History( L ), History( B ) );

        return newcode;
    end;

    UnrestrictedExtendedDirectSumCode := function( L, B, m )

        local SumCode, i, j, el, newcode, lastels;

        SumCode := L;
        for i in [ 2 .. m ] do
            SumCode := DirectSumCode( SumCode, L );
        od;
        lastels := [];
        for i in VectorCodeword( Elements( B ) ) do
            el := Copy( i );
            for j in [ 2 .. m ] do
                el := Concatenation( el, i );
            od;
            Append( lastels, el );
        od;
        newcode := AddedElementsCode( SumCode, lastels );
        newcode.history := MergeHistories( History( L ), History( B ) );
        newcode.name := Concatenation( String( m ),
          "-fold extended direct sum code" );

        return newcode;
    end;

    if IsLinearCode( L ) and IsLinearCode( B ) then
        # this is faster because it only uses the generator matrices
        return LinearExtendedDirectSumCode( L, B, m );
    else
        return UnrestrictedExtendedDirectSumCode( L, B, m );
    fi;
end;


########################################################################
##
#F  PiecewiseConstantCode( <partition>, <constraints> [, <field> ] )
##

PiecewiseConstantCode := function ( arg )
    
    local n, partition, constraints, field, i, j, 
          elements, constr, position, newels, addels, 
          ConstantWeightElements, el, sumels;
    
    # check the arguments
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "usage: PiecewiseConstantCode( <partition>, ",
               "<constraints> [, <field> ] )" );
    fi;
    if Length( arg ) = 3 then
        field := arg[ 3 ];
        if not IsField( field ) then
            if not IsInt( field ) then
                Error( "PiecewiseConstantCode: <field> must be a field" );
            else
                field := GF( field );
            fi;
        fi;
    else
        field := GF( 2 );
    fi;
    
    # find out the partition
    partition := arg[ 1 ];
    # allow for an integer as "partition"
    if not IsList( partition ) then
        partition := [ partition ];
    fi;
    for i in partition do
        if not IsInt( i ) then
            Error( "PiecewiseConstantCode: <partition> must be ",
                   "a list of integers" );
        fi;
        if i <= 0 then
            Error( "PiecewiseConstantCode: <partition> must be ",
                   "a list of positive integers" );
        fi;
    od;
    n := Sum( partition );
    
    # check the constraints
    constraints := arg[ 2 ];
    if not IsList( constraints ) then
        constraints := [ constraints ];
    fi;
    for i in [ 1 .. Length( constraints ) ] do
        if not IsList( constraints[ i ] ) then 
            constraints[ i ]:= [ constraints[ i ] ];
        fi;
        if Length( constraints[ i ] ) <> Length( partition ) then
            Error( "PiecewiseConstantCode: the length of constraint ", i,
                   " is not equal to the length of <partition>" );
        fi;
        for j in [ 1 .. Length( constraints[ i ] ) ] do
            if not IsInt( constraints[ i ][ j ] ) then
                Error( "PiecewiseConstantCode: entry ", j,
                       " of constraint ", i,
                       " is not an integer" );
            fi;
            if constraints[ i ][ j ] < 0 
               or constraints[ i ][ j ] > partition[ j ] then
                Error( "PiecewiseConstantCode: entry ", j,
                       " of constraint ", i,
                       " must >= 0 and <= ", partition[ j ] );
            fi;
        od;
    od;
    
    ConstantWeightElements := function( place, weight )
        local els, i, numberofzeros, zerovector;
        if weight = 0 then
            return [ NullVector( n, field ) ];
        fi;
        els := VectorCodeword( Elements( 
                       ConstantWeightSubcode( WholeSpaceCode( 
                               partition[ place ], field ), weight ) ) );
        # add zeros to the front
        if place <> 1 then
            numberofzeros := Sum( partition{ [ 1 .. place - 1 ] } );
            zerovector := NullVector( numberofzeros, field );
            for i in [ 1 .. Length( els ) ] do
                els[ i ] := Flat( [ zerovector, els[ i ] ] );
            od;
        fi;
        if place <> Length( partition ) then
            numberofzeros := Sum( partition{ [ place + 1 
                                     .. Length( partition ) ] } );
            zerovector := NullVector( numberofzeros, field );
            for i in [ 1 .. Length( els ) ] do
                els[ i ] := Flat( [ els[ i ], zerovector ] );
            od;
        fi;
        return els;
    end;
    
    # make the new code
    elements := [ ];
    for constr in constraints do
        newels := [ ];
        for i in [ 1 .. Length( partition ) ] do
            addels := ConstantWeightElements( i, constr[ i ] );
            if Length( newels ) > 0 then
                sumels := [ ];
                for el in addels do
                    Append( sumels, List( newels, x -> x + el ) );
                od;
                newels := Copy( sumels );
            else
                newels := Copy( addels );
            fi;
        od;
        Append( elements, newels );
    od;
    return ElementsCode( elements, "piecewise constant code", field );
end;


########################################################################
##
#F  GabidulinCode( );
##

GabidulinCode := function ( arg )
    
    local checkmat, m, nmat, els, i, j, w1, w2, w3,
          fieldels, binels, newels, binnewels, size, one, newcode;
    
    m := arg[ 1 ];
    if arg[ Length( arg ) ] <> false and m < 4 then
        Error( "GabidulinCode: <m> must be at least 4" );
    fi;
    
    w1 := arg[ 2 ];
    w2 := arg[ 3 ];
    
    w3 := w1 + w2;
            
    checkmat := NullMat( 5 * 2^(m-2) - 1, 2 * m - 1, GF( 2 ) );
    
    size := 2^(m-2);
    
    fieldels := SortedGaloisFieldElements( size );
    binels := BinaryRepresentation( fieldels, m-2 );
    
    one := Z(2)^0;
    
    # make matrix N
    
    for i in [ 1 .. size - 1 ] do
        for j in [ 1 .. m-2 ] do
            checkmat[ i ][ j + 1 ] := binels[ i + 1 ][ j ];
        od;
    od;
    
    # make matrix D
    
    for i in [ 1 .. size ] do
        checkmat[ size - 1 + i ][ 1 ] := one;
        
        for j in [ 1 .. m - 2 ] do
            checkmat[ size - 1 + i ][ j + 1 ] := binels[ i ][ j ];
        od;
        
        newels := Copy( fieldels );
        for j in [ 2 .. size ] do
            newels[ j ] := w1/fieldels[ j ];
        od;
        binnewels := BinaryRepresentation( newels, m-2 );
        
        for j in [ 1 .. m - 2 ] do
            checkmat[ size - 1 + i ][ j + m - 1 ] := binnewels[ i ][ j ];
        od;
    od;
    
    # make matrix Q
    
    for i in [ 1 .. size ] do
        checkmat[ 2 * size - 1 + i ][ 1 ] := one;
        checkmat[ 2 * size - 1 + i ][ 2 * m - 1 ] := one;
        
        for j in [ 1 .. m - 2 ] do
            checkmat[ 2 * size - 1 + i ][ j + 1 ] := binels[ i ][ j ];
        od;
        
        newels := Copy( fieldels );
        for j in [ 2 .. size ] do
            newels[ j ] := w2/fieldels[ j ];
        od;
        binnewels := BinaryRepresentation( newels, m-2 );
        
        for j in [ 1 .. m - 2 ] do 
            checkmat[ 2 * size - 1 + i ][ j + m - 1 ] := binnewels[ i ][ j ];
        od;
    od;
    
    # make matrix M
    
    for i in [ 1 .. size ] do
        checkmat[ 3 * size - 1 + i ][ 1 ] := one;
        checkmat[ 3 * size - 1 + i ][ 2 * m - 2 ] := one;
        
        for j in [ 1 .. m - 2 ] do
            checkmat[ 3 * size - 1 + i ][ j + 1 ] := binels[ i ][ j ];
        od;
        
        newels := Copy( fieldels );
        for j in [ 2 .. size ] do
            newels[ j ] := w3/fieldels[ j ];
        od;
        binnewels := BinaryRepresentation( newels, m-2 );
        
        for j in [ 1 .. m - 2 ] do
            checkmat[ 3 * size - 1 + i ][ j + m - 1 ] := binnewels[ i ][ j ];
        od;
    od;
    
    # make matrix G
    
    for i in [ 1 .. size ] do
        checkmat[ 4 * size - 1 + i ][ 1 ] := one;
        checkmat[ 4 * size - 1 + i ][ 2 * m - 2 ] := one;        
        checkmat[ 4 * size - 1 + i ][ 2 * m - 1 ] := one;    
        
        for j in [ 1 .. m - 2 ] do
            checkmat[ 4 * size - 1 + i ][ m - 1 + j ] := binels[ i ][ j ];
        od;
    od;
    
    checkmat := TransposedMat( checkmat );
    
    newcode := CheckMatCode( checkmat, 
                       Concatenation( "Gabidulin code (m=",String(m),")" ),
                       GF(2) );
    
#    newcode.wordLength := 5 * size - 1;
#    newcode.dimension := 2 * m - 1;
    newcode.lowerBoundMinimumDistance := 3;
    newcode.upperBoundMinimumDistance := 3;
    newcode.boundsCoveringRadius := [ 2 ];
        
    return( newcode );    
end;


########################################################################
##
#F  EnlargedGabidulinCode( );
##

EnlargedGabidulinCode := function ( arg )
    
    local checkmat, m, nmat, els, i, j, w1, w2, w3,
          fieldels, binels, newels, binnewels, size, one, newcode, bmat, binel, el, k;
    
    m := arg[ 1 ];
    if m < 4 then
        Error( "EnlargedGabidulinCode: <m> must be at least 4" );
    fi;
    
    w1 := arg[ 2 ];
    w2 := arg[ 3 ];
    el := arg[ 4 ];
    
    w3 := w1 + w2;
            
    checkmat := NullMat( 7 * 2^(m-2) - 2, 2 * m, GF( 2 ) );
    
    size := 2^(m-1);
    
    fieldels := SortedGaloisFieldElements( GF( size ) );
    binels := BinaryRepresentation( fieldels, m-1 );
    
    one := Z(2)^0;
    
    # make matrix Z
    
    bmat := TransposedMat( Copy( CheckMat( 
                    GabidulinCode( m, w1, w2, false ) ) ) );

    for i in [ 1 .. Length( bmat ) - 1 ] do
        k := i;
        if k > 2^(m-2) - 1 then
            k := k + 1;
        fi;
        for j in [ 1 .. Length( bmat[ 1 ] ) ] do
            checkmat[ i ][ j + 1 ] := bmat[ k ][ j ];
        od;
    od;
    
    # make matrix Y
    
    binel := BinaryRepresentation( el, m );
    
    for i in [ 1 .. size ] do
        checkmat[ Length( bmat ) - 1 + i ][ 1 ] := one;
        for j in [ 1 .. m-1 ] do
            checkmat[ Length( bmat ) - 1 + i ][ j + 1 ] := binels[ i ][ j ];
        od;
        for j in [ 1 .. m ] do
            checkmat[ Length( bmat ) - 1 + i ][ m + j ] := binel[ j ];
        od;
    od;
        
    checkmat := TransposedMat( checkmat );
    newcode := CheckMatCode( checkmat, 
                       Concatenation( "enlarged Gabidulin code (m="
                               ,String(m),")" ),
                       GF(2) );

    newcode.lowerBoundMinimumDistance := 3;
    newcode.upperBoundMinimumDistance := 3;
    newcode.boundsCoveringRadius := [ 2 ];
        
    return( newcode );    
end;


########################################################################
##
#F  DavydovCode( );
##

DavydovCode := function ( arg )
    
    local checkmat, r, i, j, v, field2el, binel2, fieldels, binels, 
          k, l, field1el, field2el, binel1, newcode;
    
    r := arg[ 1 ];
    if r < 4 then
        Error( "DavydovCode: <m> must be at least 5" );
    fi;
    
    v := arg[ 2 ];
    i := arg[ 3 ]; # 0 < i < 2^v ??
    j := arg[ 4 ]; # 0 < j < 2^(r-v) ??

    checkmat := NullMat( 2^v + 2^(r-v) - 3, r, GF( 2 ) );
    
    field1el := i;
    binel1 := BinaryRepresentation( field1el, v );
    field2el := j;
    binel2 := BinaryRepresentation( field2el, r - v );

    
    # make matrix K
    
    fieldels := SortedGaloisFieldElements( GF( 2^v ) ); # remove 0
    SubtractSet( fieldels, [ fieldels[ 1 ], i ] );
    binels := BinaryRepresentation( fieldels, v );
    
    
    for k in [ 1 .. 2^v - 2 ] do
        for l in [ 1 .. v ] do
            checkmat[ k ][ l ] := binels[ k ][ l ];
        od;
        for l in [ 1 .. r-v ] do
            checkmat[ k ][ v + l ] := binel2[ l ];
        od;
    od;
    
    # make matrix (column) A
    
    for l in [ 1 .. v ] do
        checkmat[ 2^v - 1 ][ l ] := binel1[ l ];
    od;
    for l in [ 1 .. r-v ] do
        checkmat[ 2^v - 1 ][ v + l ] := binel2[ l ];
    od;
    
    # make matrix S
    
    fieldels := SortedGaloisFieldElements( GF( 2^(r-v) ) );
    SubtractSet( fieldels, [ fieldels[ 1 ], j ] );
    binels := BinaryRepresentation( fieldels, r - v );
    
    for k in [ 1 .. 2^(r-v) - 2 ] do
        for l in [ 1 .. v ] do
            checkmat[ 2^v - 1 + k ][ l ] := binel1[ l ];
        od;
        for l in [ 1 .. r-v ] do
            checkmat[ 2^v - 1 + k ][ v + l ] := binels[ k ][ l ];
        od;
    od;
    
    
    checkmat := TransposedMat( checkmat );
    
    newcode := CheckMatCode( checkmat, 
                       Concatenation( "Davydov code (r=",String(r),
                               ", v=", String(v), ")" ),
                       GF(2) );
    
#    newcode.wordLength := 5 * size - 1;
#    newcode.dimension := newcode.wordLength - (2 * m);
    newcode.lowerBoundMinimumDistance := 4;
    newcode.upperBoundMinimumDistance := 4;
    newcode.boundsCoveringRadius := [ 2 ];
        
    return( newcode );    
end;



########################################################################
##
#F  TombakCode( );
##

TombakCode := function ( arg )
    
    local checkmat, one, m, i, beta, gamma, delta, w1, w2, w3, newcode,
          k, size, fieldels, binels, binel, l, sortlist, sortedlist, 
          newels, binnewels, theta, newsize, betabin, gammbin, gammabin, deltabin, w1bin, w2bin, w3bin;
    
    m := arg[ 1 ];
    if arg[ Length( arg ) ] <> false and m < 5 then
        Error( "TombakCode: <m> must be at least 5" );
    fi;
    
    beta := arg[ 3 ];
    gamma := arg[ 4 ];
    delta := beta + gamma;
    
    betabin := BinaryRepresentation( beta, m-1 );
    gammabin := BinaryRepresentation( gamma, m-1 );
    deltabin := betabin + gammabin;
    
    w1 := arg[ 5 ];
    w2 := arg[ 6 ];
    w3 := w1 + w2;
    
    w1bin := BinaryRepresentation( w1, m-3 );
    w2bin := BinaryRepresentation( w2, m-3 );
    w3bin := w1bin + w2bin;
    
    i := arg[ 2 ];

    checkmat := NullMat( 15 * 2^(m-3) - 3, 2*m, GF( 2 ) );
    
    one := Z(2)^0;
           
    # make matrix C
    
    size := 2^(m-3);
    
    fieldels := SortedGaloisFieldElements( size );
    binels := BinaryRepresentation( fieldels, m-3 );
    
    binel := betabin;
    
    for k in [ 1 .. size - 1 ] do
        for l in [ 1 .. m - 3 ] do
            checkmat[ k ][ 3 + l ] := binels[ k + 1 ][ l ];
        od;
        for l in [ 1 .. m - 1 ] do
            checkmat[ k ][ m + l ] := binel[ l ];
        od;
        checkmat[ k ][ 2 * m ] := one;
    od;
    
    # make matrix V
    
    binel := gammabin;
    
    for k in [ 1 .. 2^(m-3) ] do
        checkmat[ size - 1 + k ][ 2 ] := one;
        for l in [ 1 .. m - 3 ] do
            checkmat[ size - 1 + k ][ 3 + l ] := binels[ k ][ l ];
        od;
        for l in [ 1 .. m - 1 ] do
            checkmat[ size - 1 + k ][ m + l ] := binel[ l ];
        od;
        checkmat[ size - 1 + k ][ 2 * m ] := one;
    od;
    
    # make matrix X
    
    binel := deltabin;
    
    for k in [ 1 .. size ] do
        checkmat[ 2 * size - 1 + k ][ 2 ] := one;
        checkmat[ 2 * size - 1 + k ][ 3 ] := one;
        for l in [ 1 .. m - 3 ] do
            checkmat[ 2 * size - 1 + k ][ 3 + l ] := binels[ k ][ l ];
        od;
        for l in [ 1 .. m - 1 ] do
            checkmat[ 2 * size - 1 + k ][ m + l ] := binel[ l ];
        od;
    od;
    
    # make sub-matrix Theta
    theta := NullMat( 4*size - 1, 2 * (m-3) + 2, GF( 2 ) );
    
    for k in [ 1 .. size - 1 ] do
        for l in [ 1 .. m-3 ] do
            theta[ k ][ l ] := binels[ k + 1 ][ l ];
        od;
        
        newels := Copy( fieldels );
        for l in [ 2 .. size ] do 
            newels[ l ] := w1/fieldels[ l ];
        od;
        binnewels := BinaryRepresentation( newels, m-3 );
        
        for l in [ 1 .. m-3 ] do
            theta[ k ][ m-3 + l ] := binnewels[ k + 1 ][ l ];
        od;
    od;
    
    for k in [ 1 .. size ] do
        for l in [ 1 .. m-3 ] do
            theta[ size - 1 + k ][ l ] := binels[ k ][ l ];
        od;
        
        newels := Copy( fieldels );
        for l in [ 2 .. size ] do
            newels[ l ] := w2/fieldels[ l ];
        od;
        binnewels := BinaryRepresentation( newels, m-3 );
        
        for l in [ 1 .. m-3 ] do
            theta[ size - 1 + k ][ m-3 + l ] := binnewels[ k ][ l ];
        od;
        theta[ size - 1 + k ][ 2*(m-3) + 2 ] := one;
    od;
    
    for k in [ 1 .. size ] do
        for l in [ 1 .. m-3 ] do
            theta[ 2 * size - 1 + k ][ l ] := binels[ k ][ l ];
        od;
        
        newels := Copy( fieldels );
        for l in [ 2 .. size ] do
            newels[ l ] := w3/fieldels[ l ];
        od;
        binnewels := BinaryRepresentation( newels, m-3 );
        
        for l in [ 1 .. m-3 ] do
            theta[ 2 * size - 1 + k ][ m-3 + l ] := binnewels[ k ][ l ];
        od;
        theta[ 2 * size - 1 + k ][ 2*(m-3) + 1 ] := one;
    od;
    
    for k in [ 1 .. size ] do
        for l in [ 1 .. m-3 ] do
            theta[ 3*size - 1 + k ][ m-3 + l ] := binels[ k ][ l ];
        od;
        
        theta[ 3 * size - 1 + k ][ 2*(m-3) + 1 ] := one;
        theta[ 3 * size - 1 + k ][ 2*(m-3) + 2 ] := one;
    od;
    
    # make matrix Phi
    
    for k in [ 1 .. 4*size - 1 ] do
        checkmat[ 3*size - 1 + k ][ 3 ] := one;
        for l in [ 1 .. 2*m - 4 ] do
            checkmat[ 3*size - 1 + k ][ 3 + l ] := theta[ k ][ l ];
        od;
    od;
    
    # make matrix Lambda
    
    binel := betabin;
    
    for k in [ 1 .. 4*size - 1 ] do
        checkmat[ 7 * size - 2 + k ][ 3 ] := one;
        for l in [ 1 .. m - 3 ] do
            checkmat[ 7 * size - 2 + k ][ 3 + l ] := theta[ k ][ l ];
        od;
        for l in [ 1 .. m - 1 ] do
            checkmat[ 7 * size - 2 + k ][ m-1 + 3 + l ] := 
              theta[ k ][ m-3 + l ] + binel[ l ];
        od;
        checkmat[ 7 * size - 2 + k ][ 2*m ] := one;
    od;
    
    # make matrix Y
    
    newsize := 2^(m-1);
    fieldels := SortedGaloisFieldElements( newsize );
    binels := BinaryRepresentation( fieldels, m-1 );
    
    binel := BinaryRepresentation( i, m );
    
    for k in [ 1 .. newsize ] do
        checkmat[ 11*size - 3 + k ][ 1 ] := one;
        for l in [ 1 .. m-1 ] do
            checkmat[ 11*size - 3 + k ][ 1 + l ] := binels[ k ][ l ];
        od;
        for l in [ 1 .. m ] do
            checkmat[ 11*size - 3 + k ][ m + l ] := binel[ l ];
        od;
    od;
    
    checkmat := TransposedMat( checkmat );
    
    newcode := CheckMatCode( checkmat, 
                       Concatenation( "Tombak code (m=",String(m),")" ),
                       GF(2) );
    
#    newcode.wordLength := 5 * size - 1;
#    newcode.dimension := newcode.wordLength - (2 * m);
    newcode.lowerBoundMinimumDistance := 4;
    newcode.upperBoundMinimumDistance := 4;
    newcode.boundsCoveringRadius := [ 2 ];
        
    return( newcode );    
end;

########################################################################
##
#F  EnlargedTombakCode( );
##

EnlargedTombakCode := function ( arg )
    
    local checkmat, m, nmat, els, i, j, w1, w2, w3,
          fieldels, binels, newels, binnewels, size, one, newcode, bmat, binel, el, k, l, beta, gamma, delta, w1, w2, w3, i, u, betabin, gammbin, deltabin, w1bin, w2bin, w3bin;
    
    m := arg[ 1 ];
    
    if m < 6 then
        Error( "EnlargedTombakCode: <m> must be at least 6" );
    fi;
    
    beta := arg[ 3 ];
    gamma := arg[ 4 ];
    
    w1 := arg[ 5 ];
    w2 := arg[ 6 ];
    
    i := arg[ 2 ];

    u := arg[ 7 ];
    
    checkmat := NullMat( 23 * 2^(m-4) - 3, 2 * m - 1, GF( 2 ) );
    
    size := 2^(m-1);
    
    fieldels := SortedGaloisFieldElements( GF( size ) );
    binels := BinaryRepresentation( fieldels, m-1 );
    
    one := Z(2)^0;
    
    # make matrix pi
    
    bmat := TransposedMat( Copy( CheckMat( 
                    TombakCode( m-1, i, beta, gamma, w1, w2, false ) ) ) );

    for k in [ 1 .. Length( bmat ) ] do
        for l in [ 1 .. Length( bmat[ 1 ] ) ] do
            checkmat[ k ][ l + 1 ] := bmat[ k ][ l ];
        od;
    od;
    
    # make matrix J
    
    binel := BinaryRepresentation( u, m-1 );
    
    for k in [ 1 .. size ] do
        checkmat[ Length( bmat ) + k ][ 1 ] := one;
        for l in [ 1 .. m-1 ] do
            checkmat[ Length( bmat ) + k ][ 1 + l ] := binel[ l ];
        od;
        for l in [ 1 .. m-1 ] do
            checkmat[ Length( bmat ) + k ][ m + l ] := binels[ k ][ l ];
        od;
    od;
        
    checkmat := TransposedMat( checkmat );
    
    newcode := CheckMatCode( checkmat, 
                       Concatenation( "enlarged Tombak code (m="
                               ,String(m),")" ),
                       GF(2) );
    
#    newcode.wordLength := 5 * size - 1;
#    newcode.dimension := newcode.wordLength - (2 * m);
    newcode.lowerBoundMinimumDistance := 4;
    newcode.upperBoundMinimumDistance := 4;
    newcode.boundsCoveringRadius := [ 2 ];
        
    return( newcode );    
end;
        
