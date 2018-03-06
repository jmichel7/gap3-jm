########################################################################
##
#F  AmalgatedDirectSum( <C>, <D>, [, <check> ] )
##
##  Return the amalgated direct sum code of C en D.
##  
##  This construction is derived from the direct sum construction,
##  but it saves one coordinate over the direct sum.
##  
##  The amalgated direct sum construction goes as follows:

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

AmalgatedDirectSum := function ( arg )
    
    local C, D, check,
          G, H, Cstandard, Dstandard,
          NewG, NewH, Nulmat, NewC,
          i;
    
    # check the arguments
    if Length( arg ) < 2 or Length( arg ) > 3 then
        Error( "Usage: AmalgatedDirectSum( <C>, <D>, [, <check> ] )" );
    fi;
    
    if Length( arg ) = 2 then
        check := true;
    else
        check := arg[ 3 ];
        if check <> true and check <> false then
            Error( "AmalgatedDirectSum: <check> must be true or false" );
        fi;
    fi;
    
    C := arg[ 1 ];
    D := arg[ 2 ];
    if not IsCode( C ) or not IsCode( D ) then
        Error( "AmalgatedDirectSum: <C> and <D> must be codes" );
    fi;
    
    if not IsLinearCode( C ) or not IsLinearCode( D ) then
        Error( "AmalgatedDirectSum: <C> and <D> must be linear." );
    fi;
    
    if Field( C ) <> Field( D ) then
        Error( "AmalgatedDirectSum: <C> and <D> must be in the same field.");
    fi;
    
    # Copy, otherwise strange things may happen
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
        if not CoordinateAcceptable( Cstandard, WordLength( Cstandard ) ) then
            Error( "AmalgatedDirectSum: Standard form of <C> is not", 
                   "acceptable in the last coordinate.");
        fi;
        # is the last coordinate of the standardcode D' acceptable ?
        if not CoordinateAcceptable( Dstandard, 1 ) then
            Error( "AmalgatedDirectSum: Standard form of <D> is not", 
                   "acceptable in first coordinate.");
        fi;
    fi;
    
    NewG := G;
    # build upper part of the new generator matrix
    # append n(D)-1 zeroes to all rows of G'
    for i in NewG do 
        Append( i, List( [ 1..WordLength( D ) - 1 ], x -> Field(C).zero ) );
    od;
    
    # concatenate the last row of G' and the first row of H'
    for i in [2..WordLength(D)] do
        NewG[Dimension(C)][WordLength(C)+i-1] := H[1][i];
    od;
    
    # throw away the first row of H' (it is already appended)
    NewH := List( [ 2..Length( H ) ], x -> H[ x ] );
    
    # throw away the first column of H' (it is (1, 0, ..., 0), 
    # the one is already present in the new generator matrix)
    NewH := List( [ 1..Length( NewH ) ], 
                  x -> List( [ 2..WordLength( D ) ], y -> NewH[x][y] ) );
    
    # fill the lower left part with zeroes
    Nulmat := List( [ 1..Dimension( D ) - 1 ],
                    x -> NullVector( WordLength( C ) , Field(C) ) );
    
    # and put the remainmants of H' in the lower right corner
    for i in [1..Length(Nulmat)] do
        Append( Nulmat[i], NewH[i] );
    od;
    
    # paste the lower and the upper part
    Append(NewG, Nulmat);
    
    # construct the new code with the new generator matrix
    NewC := GeneratorMatCode( NewG, "amalgated sum code", Field( C ) );
    
    # write history
    NewC.history := MergeHistories( History( C ), History( D ) );
    
    return NewC;
    
end;

########################################################################
##
#F  BlockwiseDirectSumCode( <C1>, <L1>, <C2>, <L2> )
##
##  Return the blockwise direct sum of C1 and C2 with respect to 
##  the cosets defined by the codewords in L1 and L2.
##
##  This function should be generalised, as in some
##  articles from Sloane. (I think)

BlockwiseDirectSumCode := function ( arg )
    
    local C1, C2, L1, L2, 
          k, i,
          subcode1, subcode2, 
          sum, newcode,
          d11, d22, ds1, ds2,
          newels;
    
    # check the arguments
    if Length( arg ) <> 4 then
        Error( "Usage: BlockwiseDirectSumCode( <C1>, <L1>, <C2>, <L2> )" );
    fi;
    
    C1 := arg[ 1 ]; C2 := arg[ 3 ];
    L1 := arg[ 2 ]; L2 := arg[ 4 ];
    
    if not IsCode( C1 ) or not IsCode( C2 ) then
        Error( "BlockwiseDirectSumCode: <C1> and <C2> must be codes" );
    fi;
    
    if Field(C1) <> Field(C2) then
        Error( "BlockwiseDirectSumCode: <C1> and <C2>",
               "must be in the same field" );
    fi;
    
    if Length(L1) <> Length(L2) then
        Error( "BlockwiseDirectSumCode: <L1> and <L2>",
               "must have equal lengths" );
    fi;
    
    # take very large minimum distances
    d11:=10^10000;
    d22:=10^10000;
    
    k := Length(L1);
    
    # make the elements of the new code
    newels := [];
    for i in [1..k] do
        # build subcode 1, using the GUAVA-function CosetCode,
        # which adds the word L1[i] to all codewords of C1
        subcode1 := CosetCode(C1, L1[i]);
        ds1 := MinimumDistance( subcode1 );
        # is the new minimumdistance smaller ?
        if ds1 < d11 then 
            d11 := ds1;
        fi;
        # same for subcode 2
        subcode2 := CosetCode(C2, L2[i]);
        ds2 := MinimumDistance( subcode2 );
        if ds2 < d22 then
            d22 := ds2;
        fi;
        # now build the the direct sum of the two subcodes
        sum := DirectSumCode( subcode1, subcode2 );
        # add the elements of the direct sum-code to
        # the elements we already found 
        # in the previous steps
        newels := Union( newels, Elements(sum) );
    od;
    # finally build the new code with the computed elements
    newcode := ElementsCode( newels, "blockwise direct sum code",
                       Field( C1 ) );
    # compute the lowerbound for the minimumdistance
    newcode.lowerBoundMinimumDistance := Minimum( d11, d22,
                          MinimumDistance(C1) + MinimumDistance(C2) );
    # write history
    newcode.history := MergeHistories( History( C1 ), History( C2 ) );
    return newcode;
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
    local i, G, GL, NewG, j, firstzeros, lastzeros, GB;
    # check the arguments
    # are we dealing with codes ?
    if not IsCode( L ) or not IsCode( B ) then
        Error( "L and B must be codes" );
    fi;
    # are they linear ?
    if not IsLinearCode( L ) or not IsLinearCode( B ) then
        Error( "L and B must be linear code" );
    fi;
    # do they have the same length
    if WordLength( L ) <> WordLength( B ) then
        Error( "L and B must have the same length" );
    fi;
    # m is the number of copies, must be integer >=1
    if not IsInt( m ) or m < 1 then
        Error( "m must be an positive integer" );
    fi;
    
    # get the generator matrices of L en B
    # Copy, so we don't accidentally mess with them
    GL := Copy( GeneratorMat( L ) );
    GB := Copy( GeneratorMat( B ) );
    
    # the new generator matrix, fill with zeros first
    NewG := List( [ 1..Dimension( L ) * ( m + 1) ],
                  x -> NullWord( WordLength( L ) * m ) );
    # first m * Dimension(L) rows,
    # form: [ GL 0  0  0 ... 0  ]
    #       [ 0  GL 0  0 ... 0  ]
    #       [ ................  ]
    #       [ 0  0  0  0     GL ]
    
    for i in [ 1..m ] do
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
    return GeneratorMatCode( NewG, Field( L ) );
end;

