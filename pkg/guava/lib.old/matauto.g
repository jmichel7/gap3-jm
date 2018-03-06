#############################################################################
##
#A  matauto.g                   GUAVA library                  Heiko Thei"sen
##
#H  @(#)$Id: matauto.g,v 1.1 1997/01/20 15:10:50 werner Exp $
##
##  This file  contains    functions  that perform  backtrack   searches  for
##  automorphisms resp. isomorphisms of matrices over GF(2).
##
#H  $Log: matauto.g,v $
#H  Revision 1.1  1997/01/20 15:10:50  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
##

#############################################################################
##

#V  BmatrixOps  . . . . . . . . . . .  operations record for boolean matrices
##
BmatrixOps := OperationsRecord( "BmatrixOps", GroupElementOps );

#############################################################################
##
#F  BmatrixMatrix( <A> )  . . . . . . . . . . . . . convert matrix to bmatrix
##
BmatrixMatrix := function( A )
    local  B,  i,  row,  bit;
    
    B := rec( rows := Length( A ),
           columns := Length( A[ 1 ] ),
        operations := BmatrixOps );
    B.blist := BlistList( [ 1 .. B.rows * B.columns ], [  ] );
    i := 0;
    for row  in A  do
        for bit  in row  do
            i := i + 1;
            if bit <> 0 * bit  then
                B.blist[ i ] := true;
            fi;
        od;
    od;
    return B;
end;

#############################################################################
##
#F  MatrixBmatrix( <A> )  . . . . . . . . . . . . . convert bmatrix to matrix
##
MatrixBmatrix := function( A )
    local  B,  one,  pos;
    
    B := NullMat( A.rows, A.column, GF(2) );
    one := Z(2);
    pos := Position( A.blist, true );
    while pos <> false  do
        B[ QuoInt( pos - 1, A.rows ) + 1 ][ ( pos - 1 ) mod A.rows + 1 ]
          := one;
        pos := Position( A.blist, true, pos );
    od;
    return B;
end;
        
BmatrixOps.( "=" ) := function( A, B )
    return     A.rows    = B.rows
           and A.columns = B.columns
           and A.blist   = B.blist;
end;

BmatrixOps.( "<" ) := function( A, B )
    return   A.rows    < B.rows or
           ( A.rows    = B.rows and
           ( A.columns < B.columns or
           ( A.columns = B.columns and
             A.blist   < B.blist ) ) );
end;

BmatrixOps.( "*" ) := function( A, B )
    local  P,  col,  row,  i,  r,  R,  c,  C;
    
    P := rec( rows := A.rows,
           columns := B.columns,
             blist := BlistList( [ 1 .. A.rows * B.columns ], [  ] ),
        operations := BmatrixOps );
    col := [ 0, B.rows .. ( B.rows - 1 ) * B.columns ];
    row := [ -A.rows + 1 .. 0 ];
    i := 0;
    for r  in [ 1 .. A.rows ]  do
        R := A.blist{ row + r * A.rows };
        for c  in [ 1 .. B.columns ]  do
            C := B.blist{ col + c };
            i := i + 1;
            IntersectBlist( C, R );
            if SizeBlist( C ) mod 2 = 1  then
                P[ i ] := true;
            fi;
        od;
    od;
    return P;
end;

BmatrixOps.( "^" ) := function( A, B )
    if not IsInt( B )  then
        B := MatrixBmatrix( B );
    fi;
    return BmatrixMatrix( MatrixBmatrix( A ) ^ B );
end;

BmatrixOps.Print := function( A )
    local  i,  row,  bit;
    
    i := 0;
    for row  in [ 1 .. A.rows ]  do
        for bit  in [ 1 .. A.columns ]  do
            i := i + 1 ;
            if A.blist[ i ]  then
                Print( "1" );
            else
                Print( "." );
            fi;
        od;
        Print( "\n" );
    od;
    Print( "\n" );
end;

#############################################################################
##
#F  ConjugateBmatrix( <A>, <g> )  . . . . . . . . . .  conjugate of a bmatrix
##
ConjugateBmatrix := function( A, g )
    local  R,  r,  i;
    
    R := rec( rows := A.rows,
           columns := A.columns,
             blist := BlistList( [ 1 .. Length( A.blist ) ], [  ] ),
        operations := A.operations );
    for r  in [ 1 .. A.rows ]  do
        i := ( r + A.columns ) ^ g - A.columns;
        R.blist{ [ 1 .. A.columns ] + ( i - 1 ) * A.columns } :=
          Permuted( A.blist
                  { [ 1 .. A.columns ] + ( r - 1 ) * A.columns }, g );
    od;
    return R;
end;

#############################################################################
##
#F  BmatrixPartition( <A>, <CR>, <i> )  . . . . . . . . . . . . . . . . local
##
BmatrixPartition := function( A, CR, i )
    local  rows,  P,  Q,  c,  num,  p;
    
    rows := Cell( CR, i );
    P := [  ];
    if rows[ 1 ] > A.columns  then
        rows := A.columns * ( rows - A.columns - 1 );
        for c  in [ 1 .. A.columns ]  do
            num := SizeBlist( A.blist{ c + rows } ) + 1;
            if not IsBound( P[ num ] )  then
                P[ num ] := [  ];
            fi;
            Add( P[ num ], c );
        od;
        Q := [  ];
        for p  in P  do
            Add( Q, p );
        od;
        Add( Q, [ A.columns + 1 .. A.columns + A.rows ] );
    else
        for c  in [ 1 .. A.rows ]  do
            num := SizeBlist( A.blist{ (c-1) * A.columns + rows } ) + 1;
            if not IsBound( P[ num ] )  then
                P[ num ] := [  ];
            fi;
            Add( P[ num ], c );
        od;
        Q := [ [ 1 .. A.columns ] ];
        for p  in P  do
            Add( Q, p + A.columns );
        od;
    fi;
    return Partition( Q );
end;

#############################################################################
##
#F  Refinements.Bmatrix( <P>, <rbase>, <Rf> ) . . . .  refinement for bmatrix
##
Refinements.Bmatrix := function( P, rbase, Rf )
    local  R;
    
    R := BmatrixPartition( Rf[ 3 ], P, Rf[ 4 ] );
    return MeetPartitionStrat( P, R, Rf[ 2 ] );
end;

#############################################################################
##
#F  RepOpMatricesCodes( <repr>, <A>, <B>, <Pr> )   repop on matrices or codes
##
RepOpMatricesCodes := function( repr, A, B, Pr )
    local  res,        # representative or automorphism group: the result
           G,          # the full symmetric group on rows and columns
           Omega,      # operation domain of <G>
           n,          # maximum of <Omega>
           rbase,      # the R-base for the backtrack algorithm
           K,          # loop variable along stabilizer chain of <G>
           num,        # number of cells before refinement
           R,          # partition to meet with <P>
           P,          # partition refined during construction of <rbase>
           cell,       # cell from which next R-base point is taken
           strat,      # strategy for meeting partitions
           i;          # loop variable
    
    G := SymmetricGroup( A.columns + A.rows );
    Omega := [ 1 .. A.columns + A.rows ];
    P := Partition( [ [             1 .. A.columns ],
                      [ A.columns + 1 .. A.columns + A.rows ] ] );
    
    # Construct an R-base.
    rbase := EmptyRBase( P, P );
    K := InitializeK( G, rbase );
    while K <> false  do
        cell := NextRBasePoint( G, K, Omega, rbase, P );
        if cell <> false  then
            repeat
                num := NumberCells( P );
                for i  in [ 1 .. num ]  do
                    R := BmatrixPartition( A, P, i );
                    strat := StratMeetPartition( P, R );
                    if strat <> false  then
                        Add( rbase.RRfm[ Length( rbase.RRfm ) ],
                             [ "Bmatrix", strat, B, i ] );
                    fi;
                od;
            until NumberCells( P ) = num;
        fi;
        K := StepDownK( G, K, rbase );
    od;
    Add( rbase.RRfd, FingerprintPartition( P ) );
    
    res := PartitionBacktrack( G,
                               Omega,
                               Pr,
                               repr,
                               rbase );
    if res <> false  then
        if repr  then
            return RestrictedPerm( res, [ 1 .. A.columns ] );
        else
            return Operation( res, [ 1 .. A.columns ] );
        fi;
    fi;
end;

#############################################################################
##
#F  RepOpMatrices( <arg> )  . . . . . . . . . . . . . . . . repop on matrices
##
RepOpMatrices := function( arg )
    local  A,  B,  repr;
    
    A := BmatrixMatrix( arg[ 1 ] );
    if Length( arg ) > 1  then
        B := BmatrixMatrix( arg[ 2 ] );
        repr := true;
    else
        B := A;
        repr := false;
    fi;
    return RepOpMatricesCodes( repr, A, B,
                   g -> ConjugateBmatrix( A, g ) = B );
end;

#############################################################################
##
#F  AutomorphismGroupMatrix( <M> )  . . . . .  automorphism group of a matrix
##
AutomorphismGroupMatrix := RepOpMatrices;

#############################################################################
##
#F  PrBinaryLinearCodes( <c1>, <g>, <c2> )  . . . . . . . . . test if c1^g=c2
##
PrBinaryLinearCodes := function( c1, g, c2 )
    local  m1,  m2;
    
    m1 := PermutedCols( c1.generatorMat, g );
    m2 := c2.generatorMat;
    return BaseMat( m1 ) = BaseMat( m2 );
end;
    
#############################################################################
##
#F  RepOpBinaryLinearCodes( <arg> ) . . . . . .  repop on binary linear codes
##
RepOpBinaryLinearCodes := function( arg )
    local  c1,  A,  c2,  B,  repr;
    
    c1 := arg[ 1 ];
    A := BmatrixMatrix( List( MinimumWeightWords( c1 ), v -> v.vector ) );
    if Length( arg ) > 1  then
        c2 := arg[ 2 ];
        B := BmatrixMatrix( List( MinimumWeightWords( c2 ), v->v.vector ) );
        repr := true;
    else
        c2 := c1;
        B := A;
        repr := false;
    fi;
    return RepOpMatricesCodes( repr, A, B,
                   g -> PrBinaryLinearCodes( c1, g, c2 ) );
end;

#############################################################################
##
#F  AutomorphismGroupBinaryLinearCode( <C> )  .  automorphism group of a code
##
AutomorphismGroupBinaryLinearCode := RepOpBinaryLinearCodes;

#############################################################################
##

#E  Emacs variables . . . . . . . . . . . . . . local variables for this file
##  Local Variables:
##  mode:             outline-minor
##  outline-regexp:   "#[AEFTV]"
##  fill-column:      77
##  End:
#############################################################################

