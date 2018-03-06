#############################################################################
##
#A  combinat.tst                GAP tests                    Martin Schoenert
##
#A  @(#)$Id: combinat.tst,v 1.1.1.1 1996/12/11 12:44:37 werner Exp $
##
#Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  tests  the functions that  mainly  deal  with  combinatorics.
##
#H  $Log: combinat.tst,v $
#H  Revision 1.1.1.1  1996/12/11 12:44:37  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.2  1992/11/25  14:49:10  martin
#H  changed the format of the message
#H
#H  Revision 3.1  1992/04/29  08:53:40  martin
#H  changed the output (because of the new printer)
#H
#H  Revision 3.0  1991/07/22  18:13:03  martin
#H  initial revision under RCS
#H
##
SizeScreen([80,24]);
#>[ 80, 24 ]
#F  Factorial( <n> )  . . . . . . . . . . . . . . . . factorial of an integer
List( [0..10], Factorial );
#>[ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 ]
Factorial( 50 );
#>30414093201713378043612608166064768844377641568960512000000000000

#F  Binomial( <n>, <k> )  . . . . . . . . .  binomial coefficient of integers
List( [-8..8], k -> Binomial( 0, k ) );
#>[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ]
List( [-8..8], n -> Binomial( n, 0 ) );
#>[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
ForAll( [-8..8], n -> ForAll( [-2..8], k ->
    Binomial(n,k) = Binomial(n-1,k) + Binomial(n-1,k-1) ) );
#>true
Binomial( 400, 50 );
#>17035900270730601418919867558071677342938596450600561760371485120

#F  Bell( <n> ) . . . . . . . . . . . . . . . . .  value of the Bell sequence
List( [0..10], n -> Bell(n) );
#>[ 1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 ]
List( [0..10], n -> Sum( [0..n], k -> Stirling2( n, k ) ) );
#>[ 1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 ]
Bell( 60 );
#>976939307467007552986994066961675455550246347757474482558637

#F  Stirling1( <n>, <k> ) . . . . . . . . . Stirling number of the first kind
List( [-8..8], k -> Stirling1( 0, k ) );
#>[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ]
List( [-8..8], n -> Stirling1( n, 0 ) );
#>[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ]
ForAll( [-8..8], n -> ForAll( [-8..8], k ->
    Stirling1(n,k) = (n-1) * Stirling1(n-1,k) + Stirling1(n-1,k-1) ) );
#>true
Stirling1( 60, 20 );
#>568611292461582075463109862277030309493811818619783570055397018154658816

#F  Stirling2( <n>, <k> ) . . . . . . . .  Stirling number of the second kind
List( [-8..8], k -> Stirling2( 0, k ) );
#>[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ]
List( [-8..8], n -> Stirling2( n, 0 ) );
#>[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ]
ForAll( [-8..8], n -> ForAll( [-8..8], k ->
    Stirling2(n,k) = k * Stirling2(n-1,k) + Stirling2(n-1,k-1) ) );
#>true
Stirling2( 60, 20 );
#>170886257768137628374668205554120607567311094075812403938286

#F  Combinations( <mset>, <k> ) . . . .  set of sorted sublists of a multiset
Combinations( [] );
#>[ [  ] ]
List( [0..1], k -> Combinations( [], k ) );
#>[ [ [  ] ], [  ] ]
Combinations( [1..4] );
#>[ [  ], [ 1 ], [ 1, 2 ], [ 1, 2, 3 ], [ 1, 2, 3, 4 ], [ 1, 2, 4 ], [ 1, 3 ], 
#>  [ 1, 3, 4 ], [ 1, 4 ], [ 2 ], [ 2, 3 ], [ 2, 3, 4 ], [ 2, 4 ], [ 3 ], 
#>  [ 3, 4 ], [ 4 ] ]
List( [0..5], k -> Combinations( [1..4], k ) );
#>[ [ [  ] ], [ [ 1 ], [ 2 ], [ 3 ], [ 4 ] ], 
#>  [ [ 1, 2 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ], [ 3, 4 ] ], 
#>  [ [ 1, 2, 3 ], [ 1, 2, 4 ], [ 1, 3, 4 ], [ 2, 3, 4 ] ], [ [ 1, 2, 3, 4 ] ], 
#>  [  ] ]
Combinations( [1,2,2,3] );
#>[ [  ], [ 1 ], [ 1, 2 ], [ 1, 2, 2 ], [ 1, 2, 2, 3 ], [ 1, 2, 3 ], [ 1, 3 ], 
#>  [ 2 ], [ 2, 2 ], [ 2, 2, 3 ], [ 2, 3 ], [ 3 ] ]
List( [0..5], k -> Combinations( [1,2,2,3], k ) );
#>[ [ [  ] ], [ [ 1 ], [ 2 ], [ 3 ] ], 
#>  [ [ 1, 2 ], [ 1, 3 ], [ 2, 2 ], [ 2, 3 ] ], 
#>  [ [ 1, 2, 2 ], [ 1, 2, 3 ], [ 2, 2, 3 ] ], [ [ 1, 2, 2, 3 ] ], [  ] ]
Combinations( [1..12] )[4039];
#>[ 7, 8, 9, 10, 11, 12 ]
Combinations( [1..16], 4 )[266];
#>[ 1, 5, 9, 13 ]
Combinations( [1,2,3,3,4,4,5,5,5,6,6,6,7,7,7,7] )[378];
#>[ 1, 2, 3, 4, 5, 6, 7 ]
Combinations( [1,2,3,3,4,4,5,5,5,6,6,6,7,7,7,7,8,8,8,8], 8 )[97];
#>[ 1, 2, 3, 4, 5, 6, 7, 8 ]

#F  NrCombinations( <mset>, <k> ) . . number of sorted sublists of a multiset
NrCombinations( [] );
#>1
List( [0..1], k -> NrCombinations( [], k ) );
#>[ 1, 0 ]
NrCombinations( [1..4] );
#>16
List( [0..5], k -> NrCombinations( [1..4], k ) );
#>[ 1, 4, 6, 4, 1, 0 ]
NrCombinations( [1,2,2,3] );
#>12
List( [0..5], k -> NrCombinations( [1,2,2,3], k ) );
#>[ 1, 3, 4, 3, 1, 0 ]
NrCombinations( [1..12] );
#>4096
NrCombinations( [1..16], 4 );
#>1820
NrCombinations( [1,2,3,3,4,4,5,5,5,6,6,6,7,7,7,7] );
#>2880
NrCombinations( [1,2,3,3,4,4,5,5,5,6,6,6,7,7,7,7,8,8,8,8], 8 );
#>1558

#F  Arrangements( <mset> )  . . . . set of ordered combinations of a multiset
Arrangements( [] );
#>[ [  ] ]
List( [0..1], k -> Arrangements( [], k ) );
#>[ [ [  ] ], [  ] ]
Arrangements( [1..3] );
#>[ [  ], [ 1 ], [ 1, 2 ], [ 1, 2, 3 ], [ 1, 3 ], [ 1, 3, 2 ], [ 2 ], [ 2, 1 ], 
#>  [ 2, 1, 3 ], [ 2, 3 ], [ 2, 3, 1 ], [ 3 ], [ 3, 1 ], [ 3, 1, 2 ], [ 3, 2 ], 
#>  [ 3, 2, 1 ] ]
List( [0..4], k -> Arrangements( [1..3], k ) );
#>[ [ [  ] ], [ [ 1 ], [ 2 ], [ 3 ] ], 
#>  [ [ 1, 2 ], [ 1, 3 ], [ 2, 1 ], [ 2, 3 ], [ 3, 1 ], [ 3, 2 ] ], 
#>  [ [ 1, 2, 3 ], [ 1, 3, 2 ], [ 2, 1, 3 ], [ 2, 3, 1 ], [ 3, 1, 2 ], 
#>      [ 3, 2, 1 ] ], [  ] ]
Arrangements( [1,2,2,3] );
#>[ [  ], [ 1 ], [ 1, 2 ], [ 1, 2, 2 ], [ 1, 2, 2, 3 ], [ 1, 2, 3 ], 
#>  [ 1, 2, 3, 2 ], [ 1, 3 ], [ 1, 3, 2 ], [ 1, 3, 2, 2 ], [ 2 ], [ 2, 1 ], 
#>  [ 2, 1, 2 ], [ 2, 1, 2, 3 ], [ 2, 1, 3 ], [ 2, 1, 3, 2 ], [ 2, 2 ], 
#>  [ 2, 2, 1 ], [ 2, 2, 1, 3 ], [ 2, 2, 3 ], [ 2, 2, 3, 1 ], [ 2, 3 ], 
#>  [ 2, 3, 1 ], [ 2, 3, 1, 2 ], [ 2, 3, 2 ], [ 2, 3, 2, 1 ], [ 3 ], [ 3, 1 ], 
#>  [ 3, 1, 2 ], [ 3, 1, 2, 2 ], [ 3, 2 ], [ 3, 2, 1 ], [ 3, 2, 1, 2 ], 
#>  [ 3, 2, 2 ], [ 3, 2, 2, 1 ] ]
List( [0..5], k -> Arrangements( [1,2,2,3], k ) );
#>[ [ [  ] ], [ [ 1 ], [ 2 ], [ 3 ] ], 
#>  [ [ 1, 2 ], [ 1, 3 ], [ 2, 1 ], [ 2, 2 ], [ 2, 3 ], [ 3, 1 ], [ 3, 2 ] ], 
#>  [ [ 1, 2, 2 ], [ 1, 2, 3 ], [ 1, 3, 2 ], [ 2, 1, 2 ], [ 2, 1, 3 ], 
#>      [ 2, 2, 1 ], [ 2, 2, 3 ], [ 2, 3, 1 ], [ 2, 3, 2 ], [ 3, 1, 2 ], 
#>      [ 3, 2, 1 ], [ 3, 2, 2 ] ], 
#>  [ [ 1, 2, 2, 3 ], [ 1, 2, 3, 2 ], [ 1, 3, 2, 2 ], [ 2, 1, 2, 3 ], 
#>      [ 2, 1, 3, 2 ], [ 2, 2, 1, 3 ], [ 2, 2, 3, 1 ], [ 2, 3, 1, 2 ], 
#>      [ 2, 3, 2, 1 ], [ 3, 1, 2, 2 ], [ 3, 2, 1, 2 ], [ 3, 2, 2, 1 ] ], [  ] ]
Arrangements( [1..6] )[736];
#>[ 3, 2, 1, 6, 5, 4 ]
Arrangements( [1..8], 4 )[443];
#>[ 3, 1, 7, 5 ]
Arrangements( [1,2,3,3,4,4,5] )[3511];
#>[ 5, 4, 3, 2, 1 ]
Arrangements( [1,2,3,4,4,5,5,6,6], 5 )[424];
#>[ 2, 3, 4, 5, 6 ]

#F  NrArrangements( <mset>, <k> ) . . number of sorted sublists of a multiset
NrArrangements( [] );
#>1
List( [0..1], k -> NrArrangements( [], k ) );
#>[ 1, 0 ]
NrArrangements( [1..3] );
#>16
List( [0..4], k -> NrArrangements( [1..3], k ) );
#>[ 1, 3, 6, 6, 0 ]
NrArrangements( [1,2,2,3] );
#>35
List( [0..5], k -> NrArrangements( [1,2,2,3], k ) );
#>[ 1, 3, 7, 12, 12, 0 ]
NrArrangements( [1..6] );
#>1957
NrArrangements( [1..8], 4 );
#>1680
NrArrangements( [1,2,3,3,4,4,5] );
#>3592
NrArrangements( [1,2,3,4,4,5,5,6,6], 5 );
#>2880

#F  UnorderedTuples( <set>, <k> ) . . . .  set of unordered tuples from a set
List( [0..1], k -> UnorderedTuples( [], k ) );
#>[ [ [  ] ], [  ] ]
List( [0..4], k -> UnorderedTuples( [1..3], k ) );
#>[ [ [  ] ], [ [ 1 ], [ 2 ], [ 3 ] ], 
#>  [ [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 2 ], [ 2, 3 ], [ 3, 3 ] ], 
#>  [ [ 1, 1, 1 ], [ 1, 1, 2 ], [ 1, 1, 3 ], [ 1, 2, 2 ], [ 1, 2, 3 ], 
#>      [ 1, 3, 3 ], [ 2, 2, 2 ], [ 2, 2, 3 ], [ 2, 3, 3 ], [ 3, 3, 3 ] ], 
#>  [ [ 1, 1, 1, 1 ], [ 1, 1, 1, 2 ], [ 1, 1, 1, 3 ], [ 1, 1, 2, 2 ], 
#>      [ 1, 1, 2, 3 ], [ 1, 1, 3, 3 ], [ 1, 2, 2, 2 ], [ 1, 2, 2, 3 ], 
#>      [ 1, 2, 3, 3 ], [ 1, 3, 3, 3 ], [ 2, 2, 2, 2 ], [ 2, 2, 2, 3 ], 
#>      [ 2, 2, 3, 3 ], [ 2, 3, 3, 3 ], [ 3, 3, 3, 3 ] ] ]
UnorderedTuples( [1..10], 6 )[1459];
#>[ 1, 3, 5, 7, 9, 10 ]

#F  NrUnorderedTuples( <set>, <k> ) . . number unordered of tuples from a set
List( [0..1], k -> NrUnorderedTuples( [], k ) );
#>[ 1, 0 ]
List( [0..4], k -> NrUnorderedTuples( [1..3], k ) );
#>[ 1, 3, 6, 10, 15 ]
NrUnorderedTuples( [1..10], 6 );
#>5005

#F  Tuples( <set>, <k> )  . . . . . . . . .  set of ordered tuples from a set
List( [0..1], k -> Tuples( [], k ) );
#>[ [ [  ] ], [  ] ]
List( [0..3], k -> Tuples( [1..3], k ) );
#>[ [ [  ] ], [ [ 1 ], [ 2 ], [ 3 ] ], 
#>  [ [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 1 ], [ 2, 2 ], [ 2, 3 ], [ 3, 1 ], 
#>      [ 3, 2 ], [ 3, 3 ] ], 
#>  [ [ 1, 1, 1 ], [ 1, 1, 2 ], [ 1, 1, 3 ], [ 1, 2, 1 ], [ 1, 2, 2 ], 
#>      [ 1, 2, 3 ], [ 1, 3, 1 ], [ 1, 3, 2 ], [ 1, 3, 3 ], [ 2, 1, 1 ], 
#>      [ 2, 1, 2 ], [ 2, 1, 3 ], [ 2, 2, 1 ], [ 2, 2, 2 ], [ 2, 2, 3 ], 
#>      [ 2, 3, 1 ], [ 2, 3, 2 ], [ 2, 3, 3 ], [ 3, 1, 1 ], [ 3, 1, 2 ], 
#>      [ 3, 1, 3 ], [ 3, 2, 1 ], [ 3, 2, 2 ], [ 3, 2, 3 ], [ 3, 3, 1 ], 
#>      [ 3, 3, 2 ], [ 3, 3, 3 ] ] ]
Tuples( [1..8], 4 )[167];
#>[ 1, 3, 5, 7 ]

#F  NrTuples( <set>, <k> )  . . . . . . . number of ordered tuples from a set
List( [0..1], k -> NrTuples( [], k ) );
#>[ 1, 0 ]
List( [0..3], k -> NrTuples( [1..3], k ) );
#>[ 1, 3, 9, 27 ]
NrTuples( [1..8], 4 );
#>4096

#F  PermutationsList( <mset> )  . . . . . . set of permutations of a multiset
PermutationsList( [] );
#>[ [  ] ]
PermutationsList( [1..4] );
#>[ [ 1, 2, 3, 4 ], [ 1, 2, 4, 3 ], [ 1, 3, 2, 4 ], [ 1, 3, 4, 2 ], 
#>  [ 1, 4, 2, 3 ], [ 1, 4, 3, 2 ], [ 2, 1, 3, 4 ], [ 2, 1, 4, 3 ], 
#>  [ 2, 3, 1, 4 ], [ 2, 3, 4, 1 ], [ 2, 4, 1, 3 ], [ 2, 4, 3, 1 ], 
#>  [ 3, 1, 2, 4 ], [ 3, 1, 4, 2 ], [ 3, 2, 1, 4 ], [ 3, 2, 4, 1 ], 
#>  [ 3, 4, 1, 2 ], [ 3, 4, 2, 1 ], [ 4, 1, 2, 3 ], [ 4, 1, 3, 2 ], 
#>  [ 4, 2, 1, 3 ], [ 4, 2, 3, 1 ], [ 4, 3, 1, 2 ], [ 4, 3, 2, 1 ] ]
PermutationsList( [1,2,2,3,] );
#>[ [ 1, 2, 2, 3 ], [ 1, 2, 3, 2 ], [ 1, 3, 2, 2 ], [ 2, 1, 2, 3 ], 
#>  [ 2, 1, 3, 2 ], [ 2, 2, 1, 3 ], [ 2, 2, 3, 1 ], [ 2, 3, 1, 2 ], 
#>  [ 2, 3, 2, 1 ], [ 3, 1, 2, 2 ], [ 3, 2, 1, 2 ], [ 3, 2, 2, 1 ] ]
PermutationsList( [1..6] )[ 128 ];
#>[ 2, 1, 4, 3, 6, 5 ]
PermutationsList( [1,2,2,3,3,4,4,4] )[1359];
#>[ 4, 3, 2, 1, 4, 3, 2, 4 ]

#F  NrPermutationsList( <mset> )  . . .  number of permutations of a multiset
NrPermutationsList( [] );
#>1
NrPermutationsList( [1..4] );
#>24
NrPermutationsList( [1,2,2,3] );
#>12
NrPermutationsList( [1..6] );
#>720
NrPermutationsList( [1,2,2,3,3,4,4,4] );
#>1680

#F  Derangements( <list> ) . . . . set of fixpointfree permutations of a list
Derangements( [] );
#>[ [  ] ]
Derangements( [1..4] );
#>[ [ 2, 1, 4, 3 ], [ 2, 3, 4, 1 ], [ 2, 4, 1, 3 ], [ 3, 1, 4, 2 ], 
#>  [ 3, 4, 1, 2 ], [ 3, 4, 2, 1 ], [ 4, 1, 2, 3 ], [ 4, 3, 1, 2 ], 
#>  [ 4, 3, 2, 1 ] ]
Derangements( [1..6] )[ 128 ];
#>[ 4, 3, 6, 1, 2, 5 ]
Derangements( [1,2,2,3,3,4,4,4] )[64];
#>[ 4, 1, 4, 2, 4, 2, 3, 3 ]

#F  NrDerangements( <list> ) .  number of fixpointfree permutations of a list
NrDerangements( [] );
#>1
NrDerangements( [1..4] );
#>9
NrDerangements( [1..6] );
#>265
NrDerangements( [1,2,2,3,3,4,4,4] );
#>126

#F  Permanent( <mat> )  . . . . . . . . . . . . . . . . permanent of a matrix
Permanent( [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]] );
#>9
Permanent( [[1,1,0,1,0,0,0],[0,1,1,0,1,0,0],[0,0,1,1,0,1,0],[0,0,0,1,1,0,1],
            [1,0,0,0,1,1,0],[0,1,0,0,0,1,1],[1,0,1,0,0,0,1]] );
#>24

#F  PartitionsSet( <set> )  . . . . . . . . . . .  set of partitions of a set
PartitionsSet( [] );
#>[ [  ] ]
List( [0..1], k -> PartitionsSet( [], k ) );
#>[ [ [  ] ], [  ] ]
PartitionsSet( [1..4] );
#>[ [ [ 1 ], [ 2 ], [ 3 ], [ 4 ] ], [ [ 1 ], [ 2 ], [ 3, 4 ] ], 
#>  [ [ 1 ], [ 2, 3 ], [ 4 ] ], [ [ 1 ], [ 2, 3, 4 ] ], 
#>  [ [ 1 ], [ 2, 4 ], [ 3 ] ], [ [ 1, 2 ], [ 3 ], [ 4 ] ], 
#>  [ [ 1, 2 ], [ 3, 4 ] ], [ [ 1, 2, 3 ], [ 4 ] ], [ [ 1, 2, 3, 4 ] ], 
#>  [ [ 1, 2, 4 ], [ 3 ] ], [ [ 1, 3 ], [ 2 ], [ 4 ] ], [ [ 1, 3 ], [ 2, 4 ] ], 
#>  [ [ 1, 3, 4 ], [ 2 ] ], [ [ 1, 4 ], [ 2 ], [ 3 ] ], [ [ 1, 4 ], [ 2, 3 ] ] ]
List( [0..4], k -> PartitionsSet( [1..3], k ) );
#>[ [  ], [ [ [ 1, 2, 3 ] ] ], 
#>  [ [ [ 1 ], [ 2, 3 ] ], [ [ 1, 2 ], [ 3 ] ], [ [ 1, 3 ], [ 2 ] ] ], 
#>  [ [ [ 1 ], [ 2 ], [ 3 ] ] ], [  ] ]
PartitionsSet( [1..7] )[521];
#>[ [ 1, 3, 5, 7 ], [ 2, 4, 6 ] ]
PartitionsSet( [1..8], 3 )[96];
#>[ [ 1, 2, 3 ], [ 4, 5 ], [ 6, 7, 8 ] ]

#F  NrPartitionsSet( <set> )  . . . . . . . . . number of partitions of a set
NrPartitionsSet( [] );
#>1
List( [0..1], k -> NrPartitionsSet( [], k ) );
#>[ 1, 0 ]
NrPartitionsSet( [1..4] );
#>15
List( [0..4], k -> NrPartitionsSet( [1,2,3], k ) );
#>[ 0, 1, 3, 1, 0 ]
NrPartitionsSet( [1..8] );
#>4140
NrPartitionsSet( [1..9], 3 );
#>3025

#F  Partitions( <n> ) . . . . . . . . . . . . set of partitions of an integer
Partitions( 0 );
#>[ [  ] ]
List( [0..1], k -> Partitions( 0, k ) );
#>[ [ [  ] ], [  ] ]
Partitions( 6 );
#>[ [ 1, 1, 1, 1, 1, 1 ], [ 2, 1, 1, 1, 1 ], [ 2, 2, 1, 1 ], [ 2, 2, 2 ], 
#>  [ 3, 1, 1, 1 ], [ 3, 2, 1 ], [ 3, 3 ], [ 4, 1, 1 ], [ 4, 2 ], [ 5, 1 ], 
#>  [ 6 ] ]
List( [0..7], k -> Partitions( 6, k ) );
#>[ [  ], [ [ 6 ] ], [ [ 3, 3 ], [ 4, 2 ], [ 5, 1 ] ], 
#>  [ [ 2, 2, 2 ], [ 3, 2, 1 ], [ 4, 1, 1 ] ], 
#>  [ [ 2, 2, 1, 1 ], [ 3, 1, 1, 1 ] ], [ [ 2, 1, 1, 1, 1 ] ], 
#>  [ [ 1, 1, 1, 1, 1, 1 ] ], [  ] ]
Partitions( 20 )[314];
#>[ 7, 4, 3, 3, 2, 1 ]
Partitions( 20, 10 )[17];
#>[ 5, 3, 3, 2, 2, 1, 1, 1, 1, 1 ]

#F  NrPartitions( <n> ) . . . . . . . . .  number of partitions of an integer
NrPartitions( 0 );
#>1
List( [0..1], k -> NrPartitions( 0, k ) );
#>[ 1, 0 ]
NrPartitions( 6 );
#>11
List( [0..7], k -> NrPartitions( 6, k ) );
#>[ 0, 1, 3, 3, 2, 1, 1, 0 ]
NrPartitions( 100 );
#>190569292
NrPartitions( 100, 10 );
#>2977866

#F  OrderedPartitions( <n> ) . . . .  set of ordered partitions of an integer
OrderedPartitions( 0 );
#>[ [  ] ]
List( [0..1], k -> OrderedPartitions( 0, k ) );
#>[ [ [  ] ], [  ] ]
OrderedPartitions( 5 );
#>[ [ 1, 1, 1, 1, 1 ], [ 1, 1, 1, 2 ], [ 1, 1, 2, 1 ], [ 1, 1, 3 ], 
#>  [ 1, 2, 1, 1 ], [ 1, 2, 2 ], [ 1, 3, 1 ], [ 1, 4 ], [ 2, 1, 1, 1 ], 
#>  [ 2, 1, 2 ], [ 2, 2, 1 ], [ 2, 3 ], [ 3, 1, 1 ], [ 3, 2 ], [ 4, 1 ], [ 5 ] ]
List( [0..6], k -> OrderedPartitions( 5, k ) );
#>[ [  ], [ [ 5 ] ], [ [ 1, 4 ], [ 2, 3 ], [ 3, 2 ], [ 4, 1 ] ], 
#>  [ [ 1, 1, 3 ], [ 1, 2, 2 ], [ 1, 3, 1 ], [ 2, 1, 2 ], [ 2, 2, 1 ], 
#>      [ 3, 1, 1 ] ], 
#>  [ [ 1, 1, 1, 2 ], [ 1, 1, 2, 1 ], [ 1, 2, 1, 1 ], [ 2, 1, 1, 1 ] ], 
#>  [ [ 1, 1, 1, 1, 1 ] ], [  ] ]
OrderedPartitions( 13 )[2048];
#>[ 1, 12 ]
OrderedPartitions( 16, 6 )[1001];
#>[ 1, 11, 1, 1, 1, 1 ]

#F  NrOrderedPartitions( <n> ) . . number of ordered partitions of an integer
NrOrderedPartitions( 0 );
#>1
List( [0..1], k -> NrOrderedPartitions( 0, k ) );
#>[ 1, 0 ]
NrOrderedPartitions( 5 );
#>16
List( [0..6], k -> NrOrderedPartitions( 5, k ) );
#>[ 0, 1, 4, 6, 4, 1, 0 ]
NrOrderedPartitions( 13 );
#>4096
NrOrderedPartitions( 16, 6 );
#>3003

#F  RestrictedPartitions( <n>, <set> )  . restricted partitions of an integer
RestrictedPartitions( 0, [1..10] );
#>[ [  ] ]
List( [0..1], k -> RestrictedPartitions( 0, [1..10], k ) );
#>[ [ [  ] ], [  ] ]
RestrictedPartitions( 10, [1,2,5,10] );
#>[ [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ], [ 2, 1, 1, 1, 1, 1, 1, 1, 1 ], 
#>  [ 2, 2, 1, 1, 1, 1, 1, 1 ], [ 2, 2, 2, 1, 1, 1, 1 ], [ 2, 2, 2, 2, 1, 1 ], 
#>  [ 2, 2, 2, 2, 2 ], [ 5, 1, 1, 1, 1, 1 ], [ 5, 2, 1, 1, 1 ], [ 5, 2, 2, 1 ], 
#>  [ 5, 5 ], [ 10 ] ]
List( [1..10], k -> RestrictedPartitions( 10, [1,2,5,10], k ) );
#>[ [ [ 10 ] ], [ [ 5, 5 ] ], [  ], [ [ 5, 2, 2, 1 ] ], 
#>  [ [ 2, 2, 2, 2, 2 ], [ 5, 2, 1, 1, 1 ] ], 
#>  [ [ 2, 2, 2, 2, 1, 1 ], [ 5, 1, 1, 1, 1, 1 ] ], [ [ 2, 2, 2, 1, 1, 1, 1 ] ],
#>  [ [ 2, 2, 1, 1, 1, 1, 1, 1 ] ], [ [ 2, 1, 1, 1, 1, 1, 1, 1, 1 ] ], 
#>  [ [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ] ] ]
RestrictedPartitions( 20, [2,5,10] );
#>[ [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ], [ 5, 5, 2, 2, 2, 2, 2 ], [ 5, 5, 5, 5 ], 
#>  [ 10, 2, 2, 2, 2, 2 ], [ 10, 5, 5 ], [ 10, 10 ] ]
List( [1..20], k -> RestrictedPartitions( 20, [2,5,10], k ) );
#>[ [  ], [ [ 10, 10 ] ], [ [ 10, 5, 5 ] ], [ [ 5, 5, 5, 5 ] ], [  ], 
#>  [ [ 10, 2, 2, 2, 2, 2 ] ], [ [ 5, 5, 2, 2, 2, 2, 2 ] ], [  ], [  ], 
#>  [ [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ] ], [  ], [  ], [  ], [  ], [  ], [  ], 
#>  [  ], [  ], [  ], [  ] ]
RestrictedPartitions( 60, [2,3,5,7,11,13,17] )[600];
#>[ 13, 7, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ]
RestrictedPartitions( 100, [2,3,5,7,11,13,17], 10 )[75];
#>[ 17, 17, 13, 13, 13, 7, 5, 5, 5, 5 ]

#F  NrRestrictedPartitions(<n>,<set>) . . . . number of restricted partitions
NrRestrictedPartitions( 0, [1..10] );
#>1
List( [0..1], k -> NrRestrictedPartitions( 0, [1..10], k ) );
#>[ 1, 0 ]
NrRestrictedPartitions( 50, [1,2,5,10] );
#>341
List( [1..50], k -> NrRestrictedPartitions( 50, [1,2,5,10], k ) );
#>[ 0, 0, 0, 0, 1, 1, 1, 2, 4, 6, 6, 8, 10, 11, 11, 12, 13, 14, 14, 14, 15, 15, 
#>  14, 14, 14, 13, 12, 12, 11, 10, 9, 9, 8, 7, 6, 6, 6, 5, 4, 4, 4, 3, 2, 2, 
#>  2, 2, 1, 1, 1, 1 ]
NrRestrictedPartitions( 50, [2,5,10] );
#>21
List( [1..50], k -> NrRestrictedPartitions( 50, [2,5,10], k ) );
#>[ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 
#>  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
NrRestrictedPartitions( 60, [2,3,5,7,11,13,17] );
#>1213
NrRestrictedPartitions( 100, [2,3,5,7,11,13,17], 10 );
#>125

#F  Lucas(<P>,<Q>,<k>)  . . . . . . . . . . . . . . value of a lucas sequence
List( [0..10], i->Lucas(1,-2,i)[1] );
#>[ 0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341 ]
List( [0..10], i->Lucas(1,-2,i)[2] );
#>[ 2, 1, 5, 7, 17, 31, 65, 127, 257, 511, 1025 ]
List( [0..10], i->Lucas(1,-1,i)[1] );
#>[ 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55 ]
List( [0..10], i->Lucas(2,1,i)[1] );
#>[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
Lucas( 0, -4, 100 ) = [ 0, 2^101, 4^100 ];
#>true

#F  Fibonacci( <n> )  . . . . . . . . . . . . value of the Fibonacci sequence
List( [0..17], Fibonacci );
#>[ 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597 ]
Fibonacci( 333 );
#>1751455877444438095408940282208383549115781784912085789506677971125378

#F  Bernoulli( <n> )  . . . . . . . . . . . . value of the Bernoulli sequence
List( [0..14], Bernoulli );
#>[ 1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66, 0, -691/2730, 0, 7/6 ]
Bernoulli( 80 );
#>-4603784299479457646935574969019046849794257872751288919656867/230010

# thats it for the combinatorical package  ##################################
Print( "$Id: combinat.tst,v 1.1.1.1 1996/12/11 12:44:37 werner Exp $  ",
       QuoInt(700000000,time), " GAPstones\n" );
if IsBound( GAPSTONES )  then Add( GAPSTONES, QuoInt(700000000,time) );  fi;
