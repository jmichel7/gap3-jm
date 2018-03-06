SizeScreen( [ 70, 24 ] );
#>[ 70, 24 ]

InfoRead1 := Ignore;
#>function (...) internal; end

RequirePackage( "anupq" );
F := FreeGroup(2);
#>Group( f.1, f.2 )

F.relators := [ F.1^2, F.2^2, Comm(F.1,F.2) ];
#>[ f.1^2, f.2^2, f.1^-1*f.2^-1*f.1*f.2 ]

G := AgGroupFpGroup(F);
#>Group( f.1, f.2 )

a1 := GroupHomomorphismByImages( G, G, [G.1, G.2], [G.2, G.1 * G.2] );
#>GroupHomomorphismByImages( Group( f.1, f.2 ), Group( f.1, f.2 ), 
#>[ f.1, f.2 ], [ f.2, f.1*f.2 ] )

a2 := GroupHomomorphismByImages( G, G, [G.1, G.2], [G.2, G.1] );
#>GroupHomomorphismByImages( Group( f.1, f.2 ), Group( f.1, f.2 ), 
#>[ f.1, f.2 ], [ f.2, f.1 ] )

G.automorphisms := [ a1, a2 ];
#>[ GroupHomomorphismByImages( Group( f.1, f.2 ), Group( f.1, f.2 ), 
#>    [ f.1, f.2 ], [ f.2, f.1*f.2 ] ), 
#>  GroupHomomorphismByImages( Group( f.1, f.2 ), Group( f.1, f.2 ), 
#>    [ f.1, f.2 ], [ f.2, f.1 ] ) ]

L := PqDescendants( G, "OrderBound", 4, "ClassBound", 4 );
#>[ Group( G.1, G.2, G.3 ), Group( G.1, G.2, G.3 ), 
#>  Group( G.1, G.2, G.3, G.4 ), Group( G.1, G.2, G.3, G.4 ), 
#>  Group( G.1, G.2, G.3, G.4 ), Group( G.1, G.2, G.3, G.4 ), 
#>  Group( G.1, G.2, G.3, G.4 ) ]

List( L, x -> x.relators );
#>[ [ G.1^2*G.3^-1, G.2^2, G.3^2 ], 
#>  [ G.1^2, G.2^2, G.3^2, G.2^-1*G.1^-1*G.2*G.1*G.3^-1 ], 
#>  [ G.1^2*G.3^-1, G.2^2*G.4^-1, G.3^2, G.4^2 ], 
#>  [ G.1^2*G.4^-1, G.2^2, G.3^2, G.4^2, G.2^-1*G.1^-1*G.2*G.1*G.3^-1 
#>     ], 
#>  [ G.1^2*G.4^-1, G.2^2*G.3^-1, G.3^2, G.4^2, G.2^-1*G.1^-1*G.2*G.1*\
#>G.3^-1 ], [ G.1^2*G.3^-1, G.2^2, G.3^2*G.4^-1, G.4^2 ], 
#>  [ G.1^2, G.2^2, G.3^2*G.4^-1, G.4^2, G.2^-1*G.1^-1*G.2*G.1*G.3^-1,
#>      G.3^-1*G.1^-1*G.3*G.1*G.4^-1, G.3^-1*G.2^-1*G.3*G.2*G.4^-1 ] ]

List( L, x -> x.automorphisms );
#>[ [ GroupHomomorphismByImages( Group( G.1, G.2, G.3 ), Group( G.1, 
#>        G.2, G.3 ), [ G.1, G.2, G.3 ], [ G.1*G.3, G.2, G.3 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3 ), Group( 
#>        G.1, G.2, G.3 ), [ G.1, G.2, G.3 ], [ G.1, G.2*G.3, G.3 ] ),
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3 ), Group( 
#>        G.1, G.2, G.3 ), [ G.1, G.2, G.3 ], [ G.1*G.2, G.2, G.3 ] ) 
#>     ], 
#>  [ GroupHomomorphismByImages( Group( G.1, G.2, G.3 ), Group( G.1, 
#>        G.2, G.3 ), [ G.1, G.2, G.3 ], [ G.2, G.1, G.3 ] ) ], 
#>  [ GroupHomomorphismByImages( Group( G.1, G.2, G.3, G.4 ), Group( 
#>        G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.3, G.2, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.4, G.2, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1, G.2*G.3, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1, G.2*G.4, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.2, G.1*G.2, G.4, G.3*G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.2, G.1, G.4, G.3 ] ) ], 
#>  [ GroupHomomorphismByImages( Group( G.1, G.2, G.3, G.4 ), Group( 
#>        G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.4, G.2, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1, G.2*G.4, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.2, G.2, G.3, G.3*G.4 ] ) ], 
#>  [ GroupHomomorphismByImages( Group( G.1, G.2, G.3, G.4 ), Group( 
#>        G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.4, G.2, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1, G.2*G.4, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.2, G.2, G.3, G.4 ] ) ], 
#>  [ GroupHomomorphismByImages( Group( G.1, G.2, G.3, G.4 ), Group( 
#>        G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.4, G.2, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1, G.2*G.4, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.3, G.2, G.3*G.4, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1*G.2, G.2, G.3, G.4 ] ) ], 
#>  [ GroupHomomorphismByImages( Group( G.1, G.2, G.3, G.4 ), Group( 
#>        G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.1, G.2*G.4, G.3, G.4 ] ), 
#>      GroupHomomorphismByImages( Group( G.1, G.2, G.3, 
#>        G.4 ), Group( G.1, G.2, G.3, G.4 ), [ G.1, G.2, G.3, G.4 ], 
#>        [ G.2, G.1, G.3*G.4, G.4 ] ) ] ]

