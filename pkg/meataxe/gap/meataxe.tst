LogTo("testrun.out");
InfoRead1:= Print;
RequirePackage( "meataxe" );
a5:= Group( (1,2,3,4,5), (1,2,3) );;
Size( a5 );
f:= GF(2);;
m1:= MeatAxeMat( a5.1, f, [5,5] );
m2:= MeatAxeMat( a5.2, f, [5,5] );
Display( m1 );
a:= Algebra( f, [ m1, m2 ] );
nat:= NaturalModule( a );
fix:= FixedSubmodule( nat );
Dimension( fix );
quot:= nat / fix; 
op:= Operation( a, quot );
m:= NaturalModule( op );
IsIrreducible( m );
IsAbsolutelyIrreducible( m );
deg4:= m.ring;
tens:= KroneckerProduct( m, m );
comp:= CompositionFactors( tens );
IsIrreducible( comp[3] );
IsAbsolutelyIrreducible( comp[3] );
sf:= SplittingField( comp[3] );

comp[3].ring;
List( comp[3].ring.generators, x -> x.abstract );
comp[3].ring.1;
comp[3].ring.2;

new:= Algebra( sf, [ comp[3].ring.1, comp[3].ring.2 ] );

nat:= NaturalModule( new );
comp:= CompositionFactors( nat );
deg2:= List( comp, x -> x.ring ); 
repres:= [ a.1^0, a.1* a.2 * a.1^3, a.1, a.1^2 ];;
List( repres, OrderMeatAxeMat );
abstracts:= List( repres, x -> x.abstract );
mapped:= List( [ 1 .. 4 ],
  x-> MappedExpression( abstracts[x],
      a.freeAlgebra.generators, deg4.generators ) );
List( mapped, OrderMeatAxeMat );
List( mapped, BrauerCharacterValue );

deg2[1].generators;
List( deg2[1].generators, GapObject );

mapped:= List( [ 1 .. 4 ],
  x-> MappedExpression( abstracts[x],
      a.freeAlgebra.generators, deg2[1].generators ) );
List( mapped, BrauerCharacterValue );
FFList( GF( 8 ) );
f:= GF(2);;
m:= [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] * f.one;;
m1:= MeatAxeMat( m, "file2" );
p:= (1,2,3);;
m2:= MeatAxeMat( p, f, [ 3, 3 ], "file" );
Display( m2 );
n:= MeatAxeMat( "file", f, [ 3, 3 ] );; # just notify a matrix 
m:= [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ] * f.one;;
m1:= MeatAxeMat( m, "file2" );;
GapObject( m1 );
p:= (1,2,3);;
m2:= MeatAxeMat( p, f, [ 3, 3 ], "file" );
Display( m2 );
comp:= CompositionFactors( tens );;
f:= GF(2);;
a:= Algebra( f, [ MeatAxeMat( (1,2,3,4,5), f, [ 5, 5 ] ),
                   MeatAxeMat( (1,2)      , f, [ 5, 5 ] ) ] );;
Fingerprint( a );
f:= GF(2);;
a:= Algebra( f, [ MeatAxeMat( (1,2,3,4,5), f, [ 5, 5 ] ),
                   MeatAxeMat( (1,2)      , f, [ 5, 5 ] ) ] );;
nat:= NaturalModule( a );;
fix:= FixedSubmodule( nat );;
Dimension( fix );
g:= MeatAxeMat( (1,2,3,4,5), GF(2), [ 5, 5 ] );;
BrauerCharacterValue( g );

perm:= MeatAxePerm( (1,2,3) );
GapObject( perm );

perm:= MeatAxePerm( () );
GapObject( perm );


# MeatAxe.Unbind(); 
