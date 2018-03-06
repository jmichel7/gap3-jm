SizeScreen( [ 72, ] );;
InfoChevie:=Ignore;;
W := CoxeterGroup( "A", 4 );
#>CoxeterGroup("A",4)
w := Braid( W )( 1, 2, 3, 4 );
#>1234
w ^ 3;
#>121321432.343
w ^ 4;
#>w0.232432
w ^ -1;
#>(1234)^-1
CHEVIE.PrintGarside;
#>rec(
#>   )
CHEVIE.PrintGarside := rec(GAP:=true);
#>rec(
#>  GAP := true )
w;
#>B(1,2,3,4)
w ^ 3;
#>B(1,2,1,3,2,1,4,3,2,3,4,3)
w ^ -1;
#>B(1,2,3,4)^-1
CHEVIE.PrintGarside := rec();
#>rec(
#>   )
w ^ -1;
#>(1234)^-1
CHEVIE.PrintGarside := rec(Greedy:=true);;
w^-1;                                    
#>w0^-1.232432
W := CoxeterGroup( "A", 3 );;
B := Braid( W );
#>function ( arg ) ... end
B( W.generators[1] );
#>1
B( 2, 1, 2, 1, 1 );
#>121.1.1
B( [ 2, 1, 2, 1, 1 ], -1 );
#>w0^-1.121.1.1
W := CoxeterGroup( "A", 2 );;
a := Braid( W )( [1] );
#>1
b := Braid( W )( [2] );
#>2
a * b;
#>12
( a * b ) ^ 4;
#>w0^2.12
( a * b ) ^ -1;
#>w0^-1.2
a ^ b;
#>w0^-1.21.12
a / b;
#>w0^-1.2.21
a.monoid;
#>BraidMonoid(CoxeterGroup("A",2))
CHEVIE.PrintGarside := rec(GAP:=true);;
( a * b ) ^ -1;
#>B(1,2)^-1
CHEVIE.PrintGarside := rec();;
( a * b ) ^ -1;
#>(12)^-1
CHEVIE.PrintGarside := rec(Greedy:=true);;
( a * b ) ^ -1;
#>w0^-1.2
W := CoxeterGroup( "A", 3 );;
b := Braid( W )( [ 2, 1, 2, 1, 1 ] );
#>121.1.1
p := EltBraid( b );
#>( 1, 8)( 2, 7)( 3, 6)( 4,10)( 9,12)
CoxeterWord( W, p );
#>[ 1, 2, 1 ]
AsWord( b );
#>[ 1, 2, 1, 1, 1 ]
W := CoxeterGroup( "F", 4 );;
w:=[ 2, 3, 2, 3, 4, 3, 2, 1, 3, 4 ];;
GoodCoxeterWord( W, w );
#>[ [ [ 1, 2, 3, 4 ], 2 ], [ [ 3, 4 ], 4 ] ]
OrderPerm( EltWord( W, w ) );
#>6
Braid( W )( w ) ^ 6;
#>w0^2.343.343.343.343
GoodCoxeterWord( W, [ 3, 2, 3, 4, 3, 2, 1, 3, 4, 2 ] );
#>false
W := CoxeterGroup( [ [ 1, -1, 0 ], [ 0, 1, -1 ] ],
                   [ [ 1, -1, 0 ], [ 0, 1, -1 ] ] );;
gu3 := CoxeterCoset( W, -IdentityMat( 3 ) );
#>2A2.(q+1)
F4 := CoxeterGroup( "F", 4 );;
D4 := ReflectionSubgroup( F4, [ 1, 2, 16, 48 ] );
#>ReflectionSubgroup(CoxeterGroup("F",4), [ 1, 2, 9, 16 ])
PrintDiagram( D4 );
#>D4 9
#>    \
#>     1 - 16
#>    /
#>   2
3D4 := CoxeterCoset( D4, MatXPerm(D4, (2,9,16)) );
#>3D4<9,16,1,2>
f4coset := CoxeterCoset( CoxeterGroup( "F", 4 ) );
#>F4
w:=RepresentativeOperation(F4,[1,2,9,16],[1,9,16,2],OnTuples);
#>( 2, 9,16)( 3, 4,31)( 5,11,18)( 6,13,10)( 7,27,28)( 8,15,12)(14,22,20)
#>(17,19,21)(26,33,40)(29,35,42)(30,37,34)(32,39,36)(38,46,44)(41,43,45)
3d4again := CoxeterSubCoset( f4coset, [ 1, 2, 9, 16], w );
#>3D4<9,16,1,2>
PrintDiagram( 3d4again );
#>phi acts as ( 2, 9,16) on the component below
#>D4 9
#>    \
#>     1 - 2
#>    /
#>   16
WF := CoxeterCoset( CoxeterGroup( "A", 4 ), (1,4)(2,3) );
#>2A4
Display( InductionTable( CoxeterSubCoset( WF, [ 2, 3 ] ), WF ) );
#>Induction from 2A2<2,3>.(q-1)(q+1) to 2A4
#>      |111 21 3
#>________________
#>11111 |  1  . .
#>2111  |  .  1 .
#>221   |  1  . .
#>311   |  1  . 1
#>32    |  .  . 1
#>41    |  .  1 .
#>5     |  .  . 1
HF := CoxeterSubCoset( WF, [1, 2],
                LongestCoxeterElement( CoxeterGroup( WF ) ) );
#>2A2.(q+1)^2
Display( InductionTable( HF, WF ) );
#>Induction from 2A2.(q+1)^2 to 2A4
#>      |111 21  3
#>_________________
#>11111 | -1  .  .
#>2111  | -2 -1  .
#>221   | -1 -2  .
#>311   |  1  2 -1
#>32    |  . -2  1
#>41    |  .  1 -2
#>5     |  .  .  1
W := CoxeterCoset( CoxeterGroup( "A", 2, "A", 2 ), (1,3)(2,4) );
#>(A2xA2)
Display( InductionTable( CoxeterSubCoset( W, [ 1, 3 ] ), W ) );
#>Induction from (A1xA1)<1,3>.(q-1)(q+1) to (A2xA2)
#>    |11 2
#>__________
#>111 | 1 .
#>21  | 1 1
#>3   | . 1
W := CoxeterCoset( CoxeterGroup( "A", 2, "G", 2, "A", 2 ),(1,5,2,6) );
#>2(A2xA2)<1,2,5,6>xG2<3,4>
ReflectionName( W );
#>"2(A2xA2)<1,2,5,6>xG2<3,4>"
W := CoxeterCoset( CoxeterGroup( "A", 2, "A", 2 ), (1,3,2,4) );
#>2(A2xA2)
PrintDiagram( W );
#>phi permutes the next 2 components
#>phi^2 acts as (1,2) on the component below
#>A2 1 - 2
#>A2 3 - 4
W := CoxeterCoset( CoxeterGroup( "A", 2, "A", 2 ), (1,3,2,4) );
#>2(A2xA2)
ReflectionType( W );
#>[ rec(orbit := [ rec(rank    := 2,
#>      series  := "A",
#>      indices := [ 1, 2 ]), rec(rank    := 2,
#>      series  := "A",
#>      indices := [ 3, 4 ]) ],
#>      twist := (1,2)) ]
W := CoxeterGroup( "D", 4 );;
ChevieClassInfo( CoxeterCoset( W, (1,2,4) ) );
#>rec(
#>  classtext := [ [ 1 ], [  ], [ 1, 2, 3, 1, 2, 3 ], [ 3 ], [ 1, 3 ], 
#>      [ 1, 2, 3, 1, 2, 4, 3, 2 ], [ 1, 2, 3, 2 ] ],
#>  classnames := [ "C_3", "\\tilde A_2", "C_3+A_1", "\\tilde A_2+A_1", 
#>      "F_4", "\\tilde A_2+A_2", "F_4(a_1)" ],
#>  orders := [ 6, 3, 6, 6, 12, 3, 6 ],
#>  classes := [ 48, 16, 16, 48, 48, 8, 8 ],
#>  classparams := 
#>   [ [ "C_3" ], [ "\\tilde A_2" ], [ "C_3+A_1" ], [ "\\tilde A_2+A_1" 
#>         ], [ "F_4" ], [ "\\tilde A_2+A_2" ], [ "F_4(a_1)" ] ] )
W := CoxeterCoset( CoxeterGroup( "D", 4 ), (1,2,4) );
#>3D4
Display( CharTable( W ) );
#>3D4
#>
#>       2  2   2     2      2  2      3      3
#>       3  1   1     1      .  .      1      1
#>
#>         C3 ~A2 C3+A1 ~A2+A1 F4 ~A2+A2 F4(a1)
#>
#>.4        1   1     1      1  1      1      1
#>.1111    -1   1     1     -1  1      1      1
#>.22       .   2     2      . -1     -1     -1
#>11.2      .   .     .      . -1      3      3
#>1.3       1   1    -1     -1  .     -2      2
#>1.111    -1   1    -1      1  .     -2      2
#>1.21      .   2    -2      .  .      2     -2
#>
W := CoxeterGroup( "E", 6 );; WF := CoxeterCoset( W );
#>E6
phi := EltWord( W,
       [ 6, 5, 4, 2, 3, 1, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1 ] );;
HF := CoxeterSubCoset( WF, [ 2..5 ], phi );
#>3D4<2,3,4,5>.(q^2+q+1)
PrintDiagram( HF );
#>phi acts as (2,3,5) on the component below
#>D4 2
#>    \
#>     4 - 5
#>    /
#>   3
ReflectionDegrees( HF );
#>[ [ 1, E(3) ], [ 1, E(3)^2 ], [ 2, 1 ], [ 4, E(3) ], [ 6, 1 ], 
#>  [ 4, E(3)^2 ] ]
ReflectionDegrees( CoxeterGroup( HF ) );
#>[ 1, 1, 2, 4, 6, 4 ]
W := CoxeterGroup( "A", 3 );;
Display( CharTable( W ));
#>A3
#>
#>      2    3    2    3    .  2
#>      3    1    .    .    1  .
#>
#>        1111  211   22   31  4
#>     2P 1111 1111 1111   31 22
#>     3P 1111  211   22 1111  4
#>
#>1111       1   -1    1    1 -1
#>211        3   -1   -1    .  1
#>22         2    .    2   -1  .
#>31         3    1   -1    . -1
#>4          1    1    1    1  1
#>
W := CoxeterGroup( "G", 2);;
ct := CharTable( W );
#>CharTable( "G2" )
ct.classtext;
#>[ [  ], [ 2 ], [ 1 ], [ 1, 2 ], [ 1, 2, 1, 2 ], [ 1, 2, 1, 2, 1, 2 ] ]
ct.classnames;
#>[ "A0", "~A1", "A1", "G2", "A2", "A1+~A1" ]
ct.irredinfo;
#>[ rec(
#>      charparam := [ [ 1, 0 ] ],
#>      charname := "phi{1,0}" ), rec(
#>      charparam := [ [ 1, 6 ] ],
#>      charname := "phi{1,6}" ), rec(
#>      charparam := [ [ 1, 3, 1 ] ],
#>      charname := "phi{1,3}'" ), rec(
#>      charparam := [ [ 1, 3, 2 ] ],
#>      charname := "phi{1,3}''" ), rec(
#>      charparam := [ [ 2, 1 ] ],
#>      charname := "phi{2,1}" ), rec(
#>      charparam := [ [ 2, 2 ] ],
#>      charname := "phi{2,2}" ) ]
W := CoxeterGroup( "A", 3 );
#>CoxeterGroup("A",3)
List( Elements( W ), x -> ReflectionCharValue( W, x ) );
#>[ 3, 1, 1, 1, 0, 0, 0, -1, 0, -1, -1, 1, 1, -1, -1, -1, 0, 0, 0, 0, 
#>  -1, -1, 1, -1 ]
W := CoxeterGroup( "H", 4 );
#>CoxeterGroup("H",4)
ReflectionDegrees( W );
#>[ 2, 12, 20, 30 ]
Size( W );
#>14400
q := X( Rationals );; q.name := "q";;
FakeDegrees( CoxeterGroup( "A", 2 ), q );
#>[ q^3, q^2 + q, q^0 ]
FakeDegree( CoxeterGroup( "A", 2 ), [ [ 2, 1 ] ], q );
#>q^2 + q
LowestPowerFakeDegrees( CoxeterGroup( "D", 4 ) );
#>[ 6, 6, 7, 12, 4, 3, 6, 2, 2, 4, 1, 2, 0 ]
HighestPowerFakeDegrees( CoxeterGroup( "D", 4 ) );
#>[ 10, 10, 11, 12, 8, 9, 10, 6, 6, 8, 5, 6, 0 ]
LowestPowerGenericDegrees( CoxeterGroup( "D", 4 ) );
#>[ 6, 6, 7, 12, 3, 3, 6, 2, 2, 3, 1, 2, 0 ]
HighestPowerGenericDegrees( CoxeterGroup( "D", 4 ) );
#>[ 10, 10, 11, 12, 9, 9, 10, 6, 6, 9, 5, 6, 0 ]
 ChevieCharInfo( CoxeterGroup( "G", 2 ) );
#>rec(
#>  charparams := 
#>   [ [ [ 1, 0 ] ], [ [ 1, 6 ] ], [ [ 1, 3, 1 ] ], [ [ 1, 3, 2 ] ], 
#>      [ [ 2, 1 ] ], [ [ 2, 2 ] ] ],
#>  extRefl := [ 1, 5, 2 ],
#>  a := [ 0, 6, 1, 1, 1, 1 ],
#>  A := [ 0, 6, 5, 5, 5, 5 ],
#>  b := [ 0, 6, 3, 3, 1, 2 ],
#>  spaltenstein := 
#>   [ "1", "\\varepsilon", "\\varepsilon_l", "\\varepsilon_c", 
#>      "\\theta'", "\\theta''" ],
#>  charnames := [ "phi{1,0}", "phi{1,6}", "phi{1,3}'", "phi{1,3}''", 
#>      "phi{2,1}", "phi{2,2}" ],
#>  positionId := 1,
#>  positionDet := 2,
#>  B := [ 0, 6, 3, 3, 5, 4 ] )
W := CoxeterGroup( "E", 6 );;
ChevieCharInfo( W ).frame;
#>[ "1_p", "1_p'", "10_s", "6_p", "6_p'", "20_s", "15_p", "15_p'", 
#>  "15_q", "15_q'", "20_p", "20_p'", "24_p", "24_p'", "30_p", "30_p'", 
#>  "60_s", "80_s", "90_s", "60_p", "60_p'", "64_p", "64_p'", "81_p", 
#>  "81_p'" ]
G := ComplexReflectionGroup( 4 );
#>ComplexReflectionGroup(4)
ReflectionDegrees( G );
#>[ 4, 6 ]
Size( G );
#>24
q := X( Cyclotomics );; q.name := "q";;
FakeDegrees( G, q );
#>[ q^0, q^4, q^8, q^7 + q^5, q^5 + q^3, q^3 + q, q^6 + q^4 + q^2 ]
G := ComplexReflectionGroup( 4, 2, 3 );
#>ComplexReflectionGroup(4,2,3)
v := X( Cyclotomics );; v.name := "v";;
CH := Hecke( G, v );
#>Hecke(G423,v)
G := ComplexReflectionGroup( 4 );
#>ComplexReflectionGroup(4)
v := X( Cyclotomics );; v.name := "v";;
CH := Hecke( G, v );
#>Hecke(G4,v)
Display( CharTable( CH ) );
#>H(G4)
#>
#>          2 3     3   2     1        1      1            1
#>          3 1     1   .     1        1      1            1
#>
#>            .     z 212    12      z12      1           1z
#>         2P .     .   z     1        1    z12          z12
#>         3P .     z 212     z        .      .            z
#>         5P .     z 212    1z        1    z12           12
#>
#>phi{1,0}    1   v^6 v^3   v^2      v^8      v          v^7
#>phi{1,4}    1     1   1  E3^2     E3^2     E3           E3
#>phi{1,8}    1     1   1    E3       E3   E3^2         E3^2
#>phi{2,5}    2    -2   .     1       -1     -1            1
#>phi{2,3}    2 -2v^3   . E3^2v -E3^2v^4 v+E3^2 -v^4-E3^2v^3
#>phi{2,1}    2 -2v^3   .   E3v   -E3v^4   v+E3   -v^4-E3v^3
#>phi{3,2}    3  3v^2  -v     .        .    v-1      v^3-v^2
#>
W := CoxeterGroup( "D", 4 );
#>CoxeterGroup("D",4)
PrintArray( W.cartan );
#>[[ 2,  0, -1,  0],
#> [ 0,  2, -1,  0],
#> [-1, -1,  2, -1],
#> [ 0,  0, -1,  2]]
W := CoxeterGroup( "A", 2, "B", 2 );; PrintArray( W.cartan );
#>[[ 2, -1,  0,  0],
#> [-1,  2,  0,  0],
#> [ 0,  0,  2, -2],
#> [ 0,  0, -1,  2]]
W := CoxeterGroup( [ [ -1, 1, 0], [ 0, -1, 1 ] ],
                          [ [ -1, 1, 0], [ 0, -1, 1 ] ] );
#>CoxeterGroup([[-1,1,0],[0,-1,1]],[[-1,1,0],[0,-1,1]])
MatXPerm( W, W.generators[1] );
#>[ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ]
C := CartanMat( "F", 4 );;
PrintArray( C );
#>[[ 2, -1,  0,  0],
#> [-1,  2, -1,  0],
#> [ 0, -2,  2, -1],
#> [ 0,  0, -1,  2]]
CartanMat( "I", 2, 5 );
#>[ [ 2, E(5)^2+E(5)^3 ], [ E(5)^2+E(5)^3, 2 ] ]
C := [ [ 2, 0, -1 ], [ 0, 2, 0 ], [ -1, 0, 2 ] ];;
t:= ReflectionType( C );
#>[ rec(rank    := 2,
#>      series  := "A",
#>      indices := [ 1, 3 ]), rec(rank    := 1,
#>      series  := "A",
#>      indices := [ 2 ]) ]
PrintDiagram(t);
#>A2 1 - 3
#>A1 2
ReflectionName(t);
#>"A2xA1"
ReflectionName( ReflectionType( CartanMat( "I", 2, 7 ) ) );
#>"I2(7)"
PrintDiagram( ReflectionType( CartanMat( "E", 8) ) );
#>E8      2
#>        |
#>1 - 3 - 4 - 5 - 6 - 7 - 8
W := CoxeterGroup( "A", 4 );;
PrintArray( W.cartan );
#>[[ 2, -1,  0,  0],
#> [-1,  2, -1,  0],
#> [ 0, -1,  2, -1],
#> [ 0,  0, -1,  2]]
W.matgens;
#>[ [ [ -1, 0, 0, 0 ], [ 1, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ],
#>  [ [ 1, 1, 0, 0 ], [ 0, -1, 0, 0 ], [ 0, 1, 1, 0 ], [ 0, 0, 0, 1 ] ],
#>  [ [ 1, 0, 0, 0 ], [ 0, 1, 1, 0 ], [ 0, 0, -1, 0 ], [ 0, 0, 1, 1 ] ],
#>  [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 1 ], [ 0, 0, 0, -1 ] ] 
#> ]
W.roots;
#>[ [ 1, 0, 0, 0 ], [ 0, 1, 0, 0 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ], 
#>  [ 1, 1, 0, 0 ], [ 0, 1, 1, 0 ], [ 0, 0, 1, 1 ], [ 1, 1, 1, 0 ], 
#>  [ 0, 1, 1, 1 ], [ 1, 1, 1, 1 ], [ -1, 0, 0, 0 ], [ 0, -1, 0, 0 ], 
#>  [ 0, 0, -1, 0 ], [ 0, 0, 0, -1 ], [ -1, -1, 0, 0 ], 
#>  [ 0, -1, -1, 0 ], [ 0, 0, -1, -1 ], [ -1, -1, -1, 0 ], 
#>  [ 0, -1, -1, -1 ], [ -1, -1, -1, -1 ] ]
W := CoxeterGroup( "D", 4 );;
p := EltWord( W, [ 1, 3, 2, 1, 3 ] );
#>( 1,14,13, 2)( 3,17, 8,18)( 4,12)( 5,20, 6,15)( 7,10,11, 9)(16,24)
#>(19,22,23,21)
CoxeterWord( W, p );
#>[ 1, 3, 1, 2, 3 ]
LongestCoxeterWord( W );
#>[ 1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4 ]
w0 := LongestCoxeterElement( W );
#>( 1,13)( 2,14)( 3,15)( 4,16)( 5,17)( 6,18)( 7,19)( 8,20)( 9,21)(10,22)
#>(11,23)(12,24)
CoxeterLength( W, w0 );
#>12
List( Reflections( W ), i -> CoxeterWord( W, i ) );
#>[ [ 1 ], [ 2 ], [ 3 ], [ 4 ], [ 1, 3, 1 ], [ 2, 3, 2 ], [ 3, 4, 3 ], 
#>  [ 1, 2, 3, 1, 2 ], [ 1, 3, 4, 3, 1 ], [ 2, 3, 4, 3, 2 ], 
#>  [ 1, 2, 3, 4, 3, 1, 2 ], [ 3, 1, 2, 3, 4, 3, 1, 2, 3 ] ]
l := List( [1 .. W.N+1], x -> CoxeterElements( W, x-1 ) );;
List( l, Length );
#>[ 1, 4, 9, 16, 23, 28, 30, 28, 23, 16, 9, 4, 1 ]
W := CoxeterGroup( "G", 2 );;
EltWord( W, [1,2,1] );
#>( 1,12)( 2, 4)( 3, 9)( 6, 7)( 8,10)
W := CoxeterGroup( "A", 3 );;
w := ( 1,11)( 3,10)( 4, 9)( 5, 7)( 6,12);;
w in W;
#>true
CoxeterWord( W, w );
#>[ 1, 2, 3, 2, 1 ]
BrieskornNormalForm( W, w );
#>[ [ 1, 3 ], [ 2 ], [ 1, 3 ] ]
W := CoxeterGroup( "F", 4 );;
p := EltWord( W, [ 1, 2, 3, 4, 2 ] );
#>( 1,44,38,25,20,14)( 2, 5,40,47,48,35)( 3, 7,13,21,19,15)
#>( 4, 6,12,28,30,36)( 8,34,41,32,10,17)( 9,18)(11,26,29,16,23,24)
#>(27,31,37,45,43,39)(33,42)
CoxeterLength( W, p );
#>5
CoxeterWord( W, p );
#>[ 1, 2, 3, 2, 4 ]
W := CoxeterGroup( "E", 6 );;
ReducedCoxeterWord( W, [ 1, 1, 1, 1, 1, 2, 2, 2, 3 ] );
#>[ 1, 2, 3 ]
W := CoxeterGroup( "A", 2 );;
w := EltWord( W, [ 1, 2 ] );;
LeftDescentSet( W, w );
#>[ 1 ]
RightDescentSet( W, w );
#>[ 2 ]
W := CoxeterGroup( "B", 2 );; W.roots;
#>[ [ 1, 0 ], [ 0, 1 ], [ 1, 1 ], [ 2, 1 ], [ -1, 0 ], [ 0, -1 ], 
#>  [ -1, -1 ], [ -2, -1 ] ]
Reflections( W );
#>[ (1,5)(2,4)(6,8), (1,3)(2,6)(5,7), (2,8)(3,7)(4,6), (1,7)(3,5)(4,8) ]
LongestCoxeterElement( CoxeterGroup( "A", 4 ) );
#>( 1,14)( 2,13)( 3,12)( 4,11)( 5,17)( 6,16)( 7,15)( 8,19)( 9,18)(10,20)
LongestCoxeterWord( CoxeterGroup( "A", 4 ) );
#>[ 1, 2, 1, 3, 2, 1, 4, 3, 2, 1 ]
W := CoxeterGroup( "G", 2 );;  W.roots;
#>[ [ 1, 0 ], [ 0, 1 ], [ 1, 1 ], [ 1, 2 ], [ 1, 3 ], [ 2, 3 ], 
#>  [ -1, 0 ], [ 0, -1 ], [ -1, -1 ], [ -1, -2 ], [ -1, -3 ], 
#>  [ -2, -3 ] ]
HighestShortRoot( W );
#>4
W1 := CoxeterGroup( "A", 1, "B", 3 );;
W := CoxeterGroup( "G", 2 );;
e := CoxeterElements( W, 6 );
#>[ ( 1, 7)( 2, 8)( 3, 9)( 4,10)( 5,11)( 6,12) ]
e[1] = LongestCoxeterElement( W );
#>true
CoxeterWords( W );
#>[ [  ], [ 2 ], [ 1 ], [ 2, 1 ], [ 1, 2 ], [ 2, 1, 2 ], [ 1, 2, 1 ], 
#>  [ 2, 1, 2, 1 ], [ 1, 2, 1, 2 ], [ 2, 1, 2, 1, 2 ], 
#>  [ 1, 2, 1, 2, 1 ], [ 1, 2, 1, 2, 1, 2 ] ]
WordsClassRepresentatives( CoxeterGroup( "F", 4 ) );
#>[ [  ], 
#>  [ 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 
#>      3, 4 ], [ 2, 3, 2, 3 ], [ 2, 1 ], 
#>  [ 1, 2, 3, 4, 2, 3, 2, 3, 4, 3 ], 
#>  [ 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 ], [ 4, 3 ], 
#>  [ 1, 2, 1, 3, 2, 3, 1, 2, 3, 4 ], 
#>  [ 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 ], 
#>  [ 1, 2, 3, 4, 1, 2, 3, 4 ], [ 1, 2, 3, 4 ], [ 1 ], 
#>  [ 2, 3, 2, 3, 4, 3, 2, 3, 4 ], [ 1, 4, 3 ], [ 4, 3, 2 ], 
#>  [ 2, 3, 2, 1, 3 ], [ 3 ], [ 1, 2, 1, 3, 2, 1, 3, 2, 3 ], 
#>  [ 2, 1, 4 ], [ 3, 2, 1 ], [ 2, 4, 3, 2, 3 ], [ 1, 3 ], [ 3, 2 ], 
#>  [ 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3 ], [ 1, 2, 3, 4, 2, 3 ] ]
W := CoxeterGroup( "D", 4 );;
ChevieClassInfo( W );
#>rec(
#>  classtext := 
#>   [ [  ], [ 1, 2 ], [ 1, 2, 3, 1, 2, 3, 4, 3, 1, 2, 3, 4 ], [ 1 ], 
#>      [ 1, 2, 3 ], [ 1, 2, 4 ], [ 1, 4 ], [ 2, 4 ], 
#>      [ 1, 3, 1, 2, 3, 4 ], [ 1, 3 ], [ 1, 2, 3, 4 ], [ 1, 4, 3 ], 
#>      [ 2, 4, 3 ] ],
#>  classparams := 
#>   [ [ [ [ 1, 1, 1, 1 ], [  ] ] ], [ [ [ 1, 1 ], [ 1, 1 ] ] ], 
#>      [ [ [  ], [ 1, 1, 1, 1 ] ] ], [ [ [ 2, 1, 1 ], [  ] ] ], 
#>      [ [ [ 1 ], [ 2, 1 ] ] ], [ [ [ 2 ], [ 1, 1 ] ] ], 
#>      [ [ [ 2, 2 ], '+' ] ], [ [ [ 2, 2 ], '-' ] ], 
#>      [ [ [  ], [ 2, 2 ] ] ], [ [ [ 3, 1 ], [  ] ] ], 
#>      [ [ [  ], [ 3, 1 ] ] ], [ [ [ 4 ], '+' ] ], [ [ [ 4 ], '-' ] ] ]
#>   ,
#>  classnames := [ "1111.", "11.11", ".1111", "211.", "1.21", "2.11", 
#>      "22.+", "22.-", ".22", "31.", ".31", "4.+", "4.-" ],
#>  orders := [ 1, 2, 2, 2, 4, 2, 2, 2, 4, 3, 6, 4, 4 ],
#>  centralizers := [ 192, 32, 192, 16, 8, 16, 32, 32, 16, 6, 6, 8, 8 ],
#>  classes := [ 1, 6, 1, 12, 24, 12, 6, 6, 12, 32, 32, 24, 24 ] )
W := CoxeterGroup( "H", 3 );;
w := EltWord( W, [ 1, 2, 1, 3 ] );;
b := Filtered( Elements( W ), i -> Bruhat( W, i, w) );;
List( b, x -> CoxeterWord( W, x ) );
#>[ [  ], [ 3 ], [ 2 ], [ 1 ], [ 2, 1 ], [ 2, 3 ], [ 1, 3 ], [ 1, 2 ], 
#>  [ 2, 1, 3 ], [ 1, 2, 1 ], [ 1, 2, 3 ], [ 1, 2, 1, 3 ] ]
W := CoxeterGroup(
    [ [ 2, 0,-1, 0, 0, 0, 1 ], [ 0, 2, 0,-1, 0, 0, 0 ],
      [-1, 0, 2,-1, 0, 0,-1 ], [ 0,-1,-1, 2,-1, 0, 0 ],
      [ 0, 0, 0,-1, 2,-1, 0 ], [ 0, 0, 0, 0,-1, 2, 0 ] ],
    [ [ 1, 0, 0, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0, 0, 0 ],
      [ 0, 0, 1, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0, 0, 0 ],
      [ 0, 0, 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 0, 1, 0 ] ] );;
w0 := LongestCoxeterElement( W );;
mx := MatXPerm( W, w0 );
#>[ [ 0, 0, 0, 0, 0, -1, 1 ], [ 0, -1, 0, 0, 0, 0, 2 ], 
#>  [ 0, 0, 0, 0, -1, 0, 3 ], [ 0, 0, 0, -1, 0, 0, 4 ], 
#>  [ 0, 0, -1, 0, 0, 0, 3 ], [ -1, 0, 0, 0, 0, 0, 1 ], 
#>  [ 0, 0, 0, 0, 0, 0, 1 ] ]
my := MatYPerm( W, w0 );
#>[ [ 0, 0, 0, 0, 0, -1, 0 ], [ 0, -1, 0, 0, 0, 0, 0 ], 
#>  [ 0, 0, 0, 0, -1, 0, 0 ], [ 0, 0, 0, -1, 0, 0, 0 ], 
#>  [ 0, 0, -1, 0, 0, 0, 0 ], [ -1, 0, 0, 0, 0, 0, 0 ], 
#>  [ 1, 2, 3, 4, 3, 1, 1 ] ]
PermMatX( W, mx ) = w0;
#>true
PermMatY( W, my ) = w0;
#>true
W := CoxeterGroup( "A", 2 );;
q := X( Rationals );; q.name := "q";;
HF := Hecke( CoxeterCoset( W, (1,2) ), q^2, q );
#>Hecke(2A2,q^2,q)
Display( CharTable( HF ) );
#>H(2A2)
#>
#>     2     1   1   .
#>     3     1   .   1
#>
#>         111  21   3
#>    2P   111 111   3
#>    3P   111  21 111
#>
#>111       -1   1  -1
#>21     -2q^3   .   q
#>3        q^6   1 q^2
#>
W := CoxeterGroup( "A", 2, "A", 3 );;
HF := Hecke( CoxeterCoset( W, (1,2) ), q^2, q );
#>Hecke(2A2xA3,q^2,q)
Display( CharTable( HF ) );
#>H(2A2xA3)
#>
#>          2        4          3          4          1      3        4
#>          3        2          1          1          2      1        1
#>
#>            111,1111    111,211     111,22     111,31  111,4  21,1111
#>         2P 111,1111   111,1111   111,1111     111,31 111,22 111,1111
#>         3P 111,1111    111,211     111,22   111,1111  111,4  21,1111
#>
#>111,1111          -1          1         -1         -1      1        1
#>111,211           -3     -q^2+2     2q^2-1      q^2-1   -q^2        3
#>111,22            -2     -q^2+1     -q^4-1        q^2      .        2
#>111,31            -3    -2q^2+1  -q^4+2q^2   -q^4+q^2    q^4        3
#>111,4             -1       -q^2       -q^4       -q^4   -q^6        1
#>21,1111        -2q^3       2q^3      -2q^3      -2q^3   2q^3        .
#>21,211         -6q^3 -2q^5+4q^3  4q^5-2q^3  2q^5-2q^3  -2q^5        .
#>21,22          -4q^3 -2q^5+2q^3 -2q^7-2q^3       2q^5      .        .
#>21,31          -6q^3 -4q^5+2q^3 -2q^7+4q^5 -2q^7+2q^5   2q^7        .
#>21,4           -2q^3      -2q^5      -2q^7      -2q^7  -2q^9        .
#>3,1111           q^6       -q^6        q^6        q^6   -q^6        1
#>3,211           3q^6   q^8-2q^6  -2q^8+q^6   -q^8+q^6    q^8        3
#>3,22            2q^6    q^8-q^6   q^10+q^6       -q^8      .        2
#>3,31            3q^6   2q^8-q^6  q^10-2q^8   q^10-q^8  -q^10        3
#>3,4              q^6        q^8       q^10       q^10   q^12        1
#>
#>          2        3        4       1      3        3        2
#>          3        .        .       1      .        2        1
#>
#>              21,211    21,22   21,31   21,4   3,1111    3,211
#>         2P 111,1111 111,1111  111,31 111,22   3,1111   3,1111
#>         3P   21,211    21,22 21,1111   21,4 111,1111  111,211
#>
#>111,1111          -1        1       1     -1       -1        1
#>111,211        q^2-2  -2q^2+1  -q^2+1    q^2       -3   -q^2+2
#>111,22         q^2-1    q^4+1    -q^2      .       -2   -q^2+1
#>111,31        2q^2-1 q^4-2q^2 q^4-q^2   -q^4       -3  -2q^2+1
#>111,4            q^2      q^4     q^4    q^6       -1     -q^2
#>21,1111            .        .       .      .        q       -q
#>21,211             .        .       .      .       3q   q^3-2q
#>21,22              .        .       .      .       2q    q^3-q
#>21,31              .        .       .      .       3q   2q^3-q
#>21,4               .        .       .      .        q      q^3
#>3,1111            -1        1       1     -1      q^2     -q^2
#>3,211          q^2-2  -2q^2+1  -q^2+1    q^2     3q^2 q^4-2q^2
#>3,22           q^2-1    q^4+1    -q^2      .     2q^2  q^4-q^2
#>3,31          2q^2-1 q^4-2q^2 q^4-q^2   -q^4     3q^2 2q^4-q^2
#>3,4              q^2      q^4     q^4    q^6      q^2      q^4
#>
#>          2         3        .     2
#>          3         1        2     1
#>
#>                 3,22     3,31   3,4
#>         2P    3,1111     3,31  3,22
#>         3P    111,22 111,1111 111,4
#>
#>111,1111           -1       -1     1
#>111,211        2q^2-1    q^2-1  -q^2
#>111,22         -q^4-1      q^2     .
#>111,31      -q^4+2q^2 -q^4+q^2   q^4
#>111,4            -q^4     -q^4  -q^6
#>21,1111             q        q    -q
#>21,211        -2q^3+q   -q^3+q   q^3
#>21,22           q^5+q     -q^3     .
#>21,31        q^5-2q^3  q^5-q^3  -q^5
#>21,4              q^5      q^5   q^7
#>3,1111            q^2      q^2  -q^2
#>3,211       -2q^4+q^2 -q^4+q^2   q^4
#>3,22          q^6+q^2     -q^4     .
#>3,31         q^6-2q^4  q^6-q^4  -q^6
#>3,4               q^6      q^6   q^8
#>
W := CoxeterGroup( "A", 2 );
#>CoxeterGroup("A",2)
H := Hecke( W, 0 );
#>Hecke(A2,0)
T := Basis( H, "T" );
#>function ( arg ) ... end
el := CoxeterWords( W );
#>[ [  ], [ 2 ], [ 1 ], [ 2, 1 ], [ 1, 2 ], [ 1, 2, 1 ] ]
PrintArray(List(el,i->List(el,j->T(i)*T(j))));
#>[[      T(),      T(2),      T(1),    T(2,1),    T(1,2),  T(1,2,1)],
#> [     T(2),     -T(2),    T(2,1),   -T(2,1),  T(1,2,1), -T(1,2,1)],
#> [     T(1),    T(1,2),     -T(1),  T(1,2,1),   -T(1,2), -T(1,2,1)],
#> [   T(2,1),  T(1,2,1),   -T(2,1), -T(1,2,1), -T(1,2,1),  T(1,2,1)],
#> [   T(1,2),   -T(1,2),  T(1,2,1), -T(1,2,1), -T(1,2,1),  T(1,2,1)],
#> [ T(1,2,1), -T(1,2,1), -T(1,2,1),  T(1,2,1),  T(1,2,1), -T(1,2,1)]]
W := CoxeterGroup( "B", 3 );
#>CoxeterGroup("B",3)
u := X( Rationals );; u.name := "u";;
H := Hecke( W, u );
#>Hecke(B3,u)
H := Hecke( W, u^2, u );
#>Hecke(B3,u^2,u)
H := Hecke( W, [ u^6, u^4, u^4 ], [ u^3, -u^2, -u^2 ] );
#>Hecke(B3,[u^6,u^4,u^4],[u^3,-u^2,-u^2])
H := Hecke( W, 9, 3 );
#>Hecke(B3,9,3)
H := Hecke( CoxeterGroup( "B", 2 ), u );
#>Hecke(B2,u)
HeckeSubAlgebra( H, [ 1, 4 ] );
#>Hecke(B2,u)
T := Basis( H, "T" );;
W := CoxeterGroup( "B", 3 );;
H := Hecke( W, u );;
T := Basis( H, "T" );;
T( 1, 2 ) = T( [ 1, 2 ] );
#>true
T( 1, 2 ) = T( EltWord( W, [ 1, 2 ] ) );
#>true
l := [ [], [ 1, 2, 3 ], [ 1 ], [ 2 ], [ 3 ] ];;
pl := List( l, i -> EltWord( W, i ) );;
h := T( pl, [ u^100, 1/u^20, 1, -5, 0 ] );
#>u^100T()+T(1)-5T(2)+u^-20T(1,2,3)
h.elm;
#>[ (), ( 1, 4)( 2,11)( 3, 5)( 8, 9)(10,13)(12,14)(17,18), 
#>  ( 1,10)( 2, 6)( 5, 8)(11,15)(14,17), 
#>  ( 1,16,13,10, 7, 4)( 2, 8,12,11,17, 3)( 5, 9, 6,14,18,15) ]
h.coeff;
#>[ u^100, -5, 1, u^(-20) ]
CHEVIE.PrintHecke := rec(GAP:=true);
#>rec(
#>  GAP := true )
T( pl, [ u^100, 1/u^20, 1, -5, 0 ] );
#>u^100*T()+T(1)+u^-20*T(2)-5*T(1,2,3)
CHEVIE.PrintHecke := rec();;
T( pl, [ u^100, 1/u^20, 1, -5, 0 ] );
#>u^100T()+T(1)+u^-20T(2)-5T(1,2,3)
q := X( Rationals );; q.name := "q";;
H := Hecke( CoxeterGroup( "A", 2 ), q );
#>Hecke(A2,q)
T := Basis( H, "T" );;
( T() + T( 1 ) ) * ( T() + T( 2 ) );
#>T()+T(1)+T(2)+T(1,2)
T( 1 ) * T( 1 );
#>qT()+(q-1)T(1)
T( 1, 1 );
#>qT()+(q-1)T(1)
( q * T( 1, 2 ) ) ^ -1;
#>(q^-1-2q^-2+q^-3)T()+(-q^-2+q^-3)T(1)+(-q^-2+q^-3)T(2)+q^-3T(2,1)
( T( 1 ) + T( 2 ) ) ^ 2;
#>2qT()+(q-1)T(1)+(q-1)T(2)+T(1,2)+T(2,1)
T( 1 ) + T();
#>T()+T(1)
T( 1 ) - T( 1 );
#>0
h := T( [ EltWord( CoxeterGroup( H ), [ 1 ] ),() ], [ 0, q^100 ] );
#>q^100T()
AlphaInvolution( T( 1, 2 ) );
#>T(2,1)
CreateHeckeBasis( "t", rec(
    T := h->Basis( Hecke(h), "T" )( h.elm, List( [1 .. Length( h.elm )],
     i->Hecke(h).rootParameter[1]^CoxeterLength(
                 CoxeterGroup( Hecke(h) ), h.elm[i] ) *  h.coeff[i] ) ),
    t := h->Basis( Hecke(h), "t" )( h.elm, List( [1 .. Length( h.elm )],
     i->Hecke(h).rootParameter[1]^-CoxeterLength(
                  CoxeterGroup( Hecke(h) ), h.elm[i] ) * h.coeff[i] ) ),
    BetaInvolution := h->Basis( Hecke( h ),"t")(
                        HeckeAlgebraOps.T.BetaInvolution(
                                Basis( Hecke( h ), "T" )( h ) ) ) ) );
v := X( Rationals );; v.name := "v";;
H := Hecke( CoxeterGroup( "A", 3 ), v ^ 2, v );;
h := Basis( H, "t" )( 3, 1, 2 );
#>t(1,3,2)
h1 := Basis( H, "T")( h );
#>v^3T(1,3,2)
h2 := Basis( H, "t" )( h1 );
#>t(1,3,2)
BetaInvolution( h2 );
#>v^-12t(2,1,3)
H := Hecke(CoxeterGroup( "B", 2) , v^2, v);
#>Hecke(B2,v^2,v)
ref:= HeckeReflectionRepresentation( H );
#>[ [ [ -v^0, 0*v^0 ], [ -v^2, v^2 ] ], 
#>  [ [ v^2, -2*v^0 ], [ 0*v^0, -v^0 ] ] ]
H := Hecke( CoxeterGroup( "H", 3 ));;
HeckeReflectionRepresentation( H );
#>[ [ [ -1, 0, 0 ], [ -1, 1, 0 ], [ 0, 0, 1 ] ], 
#>  [ [ 1, E(5)+2*E(5)^2+2*E(5)^3+E(5)^4, 0 ], [ 0, -1, 0 ], 
#>      [ 0, -1, 1 ] ], [ [ 1, 0, 0 ], [ 0, 1, -1 ], [ 0, 0, -1 ] ] ]
H := Hecke(CoxeterGroup( "F", 4 ));;
r := HeckeReflectionRepresentation( H );;
CheckHeckeDefiningRelations( H, r );
#>true
W := CoxeterGroup( "G", 2 );;
u := X( Rationals );;  u.name := "u";;
v := X( LaurentPolynomialRing( Rationals ) );; v.name := "v";;
u := u * v^0;;
H := Hecke( W, [ u^2, v^2 ], [ u, v ] );
#>Hecke(G2,[u^2,v^2],[u,v])
Display( CharTable( H ) );
#>H(G2)
#>
#>            2  2     2     2      1       1        2
#>            3  1     .     .      1       1        1
#>
#>              A0   ~A1    A1     G2      A2   A1+~A1
#>           2P A0    A0    A0     A2      A2       A0
#>           3P A0   ~A1    A1 A1+~A1      A0   A1+~A1
#>
#>phi{1,0}       1   v^2   u^2 u^2v^2  u^4v^4   u^6v^6
#>phi{1,6}       1    -1    -1      1       1        1
#>phi{1,3}'      1   v^2    -1   -v^2     v^4     -v^6
#>phi{1,3}''     1    -1   u^2   -u^2     u^4     -u^6
#>phi{2,1}       2 v^2-1 u^2-1     uv -u^2v^2 -2u^3v^3
#>phi{2,2}       2 v^2-1 u^2-1    -uv -u^2v^2  2u^3v^3
#>
H1 := Hecke( W, [ E(6)^2, E(6)^4 ],[ E(6), E(6)^2 ] );
#>Hecke(G2,[E3,E3^2],[-E3^2,E3])
ct := CharTable( H1 );
#>CharTable( "H(G2)" )
Display( ct );
#>H(G2)
#>
#>            2  2             2             2      1    1      2
#>            3  1             .             .      1    1      1
#>
#>              A0           ~A1            A1     G2   A2 A1+~A1
#>           2P A0            A0            A0     A2   A2     A0
#>           3P A0           ~A1            A1 A1+~A1   A0 A1+~A1
#>
#>phi{1,0}       1          E3^2            E3      1    1      1
#>phi{1,6}       1            -1            -1      1    1      1
#>phi{1,3}'      1          E3^2            -1  -E3^2   E3     -1
#>phi{1,3}''     1            -1            E3    -E3 E3^2     -1
#>phi{2,1}       2 (-3-ER(-3))/2 (-3+ER(-3))/2     -1   -1      2
#>phi{2,2}       2 (-3-ER(-3))/2 (-3+ER(-3))/2      1   -1     -2
#>
RankMat( ct.irreducibles );
#>5
q := X( Rationals );; q.name := "q";;
H := Hecke( CoxeterGroup( "B", 2 ), q ^ 2, q );;
HeckeCharValues( Basis( H, "C'" )( 1, 2, 1 ) );
#>[ -q - q^(-1), q + q^(-1), 0*q^0, q^3 + 2*q + 2*q^(-1) + q^(-3), 
#>  0*q^0 ]
u := X( Rationals );; u.name := "u";;
W := CoxeterGroup( "A", 3 );
#>CoxeterGroup("A",3)
H := Hecke( W, u );;
h := Basis( H, "T" )( LongestCoxeterElement( W ) );
#>T(1,2,1,3,2,1)
cp := HeckeClassPolynomials( h );
#>[ 0*u^0, 0*u^0, u^2, u^3 - 2*u^2 + u, u^3 - u^2 + u - 1 ]
CharTable( H ).irreducibles * cp;
#>[ u^0, -u^2, 2*u^3, -u^4, u^6 ]
q := X( Rationals );; q.name := "q";;
W := CoxeterGroup( "G", 2 );; H := Hecke( W, q );
#>Hecke(G2,q)
PoincarePolynomial( H );
#>q^6 + 2*q^5 + 2*q^4 + 2*q^3 + 2*q^2 + 2*q + 1
u := X( Rationals );; u.name := "u";;
v := X( LaurentPolynomialRing( Rationals ) );; v.name := "v";;
H := Hecke( CoxeterGroup( "G", 2 ), [ u ^ 2, v ^ 2 ], [ u, v ] );
#>Hecke(G2,[u^2,v^2],[u,v])
schur := SchurElements( H);
#>[ (u^6 + u^4)*v^6 + (u^6 + 2*u^4 + u^2)*v^4 + (u^4 + 2*u^2 + 1)*v^
#>    2 + (u^2 + 1), (1 + u^(-2)) + (1 + 2*u^(-2) + u^(-4))*v^(
#>    -2) + (u^(-2) + 2*u^(-4) + u^(-6))*v^(-4) + (u^(-4) + u^(-6))*v^(
#>    -6), (u^(-4) + u^(-6))*v^6 + (u^(-2) + 2*u^(-4) + u^(-6))*v^4 + (
#>    1 + 2*u^(-2) + u^(-4))*v^2 + (1 + u^(-2)), 
#>  (u^2 + 1) + (u^4 + 2*u^2 + 1)*v^(-2) + (u^6 + 2*u^4 + u^2)*v^(
#>    -4) + (u^6 + u^4)*v^(-6), (2*u^0)*v^2 + (2*u - 2*u^(-1))*v + (2*u^
#>    2 - 2 + 2*u^(-2)) + (-2*u + 2*u^(-1))*v^(-1) + (2*u^0)*v^(-2), 
#>  (2*u^0)*v^2 + (-2*u + 2*u^(-1))*v + (2*u^2 - 2 + 2*u^(-2)) + (2*u - 
#>    2*u^(-1))*v^(-1) + (2*u^0)*v^(-2) ]
schur[1];
#>(u^6 + u^4)*v^6 + (u^6 + 2*u^4 + u^2)*v^4 + (u^4 + 2*u^2 + 1)*v^
#>2 + (u^2 + 1)
SchurElement( H, [ [ 1, 3, 1 ] ] );
#>(u^(-4) + u^(-6))*v^6 + (u^(-2) + 2*u^(-4) + u^(-6))*v^4 + (1 + 2*u^(
#>-2) + u^(-4))*v^2 + (1 + u^(-2))
v := X( Cyclotomics );; v.name := "v";;
H := Hecke( CoxeterGroup( "H", 3 ), v ^ 2, v );;
HeckeCentralMonomials( H );
#>[ v^0, v^60, v^24, v^36, v^20, v^20, v^40, v^40, v^30, v^30 ]
HeckeCharValuesGood( H, [ 1, 2, 3 ] );
#>[ v^0, v^60, 5*v^24, 5*v^36, 3*v^20, 3*v^20, 3*v^40, 3*v^40, 4*v^30, 
#>  4*v^30 ]
klpol := function( W, x, y)
      return KazhdanLusztigPolynomial( W, EltWord( W, x ), EltWord( W, y ));
    end;
#>function ( W, x, y ) ... end
q := X( Rationals );; q.name := "q";;
W := CoxeterGroup( "B", 3 );;
el := CoxeterWords( W );
#>[ [  ], [ 3 ], [ 2 ], [ 1 ], [ 3, 2 ], [ 2, 1 ], [ 2, 3 ], [ 1, 3 ], 
#>  [ 1, 2 ], [ 2, 1, 2 ], [ 3, 2, 1 ], [ 2, 3, 2 ], [ 2, 1, 3 ], 
#>  [ 1, 2, 1 ], [ 1, 3, 2 ], [ 1, 2, 3 ], [ 3, 2, 1, 2 ], 
#>  [ 2, 1, 2, 3 ], [ 2, 3, 2, 1 ], [ 2, 1, 3, 2 ], [ 1, 2, 1, 2 ], 
#>  [ 1, 3, 2, 1 ], [ 1, 2, 1, 3 ], [ 1, 2, 3, 2 ], [ 3, 2, 1, 2, 3 ], 
#>  [ 2, 1, 2, 3, 2 ], [ 2, 3, 2, 1, 2 ], [ 2, 1, 3, 2, 1 ], 
#>  [ 1, 3, 2, 1, 2 ], [ 1, 2, 1, 2, 3 ], [ 1, 2, 1, 3, 2 ], 
#>  [ 1, 2, 3, 2, 1 ], [ 2, 3, 2, 1, 2, 3 ], [ 2, 1, 2, 3, 2, 1 ], 
#>  [ 2, 1, 3, 2, 1, 2 ], [ 1, 3, 2, 1, 2, 3 ], [ 1, 2, 1, 2, 3, 2 ], 
#>  [ 1, 2, 1, 3, 2, 1 ], [ 1, 2, 3, 2, 1, 2 ], [ 2, 1, 2, 3, 2, 1, 2 ],
#>  [ 2, 1, 3, 2, 1, 2, 3 ], [ 1, 2, 3, 2, 1, 2, 3 ], 
#>  [ 1, 2, 1, 2, 3, 2, 1 ], [ 1, 2, 1, 3, 2, 1, 2 ], 
#>  [ 2, 1, 2, 3, 2, 1, 2, 3 ], [ 1, 2, 1, 2, 3, 2, 1, 2 ], 
#>  [ 1, 2, 1, 3, 2, 1, 2, 3 ], [ 1, 2, 1, 2, 3, 2, 1, 2, 3 ] ]
List( el, w -> klpol( W, [], w ) );
#>[ [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], 
#>  [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], 
#>  [ 1 ], [ 1, 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1 ], [ 1, 1 ], [ 1 ], 
#>  [ 1 ], [ 1, 1 ], [ 1 ], [ 1 ], [ 1, 1 ], [ 1, 1 ], [ 1 ], [ 1, 1 ], 
#>  [ 1 ], [ 1, 1 ], [ 1 ], [ 1, 0, 1 ], [ 1, 1 ], [ 1, 1, 1 ], 
#>  [ 1, 1 ], [ 1, 1 ], [ 1 ], [ 1 ], [ 1, 0, 1 ], [ 1 ], [ 1, 1 ], 
#>  [ 1 ] ]
W.klpol;
#>Dictionary with 19 entries
W := CoxeterGroup( "F", 4 );;
w := LongestCoxeterElement( W ) * W.generators[1];;
CoxeterLength( W, w );
#>23
y := EltWord( W, [ 1, 2, 3, 4 ] );;
cr := CriticalPair( W, y, w);;
CoxeterWord( W, cr);
#>[ 2, 3, 2, 1, 3, 4, 3, 2, 1, 3, 2, 3, 4, 3, 2, 3 ]
KazhdanLusztigPolynomial( W, y, w );
#>[ 1, 0, 0, 1 ]
KazhdanLusztigPolynomial( W, cr, w);
#>[ 1, 0, 0, 1 ]
W := CoxeterGroup( "B", 4 );;
y := [ 1, 2, 3, 4, 3, 2, 1 ];;
py := EltWord( W, y );
#>( 1,28)( 2,15)( 4,27)( 6,16)( 7,24)( 8,23)(11,20)(12,17)(14,30)(18,31)
#>(22,32)
x := [ 1 ];;
px := EltWord( W, x );
#>( 1,17)( 2, 8)( 6,11)(10,14)(18,24)(22,27)(26,30)
Bruhat( W, px, py );
#>true
List([0..3],i->KazhdanLusztigCoefficient( W, px, py, i ) );
#>[ 1, 2, 1, 0 ]
W := CoxeterGroup( "G", 2 );;
LeftCells(W);
#>[ LeftCell<G2: duflo= character=phi{1,0}>, 
#>  LeftCell<G2: duflo=1,2 character=phi{1,6}>, 
#>  LeftCell<G2: duflo=2 character=phi{2,1}+phi{1,3}'+phi{2,2}>, 
#>  LeftCell<G2: duflo=1 character=phi{2,1}+phi{1,3}''+phi{2,2}> ]
v := X( Cyclotomics ) ;; v.name := "v";;
H := Hecke(CoxeterGroup( "H", 3), v^2, v );
#>Hecke(H3,v^2,v)
c := LeftCells( CoxeterGroup( H ) );;
List( c, Size);
#>[ 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 8, 8, 8, 8 ]
Representation(c[2],H);
#>[ [ [ v^2 ] ], [ [ v^2 ] ], [ [ v^2 ] ] ]
W := CoxeterGroup( "B", 3 );;
v := X( Rationals );; v.name := "v";;
H := Hecke( W, v^2, v );
#>Hecke(B3,v^2,v)
T := Basis( H, "T" );;
C := Basis( H, "C" );;
T( C( 1 ) );
#>-vT()+v^-1T(1)
C( T( 1 ) );
#>v^2C()+vC(1)
ref := HeckeReflectionRepresentation( H );;
c := CharRepresentationWords( ref, WordsClassRepresentatives( W ) );
#>[ 3*v^0, 2*v^2 - 1, v^8 - 2*v^4, -3*v^12, 2*v^2 - 1, v^4, 
#>  v^4 - 2*v^2, -v^6, v^4 - v^2, 0*v^0 ]
List( ChevieClassInfo( W ).classtext, i ->
                                HeckeCharValues( C( i ), c ) );
#>[ 3*v^0, -v - v^(-1), 0*v^0, 0*v^0, -v - v^(-1), 2*v^0, 0*v^0, 0*v^0, 
#>  v^0, 0*v^0 ]
v := X( Rationals );; v.name := "v";;
H := Hecke( CoxeterGroup( "B", 2 ), v ^ 2, v );;
h := Basis( H, "C'" )( 1 );
#>C'(1)
h2 := h * h;
#>(v+v^-1)C'(1)
Basis( H, "T" )( h2 );
#>(1+v^-2)T()+(1+v^-2)T(1)
W := CoxeterGroup( "B", 2 );;
v := X( Rationals );; v.name := "v";;
H := Hecke( W, v^2, v );
#>Hecke(B2,v^2,v)
T := Basis( H, "T" );;
D := Basis( H, "D" );;
D( T( 1 ) );
#>vD(1)-v^2D(1,2)-v^2D(2,1)+v^3D(1,2,1)+v^3D(2,1,2)-v^4D(1,2,1,2)
BetaInvolution( D( 1 ) );
#>C'(2,1,2)
Dp := Basis( H, "D'" );;
AltInvolution( Dp( 1 ) );
#>D(1)
Dp( 1 )^3;
#>(v+2v^-1-5v^-5-9v^-7-8v^-9-4v^-11-v^-13)D'()+(v^2+2+v^-2)D'(1)
W := CoxeterGroup( "G", 2 );;
W.roots[4];
#>[ 1, 2 ]
R := ReflectionSubgroup( W, [ 2, 4 ] );
#>ReflectionSubgroup(CoxeterGroup("G",2), [ 2, 3 ])
PrintDiagram( R );
#>~A2 2 - 3
R.rootInclusion;
#>[ 2, 3, 4, 8, 9, 10 ]
CoxeterLength( W, R.generators[2] );
#>3
CoxeterLength( R, R.generators[2] );
#>1
el := CoxeterWords( R );
#>[ [  ], [ 2 ], [ 3 ], [ 2, 3 ], [ 3, 2 ], [ 2, 3, 2 ] ]
el1 := List( el, x -> R.rootRestriction{ x } );
#>[ [  ], [ 1 ], [ 2 ], [ 1, 2 ], [ 2, 1 ], [ 1, 2, 1 ] ]
el2 := List( el, x -> EltWord( R, x ) );
#>[ (), ( 1, 5)( 2, 8)( 3, 4)( 7,11)( 9,10), 
#>  ( 1,12)( 2, 4)( 3, 9)( 6, 7)( 8,10), 
#>  ( 1, 5,12)( 2,10, 3)( 4, 9, 8)( 6, 7,11), 
#>  ( 1,12, 5)( 2, 3,10)( 4, 8, 9)( 6,11, 7), 
#>  ( 2, 9)( 3, 8)( 4,10)( 5,12)( 6,11) ]
List( el2, x -> CoxeterWord( W, x ) );
#>[ [  ], [ 2 ], [ 1, 2, 1 ], [ 2, 1, 2, 1 ], [ 1, 2, 1, 2 ], 
#>  [ 2, 1, 2, 1, 2 ] ]
W := CoxeterGroup( "G", 2 );;
W.roots[4];
#>[ 1, 2 ]
R := ReflectionSubgroup( W, [ 2, 4 ] );;
Display( InductionTable( R, W ) );
#>Induction from ~A2<2,3> to G2
#>           |111 21 3
#>_____________________
#>phi{1,0}   |  .  . 1
#>phi{1,6}   |  1  . .
#>phi{1,3}'  |  .  . 1
#>phi{1,3}'' |  1  . .
#>phi{2,1}   |  .  1 .
#>phi{2,2}   |  .  1 .
W := CoxeterGroup( "F", 4 );;
H := ReflectionSubgroup( W, [ 1, 2, 11, 20 ] );
#>ReflectionSubgroup(CoxeterGroup("F",4), [ 1, 2, 9, 16 ])
ReflectionName( H );
#>"D4<9,2,1,16>"
PrintDiagram( H );
#>D4 9
#>    \
#>     1 - 16
#>    /
#>   2
H.rootRestriction;
#>[ 1, 2,,, 5,,,, 3,, 6,,, 8,, 4,, 7,, 9,, 10, 11, 12, 13, 14,,, 17,,,, 
#>  15,, 18,,, 20,, 16,, 19,, 21,, 22, 23, 24 ]
w := EltWord( W, [ 3, 2, 3, 4, 3, 2 ] );;
f := ReducedInRightCoset( H, w );;
CoxeterWord( W, f );
#>[ 4, 3 ]
H.rootInclusion{[ 1 ..4 ]};
#>[ 1, 2, 9, 16 ]
OnTuples( H.rootInclusion{[ 1 .. 4 ]}, f );
#>[ 1, 9, 16, 2 ]
ReflectionSubgroup( H, [ 1, 2, 6 ] );
#>ReflectionSubgroup(CoxeterGroup("F",4), [ 1, 2, 3 ])
W := CoxeterGroup( "F", 4 );;
H := ReflectionSubgroup( W, [ 10, 11, 12 ] );;
PrintDiagram( H );
#>B2 10 <=< 11
#>~A1 12
LongestCoxeterWord( H );
#>[ 10, 11, 10, 11, 12 ]
W := CoxeterGroup( "B", 3 );;
H := ReflectionSubgroup(W, [ 2, 3 ]);;
List( ReducedRightCosetRepresentatives( W, H ),
                                   x-> CoxeterWord( W, x ) );
#>[ [  ], [ 1 ], [ 1, 2 ], [ 1, 2, 1 ], [ 1, 2, 3 ], [ 1, 2, 1, 3 ], 
#>  [ 1, 2, 1, 3, 2 ], [ 1, 2, 1, 3, 2, 1 ] ]
W := CoxeterGroup( "F", 4 );;
PermCosetsSubgroup( W, ReflectionSubgroup( W, [ 1, 2, 3 ] ) );
#>[ ( 4, 5)( 6, 7)( 8,10)(16,18)(17,20)(19,21), 
#>  ( 3, 4)( 7, 9)(10,12)(14,16)(15,17)(21,22), 
#>  ( 2, 3)( 4, 6)( 5, 7)( 9,11)(12,14)(13,15)(17,19)(20,21)(22,23), 
#>  ( 1, 2)( 6, 8)( 7,10)( 9,12)(11,13)(14,15)(16,17)(18,20)(23,24) ]
W := CoxeterGroup( "D", 4);;
H := ReflectionSubgroup( W, [ 1, 3 ] );;
Display( jInductionTable( H, W ) );
#>j-Induction from A2<1,3>.(q-1)^2 to D4
#>      |111 21 3
#>________________
#>11+   |  .  . .
#>11-   |  .  . .
#>1.111 |  .  . .
#>.1111 |  .  . .
#>11.2  |  .  . .
#>1.21  |  1  . .
#>.211  |  .  . .
#>2+    |  .  . .
#>2-    |  .  . .
#>.22   |  .  . .
#>1.3   |  .  1 .
#>.31   |  .  . .
#>.4    |  .  . 1
W := CoxeterGroup( "D", 4 );;
H := ReflectionSubgroup( W, [ 1, 3 ] );;
Display( JInductionTable( H, W ) );
#>J-Induction from A2<1,3>.(q-1)^2 to D4
#>      |111 21 3
#>________________
#>11+   |  .  . .
#>11-   |  .  . .
#>1.111 |  .  . .
#>.1111 |  .  . .
#>11.2  |  1  . .
#>1.21  |  1  . .
#>.211  |  .  . .
#>2+    |  .  . .
#>2-    |  .  . .
#>.22   |  .  . .
#>1.3   |  .  1 .
#>.31   |  .  . .
#>.4    |  .  . 1
W := Group( [ (1,2), (2,3), (3,4) ], () );
#>Group( (1,2), (2,3), (3,4) )
H:=Subgroup( W, [ (1,2), (3,4) ] );
#>Subgroup( Group( (1,2), (2,3), (3,4) ), [ (1,2), (3,4) ] )
W.name := "W";; H.name := "H";;
Display( InductionTable( H, W ) );
#>Induction from H to W
#>    |X.1 X.2 X.3 X.4
#>_____________________
#>X.1 |  1   .   .   .
#>X.2 |  .   .   .   1
#>X.3 |  1   .   .   1
#>X.4 |  .   1   1   1
#>X.5 |  1   1   1   .
W := CoxeterGroup( "G", 2 );;
H := ReflectionSubgroup( W, [ 1, 4 ] );
#>ReflectionSubgroup(CoxeterGroup("G",2), [ 1, 4 ])
ReflectionName( H );
#>"A1x~A1<4>"
Display(InductionTable( H, W ) );
#>Induction from A1x~A1<4> to G2
#>           |11,11 11,2 2,11 2,2
#>________________________________
#>phi{1,0}   |    .    .    .   1
#>phi{1,6}   |    1    .    .   .
#>phi{1,3}'  |    .    1    .   .
#>phi{1,3}'' |    .    .    1   .
#>phi{2,1}   |    .    1    1   .
#>phi{2,2}   |    1    .    .   1
H := Hecke(CoxeterGroup( "F", 4 ));;
r := ChevieClassInfo( Group( H ) ).classtext;;
t := HeckeReflectionRepresentation( H );;
CharRepresentationWords( t, r );
#>[ 4, -4, 0, 1, -1, 0, 1, -1, -2, 2, 0, 2, -2, -1, 1, 0, 2, -2, -1, 1, 
#>  0, 0, 2, -2, 0 ]
G := Group( (1,2)(3,4), (1,2,3,4,5) );;
ConjugacyClasses( G );
#>[ ConjugacyClass( Group( (1,2)(3,4), (1,2,3,4,5) ), () ), 
#>  ConjugacyClass( Group( (1,2)(3,4), (1,2,3,4,5) ), (3,4,5) ), 
#>  ConjugacyClass( Group( (1,2)(3,4), (1,2,3,4,5) ), (2,3)(4,5) ), 
#>  ConjugacyClass( Group( (1,2)(3,4), (1,2,3,4,5) ), (1,2,3,4,5) ), 
#>  ConjugacyClass( Group( (1,2)(3,4), (1,2,3,4,5) ), (1,2,3,5,4) ) ]
g:=(1,5,3,2,4);
#>(1,5,3,2,4)
PositionClass( G, g );
#>4
G := Group( (1,7)(2,3)(5,6)(8,9)(11,12),
                   (1,5)(2,8)(3,4)(7,11)(9,10) );;
PointsAndRepresentativesOrbits( G );
#>[ [ [ 1, 7, 5, 11, 6, 12 ], [ 2, 3, 8, 4, 9, 10 ] ], 
#>  [ [ (), ( 1, 7)( 2, 3)( 5, 6)( 8, 9)(11,12), 
#>          ( 1, 5)( 2, 8)( 3, 4)( 7,11)( 9,10), 
#>          ( 1,11,12, 7, 5, 6)( 2, 4, 3, 8,10, 9), 
#>          ( 1, 6, 5, 7,12,11)( 2, 9,10, 8, 3, 4), 
#>          ( 1,12)( 2, 4)( 3, 9)( 6, 7)( 8,10) ], 
#>      [ (), ( 1, 7)( 2, 3)( 5, 6)( 8, 9)(11,12), 
#>          ( 1, 5)( 2, 8)( 3, 4)( 7,11)( 9,10), 
#>          ( 1,11,12, 7, 5, 6)( 2, 4, 3, 8,10, 9), 
#>          ( 1, 6, 5, 7,12,11)( 2, 9,10, 8, 3, 4), 
#>          ( 1, 6)( 2,10)( 4, 8)( 5,11)( 7,12) ] ] ]
C1 := [ [   2,  -1,   0,   0 ],
               [  -1,   2,  -1,   0 ],
               [   0,  -1,   2,  -1 ],
               [   0,   0,  -1,   2 ] ];;
C2 := [ [   2,   0,  -1,   0 ],
               [   0,   2,  -1,   0 ],
               [  -1,  -1,   2,  -1 ],
               [   0,   0,  -1,   2 ] ];;
d := DiagonalMat( C1, C2 );;
PrintArray( d );
#>[[ 2, -1,  0,  0,  0,  0,  0,  0],
#> [-1,  2, -1,  0,  0,  0,  0,  0],
#> [ 0, -1,  2, -1,  0,  0,  0,  0],
#> [ 0,  0, -1,  2,  0,  0,  0,  0],
#> [ 0,  0,  0,  0,  2,  0, -1,  0],
#> [ 0,  0,  0,  0,  0,  2, -1,  0],
#> [ 0,  0,  0,  0, -1, -1,  2, -1],
#> [ 0,  0,  0,  0,  0,  0, -1,  2]]
m := [ [  0,  0,  0,  1 ],
       [  0,  0,  1,  0 ],
       [  0,  1,  0,  0 ],
       [  1,  0,  0,  0 ] ];;
DecomposedMat( m );
#>[ [ 1, 4 ], [ 2, 3 ] ]
PrintArray( m{[ 1, 4 ]}{[ 1, 4 ]});
#>[[0, 1],
#> [1, 0]]
a := [ [ 1, 2 ], [ 3, 1 ] ];;
IsDiagonalMat( a );
#>false
IsLowerTriangularMat( a );
#>false
a[1][2] := 0;;
IsLowerTriangularMat( a );
#>true
a := [ [ 1, 2 ], [ 3, 1 ] ];;
l := [ [ 1, 0 ], [ 0, 1 ], [ 1, 1 ], [ 0, 0 ] ];;
l * a;
#>[ [ 1, 2 ], [ 3, 1 ], [ 4, 3 ], [ 0, 0 ] ]
IsNormalizing( l, a );
#>false
l := [ 1, , 2, , , 3 ];;
a := [ 1, 2, 3, 4 ];;
b := [ 2, 3, 1, 4 ];;
p := PermListList( a, b );
#>(1,2,3)
Permuted( b, p );
#>[ 1, 2, 3, 4 ]
PermListList( a, [ 3, 2, 1, 5 ] );
#>false
IntListToString( [ 4, 2, 2, 1, 1 ] );
#>"42211"
IntListToString( [ 14, 2, 2, 1, 1 ] );
#>"(14,2,2,1,1)"
IntListToString( [ 14, 2, 2, 1, 1 ], "{}" );
#>"{14,2,2,1,1}"
d := PartitionTuples( 3,2 );;
Print(Join(List(d,PartitionTupleToString),"   "),"\n");
#>111.   11.1   1.11   .111   21.   1.2   2.1   .21   3.   .3
W:=CoxeterGroupByCoxeterMatrix([[1,5],[5,1]]);
#>CoxeterGroupByCoxeterMatrix([ [ 1, 5 ], [ 5, 1 ] ])
PrintDiagram(W);
#>I2(5) 1 -5- 2
W:=ComplexReflectionGroup(24);;
SchurElements(Hecke(W,1));
#>[ 336, 336, 112, 112, 112, 112, 56, 56, 48, 48, 42, 42 ]
Print( "$Id: chevie.tst 2016/12/11 ",QuoInt(4000000000,time)," GAPstones\n");
