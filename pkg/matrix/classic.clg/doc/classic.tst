ResetRandom := function()
  R_X :=
  [ 264044297, 183552375, 186557741, 46008580, 91966775, 40228334, 111587641, 
  262452744, 44854923, 39845682, 206581571, 92903376, 231314873, 66408590, 
  258714463, 184978940, 247212781, 95976059, 76757457, 166475183, 209286133, 
  155943299, 70960569, 13942903, 202135272, 242522493, 239638894, 42966699, 
  154546372, 38965769, 121766026, 134939415, 151605824, 75410773, 197550847, 
  138731753, 64831979, 175307525, 190299479, 107235201, 53348259, 103934810, 
  207032077, 145748420, 124215055, 144100830, 6088801, 180843272, 27153689, 
  259704742, 66996111, 23305316, 136497461, 203183570, 153356235 ];
  R_N := 16;
end;;

#############################################################################
##
##  GeneralOrthogonalGroup, example 1
##
ResetRandom(); Print( "#T  GeneralOrthogonalGroup, example 1\n" );

g := GeneralOrthogonalGroup(0,5,3);
#>O(0,5,3)
Size( g );
#>103680
Size( SP(4,3) ); 
#>51840
DeterminantMat(g.1);
#>Z(3)^0
DeterminantMat(g.2);
#>Z(3)
DisplayMat( g.symmetricForm );    
#> . 1 . . .
#> 1 . . . .
#> . . 2 . .
#> . . . 2 .
#> . . . . 2
DisplayMat( g.quadraticForm );
#> . 1 . . .
#> . . . . .
#> . . 1 . .
#> . . . 1 .
#> . . . . 1
v1 := [1,2,0,1,2] * Z(3);     
#>[ Z(3), Z(3)^0, 0*Z(3), Z(3), Z(3)^0 ]
v1 * g.quadraticForm * v1;
#>Z(3)^0
v1 * g.symmetricForm * v1;
#>Z(3)


#############################################################################
##
##  SpinorNorm, example 1
##
ResetRandom(); Print( "#T  SpinorNorm, example 1\n" );

z  := GF(9).root;;
m1 := [[0,1,0,0,0,0,0,0,0],[1,2,2,0,0,0,0,0,0],
 [0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0],
  [0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],
  [0,2,1,0,0,0,0,0,0]]*z^0;;
m2 := [[z,0,0,0,0,0,0,0,0],[0,z^7,0,0,0,0,0,0,0],
  [0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],
  [0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],
  [0,0,0,0,0,0,0,0,1]]*z^0;;
form := IdentityMat( 9, GF(9) );;
form{[1,2]}{[1,2]} := [[0,2],[2,0]] * z^0;;
m1 * form * TransposedMat(m1) = form;
#>true
m2 * form * TransposedMat(m2) = form; 
#>true
SpinorNorm( form, m1 );
#>Z(3)^0
SpinorNorm( form, m2 );
#>Z(3^2)^5


#############################################################################
##
##  OrderMat, example 1
##
ResetRandom(); Print( "#T  OrderMat, example 1\n" );

OrderMat( [ [ Z(17)^4, Z(17)^12, Z(17)^4, Z(17)^7 ], 
   [ Z(17)^10, Z(17), Z(17)^11, 0*Z(17) ], 
   [ Z(17)^8, Z(17)^13, Z(17)^0, Z(17)^14 ], 
   [ Z(17)^14, Z(17)^10, Z(17), Z(17)^10 ] ] );
#>5220


#############################################################################
##
##  RecognizeSL, example 1
##
ResetRandom(); Print( "#T  RecognizeSL, example 1\n" );

m1 := [ [ 0*Z(17), Z(17), Z(17)^10, Z(17)^12, Z(17)^2 ], 
      [ Z(17)^13, Z(17)^10, Z(17)^15, Z(17)^8, Z(17)^0 ], 
      [ Z(17)^10, Z(17)^6, Z(17)^9, Z(17)^8, Z(17)^10 ], 
      [ Z(17)^13, Z(17)^5, Z(17)^0, Z(17)^12, Z(17)^5 ], 
      [ Z(17)^14, Z(17)^13, Z(17)^5, Z(17)^10, Z(17)^0 ] ];;
m2 := [ [ 0*Z(17), Z(17)^10, Z(17)^2, 0*Z(17), Z(17)^10 ], 
      [ 0*Z(17), Z(17)^6, Z(17)^0, Z(17)^4, Z(17)^15 ], 
      [ Z(17)^7, Z(17)^6, Z(17)^10, Z(17), Z(17)^2 ], 
      [ Z(17)^3, Z(17)^10, Z(17)^5, Z(17)^4, Z(17)^6 ], 
      [ Z(17)^0, Z(17)^8, Z(17)^0, Z(17)^5, Z(17) ] ];;
g := Group( m1, m2 );;

sl := RecognizeSL( g, 1 );
#>#I  field: 17, dimension: 5, number of generators: 2
#>#I  <G> could be almost simple
#>#I  <G> could be reducible
#>#I  <G> could still contain SL( 5, 17 )
#><< SL recognition record >>

sl.containsSL;
#>false

sl := RecognizeSL( sl, 20 );
#>#I  field: 17, dimension: 5, number of generators: 2
#>#I  <G> is SL( 5, 17 )
#><< SL recognition record >>

sl.containsSL;
#>true


#############################################################################
##
##  RecognizeSL, example 2
##
ResetRandom(); Print( "#T  RecognizeSL, example 2\n" );

sl := RecognizeSL( SL(10,4), 0 );
#>#I  field: 4, dimension: 10, number of generators: 2
#>#I  <G> could be almost simple
#>#I  <G> could be an almost sporadic group
#>#I  <G> could be an alternating group
#>#I  <G> could be a classical group
#>#I  <G> could be definable over a larger field
#>#I  <G> could be definable over GF(2)
#>#I  <G> could be imprimitive
#>#I  <G> could be reducible
#>#I  <G> could be a tensor product
#>#I  <G> could still contain SL( 10, 4 )
#><< SL recognition record >>

SetPrintLevel(sl,2); sl;         
#>#I  field: 4, dimension: 10, number of generators: 2
#>#I  <G> could be almost simple: A_1(2^2) A_1(2^3) A_1(2^4) A_1(2^5) A_1(2^
#>6) A_1(2^7) A_1(2^8) A_1(2^9) A_1(2^10) A_1(2^11) A_1(2^12) A_1(2^13) A_1(2^
#>14) A_1(2^15) A_1(2^16) A_1(2^17) A_1(2^18) A_1(2^19) A_1(2^20) A_2(2^1) A_2(
#>2^2) A_2(2^3) A_2(2^4) A_2(2^5) A_2(2^6) A_2(2^7) A_3(2^1) A_3(2^2) A_3(2^
#>3) A_3(2^4) A_4(2^1) A_4(2^2) A_5(2^1) A_6(2^1) A_2(3^1) A_1(3^2) A_1(5^1) A_
#>1(7^1) A_1(11^1) A_1(13^1) A_1(17^1) A_1(19^1) C_2(2^2) C_2(2^3) C_2(2^4) C_2(
#>2^5) C_2(2^6) C_3(2^1) C_3(2^2) C_4(2^1) C_5(2^1) C_2(3^1) D_4(2^1) D_4(2^
#>2) D_5(2^1) G_2(2^1) G_2(2^2) G_2(2^3) G_2(2^4) F_4(2^1) 2A_2(2^2) 2A_2(2^
#>3) 2A_2(2^4) 2A_2(2^5) 2A_2(2^6) 2A_2(2^7) 2A_3(2^2) 2A_3(2^3) 2A_3(2^4) 2A_4(
#>2^2) 2A_2(3^1) 2A_3(3^1) 2B_2(2^3) 2B_2(2^5) 2B_2(2^7) 2B_2(2^9) 2B_2(2^
#>11) 2F_4(2^1) 2D_4(2^1) 2D_4(2^2) 2D_5(2^1) 3D_4(2^1) 3D_4(2^2) 
#>#I  <G> could be an almost sporadic group: M11 M12 M22 J2 J3 
#>#I  <G> could be an alternating group: A5 A6 A7 A8 A9 A10 A11 A12 
#>#I  <G> could be a classical group: SU SP Omega 
#>#I  <G> could be definable over a larger field
#>#I  <G> could be definable over GF(2)
#>#I  <G> could be imprimitive
#>#I  <G> could be reducible: 1 .. 9
#>#I  <G> could be a tensor product
#>#I  <G> could still contain SL( 10, 4 )
#><< SL recognition record >>


#############################################################################
##
##  RecognizeSL, example 3
##
ResetRandom(); Print( "#T  RecognizeSL, example 3\n" );

m1 := [[0,9,10,10,5,9,10],[9,0,6,2,7,7,2],[10,6,1,4,7,3,1],
       [10,2,4,5,6,3,0],[5,7,7,6,4,10,6],[9,7,3,3,10,6,7],
       [10,2,1,0,6,7,0]]*Z(11);;
m2 := [[8,9,6,4,8,4,1],[5,2,5,0,7,7,4],[8,5,4,10,6,1,6],
       [1,7,7,3,6,1,10],[9,2,9,1,5,9,7],[7,10,0,2,1,7,8],
       [7,2,7,7,10,2,10]]*Z(11);;
g := Group( m1, m2 );;

sl := RecognizeSL( g, 20 );; SetPrintLevel( sl, 2 ); sl;
#>#I  field: 11, dimension: 7, number of generators: 2
#>#I  <G> could be almost simple: A_1(11^3) A_2(11^2) A_3(11^1) B_3(11^1) C_3(
#>11^1) G_2(11^1) 
#>#I  <G> could be an almost sporadic group: J1 
#>#I  <G> could be a classical group: Omega 
#>#I  <G> could be reducible: 0 1 3 4 6 
#>#I  <G> could still contain SL( 7, 11 )
#><< SL recognition record >>

so := RecognizeSO( sl, 0 );; SetPrintLevel( so, 2 );  so;
#>#W  Warning: group must act absolutely irreducible
#>#I  field: 11, dimension: 7, number of generators: 2
#>#I  symmetric form is known
#>#I  quadratic form is known
#>#I  <G> could be almost simple: A_1(11^3) A_2(11^2) A_3(11^1) G_2(11^1) 
#>#I  <G> could be an almost sporadic group: J1 
#>#I  <G> could still be an orthogonal group
#><< SO recognition record >>


#############################################################################
##
##  RecognizeSO, example 1
##
ResetRandom(); Print( "#T  RecognizeSO, example 1\n" );

z  := GF(9).root;;
m1 := [[0,1,0,0,0,0,0,0,0],[1,2,2,0,0,0,0,0,0],
  [0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0,0],
  [0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,1],
  [0,2,1,0,0,0,0,0,0]]*z^0;;
m2 := [[z,0,0,0,0,0,0,0,0],[0,z^7,0,0,0,0,0,0,0],
  [0,0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0,0],
  [0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1,0],
  [0,0,0,0,0,0,0,0,1]]*z^0;;
g := Group( m1, m2 );;

sl := RecognizeSL( g, 10 );; SetPrintLevel(sl,2); sl;
#>#I  field: 9, dimension: 9, number of generators: 2
#>#I  <G> could be a classical group: Omega 
#>#I  <G> could be reducible: 0 1 8 
#>#I  <G> could still contain SL( 9, 9 )
#><< SL recognition record >>

so := RecognizeSO( sl, 0 );
#>#W  Warning: group must act absolutely irreducible
#>#I  field: 9, dimension: 9, number of generators: 2
#>#I  symmetric form is known
#>#I  quadratic form is known
#>#I  <G> contains Omega0( 9, 9 )
#><< SO recognition record >>

DisplayMat( so.quadraticForm );
#> . 2 . . . . . . .
#> . . . . . . . . .
#> . . 2 . . . . . .
#> . . . 2 . . . . .
#> . . . . 2 . . . .
#> . . . . . 2 . . .
#> . . . . . . 2 . .
#> . . . . . . . 2 .
#> . . . . . . . . 2
DisplayMat( m1 * so.quadraticForm * TransposedMat(m1) );
#> . . . . . . . . .
#> 2 . . . . . . . 2
#> . . 2 . . . . . .
#> . . . 2 . . . . .
#> . . . . 2 . . . .
#> . . . . . 2 . . .
#> . . . . . . 2 . .
#> . . . . . . . 2 .
#> . 1 . . . . . . 2
DisplayMat( m2 * so.quadraticForm * TransposedMat(m2) );
#> . 2 . . . . . . .
#> . . . . . . . . .
#> . . 2 . . . . . .
#> . . . 2 . . . . .
#> . . . . 2 . . . .
#> . . . . . 2 . . .
#> . . . . . . 2 . .
#> . . . . . . . 2 .
#> . . . . . . . . 2

so.gScalars;
#>[ Z(3)^0, Z(3)^0 ]

List( [m1,m2], DeterminantMat );
#>[ Z(3), Z(3)^0 ]

List( [m1,m2], x -> SpinorNorm(so.symmetricForm,x) );
#>[ Z(3)^0, Z(3^2)^5 ]


#############################################################################
##
##  CRecognizeSL, example 2
##
ResetRandom(); Print( "#T  CRecognizeSL, example 2\n" );

m1 := [ [  6, 10,  6, 10,  7,  6 ],
        [  8,  6,  3,  7,  3,  2 ],
        [  7,  5,  0,  8,  3,  2 ],
        [  7,  6, 10,  4, 10,  5 ],
        [  6,  4,  1,  2,  8,  0 ],
        [  4,  5, 10, 10,  1,  6 ] ] * Z(11)^0;;
m2 := [ [  4,  3,  7,  0,  3,  2 ],
        [  3,  5,  4,  8,  7,  2 ],
        [  7,  2,  3,  6,  9,  0 ],
        [  4,  8,  8,  1,  1,  6 ],
        [  0, 10,  3,  6, 10, 10 ],
        [  7,  4,  6, 10,  2,  7 ] ] * Z(11)^0;;

g := Group( m1, m2 );;
sl := RecognizeSL( g, 20 );
#>#I  field: 11, dimension: 6, number of generators: 2
#>#I  <G> could be almost simple
#>#I  <G> could be definable over a larger field
#>#I  <G> could be imprimitive
#>#I  <G> could be reducible
#>#I  <G> could be a tensor product
#>#I  <G> could still contain SL( 6, 11 )
#><< SL recognition record >>

split := SplitMatGroup(g);;
quo := Group( split.quotient, split.quotient[1]^0 );;
RecognizeSL(quo,10);
#>#I  field: 11, dimension: 3, number of generators: 2
#>#I  <G> is SL( 3, 11 )
#><< SL recognition record >>

csl := CRecognizeSL( quo, quo.generators );
#>#I  <G> is SL( 3, 11 )
#><< constructive SL recognition record >>
w := quo.1 * quo.2 * quo.1;;
t := CRecSL.Rewrite( csl, w );;
c := Value( t, [m1,m2] );;
k1 := (m1*m2*m1/c) ^ split.basis;;
DisplayMat(k1);
#>  1  .  .  2  4  3
#>  .  1  .  .  2  9
#>  .  .  1  3 10 10
#>  .  .  .  3  8  7
#>  .  .  .  3  2  4
#>  .  .  .  9  5  9


w := quo.2 * quo.1;;
c := Value( CRecSL.Rewrite(csl,w), [m1,m2] );;
k2 := (m2*m1/c) ^ split.basis;;
DisplayMat(k2);
#>  1  .  .  3  3  1
#>  .  1  .  1  4  8
#>  .  .  1  .  4  5
#>  .  .  .  4  5  3
#>  .  .  .  2  1  9
#>  .  .  .  2  .  1
sub := Group( k1{[4..6]}{[4..6]}, k2{[4..6]}{[4..6]} );;
csl := CRecognizeSL( sub, sub.generators );
#>#I  <G> is SL( 3, 11 )
#><< constructive SL recognition record >>

w := sub.1 * sub.2;;
c := Value( CRecSL.Rewrite(csl,w), [k1,k2] );;
DisplayMat( k1*k2/c );
#>  1  .  .  .  .  .
#>  .  1  .  .  .  .
#>  .  .  1  .  .  .
#>  .  .  .  1  .  .
#>  .  .  .  .  1  .
#>  .  .  .  .  .  1


#############################################################################
##
##  CRecognizeSP, example 1
##
ResetRandom(); Print( "#T  CRecognizeSP, example 1\n" );

g := SPwithForm( 6, 11 );
#>SP(6,11)
sp := CRecognizeSP( g, g.generators, g.symplecticForm );           
#>#I  split element is known
#>#I  upper right corner is known
#>#I  upper left GL is known
#>#I  lower left corner is known
#>#I  <G> is SP( 6, 11 )
#><< constructive SP recognition record >>
a := RecSL.Random(g);;
t := Rewrite( sp, a );;
a = Value( t, g.generators );
#>true
