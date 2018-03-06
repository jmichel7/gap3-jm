# File:  xmod.t,  version 13/ 1/97

# contains all the commands in the xmod manual chapter

c5 := CyclicGroup(5);; c5.name :="c5";;
X1 := AutomorphismXMod( c5 );
Print( X1, "\n" );
XModPrint( X1 );

shom := GroupHomomorphismByImages( c5, c5, [(1,2,3,4,5)], [(1,5,4,3,2)] );
rhom := InclusionMorphism( X1.range, X1.range );
mor1 := XModMorphism( X1, X1, [ shom, rhom ] );
Print( mor1, "\n", IsXModMorphism( mor1 ), "\n" );
XModMorphismPrint( mor1 );

s3c4 := Group( (1,2), (2,3), (4,5,6,7) );
s3c4.name := "s3c4";
s3 := Subgroup( s3c4, [ (1,2), (2,3) ] );
s3.name := "s3";
gen := s3c4.generators;
imb := [ (1,2), (2,3), () ];
bX := GroupHomomorphismByImages( s3c4, s3, gen, imb );
im1 := List( gen, g -> g^(1,2) );
a12 := GroupHomomorphismByImages( s3c4, s3c4, gen, im1 );
im2 := List( gen, g -> g^(2,3) );
a23 := GroupHomomorphismByImages( s3c4, s3c4, gen, im2 );
autX := Group( a12, a23 );
aX := GroupHomomorphismByImages( s3, autX, [(1,2),(2,3)], [a12,a23] );
X := XMod( bX, aX );
Print( "\n", X, "\n", IsXMod( X ), "\n" );
XModPrint( X );

s4 := Group( (1,2,3,4), (1,2) );; s4.name := "s4";;
a4 := Subgroup( s4, [ (1,2,3), (2,3,4) ] );; a4.name :="a4";;
k4 := Subgroup( a4, [ (1,2)(3,4), (1,3)(2,4) ] );; k4.name := "k4";;
CX := ConjugationXMod( a4, k4 );
Print( CX, "\n" );
k4.name := "v4";
XModName( CX );
Print( CX.name, "\n" );

d8 := Subgroup( s4, [ (1,2,3,4), (1,3) ] );  d8.name := "d8";
gend8 := d8.generators;
genk4 := k4.generators;
f := GroupHomomorphismByImages( d8, k4, gend8, genk4 );
EX := CentralExtensionXMod( f );
Print( EX, "\n" );
XModPrint( EX );

q8 := Group( (1,2,3,4)(5,8,7,6), (1,5,3,7)(2,6,4,8) );
q8.name := "q8";  genq8 := q8.generators;
iaq8 := InnerAutomorphismGroup( q8 );
aq8 := GroupHomomorphismByImages( q8, q8, genq8,
             [(1,5,3,7)(2,6,4,8),(1,4,3,2)(5,6,7,8)]);
genA := Concatenation( iaq8.generators, [aq8] );
Print( genA, "\n" );
id := IdentityMapping( q8 );
A := Group( genA, id );
AX := AutomorphismXMod( q8, A );
Print( AX , "\n" );
Print( RecFields( AX ), "\n" );
IX := InnerAutomorphismXMod( q8 );
Print( IX, "\n" );

imf := [ (1,3)(2,4), (1,3)(2,4) ];
f := GroupHomomorphismByImages( k4, d8, genk4, imf );
TX := TrivialActionXMod( f );
Print( TX, "\n" );
XModPrint( TX );

k4gen := k4.generators;
k4im := [ (1,3)(2,4), (1,4)(2,3) ];
ak4 := GroupHomomorphismByImages( k4, k4, k4gen, k4im );
Ak4 := Group( ak4 );
R := rec( );;
R.module := k4;;
R.auto := Ak4;;
Print( IsRModule( R ), "\n", RecFields( R ), "\n", R.perm, "\n" );
RX := RModuleXMod( R );
Print( RX, "\n" );
XModPrint( RX );

SX := XModSelect( 24, 14, "conj", 3 );
Print( SX, "\n" );
XModPrint( SX );

Print( RecFields( XModOps ), "\n" );
XModOps.Print( CX );
Print( "\n", Size( CX ), "\n", Elements( CX ), "\n" );
Print( IsConjugation( CX ), "\n" );
Print( IsAspherical( CX ), "\n" );
Print( IsSimplyConnected( EX ), "\n" );
Print( IsCentralExtension( EX ), "\n" );
Print( IsAutomorphismXMod( AX ), "\n" );
Print( IsTrivialAction( TX ), "\n" );
Print( IsZeroBoundary( EX ), "\n" );
Print( IsRModule( RX ), "\n" );
Print( WhatTypeXMod( EX ), "\n" );

DX := XModOps.DirectProduct( CX, CX );
Print( DX, "\n" );
XModPrint( DX );

smor := GroupHomomorphismByImages( q8, k4, genq8, genk4 );
Print( smor, "\n", IsHomomorphism( smor ), "\n" );
sl23 := SX.range;
gensl23 := sl23.generators;
images := [ (1,2)(3,4), (1,3)(2,4), (2,3,4) ];
Print( gensl23, "\n", images, "\n" );
rmor := GroupHomomorphismByImages( sl23, a4, gensl23, images );
Print( IsHomomorphism( rmor ), "\n" );
mor := XModMorphism( SX, CX, [ smor, rmor ] );
Print( mor, "\n" );

XModPrintLevel := 3;
Print( IsXModMorphism( mor ), "\n" );
XModPrintLevel := 1;

XModMorphismPrint( mor );
k4.name := "k4"; XModName( CX ); XModMorphismName( mor );
Print( mor.name, "\n" );
Print( IsMonomorphism( mor ), "\n" );
Print( IsEpimorphism( mor ), "\n" );
Print( IsIsomorphism( mor ), "\n" );
Print( IsEndomorphism( mor ), "\n" );
Print( IsAutomorphism( mor ), "\n" );

Print( IdentitySubXMod( CX ), "\n" );
q8 := SX.source;  genq8 := q8.generators;
q8.name := "q8";  XModName( SX );
c4 := Subgroup( q8, [ genq8[1] ] );
Print( c4, "\n" );
c4.name := "c4";
subSX := SubXMod( SX, c4, q8 );
Print( c4, "\n", subSX, "\n" );
XModPrint( subSX );
Print( IsSubXMod( SX, subSX ), "\n" );
inc := InclusionMorphism( subSX, SX );
Print( inc, "\n", IsXModMorphism( inc ), "\n" );
XModMorphismPrint( inc );

Print( IsNormalSubXMod( SX, subSX ), "\n" );
NSX := NormalSubXMods( SX );
Print( NSX, "\n" );
Print( Size( NSX[3] ), "\n" );
FX := FactorXMod( SX, NSX[3] );
Print( FX, "\n", Size( FX ), "\n" );
KX := Kernel( mor );
Print( KX, "\n" );
XModPrint( KX );
Print( IsNormalSubXMod( SX, KX ), "\n" );
JX := ImageXModMorphism( mor, subSX );
Print( subSX, "\n", JX, "\n" );
Print( RecFields( mor ), "\n" );
XModPrint( JX );

g := Elements( q8 )[8];
Print( g, "\n" );
psi := XModOps.InnerAutomorphism( subSX, g );
Print( psi, "\n" );
XModMorphismPrint( psi );
Print( XModMorphismOps.Order( psi ), "\n" );
xcomp := XModMorphismOps.CompositeMorphism( psi, inc );
Print( xcomp, "\n" );
XModMorphismPrint( xcomp );

c2 := Subgroup( q8, [ genq8[1]^2 ] );
Print( c2, "\n" );
c2.name := "c2";
sub2 := SubXMod( subSX, c2, q8 );
Print( sub2, "\n" );
inc2 := InclusionMorphism( sub2, subSX );
Print(inc2, "\n" );
PX := SourceXModXPModMorphism( inc2 );
Print( PX, "\n", IsConjugation( PX ), "\n" );

Print( "\nEnd of section: CROSSED MODULES \n" );

#######################################################################

s3c4gen := s3c4.generators;
t1 := GroupHomomorphismByImages( s3c4, s3, s3c4gen, [(1,2),(2,3),()] );
C1 := Cat1( s3c4, t1, t1 );
Print( s3c4gen, "\n", C1, "\n" );
Cat1Print( C1 );
CX1 := Cat1XMod( X1 );
Print( X1, "\n", CX1, "\n", CX1.source.generators, "\n" );
XC1 := XModCat1( C1 );
Print( XC1, "\n", WhatTypeXMod( XC1 ), "\n\n" );

h20 := Group( (1,2,3,4,5), (2,3,5,4) );
h20.name := "h20";
genh20 := h20.generators;
imh20 := [ (), (2,3,5,4) ];
h := GroupHomomorphismByImages( h20, h20, genh20, imh20 );
t := h;
C := Cat1( h20, t, h );
Print( C, "\n", IsCat1( C ), "\n" );
Cat1Print( C );
C.range.name := "c4";
Cat1Name( C );
Print( C.name, "\n" );

c5 := Subgroup( h20, [ (1,2,3,4,5) ] );
c5.name := "c5";
CC := ConjugationCat1( h20, c5 );
Print( CC, "\n" );
Cat1Print( CC );
ct := CC.tail;
ch := CC.head;
CG := CC.source;
genCG := CG.generators;
x := genCG[2] * genCG[3];
tx := Image( ct,x );
hx := Image( ch, x );
Print( x, "\n", tx, "\n", hx, "\n" );
Print( RecFields( CC ), "\n" );

Print( RecFields( Cat1Ops ), "\n" );
Print( Size( C ), "\n", Elements( C ), "\n" );
XC := XModCat1( C );
Print( XC, "\n" );
XModPrint( XC );

CX := ConjugationXMod( a4, k4 );
CCX := Cat1XMod( CX );
Print( CX, "\n", CCX, "\n" );
Cat1Print( CCX );
Unbind( CX.cat1 );
SCX := SemidirectCat1XMod( CX );
Print( SCX, "\n" );
Cat1Print( SCX );

Print( Length( Cat1List ), "\n", Cat1List[8], "\n" );
Print( Cat1ListMaxSize, "\n" );
Print( NumbersOfIsomorphismClasses[18], "\n" );
L := Cat1Select( 18 );
Print( L, "\n" );
SC := Cat1Select( 18, 5 );
Print( SC, "\n" );
SC := Cat1Select( 18, 5, 4 );
Print( SC, "\n" );
Cat1Print( SC );
XSC := XModCat1( SC );
Print( XSC, "\n" );
AC := Cat1Select( 12, 5, 1 );
Print( AC, "\n" );

GCCX := CCX.source;
GAC := AC.source;
genGAC := GAC.generators;
im := Sublist( GCCX.generators, [1..2] );
Print( GCCX, "\n", GAC, "\n", genGAC, "\n", im, "\n" );
musrc := GroupHomomorphismByImages( GAC, GCCX, genGAC, im );
murng := InclusionMorphism( a4, a4 );
mu := Cat1Morphism( AC, CCX, [ musrc, murng ] );
Print( mu, "\n", IsCat1Morphism( mu ), "\n" );
CCX.source.name := "a4.k4";
Cat1Name( CCX );
Cat1MorphismName( mu );
Print( CCX.name, "\n", mu.name, "\n" );
Cat1MorphismPrint( mu );

Print( IsMonomorphism( mu ), "\n" );
Print( IsEpimorphism( mu ), "\n" );
Print( IsIsomorphism( mu ), "\n" );
Print( IsEndomorphism( mu ), "\n" );
Print( IsAutomorphism( mu ), "\n" );

GSC := SC.source;
homsrc := GroupHomomorphismByImages( a4, GSC,
         [(1,2,3),(2,3,4)], [(4,5,6),(4,6,5)] );
musrc := Cat1MorphismSourceHomomorphism( AC, SC, homsrc );
Print( musrc, "\n", IsCat1Morphism( musrc ), "\n" );
Cat1MorphismPrint( musrc );

revCC := ReverseCat1( CC );
Print( revCC, "\n" );
revmu := ReverseIsomorphismCat1( CC );
Print( revmu, "\n", IsCat1Morphism( revmu ), "\n" );

CX.cat1 := CCX;
CSX := Cat1XMod( SX );
Print( CSX, "\n", mor, "\n" );
catmor := Cat1MorphismXModMorphism( mor );
Print( catmor, "\n", IsCat1Morphism( catmor ), "\n" );
Cat1MorphismPrint( catmor );
xmu := XModMorphismCat1Morphism( mu );
Print( mu, "\n", xmu, "\n" );

Print( psi, "\n", inc, "\n" );
mupsi := Cat1MorphismXModMorphism( psi );
muinc := Cat1MorphismXModMorphism( inc );
Print( mupsi, "\n", muinc, "\n" );
mucomp := Cat1MorphismOps.CompositeMorphism( mupsi, muinc );
muxcomp := Cat1MorphismXModMorphism( xcomp );
Print( mucomp, "\n", muxcomp, "\n" );
Print( "mucomp = muxcomp ? ", mucomp = muxcomp, "\n" );

Print( IdentitySubCat1( SC ), "\n" );
d20 := Subgroup( h20, [ (1,2,3,4,5), (2,5)(3,4) ] );
subC := SubCat1( C, d20 );
Print( subC, "\n" );
Cat1Print( subC );
Print( InclusionMorphism( subC, C ), "\n" );
Print( NormalSubCat1s( SC ), "\n" );

all := AllCat1s( a4 );
Print( all, "\n" );

Print( "\nEnd of section: CAT1-GROUPS \n" );

#########################################################################

Print( X1, "\n" );
chi1 := XModDerivationByImages( X1, [ () ] );
Print( chi1, "\n", IsDerivation( chi1 ), "\n", RecFields( chi1 ), "\n" );
xi1 := SectionDerivation( chi1 );
Print( xi1, "\n", xi1.cat1, "\n\n" );

Print( RegularDerivations( X1 ), "\n", AllDerivations( X1 ), "\n" );
Print( DerivationsSorted( X1 ), "\n" );
imder1 := X1.derivations.genimageList;
Print( imder1, "\n" );

CX1 := Cat1XMod( X1 );
Print( CX1, "\n");
CX1.source.name := "Hol(c5)";
Cat1Name( CX1 ); 
Print( RegularSections( CX1 ), "\n" );
Print( CX1.sections.genimageList, "\n\n" );

chi2 := XModDerivationByImages( X1, imder1[2] );
Print( chi2, "\n", DerivationImage( chi2, (1,4)(2,3) ), "\n" );
Print( DerivationImages( chi2 ), "\n\n" );
PrintList( DerivationTable( X1) );
PrintList( WhiteheadGroupTable( X1 ) );

sigma2 := SourceEndomorphismDerivation( chi2 );
rho2 := RangeEndomorphismDerivation( chi2 );
xi2 := SectionDerivation( chi2 );;
gamma2 := SourceEndomorphismSection( xi2 );
Print( sigma2, "\n", rho2, "\n", gamma2, "\n" );
mor2 := XModMorphism( X1, X1, [sigma2,rho2] );
mu2 := Cat1Morphism( CX1, CX1, [gamma2,rho2] );
Print( mor2, "\n", mu2, "\n" );

Print( "\nEnd of section: About derivations and sections.\n\n" );

Print( XSC, "\n" );
imchi := [ (1,2,3)(4,6,5), (1,2,3)(4,6,5) ];;
chi := XModDerivationByImages( XSC, imchi );
Print( chi, "\n" );
im0 := [ (1,3,2)(4,5,6), () ];;
Print( IsDerivation( XSC, im0 ), "\n" );
Print( DerivationImage( chi, (4,6,5) ), "\n" );
Print( XSC.source.elements, "\n" );
Print( DerivationImages( chi ), "\n" );
Print( InnerDerivation( XSC, (1,2,3)(4,6,5) ), "\n" );
PrintList( ListInnerDerivations( XSC ) );
Print( RecFields( chi.operations ), "\n" );

Print( "\n", SC, "\n" );
imxi := [ (1,2,3), (1,2)(4,6) ];
xi := Cat1SectionByImages( SC, imxi );
Print( xi, "\n" );
im0 := [ (1,2,3), (2,3)(4,5)];;
Print( IsSection( SC, im0 ), "\n" ); 
Print( IsRegular( chi ), "\n", IsRegular( xi ), "\n\n" );
Print( RecFields( xi.operations ), "\n" );

regXSC := RegularDerivations( XSC );
Print( regXSC, "\n" );
PrintList( regXSC.genimageList );
Print( RecFields( regXSC ), "\n" );
allXSC := AllDerivations( XSC );
Print( allXSC, "\n", DerivationsSorted( allXSC ), "\n" );
PrintList( allXSC.genimageList );
PrintList( DerivationTable( allXSC ) );
Print( AreDerivations( regXSC ), "\n" );

Unbind( XSC.derivations );
regSC := RegularSections( SC );
Print( regSC, "\n" );
allSC := AllSections( SC );
Print( allSC, "\n", RecFields( allSC ), "\n" );
PrintList( allSC.genimageList );
allXSC := AllDerivations( XSC, "cat1" );
Print( allXSC, "\n\n" );
Print( AreSections( allSC ), "\n" );

chi8 := XModDerivationByImages( XSC, allXSC.genimageList[8] );
Print( chi8, "\n" );
xi8 := SectionDerivation( chi8 );
Print( xi8, "\n\n" );
xi4 := Cat1SectionByImages( SC, allSC.genimageList[4] );
Print( xi4, "\n" );
chi4 := DerivationSection( xi4 );
Print( chi4, "\n" );

chi48 := CompositeDerivation( chi4, chi8 );
xi48 := CompositeSection( xi4, xi8 );
Print( chi48, "\n", xi48, "\n" );
Print( ( SectionDerivation( chi48 ) = xi48 ), "\n" );
WGT := WhiteheadGroupTable( XSC );
PrintList( WGT );
WMT := WhiteheadMonoidTable( XSC );
PrintList( WMT );

inv4 := InverseDerivations( chi4 );
inv8 := InverseDerivations( chi8 );
Print( inv4, "\n", inv8, "\n" );
inv := ListInverseDerivations( XSC );
Print( "inv = ", inv, "\n" );

sigma8 := SourceEndomorphismDerivation( chi8 );
sigma4 := SourceEndomorphismDerivation( chi4 );
Print( sigma8, "\n", sigma4, "\n" );
TSE := TableSourceEndomorphismDerivations( XSC );
PrintList( TSE );
rho8 := RangeEndomorphismDerivation( chi8 );
rho4 := RangeEndomorphismDerivation( chi4 );
Print( rho8, "\n", rho4, "\n" );
TRE := TableRangeEndomorphismDerivations( XSC );
PrintList( TRE );
phi4 := XModEndomorphismDerivation( chi4 );
Print( phi4, "\n" );
gamma4 := SourceEndomorphismSection( xi4 );
rho4 := RangeEndomorphismSection( xi4 );
psi4 := Cat1EndomorphismSection( xi4 );
Print( gamma4, "\n", rho4, "\n", psi4, "\n" );

Print( "End of section: DERIVATIONS & SECTIONS \n" );

#######################################################################

Print( X1, "\n" );
WGX1 := WhiteheadPermGroup( X1 );
Print( WGX1, "\n", WGX1.generators, "\n" );
AX1 := AutomorphismPermGroup( X1 );
Print( AX1, "\n", AX1.generators, "\n" );
Print( XModMorphismAutoPerm( X1, AX1.generators[1] ), "\n\n" );
WX1 := Whitehead( X1 );
Print( WX1,"\n" );
NX1 := Norrie( X1 );
Print( NX1,"\n" );
LX1 := Lue( X1 );
Print( LX1, "\n" );
ActX1 := Actor( X1 );
XModPrint( ActX1 );
InActX1 := InnerActor( X1 );
Print( InActX1, "\n", InActX1 = ActX1, "\n" );
Print( InnerMorphism( X1 ), "\n", Centre( X1 ), "\n" );

Print( ActorSquareRecord( X1 ), "\n" );
WG := WhiteheadPermGroup( XSC );
Print( WG, "\n", XSC.derivations.genpos, "\n" );
Print( Elements( WG ), "\n" );

WXSC := Whitehead( XSC );
Print( WXSC, "\n" );
XModPrint( WXSC );

autXSC := AutomorphismPermGroup( XSC );
Print( autXSC, "\n" );
Print( autXSC.projsrc, "\n", autXSC.projrng, "\n" );
Print( autXSC.embedSourceAuto, "\n", autXSC.embedRangeAuto, "\n" );
Print( autXSC.autogens, "\n" );
Print( XModMorphismAutoPerm( XSC, (1,2)(3,4)(6,7) ), "\n" );

chi8im := ImageAutomorphismDerivation( phi4, chi8 );
Print( chi8im, "\n" );
Print( Position( allXSC.genimageList, chi8im.genimages ), "\n" );

NXSC := Norrie( XSC );
Print( NXSC, "\n" );
XModPrint( NXSC );
LXSC := Lue( XSC );
Print( LXSC, "\n" );
XModPrint( LXSC );
ActXSC := Actor( XSC );
Print( ActXSC, "\n" );
XModPrint( ActXSC );

innXSC := InnerMorphism( XSC );
Print( innXSC, "\n" );
XModMorphismPrint( innXSC );
ZXSC := Centre( XSC );
Print( ZXSC, "\n" );
InnActXSC := InnerActor( XSC );
Print( InnActXSC, "\n" );
XModPrint( InnActXSC );

ActSC := Actor( SC );
Cat1Print( ActSC );

Print( "\nEnd of section: Actor Squares \n" );

#######################################################################

d16 := DihedralGroup( 16 );
Print( d16, "\n" );  d16.name := "d16";
d8 := Subgroup( d16, [ (1,3,5,7)(2,4,6,8), (1,3)(4,8)(5,7) ] );
d8.name := "d8";
c4 := Subgroup( d8, [ (1,3,5,7)(2,4,6,8) ] );
c4.name := "c4";
DX := ConjugationXMod( d8, c4 );
iota := InclusionMorphism( d8,d16 );
Print( DX, "\n" );
IDXincl := InducedXMod( DX, iota );
XModPrint( IDXincl ); 

d8gen := d8.generators;
k4gen := k4.generators;
Print( d8gen, "\n", k4gen, "\n", DX, "\n" );
iota := GroupHomomorphismByImages( d8, k4, d8gen, k4gen );
IDXsurj := InducedXMod( DX, iota );
XModPrint( IDXsurj );

s3 := Subgroup( s4, [ (2,3), (1,2,3) ] );
c3 := Subgroup( s3, [ (1,2,3) ] );
s3.name := "s3"; c3.name := "c3";
InducedXMod( s4, s3, c3 );
Print( last, "\n" );

Print( AllInducedXMods( d8 ), "\n" );

CDX := Cat1XMod( DX );
Print( CDX, "\n" );
inc := InclusionMorphism( d8, d16 );
Print( inc, "\n" );
ICDX := InducedCat1( CDX, inc );
XICDX := XModCat1( ICDX );
Print( ICDX, "\n", XICDX, "\n", AbelianInvariants( XICDX.source ), "\n" );

Print( "\nEnd of section: INDUCED CONSTRUCTIONS \n" );

##########################################################################

incs3 := InclusionMorphism( s3, s3 );
Print( incs3, "\n", incs3.genimages, "\n" );
end8 := EndomorphismClasses( d8 );;
Print( RecFields( end8 ), "\n" );
Print( Length( end8.classes ), "\n" );
Print( end8.classes[8], "\n" );
innd8 := InnerAutomorphismGroup( d8 );
Print( innd8, "\n", innd8.generators, "\n" );
Print( IsAutomorphismGroup( innd8 ), "\n" );
Apair := AutomorphismPair( innd8 );
Print( Apair, "\n", IsAutomorphismPair( Apair ), "\n" );

L := [ [1,4], [1,2], [2,3], [1,3], [5] ];
Print( L, "\n", DistinctRepresentatives( L ), "\n" );
M := [ [2,5], [3,5], [4,5], [1,2,3], [1,2,3] ];
Print( M, "\n", CommonRepresentatives( L, M ), "\n" );

s4 := Group( (1,2,3,4), (1,2) );; s4.name := "s4";;
a4 := Subgroup( s4, [ (1,2,3), (2,3,4) ] );; a4.name :="a4";;
inc := InclusionMorphism( a4, s4 );
Print( inc, "\n" );
zero := ZeroMorphism( s4, a4 );
Print( zero, "\n" );

Ea4 := EndomorphismClasses( a4, 7 );
Ea4 := EndomorphismClasses( a4 );
Print( Ea4, "\n" );
Print( EndomorphismImages( a4 ), "\n" );
IdempotentImages( a4, 7 );
Print( IdempotentImages( a4, 2 ), "\n" );
Print( IdempotentImages( a4, 3 ), "\n" );

inna4 := InnerAutomorphismGroup( a4 );
Print( inna4, "\n", inna4.generators, "\n" );
Print( IsAutomorphismGroup( inna4 ), "\n" );

c3 := Subgroup( a4, [(1,2,3)] );;  c3.name := "c3";;
ac3 := AutomorphismGroup( c3 );
Print( ac3, "\n" );
pairc3 := AutomorphismPair( ac3 );
Print( pairc3, "\n", IsAutomorphismPair( pairc3 ), "\n" );
pc3 := pairc3.perm;
Print( pc3, "\n" );

P := AutomorphismPermGroup( a4 );
Print( P, "\n", P.generators, "\n" );

f := FreeGroup( 2 );
rels := [ f.1^3, f.2^3, (f.1*f.2)^2 ];
g := f/rels;
pairg := FpPair( g );
Print( pairg, "\n" );
Print( h20.generators, "\n" );
pairh := FpPair( h20 );
Print( pairh, "\n", pairh.fp.relators, "\n", IsFpPair( pairh ), "\n" );


agen := ac3.generators;
pgen := pc3.generators;
a := GroupHomomorphismByImages( pc3, ac3, pgen, agen );
Print( a, "\n" );
G := SemidirectProduct( pc3, a, c3 );
Print( G, "\n" );
G.name := "G";;
PG := SemidirectPair( G );
Print( PG, "\n", IsSemidirectPair( PG ), "\n" );

J := [ [1,2,3], [3,4], [3,4], [1,2,4] ];;
PrintList( J );
Print( DistinctRepresentatives( J ), "\n" );
K := [ [3,4], [1,2], [2,3], [2,3,4] ];
Print( K, "\n", CommonRepresentatives( J, K ), "\n" );
T := CommonTransversal( a4, c3 );
Print( T, "\n", IsCommonTransversal( a4, c3, T ), "\n" );

Print( "\nEnd of section: UTILITIES \n" );
