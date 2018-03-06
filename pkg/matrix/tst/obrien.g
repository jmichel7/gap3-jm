#tests for matrix package 

RequirePackage ("matrix");
InfoTensor := Print;
InfoSmash := Print;
InfoPrim := Print;
InfoSemiLinear := Print;
InfoExtraSpecial := Print;

#not a tensor product 
ReadDataPkg( "matrix", "data", "j1d14.gap" );
x:=IsTensor(G);
Print ("Group tensor product? ", x[1], "\n");

#semilinear 
ReadDataPkg( "matrix", "data", "onand18.gap" );
x := IsSemiLinear (G);
Print ("Group semilinear? ", x[1], "\n");

#tensor product 
ReadDataPkg( "matrix", "data", "onand27.gap" );
x:=IsTensor(G);
Print ("Group tensor product? ", x[1], "\n");
Print ("Tensor decomposition is ", x[2], "\n");

#tensor product -- this examples requires Meataxe package 
#ReadDataPkg( "matrix", "data", "a6d16.gap" );
#x:=IsTensor(G);
#Print ("Group tensor product? ", x[1], "\n");
#Print ("Tensor decomposition is ", x[2], "\n");

#imprimitive 
ReadDataPkg( "matrix", "data", "sl25d6.gap" );
x:=IsPrimitive (G, [[6, 1]]);
Print ("Group primitive? ", x[1], "\n");
B := BlockSystemFlag (x[2]);
Print ("Block system is ", B, "\n");

#imprimitive; compute kernel
ReadDataPkg( "matrix", "data", "l217.gap" );
x:=IsPrimitive (G, [[18, 1]]);
Print ("Group primitive? ", x[1], "\n");
B := BlockSystemFlag (x[2]);
Print ("Permutation representation is  ", PermGroupFlag (B), "\n");
P:=PermGroupFlag (B);
maps:=MapsFlag (B); #what are images of generators of G? 
K:=ApproximateKernel (G, P, 100, 25, maps);; 

#SmashGModule decomposition  
ReadDataPkg( "matrix", "data", "esd4.gap" );
M := GModule (G);
S:=[G.4];
x:=SmashGModule (M, S);
Print ("Group has extraspecial decomposition? ", ExtraSpecialFlag (M), "\n");
P := ExtraSpecialGroupFlag (M);
Print ("G normalises the  p-group\n", P, "\n");

#symmetric tensor power decomposition 
MG:=GL(3,2);;
PG:=SymmetricGroup (3);;
G:=PowerWreathProduct(MG,PG);;
gens:=Generators (G);;
M := GModule (G);
S:=[gens[1]];;
x:=SmashGModule (M, S);
Print ("Group has symmetric tensor product decomposition? ", 
SymTensorProductFlag (M), "\n");
Print ("The permutations are ", SymTensorPermFlag (M), "\n");

#imprimitive decomposition -- also compute kernel 
MG:=GL(5,5);
PG:=SymmetricGroup (3);
G:= ImprimitiveWreathProduct(MG,PG);
S:=[G.1, G.2];
M:=GModule (G);
x:= SmashGModule (M, S);
Print ("Group has imprimitive decomposition? ", ImprimitiveFlag (M), "\n");
#now calculate kernel 
B:=BlockSystemFlag(M);;
P:=PermGroupFlag (B);
m:=MapsFlag (B); #what are images of generators of G? 
K:=ApproximateKernel (G, P, 100, 25, m);; 

#tensor product decomposition 
MG1:=GL(5,5);
MG2:=GL(3,5);
g1:=MG1.1; g2:=MG1.2; g3:=Identity(MG1);
h1:=MG2.1; h2:=MG2.2; h3:=Identity(MG2);
w:=KroneckerProduct(g1,h3);
x:=KroneckerProduct(g2,h3);
y:=KroneckerProduct(g3,h1);
z:=KroneckerProduct(g3,h2);
gens := [ w,x,y,z ];
G := Group (gens, gens[1]^0);
S:=[gens[1],gens[2]];
M:=GModule (G);
x:=SmashGModule (M, S);
Print ("Group has tensor decomposition? ", TensorProductFlag (M), "\n");
Print ("Tensor factors are ", TensorFactorsFlag (M), "\n");

Print ("Now the manual example\n");
#manual example 
# First set up the natural permutation module for the alternating group
# $A_5$ over the field $GF(2)$.
P := Group ((1,2,3), (3,4,5));;
M := PermGModule (P, GF(2));
# Now test for irreducibility, and calculate a proper submodule.
IsIrreducible (M);
SM := SubGModule (M, SubbasisFlag (M));;
DimensionFlag (SM);
DSM := DualGModule (SM);;
# We must prove irreducibility first.
# before we can test to see if SM is self-dual.
IsIrreducible (SM);
IsAbsolutelyIrreducible (SM);
IsomorphismGModule (SM, DSM);
# This is an explicit isomorphism.
# Now form a tensor product and decompose it into composition factors.
TM := TensorProductGModule (SM, SM);;
cf := CompositionFactors (TM);;
Length (cf);

#DimensionFlag(cf[1][1]); cf[1][2];
#DimensionFlag(cf[2][1]); cf[2][2];
#DimensionFlag(cf[3][1]); cf[3][2];
# This tells us that TM has three composition factors, of dimensions
# 1, 4 and 4, with multiplicities 4, 2 and 1, respectively.
# Is one of the 4-dimensional factors isomorphic to TM?
#IsomorphismGModule (SM, cf[2][1]);
#IsomorphismGModule (SM, cf[3][1]);
#IsAbsolutelyIrreducible (cf[2][1]);
#DegreeFieldExtFlag(cf[2][1]);
# If we extend the field of  cf[2][1] to $GF(4)$, it should
# become reducible.  
#MM := GModule (GeneratorsFlag (cf[2][1]), GF(4));;
#CF2 := CompositionFactors (MM);;
#Length (CF2);
#DimensionFlag (CF2[1][1]); CF2[1][2];
#DimensionFlag (CF2[2][1]); CF2[2][2];
# It reduces into two non-isomorphic 2-dimensional factors.

a := [
[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]] * Z(2)^0;;
b := [
[1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]] * Z(2)^0;;
c := [
[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]] * Z(2)^0;;
gens := [a, b, c];;
# Next we define the module.
M := GModule (gens);;
# So far only the basic components have been set.
RecFields (M);
# First we check for irreducibility and absolute irreducibility.
IsIrreducible (M);
IsAbsolutelyIrreducible (M);
# A few more components have been set during these two function calls.
RecFields(M);
# The function Commutators forms the list of commutators of generators.
S := Commutators(gens);;
InfoSmash := Print;;
# Setting InfoSmash to Print means that SmashGModule prints out
# intermediate output to tell us what it is doing. If we
# read this output it tells us what kind of decomposition SmashGModule
# has found. Otherwise the output is only a true or false.
# All the relevant information is also contained in the components
# of M which are set by this function.

SmashGModule (M, S);

# Additional components are set during the call to SmashGModule 
RecFields(M);
# We can access the components either directly or through flag functions.
SemiLinearFlag (M);
# This flag tells us G that acts semilinearly.

DegreeFieldExtFlag (M);
#This flag tells us the relevant extension field is GF(2\^3)
Length (LinearPartFlag (M));

# LinearPartFlag (M) is a set of normal subgroup generators for the 
# intersection of G with GL(4, GF(2\^3)). It is also the value of the set S 
# at the end of the call to SmashGModule is bigger than the set S which was
# input (which had 3 elements) since conjugates have been added.

FrobeniusAutomorphismsFlag (M);
# The first two generators of G act linearly, the last induces the field
# automorphism which maps x to x\^2 (= x\^(2\^1)) on GF(2\^3)

ReadDataPkg( "matrix", "data", "l217.gap" );
# Initialise a seed for random element generation.
InitPseudoRandom (G, 10, 100);;
# Now select a random element.
g := PseudoRandom (G);;
OrderMat (g);
h := ElementOfOrder (G, 8, 10);;
OrderMat (h);
#Is the group primitive?
R := IsPrimitive(G);;
#Examine the boolean returned.
R[1];
M := R[2];;
#What is the block system found?
B := BlockSystemFlag (M);
v := BlockFlag (B);
#Illustrate use of MinBlocks
B := MinBlocks (M, v);;
B;






