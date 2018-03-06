f:=GaloisField(2,2);
w:=f.root; z:=f.zero; o:=f.one;
g1:=[
[ z, o, z, z, z, z, z, z, z, z, z, z],
[ z, z, o, z, z, z, z, z, z, z, z, z],
[ o, z, z, z, z, z, z, z, z, z, z, z],
[ z, z, z, o, z, z, z, z, z, z, z, z],
[ z, z, z, z, o, z, z, z, z, z, z, z],
[ z, z, z, z, z, o, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, o, z, z, z, z],
[ z, z, z, z, z, z, z, z, o, z, z, z],
[ z, z, z, z, z, z, o, z, z, z, z, z],
[ z, z, z, z, z, z, z, z, z, z, o, z],
[ z, z, z, z, z, z, z, z, z, z, z, o],
[ z, z, z, z, z, z, z, z, z, o, z, z]
];
g2:=[
[ o, z, z, z, z, z, z, z, z, z, z, z],
[ z, o, z, z, z, z, z, z, z, z, z, z],
[ z, z, o, z, z, z, z, z, z, z, z, z],
[ z, z, z, z, o, z, z, z, z, z, z, z],
[ z, z, z, z, z, o, z, z, z, z, z, z],
[ z, z, z, o, z, z, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, o, z, z, z, z],
[ z, z, z, z, z, z, z, z, o, z, z, z],
[ z, z, z, z, z, z, o, z, z, z, z, z],
[ z, z, z, z, z, z, z, z, z, z, z, o],
[ z, z, z, z, z, z, z, z, z, o, z, z],
[ z, z, z, z, z, z, z, z, z, z, o, z]
];
g3:=[
[ z, o, z, z, z, z, z, z, z, z, z, z],
[ z, z, z, z, z, o, z, z, z, z, z, z],
[ z, z, o, z, z, z, z, z, z, z, z, z],
[ z, z, z, o, z, z, z, z, z, z, z, z],
[ z, z, z, z, o, z, z, z, z, z, z, z],
[ o, z, z, z, z, z, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, z, o, z, z, z],
[ z, z, z, z, z, z, z, z, z, z, z, o],
[ z, z, z, z, z, z, z, z, z, o, z, z],
[ z, z, z, z, z, z, o, z, z, z, z, z],
[ z, z, z, z, z, z, z, o, z, z, z, z],
[ z, z, z, z, z, z, z, z, z, z, o, z]
];
g4:=[
[ z, o, z, z, z, z, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, o, z, z, z, z],
[ z, z, o, z, z, z, z, z, z, z, z, z],
[ z, z, z, z, o, z, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, z, z, o, z, z],
[ z, z, z, z, z, z, z, z, z, z, o, z],
[ z, z, z, z, z, z, o, z, z, z, z, z],
[ o, z, z, z, z, z, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, z, o, z, z, z],
[ z, z, z, o, z, z, z, z, z, z, z, z],
[ z, z, z, z, z, z, z, z, z, z, z, o],
[ z, z, z, z, z, o, z, z, z, z, z, z]
];
g5:=[
[ o, o, o, z, z, z, z, z, z, z, z, z],
[ o, w, w^2, z, z, z, z, z, z, z, z, z],
[ o, w^2, w, z, z, z, z, z, z, z, z, z],
[ z, z, z, o, o, o, z, z, z, z, z, z],
[ z, z, z, o, w^2, w, z, z, z, z, z, z],
[ z, z, z, o, w, w^2, z, z, z, z, z, z],
[ z, z, z, z, z, z, o, o, o, z, z, z],
[ z, z, z, z, z, z, o, w^2, w, z, z, z],
[ z, z, z, z, z, z, o, w, w^2, z, z, z],
[ z, z, z, z, z, z, z, z, z, o, o, o],
[ z, z, z, z, z, z, z, z, z, o, w, w^2],
[ z, z, z, z, z, z, z, z, z, o, w^2, w]
];
g6:=[
[w, z, z, z, z, z, z, z, z, z, z, z],
[ z, w, z, z, z, z, z, z, z, z, z, z],
[ z, z, w, z, z, z, z, z, z, z, z, z],
[ z, z, z, w^2, z, z, z, z, z, z, z, z],
[ z, z, z, z, w^2, z, z, z, z, z, z, z],
[ z, z, z, z, z, w^2, z, z, z, z, z, z],
[ z, z, z, z, z, z, w, z, z, z, z, z],
[ z, z, z, z, z, z, z, w, z, z, z, z],
[ z, z, z, z, z, z, z, z, w, z, z, z],
[ z, z, z, z, z, z, z, z, z, w^2, z, z],
[ z, z, z, z, z, z, z, z, z, z, w^2, z],
[ z, z, z, z, z, z, z, z, z, z, z, w^2]
];
gens:=[g1,g2,g3,g4,g5,g6];
G := Group (gens, gens[1]^0); 
