# construct the group

C6 := CyclicGroup(AgWords, 6);
S3 := SymmetricGroup(AgWords, 3);
H := WreathProduct(C6,S3);
G := SpecialAgGroup(H / Centre(H));
G.name := "G";

# compute the automorphism group

A := AutGroupSagGroup(G);
Size(A);

# compute the size of the outer automorphism group

innsize := Size(G) / Size(Centre(G));
outsize := Size(A) / innsize;

# examine the structure

AutGroupStructure(A);

