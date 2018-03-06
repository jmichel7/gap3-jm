F:=GaloisField(2);;
o:=F.one; z:=F.zero;;
gens:=[[[z,z,o,z],[z,z,z,o],[o,z,z,z],[z,o,z,z,]],
[[z,o,z,o],[o,o,o,o],[z,z,o,o],[z,z,o,z]]];;
G := Group (gens, gens[1]^0);;
