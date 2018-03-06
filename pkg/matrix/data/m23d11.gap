F:=GaloisField(2);;
o:=F.one;; z:=F.zero;;
a:=[
[z, o, z, z, z, z, z, z, z, z, z],
[o, z, z, o, z, z, z, z, z, z, z],
[z, z, o, z, o, z, z, z, z, z, z],
[o, o, z, z, z, z, o, z, z, z, z],
[z, o, z, z, z, z, z, z, o, z, z],
[o, z, z, o, o, o, z, z, z, o, z],
[o, o, z, z, z, z, z, z, z, z, z],
[o, o, z, z, o, z, o, o, o, o, o],
[z, o, z, o, z, z, z, z, z, z, o],
[z, z, z, z, z, z, o, z, z, o, z],
[z, o, z, o, o, z, z, z, o, z, o]];;

b:=[
[z, o, o, z, z, z, z, z, z, z, z],
[z, z, o, z, z, z, z, z, z, z, z],
[o, z, o, o, o, o, z, z, z, z, z],
[z, z, z, z, o, z, z, o, z, z, z],
[z, z, z, o, o, o, o, z, o, o, z],
[o, z, z, o, z, o, z, z, o, z, o],
[z, o, o, o, o, z, z, o, o, z, z],
[z, z, z, z, z, o, z, z, z, o, z],
[z, z, o, o, z, z, o, o, o, o, z],
[z, o, z, o, z, o, o, z, z, o, z],
[z, o, z, o, z, o, z, z, z, o, o]];;
gens:=[a,b];;

G := Group (gens,gens[1]^0);;

