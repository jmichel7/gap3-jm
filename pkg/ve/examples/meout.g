field := GF(3);
VErunTime := 16;
gen_gena := field.one*[[1,0,0,0],
[0,0,1,0],
[0,0,1,0],
[0,0,0,1]];
gen_genb := field.one*[[0,1,0,0],
[0,1,0,0],
[0,0,0,1],
[0,0,0,1]];
gens := [gen_gena,gen_genb];
