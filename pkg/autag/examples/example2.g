
# construct the group

f := FreeGroup(List([1..9],i->Concatenation("g",String(i))));

P := f / [ f.1^2/(f.2*f.4),  f.2^f.1/(f.2*f.5),  f.2^2/f.3,  
           f.3^f.1/(f.3*f.6),  f.3^f.2/f.3 ,  f.3^2/f.7 ,  
           f.4^f.1/(f.4*f.5^2*f.8),  f.4^f.2/(f.4*f.5*f.6^2*f.7*f.8^2),
           f.4^f.3/(f.4*f.5^2*f.6^2*f.8^2),  f.4^3/f.8^2,
           f.5^f.1/(f.5*f.6*f.7^2*f.9) ,  f.5^f.2/(f.5^2*f.6*f.8^2*f.9),
           f.5^f.3/(f.5^2*f.7^2*f.9),  f.5^f.4/(f.5*f.8),  f.5^3/f.8,  
           f.6^f.1/(f.5*f.8),  f.6^f.2/(f.5*f.6*f.7^2*f.9^2),  
           f.6^f.3/(f.6^2*f.7*f.8*f.9),  f.6^f.4/(f.6*f.8^2),  
           f.6^f.5/(f.6*f.9^2),  f.6^3/f.8^2,  f.7^f.1/(f.7^2*f.9), 
           f.7^f.2/f.7,  f.7^f.3/f.7,  f.7^f.4/(f.7*f.8*f.9),  f.7^f.5/f.7,  
           f.7^f.6/f.7,  f.7^3/f.8,  f.8^f.1/f.8^2,  f.8^f.2/f.8,  
           f.8^f.3/f.8,  f.8^f.4/f.8,  f.8^f.5/f.8,  
           f.8^f.6/f.8,  f.8^f.7/f.8,  f.8^3,  
           f.9^f.1/f.9^2,  f.9^f.2/f.9,  f.9^f.3/f.9,  
           f.9^f.4/f.9,  f.9^f.5/f.9,  f.9^f.6/f.9,  
           f.9^f.7/f.9,  f.9^f.8/f.9,  f.9^3 ];
G := SpecialAgGroup(AgGroupFpGroup(P));
G.name := "G";

# compute the automorphism group

A := AutGroupSagGroup(G);

# examine the structure of the automorphism group

AutGroupStructure(A);

# examine the 2nd factor group in more detail

factors := AutGroupFactors(A);;
IsAbelian(factors[2]);
IsCyclic(factors[2]); 
# so the matrix group is cyclic of order 8


# check whether conjugation by G.1*G.3*G.5 has weight 12 or more

series := AutGroupSeries(A);;
a := InnerAut(G, G.1*G.3*G.5);
series[4].weight;   # so series[4] is the weight 12 subgroup of <A>
a in series[4];     # returns false, so <a> has weight less than 12

