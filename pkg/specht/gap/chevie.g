#######################################################################
##  moduleElements.g - defines a "module" record structure           ##
##                                                                   ##
##     Andrew Mathas          mathas@maths.usyd.edu.au               ##
##     University of Sydney   Sydney, 2006                           ##
##                                                                   ##
#######################################################################

#D Combinatorial functions which use or interface with chevie

# these functions all require chevie
RequirePackage("chevie");

# return a Chevie permation corresponding to the tableau <t>
CoxeterPermTableau:=function(W,t)
  return EltWord(W,CoxeterWordTableau(t));
end;

# Given a Coxeter permutation in Sym(n) return the corresponding 
# permutation in Sym(n).
SymmPermCoxeterPerm:=function(W,w) local i;
  return Product(CoxeterWord(W,w), i->(i,i+1));
end;

# return the row semistandard mu-tableau of type nu corresponding to
# the permutation d.
RowStandardTableauDoubleCosetRep:=function(mu, nu, d)
  local t, tnu, row, col, r;
  t:=ActOnTableau(InitialTableau(mu),d); # t^mu*d
  tnu:=InitialTableau(nu);       # T*nu
  for row in [1..Length(t)] do
    for col in [1..Length(t[row])] do
      t[row][col]:=First([1..Length(tnu)], r->t[row][col] in tnu[r]);
    od;
  od;
  return t;
end;

# return the Young subgroup Sym(mu) as a Coxeter subgroup of <W>.
YoungSubgroup:=function(W,mu)
  return ReflectionSubgroup(W,GeneratorsYoungSubgroupComposition(mu));
end;


#U CompositionYoungSubgroupGenerators(W,J)
#U CompositionYoungSubgroupGenerators(n,J)
#M Return the composition mu such that W_J=Sym(mu), where J is a subset
#M of the simple roots [1..rk(W)] for W. In the second form n is the
#M degree of the the symmetric group containing W_J.
CompositionYoungSubgroupGenerators:=function(W, J) local rank, carry, mu, i;
  if IsGroup(W) then rank:=W.rank; else rank:=W-1; fi;
  if J=[] then return [1..rank+1]*0+1; fi;
  carry:=1;
  if J[1]>1 then mu:=[1..J[1]-1]*0+1; else mu:=[]; fi;
  for i in [2..Length(J)] do
    if J[i]=J[i-1]+1 then carry:=carry+1;
    else 
      Add(mu, carry+1); 
      Append(mu, [1..J[i]-J[i-1]-2]*0+1);
      carry:=1;
    fi;
  od;
  Add(mu, carry+1);
  if Sum(mu)<rank+1 then Append(mu, [1..rank+1-Sum(mu)]*0+1); fi;
  return mu;
end;

# Return the set of super (H,K)-coset representatives. That is, the set   
# of distinguished (H,K)-double representatives which also satisfy the    
# property that alpha_K d\cup beta_H=\0 and beta_K d\cup \alpha_H=\0.     
DistinguishedDoubleCosetRepresentatives:=function(W, H,K)
  local reps, sd, d, i;
  reps:=[];
  for d in ReducedRightCosetRepresentatives(W, H) do
    if d^-1=ReducedInCoxeterCoset(K, d^-1) then # (H,K)-double coset representative 
      Add(reps, d);
    fi;
  od;
  return reps;
end;

DistinguishedRightCosetRepresentatives:=ReducedRightCosetRepresentatives;
DistinguishedRightCosetPermutations:=function(W,K) local w;
  return List(DistinguishedRightCosetRepresentatives(W,K),
        w->SymmPermCoxeterPerm(W,w));
end;
DistinguishedLeftCosetRepresentatives:=function(W,H) local w;
  return List(ReducedRightCosetRepresentatives(W,H),w->w^-1);
end;

#U CoxeterPermSymmPerm(x)
#M Given a permutation <x> in Sym(n) return the corresponding CoxeterWord
#M (we assume that <x> does indeed belong to W=Sym(n) (certainly it belongs
#M in some symmetric group).
CoxeterPermSymmPerm:=function(W,x) local w, i;
  w:=[];
  while x<>() do
    i:=1; while i^x<(i+1)^x do i:=i+1; od;
    Add(w,i);
    x:=(i,i+1)*x;
  od;
  return EltWord(W,w);
end;

# Given a Coxeter word in Sym(n) return the corresponding permutation 
# in Sym(n).
SymmPermCoxeterWord:=function(w) local i;
  if w=[] then return ();
  else return Product(w, i->(i,i+1));
  fi;
end;

#U Multicombinations(set, size)............................................
#M Return the list of all unordered partitions of the set <set> into pieces
#M determined by the set <size>.                                           
Multicombinations:=function(set,size) local W, Y, d, i, comb;
  # first create the group which generates these partitions
  if Length(set)=0 then return List(size, i->[]); fi;
  comb:=[];
  d:=1;
  for i in size do
    Add(comb,[d..d+i-1]); 
    d:=d+i; 
  od;
  if Length(set)=1 then return [comb+set[1]-1]; fi;
  W:=CoxeterGroup("A",Length(set)-1);
  Y:=YoungSubgroup(W,size);
  return List(DistinguishedRightCosetPermutations(W,Y), 
                   w->List(comb,i->set{List(i,d->d^w)}));
end;

