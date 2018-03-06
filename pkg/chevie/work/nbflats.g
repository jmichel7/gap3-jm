# Number of flats in the Hypeprplane arrangement of W
numberflats:=function(W) local l;
  l:=ParabolicRepresentatives(W);
  return List(l, x->Size(W)/Size(Stabilizer(W,x,OnSets))/
    Size(ReflectionSubgroup(W,x)));
end;  
