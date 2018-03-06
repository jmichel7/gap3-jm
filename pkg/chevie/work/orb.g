orbit:=function(M,v)local res,n,new;
  res:=[]; new:=[v];
  repeat
    n:=Set(Concatenation(List(M,g->new*g)));
    UniteSet(res,new);
    Print(Length(res),"\n");
    new:=Difference(n,res);
  until Length(new)=0;
  return res;
end;

orbits:=function(W,l)local res,o;
  res:=[];
  while Length(l)>0 do
    o:=orbit(W.matgens,l[1]);
    Add(res,o);
    l:=Difference(l,o);
    Print(Length(l)," remaining\n");
  od;
  return res;
end;

# points fixes communs a tous les elements du groupe  H
FixSpace:=function(W,H)local l;
  l:=H.generators;
  if Length(l)=0 then return 
    VectorSpace(W.roots{W.generatingReflections},Cyclotomics);fi;
  l:=List(l,function(w)local m;
    m:=MatXPerm(W,w);
    return NullspaceMat(m-m^0);end);
  l:=List(l,x->VectorSpace(x,Cyclotomics));
  return Intersection(l);
end;

IsStandard:=H->ForAny(Arrangements(H.generators,Length(H.generators)),
   a->CheckRelations(ReflectionType(H),a));
