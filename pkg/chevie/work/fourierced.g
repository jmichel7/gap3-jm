
TailleFourierSousgroupe:=function(G,H)
  return Sum(ConjugacyClasses(G),function(i)local res,C,l,ct;
    i:=Representative(i);
    C:=Centralizer(G,i);
    l:=List(DoubleCosets(G,C,H),Representative);
    l:=Filtered(l,j->j*i*j^-1 in H);
    l:=List(l,j->PermutationCharacter(C,Intersection(C,ConjugateSubgroup(H,j))));
    ct:=CharTable(C);
    return Number(ct.irreducibles,k->ForAny(l,j->ScalarProduct(ct,j,k)>0));
  end);
end;

TraceFourierChar:=function(G)local conj;
  conj:=List(ConjugacyClasses(G),Representative);
  return List(conj,function(i)local C;
    C:=Centralizer(G,i);
    return Size(C)*Number(List(DoubleCosets(G,C,C),Representative),
       j->i^j in C and j in DoubleCoset(C,j^-1,C));
  end);
end;

TraceFourier:=G->ScalarProduct(CharTable(G),TraceFourierChar(G),
  [1..Length(ConjugacyClasses(G))]*0+1);
