FindComplexConjugation:=function(W,i)local q,ud,p,i,j,p1,f,ambig,nambig,j;
  q:=X(Cyclotomics);q.name:="q";
  f:=UnipotentCharacters(W).families[i];
  ud:=UnipotentDegrees(W,q){f.charNumbers};
  p:=List(ud,p->PositionsSgn(ud,ComplexConjugate(p)));
  p1:=ComplexConjugate(f.eigenvalues);
  if Permuted(p1,SortingPerm(p1))<>
     Permuted(f.eigenvalues,SortingPerm(f.eigenvalues)) then
    InfoChevie("#I WARNING! the family ",i," is not real\n");
    return false;
  fi;
  p1:=List(f.eigenvalues,e->Positions(f.eigenvalues,ComplexConjugate(e)));
  for j in [1..Length(p)] do p[j]:=Filtered(p[j],x->AbsInt(x) in p1[j]);od;
  if ForAny(p,x->Length(x)=0) then
    InfoChevie("#I WARNING! no solution for family ",i," \n");
    return false;
  fi;
  ambig:=Filtered([1..Length(p)],i->Length(p[i])>1);
  nambig:=Difference([1..Length(p)],ambig);
  for i in nambig do p[i]:=p[i][1];od;
  for i in ambig do
    for j in nambig do
      p[i]:=Filtered(p[i],n->SignInt(n)*SignInt(p[j])*f.fourierMat[AbsInt(n)]
       [AbsInt(p[j])]=ComplexConjugate(f.fourierMat[i][j]));
    od;
  od;
  for i in ambig do if Length(p[i])=1 then p[i]:=p[i][1];fi;od;
  return p;
end;

All:=W->List([1..Length(UnipotentCharacters(W).families)],i->FindComplexConjugation(W,i));
