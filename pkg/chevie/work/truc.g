# get relations  from sq(r) and matrix O of eigenvalues of family
getrels:=function(M,O)local N;
   N:=Concatenation(M*ComplexConjugate(O)*M-O*M*O);
   Append(N,Concatenation(M^2*O-O*M^2));
   return Set(N);
end;

detfam2:=function(W,j)local M,f,O,res,p,v;
  M:=sq(detfam(W,j));
  f:=UnipotentCharacters(W).families[j];
  O:=DiagonalMat(f.eigenvalues);
  res:=rec(M:=M,f:=f.fourierMat,O:=O,rel:=sshr(getrels(M,O)));
  while true do
    p:=PositionProperty(res.rel,x->linear(x)<>false);
    if p=false then return res;fi;
    v:=linear(res.rel[p]);
    res.M:=Value(res.M,v);
    res.rel:=sshr(Value(res.rel,v));
  od;
end;

# tuples of unip chars with same degree and eigenvalue
db:=function(W)local uc,deg,eig;
  uc:=UnipotentCharacters(W);
  deg:=CycPolUnipotentDegrees(W);
  eig:=Eigenvalues(uc,[1..Length(deg)]);
  uc:=List([1..Length(deg)],i->[deg[i],eig[i],
    PositionProperty(uc.families,f->i in f.charNumbers),uc.TeXCharNames[i]]);
  Sort(uc);
  eig:=Filtered([1..Length(uc)-1],i->uc[i]{[1..3]}=uc[i+1]{[1..3]});
  eig:=Union(eig,eig+1);
  return uc{eig};
end;

sum:=function(M,i,j,k,special)
 return Sum([1..Length(M)],l->M[i][l]*M[j][l]*M[k][l]/M[special][l]);
end;
