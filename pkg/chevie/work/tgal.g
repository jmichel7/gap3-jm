# test galois actions

# test unipotent degrees permuted-with-signs by Gal(K/Q) by a permutation
# which fixes eigenvalues
testgal:=function(W)local K,e,ud,ei,x,res,l;
  K:=Field(Flat(CartanMat(W)));
  e:=Filtered(Elements(GaloisGroup(K)),x->x<>x^0);
  ud:=UnipotentDegrees(W,Mvp("x"));
  ei:=CollectBy([1..Length(ud)],Eigenvalues(UnipotentCharacters(W)));
  return List(e,function(x)local res; res:=[];
    for l in ei do res{l}:=SignPermuted(l,SignedPermListList(GaloisAction(ud{l},x),ud{l}));od;
    return res;
  end);
end;

# Field(W), families such that Field(Fourier(f))<>Field(W)
ff:=function(W)local K,KF,uc;
  K:=Field(Flat(CartanMat(W)));
  uc:=UnipotentCharacters(W);
  KF:=List(uc.families,f->Field(Flat(Fourier(f))));
  KF:=TransposedMat([KF,[1..Length(KF)]]);
  KF:=Filtered(KF,x->not ForAll(x[1].generators,y->y in K));
  if KF<>[] then Print(K,KF,"\n"); fi;
end;

# check parameters of zeta^-1 series complex conjugate of those of zeta-series
# checkzeta(W[,d])
checkzeta:=function(arg)local W,d,p,i,s,s1,r;
  W:=arg[1];
  if Length(arg)=2 then r:=arg[2];else r:=RegularEigenvalues(W);fi;
  for d in r do
    p:=Filtered(PrimeResidues(d),i->2*i<d);
    for i in p do
      s:=ComplexConjugate(Hecke(PrincipalSeries(W,i/d)).parameter);
      s:=List(s,x->Permuted(x,SortingPerm(x)));
      s1:=Hecke(PrincipalSeries(W,1-i/d)).parameter;
      s1:=List(s,x->Permuted(x,SortingPerm(x)));
      if s<>s1 then Error();fi;
    od;
  od;
end;

# build S^2-correspondence by zeta-series to zeta^-1-series
checkS2:=function(arg)local W,d,r,uc,s2;
  uc:=UnipotentCharacters(W);
  W:=arg[1];
  if Length(arg)=2 then r:=arg[2];else r:=RegularEigenvalues(W);fi;
  s2:=[];
  for d in r do
    uc:=1;
  od;
end;
