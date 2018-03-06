spieces:=function(W)local spec,H,b,tbl,q,fd,conjdet,O,irr,uc,l,a,p;
  q:=X(Rationals);q.name:="q";
  irr:=CharTable(W).irreducibles;
  p:=PermListList(irr,ComplexConjugate(irr));
  uc:=UnipotentCharacters(W);
  if uc<>false then
    spec:=List(uc.families,f->Position(uc.harishChandra[1].charNumbers,
			     f.charNumbers[f.special]));
    a:=ChevieCharInfo(W).a;
  else
    H:=Hecke(W,Mvp("x"));
    if IsBound(ChevieCharInfo(W).b) then b:=ChevieCharInfo(W).b;
    else b:=List(FakeDegrees(W,q),Valuation);
    fi;
    a:=-List(SchurElements(H),Valuation);
    spec:=List(RouquierBlocks(H),x->First(x,y->b[y]=Minimum(b{x})));
  fi;
  l:=[1..NrConjugacyClasses(W)];
  b:=LowestPowerFakeDegrees(W);
  conjdet:=PositionDet(W)^p;
  fd:=FakeDegrees(W,q);
  O:=q^Sum(ReflectionCoDegrees(W)+1)*List(l,
    i->List(l,j->q^(-b[i]-b[j])*fd*DecomposeTensor(W,i,j,conjdet)));
  PrintArray(O);
  tbl:=BigCellDecomposition(O,CollectBy(l,x->-a[x]));
  PrintArray(tbl[1]);
  l:=DiagonalOfMat(tbl[2]){spec};
  return [Sum(l),l];
end;

checkShoji:=function(W)local uc,a,t,iszero,i,j,n;
  iszero:=x->x=x*0;
  uc:=UnipotentClasses(W);
  a:=ChevieCharInfo(W).a;
  t:=ICCTable(uc);
  n:=CharNames(W);
  if t.locsys<>uc.springerSeries[1].locsys then Error("pas prevu!");fi;
  for i in [1..Length(a)] do
    for j in [1..Length(a)] do
      if not iszero(t.scalar[i][j]) and a[j]>a[i] then
        Print("for ",n[i],",a=",a[i]," coeff<>0 on ",n[j],",a=",a[j],"\n");
      fi;
    od;
  od;
end;

checkShojiA:=function(W)local uc,a,t,iszero,i,j,n;
  iszero:=x->x=x*0;
  uc:=UnipotentClasses(W);
  a:=ChevieCharInfo(W).A;
  t:=ICCTable(uc);
  n:=CharNames(W);
  if t.locsys<>uc.springerSeries[1].locsys then Error("pas prevu!");fi;
  for i in [1..Length(a)] do
    for j in [1..Length(a)] do
      if not iszero(t.scalar[i][j]) and a[j]>a[i] then
        Print("for ",n[i],",A=",a[i]," coeff<>0 on ",n[j],",A=",a[j],"\n");
      fi;
    od;
  od;
end;
