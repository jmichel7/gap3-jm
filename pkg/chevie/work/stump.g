test:=function(W)local exps,h,ir,r,sr,i,nr;
  exps:=p->Set(List(p.elm,x->x.coeff[1]));
  h:=Maximum(ReflectionDegrees(W));
  ir:=CharTable(W).irreducibles;
  r:=ChevieCharInfo(W).extRefl[2];
  sr:=exps(FakeDegrees(W,Mvp("x"))[r]);
  for i in PrimeResidues(h) do
    nr:=Position(ir,List(ir[r],x->GaloisCyc(x,i)));
    if Set(List(sr*i,x->x mod h))=exps(FakeDegrees(W,Mvp("x"))[nr]) then
      Print(i," OK\n");
    elif Set(List(sr*i,x->x mod h))=sr then Print(i," ID\n");
    else Print(i," NO\n");Error();
    fi;
  od;
end;
