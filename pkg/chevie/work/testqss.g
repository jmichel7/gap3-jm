sp:=n->RootDatum("spin",2*n);
hf:=n->RootDatum("halfspin",2*n);
so:=n->RootDatum("so",2*n);
pso:=n->RootDatum("pso",2*n);
d:=function(x)Display(UnipotentClasses(x),rec(order:=false));end;
g:=function(x)Print(FundamentalGroup(x),AlgebraicCentre(x).descAZ,"\n");end;
D:=function(n)Print(CHEVIE.D.FundamentalGroup(n),
       CHEVIE.D.CenterSimplyConnected(n),"\n");end;
p:=function(n,i)local W;W:=CoxeterGroup("D",n);
  return RestrictedPerm(LongestCoxeterElement(W)/
    LongestCoxeterElement(W,Difference([1..n],[i])),
    Concatenation([1..n],[2*W.N]));end;
w:=function(W)local g;
  g:=ApplyFunc(Group,AdjointFundamentalGroup(W));
  return List(FundamentalGroup(W).generators,x->GetWord(g,x));
end;
