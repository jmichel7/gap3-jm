# Invariant form stupidly computed
Ibete:=function(W)local h;
  h:=x->x*ComplexConjugate(TransposedMat(x));
  return Sum(Elements(W),x->h(MatXPerm(W,x)))/Size(W);
end;

# Invariant form computed from Lehrer-Taylor: if (,) invariant form then
#  C[i][j]/C[i][i]=(r_j,r_i)/(r_i,r_i)
InvariantHermitianForm:=function(W)local C,i,j,M,b,k,N,next;
  i:=IndependentLines(W.roots{W.generatingReflections});
  C:=CartanMat(W); C:=C{i}{i}; N:=W.roots{i}^-1;
  M:=C*0;
  for b in DecomposedMat(C) do
    # first fill in the diagonal terms
    M[b[1]][b[1]]:=1;next:=[1];
    while next<>[] do
      i:=b[next[1]];next:=next{[2..Length(next)]};
      for k in Filtered([1..Length(b)],k->M[b[k]][b[k]]=0 and C[i][b[k]]<>0) do
	j:=b[k];Add(next,k);
        M[j][j]:=M[i][i]*C[i][j]*ComplexConjugate(C[j][j])/
	     ComplexConjugate(C[j][i])/C[i][i];
      od;
    od;
    # then fill in the rest
    for i in b do for j in b do
      if i<>j then M[j][i]:=C[i][j]*M[i][i]/C[i][i];fi;
    od;od;
  od;
  M:=N*M*TransposedMat(ComplexConjugate(N));
  return M/M[1][1];
end;

# Check form M is invariant
checkform:=function(W,M)local rr,mm;
  rr:=W.simpleRoots{W.generatingReflections};mm:=W.matgens;
  if not ForAll(rr,r->ForAll(rr,s->ForAll(mm,m->
    (r*m)*M*ComplexConjugate(s*m)=r*M*ComplexConjugate(s)))) then
    Error();fi;
end;

checkLT:=function(W)local M,C,I,h;
  h:=function(x,M,y)return x*M*ComplexConjugate(y);end;
  M:=Ibete(W);C:=CartanMat(W);I:=W.generatingReflections;
  if not ForAll(I,i->ForAll(I,j->C[i][j]*h(W.roots[i],M,W.roots[i])=
    C[i][i]*h(W.roots[j],M,W.roots[i]))) then Error();fi;
end;

checkCartan:=function(W)local I,C,mm;
  I:=W.generatingReflections;
  C:=CartanMat(W);mm:=W.matgens;
  if not ForAll(I,i->W.roots[i]*mm[i]=
     E(1/W.EigenvaluesGeneratingReflections[i])*W.roots[i]) then Error();fi;
  if not ForAll(I,i->ForAll(I,
    j->W.roots[j]*mm[i]=W.roots[j]-C[i][j]*W.roots[i])) then Error();fi;
end;

# check invariant form agrees with stupid
checkI:=function(W)local a,b;
 a:=Ibete(W);a:=a/a[1][1];b:=InvariantHermitianForm(W);
 if a<>b then Error();fi;
end;
