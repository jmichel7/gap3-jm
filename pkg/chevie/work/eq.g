W:=ComplexReflectionGroup(28);
uc:=UnipotentCharacters(W);
f:=uc.families[11];
ud:=UnipotentDegrees(W,x){f.charNumbers};
fd:=FakeDegrees(uc,x){f.charNumbers};
S0:=Fourier(f);
O:=DiagonalMat(Eigenvalues(f));
ReadChv("work/mvptools");
S:=gmat(Length(O),"a");
l:=S*ud-fd;Append(l,TransposedMat(S)*fd-ComplexConjugate(ud));
Append(l,S*O^-1*ud-O*ud);
Append(l,S*O*ud-O^-1*S*fd);
# car S^2O^-1=O^-1S^2 apply to ud and use previous
# Append(l,O*ComplexConjugate(ud)-S*O^-1*ComplexConjugate(ud));
# comes from sh^2 preserves ud. Seems no new info
l:=Concatenation(List(l,x->Coefficients(x,"x")));
Append(l,Flat(S-TransposedMat(S)));
l:=sortbylg(l);
Print("f=",f,"\n");
# S*O^-1*S=O*S*O  (*)
# car SOSOSO=1 => SOS=O^-1 S^3 O^-1 => SOS=O^-1SO^-1S^2=> SO=O^-1 SO^-1S=>qed

t:=function(j,l0) # j an index in l0
  return Concatenation(List(l0,k->List([1..Length(S)],
    i->Sum(Zip(S[i],S[j],S[k],S[f.special],
   function(a,b,c,d)return a*b*ComplexConjugate(c)/d;end)))));end;

# g14.2 remains 105 -- 0 left by Ennola
# g14.5 two sols for last remaining variable
# G26.7 21 vars -- 0 left by Ennola
# g27.2 finish 3 vars by S^2O-OS^2 and (SO)^3
# g27.3 2 vars use above and also S^4=1
# g27.6 remains 15 hard variables =>2 by OSO and 4 sols (0 by Ennola)
# g27.9 like g27.3
# g27.10 remains 3 variables (solve by OSO)
# g28.11 remains 10 variables
# g29.8 two sols for last remaining variable (3 apres Ennola =>0)
# g30.13 
# g32.8 two sols for last remaining variable
# g32.12 Ennola=>1
# g34.9 two sols for last remaining variable
# g34.20 remains 55 (10 after Ennola->3 using quadratic=>2 sols)
# g34.32 two sols for last remaining variable
