q:=X(Rationals);q.name:="q";

W:=ComplexReflectionGroup(28);
H:=Hecke(W,q^2);T:=Basis(H,"T");C:=Basis(H,"C'");
w0:=T(LongestCoxeterElement(W));

e:=Elements(W);; List(e,x->T(C(x)));;time;

W:=ComplexReflectionGroup(2,1,3);
H:=Hecke(W,q^2); C:=Basis(H,"C'");T:=Basis(H,"T");C(T(1));
p:=List([1..1000],i->[Random(W),Random(W)]);

a:=function(p)return ValuePol(KazhdanLusztigPolynomial(W,p[1],p[2]),q^2);end;
b:=function(p)return q^CoxeterLength(W,p[2])*Coefficient(
  H.operations.getCp(H,p[2]),p[1]);end;
