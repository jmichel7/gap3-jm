# tests the  word w in  p.generators is  central in the  presentation p.
# Actually  returns  the presentation  with  the  relations that  it  is
# central added.  It should  be checked  by usage  of ShrinkPresentation
# that the new presentation is the same.

TestCentral:=function(p,w)local g;
  p:=Copy(p);
  for g in p.generators do AddRelator(p,Comm(w,g)); od;
# ShrinkPresentation(p);
  return p;
end;

test:=function(n,p)local W,m;
  W:=ComplexReflectionGroup(n);
  m:=ReflectionDegrees(W);
  m:=m[Length(m)]/Gcd(m);
  Print("power to try:",m,"\n");
  n:=Length(W.generators);
  return List(Arrangements([1..n],n),a->TestCentral(p,
    Product(p.generators{a})^m));
end;
