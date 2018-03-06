# compute .cycpolfakedegrees from List(Fakedegrees(W,x),CycPol)
Computecycpolfakedegrees:=function(p)local coeff,l;
  coeff:=p.coeff;
  l:=Concatenation(List(CycPolOps.decompose(p),x->[1..x[1]]*0+x[2]));
  if IsPolynomial(coeff) then
    coeff:=ScalMvp(Value(Value(coeff,x^(1/2)),["x",X(Rationals)]));
    if coeff.valuation<>0 then Error();fi;
    coeff:=coeff.coefficients;
  fi;
  return Concatenation([coeff,p.valuation],l);
end;
