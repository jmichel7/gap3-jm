

#U ErdmannsPartitionMap(S,mu,nu)                                          
#M Erdmann has proved that if e-p>0 (i.e. the classical case, q=1) then   
#M              [Delta(mu):L(nu)]=[S(t(mu)):D(t(nu))],                    
#M where t(tau) = p*tau' + (p-1)*(n-1,n-2,...,1,0). This function returns 
#M t(mu). Note: Erdmann defines t(tau) to be our t(tau')...               
ErdmannPartitionMap:=function(S,mu)
  local n, i;
  if S.e=S.p and S.p>0 then
    n:=Sum(mu);
    mu:=Copy(mu);
    for i in [1..n] do
      if i<=Length(mu) then mu[i]:=S.p*mu[i]+n-i;
      elif i<n then mu[i]:=n-i;
      fi;
    od;
    return mu;
  else return false;
  fi;
end;

ErdmannDecompositionNumber:=function(S,mu,nu)
  local tmu, tnu;
  tmu:=ErdmannPartitionMap(S,mu);
  tnu:=ErdmannPartitionMap(S,nu);
  if IsInRouquierBlock(S,tmu) and IsInRouquierBlock(S,tnu) then 
    return DecompositionNumberRouquierBlock(S,tmu,tnu);
  else return ErdmannDecompositionNumber(S,tmu,tnu);
  fi;
end;

ErdmannDecompositionMatrix:=function(S,n)
  local rows, cols, dmat, decnums, d, mu, nu;
  rows:=Partitions(n);
  if S.IsSpecht then cols:=Filtered(rows, mu->IsERegular(S,mu));
  else cols:=rows;
  fi;
  dmat:=EmptyDecompositionMatrix(S,rows,cols);
  for mu in [1..Length(cols)] do
    decnums:=rec(coeffs:=[], parts:=[]);
    for nu in [1..Length(rows)] do
      d:=ErdmannDecompositionNumber(S, rows[nu], cols[mu]);
      if d<>0 then # add into b
        Add(decnums.coeffs, d);
        Add(decnums.parts, nu);
      fi;
    od;
    Add(dmat.d, decnums);
  od;
  return dmat;
end;
    

