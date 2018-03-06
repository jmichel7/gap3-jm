e:=3;
t:=Mvp("t");q:=Mvp("q");
Id:=Sum([0..e-1],i->H(i,0))*(q-1)/(q^e-1);
H:=function(i,j)
  if i>=e then return H(i-e,j)+t^j*q^(i-e)*(q^e-1)*Id;fi;
  if j>=e then return H(i,j-e)+t^(j-e)*q^i*(t^e-1)*Id;fi;
  return Mvp(SPrint("H",i,j));end;

relE:=j->Sum([0..e-1],i->E(e)^(-i*j)*H(i,j))-t^j*(q^e-1)/(E(e)^-j*q-1)*Id;

relE1:=function(j,k)return Sum([0..e-1],i->E(e)^(-i*k)*H(i,j)-
t^(j-k)*E(e)^(-i*k)*H(i,k)+E(e)^(-i*j)*H(i,k)-t^(k-j)*E(e)^(-i*j)*H(i,j));end;

vars:=Concatenation(List([0..e-1],i->List([0..e-1],j->Mvp(SPrint("H",i,j)))));

tovec:=function(p)
  if IsRatFrac(p) then return tovec(p.num)/p.den;fi;
  return Concatenation(List([0..e-1],i->List([0..e-1],
    function(j)local v;
      v:=Coefficients(p,SPrint("H",i,j));
      if Length(v)>=2 then return v[2]*H(0,0)^0;else return 0*H(0,0);fi;end)));
end;

Sh:=function(p)local s;
  p:=Value(p,["q",q*t]);
  if IsRatFrac(p) then return Sh(p.num)/p.den;fi;
  s:=Concatenation(List([0..e-1],i->List([0..e-1],j->H(i,i+j))));
  p:=tovec(p)*s;
  return p;
end;

Om:=function(p)
  p:=Value(p,["t",q*t]);
  if IsRatFrac(p) then return Om(p.num)/p.den;fi;
  p:=tovec(p)*Concatenation(List([0..e-1],i->List([0..e-1],j->H(i+j,j))));
  return p;
end;

rels:=function()local v,i,j;
  v:=[];for i in [1..e-1] do Add(v,relE(i));od;
  for i in [0..e-1] do for j in [0..i-1] do Add(v,relE1(i,j));od;od;
  return v;
end;
