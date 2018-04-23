##############################################################################
##
#A  mvrf.g       VKCURVE package         Gwenaelle Genet and Jean Michel
##  Multivariate Rational Fractions
##
#Y  Copyright (C) 2002  University Paris VII, France.
##
##  This file was initially written by Gwenaelle Genet in January 2002 
## 
##  Rational Fractions are represented as
##  rec(num:=numerator (an Mvp which is a true polynomial)
##      den:=denominator (an Mvp which is a true polynomial)
##      operations:=FracRatOps)
##
##  Computations are done in the field K(x1,x2,...,xn)
##  where K is the (presumed field) of coefficients of the Mvps.
##  For the moment K is presumed to be a subfield of the cyclotomics.
#############################################################################

RatFracOps:=OperationsRecord("RatFracOps");

IsRatFrac:=x->IsRec(x) and IsBound(x.operations) and RatFracOps=x.operations;

# Given a list arg of Laurent polynomials, returns the unique monomial 
# m such that m*arg is a list of true polynomials, one of them with a
# nonzero constant term
if VKCURVE.mvp=1 then
LaurentDenominator := function(arg)local res, v;
  res:=MonomialGcd(Set(Concatenation(List(arg,function(x)
    if IsRec(x) then return x.elm;else return [];fi;end))));
  res.coeff:=List(res.coeff,
    function(x)if x>=0 then return x-Int(x);else return x;fi;end);
  v:=Filtered([1..Length(res.coeff)],i->res.coeff[i]<>0);
  return rec(elm:=[rec(coeff:=-res.coeff{v},elm:=res.elm{v})],
             coeff:=[1],operations:=MvpOps);
end;
elif VKCURVE.mvp=2 then
LaurentDenominator := function(arg)local res, f;
  res:=MonomialGcd(Set(Concatenation(List(arg,function(x)
    if IsRec(x) then return x.elm;else return [];fi;end))));
  for f in RecFields(res) do
    if res.(f)>=0 then res.(f):=Int(res.(f))-res.(f);
      if res.(f)=0 then Unbind(res.(f));fi;
    else res.(f):=-res.(f);
    fi;
  od;
  return MvpOps.Mvp([res],[1]);
end;
else
LaurentDenominator := function(arg)local res, v;
  res:=MonomialGcd(Set(Concatenation(List(arg,function(x)
    if IsRec(x) then return List(x.pairs,y->y[1]);else return [];fi;end))));
  for p in res.pairs do if p[2]>=0 then p[2]:=Int(p[2])-p[2];
                                   else p[2]:=-p[2];fi;od;
  res.pairs:=Filtered(res.pairs,x->x[2]<>0);
  return Mvp([[res,1]]);
end;
fi;
    
# Returns a Lcm of the arguments that are PolMvps or scalars.
MvpLcm := function(arg) local lcm,i;
  if ForAny(arg,x->x=0*x) then return Mvp(0);fi;
  arg:=Filtered(arg,x->ScalMvp(x)<>1);
  if Length(arg)=0 then return Mvp(1); fi;
  lcm:=arg[1];
  for i in [2..Length(arg)] do 
#   Print("ExactDiv(",arg[i],",",MvpGcd(lcm,arg[i]),")\n");
    lcm:=lcm*MvpOps.ExactDiv(arg[i],MvpGcd(lcm,arg[i]));
  od;
  return lcm/lcm.coeff[Length(lcm.elm)];
end;

RatFracOps.RatFrac:=function(n,d)
  if Length(d.coeff)=0 then Error("division by zero");fi;
  if d.coeff[1]<>1 then
#   Print("Making ",n,"/",d,"\n");
#   if Valuation(n/d.coeff[1])<0 then Error();fi;
#   if Valuation(d/d.coeff[1])<0 then Error();fi;
    return rec(num:=n/d.coeff[1],den:=d/d.coeff[1],operations:=RatFracOps);
  else return rec(num:=n,den:=d,operations:=RatFracOps);
  fi;
end;
    
# takes one or two Mvps, returning arg[1] as a RatFrac (resp arg[1] / arg[2])
RatFrac:=function(arg) local p,n;
  if Length(arg) = 1 then 
    n:=arg[1];
    if IsRatFrac(n) then return n;
    elif IsMvp(n) then p:=LaurentDenominator(n);
#     Print("calling RatFrac(",p*n,",",p,")\n");
      return RatFracOps.RatFrac(p*n,p);
    else # assume n is a 'scalar'
      return rec(num:=Mvp(n),den:=Mvp(1),operations:=RatFracOps);
    fi;
  elif Length(arg)=2 then 
    if arg[2]=0*arg[2] then Error("denominator is ",arg[2]);fi;
    arg:=arg*LaurentDenominator(arg[1],arg[2]);
    p:=MvpGcd(arg[1],arg[2]);
#     Print("calling2 RatFrac(",arg[1],",",p,")\n");
#     Print("ExactDiv2(",arg[1],",",p,")\n");
#     Print("ExactDiv2b(",arg[2],",",p,")\n");
    if ScalMvp(p)<>false then return RatFracOps.RatFrac(arg[1],arg[2]);fi;
    return RatFracOps.RatFrac(MvpOps.ExactDiv(arg[1],p),
      MvpOps.ExactDiv(arg[2],p));
  else Error("The number of arguments must be 1 or 2. \n");
  fi;
end;

RatFracOps.\*:=function(x,y) local d1,d2;
  if IsRatFrac(x) then 
     if IsRatFrac(y) then 
       d1:=MvpGcd(x.num,y.den);d2:=MvpGcd(y.num,x.den);
#      Print("calling3 RatFrac(",x,",",y,")\n");
#      Print("ExactDiv3(",x.num,",",d1,")\n");
#      Print("ExactDiv3b(",y.num,",",d2,")\n");
#      Print("ExactDiv3c(",y.den,",",d2,")\n");
#      Print("ExactDiv3d(",y.den,",",d1,")\n");
       return RatFracOps.RatFrac(MvpOps.ExactDiv(x.num,d1)*
         MvpOps.ExactDiv(y.num,d2),
         MvpOps.ExactDiv(x.den,d2)*MvpOps.ExactDiv(y.den,d1));
     elif IsList(y) then return List(y, z -> x*z);
     elif IsMvp(y) then return RatFrac(x.num*y,x.den);
     else 
#       Print("calling4 RatFrac(",x.num*y,",",x.den,")\n");
        return RatFracOps.RatFrac(x.num*y,x.den);# assume y is a 'scalar'
     fi;
  elif IsList(x) then return List(x, z -> z*y);
  elif IsMvp(x) then return RatFrac(x*y.num,y.den);
  else 
#    Print("calling5 RatFrac(",y.num*x,",",y.den,")\n");
     return RatFracOps.RatFrac(y.num*x,y.den);# assume x is a 'scalar'
  fi;
end;           

RatFracOps.\/:=function(x,y)local d1,d2;
  if IsRatFrac(y) then 
    if IsRatFrac(x) then 
      d1:=MvpGcd(x.num,y.num);d2:=MvpGcd(x.den,y.den);
#     Print("calling6 RatFrac(",d1,",",d2,")\n");
      return RatFracOps.RatFrac(MvpOps.ExactDiv(x.num,d1)*
         MvpOps.ExactDiv(y.den,d2),
         MvpOps.ExactDiv(x.den,d2)*MvpOps.ExactDiv(y.num,d1));
    elif IsMvp(x) then return RatFrac(x)*RatFracOps.RatFrac(y.den,y.num);
    elif IsList(x) then return List(x, e->e/y);
    else return RatFracOps.RatFrac(y.den*x,y.num); # assume x is a 'scalar'
    fi;
  else return RatFracOps.RatFrac(x.num/y,x.den);# assume y is a 'scalar'
  fi;
end;  

RatFracOps.\+:=function(x,y) local d,yd;
  if IsRatFrac(y) then 
    if IsMvp(x) then return RatFrac(x)+y;
    elif IsRatFrac(x) then 
      d:=MvpGcd(x.den,y.den);yd:=MvpOps.ExactDiv(y.den,d);
      return RatFrac(x.num*yd+MvpOps.ExactDiv(x.den,d)*y.num,x.den*yd);
    elif IsList(x) then return List(x, e->e+y); 
    else return RatFracOps.RatFrac(x*y.den+y.num,y.den);# assume x is a 'scalar'
    fi;
  elif IsList(y) then return List(y,z->x+z);
  else return RatFracOps.RatFrac(x.den*y+x.num,x.den);# assume y is a 'scalar'
  fi;
end;       

RatFracOps.\-:=function(x,y) 
  return x + (-1)*y; end; 

############          Integral powers           ############
RatFracOps.\^:=function(x,i) local y;
  if i=0 then return RatFrac(x.num^0);fi;
  if not IsInt(i) then Error("only integer powers implemented"); fi;
  if i<0 then return RatFracOps.RatFrac(x.den,x.num)^-i;fi;
  y:=1;
  while i>0 do 
    if i mod 2 <> 0 then y:=y*x;fi;
    if i>=2 then x:=x*x;fi;
    i:=QuoInt(i,2); 
  od;
  return y; 
end;

###########     Functions String and Print     ###########
RatFracOps.Format:=function(f,o) local n,d;
  n:=Format(f.num,o); d:=Format(f.den,o);
  if d="1" then return n; fi;
  if Length(Intersection("+-/*",d{[2..Length(d)]}))<>0  then
    d:=Concatenation("(",d,")"); fi;
  if Length(Intersection("+-/",n{[2..Length(n)]}))<>0  then
    n:=Concatenation("(",n,")"); fi;
  return Concatenation(n,"/",d);
end;

RatFracOps.String:=f->Format(f);

RatFracOps.Print:=function(f) Print(String(f));end;

RatFracOps.Variables:=x->Union(Variables(x.num),Variables(x.den));

MvpGcd:=function(arg)
  local i,gcd,coef,cont,reseuc,listvar,l,VecGcd,Vecmod,v,nvar;
  arg:=Filtered(arg,x->ScalMvp(x)<>0*x);
  if Length(arg)=0 then return Mvp(0); 
  elif Length(arg)=1 then return arg[1];
  elif Length(arg)>2 then
    gcd:=arg[1];
    for i in [2..Length(arg)] do gcd:=MvpGcd(gcd,arg[i]); od;
    return gcd; 
  fi;
  if ForAny(arg,x->ScalMvp(x)<>false) then return arg[1]^0; fi;

  if ForAny(arg,x->Length(x.coeff)=1) then
  # first find highest monomial which divides both
    return rec(coeff:=[1],elm:=[MonomialGcd(Concatenation(arg[1].elm,arg[2].elm))],
    operations:=MvpOps);
  fi;
# Print("arg=",arg,"\n");
  listvar:=List(arg,Variables);
  v:=ApplyFunc(SymmetricDifference,listvar);
  listvar:=Union(listvar);
  nvar:=Length(listvar);
# if nvar=2 then Print("v=",v," listvar=",listvar,"\n");fi;
  if Length(v)>0 then v:=v[1];
  else v:=listvar[1];
  fi;
  coef:=List(arg,x->Coefficients(x,v));
# if nvar=2 then Print("coef=",coef,"\n");fi;

  Vecmod := function(p, q)local lp,lq,plp,qlq,f;
    lq:=Length(q);if lq=1 then return [];fi;
    qlq:=q[lq];lp:=Length(p);p:=ShallowCopy(p);
    f:=ScalMvp(qlq);q:=q{[1..lq-1]};
    if f<>false then q:=q/f;fi;
    while lp>=lq do
      plp:=p[lp]; lp:=lp-1; 
      if f=false then p:=p*qlq;fi;
      p{[lp-lq+2..lp]}:=p{[lp-lq+2..lp]}*plp^0-plp*q;
      while lp>0 and p[lp]=0*p[lp] do lp:=lp-1;od;
    od;
    p:=p{[1..lp]};
    if lp=0 then return p;fi;
    plp:=ScalMvp(p[lp]);
    if plp=false then return p/ApplyFunc(MvpGcd,p);fi;
    return p/plp;
  end;

  VecGcd:=function(p,q)local tmp;
    while Length(q)>0 do tmp:=Vecmod(p,q);p:=q;q:=tmp;od;
    return p/p[Length(p)];
  end;
  
  if nvar=1 then return Mvp(v,VecGcd(coef[1],coef[2]));fi; # faster
  cont:=List(coef,x->ApplyFunc(MvpGcd,x));
  if ScalMvp(cont[1]=false) then
    coef[1]:=List(coef[1],x->MvpOps.ExactDiv(x,cont[1]));
  else coef[1]:=List(coef[1],Mvp);cont[1]:=Mvp(cont[1]);
  fi;
  if ScalMvp(cont[2]=false) then
    coef[2]:=List(coef[2],x->MvpOps.ExactDiv(x,cont[2]));
  else coef[2]:=List(coef[2],Mvp);cont[2]:=Mvp(cont[2]);
  fi;
  reseuc:=VecGcd(coef[1],coef[2]);
  l:=ApplyFunc(MvpLcm,List(reseuc,x->RatFrac(x).den));
  reseuc:=List(reseuc,function(p)
      if IsRatFrac(p) then return p.num*MvpOps.ExactDiv(l,p.den);
      else return p*l;fi;end);
# if not ForAll(reseuc,IsMvp) then Error();fi;
# if not IsMvp(MvpGcd(cont[1],cont[2])) then Error();fi;
  return MvpGcd(cont[1],cont[2])* ValuePol(reseuc,Mvp(v));
end;

RatFracOps.Value:=function(f,v)local i;
  if Length(v)=0 then return f;fi;
  if not IsRatFrac(v[2]) then f:=Value(f.num,v{[1,2]})/
                                            Value(f.den,v{[1,2]});
  else f:=ValuePol(Coefficients(f.num,v[1]),v[2])/
               ValuePol(Coefficients(f.den,v[1]),v[2]);
  fi;
  return Value(f,v{[3..Length(v)]});
end;

RatFracOps.ComplexConjugate:=f->RatFrac(ComplexConjugate(f.num),ComplexConjugate(f.den));

RatFracOps.\=:=function(x,y)
  if IsCyc(y) then return x=RatFrac(y);
  elif IsCyc(x) then return y=RatFrac(x);
  elif IsRatFrac(x) then return IsRatFrac(y) and 
      x.num=y.num and  (x.num=0*x.num or x.den=y.den);
  elif IsMvp(x) then return y=RatFrac(x);
  else return y=x;
  fi;
end;

RatFracOps.CycPol:=p->CycPol(p.num)/CycPol(p.den);

RatFracOps.Galois:=function(x,e)
  x:=ShallowCopy(x);x.num:=Galois(x.num,e);x.den:=Galois(x.den,e);
  return x;
end;
