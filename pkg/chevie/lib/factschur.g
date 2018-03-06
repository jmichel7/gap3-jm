#############################################################################
##
#A  factschur.g          CHEVIE library        Maria Chlouveraki, Jean Michel
##
#Y  Copyright (C) 2010 - 2016   University  Paris VII.
##
##  This  file  contains  routines  to compute and work with factorized
##  Schur elements according to Maria's Ph. D. thesis.
##
FactorizedSchurElementsOps:=OperationsRecord("FactorizedSchurElementsOps");

IsFactorizedSchurElement:=x->IsRec(x) and IsBound(x.operations) and
  x.operations=FactorizedSchurElementsOps;

FactorizedSchurElementsOps.Format:=function(x,options)local v,res;
  v:=List(x.vcyc,function(l)
    if IsBound(l.cyc) then 
      if IsBound(options.Maple) then
        return SPrint("(",Format(Value(l.cyc,l.monomial),options),")");
      else return SPrint("P",l.cyc,"(",Format(l.monomial,options),")");
      fi;
    else 
      if IsBound(options.Maple) then
        return SPrint("(",Format(Value(l.pol,l.monomial),options),")");
      else
      return SPrint(Format(l.pol,options),
       "(",Format(l.monomial,options),")");
      fi;
    fi;
  end);
  if IsBound(options.GAP) or IsBound(options.Maple) then v:=Join(v,"*");
  else v:=Join(v,""); fi;
  return FormatCoefficient(x.factor,v,options);
end;

FactorizedSchurElementsOps.String:=x->Format(x,rec());

FactorizedSchurElementsOps.Print:=function(x)Print(String(x));end;

FactorizedSchurElementsOps.Expand:=function(x)return 
  x.factor*Product(x.vcyc,v->Value(v.pol,v.monomial));end;

FactorizedSchurElementsOps.\*:=function(a,b)local t,p;
  if IsFactorizedSchurElement(a) then
    a:=Copy(a);
    if IsFactorizedSchurElement(b) then
      for t in b.vcyc do
        p:=PositionProperty(a.vcyc,x->x.monomial=t.monomial);
        if p=false then Add(a.vcyc,t);
        else a.vcyc[p].pol:=a.vcyc[p].pol*t.pol;
        fi;
      od;
      a.vcyc:=Filtered(a.vcyc,x->Degree(x.pol)>0);
      a.factor:=a.factor*b.factor;
    else a.factor:=a.factor*b;
    fi;
  else
    b:=ShallowCopy(b);b.factor:=a*b.factor;a:=b;
  fi;
  return a;
end;

FactorizedSchurElementsOps.\/:=function(a,b)local t,p;
  if IsFactorizedSchurElement(a) then
    a:=Copy(a);
    if IsFactorizedSchurElement(b) then
      for t in b.vcyc do
        p:=PositionProperty(a.vcyc,x->x.monomial=t.monomial);
        if p=false then Add(a.vcyc,rec(monomial:=t.monomial,pol:=1/t.pol));
        else a.vcyc[p].pol:=a.vcyc[p].pol/t.pol;
        fi;
      od;
      a.vcyc:=Filtered(a.vcyc,x->Degree(x.pol)>0);
      a.factor:=a.factor/b.factor;
    else a.factor:=a.factor/b;
    fi;
  else
    b:=ShallowCopy(b);b.factor:=a/b.factor;
    a:=b;a.vcyc:=List(a.vcyc,t->rec(monomial:=t.monomial,pol:=1/t.pol));
  fi;
  return a;
end;

FactorizedSchurElementsOps.Value:=function(x,y)
   x:=ShallowCopy(x);
   x.factor:=Value(x.factor,y);
   x.vcyc:=List(x.vcyc,p->rec(pol:=p.pol,monomial:=Value(p.monomial,y)));
   return FactorizedSchurElementsOps.Simplify(x);
end;

# FactorizedSchurElementsOps.Simplify:=x->x  # to stop simplifying
FactorizedSchurElementsOps.Simplify:=function(res)local simplify;
# Print(res,"\n=>");
  res.vcyc:=List(res.vcyc,function(r)local k,c,v,n;
    k:=ScalMvp(r.monomial);
    if k<>false then res.factor:=res.factor*Value(r.pol,k); return false;fi;
# we arrange that the first variable appears to a positive power
    k:=r.monomial.elm[1].coeff;
    v:=rec(pol:=r.pol,monomial:=r.monomial);
    if k[1]<0 then
      v.pol:=CycPolOps.DescentOfScalars(v.pol,-1);
      k:=-k;v.monomial:=1/v.monomial;
      res.factor:=res.factor*v.pol.coeff*r.monomial^-v.pol.valuation;
      v.pol.coeff:=1;v.pol.valuation:=0;
    fi;
#     Print("m=",m," pol=",r.pol);
# If the coefficient of v.monomial is n*zeta where n^2 is a rational and zeta
# a root of unity we replace pol by pol(x*zeta) and divide v.monomial by zeta
    c:=v.monomial.coeff[1];n:=c*ComplexConjugate(c);
    if IsRat(n) then n:=GetRoot(n)/c;
      if n<>1 then
        v.monomial:=v.monomial*n;v.pol:=CycPolOps.EnnolaTwist(v.pol,1/n);
        res.factor:=res.factor*v.pol.coeff;v.pol:=v.pol/v.pol.coeff;
      fi;
    fi;
# We find the unique positive rational pow such that in r.monomial^pow
# - Each variable appears to an integral power.
# - The gcd of these powers is 1.
#     Print(" => m=",v.coeff," pol=",v.pol,"\n");
    v.power:=AbsInt(Gcd(List(k,Numerator)))/Lcm(List(k,Denominator));
    v.monomial:=v.monomial^(1/v.power);
    return v;
  end);
  res.vcyc:=Filtered(res.vcyc,x->x<>false);
# Print(res,"\n");
# the following function uses CycPol in order to obtain the factorization
# of the product of cyclotomic polynomials evaluated on different powers of the
# same monomial. The argument is a list of records with fields
# (coeff, pol, monomial, power) representing pol(coeff*(monomial)^power)
  simplify:=function(fil)local P,f,D,v;
#   Print(fil[1].monomial,":",Join(List(fil,x->Format(x.pol))),
#      " pow:",FormatGAP(List(fil,x->x.power)));
    D:=Lcm(List(fil,x->Denominator(x.power)));
  # We use an indeterminate representing m^(1/D)
    P:=Product(fil,x->CycPolOps.DescentOfScalars(x.pol,D*x.power));
# If all the terms appearing in P (which has valuation 0) are integral
# powers of some x^f, where f|D, then we change the indeterminate to m^(f/D)
    P:=Value(P,Indeterminate(Cyclotomics));
    f:=Filtered([1..Length(P.coefficients)],i->P.coefficients[i]<>0)-1;
    f:=Gcd(Gcd(f),D);
    if f>1 then
      P.coefficients:=P.coefficients{[1,1+f..Length(P.coefficients)]};
      D:=D/f;
    fi;
    v:=rec(pol:=CycPol(P),monomial:=fil[1].monomial^(1/D));
#   Print("=>",v.pol,"(",v.monomial,")\n");
    return v;
  end;
  res.vcyc:=List(CollectBy(res.vcyc,x->x.monomial),simplify);
  return res;
end;

FactorizedSchurElementsOps.Lcm:=function(l)
  l:=Concatenation(List(l,x->x.vcyc));
  l:=CollectBy(l,x->x.monomial);
  l:=List(l,x->[x[1].monomial,ApplyFunc(LcmCycPol,List(x,y->y.pol))]);
  return rec(factor:=1,vcyc:=List(l,x->rec(monomial:=x[1],pol:=x[2])),
    operations:=FactorizedSchurElementsOps);
end;
