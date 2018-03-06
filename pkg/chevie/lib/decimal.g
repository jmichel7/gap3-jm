#############################################################################
##
#A  decimal.g           CHEVIE library                            Jean Michel
##
#A  $Id: util.g,v 1.1 1997/01/21 13:46:35 gap Exp $
##
#Y  Copyright (C) may 1997 - 2001  University  Paris VII.
##
#  This file defines evalf(), Exp() and decimal numbers. 
#
# The main function defined in this  file is a function evalf similar to
# the Maple function with the same name.
#
#   gap> evalf(1/3);
#   0.3333333333
#
# As one  can see, the  resulting 'decimal numbers' have  an appropriate
# Print method (which uses the  String method). The precision (number of
# digits after the decimal point) in  which the conversion of numbers to
# decimal  approximations is  effected can  be  changed in  one call  by
# giving a second argument to evalf, or globally by calling the function
# SetDecimalPrecision.
#
#   gap> evalf(1/3,20);
#   0.33333333333333333333
#   gap> evalf(1/3);
#   0.3333333333
#   gap> SetDecimalPrecision(20);
#   gap> evalf(1/3);
#   0.33333333333333333333
#
# The result  is a ``decimal  number'', which is a  fixed-precision real
# number for which the  operations .+,-,*,/,^,<,GetRoot are defined. The
# precision of the result of an operation will that of the least precise
# number used.
#
#   gap> evalf(1/3)+1;
#   1.3333333333
#   gap> last^3;
#   2.3703703704
#
# The functions Exp, Log and Pi() are also defined by this file. Their code
# provides an example of how one can use decimal numbers.
#
#   gap> Exp(evalf(1));
#   2.7182818285
#
# It is  also possible  to evalf() cyclotomic  numbers. For  this reason
# this file needs also the file complex.g.
#
#   gap> evalf(ER(2));
#   1.4142135623
#   gap> evalf(E(3)); 
#   -0.5+0.8660254038I
#   gap> last^3;
#   1
#   gap> Exp(Pi()*evalf(E(4)));
#   -1
#
# To  make  evaluation  of  cyclotomics  relatively  fast,  a  cache  of
# previously computed primitive roots of unity is maintained (this cache
# is kept for each different precision ever used).
#
#  evalf() works  also  for  strings (the  result  is  truncated if  too
# precise):
#
#  gap> evalf(".3450000000000000000001");
#  0.345 # precision is 10
#
# and for lists (it is applied recursively to each element).
#
# There is a last function defined in this file we should mention.
# IsDecimal(x) returns true iff x is a decimal number.
#########################################################################

DecimalOps:=OperationsRecord("DecimalOps");

SetDecimalPrecision:=function(n)
  DecimalOps.guard:=10^(n+1);
  DecimalOps.precision:=n;
end;

SetDecimalPrecision(10);

IsDecimal:=n->IsRec(n) and IsBound(n.operations) and n.operations=DecimalOps;

Decimal:=function(x,guard)
  return rec(mantissa:=x*guard,guard:=guard,operations:=DecimalOps);
end;

Rational:=function(a)
  if IsRat( a )  then return a;
  elif IsDecimal(a) then return a.mantissa / a.guard;
  else Error("no method for Rational(",a,")");
  fi;
end;

# evalf(x [,precision])
evalf:=function(arg)local x,p,q,res,N,v,i,sign;
  x:=arg[1];
  if Length(arg)=1 then 
    if IsDecimal(x) then q:=x.guard;
    else p:=DecimalOps.precision;q:=DecimalOps.guard;
    fi;
  else p:=Maximum(arg[2],0);q:=10^(p+1);
  fi;
  if IsDecimal(x) then 
    x:=ShallowCopy(x);
    if q<x.guard then x.mantissa:=QuoInt(x.mantissa+x.guard/2/q,x.guard/q);
    elif q>x.guard then x.mantissa:=x.mantissa*q/x.guard;
    fi;
    x.guard:=q;
    return x;
  elif IsInt(x) then return Decimal(x,q);
  elif IsRat(x) then return  Numerator(x)/Decimal(Denominator(x),q);
  elif IsCyc(x) then# if you don't want to use complex numbers comment this out
    N:=NofCyc(x);
    res:=ValuePol(evalf(CoeffsCyc(x,N),p),DecimalOps.getRootOfOne(N,p));
    if GaloisCyc(x,-1)=x then res:=res.r;fi;
    return res;
  elif IsString(x) then
    if x[Length(x)]='I' then
      i:=PositionProperty(x{[2..Length(x)]},y->y in "+-");
      res:= Complex(evalf(x{[1..i]},p),evalf(x{[i+1..Length(x)-1]},p));
    fi;
    if x[1]='-' then sign:=-1;x:=x{[2..Length(x)]};
    elif x[1]='+' then sign:=1;x:=x{[2..Length(x)]};
    else sign:=1;fi;
    i:=Position(x,'.');
    if i=false then v:=1;
    else x:=Filtered(x,y->y<>'.');v:=10^(1+Length(x)-i);
    fi;
    x:=sign*ValuePol(List(Reversed(x),y->Position("0123456789",y))-1,10);
    if v<=q then return Decimal(x/v,q);
    else return rec(mantissa:=QuoInt(x,v/q),guard:=q,operations:=DecimalOps);
    fi;
  elif IsList(x) then return List(x,y->evalf(y,p));
  elif IsRec(x) and IsBound(x.operations) and IsBound(x.operations.evalf) then
    return  x.operations.evalf(x,p);
  else return x;
  fi;
end;

# check a and b are decimals and convert them to same (minimum) precision
DecimalOps.adjust:=function(a,b)local p,res;
  if IsDecimal(a) then p:=a.guard;
    if IsDecimal(b) and b.guard<p then p:=b.guard; fi;
  elif IsDecimal(b) then p:=b.guard;
  else p:=DecimalOps.guard;
  fi;
  if p=DecimalOps.guard then p:=DecimalOps.precision;
  else res:=0; while p>1 do res:=res+1;p:=p/10;od;p:=res-1;
  fi;
  res:=[evalf(a,p),evalf(b,p)];
  if IsComplex(res[1]) then res[2]:=Complex(res[2]);fi;
  if IsComplex(res[2]) then res[1]:=Complex(res[1]);fi;
  return res;
end;

DecimalOps.\<:=function(a,b)local res;
  res:=DecimalOps.adjust(a,b);
  return res[1].mantissa<res[2].mantissa;
end;

DecimalOps.\=:=function(a,b)local res;
  res:=DecimalOps.adjust(a,b);
  return ForAll(res,IsRec) and res[1].mantissa=res[2].mantissa;
end;

DecimalOps.\+:=function(a,b)local res;
  if IsList(a) then return List(a,x->x+b);
  elif IsList(b) then return List(b,x->a+x);fi;
  res:=DecimalOps.adjust(a,b);
  if IsComplex(res[1]) then return res[1]+res[2];fi;
  res[1].mantissa:=res[1].mantissa+res[2].mantissa;return res[1];
end;

DecimalOps.\-:=function(a,b)local res;
  if IsList(a) then return List(a,x->x-b);
  elif IsList(b) then return List(b,x->a-x);fi;
  res:=DecimalOps.adjust(a,b);
  res[1].mantissa:=res[1].mantissa-res[2].mantissa;return res[1];
end;

DecimalOps.roundquotient:=function(a,b)
  if a*b<0 then return QuoInt(2*a-b,2*b);else return QuoInt(2*a+b,2*b);fi;
end;

DecimalOps.\*:=function(a,b)local res;
  if IsList(a) then return List(a,x->x*b);
  elif IsList(b) then return List(b,x->a*x);fi;
  res:=DecimalOps.adjust(a,b);
  if IsComplex(res[1]) then return res[1]*res[2];fi;
  res[1].mantissa:=DecimalOps.roundquotient(res[1].mantissa*res[2].mantissa,
                   res[1].guard);
  return res[1];
end;

DecimalOps.\/:=function(a,b)local res;
  if IsList(a) then return List(a,x->x*b);fi;
  res:=DecimalOps.adjust(a,b);
  res[1].mantissa:=DecimalOps.roundquotient(res[1].mantissa*res[1].guard,
                                         res[2].mantissa);
  return res[1];
end;

DecimalOps.\^:=function(h,n)local p;
  p:=1;
  if not IsInt(n) then return GetRoot(h,Denominator(n))^Numerator(n);fi;
  if n<0 then h:=1/h;n:=-n;fi;
  while n>0 do
    if n mod 2 <> 0 then p:=p*h;fi;
    h:=h*h;
    n:=QuoInt(n,2);
  od;
  return p;
end;

DecimalOps.String:=function(a)local res,as,i,d;
  res:="";as:=a.mantissa;
  if as<0 then res:="-";as:=-as;fi;
  as:=DecimalOps.roundquotient(as,10);
  if a.guard<=1 then Error("precision should be at least 1");return;fi;
  d:=a.guard/10;
  Append(res,String(QuoInt(as,d)));as:=String(as mod d);
  as:=Concatenation(String(d/10^Length(as)),as);
  as:=as{[1..First([Length(as),Length(as)-1..1],x->as[x]<>'0')]};
  if Length(as)>1 then PrintToString(res,".",as{[2..Length(as)]});fi;
  if res="-0" then res:="0";fi;
  return String(res);
end;

DecimalOps.Print:=function(a)Print(String(a));end;

Exp:=function(x)local res,i,p,z;
  if IsCyc(x) then x:=evalf(x);fi;
  z:=0*x;res:=z;p:=1;i:=1;
  while p<>z do res:=p+res;p:=(1/i)*p*x;i:=i+1;od;
  return res;
end;

Pi:=function(arg)local a,i,r,pi,p;
  if Length(arg)=1 then p:=arg[1];else p:=DecimalOps.precision;fi;
  if not IsBound(DecimalOps.pi) or DecimalOps.pi.guard<10^(p+1) then
    pi:=0; a:=[0,2,2,1];i:=1;
    repeat
      r:=evalf(a[1+i mod 4]/((-4)^QuoInt(i,4)*i),p+1);
      pi:=pi+r;
      i:=i+1;
    until i mod 4<>1 and r=0;
    DecimalOps.pi:=evalf(pi,p);
  fi;
  return evalf(DecimalOps.pi,p);
end;

Log:=function(x)local e,res,i;
  if not IsDecimal(x) then x:=evalf(x);fi;
  if x<0 then Error(x," should be positive");fi;
  res:=Decimal(0,x.guard);
  e:=Exp(Decimal(1,x.guard));
  while x>1 do x:=x/e;res:=res+1;od;
  x:=1-x;e:=evalf(x);i:=1; # ln(1-x)=-\sum([1..],x^i/i)
  while e.mantissa>1 do res:=res-e/i;e:=e*x;i:=i+1;od;
  return res;
end;

DecimalOps.getRootOfOne:=function(N,p)local r;
  if not IsBound(DecimalOps.rootsOfOne) then DecimalOps.rootsOfOne:=rec();fi;
  r:=DecimalOps.rootsOfOne;
  if not IsBound(r.(N)) or r.(N).r.guard<10^(p+1) then
    r.(N):=Exp(2*Pi(p)*evalf(Complex(0,1),p)/N);
  fi;
  return evalf(r.(N),p);
end;

DecimalOps.GetRoot:=function(x,n)local z,new;
  z:=x;
# if not IsDecimal(z) then z:=evalf(z);fi;
  if n mod 2=0 and x<0 then Error("illegal: odd root of negative decimal");fi;
  if z=0*z then return z;fi;
  while true do
    new:=((n-1)*z+x/z^(n-1))/n; 
    if AbsInt(z-new)<10/z.guard then return new;fi;
    z:=new;
  od;
end;

DecimalOps.Format:=function(x,option)local p,res;
  if IsBound(option.GAP) then 
    p:=x.guard;res:=0; while p>1 do res:=res+1;p:=p/10;od;p:=res-1;
    return SPrint("evalf(",x.mantissa/x.guard,",",p,")");
  else return String(x);
  fi;
end;
