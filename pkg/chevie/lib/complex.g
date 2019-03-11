#############################################################################
##
#A  complex.g           CHEVIE library                            Jean Michel
##
#A  $Id: util.g,v 1.1 1997/01/21 13:46:35 gap Exp $
##
#Y  Copyright (C) may 1997 - 2001  University  Paris VII.
##
#  This file defines complex numbers on an arbitrary ring.
#
#   gap> Complex(0,1);
#   I
#   gap> last+1;
#   1+I
#   gap> last^2;
#   2I
#   gap> last^2;
#   -4
#   gap> Complex(E(3));
#   -1/2+(-1/2*E(12)^7+1/2*E(12)^11)I
#   gap> Complex(E(3)^2);
#   -1/2+(1/2*E(12)^7-1/2*E(12)^11)I
#   gap> last+last2;
#   -1
# 
# It can be combined with the file ``decimal.g''
# 
#   gap> Complex(E(3));
#   -1/2+(-1/2*E(12)^7+1/2*E(12)^11)I
#   gap> evalf(last);
#   0.5+0.8660254039I
# 
# In addition to the function Complex, this file defines two other
# functions: ComplexConjugate and IsComplex
#   
#   gap> x:=X(Rationals);;x.name:="x";;Complex(0,x);
#   xI
#   gap> last^2;
#   -x^2
#   gap> IsComplex(last);
#   true
#   gap> ComplexConjugate(last2);
#   -x^2
#   
##########################################################################
ComplexOps:=OperationsRecord("ComplexOps");

IsComplex:=n->IsRec(n) and IsBound(n.operations) and n.operations=ComplexOps;

Complex:=function(arg)local x;
  x:=arg[1];
  if Length(arg)=2 then return rec(r:=x,i:=arg[2],operations:=ComplexOps);
  elif IsComplex(x) then return x;
  elif IsCyc(x) then
       return rec(r:=(x+GaloisCyc(x,-1))/2,i:=(x-GaloisCyc(x,-1))/(2*E(4)),
		  operations:=ComplexOps);
  else return rec(r:=x,i:=0*x,operations:=ComplexOps);fi;
end;

ComplexOps.\+:=function(a,b)
  if IsList(a) then return List(a,x->x+b);
  elif IsList(b) then return List(b,x->a+x);fi;
  a:=Complex(a);b:=Complex(b);return Complex(a.r+b.r,a.i+b.i);
end;

ComplexOps.\-:=function(a,b)
  if IsList(a) then return List(a,x->x-b);
  elif IsList(b) then return List(b,x->a-x);fi;
  a:=Complex(a);b:=Complex(b);return Complex(a.r-b.r,a.i-b.i);
end;

ComplexOps.\*:=function(a,b)local res;
  if IsList(a) then return List(a,x->x*b);
  elif IsList(b) then return List(b,x->a*x);fi;
  a:=Complex(a);b:=Complex(b);return Complex(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);
end;

ComplexOps.ComplexConjugate:=a->Complex(a.r,-a.i);

#complex conjugation of lists, matrices, polynomials...
ComplexConjugate:=function(x)
  if IsCyc(x) then return GaloisCyc(x,-1);
  elif IsList(x) then return List(x,ComplexConjugate);
  elif IsPolynomial(x) then 
    x:=Copy(x);x.coefficients:=ComplexConjugate(x.coefficients);return x;
  elif IsUnknown(x) then return x;
  elif IsRec(x) and IsBound(x.operations) and 
    IsBound (x.operations.ComplexConjugate)then 
    return x.operations.ComplexConjugate(x);
  else Error("no method for ComplexConjugate(",x,")");
  fi;
end;

# try to covert x to a cyclotomic
Cyclotomic:=function(x)
  if IsCyc(x) then return x;
  elif IsComplex(x) then 
    if IsCyc(x.r) and IsCyc(x.i) then return x.r+E(4)*x.i;
    elif IsDecimal(x.r) and IsDecimal(x.i) then
      return Rational(x.r)+E(4)*Rational(x.i);
    else Error(x," is not a cyclotomic");
    fi;
  elif IsDecimal(x) then return Rational(x);
  else Error("no method for Cyclotomic(",x,")");
  fi;
end;

ComplexOps.Norm:=x->x.i*x.i+x.r*x.r;

ComplexOps.\/:=function(a,b)
  if IsList(a) then return List(a,x->x/b);fi;
  a:=Complex(a);b:=Complex(b);
  a:=a*ComplexConjugate(b);b:=ComplexOps.Norm(b);
  return Complex(a.r/b,a.i/b);
end;

ComplexOps.\^:=function(h,n)local p,N;
   if h=0*h then
     if n=0 then return Complex(h.r^0,h.i);
     else return h;fi;
   fi;
   p:=1;
   if n<0 then h:=1/h; n:=-n; fi;
   if n=1 then return h;fi;
   while n>0 do
    if n mod 2 <> 0 then p:=p*h;fi;
    h:=h*h;
    n:=QuoInt(n,2);
   od;
   return p;
end;

ComplexOps.Format:=function(a,option)local r,i;
  if IsBound(option.GAP) then return ConcatenationString("Complex(",
    Format(a.r,option),",",Format(a.i,option),")");
  else 
    r:=ShallowCopy(Format(a.r,option)); 
    if r="0" then r:="";fi;
    i:=Format(a.i,option);
    if i<>"0" then
      if '+' in i{[2..Length(i)]} or '-' in i{[2..Length(i)]} then
	if r<>"" then Add(r,'+');fi;
	Append(r,"("); Append(r,i);Add(r,')');
      else
	 if i[1]<>'-' and r<>"" then Add(r,'+');fi;
	 if i="-1" then Append(r,"-");
	 elif i<>"1" then Append(r,i);fi;
      fi;
      Add(r,'I');
    fi;
    if r="" then r:="0";fi;
    return r;
  fi;
end;

ComplexOps.String:=a->Format(a);

ComplexOps.Print:=function(a)Print(String(a));end;

ComplexOps.evalf:=function(arg)local x,p;
  x:=arg[1];
  if Length(arg)=2 then p:=arg[2];else p:=DecimalOps.prec;fi;
  x:=Complex(evalf(x.r,p),evalf(x.i,p));
  if IsComplex(x.r) then x.r:=x.r.r;fi;
  if IsComplex(x.i) then x.i:=x.i.r;fi;
  return x;
end;
