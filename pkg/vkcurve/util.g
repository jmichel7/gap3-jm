##############################################################################
##
#A  util.g       VKCURVE package         Jean Michel
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds various utility functions.
## 
#############################################################################
BigNorm:=x->AbsInt(x.r)+AbsInt(x.i); # alternative to ComplexOps.Norm(x)

SmallNorm:=x->Maximum(AbsInt(x.r),AbsInt(x.i)); # another alternative

# 'vraie' norme d'un nombre complexe
Rho:=function(a)
  a:=ComplexOps.Norm(a);
  if IsRat(a) then a:=evalf(a,-DecimalLog(a)+5);fi;
  return GetRoot(a,2);
end;

# For a non-zero rational x, returns k such that 10^k<|x|<=10^(k+1) 
#
DecimalLog:=function(x)local d;
  if x=0*x then Error("trying to take decimal log of 0");fi;
  d:=QuoInt(Denominator(x),Numerator(x));
  if d=0 then d:=QuoInt(Numerator(x),Denominator(x));
     return Length(String(AbsInt(d)))-1;
  else return -Length(String(AbsInt(d)));fi;
end;

ComplexRational:=function(a)
  if IsCyc(a) then a:=Complex(a);
    if not IsRat(a.r) then a.r:=evalf(a.r);fi;
    if not IsRat(a.i) then a.i:=evalf(a.i);fi;
  fi;
  return Complex(Rational(a.r),Rational(a.i));
end;

ComplexRatToGaussian:=a->a.r+E(4)*a.i;

# returns [minimum of l, position of that minimum in l]
#
MinPos:=function(l)local min,pos,i; pos:=1;min:=l[1];
  for i in [2..Length(l)] do if l[i]<min then min:=l[i];pos:=i;fi;od;
  return [min,pos];
end;

# minimum distance between 2 points of v
# if the minimum m is realized between p1 and p2,
# returns [m,[index of p1 in v, index of p2 in v]]
#
Dispersal:=function(v)local l,p,minpos;
  l:=Concatenation(List([1..Length(v)],i->List([1..i-1],j->[j,i])));
  p:=MinPos(List(l,x->BigNorm(v[x[1]]-v[x[2]])));
  return [p[1],l[p[2]]];
end;

# distance of z to segment [a,b]
DistSeg:=function(z,a,b) local r;
  b:=b-a;z:=z-a;r:=Rho(b);z:=r/b*z;
  if z.r<0 then return Rho(z);
  elif z.r>r then return Rho(z-b);
  elif z.i>0 then return z.i;
  else return -z.i;
  fi;
end;

PrimeCoeffs:=function(p)local F,q;
  F:=DefaultField(p.coeff);
  q:=List(Coefficients(p,"y"),q->Polynomial(F,ScalMvp(Coefficients(q,"x"))));
  q:=Gcd(q);
  if Degree(q)>0 then
    Error("coefficients of ",p," have gcd=",Value(q,Mvp("x")),"\n");
  fi;
end;

# discriminant with respect to x of an Mvp in x,y with rational coeffs
Discy:=function(p)local n,v;
  n:=2*Length(Coefficients(p,"x"))*Length(Coefficients(p,"y"));
  v:=List([1..n],function(i)local q;
    q:=ScalMvp(Coefficients(Value(p,["y",i]),"x"));
    if Length(q)=0 then return 0;
    elif Length(q)=1 then return q[1];
    else return DeterminantMat(ResultantMat(q,Derivative(q)));
    fi;
    end);
  v:=InterpolatedPolynomial(DefaultField(v),[1..n],v);
  return Value(v,Mvp("y"));
end;

# Resultant matrix of polynomials with vector-of-coeffs v, w
ResultantMat:=function(v,w)local m,i;
  v:=Reversed(v);
  w:=Reversed(w);
  m:=List([1..Length(v)+Length(w)-2],x->[1..Length(v)+Length(w)-2]*0*v[1]);
  for i in [1..Length(w)-1] do m[i]{[i..i+Length(v)-1]}:=v;od;
  for i in [1..Length(v)-1] do m[i+Length(w)-1]{[i..i+Length(w)-1]}:=w;od;
  return m;
end;

# Cut(string[,opt])
# opt: option record with fields
#    .width  cutting width  [default SizeScreen()[1]-2 ] 
#    .after  cutting after these chars [default ","]
#    .before cutting before these chars [default ""]
#    .file   where to print result [default stdout]
Cut:=function(arg)local opt,s,a,pa,pb,wr,pos,res;
  if Length(arg)>1 then opt:=ShallowCopy(arg[2]);else opt:=rec();fi;
  if not IsBound(opt.width) then opt.width:=SizeScreen()[1]-2;fi;
  if IsBound(opt.places) then
    Print("option .places is deprecated: use .after");
    opt.after:=opt.places;
  fi;
  if not IsBound(opt.after) then opt.after:=",";fi;
  if not IsBound(opt.before) then opt.before:="";fi;
  res:="";
  a:=Split(arg[1],'\n');
  if a[Length(a)]="" then a:=a{[1..Length(a)-1]};fi;
  for s in a do
    if s=[] then s:="";fi; # fix a bug in GAP3
    pos:=0;
    while Length(s)>pos+opt.width do
      pa:=PositionProperty(Reversed(pos+[1..opt.width]),x->s[x] in opt.after);
      pb:=PositionProperty(Reversed(pos+[1..opt.width+1]),
         x->s[x] in opt.before);
      if pa<>false and (pb=false or pa<=pb) then
        PrintToString(res,s{pos+[1..opt.width+1-pa]},"\n");
        pos:=pos+opt.width+1-pa;
      elif pb<>false then
        PrintToString(res,s{pos+[1..opt.width-pb+1]},"\n");
        pos:=pos+opt.width-pb+1;
      else Error("could not cut ",s{pos+[1..opt.width+1]});fi;
    od;
    PrintToString(res,s{[pos+1..Length(s)]},"\n");
  od;
  if IsBound(opt.file) then AppendTo(opt.file,res);else Print(res);fi;
end;
