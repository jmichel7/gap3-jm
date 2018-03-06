##############################################################################
##
#A  polyroot.g       VKCURVE package         Jean Michel
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds various root-finding functions for polynomials.
## 
#############################################################################
# p=vector-polynomial of complex rationals
# z=initial guess (complex rational)
# precision=rational (asked-for precision)
# returns [zero,guaranteed precision]
# where guaranteed precision <= asked-for precision (but may be much better)
# and where zero is a rational complex with 'as rounded denominator as possible'
# Returns false if does not converge after VKCURVE.NewtonLim iterations
NewtonRoot:=function(p,z,precision)local a,b,c,err,d,cnt,deriv,pr;
  deriv:=Derivative(p);
  for cnt in [1..VKCURVE.NewtonLim] do
    a:=ValuePol(p,z); b:=ValuePol(deriv,z);
    if b=0*b then c:=a;else c:=a*b^-1;fi;
    err:=AbsInt(BigNorm(c)); # AbsInt actually works for Rationals
    if err=0 then err:=precision/100/Length(p);fi;
    if err>precision then err:=precision;fi;
    pr:=Maximum(0,-1-DecimalLog(err)); # 10^(-pr-1)<=err<=10^-pr
    z:=ComplexRational(evalf(z-c,pr+1));
    if VKCURVE.showallnewton then Print(cnt,": ",evalf(z,pr),"\n");fi;
    if 10^(-pr)*(Length(p)-1)<=precision then 
    # note that |x-evalf(x,pr)|<=10^(-pr-1)
      if VKCURVE.showNewton then Print(cnt,".\c");fi;
      return [z,10^(-pr)];
    fi;
  od;
  if VKCURVE.showNewton then 
    Print("\n\n****** Non-Convergent Newton after ",VKCURVE.NewtonLim,
	  " iterations ******\n\n");
    Print("p=",ValuePol(p,Mvp("x"))," initial=",z," prec=",precision,"\n");
  fi;
  return false;
end;

SeparateRootsInitialGuess:=function(p,v,safety)local radv;
  if Length(p)=2 then return [-p[1]/p[2]];fi;
  radv:=Rational(Dispersal(v)[1]/safety/2);
  v:=List(v,e->NewtonRoot(p,ComplexRational(e),radv));
  if ForAny(v,x->x=false) or
   Rational(Dispersal(List(v,x->x[1]))[1]/safety/2)<Maximum(List(v,x->x[2]))
  then return false;
  else return List(v,x->x[1]);
  fi;
end;

# p is a vector of rationals or complex rationals representing a polynomial
# or an univariate Mvp with rational or complex rational coeffs
# it is assumed p has no multiple roots
#
# output: 
#  if d=1/2 min. dist. between 2 elts of output then it is guaranteed
#  that for any z in output there is a zero of p within radius d/safety
#
SeparateRoots:=function(p,safety)local v,subtractroot,e,cnt;
  subtractroot:=function(p,r)local d,res,i;
    d:=Length(p);
    res:=[];res[d-1]:=p[d];
    for i in [d-1,d-2..2] do res[i-1]:=p[i]+res[i]*r;od;
    return res;
  end;
  if IsMvp(p) then v:=Variables(p);
    if Length(v)>1 then Error(p," is not univariate");
    elif Length(v)=0 then v:=["x"];
    fi;
    p:=ScalMvp(Coefficients(p,v[1]));
  fi;
  p:=List(p,Complex);
  if Length(p)<=1 then return [];
  elif Length(p)=2 then return [-p[1]/p[2]];fi;
  p:=p/p[Length(p)];
  e:=evalf(E(7),10);v:=false;cnt:=0;
  while v=false and cnt<2*Length(p) do
    if VKCURVE.showNewton and cnt>0 then 
      Print("****** ",cnt," guess failed for p degree ",Length(p)-1,"\n");
    fi;
    v:=NewtonRoot(p,ComplexRational(e),
            10^(-(Length(p)+3+DecimalLog(safety))));
    e:=e*evalf((5/4)*E(2*Length(p)),10);
    cnt:=cnt+1;
  od;
  if cnt>=2*Length(p) then Error("no good initial guess");fi;
  v:=[v[1]];
  Append(v,SeparateRoots(subtractroot(p,ComplexRational(v[1])),safety));
  if safety=0 then return v;
  else return SeparateRootsInitialGuess(p,v,safety);
  fi;
end;

FindRoots:=function(p,prec)local v,subtractroot,e;
  subtractroot:=function(p,r)local d,res,i;
    d:=Length(p);
    res:=[];res[d-1]:=p[d];
    for i in [d-1,d-2..2] do res[i-1]:=p[i]+res[i]*r;od;
    return res;
  end;
  if IsMvp(p) then v:=Variables(p);
    if Length(v)>1 then Error(p," is not univariate");
    elif Length(v)=0 then v:=["x"];
    fi;
    p:=ScalMvp(Coefficients(p,v[1]));
  fi;
  if Length(p)=2 then return [-p[1]/p[2]];fi;
  p:=List(p,ComplexRational);
  e:=evalf(E(7),10);v:=false;
  while v=false do
    v:=NewtonRoot(p,ComplexRational(e),10^(-Length(p)-3))[1];
    e:=e*evalf(E(Length(p)),10);
  od;
  v:=Concatenation([v],FindRoots(subtractroot(p,ComplexRational(v)),prec));
  return List(v,e->NewtonRoot(p,ComplexRational(e),prec)[1]);
end;

# p=vector-polynomial of complex decimals
# initial guess=complex decimals
# pr=asked-for precision is 10^-pr
# returns [zero precise to its decimal precision]
# where zero's precision >= asked-for precision (but may be much better)
DecimalNewtonRoot:=function(p,z,pr)local a,b,c,err,d,cnt,deriv,fin,test,pp;
  d:=Length(p)-1; deriv:=Derivative(p);
  pp:=Minimum(List(p,x->x.r.precision));pp:=Minimum(pp,z.r.precision);
  test:=function()
    if ComplexOps.Norm(z)^d/pp>10^(-pr) then
      Error("unsufficient precision on coeffs");
    fi;
  end;
  fin:=function() test();if VKCURVE.showNewton then Print(cnt,"/\c");fi;end;
  test();
  for cnt in [1..VKCURVE.NewtonLim] do
    a:=ValuePol(p,z); b:=ValuePol(deriv,z);
    if b=0*b then 
       if a<>0*a then Error("unsufficient precision on coeffs");fi;
       fin();return z;
    fi;
    c:=a*ComplexConjugate(b);
    b:=ComplexOps.Norm(b);
    if b=0*b then 
      if c<>0*c then Error("unsufficient precision on coeffs");fi;
      fin();return z;
    fi;
    c:=Complex(c.r/b,c.i/b); 
    err:=AbsInt(BigNorm(c)); # AbsInt actually works for Rationals
    if err=0*err then fin();return z-c;fi;
    if err*d<=10^(-pr) then fin();return z-c;fi;
    z:=z-c;
    if VKCURVE.showallnewton then Print(cnt,": ",z,"\n");fi;
  od;
  if VKCURVE.showNewton then 
    Print("\n\n****** Non-Convergent Newton after ",VKCURVE.NewtonLim,
	  " iterations ******\n\n");
    Print("p=",ValuePol(p,Mvp("x"))," initial=",z," prec=",pr,"\n");
  fi;
  return [false];
end;

# p is a vector of decimals or complex decimals representing a polynomial
# or an Mvp in x with decimal or complex decimal coeffs
# pr is sought precision
#
# output: 
#  roots with desired precision
#  or error if precision on coeffs was unsufficient
#
DecimalRoots:=function(p,pr)local v,subtractroot,e,pp;
  subtractroot:=function(p,r)local d,res,i;
    d:=Length(p);
    res:=[];res[d-1]:=p[d];
    for i in [d-1,d-2..2] do res[i-1]:=p[i]+res[i]*r;od;
    return res;
  end;
  if IsMvp(p) then p:=ScalMvp(Coefficients(p,"x"));fi;
  p:=List(p,Complex);
  if Length(p)=2 then return [-p[1]/p[2]];fi;
  pp:=Minimum(List(p,x->x.r.precision));
  e:=evalf(E(7),-DecimalLog(1/pp));v:=false;
  while v=false do
    v:=DecimalNewtonRoot(p,e,Length(p)+3);
    e:=e*evalf(E(Length(p)));
  od;
  v:=Concatenation([v],DecimalRoots(subtractroot(p,v),pr));
  return List(v,e->DecimalNewtonRoot(p,e,pr));
end;
