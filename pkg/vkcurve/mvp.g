##############################################################################
##
#A  mvp.g       VKCURVE package         Jean Michel
##  Multivariate Laurent and Puiseux polynomials
##
#Y  Copyright (C) 1997 - 2002  University Paris VII, France.
##
##  The  first version  was written  in march  1997 during  my visit at the
##  Newton institute.
## 
#############################################################################
#
# example of use:
#  
#  gap> x1:=Mvp("x1");x2:=Mvp("x2");
#  gap>  (x1+1/x2)^5;
#  5x1x2^-4+10x1^2x2^-3+10x1^3x2^-2+5x1^4x2^-1+x1^5+x2^-5
#  gap>  Format(last,rec(GAP:=1));
#  "5*x1*x2^-4+10*x1^2*x2^-3+10*x1^3*x2^-2+5*x1^4*x2^-1+x1^5+x2^-5"
#  
#  Internal representation: 
#     rec(elm:=vector of monomials,
#         coeff:=vector of corresponding coefficients,
#         operations:=MvpOps)
#  
#  where a monomial itself is a rec(elm:=vector of strings (variables), 
#                                   coeff:=vector of corresponding powers)
#  in a monomial or an Mvp .elm is sorted.
VKCURVE.mvp:=1;

MvpOps:=OperationsRecord("MvpOps");

IsMvp:=x->IsRec(x) and IsBound(x.operations) and MvpOps=x.operations;
MvpOps.Mvp:=function(e,c)return rec(elm:=e,coeff:=c,operations:=MvpOps);end;

# syntax
# Mvp("var"[,coeffs[,valuation]]])
# Mvp(elm,coeff)
# Mvp(polynomial)
# Mvp(Cyclotomic)
Mvp:=function(arg)local x,p,l;
  p:=arg[1];
  if IsString(p) and Length(p)>0 then 
    if Length(arg)=1 then
      return MvpOps.Mvp([rec(elm:=[p],coeff:=[1])],[1]);
    else l:=Filtered([1..Length(arg[2])],i->arg[2][i]<>0*arg[2][i]);
      if Length(arg)=3 then x:=arg[3];else x:=0;fi;
      return MvpOps.Mvp(List(l+x-1,function(i)
        if i=0 then return rec(coeff:=[],elm:=[]);
        else return rec(coeff:=[i],elm:=[p]);fi;end),arg[2]{l});
    fi;
  elif IsList(p) then 
    p:=MvpOps.Mvp(List(p,function(e)CollectCoefficients(e);return e;end),arg[2]);
    CollectCoefficients(p);return p;
  elif IsPolynomial(p) then
    x:=Indeterminate(p.baseRing);
    if not IsBound(x.name) then Error(x," should have .name bound");fi;
    return Mvp(x.name,p.coefficients,p.valuation);
  elif IsMvp(p) then return p;
  elif IsRatFrac(p) then
    if Length(p.den.coeff)=1 then return p.num/p.den;
    else return false;
    fi;
  elif p=0*p then return MvpOps.Mvp([],[]);
  else return MvpOps.Mvp([rec(elm:=[],coeff:=[])],[p]);
    # assume p is a non-zero 'scalar'
  fi;
end;

x:=Mvp("x");y:=Mvp("y");

# Assumes a and b have .elm=Set(.elm);  the result then also.
MvpOps.Merge:=function(a,b)local i,j,r,c,la,lb;
  i:=1;j:=1;r:=rec(elm:=[],coeff:=[]);
  la:=Length(a.elm);lb:=Length(b.elm);
  while i<=la or j<=lb do
    if j>lb then Add(r.elm,a.elm[i]);Add(r.coeff,a.coeff[i]);i:=i+1;
    elif i>la then Add(r.elm,b.elm[j]);Add(r.coeff,b.coeff[j]);j:=j+1;
    elif a.elm[i]<b.elm[j] then
      Add(r.elm,a.elm[i]);Add(r.coeff,a.coeff[i]);i:=i+1;
    elif a.elm[i]>b.elm[j] then
      Add(r.elm,b.elm[j]);Add(r.coeff,b.coeff[j]);j:=j+1;
    else
      c:=a.coeff[i]+b.coeff[j];
      if c<>0*c then Add(r.elm,a.elm[i]);Add(r.coeff,c);fi;
      i:=i+1;j:=j+1;
    fi;
  od;
  return r;
end;

MvpOps.\+:=function(x,y)local tmp;
  if IsList(y) then return List(y,z->x+z);
  elif IsMvp(y) then
    if IsList(x) then return List(x,z->z+y);
    elif IsMvp(x) then ;
    elif IsRatFrac(x) then return x+RatFrac(y);
    else x:=Mvp(x);# assume x is a 'scalar'
    fi;
  else y:=Mvp(y);# assume y is a 'scalar'
  fi;
  tmp:=MvpOps.Merge(x,y); tmp.operations:=MvpOps; return tmp;
end;

MvpOps.\-:=function(x,y)return x+(-1)*y;end;

MvpOps.Format:=function(h,option)local i,GAP,elm,coeff;
  if h.elm=[] then return "0";fi;
  GAP:=IsBound(option.Maple) or IsBound(option.GAP);
  if IsBound(option.reverse) then 
    elm:=Reversed(h.elm); coeff:=Reversed(h.coeff);
  else  elm:=h.elm; coeff:=h.coeff;
  fi;
  h:=Zip(elm,coeff,function(e,c)
    e:=Zip(e.elm,e.coeff,function(e,c)local res;
      res:=FormatMonomial(e,c,option);
      if res="1" then return "";else return res;fi;
      end);
    if GAP then e:=Join(e,"*");else e:=Concatenation(e);fi;
    return FormatCoefficient(c,e,option);
    end);
  for i in [2..Length(h)] do if h[i][1]<>'-' then h[i]:=SPrint("+",h[i]);fi;od;
  return String(Concatenation(h));
end;

MvpOps.String:=h->MvpOps.Format(h,rec());

MvpOps.Print:=function(o)Print(String(o));end;

MvpOps.\*:=function(x,y)local res,i;
  if IsMvp(y) then
    if IsMvp(x) then res:=Mvp(0);
      for i in [1..Length(x.elm)] do
	Append(res.elm,List(y.elm,j->MvpOps.Merge(x.elm[i],j)));
	Append(res.coeff,x.coeff[i]*y.coeff);
      od;
      CollectCoefficients(res);return res;
    elif IsList(x) then return List(x,z->z*y);
    elif IsCyc(x) then # frequent case (e.g. see subtraction)
      if x=0 then return Mvp(0);else return MvpOps.Mvp(y.elm,y.coeff*x);fi;
    elif IsRatFrac(x) then return x*RatFrac(y);
    else # assume x is a 'scalar'
	if x=0*x then return Mvp(0);else return MvpOps.Mvp(y.elm,y.coeff*x);fi;
    fi;
  elif IsList(y) then return List(y,z->x*z);
  else # assume y is a 'scalar'
    if y=0*y then return Mvp(0);else return MvpOps.Mvp(x.elm,x.coeff*y);fi;
  fi;
end;

MvpOps.\^:=function(x,i)local y;
  if i=0 then 
    if Length(x.coeff)=0 then return 1;else return Mvp(x.coeff[1]^0);fi;
  elif Length(x.elm)=1 then
    if IsInt(i) then y:=x.coeff[1]^i;
    else y:=GetRoot(x.coeff[1],Denominator(i))^Numerator(i);
    fi;
    return Mvp([rec(elm:=x.elm[1].elm,coeff:=x.elm[1].coeff*i)],[y]);
  elif i>0 then y:=Mvp(1);
    while i>0 do
     if i mod 2 <> 0 then y:=y*x;fi;
     if i>=2 then x:=x*x;fi;
     i:=QuoInt(i,2);
    od;
    return y;
  else return (1/x)^-i;
  fi;
end;

MvpOps.GetRoot:=function(arg)local res,x,n,msg;
  x:=arg[1];n:=arg[2];msg:=arg{[3..Length(arg)]};
  if n=1 then return x;
  elif Length(x.elm)<>1 then
    ApplyFunc(Error,Concatenation(["unable to compute ",n,
      "-th root of non-monomial ",x,":\n"], msg));
    return false;
  fi;
  res:=Copy(x);
  res.elm[1].coeff:=res.elm[1].coeff/n;
  if res.coeff[1]<>res.coeff[1]^0 then 
    res.coeff[1]:=ApplyFunc(GetRoot,Concatenation([res.coeff[1],n],msg));
  fi;
  return res;
end;

# usage: Value(f,[var1,value1,var2,value2,...])
# means specialize each of var_i to value_i
MvpOps.Value:=function(arg)local x,f,res,i,m,elm,coeff,p,j,r,vars,values,d;
  f:=arg[1];
  if Length(arg)=2 then x:=arg[2];else x:=arg{[2..Length(arg)]};fi;
  res:=Mvp(0);
  vars:=x{[1,3..Length(x)-1]};values:=x{[2,4..Length(x)]};
  if Length(f.coeff)=0 then return f;fi;
  p:=List(vars,function(s)local pp,i,j;pp:=[];
    for i in [1..Length(f.elm)] do j:=Position(f.elm[i].elm,s);
      if j<>false then Add(pp,[i,j]);fi;
    od;
    return pp;end);
  r:=Filtered([1..Length(p)],i->Length(p[i])>0);
  vars:=vars{r};values:=values{r};p:=p{r};
  if p=[] then return f;fi;
  d:=List(p,pp->Lcm(List(pp,u->Denominator(f.elm[u[1]].coeff[u[2]]))));
  for i in [1..Length(vars)] do 
   if d[i]<>1 then values[i]:=GetRoot(values[i],d[i]);fi;
  od;
  for i in [1..Length(f.coeff)] do
    elm:=rec(elm:=[],coeff:=[]); 
    coeff:=f.coeff[i];
    m:=f.elm[i];
    for j in [1..Length(m.elm)] do
      p:=Position(vars,m.elm[j]);
      if p=false then Add(elm.elm,m.elm[j]);Add(elm.coeff,m.coeff[j]);
      else coeff:=coeff*values[p]^(d[p]*m.coeff[j]);
      fi;
    od;
    if IsMvp(coeff) then
      Append(res.elm,List(coeff.elm,x->MvpOps.Merge(x,elm)));
      Append(res.coeff,coeff.coeff);
    elif IsRatFrac(coeff) then return Value(RatFrac(f),x);
    else Add(res.elm,elm); Add(res.coeff,coeff);
    fi;
  od;
  CollectCoefficients(res);return res;
end;

MvpOps.EvalPolRoot:=function(pol,x,n,p)local l,ic,pc,i,j,r,root; 
  pol:=Value(pol,["foo",Mvp("foo")^(1/n)*p]);
  return Value(pol,["foo",x]);
end;

MvpOps.ComplexConjugate:=function(x)
  x:=ShallowCopy(x);x.coeff:=ComplexConjugate(x.coeff);return x;
end;

MvpOps.Degree:=function(arg)local p;
  p:=arg[1];
  if Length(p.elm)=0 then return -1;fi;
  if Length(arg)=1 then return Maximum(List(p.elm,y->Sum(y.coeff)));fi;
  return Maximum(List(p.elm,function(e)local pv;
    pv:=Position(e.elm,arg[2]);
    if pv=false then return 0;
    else return e.coeff[pv];fi;end));
end;

MvpOps.Valuation:=function(arg)local p;
  p:=arg[1];
  if Length(p.elm)=0 then return -1;fi;
  if Length(arg)=1 then return Minimum(List(p.elm,y->Sum(y.coeff)));fi;
  return Minimum(List(p.elm,function(e)local pv;
    pv:=Position(e.elm,arg[2]);
    if pv=false then return 0;
    else return e.coeff[pv];fi;end));
end;

MvpOps.Derivative:=function(arg)local x,variable,res,i,p,elm;
  x:=arg[1];
  if Length(arg)=2 then variable:=arg[2];
  else variable:=Variables(x)[1];
  fi;
  res:=Mvp(0);
  for i in [1..Length(x.coeff)] do
    p:=PositionProperty(x.elm[i].elm,z->z=variable);
    if p<>false then
      Add(res.coeff,x.coeff[i]*x.elm[i].coeff[p]);
      elm:=Copy(x.elm[i]);elm.coeff[p]:=elm.coeff[p]-1;
      CollectCoefficients(elm);
      Add(res.elm,elm);
    fi;
  od;
  CollectCoefficients(res);return res;
end;

MvpOps.Coefficients:=function(x,var)local res,i,p,elm,d;
  res:=[];
  i:=Variables(x);
  if i=[] then return x.coeff;
  elif Length(i)=1 and i[1]=var then
    for i in [1..Length(x.coeff)] do
      if Length(x.elm[i].elm)=0 then res[1]:=x.coeff[i];
      else d:=x.elm[i].coeff[1]+1;
        if d<=0 or not IsInt(d) then
          Error(x," is not a polynomial with respect to ",var,"\n");
	fi;
	res[d]:=x.coeff[i];
      fi;
    od;
    for i in [1..Length(res)] do if not IsBound(res[i]) then res[i]:=0;fi;od;
    return res;
  fi;
  for i in [1..Length(x.coeff)] do
    p:=Position(x.elm[i].elm,var);
    if p=false then  
      if not IsBound(res[1]) then res[1]:=Mvp(x.elm{[i]},x.coeff{[i]});
      else Add(res[1].coeff,x.coeff[i]);Add(res[1].elm,x.elm[i]);
      fi;
    else  
      elm:=ShallowCopy(x.elm[i]);
      d:=elm.coeff[p]+1;
      if d<=0 or not IsInt(d) then
        Error(x," is not a polynomial with respect to ",var,"\n");
      fi;
      if not IsBound(res[d]) then res[d]:=Mvp(0);fi; 
      p:=Concatenation([1..p-1],[p+1..Length(elm.coeff)]);
      elm.coeff:=elm.coeff{p};elm.elm:=elm.elm{p};
      Add(res[d].coeff,x.coeff[i]);Add(res[d].elm,elm);
    fi;
  od;
  for d in [1..Length(res)] do 
    if not IsBound(res[d]) then res[d]:=Mvp(0);
    else SortParallel(res[d].elm,res[d].coeff);
    fi;
  od;
  return res;
end;

MvpOps.LeadingCoefficient:=function(x,var)local c;
  c:=Coefficients(x,var);
  return c[Length(c)];
end;

MvpOps.Variables:=p->Union(List(p.elm,x->x.elm));

MvpOps.ExactDiv:=function(p,q)local pos,var,res,v,w,lv,lw,listvar,i;
  # Print("p=",p," q=",q,"\n");
  if Length(q.coeff)=1 then return p*q^-1;
  elif Length(q.coeff)=0 then Error("cannot divide by 0");
  elif Length(p.coeff)=0 then return p;
  fi;
  listvar:=Variables(p);
  if Length(listvar)=0 then return false;fi;
  var:=listvar[1];
  v:=Coefficients(p,var);w:=Coefficients(q,var);
  if Length(listvar)=1 then
    if ScalMvp(w)<>false then 
     v:=ScalMvp(v)*p.coeff[1]^0; w:=ScalMvp(w);
    else v:=v+w[1]*0;
    fi;
  fi;
  w:=w+v[1]*0;
  lv:=Length(v);lw:=Length(w);
  res:=[];i:=lv-lw+1;
  while i>0 do
  # Print("var=",var," v=",v," w=",w,"\n");
    if v[lv]=0*v[lv] then res[i]:=0*p.coeff[1];
    else
      if Length(listvar)=1 then res[i]:=v[lv]/w[lw];
      else res[i]:=MvpOps.ExactDiv(v[lv],w[lw]);
        if false=res[i] then return false;fi;
      fi;
      if lw>1 then v{[i..lv-1]}:=v{[i..lv-1]}-res[i]*w{[1..lw-1]};fi;
    fi;
    lv:=lv-1;i:=i-1;
  od;
  if ForAll([1..lv],i->v[i]=0*v[i]) then return Mvp(ValuePol(res,Mvp(var)));
  else return false;
  fi;
end;

# returns monomial holding minimum degrees in the list of .elm of an Mvp
MonomialGcd:=function(p)local coeff,elm,m,p,j,v;
  coeff:=[];elm:=Union(List(p,x->x.elm));
  for m in p do 
    for j in [1..Length(elm)] do
      v:=Position(m.elm,elm[j]);
      if v=false then v:=0; else v:=m.coeff[v];fi;
      if not IsBound(coeff[j]) or coeff[j]>v then coeff[j]:=v;fi;
    od; 
  od;
  v:=Filtered([1..Length(coeff)],i->coeff[i]<>0);
  return rec(coeff:=coeff{v},elm:=elm{v});
end;

MvpOps.\/:=function(a,b)local res,ma,mb;
  if IsMvp(b) then
    if IsList(a) then return List(a,z->z/b);
    elif IsRatFrac(a) then return a/RatFrac(b);
    elif IsMvp(a) then ;
    else # assume a is a 'scalar'
      a:=Mvp(a);
    fi;
  else return MvpOps.Mvp(a.elm,a.coeff/b); # assume b is a 'scalar'
  fi;
  if Length(b.elm)=1 then return a*b^-1;fi;
  ma:=Mvp([MonomialGcd(a.elm)],[1]);
  mb:=Mvp([MonomialGcd(b.elm)],[1]);
  res:=[a/ma,b/mb];
  res:=MvpOps.ExactDiv(res[1],res[2]);
  if res<>false then return res*ma/mb;fi;
  res:=RatFrac(a,b);
  a:=ScalMvp(res.den);
  if a<>false and a=a^0 then res:=res.num;fi;
  return res;
end;

MvpOps.evalf:=function(x,p)
  x:=ShallowCopy(x);x.coeff:=List(x.coeff,y->evalf(y,p));return x;
end;

# matrix H_i,j := d/dx_j(d/dx_i(p))
Hessian:=function(p,varnames)
  return List(varnames,x->List(varnames,y->Derivative(Derivative(p,x),y)));
end;

# matrix J_i,j := d/dx_j(pols[i])
Jacobian:=function(pols,varnames)
  return List(pols,i->List(varnames,j->Derivative(i,j)));
end;

# return false is Mvp is not a scalar, that scalar else
ScalMvp:=function(x)
  if IsList(x) then 
    x:=List(x,ScalMvp);
    if ForAny(x,y->y=false) then return false;
    else return x;
    fi;
  elif IsMvp(x) then 
    if Length(x.coeff)=0 then return 0;
    elif Length(x.coeff)=1 and Length(x.elm[1].elm)=0 then return x.coeff[1];
    else return false;
    fi;
  elif IsRatFrac(x) then return ScalMvp(Mvp(x));
  else return x;
  fi;
end;

# factorize a quadratic form (an Mvp of degree 2) as product of 2 linear forms
FactorizeQuadraticForm:=function(p)local v,r,m,i,t,e,n,b,d;
  v:=Variables(p);r:=Length(v)+1; m:=NullMat(r);
  for i in [1..Length(p.elm)] do
    t:=p.coeff[i];e:=p.elm[i];n:=List(e.elm,x->Position(v,x));e:=e.coeff;
    if e=[1,1] then m[n[1]][n[2]]:=t/2;m[n[2]][n[1]]:=t/2;
    elif e=[2] then m[n[1]][n[1]]:=t;
    elif e=[1] then m[n[1]][r]:=t/2; m[r][n[1]]:=t/2;
    elif e=[] then m[r][r]:=t;
    else InfoChevie("# not a quadratic form");return false;
    fi;
  od;
  if Length(m)=2 then t:=m^0;
  else n:=Copy(m);
    TriangulizeMat(m);m:=Filtered(m,x->x<>0*x);
    if Length(m)>2 then return false;fi;
    t:=List(n,x->SolutionMat(m,x));
    m:=List(m,x->SolutionMat(TransposedMat(t),x));
  fi;
  v:=Concatenation(List(v,Mvp),[Mvp("x")^0])*t;
  if Length(m)=1 then return [v[1],v[1]*m[1][1]];fi;
  b:=m[1][2]+m[2][1];
  if m[1][1]=0 then return [v[2],b*v[1]+m[2][2]*v[2]];fi;
  b:=b/m[1][1];
  d:=GetRoot(b^2-4*m[2][2]/m[1][1],2,"no");
  if d=false then return false;fi;
  return [v[1]+v[2]/2*(b-d),m[1][1]*(v[1]+v[2]/2*(b+d))];
end;

# by GAP's rules, called only if y is an Mvp or y a basic type and x Mvp
MvpOps.\=:=function(x,y)
# Print("calling ",x," MvpOps.= ",y,"\n");
  if IsCyc(y) then return y=ScalMvp(x);
  elif IsCyc(x) then return x=ScalMvp(y);
  elif IsMvp(x) then return IsMvp(y) and x.coeff=y.coeff and x.elm=y.elm;
  else return y=x;
  fi;
end;

MvpOps.CycPol:=function(p)local l;
  l:=Variables(p);
  if Length(l)=0 then return CycPol(ScalMvp(p));
  elif Length(l)=1 then return CycPol(ScalMvp(Value(p,[l[1],X(Cyclotomics)])));
  else Error("Mvp should be univariate");fi;
end;

MvpOps.Mod1:=function(p)local res,i;
  if ForAll(p.elm,x->Length(x.coeff)>0) then return p;fi;
  res:=ShallowCopy(p);res.coeff:=ShallowCopy(res.coeff);
  i:=PositionProperty(res.elm,x->Length(x.coeff)=0);
  res.coeff[i]:=Mod1(res.coeff[i]);
  return res;
end;

# action(m,p [,varnames]) matrix m acts on Mvp p [on variables varnames]
OnPolynomials:=function(arg)local m,p,varnames,vars;
  m:=arg[1];p:=arg[2];
  if Length(arg)=3 then varnames:=arg[3]{[1..Length(m)]};
  else varnames:=Variables(p);
  fi;
  return Value(p,Concatenation(TransposedMat([varnames,List(varnames,Mvp)*m])));
end;

MvpOps.Galois:=function(x,e)
  x:=ShallowCopy(x);x.coeff:=Galois(x.coeff,e);return x;
end;
