SymOps:=OperationsRecord("SymOps");

IsSym:=x->IsRec(x) and IsBound(x.operations) and SymOps=x.operations;

Sym:=function(arg)local v,n;
  if Length(arg)=3 then 
    return rec(coeff:=arg[1],elm:=arg[2],nbvars:=arg[3],operations:=SymOps);
  fi;
  v:=arg[1];if Length(arg)=2 then n:=arg[2];else n:=Length(v);fi;
  Sort(v);return rec(coeff:=[1],elm:=[v],nbvars:=n,operations:=SymOps);
end;

SymOps.Normalize:=function(t)local i,elm,coeff;
  SortParallel(t.elm,t.coeff);
  elm:=[];coeff:=[];
  for i in [1..Length(t.elm)] do
    if i=1 or t.elm[i]<>t.elm[i-1] then
	 Add(elm,t.elm[i]);Add(coeff,t.coeff[i]);
    else coeff[Length(coeff)]:=coeff[Length(coeff)]+t.coeff[i];
    fi;
  od;
  i:=List(coeff,x->x<>0*x);
  t.coeff:=ListBlist(coeff,i); t.elm:=ListBlist(elm,i);
  return t;
end;

SymOps.Merge:=function(a,b)local i,j,r,c;
  i:=1;j:=1;
  if a.nbvars<>b.nbvars then Error("not same number of variables");fi;
  r:=ShallowCopy(a);r.elm:=[];r.coeff:=[];
  while i<=Length(a.elm) or j<=Length(b.elm) do
    if j>Length(b.elm) then
      Add(r.elm,a.elm[i]);Add(r.coeff,a.coeff[i]);i:=i+1;
    elif i>Length(a.elm) then
      Add(r.elm,b.elm[j]);Add(r.coeff,b.coeff[j]);j:=j+1;
    elif a.elm[i]<b.elm[j] then
      Add(r.elm,a.elm[i]);Add(r.coeff,a.coeff[i]);i:=i+1;
    elif a.elm[i]>b.elm[j] then
      Add(r.elm,b.elm[j]);Add(r.coeff,b.coeff[j]);j:=j+1;
    else
      c:=a.coeff[i]+b.coeff[j];
      if c<>0*c then Add(r.elm,a.elm[i]);Add(r.coeff,a.coeff[i]+b.coeff[j]);fi;
      i:=i+1;j:=j+1;
    fi;
  od;
  return r;
end;

SymOps.\+:=function(x,y)local tmp;
  if IsSym(y) then
    if IsList(x) then return List(x,z->z+y);
    elif IsSym(x) then return SymOps.Merge(x,y);
    fi;
  elif IsList(y) then return List(y,z->x+z);
  fi;
end;

SymOps.\-:=function(x,y)return x+(-1)*y;end;

SymOps.\*:=function(x,y)local res,a,i,j,fill,strip,cnt,vy,vx,cy;
  if IsSym(y) then
    if IsList(x) then return List(x,z->z*y);
    elif IsSym(x) then
      fill:=v->Concatenation(v,[Length(v)+1..x.nbvars]*0);
      strip:=function(v)Sort(v);v:=Filtered(v,x->x<>0);return v;end;
      cnt:=v->NrArrangements(v,Length(v));
      res:=Sym([],[],x.nbvars);
      vy:=List(y.elm,fill);cy:=List(vy,cnt);vx:=List(x.elm,fill);
      for i in [1..Length(vx)] do
        for a in Arrangements(vx[i],x.nbvars) do
	  for j in [1..Length(vy)] do
#    Print("v=",v," cnt(v)=",cnt(v)," cy[j]=",cy[j],"\n");
	    Add(res.coeff,x.coeff[i]*y.coeff[j]*cy[j]);
	    Add(res.elm,strip(a+vy[j]));
	  od;
	od;
      od;
      SymOps.Normalize(res);
      for i in [1..Length(res.elm)] do 
        res.coeff[i]:=res.coeff[i]/cnt(fill(res.elm[i]));
      od;
      return res;
    else # assume x is a 'scalar'
      if x=0*x then return Sym([],[],y.nbvars);
      else return Sym(y.coeff*x,y.elm,y.nbvars);
      fi;
    fi;
  elif IsList(y) then return List(y,z->x*z);
  else # assume y is a 'scalar'
    if y=0*y then return Sym([],[],x.nbvars);
    else return Sym(x.coeff*y,x.elm,x.nbvars);
    fi;
  fi;
end;

SymOps.\/:=function(x,y)
  x:=ShallowCopy(x);x.coeff:=x.coeff/y;return x;
end;

SymOps.\^:=function(f,n) local p;
  if n<0 then Error("negative power");
  elif n=0 then return Sym([],f.nbvars);
  else
   p:=false;
   while n>0 do
    if n mod 2 <> 0 then 
      if p=false then p:=f;
      else p:=p*f;
      fi;
    fi;
    f:=f*f;
    n:=QuoInt(n,2);
   od;
   return p;
   fi;
end;

SymOps.Format:=function(h,option)local i,res,e;
  if option in ["GAP","Maple"] then e:=List(h.elm,x->SPrint("(",Join(x),")"));
  else e:=List(h.elm,IntListToString);
  fi;
  e:=List([1..Length(h.elm)],i->FormatCoefficient(h.coeff[i],
      SPrint("m",e[i]),option));
  res:="";
  for i in [1..Length(h.elm)] do
    if Length(res)>0 and e[i][1]<>'-' then Append(res,"+");fi;
    Append(res,e[i]);
  od;
  if res="" then res:="0";fi;
  return String(res);
end;

SymOps.String:=h->SymOps.Format(h,"");

SymOps.Print:=function(o)Print(String(o));end;

# This implements the monomial symmetric functions
SymOps.Value:=function(p,v)local i,res,e,n;
  res:=0; n:=Length(v); # number of variables
  for i in [1..Length(p.elm)] do
#  Ptime(SPrint("start S",IntListToString(v)));
    e:=Concatenation(p.elm[i],[Length(p.elm[i])+1..n]*0);
    res:=res+p.coeff[i]*Sum(Arrangements(e,n),x->Product([1..n],i->v[i]^x[i]));
#  Ptime(SPrint("end S",IntListToString(v)));
  od;
  return res;
end;

# AsSymm(p) determine si le Mvp p est une somme de fonctions symetriques.
# AsSymm(p,vars) le fait par rapport aux variables de nom vars
AsSymm:=function(arg)local p,vars,res,i,j,v,m,good,bad;
  p:=arg[1];
  if Length(arg)=2 then vars:=arg[2];
  else vars:=Set(Concatenation(List(p.elm,x->x.elm)));
  fi;
  res:=rec(coeff:=[],elm:=[],nbvars:=Length(vars),operations:=SymOps);
  while Length(p.elm)>0 do
    m:=p.elm[1].elm;
    good:=Filtered([1..Length(m)],i->m[i] in vars);
    bad:=Filtered([1..Length(m)],i->not m[i] in vars);
    v:=p.elm[1].coeff{good};Sort(v);
    if v in res.elm then Error("not symmetric");fi;
    m:=Product(bad,i->Mvp(m[i])^p.elm[1].coeff[i]);
    res:=res+p.coeff[1]*m*Sym(v,Length(vars));
    p:=p-p.coeff[1]*Value(Sym(v,Length(vars)),List(vars,Mvp))*m;
  od;
  return res;
end;
