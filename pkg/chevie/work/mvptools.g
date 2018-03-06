# nxn matrix full of new variables starting with letter n
gmat:=function(n,letter)return List([1..n],i->List([1..n],j->
Mvp(SPrint(letter,i,"_",j))));end;

# convert an Mvp of the form ax+b to ["x",-b/a]
linear:=function(p)local i,var;
  if ForAny(p.elm,x->Sum(x.coeff)>1) then return false;fi;
  i:=PositionProperty(p.elm,x->Length(x.coeff)=1);
  if i=false then Error();fi;
  var:=p.elm[i].elm[1];
  return [var,Mvp(var)-p/p.coeff[i]];
end;

solvefor:=function(l,v)
  l:=Coefficients(l,v);
  if Length(l)>2 then return false;fi;
  return [v,-l[1]/l[2]];
end;

solve:=function(x)local var,v;
  for var in Variables(x) do
    v:=solvefor(x,var);
    if v<>false and (IsMvp(v[2]) or IsCyc(v[2])) then return v;fi;
  od;
  return false;
end;
  
pcomplexity:=p->Length(String(p));
pcomplexity:=p->Length(p.elm);

# sort a list of mvps by Length(String)
sortbylg:=function(arg)local n,f;
  n:=arg[1];
  if Length(arg)=2 then f:=arg[2];else f:=pcomplexity;fi;
  n:=Filtered(Set(n),x->x<>0*x);
  SortBy(n,f);
  return n;
end;

factor2:=function(p)local v,c,q; # factors 2nd-degree Mvps
  if Degree(p)<>2 then return false;fi;
  for v in Variables(p) do
    c:=Coefficients(p,v);
    if c[1]=0*c[1] then 
      Print("v=",v," c=",c,"\n");
      return [Mvp(v),p/Mvp(v)];fi;
    q:=MvpOps.ExactDiv(c[1],c[2]);
    if Degree(c[2])>0 and q<>false and ForAll(q.elm,x->ForAll(x.coeff,y->y>0)) then 
      Print("v=",v," c=",c,"\n");
      return [c[2],p/c[2]];
    fi;
  od;
  return false;
end;

# take in account linear factors of 2nd degrees mvp when shrinking n
simplify:=function(n)local i,p,f;
  f:=x->ScalMvp(x)<>false and ScalMvp(x)<>0;
  for i in [1..Length(n)] do p:=factor2(n[i]);
   if p<>false then 
     if ForAny(n,x->f(Value(x,linear(p[1])))) then
	n[i]:=n[i]/p[1];Print(p[1]," is simplifiable so n[",i,"]=>",n[i],"\n");
     elif ForAny(n,x->f(Value(x,linear(p[2])))) then
	n[i]:=n[i]/p[2];Print(p[2]," is simplifiable so n[",i,"]=>",n[i],"\n");
     else Print("n[",i,"] factors as ",p," but not simplifiable\n");
     fi;
   fi;
  od;
  return n;
end;

# shrint(<list n>[, silent])
shrint:=function(arg)local i,n,value,l,j,q,shrown,try,norm,m,p,coeff,pr;
  n:=arg[1];
  if Length(arg)=1 then pr:=Print;else pr:=Ignore;fi;
  norm:=x->x*ComplexConjugate(x);
  coeff:=function(p,v)
    p:=Coefficients(p,v);
    if Length(p)<>2 then return 0;
    else return p[2];fi;
  end;
  try:=function(i,j)local l,k,ni,nn;
    if i=0 then return false;fi;
    l:=Length(String(i));
    for k in Variables(j) do
      if coeff(j,k)<>0 then
        ni:=ScalMvp(coeff(i,k)/coeff(j,k));
	if ni<>false and IsInt(ni) then
	  nn:=i-ni*j;
	  if Length(Variables(nn))<Length(Variables(i)) then 
#    pr(i," - ",ni," * ",j,"=>",nn,"\n");
	    return nn;fi;
	fi;
      fi;
    od;
    nn:=evalf(norm(i.coeff));
    for k in [-2..2] do
      ni:=i+k*j;
      if ni=0 or ((Length(String(ni))<l or 
         (Length(String(ni))=l and evalf(norm(ni.coeff))<nn))
	  and Length(Variables(ni))<=Length(Variables(i))) then return ni;fi;
    od;
    return false;
  end;
  pr(Length(n)," mvps of total length=",Sum(n,x->Length(String(x))),"\n");
  n:=sortbylg(n);
  for i in [1..Length(n)] do 
    pr(i,":",Length(String(n[i])));
    if Length(String(-n[i]))<Length(String(n[i])) then 
      n[i]:=-n[i];
      pr("=>",Length(String(n[i])),"\c");
    fi;
    p:=PositionProperty(n[i].elm,x->x.coeff=[]);
    if p<>false and IsRat(n[i].coeff[p]) then
      n[i]:=Copy(n[i]);n[i].coeff[p]:=n[i].coeff[p]-Int(n[i].coeff[p]);
      CollectCoefficients(n[i]);
    fi;
    m:=Minimum(Length(n),i+5);
    for j in [m,m-1..1] do
      if j<>i then
      repeat
	q:=try(n[i],n[j]);
	if q<>false then
	  n[i]:=q; pr("=>",Length(String(n[i]))," by ",j,"\c");
	fi;
      until q=false;
      fi;
    od;
    pr("\n");
  od;
  n:=sortbylg(n);
  pr(Length(n)," mvps of total length=",Sum(n,x->Length(String(x))),"\n");
  return n;
end;

# shr(l[,lim])
shr:=function(arg)local n,i,value,l,j,q,lim,shrown,try;
  try:=function(i,j)local q,l,k,ni,p;
    q:=Filtered([Length(i.elm),Length(i.elm)-1..1],
       x->Position(j.elm,i.elm[x])<>false);
    l:=pcomplexity(i);
    for k in q do
      p:=Position(j.elm,i.elm[k]);
      ni:=i-j*i.coeff[k]/j.coeff[p];
      if pcomplexity(ni)<l or ni=0 then return ni;fi;
      if pcomplexity(ni/ni.coeff[Length(ni.coeff)])<l then return 
	ni/ni.coeff[Length(ni.coeff)];fi;
    od;
    return false;
  end;
  shrown:=function(i)local p;
    repeat
      p:=PositionProperty(n[i].coeff,x->l>pcomplexity(n[i]/x));
      if p<>false then
	n[i]:=n[i]/n[i].coeff[p];
	l:=pcomplexity(n[i]); Print("=>",l,"\c");
      fi;
    until p=false;
  end;
  n:=arg[1];
  Print(Length(n)," mvps of total length=",Sum(n,pcomplexity),"\n");
  n:=sortbylg(n);
  for i in [1..Length(n)] do 
    l:=pcomplexity(n[i]);
    Print(i,":",l);shrown(i);
    if Length(arg)>1 then lim:=Maximum(1,i-arg[2]);else lim:=1;fi;
    for j in [i-1,i-2..lim] do
      repeat
	q:=try(n[i],n[j]);
	if q<>false then
	  n[i]:=q; l:=pcomplexity(n[i]);Print("=>",l," by ",j,"\c");
	  shrown(i);
	fi;
      until q=false;
    od;
    Print("\n");
  od;
#  n:=simplify(n);
  n:=sortbylg(n);
  Print(Length(n)," mvps of total length=",Sum(n,pcomplexity),"\n");
  return n;
end;

# for mvp2
#shr2:=function(arg)local n,i,value,l,j,q,lim,shrown,try;
#  try:=function(p,q)local l,lc,lp,i,nl,qt,v,np;
#    if q=0 then return false;fi;
#    l:=q.elm[Length(q.elm)]; #leading term
#    lc:=q.coeff[Length(q.elm)]; #leading coeff
#    lp:=pcomplexity(p);
#    i:=Length(p.elm);
#    while i>0 do
#      if l>p.elm[i] then return false;fi;
#      nl:=MonomialProd(p.elm[i],MonomialPower(l,-1));
#      if ForAll(RecFields(nl),f->nl.(f)>=0) then
#	qt:=MvpNC([nl],[p.coeff[i]/lc]);
#	v:=MvpOps.MulMonom(q,nl,-p.coeff[i]/lc);
#	np:=p+v;
#        if pcomplexity(np)<lp then return np;fi;
#      fi;
#      i:=i-1;
#    od;
#    return false;
#  end;
#  shrown:=function(i)local p;
#    repeat
#      p:=PositionProperty(n[i].coeff,x->l>pcomplexity(n[i]/x));
#      if p<>false then
#	n[i]:=n[i]/n[i].coeff[p];
#	l:=pcomplexity(n[i]); Print("=>",l,"\c");
#      fi;
#    until p=false;
#  end;
#  n:=arg[1];
#  Print(Length(n)," mvps of total length=",Sum(n,pcomplexity),"\n");
#  n:=sortbylg(n);
#  for i in [1..Length(n)] do 
#    l:=pcomplexity(n[i]);
#    Print(i,":",l);shrown(i);
#    if Length(arg)>1 then lim:=Maximum(1,i-arg[2]);else lim:=1;fi;
#    for j in [i-1,i-2..lim] do
#      repeat
#	q:=try(n[i],n[j]);
#	if q<>false then
#	  n[i]:=q; l:=pcomplexity(n[i]);Print("=>",l," by ",j,"\c");
#	  shrown(i);
#	fi;
#      until q=false;
#    od;
#    Print("\n");
#  od;
##  n:=simplify(n);
#  n:=sortbylg(n);
#  Print(Length(n)," mvps of total length=",Sum(n,pcomplexity),"\n");
#  return n;
#end;

sshr:=function(l)local n,new;
  new:=Sum(l,x->Length(String(x)));
  repeat
    n:=new;
    l:=shr(l);
    new:=Sum(l,x->Length(String(x)));
  until new=n;
  return l;
end;

# find a better basis of ideal generated by Mvps t
elimMvps:=function(t)local m;
 m:=Set(Concatenation(List(t,x->x.elm)));
 Sort(m,function(a,b)return Sum(a.coeff)<Sum(b.coeff)
  or (Sum(a.coeff)=Sum(b.coeff) and a<b);end);
 t:=List(t,function(p)local v,i;v:=List(m,x->0);
   for i in [1..Length(p.elm)] do v[Position(m,p.elm[i])]:=p.coeff[i];od;
   return v;end);
 TriangulizeMat(t);
 return Filtered(List(t,v->Mvp(m,v)),x->not x=0*x);
end;

# RenameVariables(p) rename vars to "a".."z","A".."Z","a1" etc...
# RenameVariables(p,v) rename Variables(p)[i] to v[i]
# RenameVariables(p,s,v) rename in p s[i] to v[i]
# can be applied to list or matrix to rename all variables simultaneously
RenameVariables:=function(arg)local p,source,target;
  p:=arg[1];
  if Length(arg)<=2 then source:=Variables(p); else source:=arg[2];fi;
  if Length(arg)=1 then target:=List([1..Length(source)]-1,
      function(i)local alp;
      alp:="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
      if i<52 then return String([alp[i+1]]);
      else return SPrint([alp[(i mod 52)+1]],QuoInt(i,52));
      fi;
    end);
  else target:=arg[Length(arg)];
  fi;
  if IsList(p) then return List(p,x->RenameVariables(x,source,target));fi;
  p:=ShallowCopy(p);
  p.elm:=List(p.elm,function(r)r:=ShallowCopy(r);r.elm:=List(r.elm,
   function(s)local p;p:=Position(source,s);
    if p<>false then return target[p];else return s;fi;end);
   return r;end);
  return p;
end;

mRatFrac:=function(m)if IsList(m) then return List(m,mRatFrac);
                     else return RatFrac(m);fi;end;

mMvp:=function(m)local v;if IsList(m) then return List(m,mMvp);
  else v:=Mvp(m); if v<>false then return v;else return m;fi;fi;end;

# in x weed out a factor v (string or mvp)
clean:=function(x,v)
  if IsList(x) then return List(x,y->clean(y,v));
  elif IsString(v) and Value(x,[v,0])=0 then 
    Print(".\c"); return x/Mvp(v);
  elif IsMvp(v) and MvpOps.ExactDiv(x,v)<>false then
    Print(".\c"); return MvpOps.ExactDiv(x,v);
  else return x;fi;
end;

cleanlist:=l->l{Filtered([1..Length(l)],i->not ForAny([1..i-1],j->
    MvpOps.ExactDiv(l[i],l[j])<>false))};

# returns a monomial representing the degree of p
MvpDegree:=function(P)local v,va,e,i,p;
  v:=Variables(P);va:=List(v,x->0);
  for e in P.elm do
    for i in [1..Length(v)] do
      p:=Position(e.elm,v[i]);
      if p<>false then va[i]:=Maximum(va[i],e.coeff[p]);fi;
    od;
  od;
  v:=ListBlist(v,List(va,x->x<>0));
  va:=ListBlist(va,List(va,x->x<>0));
  return rec(coeff:=[1],elm:=[rec(coeff:=va,elm:=v)],operations:=MvpOps);
end;

# cut list l in factors until no l[i] divides l[j]
weedfacts:=function(l)local i,v,sum;
  for i in [1..Length(l)] do 
    v:=Mvp([MonomialGcd(l[i].elm)],[1]);if v<>1 then l[i]:=l[i]/v;fi;
  od;
  l:=Filtered(l,x->Degree(x)>0);
  l:=List(l,x->x/x.coeff[1]);
  l:=Set(l);SortBy(l,Degree);
  repeat
  sum:=0;
  for i in [1..Length(l)-1] do
    v:=Filtered([i+1..Length(l)],j->MvpOps.ExactDiv(l[j],l[i])<>false);
    l{v}:=List(l{v}/l[i],x->x/x.coeff[1]);sum:=sum+Length(v);
  od;
  l:=Filtered(l,x->Degree(x)>0);l:=Set(l);SortBy(l,Degree);
  until sum=0;
  SortBy(l,x->[Degree(x),Length(Format(x))]);
  return l;
end;

# factorize p using factors in l
listfactor:=function(p,l)local res,q;
  res:=[];
  q:=Mvp([MonomialGcd(p.elm)],[1]);if q<>1 then p:=p/q;Add(res,q);fi;
  if p.coeff[1]<>1 then Add(res,p.coeff[1]*Mvp("q")^0);p:=p/p.coeff[1];fi;
  for q in l do
    while MvpOps.ExactDiv(p,q)<>false do Add(res,q);p:=p/q;od;
  od;
# if p<>1 then Error();Add(res,p);fi;
  SortBy(res,Degree);
  return res;
end;

# if 2nd arg given do not rescale variables
shrinkpol:=function(arg)local pol,varnames,den,num,primes,l,p,p1,l1,v;
  pol:=arg[1];
  varnames:=Set(Concatenation(List(pol.elm,x->x.elm)));
  den:=List(pol.coeff,Denominator);
  num:=List(pol.coeff,Numerator);
  primes:=Set(List(Concatenation(List(Concatenation(den,num),Factors)),AbsInt));
  l:=Length(String(pol));
  Print("initial length=",l,"\n");
  for p in primes do
    while Length(String(pol/p))<l do
      l:=Length(String(pol/p));
      Print("p->p/",p," => length ",l,"\n"); pol:=pol/p;
    od;
    while Length(String(pol*p))<l do
      l:=Length(String(pol*p));
      Print("p->p*",p," => length ",l,"\n"); pol:=pol*p;
    od;
  od;    
  if Length(arg)<2 then
  for v in varnames do
    for p in primes do
      p1:=Value(pol,[v,Mvp(v)/p]);l1:=Length(String(p1));
      if l1<l then
	Print("p->p(",v,"/",p,") => length ",l1,"\n"); pol:=p1;l:=l1;
      fi;
      p1:=Value(pol,[v,Mvp(v)*p]);l1:=Length(String(p1));
      if l1<l then
	Print("p->p(",v,"*",p,") => length ",l1,"\n"); pol:=p1;l:=l1;
      fi;
    od;
  od;
  fi;
  return pol;
end;

# montre les coefficients d'un Mvp factorises
analpol:=function(p)local m,w,pcoll,pel,vars;

 vars:=Set(Concatenation(List(p.elm,x->x.elm)));

 pcoll:=c->Join(List(Collected(Factors(AbsInt(c))),
   p->FormatMonomial(p[1],p[2],rec()))," ");

 pel:=function(e)local res,i;
  res:=List(vars,x->0);
  for i in [1..Length(e.elm)] do res[Position(vars,e.elm[i])]:=e.coeff[i]; od;
  return res;
 end;

  Print("length=",Length(String(p)),"\n");
  m:=List(p.coeff,x->[Numerator(x),Denominator(x)]);
  m:=List(m,x->List(x,pcoll));
  w:=List(p.elm,pel);
  m:=Concatenation(TransposedMat(w),TransposedMat(m));
  m:=TransposedMat(m);
  SortParallel(List(w,Sum),m);
  m:=Concatenation([Concatenation(vars,["num","den"])],m);
  Print(FormatTable(m,rec(rowLabels:=[1..Length(m)])));
end;

#try rescale any var of m by p
rescale:=function(m,p)local l,var,nl,nm;
  l:=Length(Format(m));
  for var in Variables(m) do
    nm:=Value(m,[var,Mvp(var)*p]);
    nl:=Length(Format(nm));
    if nl<l then 
      Print(var,"->",var,"*",p,":",l,"=>",nl,"\n");
      m:=nm;l:=nl;
    fi;
    nm:=Value(m,[var,Mvp(var)/p]);
    nl:=Length(Format(nm));
    if nl<l then 
      Print(var,"->",var,"*",p,":",l,"=>",nl,"\n");
      m:=nm;l:=nl;
    fi;
  od;
  return m;
end;

contentden:=function(m)
  m:=List(Flat(m),x->List(x.coeff,Denominator));
  m:=List(Flat(m),Factors);
  return Collected(Flat(m));
end;

contentnum:=function(m)local myfact;
  myfact:=function(n)local v,p;
    v:=[];
    for p in Primes do
      while n mod p=0 do Add(v,p);n:=n/p;od;
    od;
    return v;
  end;
  m:=List(Flat(m),x->List(x.coeff,Numerator));
  m:=List(Flat(m),myfact);
  return Collected(Flat(m));
end;
