ReadVK("contrib/det");
# This files does the computations common to finding the discriminant
# for G24, G27, G29, G31, G33 and G34.
# All results are stored as various fields of a record res;

x:=Mvp("x"); y:=Mvp("y"); z:=Mvp("z"); t:=Mvp("t"); u:=Mvp("u"); v:=Mvp("v");
varnames:=[];vars:=[]; invars:=[]; #to shut down warnings

init:=function(W)local exp,j,aux,d,cd,res;
  res:=rec();
  res.degs:=ReflectionDegrees(W);d:=res.degs;
  res.rank:=Length(d);
  cd:=ReflectionCoDegrees(W);
  # finding the ordering of the codegrees making C triangular
  cd:=Filtered(Arrangements(cd,res.rank),function(a)local C;
    C:=List(d,x->List(a,y->d[1]+y-x));
    return Number([1..Length(C)],i->C[i][i]=0)>=Length(C)-1;end);
  if Length(cd)>1 then 
    Error("there should be only one possible order of codegres");
  fi;
  cd:=cd[1];
  res.codegs:=cd;
  Print("degs=",d," codegs=",cd,"\n");
  Print("degrees of C:\n",Format(List(d,x->List(cd,y->d[1]+y-x))),"\n");
  res.q:=PositionProperty(d[1]+cd-d,x->x>0);
  varnames:="xyztuv";varnames:=List(varnames{[1..res.rank]},x->[x]);
  vars:=List(varnames,Mvp);
  if res.q<>false then
    res.hesspos:=Position(d,d[1]+cd[res.q]-d[res.q]);
    Print("hessian is ",Ordinal(res.hesspos)," degree\n");
  fi;
  return res;
end;

# Computes p mod m for univariate Mvps; gives result as a vector of coeffs
Mvpmod:=function(p,m)local v,var,c,rest,i;
  v:=Variables(m);
  if Length(v)<>1 then Error(m," is not a univariate polynomial");fi;
  var:=v[1];
  c:=Coefficients(m,var);
  rest:=-Sum([2..Length(c)],i->Mvp(var)^(i-2)*c[i-1])/c[Length(c)];
  p:=Coefficients(p,var);
  while Length(p)>=Length(c) do
#   Print(" p=",p,"\n");
    v:=List([1..Length(c)-1],i->Mvp(var)^(i-1));
    Append(v,List([Length(c)..Length(p)],i->Mvp(var)^(i-Length(c))*rest));
#   Print("v=",v);
    p:=v*p;
#   Print(" p=",p,"\n");
    p:=Coefficients(p,var);
  od;
  return p;
end;
    
getC:=function(res)local A,z,h,C1,p,v,i,res,aux,j,exps,wt,m;
  Print("q=",res.q,": end of ",Ordinal(res.q)," line of C:\n");
  aux:=["a","b","c","d","e","f","g","h","i","j","k","l"];
  j:=0;
  res.C:=IdentityMat(res.rank)*x^0;
  res.C[res.q][res.q]:=vars[res.hesspos];
  wt:=List([res.q+1..res.rank],i->res.degs[1]+res.codegs[i]-res.degs[res.q]);
  exps:=List(wt,w->Filtered(Cartesian(List(res.degs,i->[0..QuoInt(w,i)])),
     v->v*res.degs=w and v[res.hesspos]=0));
  res.C[res.q]{[res.q+1..res.rank]}:=List(exps,x->Sum(List(x,y->Product([1..Length(y)],i->
     vars[i]^y[i])),function(m)j:=j+1;return m*Mvp(aux[j]);end));
  Print(Format([wt,res.C[res.q]{[res.q+1..res.rank]}]),"\n");
  z:=[1..res.rank-1]*Mvp("x")^0;Add(z,vars[res.rank]);
  C1:=Value(res.C,
       Concatenation(List([1..res.rank],i->[varnames[i],res.invar(i,z)])));
  h:=res.hessmat(z);
  A:=myCoF(h)*TransposedMat(res.jacobmat(z))*C1;
  h:=myDet(h);
  A:=Set(Concatenation(A));
  A:=Concatenation(List(A,x->Mvpmod(x,h)));
  A:=Filtered(A,x->x<>0*x)*Mvp("x")^0;
  A:=Set(List(A,x->x/x.coeff[Length(x.coeff)]));
  h:=Set(Concatenation(List(A,x->x.elm)));
  h:=List(Filtered(h,x->Length(x.elm)>0),x->x.elm[1]);
  Sort(h);
  Add(h,"");
  m:=[];
  for p in A do
    v:=[1..Length(h)]*0;
    for i in [1..Length(p.coeff)] do
      if Length(p.elm[i].elm)=0 then
        v[Position(h,"")]:=p.coeff[i];
      else v[Position(h,p.elm[i].elm[1])]:=p.coeff[i];
      fi;
    od;
    Add(m,v);
  od;
  m:=TransposedMat(m);
  v:=SolutionMat(m{[1..Length(m)-1]},-m[Length(m)]);
  v:=List([1..Length(h)-1],i->[h[i],v[i]]);
  res.C:=Value(res.C,Concatenation(v));
  res.c:=z->ScalMvp(Value(res.C,Concatenation(List([1..res.rank],
   i->[varnames[i],res.invar(i,z)]))));
end;

# fast version of ScalMvp(Value(p,Concatenation(TransposedMat([varnames,z]))))
applyMvp:=function(p,z,varnames)local res,i,j,e,cres;
  res:=0;
  for i in [1..Length(p.coeff)] do
    cres:=p.coeff[i];
    e:=p.elm[i];
    for j in [1..Length(e.elm)] do
       cres:=cres*z[Position(varnames,e.elm[j])]^e.coeff[j];
    od;
    res:=res+cres;
  od;
  return res;
end;

# 'bord(u,v)' as defined in Orlik-Terao p. 284
bord:=function(arg) local u,v,varnames,D1,D2,D3,B;
  u:=arg[1];v:=arg[2];
  if Length(arg)=2 then varnames:=["x","y","z"];else varnames:=arg[3];fi;
  D1:=p->Derivative(p,varnames[1]);
  D2:=p->Derivative(p,varnames[2]);
  D3:=p->Derivative(p,varnames[3]);
  B:= [[D1(D1(u)),D1(D2(u)),D1(D3(u)),D1(v)],
       [D2(D1(u)),D2(D2(u)),D2(D3(u)),D2(v)],
       [D3(D1(u)),D3(D2(u)),D3(D3(u)),D3(v)],
       [D1(v),    D2(v),    D3(v),    0*Mvp("x")]];
  return myDet(B);
end;

# the following global variables must be defined:
# vars,invars,varnames

# the following functions in record res must be defined:
# hessmat, jacobmat, c

basicder:=function(res) local points,cs,Ms,addpoint,
    i,j,k,e,M,wt,exps,invars,sample,matrix,vector,bit,coefs,
    lastgood,newpoint,getpoints,security,mm;
  points:=[];cs:=[];Ms:=[];invars:=[];
  security:=3;

  newpoint:=function()return List([1..res.rank],i->Random([-6..6]));end;

  addpoint:=function() local z,j,cc,h;
    z:=newpoint();
    Print("Testing new point...\c");
    h:=res.hessmat(z); 
    if Rank(h) = res.rank then
      Print("Adding new point...\c");
      Add(points,z); 
      Add(invars,List([1..res.rank],i->res.invar(i,z)));
      j:=res.jacobmat(z);
      cc:=res.c(z);  Add(cs,cc);
      Add(Ms,(res.degs[1]-1)*j*h^-1*TransposedMat(j)*cc);
    else addpoint();
    fi;
    Print("done (",Length(points)," points).",Stime(),"\n");
  end;

  getpoints:=function(n)
    while Length(points)<n do addpoint();od;
  end;

  M:=List([1..res.rank],i->[1..res.rank]);
  for i in [1..res.rank] do
    for j in [1..res.rank] do
      Print("Starting coef ",i,",",j,"\n");
      wt:=res.degs[i]+res.codegs[j];
      exps:=Filtered(Cartesian(List(res.degs,i->[0..QuoInt(wt,i)])),v->v*res.degs=wt);
      getpoints(Length(exps)+security);
      repeat 
	sample:=[];
	repeat AddSet(sample,Random([1..Length(exps)]));
	until Length(sample)=Length(exps);
	matrix:=List([1..Length(exps)],k-> List(exps,e->
	     Product(List([1..res.rank],j->invars[sample[k]][j]^e[j]))));
	vector:=Ms{[1..Length(exps)]}[i][j];
	bit:=true;
	Print("starting rank...\c");
	if RankMat(matrix)< Length(exps) then
	   Print("Bad luck...\n");
	   points:=[] ; Ms:=[]; invars:=[] ; 
	   getpoints(Length(exps)+security);
	   bit:=false;
	fi;
	Print("rank done ",Stime(),"\n");
	Print("starting solutionmat...\c");
	coefs:=SolutionMat(TransposedMat(matrix),vector);
	Print("SolutionMat done ",Stime(),"\n");
      until bit=true;
      M[i][j]:=Sum(List([1..Length(exps)],
       i->coefs[i]*Product(List([1..res.rank],j->vars[j]^exps[i][j]))));
      Print("checking extra points...\c");
      for k in [1..Length(points)] do
       if ScalMvp(Value(M[i][j],
       Concatenation(List([1..res.rank],i->[varnames[i],invars[k][i]]))))
	<> Ms[k][i][j] then Error("Not invariant"); fi; 
      od;
      Print("done ",Stime(),"\n");
      Print(M[i][j],"\n");
    od;
  od;
  return M;
end;

common24_33:=function(res)
  res.invar:=function(i,z)return applyMvp(invars[i],z,varnames);end;
  res.J:=Jacobian(invars,varnames);
  res.hessmat:=z->List(res.H,l->List(l,p->applyMvp(p,z,varnames)));
  res.jacobmat:=z->List(res.J,l->List(l,p->applyMvp(p,z,varnames)));
  getC(res);
  res.M:=basicder(res);
  res.disc:=myDet(res.M);
end;
