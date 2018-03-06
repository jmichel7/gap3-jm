#Programs written in Abu Dhabi march 2017
# Wess-Zumino-Witten algebra for I_2(e) (according to Lusztig)
WZW:=function(e)local f,A,l,B,C,ind;
  f:=Family(rec(fourierMat:=List([1..e-1],i->List([1..e-1],
    j->(E(2*e)^(i*j)-E(2*e)^(-i*j))))/ER(-2*e),
   eigenvalues:=List([1..e-1],k->-E(24)^-1*E(4*e)^(-k^2))));
  A:=FusionAlgebra(f*f);
  l:=Cartesian([1..e-1],[1..e-1]);
  ind:=PositionsProperty(List(l,Sum),i->0=(i mod 2));
  l:=l{ind};
  A:=SubAlgebra(A,A.basis{ind});
  B:=QuotientAlgebra(A,[A.basis[1]-A.basis[Length(A.basis)]]);
  C:=QuotientAlgebra(A,[A.basis[1]+A.basis[Length(A.basis)]]);
  ind:=PositionsProperty(l,x->x[1]<x[2] or (x[1]=x[2] and 2*x[1]<e));
  return Basis(C,"U",rec(value:=List(A.basis{ind},C.projection),
                         parameters:=List(l{ind},IntListToString)));
end;

# factorize a quadratic form
qfactor:=function(p)local v,r,m,i,t,e,n,b,d;
  if Degree(p)<>2 then Error("not quadratic form");fi;
  v:=Variables(p);r:=Length(v)+1;
  m:=NullMat(r);
  for i in [1..Length(p.elm)] do
    t:=p.coeff[i];e:=p.elm[i];n:=List(e.elm,x->Position(v,x));e:=e.coeff;
    if e=[1,1] then m[n[1]][n[2]]:=t/2;m[n[2]][n[1]]:=t/2;
    elif e=[2] then m[n[1]][n[1]]:=t;
    elif e=[1] then m[n[1]][r]:=t/2; m[r][n[1]]:=t/2;
    elif e=[] then m[r][r]:=t;
    else Error();
    fi;
  od;
  n:=Copy(m);
  TriangulizeMat(m);m:=Filtered(m,x->x<>0*x);
  if Length(m)>2 then return false;fi;
  t:=List(n,x->SolutionMat(m,x));
  m:=List(m,x->SolutionMat(TransposedMat(t),x));
  v:=Concatenation(List(v,Mvp),[Mvp("x")^0])*t;
  if Length(m)=1 then return [v[1],v[1]*m[1][1]];fi;
  b:=m[1][2]+m[2][1];
  if m[1][1]=0 then return [v[2],b*v[1]+m[2][2]*v[2]];fi;
  b:=b/m[1][1];
  d:=GetRoot(b^2-4*m[2][2]/m[1][1],2,"no");
  if d=false then return false;fi;
  return [v[1]+v[2]/2*(b-d),m[1][1]*(v[1]+v[2]/2*(b+d))];
end;

# Solve system of quadratic forms
MvpSolve:=function(l)local inner,b,varnames,n,funcs,sortbylg,cleanlist;
  sortbylg:=function(l)
    l:=Filtered(Set(l),x->x<>0*x);
    SortBy(l,x->Length(x.elm));
    return l;
  end;
  cleanlist:=l->l{Filtered([1..Length(l)],i->not ForAny([1..i-1],j->
                      MvpOps.ExactDiv(l[i],l[j])<>false))};
  funcs:=[
    function(l)local p,e; # linear
      p:=PositionProperty(l,e->Degree(e)=1);
      if p=false then return false;else e:=l[p];fi;
      Print(e,"=c>");
      return [solvefor(e,Variables(e)[1])];
    end,
    function(l)local p,e,i; # quadratic
      for i in [1..Length(l)] do
        if Degree(l[i])=2 then
          e:=l[i];p:=qfactor(e);
          if p<>false then 
            Print(e,"=b>");
            return [solvefor(p[1],Variables(p[1])[1]),
                    solvefor(p[2],Variables(p[2])[1])];
          fi;
        fi;
      od;
      return false;
    end,
    function(l)local p,e,var,v1,v2,ll,i,v,b,vv,I,vv1;
      l:=Filtered(l,x->Length(Variables(x))=2);
      ll:=CollectBy(l,x->Variables(x));
      for l in ll do
        v1:=Variables(l);v2:=v1[2];v1:=v1[1];
        b:=[rec(coeff:=[],elm:=[]),
            rec(coeff:=[1],elm:=[v1]),
            rec(coeff:=[1],elm:=[v2]),
            rec(coeff:=[1,1],elm:=[v1,v2]),
            rec(coeff:=[2],elm:=[v1]),
            rec(coeff:=[2],elm:=[v2])];
        I:=IdentityMat(6);
        v:=List(l,function(p)local r; r:=[1..6]*0;
          for i in [1..Length(p.coeff)] do 
            r[Position(b,p.elm[i])]:=p.coeff[i];od;
          return r;end);
        vv:=SumIntersectionMat(v,[I[5],I[6]])[2];
        if vv<>[] then
          vv:=vv[1]{[5,6]};
          e:=GetRoot(-vv[2]/vv[1],2,"no");
          if e=false then return false;fi;
          e:=e*Mvp(v2);
          Print(l,"=e>");
          return [[v1,e],[v1,-e]];
        fi;
        vv:=SumIntersectionMat(v,[I[1],I[2],I[5]])[2];
        if vv<>[] then
          e:=quadratMvp(ValuePol(vv[1]{[1,2,5]},Mvp(v1)));
          if e=false then return false;fi;
          Print(l,"=g>");
          return e;
        fi;
        vv:=SumIntersectionMat(v,[I[1],I[3],I[6]])[2];
        if vv<>[] then
          e:=quadratMvp(ValuePol(vv[1]{[1,3,6]},Mvp(v1)));
          if e=false then return false;fi;
          Print(l,"=h>");
          return e;
        fi;
        vv:=SumIntersectionMat(v,[I[3],I[5]])[2];
        vv1:=SumIntersectionMat(v,[I[2],I[6]])[2];
        if vv<>[] and vv1<>[] then
          vv:=vv[1]{[3,5]};vv1:=vv1[1]{[2,6]};
          e:=GetRoot(-vv[1]*vv1[1]^2/vv[2]/vv1[2]^2,3,"no");
          if e=false then return false;fi;
          Print(l,"=f>");
          return [[v2,0],[v2,e],[v2,e*E(3)],[v2,e*E(3)^2]];
        fi;
#       Error(cleanlist(l));
      od;
      return false;
    end];
  inner:=function(l)local vv,f,cnt;
    l:=sortbylg(l);
    if Length(l)=0 then return [[]];
    elif ForAny(l,x->Degree(x)=0) then return [];
    else 
      for cnt in [1,2] do
      for f in funcs do
      vv:=f(l);
      if vv<>false then 
        Print(vv,"\n");
        return Concatenation(List(vv,v->List(inner(Value(l,v)),
            x->Concatenation(v,x))));
      fi;
      od;
      l:=shr(l);
      od;
      Error(l,"\n"); 
      return [l];
    fi;
  end;
  varnames:=Variables(l);b:=inner(l);
  b:=List(b,function(v)local vars,l,i,r;
    vars:=v{[1,3..Length(v)-(Length(v)mod 2)-1]};
    vars:=Filtered(vars,IsString);
    l:=Length(vars);
    for i in [1..l] do
      if IsMvp(v[2*i]) then
      while Intersection(Variables(v[2*i]),vars)<>[] do
        v[2*i]:=Value(v[2*i],v{[1..2*l]});od;
      fi;
    od;
    r:=List([1..l],i->v{[2*i-1,2*i]});
    SortBy(r,x->Position(varnames,x[1]));
    r:=Concatenation(r);
    if Length(v)>2*l then Add(r,v{[2*l+1..Length(v)]});fi;
    return r;
    end);
  return b;
end;

# equations solved in LinearCharacters
leq:=function(A)local n,b,vars,l,varnames;
  b:=A.basis;
  n:=[1..A.dimension];l:=1+LogInt(A.dimension,10); 
  varnames:=List(n,i->SPrint("x",Replace(String(i,l)," ","0"))); 
  vars:=List(varnames,Mvp); 
  b:=Set(Flat(List(n,i->List(n,function(j)local r;
    r:=b[i]*b[j];return vars[i]*vars[j]-Sum(r.coefficients,
      p->p[1]*vars[p[2]]);end))));
  b:=sortbylg(b);
  return b;
end;

LinearCharacters:=function(A)local n,vars,varnames,b,l;
  b:=A.basis;
  n:=[1..A.dimension];l:=1+LogInt(A.dimension,10); 
  varnames:=List(n,i->SPrint("x",Replace(String(i,l)," ","0"))); 
  vars:=List(varnames,Mvp); 
  b:=Set(Flat(List(n,i->List(n,function(j)local r;
    r:=b[i]*b[j];return vars[i]*vars[j]-Sum(r.coefficients,
      p->p[1]*vars[p[2]]);end))));
  b:=MvpSolve(b);
  if ForAny(b,n->n{[1,3..Length(n)-1+(Length(n) mod 2)]}<>varnames) 
  then Error(b);fi;
  b:=List(b,n->n{[2,4..Length(n)]});
  if ScalMvp(b)=false then Error("not scal");fi;
  b:=ScalMvp(b);
  b:=Filtered(b,x->ForAny(x,y->y<>0));
  return b;
end;

# S-matrix of fusion algebra A
smat:=function(A)local m,d;
  m:=CharTable(A).irreducibles;
  d:=DiagonalOfMat(m*SignedPermutationMat(A.involution)*TransposedMat(m));
  d:=List(d,GetRoot);
  return TransposedMat(Zip(m,d,function(a,b)return a/b;end));
end;

chekf:=function(f)local A;
  A:=FusionAlgebra(f);
  return PermListList(CharTable(A).irreducibles,LinearCharacters(A));
end;

chekm:=function(f)local A,m,F;
  A:=FusionAlgebra(f);F:=Fourier(f);
  m:=smat(A);
  return SignedPermListList(TransposedMat(m),TransposedMat(F));
end;

# check that l is a linear character of A
IsChar:=function(A,l)local b,r,n;
  b:=A.basis;n:=[1..A.dimension];
  return ForAll(n,i->ForAll(n,function(j)
    r:=b[i]*b[j];r:=r.coefficients;
    return Sum(r,x->x[1]*l[x[2]])=l[i]*l[j];end));
end;

# multiplication table of A
mt:=function(A)Print(Format(List(A.basis,x->x*A.basis)),"\n");end;
pmat:=function(m)Print(Format(m),"\n");end;

# idempotents of fusion algebra from chartable
idemfromct:=function(A)local d,e,ct;
  ct:=CharTable(A).irreducibles;
  d:=DiagonalOfMat(ct*ComplexConjugate(TransposedMat(ct)));
  e:=DiagonalMat(d)^-1*ct*A.involution;
  return e;
end;

# find families of W with fusion algebra structure constants in 0,1,-1
coeff1:=function(W)local ff;
  ff:=UnipotentCharacters(W).families;
  return Filtered(ff,function(f)local A,b;
    if Size(f)=1 then return false;fi;
    A:=FusionAlgebra(f);
    b:=Basis(A);
    b:=Set(Flat(List(b,x->List(b,y->x*y))));
    b:=Set(Flat(List(b,x->List(x.coefficients,y->y[1]))));
    return IsSubset([-1,0,1],b);
  end);
end;

# show explanations of families
show:=function(W)local f,ff,t;
  ff:=UnipotentCharacters(W).families;
  Print("families for ",ReflectionName(W),"\n");
  t:=List(ff,x->[Size(x),TeXStrip(x.explanation)]);
  Print(FormatTable(t,rec(rowLabels:=[1..Length(ff)],
   columnLabels:=["size","explanation"])));
end;

# quantum dimensions of elements of a family
qdim:=function(f)local n,s,cc;
  n:=Size(f);
  s:=List([1..n],i->PositionProperty(f.mellin[i],j->j<>0));
  Print("s=",s);
  cc:=List([1..n],i->Sum(f.mellin[s[i]],x->x^2));
  cc:=List(cc,x->cc[1]/x);
  Print("classes=",cc);
  return List([1..n],i->f.eigenvalues[i]*f.mellin[s[i]][i]*cc[i]);
end;

# algebra with absolute value for structure constants
AbsAlgebra:=function(A0)local A;
  A:=rec(field:=Rationals,
   operations:=OperationsRecord("AbsAlgebraOps",FDAlgebraOps),
   type:="Abs Fusion algebra",basisname:="T");
  A.dimension:=Length(A0.structureconstants);
  A.parameters:=[1..A.dimension];
  A.zero:=AlgebraElement(A,[]);
  A.one:=AlgebraElement(A,[[Rationals.one,A0.one.coefficients[1][2]]]);
  A.operations.underlyingspace(A);
  A.structureconstants:=List([1..A.dimension],i->List([1..i],j->
    List(A0.structureconstants[i][j],p->[AbsInt(p[1]),p[2]])));
  A.identification:=[A.type,A.field,A.structureconstants];
  A.multiplication:=function(i,j)
    if i>=j then return A.structureconstants[i][j];
            else return A.structureconstants[j][i];fi;end;
  A.operations.Print:=function(A)Print(A.type," dim.",A.dimension);end;
  return A;
end;

sh:=function(f)local F,o;
  F:=Fourier(f);
  o:=DiagonalMat(Eigenvalues(f));
  if IsBound(f.lusztig) and f.lusztig then return o^-1*F*o;
  else return o^-1*F*o^-1;
  fi;
end;
