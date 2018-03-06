# Probleme de Cedric.
# Dans Z(\BC W) on definit $c:=\sum_{w\in c} w$ pour c une classe de
# conjugaison, et $e_F=\sum_{\chi\in F}e_\chi$ pour $F$ une famille.
# On pose A_W=\sum_F \BC e_F$ et on filtre $i(c)=\codim(\ker(w-\id))$ pour
# $w\in c$. Pour $w$ élément régulier on a une application $F\mapsto F'$ des
# familles de W vers celles de $W':=C_W(W)$, d'où une application 
# $\pi:e_F\mapsto e_F'$ de $A_W\to A_W'$. Cette application est-elle filtrée?
# ie $\pi(A_W\cap\le i)\subset A_{W'}\cap\le i$?
#
# center of group algebra
# elements are represented internally as vectors in basis e_\chi
# elements are printed in basis of class sums.

CenterElementOps:=OperationsRecord("CenterElementOps");

CenterElementOps.\*:=function(x,y)local res;
  res:=ShallowCopy(x);
  res.v:=Zip(x.v,y.v,function(x,y)return x*y;end);
  return res;
end;

CenterElementOps.\+:=function(x,y)local res;
  res:=ShallowCopy(x); res.v:=x.v+y.v; return res;
end;

CenterElementOps.\^:=function(h,n)local p;
  if n=0 then return h.algebra.Element(1); fi;
  p:=false;
  while n>0 do
    if n mod 2<>0 then if p=false then p:=h;else p:=p*h;fi;fi;
    h:=h*h;n:=QuoInt(n,2);
  od;
  return p;
end;

CenterElementOps.Format:=function(x,opt)local v,n,letter;
  if IsBound(opt.classes) then 
    v:=SolutionMat(x.algebra.classes,x.v);letter:="C";
    n:=ChevieClassInfo(x.algebra.group).classnames;
  else v:=x.v;letter:="e";
    n:=ChevieCharInfo(x.algebra.group).charnames;
  fi;
  v:=TransposedMat([[1..Length(v)],v]);
  v:=Filtered(v,x->x[2]<>0);
  if v=[] then return "0";fi;
  v:=Concatenation(List(v,function(i)local s;
    s:=FormatCoefficient(i[2],SPrint(letter,"(",n[i[1]],")"),rec());
    if s[1]<>'-' then return SPrint("+",s);else return s;fi;end));
  if v[1]='+' then v:=v{[2..Length(v)]};fi;
  return v;
end;

PrintCenterElms:=rec(classes:=true);
CenterElementOps.Print:=function(x)local v,n;
  Print(Format(x,PrintCenterElms));
end;

CenterAlgebra:=function(W)local A,ct,uc,re;
  A:=rec(group:=W);
  ct:=CharTable(W);
  A.classes:=TransposedMat(List(ct.irreducibles,x->x/x[1]));
  A.classes:=Zip(A.classes,ct.classes,function(x,y)return x*y;end);
  A.operations:=rec(
    Print:=function(r)Print("CenterAlgebra(",r.group,")");end);
  A.Element:=function(i)local e;
    e:=rec(algebra:=A,operations:=CenterElementOps);
    if IsInt(i) then e.v:=A.classes[i];else e.v:=i;fi;
    return e;
  end;
  re:=ReflectionEigenvalues(W);
  A.filtration:=List([0..Length(re[1])],
    i->List(Filtered([1..Length(re)],j->Number(re[j],y->y<>0)=i),j->
      A.Element(j)));
  # intersections of subspace m with f<=i
  A.filter:=function(m)return List([1..Length(A.filtration)],
    i->SumIntersectionMat(m,Concatenation(List([1..i],
      j->List(A.filtration[j],y->y.v))))[2]);
  end;
  return A;
end;

# test regular number d
test:=function(W)local d,A,l,A1,uc,n,s,pi,fpi,ffamA,image,i;
  A:=CenterAlgebra(W);
  uc:=UnipotentCharacters(W);
  A.families:=List(uc.families,function(f)local v;
    v:=List(uc.harishChandra[1].charNumbers,function(i)
      if i in f.charNumbers then return 1;else return 0;fi;end);
    return A.Element(v);end);
  ffamA:=A.filter(List(A.families,x->x.v));
  for d in RegularEigenvalues(W) do
    s:=PrincipalSeries(W,d);
    A1:=CenterAlgebra(RelativeGroup(s));
    Hecke(s);n:=s.charNumbers;
    A1.families:=List(uc.families,function(f)local v;
      v:=List(n,function(i)
	if i in f.charNumbers then return 1;else return 0;fi;end);
      return A1.Element(v);end);
    A1.families:=Filtered(A1.families,x->x.v<>0*x.v);
    pi:=List([1..Length(uc.families)],function(i)local v;
      v:=Intersection(uc.families[i].charNumbers,n);
      if v=[] then return A1.families[1].v*0;
      else return First(A1.families,f->f.v[Position(n,v[1])]<>0).v;
      fi;end);
    fpi:=v->A1.Element(SolutionMat(List(A.families,x->x.v),v)*pi);
    image:=List(ffamA,x->List(x,fpi));
    for i in [1..Length(A1.filtration)] do
      if Length(A1.filter(List(image[i],x->x.v))[i])<Length(image[i])
      then Error();fi;
    od;
  od;
end;
