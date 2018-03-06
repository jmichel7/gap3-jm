C3D6:=rec(group:=CoxeterGroup("D",6),levi:=[1,4,6],
  relgroup:=CoxeterGroup("C",3),
  unipclasses:=[[2,2,2,2,2,2],[3,2,2,2,2,1],[4,4,2,2],[4,4,3,1],
    [5,3,2,2],[6,6],[7,2,2,1],[7,5],[9,3],[11,1]],
  springer:=[[4],[1],[3],[2],[8],[6],[5],[7],[10],[9]]);
#C3D6.orderclasses:=OrderClassesClassical(C3D6.unipclasses);

# effet sur les caracteres de l'echange des racines dans F4
# (a ete calcule comme ci-dessous)
#eF4:=function()local W,c,ct;
#  W:=CoxeterGroup("F",4);
#  c:=List(ChevieClassInfo(W).classtext,x->PositionClass(W,EltWord(W,5-x)));
#  ct:=CharTable(W).irreducibles;
#  return PermListList(ct,List(ct,x->x{c}));
#end;
eF4:=(2,3)(5,7)(6,8)(11,12)(18,19)(21,23)(22,24);
# F4E7bis:=Copy(Springer(CoxeterGroup("E",7),2));
# F4E7bis.springer:=List(F4E7bis.springer,x->OnTuples(x,eF4));

#ordre partiel venant de l'algebre de Hecke associe a une table comme ci-dessus
HeckeAlg:=function(t)local G;
  G:=RelativeGroup(t.group,t.levi);
  return Hecke(G,List(G.parentMap,x->Mvp("q")^(2*CoxeterLength(t.group,x))));
end;

Sommers:=function(W)local l,res,u;
  l:=Concatenation(List(SemisimpleCentralizerRepresentatives(W),
  J->List(DistinguishedParabolicSubgroups(ReflectionSubgroup(W,J)),D->[J,D])));
  res:=List(l,function(p)local p,w,L;
    u:=rec(K:=p[1],dist:=p[2]);
    L:=ReflectionSubgroup(W,p[1]);
    w:=L.generatingReflections*0+2;
    w{L.rootRestriction{p[2]}}:=w{L.rootRestriction{p[2]}}*0;
    u.dynkin:=InduceRichardsonDynkin(W,L,w);
    u.name:=IsomorphismType(L,rec(TeX:=true));
    if p[2]<>[] then Append(u.name,SPrint("(a_",Length(p[2]),")"));fi;
    u.name:=Replace(u.name,"+","{+}");
    u.name:=String(u.name);
    return u;
    end);
  SortBy(res,x->x.dynkin);
  for u in CollectBy(res,x->[-Sum(x.dynkin),x.dynkin])
    do Print(IntListToString(u[1].dynkin)," ",
      Join(List(u,x->String(TeXStrip(x.name))),", "),"\n");od;
  return res;
end;

# centralizer of the unipotent part of a levi
orthlevi:=function(W,J)local L,l;
  L:=ReflectionSubgroup(W,J);
  l:=Filtered([1..W.N],i->OnTuples(J,Reflection(W,i))=J);
  l:=Filtered(l,i->ForAll(L.rootInclusion{[1..L.N]},
      k->not ((W.roots[i]+W.roots[k]) in W.roots)));
  l:=ReflectionSubgroup(W,l);
  return ApplyFunc(ReflectionGroup,l.type)*ReflectionSubgroup(CoxeterGroup("A",
  W.rank-L.semisimpleRank-l.semisimpleRank),[]);
end;

centralizerunip:=function(W)local l,uc,u,i,c,rab,isomin;
  l:=UnipotentClasses(W).classes;
  SortBy(l,x->-x.dimred-x.dimunip);
  for u in l do
    u.min:=orthlevi(W,List(u.balacarter,AbsInt));
    u.max:=ApplyFunc(ReflectionGroup,
      ReflectionSubgroup(W,Filtered(W.generatingReflections,
       i->u.dynkin[i]=0)).type);
    if u.dimred=Dimension(u.min) then u.red:=u.min;u.min:=true;fi;
    if u.dimred=2*u.max.N+u.max.semisimpleRank then
      u.red:=u.max;u.max:=true;
    fi;
    Print(String(u.name,10)," ",IntListToString(u.dynkin)," ");
    if IsBound(u.red) then 
      Print(IsomorphismType(u.red,rec(torus:=true)),"(");
      if u.min=true then Print("min");fi;
      if u.max=true then Print("max");fi;
      Print(")");
    else
      Print(IsomorphismType(u.min,rec(torus:=true)),"<",u.dimred,"<",
            IsomorphismType(u.max));
    fi;
    Print("\n");
  od;
end;

nrunip:=function(W)local g;g:=GenericOrder(W,x);
  return Sum(UnipotentClasses(W).classes,
           c->g/GenericOrder(c.red,x)/x^c.dimunip);
end;

# partition associated to nilpotent element with matrix m
NilpotentMatrixToPartition:=function(m)local p;
   p:=List([0..Length(m)],i->RankMat(m^i));
   return AssociatedPartition(List([1..Length(p)-1],i->p[i]-p[i+1]));
end;

UnipotentMatrixToPartition:=function(m)local p;
   m:=m-IdentityMat(Length(m));p:=List([0..Length(m)],i->RankMat(m^i));
   return AssociatedPartition(List([1..Length(p)-1],i->p[i]-p[i+1]));
end;

UalphaToMatrix:=function(u)local l,Ualpha;
  Ualpha:=function(r,c)local m,p;
    p:=Position(r,1);
    m:=IdentityMat(Length(r)+1);
    m[p][Number(r,x->x=1)+p]:=c;
    return m*c^0;
  end;
  l:=u.group.weylGroup;
  l:=l.roots{[1..l.N]};
  return Product(u.list,x->Ualpha(l[x[1]],x[2]));
end;

MatrixToUalpha:=function(m)local l,res,i,r,p,q,U,W;
  W:=CoxeterGroup("A",Length(m)-1);
  l:=W.roots;
  U:=UnipotentGroup(W);
  res:=U.Element(1,0);
  m:=Copy(m);
  for i in [1..Length(l)/2] do
    r:=l[i];
    p:=Position(r,1);
    q:=Number(r,x->x=1)+p;
    if m[p][q]<>0 then 
       res:=res*U.Element(i,m[p][q]);
       m[p]:=m[p]-m[q]*m[p][q];
    fi;
  od;
  return res;
end;
    
prin2:=function(j)local r,res,ptodiagram,richardson;
# richardson element of parabolic with partition p (with random elements)
richardson:=function(p)local n,s,a,dim;
  n:=Sum(p);
  s:=List([1..Length(p)],x->Sum(p{[1..x]}));
  s:=Concatenation(List([1..Length(p)],i->List([1..p[i]],y->n-s[i])));
  a:=List([1..Length(s)],i->List([1..s[i]],z->Random([1..7+i^z])));
  dim:=Length(a);
  return List([1..dim],i->Concatenation(List([1..dim-Length(a[i])],j->0),a[i]));
end;
# diagram of a parabolic p
  ptodiagram:=p->Join(List(p,x->List([1..x-1],y->'0')),"2");
  res:=List(OrderedPartitions(j),
    p->[NilpotentMatrixToPartition(richardson(p)),ptodiagram(p)]);
  res:=List(Set(List(res,x->x[1])),
    i->[i,List(Filtered(res,x->x[1]=i),y->y[2])]);
  for r in res do
    Print(String(IntListToString(r[1]),-15),Join(r[2]," "),"\n");
  od;
end;

# matrice nilpotente rentree comme la fin de ses lignes (arguments distincts)
rentre:=function(arg)local dim; dim:=Length(arg);
  return List([1..dim],i->Concatenation([1..dim-Length(arg[i])]*0,arg[i]));
end;

# support unipotent
UnipotentSupports:=function(W,check)local uc,v,s,uch,c,p,t,tr,sup;
 tr:=x->Replace(x,"lambda","l","nu","n");
 uc:=UnipotentClasses(W);
 uch:=UnipotentCharacters(W);
 sup:=List(uch.families,function(f)local l;
   l:=uc.springerSeries[1].locsys
     [Position(uch.harishChandra[1].charNumbers,f.charNumbers[f.special])];
   if l[2]<>PositionId(uc.classes[l[1]].Au) then Error();fi;
   return l[1];end);
 if check then
   v:=UnipotentValues(uc);
   s:=List([1..Length(v.scalar)],i->Filtered([1..Length(v.scalar[i])],
     j->v.scalar[i][j]<>0*v.scalar[i][j]));
   s:=List(s,x->Set(List(x,i->v.locsys[i][1])));
   s:=List(s,function(x)local m;
      m:=Minimum(List(x,i->uc.classes[i].dimBu));
      x:=Filtered(x,i->uc.classes[i].dimBu=m);
      if Length(x)>1 then Error();fi;
      return x[1];end);
   s:=TransposedMat([Set(s),CollectBy([1..Length(s)],s)]);
   if ForAny(s,function(c)local p;
     p:=PositionProperty(uch.families,f->Set(c[2])=Set(f.charNumbers));
     return c[1]<>sup[p];end) then Error();fi;
 fi;
 s:=List([1..Length(sup)],function(i)local p;
   c:=rec(class:=sup[i],family:=i);
   c.locsys:=List([1..NrConjugacyClasses(uc.classes[c.class].Au)],
     function(i)local p,n,f;
       p:=PositionProperty(uc.springerSeries,s->[c.class,i] in s.locsys);
       if p=1 then 
         n:=Position(uc.springerSeries[p].locsys,[c.class,i]);
	 n:=uch.harishChandra[1].charNumbers[n];
         p:=PositionProperty(uch.families,f->n in f.charNumbers);
	 f:=uch.families[p];
	 return SPrint(CharNames(uc.classes[c.class].Au)[i],"->",p,":",
	     tr(f.charLabels[Position(f.charNumbers,n)]));
       else return "*";
       fi;
       end);
   return c;end);
 t:=List(s,c->[uc.classes[c.class].Au,c.family,
     Length(uch.families[c.family].charNumbers),Join(c.locsys)]);
 Print(FormatTable(t,rec(rowLabels:=List(s,c->ClassName(uc.classes[c.class])),
   columnLabels:=["Au","Family","Fam#","locsys"],rowsLabel:="class")));
 # return s;
end;

DeligneLusztigAlmostCharacter:=function(W,i)local uc,f;
  uc:=UnipotentCharacters(W);
  f:=uc.operations.Fourier(uc);
  return UnipotentCharacter(W,f[i]);
end;

acTable:=function(W)local t,opt;
  t:=List([1..NrConjugacyClasses(W)],i->DeligneLusztigAlmostCharacter(W,i).v);
  opt:=rec();
  opt.rowLabels:=CharNames(W);
  opt.columnLabels:=CharNames(UnipotentCharacters(W));
  opt.screenColumns:=SizeScreen()[1];
  Print(FormatTable(t,opt));
end;

# display coeffs of p on the basis \prod_i=1^n (q-2i-1)/2(i+1)
# of polynomials taking integral values on odd numbers
decp2:=function(p)local n,x,basis,coeff,i,r,d,l;
  if Degree(p)<1 then return p;fi;
  n:=Variables(p)[1];x:=Mvp(n);
  basis:=[x^0];
  Append(basis,List([1..Degree(p)],i->Product([1..i],j->(x-2*j+1)/2/j)));
#  Print(basis,"\n");
  coeff:=[0..Degree(p)]*0;
  i:=Degree(p);
  while p<>0 do
    l:=Coefficients(p,n);d:=Length(l);
    r:=l[d]/Coefficients(basis[d],n)[d];
    coeff[d]:=coeff[d]+r;
    p:=p-r*basis[Length(l)];
  od;
  return ValuePol(coeff,Mvp("y"));
end;

gg:=function(arg)local W,ucl,t,q,g,l,res,uch,opt,s,title,i,b;
  q:=X(Rationals);W:=arg[1];
  ucl:=UnipotentClasses(W);
  uch:=UnipotentCharacters(W);
  t:=ICCTable(ucl);
  s:=DetPerm(W);
# Gamma_\iota=a_\iota\sum_\kappa P'_{\iota,\kappa}q^{b_\kappa-b_\iota}
# R_{\hat\kappa}
  res:=List([1..Length(t.scalar)],
    i->Sum([1..Length(t.scalar)],j->Value(t.scalar[i][j],q^-1)*
      q^(t.dimBu[i]-t.dimBu[j])*DeligneLusztigAlmostCharacter(W,s[j])));
  res:=List(res,x->List(x.v,y->decp2(Value(y,Mvp("x")))));
  title:="ggg ";
  opt:=rec(rows:=[1..Length(t.scalar)]);
  if Length(arg)=2 then opt.classes:=true;fi;
  if IsBound(opt.classes) then
    PrintToString(title," on unipotent classes\n");
    res:=TransposedMat(res);
    for i in [1..Length(t.uc.classes)] do
      b:=Filtered([1..Length(t.locsys)],j->t.locsys[j][1]=i);
      res{[1..Length(res)]}{b}:=res{[1..Length(res)]}{b}*
	 CharTable(t.uc.classes[i].Au).irreducibles;
    od;
    res:=TransposedMat(res);
    opt.rowLabels:=List(t.locsys,p->ClassName(t.uc.classes[p[1]],
	Inherit(rec(class:=p[2]),opt)));
  else opt.rowLabels:=List(t.locsys,p->ClassName(t.uc.classes[p[1]],
      Inherit(rec(locsys:=p[2]),opt)));
    PrintToString(title," on local systems\n");
  fi;
  SortParallel(t.dimBu,opt.rows);
  opt.columnLabels:=CharNames(uch);
  opt.screenColumns:=SizeScreen()[1];
  Print(FormatTable(res,opt));
  return [res,opt];
end;  

CurtisDuality:=function(uch)local l1,l2;
  l1:=Concatenation(List(uch.harishChandra,x->x.charNumbers));
  l2:=Concatenation(List(uch.harishChandra,x->x.charNumbers{
     DetPerm(ApplyFunc(ReflectionGroup,x.relativeType))}));
  SortParallel(l1,l2);
  return l2;
end;

classfam:=function(W,i)local f,ucl;
  f:=UnipotentCharacters(W).families[i];
  ucl:=UnipotentClasses(W);
  return ucl.springerSeries[1].locsys[f.charNumbers[f.special]][1];
end;

goodfams:=function(W)local uc,ucl,f;
  uc:=UnipotentCharacters(W);
  ucl:=UnipotentClasses(W);
  f:=[1..Length(uc.families)];
  f:=Filtered(f,i->Size(uc.families[i])>1);
  return Filtered(f,i->Size(uc.families[i])=
     NrDrinfeldDouble(ucl.classes[classfam(W,i)].Au));
end;

# ggg projected on family f
wavefront:=function(arg)local W,f,ucl,t,q,g,l,res,uch,opt,s,title,i,b,fid,
  fd,a,fch,locsys,class,p,ss,fchp,cd;
  q:=X(Rationals);W:=arg[1];f:=arg[2];
  ucl:=UnipotentClasses(W);
  uch:=UnipotentCharacters(W);
  f:=uch.families[f];
  fch:=List(f.charNumbers,x->Position(uch.harishChandra[1].charNumbers,x));
  fchp:=Filtered(fch,x->x<>false); # principal series
  fid:=fch[f.special];
  class:=ucl.springerSeries[1].locsys[fid][1];
  cd:=CurtisDuality(uch);
  fd:=Set(List(f.charNumbers,x->PositionProperty(uch.families,
    g->cd[x] in g.charNumbers)));
  if Length(fd)>1 then Error();fi;
  fd:=uch.families[fd[1]];
  p:=PositionProperty(ucl.springerSeries{[2..Length(ucl.springerSeries)]},
    s->class in List(s.locsys,x->x[1]));
  t:=ICCTable(ucl);
  locsys:=PositionsProperty(t.locsys,x->x[1]=class);
  s:=DetPerm(W);
# Gamma_\iota=a_\iota\sum_\kappa P'_{\iota,\kappa}q^{b_\kappa-b_\iota}
# R_{\hat\kappa}
  res:=List(locsys,
    i->Sum(fchp,j->Value(t.scalar[i][j],q^-1)*
      q^(t.dimBu[i]-t.dimBu[j])*DeligneLusztigAlmostCharacter(W,
        uch.harishChandra[1].charNumbers[s[j]])));
  res:=List(res,x->x.v{fd.charNumbers});
  title:="ggg ";
  opt:=rec();
  if Length(arg)=3 then opt.classes:=true;fi;
  if IsBound(opt.classes) then
    PrintToString(title," on unipotent classes\n");
    locsys:=t.locsys{locsys};
    if p<>false then 
      ss:=ucl.springerSeries[p+1];
      p:=PositionProperty(ss.locsys,x->x[1]=class);
      Add(locsys,ss.locsys[p]);
      Add(res,q^0*
        DeligneLusztigAlmostCharacter(W,ss.parameter[p]).v{fd.charNumbers});
    fi;
    SortParallel(List(locsys,x->x[2]),res);
    res:=TransposedMat(CharTable(ucl.classes[class].Au).irreducibles)*res;
    locsys:=List([1..NrConjugacyClasses(ucl.classes[class].Au)],i->[class,i]);
    opt.rowLabels:=List(locsys,p->ClassName(t.uc.classes[p[1]],
	Inherit(rec(class:=p[2]),opt)));
  else opt.rowLabels:=List(t.locsys{locsys},p->ClassName(ucl.classes[p[1]],
      Inherit(rec(locsys:=p[2]),opt)));
    PrintToString(title," on local systems\n");
  fi;
  res:=Concatenation([CharNames(uch){fd.charNumbers}],res);
  opt.rowLabels:=Concatenation(["uch"],opt.rowLabels);
  opt.columnsLabel:="symbol in fam";
  a:=List(f.charNumbers,x->[Position(uch.families,f),f.charLabels
     [Position(f.charNumbers,x)]]);
  b:=List(cd{f.charNumbers},function(x)local f,p;
     p:=PositionProperty(uch.families,f->x in f.charNumbers);
     return [p,
       uch.families[p].charLabels[Position(uch.families[p].charNumbers,x)]];
     end);
  a:=TransposedMat([a,b]);
  a:=Filtered(a,x->x[1]<>x[2]);
  if f=fd then a:=Filtered(a,x->x[1]<x[2]);fi;
  for b in a do Print(b[1][1],":",b[1][2],"<->",b[2][1],":",b[2][2],"\n");od;
# opt.rows:=[1..Length(locsys)];
  opt.columnLabels:=fd.charLabels;
  opt.screenColumns:=SizeScreen()[1];
  Print(FormatTable(res,opt));
  return [res,opt];
end;  
