# classtext of "canonical" generator of the center
WordCenter:=W->ChevieClassInfo(W).classtext[PositionRegularClass(W,
  OrderCenter(W))];
  
PositionCox:=W->PositionRegularClass(W,Maximum(ReflectionDegrees(W)));

# values on w using it is in a parabolic
HeckeCharValuesParabolic:=function(H,w)local W,WI,p,q,pw;
  W:=Group(H);
  pw:=ParabolicClosure(W,Set(w));
  if Length(pw)=Length(w) then WI:=ReflectionSubgroup(W,Set(w));
  else WI:=ReflectionSubgroup(W,pw);
  fi;
  if WI.semisimpleRank=W.semisimpleRank then return false;fi;
  Print("w=",w);
  p:=List(WI.generators,x->Position(Reflections(W),x));
  w:=List(w,x->Position(p,x));
  if false in w then return false;fi;
  w:=WI.rootInclusion{w};
  q:=Product(H.parameter,x->Product(x))^0;
  H:=HeckeCharValues(Hecke(WI,H.parameter{W.orbitRepresentative{
    WI.rootInclusion{WI.generatingReflections}}}),w);
  if false in H then return false;fi;
  return InductionTable(WI,W).scalar*(H*q);
end;

# find the specialization of Mvp parameters leading to group algebra
SpecializationToGroup:=function(H)
  return Concatenation(Set(Concatenation(List(H.parameter,
   function(p)local e;e:=Length(p);
    return List([1..e],function(i)local v;
      v:=ScalMvp(p[i]);
      if v<>false then
        if v<>E(e)^(i-1) then Error(v," is not ",E(e)^(i-1));fi;
        return [];
      fi;
      if p[i].coeff<>[1] or Length(p[i].elm[1].coeff)<>1 then 
        Error(p[i]," is not a power of a variable");fi;
      return [p[i].elm[1].elm[1],GetRoot(E(e)^(i-1),p[i].elm[1].coeff[1])];
    end);end))));
end;

# get central characters without having to know Chartable
# for algebras above group algebra
HeckeCentralCharacters:=function(H)local W,v,v1,v2,z,spec;
  W:=Group(H);
  z:=Size(Centre(W));
  v:=List(HeckeCentralMonomials(H),m->GetRoot(m,z));
  v1:=CharTable(W).irreducibles{[1..Length(v)]}[PositionRegularClass(W,z)];
  v2:=CharTable(W).irreducibles{[1..Length(v)]}[PositionClass(W,W.identity)];
  return Zip(v,v1,v2,function(x,y,z)return x*y/z/
    Value(Mvp(x),SpecializationToGroup(H));end);
end;

# cutz(l,z) for list l returns [a,l'] such that replace a times z->[] in l gives l'
cutz:=function(l,z)local l1; l1:=Replace(l,z,[]);
  return [(Length(l)-Length(l1))/Length(z),l1];
end;

# encode class txt with ChevieClassInfo given
# ie replace z,c,etc...
EncodeClass:=function(txt,info)local dic;
  dic:=Filtered([1..Length(info.classtext)],i->Length(info.classnames[i])=1
    and Length(info.classtext[i])>1);
  SortBy(dic,x->-Length(info.classtext[x]));
  return String(ApplyFunc(Replace,Concatenation([IntListToString(txt)],
    Concatenation(List(dic,i->[IntListToString(info.classtext[i]),
    info.classnames[i]])))));
end;

# checks class representatives cunningly chosen w.r.t. parabolics
checkparabolic:=function(W)
  local WI,p,ct,zt,z,i,I,j,n,k,cl,a,b,res,f,pz,nums,cc,triples,s,ci,ps,l,ls,lI,zp;
  pz:=b->SPrint(IntListToString(b[2]),List([1..b[1]],x->'z'));
  ci:=ChevieClassInfo(W);
  ct:=ci.classtext;
  zt:=WordCenter(W);
  cc:=ct[PositionCox(W)];
  if cutz(zt,cc)[2]<>[] then Error("z=",zt," is not power of c=",cc,"\n");fi;
  Print("c=",IntListToString(cc)," z=",EncodeClass(zt,ci),"\n");
  cl:=List(ConjugacyClasses(W),Representative);
  z:=cl[PositionClass(W,EltWord(W,zt))];
  n:=W.nbGeneratingReflections;
  nums:=[];
  f:=function(i)local res;AddSet(nums,i);
    res:=SPrint("class ",String(i,3),":",EncodeClass(ct[i],ci));
    if EncodeClass(ct[i],ci)<>ci.classnames[i] then
      PrintToString(res,"(",ci.classnames[i],")");fi;
    return res;
  end;
  res:=[];triples:=[];
  l:=Combinations(W.generatingReflections);
  if W.semisimpleRank<W.nbGeneratingReflections then
    ls:=List(l,x->Length(ParabolicClosure(W,x)));
    lI:=l{Filtered([1..Length(l)],i->ls[i]=Length(l[i]))};
    ps:=w->ls[Position(l,Set(w))];
  else lI:=l;ps:=w->Length(Set(w));
  fi;
  lI:=Filtered(lI,x->Number(lI,j->IsSubset(j,x))=2);
  for I in lI do
    WI:=ReflectionSubgroup(W,I);
    if IsSubset(W.generatingReflections,WI.rootInclusion{WI.generatingReflections})
    then
    p:=FusionConjugacyClasses(WI,W);zp:=Set(p);
    for j in p do
      s:=ChevieClassInfo(WI).classtext[Position(p,j)];
      for k in [0..OrderPerm(z)-1] do
        a:=PositionClass(W,cl[j]*z^k);AddSet(zp,a);
        b:=cutz(ct[a],zt);
	if not IsSubset(I,Set(b[2]))and(Length(s)<Length(b[2])or ps(s)<ps(b[2]))
	then
	  Print("#",a," ",pz(b)," =>",pz([k,s])," ps=",ps(b[2]),"=>",ps(s),"\n");
          Add(res,SPrint(f(a)," could be ",pz([k,s])));
          Add(triples,[a,EncodeClass(ct[a],ci),pz([k,s])]);
        fi;
      od;
    od;
    Print("parabolic:",ReflectionName(WI)," gives classes ",zp,"\n");
    fi;
  od;
  for i in [1..Length(ct)] do
    b:=cutz(ct[i],zt);
    if ps(b[2])<n then
      if ps(b[2])=0 then Add(res,f(i));
      else Add(res,SPrint(f(i)," in ",pz([b[1],[]]),"W_",pz([0,Set(b[2])])));
      fi;
    elif b[1]>0 then Add(res,SPrint(f(i)));
    fi;
  od;
  for i in Set(res) do Print(i,"\n");od;
  Print("bad classes:",Difference([1..Length(ct)],nums),"\n");
  return triples;
end;

showregular:=function(W)local tbl,txt,d,cl,q,already,i,cc,lpi,lb,ci;
  ci:=ChevieClassInfo(W);
  txt:=ci.classtext;
  cc:=txt[PositionCox(W)];Print("c=",IntListToString(cc),"\n");
  cl:=List(ConjugacyClasses(W),Representative);
  lpi:=Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W));
  already:=[]; tbl:=[];lb:=[];
  for d in Reversed(RegularEigenvalues(W)) do
    q:=[];
    for i in PrimeResidues(d) do
      q:=PositionRegularClass(W,i/d);
      if not q in already then 
        Add(already,q);
        Add(lb,SPrint(i/d));
        Add(tbl,[q,lpi*i/d,Length(txt[q]),EncodeClass(txt[q],ci),
	  ChevieClassInfo(W).classnames[q]]);
      fi;
    od;
  od;
  Print(FormatTable(tbl,
    rec(rowLabels:=lb,
        columnLabels:=["class","th. lg","lg","txt","name"],
	rowsLabel:="RegEig")));
end;

# character values of Hecke algebra H on words non-cuspidal * z^i
fromparabolic:=function(H)local W,reps,n,i,t,v,a,zt,c;
  W:=Group(H); reps:=ChevieClassInfo(W).classtext; n:=Length(reps);
  t:=List([1..n],x->[1..n]*0+Unknown());
  zt:=WordCenter(W);
  v:=HeckeCentralCharacters(H);
  for i in [1..n] do
    a:=cutz(reps[i],zt);
    Print("class ",i,":z^",a[1],".",IntListToString(a[2]),"\n");
    if Length(Set(a[2]))<W.nbGeneratingReflections then
      c:=HeckeCharValuesParabolic(H,a[2]);
      if c<>false then
        Print(ReflectionName(ReflectionSubgroup(W,Set(a[2]))),"\n");
        t{[1..n]}[i]:=Zip(c,v,function(x,y)return x*y^a[1];end);
      fi;
    fi;
  od;
  return t;
end;

# character values of Hecke algebra H on words in W_I<z>
fromHI:=function(H,I)local W,reps,n,i,t,v,a,zt,c,nm,R,rgens,tbl,vv,HR;
  W:=Group(H); reps:=ChevieClassInfo(W).classtext; n:=Length(reps);
  t:=List([1..n],x->[1..n]*0+Unknown());
  zt:=WordCenter(W);
  v:=HeckeCentralCharacters(H);
  R:=ReflectionSubgroup(W,I);
  rgens:=R.rootInclusion{R.generatingReflections};
  nm:=ReflectionName(R);
  HR:=Hecke(R,H.parameter{W.orbitRepresentative{rgens}});
  tbl:=InductionTable(R,W).scalar;
  for i in [1..n] do
    a:=cutz(reps[i],zt);
    if IsSubset(I,Set(a[2])) and ForAll(a[2],i->i in rgens) then
      vv:=HeckeCharValues(HR,a[2]);
      if not false in vv then
        c:=tbl*vv;
        Print("class ",i,":",nm,"\n");
        t{[1..n]}[i]:=Zip(c,v,function(x,y)return x*y^a[1];end);
      fi;
    fi;
  od;
  return t;
end;

# merge chartable t with partial table vals
merge:=function(t,vals)local i,j;
  for i in [1..Length(t)] do
    if IsBound(vals[i]) then
    for j in [1..Length(t[i])] do
      if IsBound(vals[i][j]) and not IsUnknown(vals[i][j]) then 
        if not IsUnknown(t[i][j]) and t[i][j]<>vals[i][j] then Error(i,",",j);fi;
	t[i][j]:=vals[i][j];
      fi;
    od;
    fi;
  od;
end;

# give a primitive root of pi -- fills char for its powers
fromrootpi:=function(H,w)local W,p,o,v,v1,n,reps,t,x,bad,i;
  W:=Group(H);
  o:=(Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W)))/Length(w);
  if OrderPerm(EltWord(W,w))<>o then Error(w," is not root of pi");fi;
  p:=List([1..o-1],i->PositionClass(W,EltWord(W,w)^i));
  Print("gives classes ",p,"\n");
  p:=Set(List(p,i->Position(p,i)));
  p:=List(p,x->[x,PositionClass(W,EltWord(W,w)^x)]);
  reps:=ChevieClassInfo(W).classtext;
  bad:=Filtered([1..Length(p)],
         i->reps[p[i][2]]<>Concatenation(List([1..p[i][1]],y->w)));
  if bad<>[] then
    for i in bad do
    Print("Warning!! ",p[i][2],"-th representative is ",reps[p[i][2]],
                 " instead of ",w,"^",p[i][1],"\n");
    od;
  fi;
  n:=Length(reps);
  t:=List([1..n],x->[1..n]*0+Unknown());
  for x in Drop(p,bad) do
    v:=Zip(HeckeCentralMonomials(H),CharTable(W).irreducibles,
       function(m,c)return m^(x[1]/o)*c[x[2]]/Value(Mvp(m^(x[1]/o)),
       SpecializationToGroup(H));end);
    t{[1..n]}[x[2]]:=v;
  od;
  return t;
end;

fromreps:=function(arg)local H,inds,t,W,r,cl,i,n;
  H:=arg[1];
  W:=Group(H); cl:=ChevieClassInfo(W).classtext; n:=Length(cl);
  if Length(arg)=1 then inds:=[1..n];else inds:=arg[2];fi;
  t:=List([1..n],x->[1..n]*0+Unknown());
  for i in inds do
    r:=Representations(H,i);
    if r<>false then 
      Print(i," \c");
      t[i]:=CharRepresentationWords(r,cl);fi;
  od;
  Print("\n");
  return t;
end;

myeig:=function(table,i,class)local v;
  v:=Eigenvalues(table,table.irreducibles[i],class);
  return Concatenation(List([1..Length(v)],i->List([1..v[i]],j->i/Length(v))));
end;

coxclasses:=W->Set(List(Arrangements(W.generators,Length(W.generators)),
   x->PositionClass(W,Product(x))));

refssmaller:=function(W,w)local refs, r;
  refs:=Reflections(W);
  r:=Set(List(refs,x->Position(refs,x)));
  return Filtered(r,x->ReflectionLength(W,refs[x]^-1*w)<ReflectionLength(W,w));
end;

refl2:=function(W,w)local m; m:=MatXPerm(W,w); return RankMat(m-m^0);end;

test:=function(W)local a,b;
  Stime();
  a:=Set(List(HurwitzOrbitItems(W.generators),x->Position(Reflections(W),x)));
  Print("Hurwitz:",Stime(),"\n");
  b:=refssmaller( W, Product( W.generators ) );
  Print("smaller:",Stime(),"\n");
  return a=b;
end;

#Check eigenvalues of Coxeter elements
CheckCox:=function(W)local h,c,cl,i,ct,refs,n,gens,tbl,ref;
  PrintDiagram(W);
  h:=Maximum(ReflectionDegrees(W));
  c:=PositionRegularClass(W,h);
  n:=Length(W.generators);
  cl:=coxclasses(W);
  if not c in cl then
    Print("******* Warning *****: ",h,"-regular class not product of gens\n");
  fi;
  ct:=CharTable(W);
  refs:=[1..Length(ct.irreducibles)];
  gens:=List(W.generators,i->[PositionClass(W,i),W.rank-1+E(OrderPerm(i))]);
  refs:=Filtered(refs,i->ct.irreducibles[i][1]=W.rank);
  refs:=Filtered(refs,i->ForAll(gens,v->ct.irreducibles[i][v[1]]=v[2]));
  ref:=ChevieCharInfo(W).extRefl[2];
  Print("Chevie's reflection representation is number:",ref,"\n");
  Print("Reflection reps where generators are distinguished:",refs,"\n");
  Print("classes for products of generators:\n");
  for i in cl do 
    Print("class ",i,"\n",FormatTable(List(refs,j->myeig(ct,j,i)),
      rec(rowLabels:=refs,
      columnLabels:=
       Concatenation(["rep","eigenvalues"],List([1..W.rank-1],x->"")))),"\n");
  od;
  Print("\nh=",h,"-regular class =",c,"=",ChevieClassInfo(W).classtext[c],"\n");
end;

# find for each char. x of 1-param Hecke H an e such that x is in C[q^{1/e}]
gete:=function(H)local x,ct,fp;
  fp:=function(p)
    if IsCyc(p) then return [1];fi;
    p:=Filtered(p.elm,x->Length(x.elm)=1);
    p:=List(p,x->Mod1(x.coeff[1]));
    return Set(List(Union(Set(p),[0]),Denominator));
  end;
  ct:=CharTable(H);
  if ct<>false then 
    ct:=ct.irreducibles; return List(ct,x->Lcm(Union(List(x,fp))));
  else
    return List(HeckeCentralCharacters(H),x->Lcm(fp(x)));
  fi;
end;

# check "z/e": for each family if characters of Spetsial algebra in family
# split over q^{1/e} then a+A is divisible by |ZW|/e
checkze:=function(W)local uc,e,z,f,aA,fe;
  uc:=UnipotentCharacters(W);
  e:=ChevieCharInfo(W).opdam;
  z:=OrderCenter(W);
  Print("   z=",z,"\n");
  for f in uc.families do
    aA:=uc.a[f.charNumbers[1]]+uc.A[f.charNumbers[1]];
    fe:=List(f.charNumbers,i->Position(uc.harishChandra[1].charNumbers,i));
    fe:=Filtered(fe,x->x<>false);
    fe:=Lcm(List(Orbits(Group(e),fe)));
    Print([aA,fe,z/Gcd(aA,z),aA/(z/fe)],"\n");
  od;
end;
