#############################################################################
##
#A  ucl.g      Unipotent Classes    Jean Michel
##
#Y  Copyright 2008-2016, Univ. Paris VII
##  
##  This file contains basic functions for unipotent classes of
##  reductive groups.
##   

#------------------------------------------------------------------------------
#  Unipotent class records contain the following fields:

UnipotentClassOps:=OperationsRecord("UnipotentClassOps");

# compute name of unipotent class according to options opt
UnipotentClassOps.Name:=function(arg)local cl,opt,n;
  cl:=arg[1];if Length(arg)=2 then opt:=arg[2];else opt:=rec();fi;
  if IsBound(opt.mizuno) and IsBound(cl.mizuno) then n:=cl.mizuno;
  elif IsBound(opt.shoji) and IsBound(cl.shoji) then n:=cl.shoji;
  else n:=cl.name;fi;
  if not IsBound(opt.TeX) then 
    n:=String(Replace(n,"_","","\\tilde ","~","{","","}",""));fi;
  if IsBound(opt.locsys) then
    if opt.locsys=PositionId(cl.Au) then return n;fi;
    cl:=SPrint("(",CharNames(cl.Au,opt)[opt.locsys],")");
    if IsBound(opt.TeX) then return SPrint(n,"^{",cl,"}");
    else return SPrint(n,cl);fi;
  elif IsBound(opt.class) then
    if opt.class=PositionId(cl.Au) then return n;fi;
    cl:=ChevieClassInfo(cl.Au).classnames[opt.class];
    if IsBound(opt.TeX) then return SPrint("\\hbox{$",n,"$}_{(",cl,")}");
    else return SPrint(n,"_",cl);fi;
  else return n;
  fi;
end;

UnipotentClassOps.String:=cl->SPrint("UnipotentClass(",Name(cl),")");

UnipotentClassOps.Print:=function(cl)Print(String(cl));end;

UnipotentClassOps.FormatCentralizer:=function(u,opt)local c,AuName,n;
  c:="";
  AuName:=function(u)local res,au;
    if Size(u.Au)=1 then return "";fi;
    if IsBound(u.AuAction) or (IsBound(u.dimred) and u.dimred=0) then res:=".";
    else res:="?";fi;
    au:=ReflectionName(ReflectionType(u.Au),opt);
    if IsBound(opt.TeX) then au:=
      Replace(au,"A_1","Z_2","A_2","S_3","A_3","S_4","A_4","S_5","B_2","D_8");
    else au:=Replace(au,"A1","Z2","A2","S3","A3","S4","A4","S5","B2","D8");
    fi;
    Append(res,au);
    return res;
  end;
  if IsBound(u.dimunip) then
    if u.dimunip>0 then PrintToString(c,Format(Mvp("q")^u.dimunip,opt));fi;
  else PrintToString(c,"q^?");fi;
  if IsBound(u.AuAction) then 
    if Rank(u.red)>0 then
      PrintToString(c,".");
      if Size(u.Au)=1 or Size(u.Au)=Size(ApplyFunc(Group,u.AuAction.F0s)) then
	PrintToString(c,ReflectionName(u.AuAction,opt));
      elif ForAll(u.AuAction.F0s,x->x=x^0) then
	PrintToString(c,ReflectionName(u.AuAction.group,opt),AuName(u));
      else
	PrintToString(c,ReflectionName(u.AuAction,opt),AuName(u));
      fi;
    else
      PrintToString(c,AuName(u));
    fi;
  elif IsBound(u.red) then 
    n:=ReflectionName(u.red,opt);
    if n<>"." then PrintToString(c,".",n);fi;
    PrintToString(c,AuName(u));
  else
    if IsBound(u.dimred) then
      if u.dimred>0 then PrintToString(c,"[red dim ",u.dimred,"]");fi;
    else PrintToString(c,"[red??]");
    fi;
    PrintToString(c,AuName(u));
  fi;
  if Length(c)>1 and c[1]='.' then c:=c{[2..Length(c)]};fi;
  if Length(c)>1 and c[Length(c)]='.' then c:=c{[1..Length(c)-1]};fi;
  c:=Replace(c,"()","");
  return c;
end;

UnipotentClassOps.Display:=function(u,opt)local res;
  Print(Name(u),": ");
  if IsBound(u.dynkin) then Print("D-R",IntListToString(u.dynkin)," ");fi;
  Print("C=",UnipotentClassOps.FormatCentralizer(u,opt));
  Print("\n");
end;

#------------------------------------------------------------------------------
#  Unipotent classes records contain the following fields:
#    .spets  :  reductive group G of which this represents classes
#    .p      :  characteristic of G
#    .classes:  list of class records
#    .orderClasses: Hasse diagram of class poset
#    .springerSeries
UnipotentClassesOps:=OperationsRecord("UnipotentClassesOps");

UnipotentClassesOps.String:=function(r)local res;
  res:=SPrint("UnipotentClasses( ",ReflectionName(r.spets));
  if r.p<>0 then PrintToString(res,",",r.p);fi;
  Append(res," )");return res;
end;

UnipotentClassesOps.Print:=function(r)Print(String(r));end;

# QuotientAu(Au,chars): chars is a list of indices of characters of Au.
# If  k is the common kernel of chars, QuotientAu returns a record
# rec(Au:=Au/k,
#     chars:=index of chars as characters of Au/k,
#     gens:=words in Au preimages of generators of Au/k)
# Since  GAP3  has  many  problems  with  quotient groups, we are forced to
# program an ad hoc solution which works only for Au actually occuring for
# unipotent classes of a reductive group G.
QuotientAu:=function(Au,chars)local ct,cl,k,f,q,n,d,g,AbGens,finish,e,Z,h,p,t;
  AbGens:=function(g)local l,res,t;
    res:=[];l:=g.generators;
    while l<>[] do
      SortBy(l,x->-Order(g,x));
      t:=l[1]*Elements(Subgroup(g,res));
      if ForAny(t,x->Order(g,x)<Order(g,l[1])) then
        t:=First(t,x->Order(g,x)<Order(g,l[1]));
	if Order(g,t)>1 then l[1]:=t;
	else l:=Drop(l,1);
	fi;
      else Add(res,l[1]);
      fi;
    od;
    return res;
  end;
  # q=Au/k,  ww=words in q images of Au.generators
  finish:=function(q,ww)local h,fusion,ctu,cth;
    h:=GroupHomomorphismByImages(Au,q,Au.generators,List(ww,x->EltWord(q,x)));
    fusion:=List(ChevieClassInfo(Au).classtext,
       c->PositionClass(q,Image(h,EltWord(Au,c))));
    ctu:=CharTable(Au).irreducibles;
    cth:=CharTable(q).irreducibles;
    return rec(Au:=q,chars:=List(chars,
      c->Position(cth,List([1..NrConjugacyClasses(q)],
      j->ctu[c][Position(fusion,j)]))),
      gens:=List(q.generators,
           x->GetWord(Au,First(Elements(Au),y->Image(h,y)=x))));
  end;
  Z:=n->ComplexReflectionGroup(n,1,1);
  ct:=TransposedMat(CharTable(Au).irreducibles{chars});
  cl:=Filtered([1..Length(ct)],i->ct[i]=ct[1]);
  if Length(cl)=1 then return rec(Au:=Au,chars:=chars,
    gens:=List(Au.generatingReflections,x->[x]));fi;
  ct:=TransposedMat(Set(ct));
  k:=Subgroup(Au,Union(List(cl,i->Elements(ConjugacyClasses(Au)[i]))));
  if Size(k)=Size(Au) then
    return rec(Au:=CoxeterGroup(),chars:=[1],gens:=[]);fi;
  if Au.semisimpleRank=1 then return finish(Z(Size(Au)/Size(k)),[[1]]);
  elif IsAbelian(Au/k) then
    q:=Au/k;q.generators:=AbGens(q);h:=NaturalHomomorphism(Au,q);
    f:=List(Au.generators,x->GetWord(q,x^h));
 #  Print(Product(List(q.generators,x->Z(Order(q,x))))," ",f,"\n");
    return finish(Product(List(q.generators,x->Z(Order(q,x)))),f);
  else
    p:=PositionProperty(Au.type,t->ForAll(Elements(k),x->
                  x in ReflectionSubgroup(Au,t.indices)));
    if p<>false then 
      p:=Au.type[p].indices;
      if Size(k)=Size(ReflectionSubgroup(Au,p)) then
	return finish(
	 ReflectionSubgroup(Au,Difference(Au.generatingReflections,p)),
	 List(Au.generatingReflections,function(i)
	   if i in p then return [];else return [i];fi;end));
      elif Length(p)=1 then
        t:=Copy(Au.type);
        p:=PositionProperty(Au.type,t->t.indices=p);
	t[p].p:=t[p].p/Size(k);
        return finish(ApplyFunc(ReflectionGroup,t),
	  List(Au.generatingReflections,x->[x]));
      fi;
    elif ReflectionName(Au)="A1xB2" and Size(k)=2 
      and LongestCoxeterElement(Au) in k
    then return finish(CoxeterGroup("B",2),[[1,2,1,2],[1],[2]]);
    fi;
  fi;
# Print(" Au=",ReflectionName(Au)," sub=",List(k.generators,e.Get),"\n");
  Error("not implemented ",ReflectionName(Au),chars);
# q:=Au/k; f:=FusionConjugacyClasses(Au,q); Print(" quot=",q," fusion=",f,"\n");
# return rec(Au:=Au,chars:=chars);
end;

# When some Springer series have been suppressed/weeded out, we  quotient 
# the Aus by the common  kernel of the remaining characters of the Aus. 
AdjustAu:=function(ucl)local l,i,chars,u,f,j,k,s;
  ucl:=ShallowCopy(ucl); 
  for i in [1..Length(ucl.classes)] do
    l:=List(ucl.springerSeries,s->
       Filtered([1..Length(s.locsys)],k->s.locsys[k][1]=i));
    chars:=Concatenation(List([1..Length(l)],j->
      List(ucl.springerSeries[j].locsys{l[j]},x->x[2])));
    u:=ucl.classes[i];
    f:=QuotientAu(u.Au,chars);
#   if Size(Au)<>Size(f.Au) then
#     Print("class ",i,"=",ucl.classes[i].name," ",[Au,chars],"=>",f,"\n");
#   fi;
    u.Au:=f.Au;
    if IsBound(u.AuAction) then
      if u.AuAction.group.rank=0 then u.AuAction.F0s:=List(f.gens,x->[[]]);
      else u.AuAction.F0s:=List(f.gens,x->Product(u.AuAction.F0s{x}));
      fi;
      u.AuAction.phis:=List(f.gens,x->Product(u.AuAction.phis{x}));
    fi;
    k:=1;
    for j in [1..Length(l)] do 
      ucl.springerSeries[j].locsys:=Copy(ucl.springerSeries[j].locsys);
      for s in l[j] do ucl.springerSeries[j].locsys[s][2]:=f.chars[k];k:=k+1;od;
    od;
  od;
  return ucl;
end;

# UnipotentClasses(W [,characteristic])
UnipotentClasses:=function(arg)local obj,f;obj:=arg[1];
  if not IsRec(obj) then Error(obj," has no method for UnipotentClasses");fi;
  if IsCoxeterGroup(obj) then arg[1]:=Spets(obj);fi;
  if Length(arg)=1 then arg:=[arg[1],0];fi;
  f:=SPrint("unipotentclasses",arg[2]);
  if not IsBound(obj.(f)) then obj.(f):=
    ApplyFunc(Dispatcher("UnipotentClasses"),arg);fi;
  return obj.(f);
end;

ReflTypeOpsUnipotentClasses:=function(t,p)local s,uc,u,weights;
  s:=["UnipotentClasses",t];
  if IsBound(t.orbit) then t:=t.orbit[1];fi;
  if IsBound(t.cartanType) then Add(s,t.cartanType);fi;
  Add(s,p);
  uc:=ApplyFunc(CHEVIE.Data,s);
  if uc=false then return false;fi;
  for u in uc.classes do # fill omitted fields
    if not IsBound(u.parameter) then u.parameter:=u.name;fi;
    if IsBound(u.dynkin) then
      weights:=RootsCartan(CartanMat(t))*u.dynkin;
      p:=Number(weights,x->x=0);
      if IsBound(u.dimBu) and u.dimBu<>p+Number(weights,x->x=1)/2 then
        Error("theory");
      fi;
      u.dimBu:=p+Number(weights,x->x=1)/2;
      p:=2*p-Number(weights,x->x=2);
      u.dimunip:=2*u.dimBu-p;
      u.dimred:=p+t.rank;
    elif IsBound(u.red) then
      u.dimred:=Dimension(u.red);
      u.dimunip:=2*u.dimBu+t.rank-u.dimred;
    elif IsBound(u.dimred) then
      u.dimunip:=2*u.dimBu+t.rank-u.dimred;
    fi;
  od;
  return uc;
end;

UnipotentClassesOps.ClassName:=function(arg)local opt;
  if Length(arg)=2 then opt:=rec();else opt:=arg[3];fi;
  return Name(arg[1].classes[arg[2]],opt);
end;

# h  is a  linear form  defined by  its value  on the  simple roots  of the
# reflection subgroup K. Induce it to W by extending by 0 on the orthogonal
# of K, then conjugate it so it takes >=0 values on the simple roots.
InducedLinearForm:=function(W,K,h)local v,w,r;
# Print("W=",W," K=",K," h=",h,"\n");
  if SemisimpleRank(K)=0 then return W.generatingReflections*0;fi;
  h:=ShallowCopy(h);Append(h,[1..K.rank-K.semisimpleRank]*0);
  h:=CoxeterGroupOps.BaseX(Parent(W))*CoxeterGroupOps.BaseX(K)^-1*h;
  r:=Parent(W).roots{W.rootInclusion};
  v:=r{[1..W.N]}*h;
  w:=ElementWithInversions(W,Filtered([1..W.N],i->v[i]<0));
  return List(W.rootRestriction{OnTuples(
    W.rootInclusion{W.generatingReflections},w^-1)},i->r[i]*h);
end;

DistinguishedParabolicSubgroups:=W->Filtered(Combinations(W.rootInclusion{
  W.generatingReflections}),function(J)local p;
  if J=[] then return true;fi;
  p:=W.generatingReflections*0+1;p{W.rootRestriction{J}}:=J*0;
  p:=W.roots{[1..W.N]}*p;
  return 2*Number(p,x->x=0)+W.semisimpleRank=Number(p,x->x=1);end);

BalaCarterLabels:=function(W)local l;
  l:=Concatenation(List(ParabolicRepresentatives(W),
    J->List(DistinguishedParabolicSubgroups(ReflectionSubgroup(W,J)),
    D->[J,D])));
  return List(l,function(p)local L,w,i;
    L:=ReflectionSubgroup(W,p[1]);
    w:=L.generatingReflections*0+2;
    w{L.rootRestriction{p[2]}}:=p[2]*0;
    return [InducedLinearForm(W,L,w),
      List(p[1],function(i)if i in p[2] then return -i;else return i;fi;end)];
    end);
end;

HasTypeOpsUnipotentClasses:=function(WF,p)
  local t,tf,u,uc,ucl,ll,i,j,k,f,p,g,l,chars,Au,W,s,bc;
  W:=Group(WF);
  t:=ReflectionType(W); 
  tf:=ReflectionType(WF); 
  uc:=List(t,x->UnipotentClasses(First(tf,
     y->ForAny(y.orbit,z->Set(x.indices)=Set(z.indices))),p));
  if false in uc then return false;fi;
  if t=[] then ucl:=rec(classes:=[rec(name:="",Au:=CoxeterGroup(),
      parameter:=[],dimBu:=0,dynkin:=[],balacarter:=[],
      dimunip:=0,red:=Torus(W.rank),operations:=UnipotentClassOps)]);
  else
  ucl:=rec(classes:=List(Cartesian(List(uc,x->x.classes)),
    function(v)local u,p;
      u:=rec(name:=Join(List(v,x->x.name)),Au:=Product(v,x->x.Au),
         dimBu:=Sum(v,x->x.dimBu), parameter:=List(v,x->x.parameter));
      if Length(v)=1 then Inherit(u,v[1]);fi;
      if ForAll(v,x->IsBound(x.dimred)) then 
        u.dimred:=Sum(v,x->x.dimred);
      fi;
      if ForAll(v,x->IsBound(x.dimunip)) then 
        u.dimunip:=Sum(v,x->x.dimunip);
      fi;
      if ForAll(v,x->IsBound(x.red)) then 
        u.red:=Product(v,x->x.red);
      fi;
      if ForAll(v,x->IsBound(x.AuAction)) then 
        u.AuAction:=Product(v,x->x.AuAction);
      fi;
      l:=List(t,y->y.indices);
      if ForAll(v,x->IsBound(x.dynkin)) then
	u.dynkin:=[];
	for i in [1..Length(l)] do u.dynkin{l[i]}:=v[i].dynkin;od;
      fi;
      if ForAll(v,x->IsBound(x.balacarter)) then
        u.balacarter:=Concatenation(List(List([1..Length(l)],
	   i->List(v[i].balacarter,function(j)
	     if j>0 then return l[i][j];else return -l[i][-j];fi;end))));
      fi;
    if W.rank>W.semisimpleRank and IsBound(u.red) 
    then u.red:=u.red*Torus(W.rank-W.semisimpleRank);
      if IsBound(u.AuAction) then
        u.AuAction.group:=u.AuAction.group*Torus(W.rank-W.semisimpleRank);
	u.AuAction.F0s:=List(u.AuAction.F0s,x->DiagonalMat(x,
	  IdentityMat(W.rank-W.semisimpleRank)));
      fi;
    fi;
    u.operations:=UnipotentClassOps;
    return u;end));
  fi;
  ucl.size:=Length(ucl.classes);
  if p=0 and not IsBound(ucl.classes[1].balacarter) then
    bc:=BalaCarterLabels(W);
    for u in ucl.classes do 
      u.balacarter:=bc[PositionProperty(bc,p->p[1]=u.dynkin)][2];
    od;
  fi;
  ll:=List(uc,x->Length(x.classes));
  ucl.orderClasses:=List(Cartesian(List(ll,x->[1..x])),function(v)local o;
    o:=Cartesian(List([1..Length(v)],j->
      Concatenation(uc[j].orderClasses[v[j]],[v[j]])));
    o:=List(o,x->PositionCartesian(ll,x));
    return Difference(o,[PositionCartesian(ll,v)]);
    end);
  ucl.springerSeries:=List(Cartesian(List(uc,x->x.springerSeries)),
    function(v)local s;
      if v=[] then return rec(Z:=[],levi:=[],locsys:=[[1,1]]);
      elif Length(v)=1 then return Copy(v[1]);
      fi;
      s:=rec(levi:=Concatenation(List([1..Length(v)],
         i->l[i]{v[i].levi})));
      s.Z:=Concatenation(List(v,x->x.Z));
      s.locsys:=List(Cartesian(List(v,x->x.locsys)),function(v)
	  v:=TransposedMat(v);
	  v[2]:=PositionCartesian(List([1..Length(v[1])],i->
	    NrConjugacyClasses(uc[i].classes[v[1][i]].Au)),v[2]);
	  v[1]:=PositionCartesian(ll,v[1]);
	  return v;end);
      if ForAll(v,x->IsBound(x.parameter)) then
        s.parameter:=List(v,x->x.parameter);
      fi;
  #   if Length(v)=1 then Inherit(s,v[1],Difference(RecFields(v[1]),
  #["levi","Z","locsys","parameter"]));fi;
      return s;end);
# adjust indices of levi, relativetype so they agree with Parent(Group(WF))
  for s in ucl.springerSeries do
    s.levi:=W.rootInclusion{s.levi};
  od;
  ucl.operations:=UnipotentClassesOps;
  ucl.spets:=WF;
  ucl.p:=p;
  if Length(uc)=1 then Inherit(ucl,uc[1],Difference(RecFields(uc[1]),
    ["group","operations","springerSeries","classes","orderClasses"]));fi;
  ucl.springerSeries:=Filtered(ucl.springerSeries,
    x->OnSets(Set(x.levi),WF.phi)=Set(x.levi));
# To deal with a general group intermediate between Gad and Gsc, we discard
# the  Springer series  corresponding to  a central  character which is not
# trivial on the fundamental group (seen as a subgroup of ZGsc)
# AlgebraicCentre(W).descAZ returns the generators of the fundamental group
# of  the  algebraic  group  W  as  words  in  generators  of  the absolute
# fundamental group.
  if not ForAll(ucl.springerSeries,x->Set(x.Z)=[1]) then 
    ucl.springerSeries:=Filtered(ucl.springerSeries,
       s->ForAll(AlgebraicCentre(W).descAZ,y->Product(s.Z{y})=1));
    ucl:=AdjustAu(ucl); 
  fi;
  s:=ucl.springerSeries[1];
  s.relgroup:=RelativeCoset(WF,s.levi);
  s.locsys:=s.locsys{ChevieCharInfo(s.relgroup).charRestrictions};
  l:=Filtered([1..Length(ucl.classes)],i->ForAny(s.locsys,y->i=y[1]));
  s.locsys:=List(s.locsys,y->[Position(l,y[1]),y[2]]);
  # for now only Springerseries[1] properly twisted
  for s in ucl.springerSeries{[2..Length(ucl.springerSeries)]} do
    s.relgroup:=RelativeCoset(WF,s.levi);
    s.locsys:=s.locsys{ChevieCharInfo(s.relgroup).charRestrictions};
    s.locsys:=List(s.locsys,y->[Position(l,y[1]),y[2]]);
  od;
  ucl.classes:=ucl.classes{l};
  AdjustAu(ucl);
  ucl.orderClasses:=Hasse(Restricted(Poset(ucl.orderClasses),l));
  ucl.size:=Length(l);
  return ucl;
end;

UnipotentClassesOps.Poset:=function(uc)local res;
  res:=Poset(uc.orderClasses); res.uc:=uc;
  res.label:=function(p,n,opt)return Name(p.uc.classes[n],opt);end;
  return res;
end;

UnipotentClassesOps.DisplayOptions:=rec(
 order:=true,springer:=true,centralizer:=true,balaCarter:=true);
#order:=false,springer:=false);

UnipotentClassesOps.Format:=function(uc,opt)local W,sp,tbl,res,p,TeX,c;
  opt:=Inherit(ShallowCopy(UnipotentClassesOps.DisplayOptions),opt); 
  TeX:=function(a,b)if IsBound(opt.TeX) then return a;else return b;fi;end;
  opt.rowLabels:=List(uc.classes,i->Name(i,opt));
  if opt.order then res:=Format(Poset(uc),opt);else res:="";fi;
  sp:=List(uc.springerSeries,ShallowCopy);
  if IsBound(opt.fourier) then
    for p in sp do p.locsys:=p.locsys{DetPerm(p.relgroup)};od;
  fi;
  W:=Group(uc.spets);
  tbl:=List(uc.classes,function(u)local i,res,b;
    if uc.p=0 then res:=[IntListToString(u.dynkin)];else res:=[];fi;
    Add(res,u.dimBu);
    if uc.p=0 and opt.balaCarter then 
      if IsBound(u.balacarter) then
	b:=List(W.generatingReflections,x->'.');
	for i in Filtered(u.balacarter, x->x>0) do b[i]:='2';od;
	for i in Filtered(u.balacarter, x->x<0) do b[-i]:='0';od;
      else b:=List(W.generatingReflections,x->'?');
      fi;
      Add(res,String(b));
    fi;
    if opt.centralizer then 
      Add(res,UnipotentClassOps.FormatCentralizer(u,opt));
    fi;
    if opt.springer then
      i:=Position(uc.classes,u);
      Append(res,List(sp,ss->Join(List(PositionsProperty(ss.locsys,y->y[1]=i),
	function(i)local c1,c2;
	  c1:=CharNames(u.Au,opt)[ss.locsys[i][2]];
	  c2:=CharNames(ss.relgroup,opt)[i];
	  if c1="" then return c2;else return SPrint(c1,":",c2);fi;
	  end),TeX("\\kern 0.8em "," "))));
    fi;
    return res;end);
  opt.rowsLabel:="u";opt.columnLabels:=[];
  if uc.p=0 then 
    Add(opt.columnLabels,TeX("\\hbox{Dynkin-Richardson}","D-R"));
  fi;
  Add(opt.columnLabels,TeX("\\dim{\\cal B}_u","dBu"));
  if uc.p=0 and opt.balaCarter then 
    Add(opt.columnLabels,TeX("\\hbox{Bala-Carter}","B-C"));
  fi;
  if opt.centralizer then
    Add(opt.columnLabels,TeX("C_{\\bf G}(u)","C(u)"));
  fi;
  if opt.springer then
    Append(opt.columnLabels,List(sp,function(ss)local res;
    res:=SPrint(ReflectionName(ss.relgroup,opt),
      "(",IsomorphismType(ReflectionSubgroup(W,ss.levi),opt),")");
    if not ForAll(ss.Z,x->x=1) then PrintToString(res,"/",
        Join(List(ss.Z,x->Format(x,opt))));fi;return res;end));
  fi;
  if not IsBound(opt.rows) then
    p:=SortingPerm(List(uc.classes,x->x.dimBu));
    tbl:=Permuted(tbl,p);opt.rowLabels:=Permuted(opt.rowLabels,p);
  fi;
  Append(res,FormatTable(tbl,opt));
  return res;
end;

UnipotentClassesOps.Display:=function(uc,opt)
  opt:=ShallowCopy(opt);opt.screenColumns:=SizeScreen()[1];
  Print(Format(uc,opt));
end;

# coefficients of R_\chi on unipotently supported local systems 
# ICCTable(uc[,Springer series no[,variable]]) eg (uc,1,X(Rationals))
# Works for G split.
ICCTable:=function(arg)local W,i,q,tbl,o,res,q,uc,ss,b,f,k,R,n;
  uc:=arg[1];W:=Group(uc.spets);
  if Length(arg)<2 then i:=1;else i:=arg[2];fi;
  q:=Indeterminate(Rationals);
  ss:=uc.springerSeries[i];
  res:=rec(spets:=uc.spets,relgroup:=ss.relgroup,series:=i,q:=q,p:=uc.p);
  if IsBound(ss.warning) then Print("# ",ss.warning,"\n");
    res.warning:=ss.warning;fi;
# We are going to solve the equation in "unipotent support", page 151
# $Transposed(P)\Lambda P=\omega$
# where $\Lambda_{i,j}$ is  $\sum_{g\in G^F} Y_i(g)\overline{Y_j(g)}$
# and $\Omega_{i,j}$ is equal to
# $|Z^0(G^F)|q^{-\text{semisimple rank}L}|G^F|/P(W_G(L))
#  q^{-b_i-b_j}FakeDegree(\chi_i\otimes\chi_j\otimes\sgn)$
# where $P(W_G(L))$ is the Poincare polynomial $\prod_i(q^{d_i}-1)$
# where $d_i$ are the reflection degrees of $W_G(L)$
# res.scalar is the matrix $P$
  R:=ss.relgroup; f:=FakeDegrees(R,q);k:=PositionDet(R);
  n:=Length(f);
# Partition on characters of ss.relgroup induced by poset of unipotent classes
  res.dimBu:=List(ss.locsys,x->uc.classes[x[1]].dimBu);
  res.blocks:=CollectBy([1..Length(ss.locsys)],-res.dimBu);
  tbl:=BigCellDecomposition(List([1..n],i->List([1..n],
  # q^{-b_i-b_j}*matrix of FakeDegrees(chi_i tensor chi_j tensor sgn)
     j->q^(-res.dimBu[i]-res.dimBu[j])*f*DecomposeTensor(R,i,j,k))),
   res.blocks);
  res.scalar:=tbl[1];
  res.locsys:=ss.locsys;
# res.L:=tbl[2]*GenericOrder(W,q)/Product(ReflectionDegrees(R),d->q^d-1)/
#   q^(W.semisimpleRank-R.semisimpleRank);
  res.L:=tbl[2]*q^(W.N+SemisimpleRank(R)-SemisimpleRank(W));
  res.uc:=uc;
  if IsBound(ss.parameter) then res.parameter:=ss.parameter;
  else res.parameter:=[1..Length(ss.locsys)];
  fi;
  if Length(arg)=3 and q<>arg[3] then
    q:=arg[3];
    res.scalar:=List(res.scalar,x->List(x,y->Value(y,q)));
    res.L:=List(res.L,x->List(x,y->Value(y,q)));
  fi;
  res.operations:=rec(
    String:=x->SPrint("ICCTable(",W,",",i,",",q,")"),
    Print:=function(x)Print(String(x));end,
    Format:=function(x,opt)local tbl,res;
     if IsBound(opt.TeX) then
     res:=SPrint("Coefficients of $X_\\phi$ on $Y_\\psi$ for $",
        ReflectionName(x.relgroup,opt),"$\n\\medskip\n\n");
     else
     res:=SPrint("Coefficients of X_phi on Y_psi for ",
        ReflectionName(x.relgroup,opt),"\n\n");
     fi;
     if not IsBound(opt.columns) and not IsBound(opt.rows) then 
       opt.rows:=[1..Length(x.dimBu)];
       SortBy(opt.rows,i->[x.dimBu[i],x.locsys[i]]);
       opt.columns:=opt.rows;
     fi;
     tbl:=Copy(x.scalar);
     if not IsBound(opt.CycPol) then opt.CycPol:=true;fi;
     if opt.CycPol then tbl:=List(tbl,x->List(x,
       function(p)p:=CycPol(p);p.vname:="q";return p;end));fi;
     opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
         Inherit(rec(locsys:=p[2]),opt)));
     opt.rowLabels:=List(CharNames(x.relgroup,opt),
       function(x)if IsBound(opt.TeX) then return SPrint("X_{",x,"}");
         else return SPrint("X",x);fi;end);
     PrintToString(res,FormatTable(TransposedMat(tbl),opt));
     return String(res);
     end,
    Display:=function(x,opt)opt.screenColumns:=SizeScreen()[1];
               Print(Format(x,opt));end);
  return res;
end;

# Green functions: Green(uc[,opt]) values on unipotent classes or local systems
# opt: variable (default X(Cyclotomics))
#
# Formatting: options of FormatTable + [.classes, .CycPol]
GreenTable:=function(arg)
  local opt,uc,W,res,pieces,l,m,n,p,q,greenpieces,Au,i,b;
  uc:=arg[1];W:=Group(uc.spets);
  if Length(arg)=1 then opt:=rec();else opt:=arg[2];fi;
  if not IsBound(opt.variable) then opt.variable:=X(Cyclotomics);fi;
  q:=opt.variable;
  pieces:=List([1..Length(uc.springerSeries)],i->ICCTable(uc,i,q));
  greenpieces:=List(pieces,x->List(x.scalar,
   l->Zip(l,x.dimBu,function(x,y)return x*q^y;end)));
  l:=Concatenation(List(pieces,x->x.locsys));
  p:=SortingPerm(l);
  res:=rec(
    scalar:=TransposedMat(Permuted(ApplyFunc(DiagonalMat,greenpieces),p)),
    uc:=uc,
    Y:=OnMatrices(ApplyFunc(DiagonalMat,List(pieces,p->p.L)),p),
    locsys:=Permuted(l,p),
    parameter:=Concatenation(List(pieces,x->x.parameter)),
    relgroups:=List(uc.springerSeries,x->x.relgroup));
  n:=Length(res.locsys);
  if IsBound(opt.classes) then
    res.cardClass:=[];
    for i in [1..Length(uc.classes)] do
      Au:=uc.classes[i].Au;
      b:=Filtered([1..Length(res.locsys)],j->res.locsys[j][1]=i);
      res.scalar{[1..n]}{b}:=res.scalar{[1..n]}{b}*CharTable(Au).irreducibles;
      res.cardClass{b}:=res.Y[b[PositionId(Au)]]{b}*CharTable(Au).irreducibles;
      res.cardClass{b}:=Zip(res.cardClass{b},ChevieClassInfo(Au).classes,
        function(x,y)return x*y/Size(Au);end);
    od;
    res.classes:=true;
  fi;
  res.operations:=rec();
  res.operations.String:=x->SPrint("GreenTable(",W,",rec(variable:=",q,"))");
  res.operations.Format:=function(x,opt)local res,tbl,i,b;
    res:=SPrint("Values of character sheaves on");
    opt.rowLabels:=Concatenation(List(x.relgroups,g->
      List(CharNames(g,opt),n->SPrint("Q^{",ReflectionName(g),"}_{",n,"}"))));
    opt.rowsLabel:="Q";
    tbl:=Copy(x.scalar);
    if IsBound(x.classes) then
      PrintToString(res," unipotent classes\n");
      opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
	  Inherit(rec(class:=p[2]),opt)));
    else PrintToString(res," local systems\n");
      opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
        Inherit(rec(locsys:=p[2]),opt)));
    fi;
    if not IsBound(opt.CycPol) then opt.CycPol:=true;fi;
    if opt.CycPol then tbl:=List(tbl,y->List(y,CycPol));fi;
    PrintToString(res,FormatTable(tbl,opt));
    return String(res);
  end;
  res.operations.Print:=function(x)Print(String(x));end;
  res.operations.Display:=function(x,opt)opt.screenColumns:=SizeScreen()[1];
               Print(Format(x,opt));end;
  return res;
end;

# values of unipotent characters
# UnipotentValues(uc[,opt]) values on unipotent classes or local systems
UnipotentValues:=function(arg)local uc,opt,res,q,W;
  uc:=arg[1];W:=Group(uc.spets);
  if Length(arg)=1 then opt:=rec();else opt:=arg[2];fi;
  if not IsBound(opt.variable) then opt.variable:=X(Cyclotomics);fi;
  q:=opt.variable;
  res:=Copy(GreenTable(uc,opt));
  uc:=UnipotentCharacters(W);
  res.scalar:=TransposedMat(uc.operations.Fourier(uc){res.parameter})
    *res.scalar;
  res.operations.String:=x->SPrint("UnipotentValues(",W,",",q,")");
  res.operations.Format:=function(x,opt)local s,i,b,tbl;
    s:="Values of unipotent characters for ";
    if IsBound(opt.TeX) then 
         PrintToString(s,"$",ReflectionName(W,opt),"$\\par");
    else PrintToString(s,ReflectionName(W,opt));
    fi;
    tbl:=Copy(x.scalar);
    if IsBound(x.classes) then
      PrintToString(s," on unipotent classes\n");
      opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
	  Inherit(rec(class:=p[2]),opt)));
    else opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
        Inherit(rec(locsys:=p[2]),opt)));
      PrintToString(s," on local systems\n");
    fi;
    opt.rowLabels:=CharNames(uc,opt);
    if not IsBound(opt.CycPol) then opt.CycPol:=true;fi;
    if opt.CycPol then tbl:=List(tbl,y->List(y,CycPol));fi;
    return Concatenation(s,FormatTable(tbl,opt));
  end;
  return res;
end;

# values of mellin  on unipotent classes
MellinValues:=function(arg)local uc,q,res,b,tbl,f,labels,W;
  uc:=arg[1];if Length(arg)=1 then q:=X(Cyclotomics); else q:=arg[2];fi;
  res:=Copy(GreenTable(uc,q));
  W:=Group(uc.spets);
  uc:=UnipotentCharacters(W);
  uc.mellin:=IdentityMat(Size(uc));
  labels:=[];
  for f in uc.families do
    uc.mellin{f.charNumbers}{f.charNumbers}:=f.mellin;
    labels{f.charNumbers}:=List(f.mellinLabels,
      x->SPrint(Position(uc.families,f),x));
  od;
  res.scalar:=TransposedMat(uc.operations.Fourier(uc){res.parameter})
    *res.scalar;
  res.scalar:=uc.mellin*res.scalar;
  res.operations.String:=x->SPrint("MellinValues(",W,",",q,")");
  res.operations.Format:=function(x,opt)local s,i,b,tbl;
    s:="Values of Mellin of unipotent characters for ";
    if IsBound(opt.TeX) then 
         PrintToString(s,"$",ReflectionName(W,opt),"$\\par");
    else PrintToString(s,ReflectionName(W,opt));
    fi;
    tbl:=Copy(x.scalar);
    if IsBound(opt.classes) then
      PrintToString(s," on unipotent classes\n");
      for i in [1..Length(x.uc.classes)] do
	b:=Filtered([1..Length(x.locsys)],j->x.locsys[j][1]=i);
	tbl{[1..Length(tbl)]}{b}:=tbl{[1..Length(tbl)]}{b}*
	   CharTable(x.uc.classes[i].Au).irreducibles;
      od;
      opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
	  Inherit(rec(class:=p[2]),opt)));
    else opt.columnLabels:=List(x.locsys,p->Name(x.uc.classes[p[1]],
        Inherit(rec(locsys:=p[2]),opt)));
      PrintToString(s," on local systems\n");
    fi;
    opt.rowLabels:=labels;
    if not IsBound(opt.CycPol) then opt.CycPol:=true;fi;
    if opt.CycPol then tbl:=List(tbl,y->List(y,CycPol));fi;
    return Concatenation(s,FormatTable(tbl,opt));
  end;
  return res;
end;

# returns the Special pieces of the nilpotent cone as a list of 
# lists of class numbers.
# The list is sorted by increasing piece dimension.
# Each piece is sorted by decreasing class dimension.
# The argument is a unipotent classes record for some Weyl group W.
# The special pieces were first defined in
# Spaltenstein "Classes unipotentes et sous-sgroupes de Borel" as the
# fibers of d^2 where d is the "duality map" of chapter III
SpecialPieces:=function(uc)local W,specialc,m,ch,specialch;
  W:=Group(Spets(uc));
  ch:=ChevieCharInfo(W);
  specialch:=Positions(ch.a-ch.b,0); # special characters of W
  specialc:=List(uc.springerSeries[1].locsys{specialch},x->x[1]);
  SortBy(specialc,c->-uc.classes[c].dimBu);
  m:=TransposedMat(Incidence(Poset(uc)));
  return List([1..Length(specialc)],function(i)local p,j;
    p:=m[specialc[i]];
    for j in [1..i-1] do SubtractBlist(p,m[specialc[j]]);od;
    p:=ListBlist([1..Length(p)],p);SortBy(p,c->uc.classes[c].dimBu);
    return p;
  end);
end;
