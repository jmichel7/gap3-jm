#############################################################################
##
#A  hastype.g  CHEVIE library    Jean Michel  5-1999
##
#Y  Copyright (C) 1999 Lehrstuhl D fur Mathematik, RWTH Aachen,
#Y   and   University Paris VII.
##
##  This  file contains  operations  for reflexion  groups or  reflexion
##  cosets  which are  classified  (whose  decomposition in  irreducible
##  components is known). It builds answers by combining the answers for
##  irreducible types fetched from the CHEVIE tables. Here 'irreducible'
##  includes descents of  scalars. The code for irreducible  types is in
##  refltype.g .
##
##  A group or coset is classified  when its field .type has been filled
##  This field is accessed by the function 'ReflectionType' .

# ReflectionGroup from a ReflectionType
ReflectionGroup:=function(arg)local t,g,res,o,i,r;
  if Length(arg)=0 then return CoxeterGroup();fi;
  for t in arg do
    if t.series<>"ST" then g:=CoxeterGroup(t);
    elif IsBound(t.ST) then g:=ComplexReflectionGroup(t.ST);
    else g:=ComplexReflectionGroup(t.p,t.q,t.rank);
    fi;
    if t.series="ST" and IsBound(t.cartanType)then
      g.roots:=ShallowCopy(g.roots);
      g.simpleCoroots:=ShallowCopy(g.simpleCoroots);
      o:=Set(g.orbitRepresentative);
      if Length(o)=2 then
        o:=Filtered([1..Length(g.roots)],i->g.orbitRepresentative[i]=o[2]);
	g.roots{o}:=g.roots{o}*t.cartanType;
	o:=Intersection(o,[1..Length(g.simpleCoroots)]);
	g.simpleCoroots{o}:=g.simpleCoroots{o}/t.cartanType;
      else for i in [2..Length(o)] do
        r:=Filtered([1..Length(g.roots)],j->g.orbitRepresentative[j]=o[i]);
	g.roots{r}:=g.roots{r}*t.cartanType[i-1];od;
	r:=Intersection(r,[1..Length(g.simpleCoroots)]);
	g.simpleCoroots{r}:=g.simpleCoroots{r}/t.cartanType[i-1];
      fi;
      g.matgens:=List(g.generatingReflections,
                  j->Reflection(g.roots[j],g.simpleCoroots[j]));
    fi;
    if IsBound(res) then res:=res*g;
    else res:=g;
    fi;
  od;
  return res;
end;

HasTypeOps:=OperationsRecord("HasTypeOps");

###########################################################################
##
#F  CharName(W,p[,option]) . . . . . . . . Name of character with .charparam p
#    option may be "TeX"
##
HasTypeOps.CharName:=function(W,p,option)local t;
  t:=ReflectionType(W);
  return Join(List([1..Length(t)],i->CharName(t[i],p[i],option)));
end;

HasTypeOps.CharNames:=function(W,option)
  return List(ChevieCharInfo(W).charparams,x->CharName(W,x,option));
end;

HasTypeOps.NrConjugacyClasses:=W->Product(List(ReflectionType(W),
  NrConjugacyClasses));

#############################################################################
##
#F  HasTypeOps.ConjugacyClasses( <W> ) . . . . . . . ConjugacyClasses,
#F  using the classification
##  
HasTypeOps.ConjugacyClasses:=function(W)local t,cl;
  t:=ChevieClassInfo(W);
  if IsSpets(W) then
    if t=false then return SpetsOps.ConjugacyClasses(W);
    else cl:=List(t.classtext,x->ConjugacyClass(Group(W),
                                            EltWord(Group(W),x)*W.phi));
    fi;
  else
    if t=false then return PermGroupOps.ConjugacyClasses(W);
    else cl:=List(t.classtext,w->ConjugacyClass(W,EltWord(W,w)));
    fi;
  fi;
  if IsBound(t.classes) then 
    return Zip(cl,t.classes,function(a,n)a.size:=n;return a;end);
  else return cl;
  fi;
end;

HasTypeOps.ClassName:=function(arg)local opt,n;
  n:=ChevieClassInfo(arg[1]).classnames[arg[2]];
  if Length(arg)=2 then opt:=rec();else opt:=arg[3];fi;
  if IsBound(opt.TeX) then return n;else return TeXStrip(n);fi;
end;

#############################################################################
##
#F  HasTypeOps.ChevieClassInfo( <W> ) . . . info on conjugacy classes
##
HasTypeOps.ChevieClassInfo:=function(W) local tmp,res,inc,get;
  if IsSpets(W) then inc:=Group(W).rootInclusion;else inc:=W.rootInclusion;fi;
  tmp:=List(ReflectionType(W),ChevieClassInfo);

  get:=f->Cartesian(List(tmp,x->x.(f)));
  
  if ForAny(tmp,x->x=false) then return false;fi;
  if Length(tmp)=1 then res:=tmp[1];# keep extra fields for irreducible type
  else res:=rec();
  fi; 
  res.classtext:=List(get("classtext"),x->inc{Concatenation(x)});
  res.classnames:=List(get("classnames"),Join);
  if Length(tmp)=0 then
     res.classparams:=[[]];res.orders:=[1];res.classes:=[1];
  else
    if ForAll(tmp,x->IsBound(x.classparams)) then 
       res.classparams:=get("classparams");
    fi;
    if ForAll(tmp,x->IsBound(x.orders)) then 
      res.orders:=List(get("orders"),Lcm);
    fi;
    if ForAll(tmp,x->IsBound(x.classes)) then 
      res.classes:=List(get("classes"),Product);
    fi;
  fi;
  return res;
end;

#############################################################################
##
#F  HasTypeOps.BraidRelations( <W> ) . . . see corresponding dispatcher 
#F  function
##
HasTypeOps.BraidRelations:=function(W)local res,type,i,j;
  type:=ReflectionType(W);
  res:=Concatenation(List(type,t->List(BraidRelations(t),
    x->List(x,y->W.rootInclusion{y}))));
  if Length(type)>1 then for p in Combinations(type,2) do
    for i in p[1].indices do for j in p[2].indices do
      Add(res,[W.rootInclusion{[i,j]},W.rootInclusion{[j,i]}]);od;od;od;
  fi;
  return res;
end;
  
#############################################################################
##
#F  HasTypeOps.CharTable( <W> ) . . . . . . . . . . . 'CharTable', using
#F  the classification
##
HasTypeOps.CharTable:=function(W)local t, l, tbl, cl, d;

  l:=List(ReflectionType(W),ReflTypeOps.CharTable);
  if ForAny(l,tbl->tbl=false) then 
    if IsSpets(W) then tbl:=SpetsOps.CharTable(W);
    else tbl:=PermGroupOps.CharTable(W);
    fi;
  else
    if Length(l)=0 then 
      tbl:=rec(size:=1,order:=1,centralizers:=[1],orders:=[1],
        irreducibles:=[[1]],powermap:=[],operations:=CharTableOps);
    else tbl:=l[1];
      for t in l{[2..Length(l)]} do
	tbl:=CharTableDirectProduct(tbl,t); 
	Unbind(tbl.fusionsource); Unbind(tbl.fusions);
      od;
    fi;
  fi;
  
  tbl.irredinfo:=List(ChevieCharInfo(W).charparams,
              x->rec(charparam:=x,charname:=CharName(W,x,rec(TeX:=true))));

  # for a coset, InitClassesCharTable which is called from
  # CharTableDirectProduct uses a wrong value of size (centralizers[1])
  # so we must fix the classes
  tbl.classes:=List(tbl.centralizers,x->tbl.size/x);

  t:=ChevieClassInfo(W);
  if t<>false then Inherit(tbl,t);
    tbl.classnames:=List(tbl.classnames,TeXStrip);
  fi;

  if not IsSpets(W) and not IsBound(tbl.powermap) then
    Print("warning: powermap not bound\n");
    tbl.powermap:=[];
    cl:=List(ConjugacyClasses(W),Representative);
    for d in Reversed(Set(Factors(tbl.size))) do
    tbl.powermap[d]:=List(cl,x->PositionClass(W,x^d));
    od;
  fi;

  tbl.operations.StringEntry := function(x)
    if x=0*x then return ".";else return Format(x);fi;end;

  tbl.name:=ReflectionName(W);
  tbl.identifier:=tbl.name;
  return tbl;
end;

#############################################################################
##
#F  HasTypeOps.Representation(<W>, i) . .  'Representation', using
#F  the classification
##
HasTypeOps.Representation:=function(W,i)local l,t,reps,i,res,ind;
  t:=ReflectionType(W);
  ind:=CartesianAt(List(t,NrConjugacyClasses),i);
  reps:=List([1..Length(ind)],i->ReflTypeOps.Representation(t[i],ind[i]));
  if ForAny(reps,r->r=false) then return false;fi;
  if Length(t)=1 then return reps[1];fi;
  if IsSpets(W) then
    res:=rec(gens:=[],F:=ApplyFunc(KroneckerProduct,List(reps,x->x.F)));
    for i in [1..Length(t)] do
      res.gens{Concatenation(List(t[i].orbit,x->x.indices))}:=
       List(reps[i].gens,m->ApplyFunc(KroneckerProduct,List([1..Length(t)],
       function(k)if k<>i then return reps[k].gens[1]^0;else return m;fi;end)));
    od;
  else
    res:=[];
    for i in [1..Length(t)] do
      res{t[i].indices}:=List(reps[i],m->ApplyFunc(KroneckerProduct,
       List([1..Length(t)],
	function(k)if k<>i then return reps[k][1]^0;else return m;fi;end)));
    od;
  fi;
  return res;
end;

#############################################################################
##
#F  HasTypeOps.Representations( <W> ) . .  'Representations', using
#F  the classification
##
HasTypeOps.Representations:=function(arg)local W,inds;
  W:=arg[1]; 
  if Length(arg)=2 then inds:=arg[2];else inds:=[1..NrConjugacyClasses(W)];fi;
  if IsList(inds) then return List(inds,i->Representation(W,i));
  else return Representation(W,inds);
  fi;
end;

#############################################################################
##
#F  HasTypeOps.WGraph(<W>, i) . .  using the classification
##
HasTypeOps.WGraph:=function(W,i)local t;
  t:=ReflectionType(W); 
  if Length(t)=1 then return ReflTypeOps.WGraph(t[1],i);fi;
  Error("not implemented for non-irreducible types");
end;

############################################################################
##
#F  FakeDegree( <W>, <p>, <q> ) . . .Fake Degree of char with charparam <p>
##  
##  This   returns the   polynomial  describing the   multiplicity of  the
##  character chi with .charparam  p in the graded  version of the regular
##  representation  given  by the quotient  S/I  where S  is the symmetric
##  algebra of the reflection representation and  I is the ideal generated
##  by the homogenous invariants of positive degree in S.
##  
FakeDegree:=function(W,p,q)local t;
  t:=ReflectionType(W);
  t:=List([1..Length(t)],i->ReflTypeOps.FakeDegree(t[i],p[i],q));
  if ForAny(t,x->x=false) then return false;
  else return Product(t)*q^0;
  fi;
end;

############################################################################
##
#F  FakeDegrees( <W>, <q> ) . . . . . .  Fake Degrees of group <W>
##
##  This returns the list of polynomials describing the multiplicity of
##  each character in the graded version of the regular representation given
##  by the quotient S/I where S is the symmetric algebra of the reflection
##  representation and I is the ideal generated by the homogenous invariants
##  of positive degree in S.
##  The ordering of the result corresponds to the ordering of the characters 
##  in the  CharTable.
##
#
HasTypeOps.FakeDegrees:=function(W,q)local P,p,f;
  P:=[];
  for p in CharParams(W) do
    f:=FakeDegree(W,p,q);
    if f=false then 
      if IsSpets(W) then return SpetsOps.FakeDegrees(W,q);
      else return PermRootOps.FakeDegrees(W,q);
      fi;
    fi;
    Add(P,f);
  od;
  return P;
end;

#############################################################################
##
#F  ReflectionDegrees( <W>)  . . . . . . . degrees as a reflection group
##
##  Returns the degrees of the  reflection group W.
##  
HasTypeOps.ReflectionDegrees:=function(W)local d,ir,M,WF;
  d:=Concatenation(List(ReflectionType(W),ReflectionDegrees));
  if not IsSpets(W) then 
    return Concatenation([W.semisimpleRank+1..W.rank]*0+1,d);
  fi;
  WF:=W;W:=Group(WF);
  if W.rank=W.semisimpleRank then return d;fi;
  # prepend the (factors=eigenvalues of phi) where W acts trivially
  ir:=PermRootOps.IndependentRoots(W);
  M:=W.operations.BaseX(W);
  M:=M*WF.F0Mat*M^-1;
  ir:=[W.semisimpleRank+1..W.rank];
  return Concatenation(List(EigenvaluesMat(M{ir}{ir}),x->[1,x]),d);
end;

HasTypeOps.ReflectionCoDegrees:=function(W)local WF,ir,M,d;
  d:=Concatenation(List(ReflectionType(W),ReflectionCoDegrees));
  if not IsSpets(W) then 
    return Concatenation([W.semisimpleRank+1..W.rank]*0-1,d);
  fi;
  WF:=W;W:=Group(WF);
  if W.rank=W.semisimpleRank then return d;fi;
  # prepend the (cofactors=eigenvalues of phi^-1) where W acts trivially
  ir:=PermRootOps.IndependentRoots(W);
  M:=W.operations.BaseX(W);
  M:=M*WF.F0Mat*M^-1;
  ir:=[W.semisimpleRank+1..W.rank];
  return Concatenation(List(ComplexConjugate(EigenvaluesMat(M{ir}{ir})),
     x->[-1,x]),d);
end;

HasTypeOps.Size:=function(W)
  if IsSpets(W) then return Product(ReflectionDegrees(W),x->x[1]);
  else return Product(ReflectionDegrees(W));
  fi;
end;

#############################################################################
##
#F  HasTypeOps.PrintDiagram( <W> )  
##  
##  For a Spets information 
##  about the action of <WF>.phi on the Dynkin diagram is additionally printed.
##  
HasTypeOps.PrintDiagram:=function(W)
  PrintDiagram(List(ReflectionType(W),function(t)
    t:=ShallowCopy(t);
    if IsSpets(W) then
      t.orbit:=List(t.orbit,function(a)
	a:=ShallowCopy(a); 
        if ForAll(a.indices,i->IsBound(Group(W).reflectionsLabels[i])) then
	  a.indices:=Group(W).reflectionsLabels{a.indices};
        else
	  a.indices:=Group(W).rootInclusion{a.indices};
        fi;
	return a;end);
      if t.twist<>() then
	t.twist:=RestrictedPerm(W.phi^Length(t.orbit),t.orbit[1].indices);
      fi;
    else
      if ForAll(t.indices,i->IsBound(W.reflectionsLabels[i])) then
        t.indices:=W.reflectionsLabels{t.indices};
      else
        t.indices:=W.rootInclusion{t.indices};
      fi;
    fi;
    return t;end));
end;

#############################################################################
##
#F  HasTypeOps.ReflectionName(<W>)  see dispatcher function
##  Gives a string which describes the isomorphism type of W.  
##  
##  For a Spets, An orbit of
##  phi on the  components is put in brackets  if of length x greater than
##  1, and  is  preceded by  the  order of  phi^x on it,  e.g. 2(A2xA2xA2)
##  denotes 3 components A2  permuted by phi,  and such that phi^3 induces
##  the non-trivial diagram automorphism on any of them, while 3D4 denotes
##  an orbit of  length 1 on  which phi is of order  3.  
##    If the coset is not semi-simple, its toral part is added as a
##  product of cyclotomic factors over the minimum possible field.
##  
HasTypeOps.ReflectionName:=function(W,option)local res,t,i,total,l;
  t:=ReflectionType(W); res:="";total:=0;
  for i in [1..Length(t)] do
    Append(res,ReflectionName(t[i],option));
    if IsSpets(W) then
      if ForAll(Concatenation(List(t[i].orbit,T->T.indices)),j->IsBound(
        Group(W).reflectionsLabels[j])) then
        l:=Group(W).reflectionsLabels{Concatenation(List(t[i].orbit,
             T->T.indices))};
      else
        l:=Group(W).rootInclusion{Concatenation(List(t[i].orbit,T->T.indices))};
      fi;
    elif ForAll(t[i].indices,j->IsBound(W.reflectionsLabels[j])) then
      l:=W.reflectionsLabels{t[i].indices};
    else
      l:=W.rootInclusion{t[i].indices};
    fi;
    if l<>total+[1..Length(l)] then
 # added JM 3 oct. 2000 so that different reflection subgroups have
 # different reflection names: this is necessary in 'storefusion'
      if IsBound(option.TeX) then 
           PrintToString(res,"\\langle ",Join(l),"\\rangle ");
      else PrintToString(res,"<",Join(l),">");fi;
    fi;
    if i<Length(t) then 
      if IsBound(option.TeX) then Append(res,"\\times ");
      else Append(res,"x");fi;
    fi;
    total:=total+Length(l);
  od;
  if IsSpets(W) then 
    t:=Group(W).rank-Group(W).semisimpleRank ;
    if t>0 then
      i:=Copy(option);
      if not IsBound(i.Cyc) then i.expand:=true;fi;
      if res<>"" then Append(res,".");fi;
      l:=CycPol(Concatenation([1,0],List(
        Filtered(ReflectionDegrees(W),x->x[1]=1),x->AsRootOfUnity(x[2]))));
      l.vname:="q";
      Append(res,Format(l,i));
    fi;
  else
    t:=W.rank-W.semisimpleRank ;
    if t>0 then
      if res<>"" then Append(res,".");fi;
      Append(res,"(q-1)"); 
      if t>1 then 
        if t>=10 and IsBound(option.TeX) then PrintToString(res,"^{",t,"}");
        else PrintToString(res,"^",t);fi;
      fi;
    fi;
  fi;
  if res="" then res:=".";fi;
  return String(res);
end;

HasTypeOps.ChevieCharInfo:=function(W)local res,t,p,f,n,i,gt,keep;
  t:=ReflectionType(W);
  p:=List(t,ChevieCharInfo);
  keep:=Length(t)=1 and (not IsBound(t[1].orbit) or Length(t[1].orbit)=1);
  if keep then res:=ShallowCopy(p[1]); # keep extra fields when irreducible
  else res:=rec();
  fi; 
  res.charparams:=Cartesian(List(p,x->x.charparams));
  res.charnames:=List(res.charparams,x->CharName(W,x,rec()));
  if keep then return res;fi;
# if ForAll(p,x->IsBound(x.charnames)) then 
#   res.charnames:=List(Cartesian(List(p,x->x.charnames)),Join);
# fi;
  for f in ["positionId","positionDet"] do
    if ForAll(p,x->IsBound(x.(f))) then 
      res.(f):=PositionCartesian(List(p,x->Length(x.charparams)),
                                 List(p,x->x.(f)));
    fi;
  od;
  for f in ["b","B","a","A"] do
    if ForAll(p,x->IsBound(x.(f))) then 
      res.(f):=List(Cartesian(List(p,x->x.(f))),Sum);
    fi;
  od;
  if ForAny(p,x->IsBound(x.opdam)) then
    res.opdam:=List(p,function(x)if IsBound(x.opdam) then return x.opdam;
      else return ();fi;end);
    gt:=Cartesian(List(p,x->[1..Length(x.charparams)]));
    res.opdam:=PermListList(gt,
      List(gt,t->Zip(t,res.opdam,function(x,i)return x^i;end)));
  fi;
  if IsSpets(W) then 
    gt:=List(ReflectionType(Group(W)),x->Set(x.indices));
    n:=[];
    for i in [1..Length(t)] do for f in t[i].orbit do
      n[Position(gt,Set(f.indices{[1..f.rank]}))]:=p[i].nrGroupClasses;
    od;od;
    res.charRestrictions:=List(Cartesian(List(p,x->x.charRestrictions)),
      function(y)local m,i;m:=[];
	for i in [1..Length(t)] do for f in t[i].orbit do
	   m[Position(gt,Set(f.indices{[1..f.rank]}))]:=y[i];
	od;od;
        return PositionCartesian(n,m);
      end);
    res.nrGroupClasses:=Product([1..Length(t)],i->
      p[i].nrGroupClasses^Length(t[i].orbit));
  fi;
  return res;
end;

HasTypeOps.LowestPowerGenericDegrees:=function(W)local ci;
  ci:=ChevieCharInfo(W);
  if not IsBound(ci.a) then Error("no LowestPowerGenericDegrees for ",W);fi;
  return ci.a;
end;

HasTypeOps.HighestPowerGenericDegrees:=function(W)local ci;
  ci:=ChevieCharInfo(W);
  if not IsBound(ci.A) then Error("no HighestPowerGenericDegrees for ",W);fi;
  return ci.A;
end;

HasTypeOps.Invariants:=function(W)local V,i,N;
  V:=Parent(W);
  i:=List(ReflectionType(W),function(t)local H,ir,i;
    H:=ReflectionGroup(t);
    if CartanMat(V,W.rootInclusion{t.indices})<>CartanMat(H) then
      Error("not standard Cartan matrix: invariants not implemented");
    fi;
    ir:=PermRootOps.IndependentRoots(H);
    i:=Invariants(t);
    if i=false then return false;fi;
    return List(Invariants(t),f->function(arg)
#     Print("t=",t," H=",H,"\n");
      return ApplyFunc(f,H.simpleCoroots{ir}^-1*List(ir,i->
       PermRootOps.Coroot(V,W.rootInclusion[t.indices[i]]))*arg);end);
      end);
  if false in i then return false;fi;
  i:=Concatenation(i);
  if IsCoxeterGroup(W) then N:=W.simpleRoots;
  else N:=W.roots{PermRootOps.IndependentRoots(W)};
  fi;
  if N<>[] then
    N:=NullspaceMat(TransposedMat(N));
    Append(i,List(N,v->function(arg)return v*arg;end));
  fi;
  return i;
end;

HasTypeOps.Discriminant:=function(W)local t;
  t:=W.type;
  if Length(t)=0 then return function()return Mvp(1);end;
  elif Length(t)=1 then return Discriminant(t[1]);
  else Error("not implemented for non-irreducible");
  fi;
end;
  
# Fo an irreducible type, reps contain:
# .duflo,  .reps: elements of W represented as images of simple roots
# .character: decomposition of left cell in irreducibles
HasTypeOps.KLeftCellRepresentatives:=function(W)local n,res;
  n:=[];
  res:=List(ReflectionType(W),function(t)local R,rr,f;
    R:=ReflectionGroup(t);
    Add(n, NrConjugacyClasses(R));
    rr:=KLeftCellRepresentatives(t);
    if rr=false then return false;fi;
    return List(rr,function(r)local f;r:=ShallowCopy(r);
      f:=l->W.rootInclusion{t.indices{CoxeterWord(R,
	PermListList(R.roots,R.roots*R.roots{l}))}};
      r.duflo:=f(r.duflo);r.reps:=List(r.reps,f);Add(r.reps,r.duflo);
    return r;end);end);
  if Length(res)=0 then return;fi;
  if false in res then return false;fi;
  return List(Cartesian(res), function(l)local r;
    r:=rec(operations:=LeftCellOps,isDomain:=true,group:=W);
    r.duflo:=EltWord(W,Concatenation(List(l,x->x.duflo)));
    r.reps:=List(Cartesian(List(l,x->x.reps)),v->EltWord(W,Concatenation(v)));
    r.reps:=Difference(r.reps,[r.duflo]);
    r.character:=List(Cartesian(List(l,x->x.character)),
      p->PositionCartesian(n,p));
    r.a:=ChevieCharInfo(W).a{r.character};
    if Length(Set(r.a))>1 then Error();else r.a:=r.a[1];fi;
    return r;end);
end;
  
HasTypeOps.DecompositionMatrix:=function(W,p)local t;
  t:=List(ReflectionType(W),t->ReflTypeOps.DecompositionMatrix(t,p));
  if ForAny(t,x->x=false) then
    Error("DecompositionMatrix not implemented for ",W,"\n");
    return false;
  else return List(Cartesian(t),x->List(Cartesian(x),Product));
  fi;
end;

HasTypeOps.ParabolicRepresentatives:=function(W,s)local t,res,sols;
  if IsCoxeterCoset(W) then 
    return CoxeterCosetOps.ParabolicRepresentatives(W,s);fi;
  if IsCoxeterGroup(W) then return AbsCoxOps.ParabolicRepresentatives(W,s);fi;
  if IsSpets(W) then return SpetsOps.ParabolicRepresentatives(W,s);fi;
  t:=ReflectionType(W);
  sols:=Filtered(Cartesian(List(t,x->[0..SemisimpleRank(x)])),l->Sum(l)=s);
  return Concatenation(List(sols,c->List(Cartesian(List([1..Length(c)],
    function(i)local r,R;
     r:=ReflTypeOps.ParabolicRepresentatives(t[i],c[i]);
     if r=false then
       R:=ReflectionSubgroup(W,W.rootInclusion{t[i].indices});
       return PermRootOps.ParabolicRepresentatives(R,c[i]);
     elif ForAll(r,x->ForAll(x,y->y in [1..t[i].rank])) then
       return List(r,x->W.rootInclusion{t[i].indices{x}});
     else R:=ReflectionSubgroup(W,W.rootInclusion{t[i].indices});
       return List(r,x->R.rootInclusion{x});
     fi;end)),Concatenation)));
end;

WeightToAdjointFundamentalGroupElement:=function(W,i)local t,b,l;
  t:=First(ReflectionType(W),t->i in t.indices);
  l:=W.rootInclusion{t.indices};
  b:=LongestCoxeterElement(W,l)*LongestCoxeterElement(W,
     Difference(l,[W.rootInclusion[i]]));
  Add(l,W.rootInclusion[Maximum(Filtered([1..Length(W.roots)],
    i->ForAll([1..W.semisimpleRank],j->j in t.indices or W.roots[i][j]=0)))]);
  return RestrictedPerm(b,l);
end;

# returns a record containing minuscule coweights, decompositions
# (in terms of generators of the fundamental group)
HasTypeOps.WeightInfo:=function(W)local l,res,n;
  l:=List(ReflectionType(W),function(t)local r,g,C;
    r:=WeightInfo(t);
    g:=Filtered([1..Length(r.minusculeCoweights)],
       i->Sum(r.decompositions[i])=1); # generators of fundamental group
    r.ww:=List(t.indices{r.minusculeCoweights{g}},
      x->WeightToAdjointFundamentalGroupElement(W,x));
    C:=Mod1(CartanMat(t)^-1);
    r.csi:=NullMat(Length(g),SemisimpleRank(W));
    r.csi{[1..Length(g)]}{t.indices}:=C{r.minusculeCoweights{g}};
    r.minusculeWeights:=t.indices{r.minusculeWeights};
    r.minusculeCoweights:=t.indices{r.minusculeCoweights};
    return r;
    end);
  res:=rec(minusculeWeights:=Cartesian(List(l,
    x->Concatenation(x.minusculeWeights,[0]))),
    minusculeCoweights:=Cartesian(List(l,
      x->Concatenation(x.minusculeCoweights,[0]))),
    decompositions:=List(Cartesian(List(l,x->Concatenation(x.decompositions,
      [0*x.moduli]))),Concatenation),
    moduli:=Concatenation(List(l,x->x.moduli)));
# centre of simply connected group: the generating minuscule coweights
# mod the root lattice
  res.CenterSimplyConnected:=Concatenation(List(l,r->r.csi));
  res.AdjointFundamentalGroup:=Concatenation(List(l,r->r.ww));
  n:=Length(res.decompositions)-1;
  res.minusculeWeights:=List(res.minusculeWeights{[1..n]},
    x->Filtered(x,y->y<>0));
  res.minusculeCoweights:=List(res.minusculeCoweights{[1..n]},
    x->Filtered(x,y->y<>0));
  res.decompositions:=res.decompositions{[1..n]};
  return res;
end;
