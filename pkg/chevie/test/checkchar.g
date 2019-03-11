# test if the powermap of CharTable(W) is OK
CHEVIE.AddTest("PowerMap",
function(W)local cl,ct,p,i,cmap;
  cl:=List(ConjugacyClasses(W),Representative);ct:=CharTable(W);p:=ct.powermap;
  for i in [1..Length(p)] do
   if IsBound(p[i]) then cmap:=List(cl,x->PositionClass(W,x^i));
     if cmap<>p[i] then
       ChevieErr(i,"-th power map is ",p[i],"\n","should be:",cmap,"\n");
     fi;
   fi;
  od;
end,
x->not IsSpets(x));

# test ReflectionEigenvalues, ChevieCharInfo.extRefl, .positionDet, .positionId
CHEVIE.AddTest("ExtRefl",
function(W)local ct,ci,extRefl,v,n,ref,checkfield;
  ct:=CharTable(W);
  # compute first using ReflectionEigenvalues
  n:=NrConjugacyClasses(W);
  v:=List([1..n],i->Product(ReflectionEigenvalues(W,i),n->X(Cyclotomics)+
    E(Denominator(n))^Numerator(n)));
  v:=List(v*X(Cyclotomics)^0,x->x.coefficients);
  extRefl:=Reversed(List(TransposedMat(v),x->Position(ct.irreducibles,x)));
  # check ReflectionEigenvalues using AntiSymmetricParts
  ref:=Position(ct.irreducibles,List(ChevieClassInfo(W).classtext,
		  w->ReflectionCharValue(W,EltWord (W,w))));
  if ref<>false then #irreducible group
    ref:=List([1..W.semisimpleRank],
                i->AntiSymmetricParts(ct,ct.irreducibles{[ref]},i)[1]);
    ref:=List(ref,x->Position(ct.irreducibles,x));
    if ref<>extRefl{[2..Length(extRefl)]} then 
      ChevieErr("ReflectionEigenvalues disagrees with AntiSymmetricParts");fi;
  fi;
  ci:=ChevieCharInfo(W);
  checkfield:=function(f,v)
    if not IsBound(ci.(f)) then
      ChevieErr(".",f," not bound but should be ",v,"\n");
    elif v<>ci.(f) then 
      ChevieErr(".",f," is ",ci.(f)," should be ",v,"\n");
    fi;
  end;
  if ref<>false then checkfield("extRefl",extRefl);fi;
  checkfield("positionDet",extRefl[Length(extRefl)]);
  checkfield("positionId",extRefl[1]);
end,
x->not IsSpets(x));

##  The function below uses the formula
##  prod_{g\in W}det(g-t)=prod_i(t^{d_i}-1)^{|W|/d_i}
##  Superseded by HasTypeOps.ReflectionDegrees but can recompute data
##
PermRootOps.ReflectionDegrees:=function(W)local l,c,res,i,e,p,d,j,eig,mul,n;
  c:=CharTable(W).classes;
  l:=Concatenation(List([1..Length(c)],
    i->List(ReflectionEigenvalues(W,i),e->[e,c[i]])));
  e:=CollectBy(l,x->x[1]);e:=List(e,x->[x[1][1],Sum(x,y->y[2])]);
  e:=TransposedMat(e);
  if e=[] then return [];fi;
  eig:=e[1];mul:=e[2];
  res:=[];
  for d in Reversed(Set(List(eig,Denominator))) do
   p:=Position(eig,Mod1(1/d));
   if p<>false and mul[p]>0 then
     for i in [1..d*mul[p]/Size(W)] do Add(res,d);od;
     n:=mul[p];
     for j in [0..d-1] do p:=Position(eig,j/d); mul[p]:=mul[p]-n;od;
   fi;
  od;
  W.degrees:=Reversed(res);return W.degrees;
end;

##   the function below uses the formula
##   prod_{g\in W}det(g\phi-t)=prod_i(t^{d_i}\zeta_i-1)^{|W|/d_i}
##  Superseded by HasTypeOps.ReflectionDegrees but can recompute data
##
SpetsOps.ReflectionDegrees:=function(W)local l,res,i,e,p,d,mul,searchdeg;
  l:=Zip(ReflectionEigenvalues(W),List(ConjugacyClasses(W),Size)/Size(W),
       function(v,t)return List(v,x->[x,t]);end);
  l:=CollectBy(Concatenation(l),x->x[1]);
  e:=List(l,x->[x[1][1],Sum(x,y->y[2])]);
  # here we got LHS of formula as \prod(t-eig[i])^mul[i]
  # more exactly a list of pairs [eig in Q/Z, mul/|W|]
  searchdeg:=function(e,card,degs)local res,d,f,g,p,pos,ne;
#   Print("degs=",degs," card=",card,"\n");
    if Length(degs)=0 then return [[]];fi;
    e:=Filtered(e,x->x[2]<>0);
    res:=[];
    d:=degs[1];
    f:=Filtered([1..Length(e)],i->e[i][2]>=1/d);
    g:=List(Filtered(e{f},x->x[1]<1/d),x->x[1]);
    for p in g do
      pos:=List([0..d-1]/d,i->PositionProperty(f,j->e[j][1]=p+i));
      if not false in pos then
#       Print("d=",d," p=",p,"\n");
        ne:=Copy(e);
        ne{f{pos}}[2]:=ne{f{pos}}[2]-1/d;
        for ne in searchdeg(ne,card/d,degs{[2..Length(degs)]}) do
          Add(res,Concatenation([[d,Mod1(d*p)]],ne));
        od;
      fi;
    od;
    return res;
  end;
  mul:=Set(searchdeg(e,Size(W),ReflectionDegrees(Group(W))));
  e:=List(mul,x->Permuted(x,SortingPerm(x)));
  e:=Set(e);
  if Length(e)>1 then Error();fi;
# Print("mul=",mul[1],"\n");
  return List(e[1],x->[x[1],E(Denominator(x[2]))^Numerator(x[2])]);
end;

# check ReflectionDegrees
CHEVIE.AddTest("ReflectionDegrees",function(W)local d,d1;
  d:=ReflectionDegrees(W);Sort(d);
  if IsSpets(W) then d1:=SpetsOps.ReflectionDegrees(W);
  else d1:=PermRootOps.ReflectionDegrees(W);
  fi;
  Sort(d1);
  CHEVIE.Check.EqLists(d,d1,"ReflectionDegrees","computed ReflectionDegrees");
  if IsCoxeterGroup(W) and ForAll(Concatenation(W.cartan),IsInt)  then 
    d1:=List(Collected(List(W.roots{[1..W.N]},Sum)),x->x[2]);
    d1:=1+AssociatedPartition(d1);
    d1:=Concatenation([W.semisimpleRank+1..W.rank]*0+1,d1);
    Sort(d1);
    CHEVIE.Check.EqLists(d,d1,"ReflectionDegrees",
      "dual partition of root heights");
  fi;
end);

# check FakeDegrees, ChevieCharInfo.b, .B
CHEVIE.AddTest("FakeDegrees",function(W)local fd,ffd,q;
  q:=X(Cyclotomics);
  fd:=FakeDegrees(W,q);
  if IsSpets(W) then ffd:=SpetsOps.FakeDegrees(W,q);
  else ffd:=PermRootOps.FakeDegrees(W,q);
  fi;
  CHEVIE.Check.EqLists(ffd,fd,"FakeDegrees","computed FakeDegrees");
  CHEVIE.Check.EqLists(List(fd,Valuation),LowestPowerFakeDegrees(W),
      "ChevieCharInfo.b","computed b");
  CHEVIE.Check.EqLists(List(fd,Degree),HighestPowerFakeDegrees(W),
      "ChevieCharInfo.B","computed B");
end);

# check invariants
CHEVIE.AddTest("Invariants",
function(W)local ii,i,j,vars;
  ii:=Invariants(W);
  vars:=List([1..W.rank],i->Mvp(SPrint("x",i)));
  ii:=List(ii,x->ApplyFunc(x,vars));
  InfoChevie(" #");
  for i in [1..Length(ii)] do 
    for j in [1..Length(W.matgens)] do 
      InfoChevie("W.",j,"*I",i,",\c");
      if not OnPolynomials(W.matgens[i],ii[i],vars)=ii[i] then 
         ChevieErr("not invariant\n");
      fi;
    od;
  od;
end,
x->not IsSpets(x) and Size(x)<14400); # H4 first painful client

# check whether CharTable(W) agrees with PermGroupOps.CharTable
CHEVIE.AddTest("CharTable",function(W)local ct,W1,ct1,f1,f,p,i,cl;
  ct:=CharTable(W);
  Unbind(W.charTable);
  cl:=List(ConjugacyClasses(W),Representative);
  if IsSpets(W) then ct1:=SpetsOps.CharTable(W);
    cl:=List(ct.irreducibles,x->PositionProperty(ct1.irreducibles,
      y->ProportionalityCoefficient(y,x)<>false));
    if ForAll(cl,x->x<>false) then 
      f:=List([1..Length(cl)],i->ProportionalityCoefficient(
	ct1.irreducibles[cl[i]],ct.irreducibles[i]));
      if ForAny(f,x->x<>1) then
	#ChevieErr("coefficients:",f,"\n");
	for i in [1..Length(cl)] do 
	  ct1.irreducibles[cl[i]]:=ct1.irreducibles[cl[i]]/f[i];
	od;
      fi;
    fi;
  for f in ["cartan","classparams","classparam","classtext","sqrtparameter",
    "identifier","name","parameter","ST","degrees","dim","indexchar",
    "indexclasses","reflclasses","classnames"] do Unbind(ct1.(f));od;
  else
    W1:=Group(W.generators,());W1.name:=W.name;
    cl:=List(cl,x->PositionClass(W1,x));
    if Length(Set(cl))<>Length(cl) then Error("classes repeated");fi;
    ct1:=CharTableSplitClasses(CharTable(W1),cl);
    cl:=List(ct.irreducibles,x->Position(ct1.irreducibles,x));
  fi;
  if ForAny(cl,x->x=false) or Length(Set(cl))<>Length(cl) then
    ChevieErr("irreducibles");
  else  ct1.irreducibles:=ct1.irreducibles{cl};
    if IsBound(ct1.irredinfo) then ct1.irredinfo:=ct1.irredinfo{cl};fi;
  fi;
  Unbind(ct1.identifier); Unbind(ct1.name);
  ct:=Copy(ct);
  for f in ["cartan","classparams","classparam","classtext","sqrtparameter",
    "identifier","name","parameter","ST","degrees","dim","indexchar",
    "indexclasses","reflclasses"] do Unbind(ct.(f));od;
  CHEVIE.Check.EqObj(ct,ct1,rec(na:="CharTable",nb:="GroupOps.CharTable"));
end);

# checks thet ChevieClassInfo.classes is correct
CHEVIE.AddTest("ClassSizes",function(W)
  CHEVIE.Check.EqLists(List(ConjugacyClasses(W),Size),
   ChevieClassInfo(W).classes,"sizes of classes",
                              "ChevieClassInfo.classes");
end);

# FakeDegreesInduce(W[,J])
# check that Fake degrees induce from a Levi as InductionTable
CHEVIE.AddTest("FakeDegreesInduce",function(arg)
  local W,J,L,t,hd,ud,index,q,pred,found;
  W:=arg[1];
  if Length(arg)=1 then 
    for J in ParabolicRepresentatives(W) do 
      CHEVIE.Test("FakeDegreesInduce",W,J);
    od;
    return;
  fi;
  J:=arg[2];
  q:=X(Cyclotomics);
  if not IsSpets(W) then W:=Spets(W);fi;
  for L in Twistings(W,J) do
    t:=InductionTable(L,W);
    hd:=FakeDegrees(L,q);
    ud:=FakeDegrees(W,q);
    index:=GenericSign(L)*GenericOrder(W,q)/(GenericOrder(L,q)*GenericSign(W));
    index.valuation:=0;
    pred:=hd*index;found:=ud*t.scalar;
    CHEVIE.Testing(" from",ReflectionName(L));
    InfoChevie("\n   # R^{",ReflectionName(W),"}_{",ReflectionName(L),"}");
    if pred<>found then ChevieErr("quotient ",
      List([1..Length(hd)],i->CycPol(pred[i])/CycPol(found[i])),"\n");
    fi;
  od;
end);

#  2 methods to check CharTable(Hecke(3D4,q))
CHEVIE.AddTest("CharsHecke3D4",
function(3D4)local F4,WF,T,v,m,3D4bis,q,ct,D4;
  F4:=CoxeterGroup("F",4);
  q:=Mvp("q");
  T:=Basis(Hecke(F4,[q,q,1,1]),"T");
  v:=[T(2),T(3,2,3),T(1),T(4,3,2,3,4)]; # embedding of Hecke(3D4,q)
  m:=TransposedMat(List(ChevieClassInfo(3D4).classtext,
              x->HeckeCharValues(Product(v{x})*T(4,3))));
  WF:=Spets(F4);
  3D4bis:=CoxeterSubCoset(WF,[2,9,1,16],EltWord(F4,[4,3])); # inside F4
  ct:=CharTable(Hecke(3D4,q)).irreducibles;
  CHEVIE.Check.EqLists(ct,
   m{List(TransposedMat(InductionTable(3D4bis,WF).scalar),x->Position(x,1))},
   "charTable(Hecke(3D4))","computed by coset induction");
  D4:=ReflectionSubgroup(F4,[2,9,1,16]);
  v:=TransposedMat(InductionTable(D4,F4).scalar);
  v:=v{ChevieCharInfo(3D4).charRestrictions};
  CHEVIE.Check.EqLists(ct,-m{List(v,x->Position(x,2))},
   "charTable(Hecke(3D4))","computed by subgroup induction");
end, 
W->ReflectionName(W)="3D4");

# Check ChevieClassInfo(W).indexclasses for W in G_4--G_{22}
CHEVIE.AddTest("G4_22indexClasses",
function(W)local t,O,a,b,c,e,l;
  t:=ReflectionType(W)[1];
  O:=ComplexReflectionGroup(CHEVIE.Data("Generic",t));
  e:=CHEVIE.Data("Embed",t);
  c:=List(ChevieClassInfo(W).classtext,c->Concatenation(List(c,x->e[x])));
  c:=List(c,x->PositionClass(O,EltWord(O,x)));
  l:=List(e,x->Position(Reflections(O),EltWord(O,x)));
  a:=ChevieClassInfo(W).indexclasses;
  if not false in l then
    W:=ReflectionSubgroup(O,l);
    b:=FusionConjugacyClasses(W,O);
  else b:=c;
  fi;
  if a<>b or b<>c then 
    a:=TransposedMat(Filtered(TransposedMat([a,b,c]),x->Length(Set(x))>1));
    ChevieErr("\n", FormatTable(a{[2,3]},
      rec(rowLabels:=[SPrint("fusion to ",ReflectionName(O)),"classtext"],
        columnLabels:=a[1],rowsLabel:="indexclasses",
        screenColumns:=SizeScreen()[1])));
  fi;
end,
function(W)local t;t:=ReflectionType(W);
  return Length(t)=1 and IsBound(t[1].ST) and t[1].ST in [4..22];end);

CHEVIE.AddTest("Opdam",
function(W)local ct,complex,x,fd,ci,c,i,f,opdam;
  ct:=CharTable(W).irreducibles;
  complex:=PermListList(ct,ComplexConjugate(ct));
  x:=Mvp("x");
  fd:=FakeDegrees(W,x);
  ci:=ChevieCharInfo(W);
  if IsBound(ci.opdam) then opdam:=ci.opdam;else opdam:=();fi;
  # c_\chi in Gordon-Griffeth
  c:=function(W,i)local cl;
    cl:=ChevieClassInfo(W).classes;
    return Sum(ReflectionDegrees(W)-1)-Sum(HyperplaneOrbits(W),
      h->Sum(h.classno,c->cl[c]*ct[i][c]/ct[i][1]));
  end;
  for i in [1..NrConjugacyClasses(W)] do
    f:=fd[i^complex];
    if fd[i^opdam]<>x^c(W,i)*Value(f,["x",x^-1]) then 
      ChevieErr("Opdam wrong for ",i,"\n");
    fi;
  od;
end,
x->not IsSpets(x));
