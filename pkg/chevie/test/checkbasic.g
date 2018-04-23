CHEVIE.AddTest("BraidRelations",function(arg)local W,r,gens,rel;
  W:=arg[1];
  if Length(arg)=1 then r:=BraidRelations(W);
  else r:=BraidRelations(arg[2]);
  fi;
  gens:=List(W.generatingReflections,i->Reflection(W,i));
  for rel in r do
    rel:=List(rel,x->List(W.rootRestriction{x},
      y->Position(W.generatingReflections,y)));
    CheckRelation(gens,rel,ChevieErr);
  od;
end,
W->not IsSpets(W),
"(W[,t]) check that W satisfies braid relations of its type [of type t]"
);

CHEVIE.AddTest("RootSystem",
function(W)local index,m,R,cr,rb,integrality,Coroots,try,p,replines,x;
# Check that all elements of the list of vectors l have coefficients
#  in the ring R when expressed as linear combinations of l{indices}
  integrality:=function(l,indices,R)local rb;
    rb:=List(l,r->SolutionMat(l{indices},r));
    return ForAll(Flat(rb),y->y in R);
  end;
  Coroots:=W->List([1..Length(W.roots)],i->PermRootOps.Coroot(W,i));
# Find representatives up to scalar multiple of elements of the list of
# vectors vec. Check that other elements differ from such a representative
# by a unit of the ring R
  replines:=function(vec,R)local found,res,v,p;res:=[];
    for v in vec do found:=false;
      for x in res do p:=ProportionalityCoefficient(x,v);
        if p<>false then 
	  if not (p in R and 1/p in R) then return false;fi;
	  found:=true;
        fi;
      od;
      if not found then Add(res,v);fi;
    od;
    return res;
  end;
  R:=DefaultRing(Flat(CartanMat(W)));
  try:=Combinations(W.generatingReflections,W.rank);
  p:=PositionProperty(try,v->integrality(W.roots,v,R));
  if p=false then ChevieErr("found no root basis");else;fi;
  cr:=Coroots(W);
  p:=PositionProperty(try,v->integrality(cr,v,R));
  if p=false then ChevieErr("found no coroot basis");else;fi;
  p:=replines(W.roots,R);
  if p=false then ChevieErr("not reduced");fi;
  m:=replines(cr,R);
  if m=false then ChevieErr("coroots not reduced??");fi;
  if Length(m)>0 and not ForAll(Set(Flat(p*TransposedMat(m))),x->x in R) then 
    ChevieErr("not a root system");fi;
end,
W->not IsSpets(W),
"Check that the roots of W define a distinguished root system in the sense\
 of Broue-Corran-Michel");

CHEVIE.AddTest("GaloisAutomorphisms",
function(W)local k,Wk,m,g,gm,p;
  k:=Field(Flat(W.matgens));
  Wk:=Field(Flat(CartanMat(W)));
  if k<>Wk then 
    ChevieErr("k_W=",Wk," but matrices over ",k,"\n");
  fi;
  if k=Rationals then return;fi;
  for g in GaloisGroup(k).generators do
    for m in W.matgens do
      gm:=List(m,x->OnTuples(x,g)); p:=PermMatX(W,gm);
      if not p in W or gm<>MatXPerm(W,p) then 
        ChevieErr("not Galois stable\n");fi;
    od;
  od;
end,
W->not IsSpets(W),
"check that W's reflection representation is globally invariant by Gal(k_W/Q)");

CHEVIE.AddTest("WFromBraidRelations",
function(W)local n,F,r;
  n:=Length(W.generators);
  F:=FreeGroup(n);
  r:=List(BraidRelations(W),x->List(x,y->Product(y,z->F.(z))));
  r:=List(r,x->x[1]/x[2]);
  Append(r,List([1..n],i->F.(i)^OrderPerm(W.(i))));
  return Size(W)=Size(F/r);
end,
W->not IsSpets(W) and Size(W)<64000,
"check that the abstract group defined by the braid and order relations has\
 the expected Size");

CHEVIE.AddTest("ClassRepresentatives",
function(W)local cl,i,w,l,o,wF;
  cl:=ChevieClassInfo(W);
  for i in [1..NrConjugacyClasses(W)] do
    if IsCoxeterGroup(W) then
      w:=Braid(W)(cl.classtext[i]);
      o:=cl.orders[i]; l:=BrieskornNormalForm(w^o);
    elif IsCoxeterCoset(W) then
      w:=Braid(CoxeterGroup(W))(cl.classtext[i]);
      wF:=EltWord(W,cl.classtext[i]);
      o:=PositionProperty([1..OrderPerm(wF)],j->(wF)^j=W.phi^j);
      l:=BrieskornNormalForm(TwistedPower(o,w,Frobenius(W)));
    fi;
    if Length(l) mod 2<>0 or ForAny([1..Length(l)/2],j->l[2*j-1]<>l[2*j])
      or ForAny([1..Length(l)-1],j->not IsSubset(l[j],l[j+1])) then
      ChevieErr("class ",i," is not good\n");
    fi;
    if o mod 2=0 then 
      if IsCoxeterGroup(W) then l:=BrieskornNormalForm(w^(o/2));
      elif IsCoxeterCoset(W) then
	l:=BrieskornNormalForm(TwistedPower(o,w,Frobenius(W)));
      fi;
      if ForAny([1..Length(l)-1],j->not IsSubset(l[j],l[j+1])) then
	ChevieErr("class ",i," is not very good\n"); 
      fi;
    fi;
  od;
end,
W->IsCoxeterGroup(W) or IsCoxeterCoset(W),
"check that the class representatives are very good in the sense of\
 Geck-Michel");
