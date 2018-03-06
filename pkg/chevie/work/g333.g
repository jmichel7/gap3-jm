g333sp:=function()local G,l,t;
 G:=ComplexReflectionGroup(33);
 l:=ReflectionSubgroup(G,[4,2,3]).rootInclusion;
 t:=TwistingElements(Spets(G),[4,2,3]);
 G:=ComplexReflectionGroup(3,3,3);
 return List(t,x->Spets(G,MatXPerm(G,PermListList(OnTuples(l,x),l))^-1));
end;

EigenvectorsMat:=function(m)local l,res,v;
  l:=EigenvaluesMat(m);
  res:=[];
  for v in List(Set(l),x->Filtered([1..Length(l)],i->l[i]=x)) do 
    res{v}:=NullspaceMat(m-l[v[1]]*m^0);od;
  return res;
end;

IsRegularVector:=function(W,v)
  return 0=Number(Reflections(W),s->v=v*MatXPerm(W,s));
end;

RegularEigenvalues:=function(W,m)
  return ListBlist(EigenvaluesMat(m),List(EigenvectorsMat(m),
    v->IsRegularVector(W,v)));
end;

z:=Mvp("z");
W:=ComplexReflectionGroup(3,3,3);

factors:=function(m)
  m:=List([x*y*z,x^3+y^3+z^3],p->
  OnPolynomials(m,p,["x","y","z"]));
  return List(m,p->List([x*y*z,x^3],n->Coefficient(p,n)));
end;

MvpOps.Coefficient:=function(p,m)local i;
  if Length(m.coeff)<>1 and m.coeff[1]<>1 then
    Error(m," should be a monomial");fi;
  i:=Position(p.elm,m.elm[1]);
  if i=false then return 0;
  else return p.coeff[i];
  fi;
end;

# reps of classes N_GL(W)/\pm 1

p:=[ (), ( 3,35,46)( 6,39,42)( 9,18,54)(12,31,45)(13,20,27)(14,32,15)(16,36,17)
    (19,41,26)(23,30,53)(28,51,29)(33,40,34)(37,43,38),
  ( 3,46,35)( 6,42,39)( 9,54,18)(12,45,31)(13,27,20)(14,15,32)(16,17,36)
    (19,26,41)(23,53,30)(28,29,51)(33,34,40)(37,38,43),
  ( 2,12,22,38)( 3,35,14,32)( 5,26, 8,33)( 6,39,13,20)( 7,42,24,27)
    ( 9,51,29,54)(10,15,21,46)(16,50,30,44)(17,53,23,36)(18,47,28,49)
    (19,41,40,34)(31,45,37,43) ];
N:=Group(p,());e:=Elements(N);

# note; p{[2,4]} already generates an N=A4.2 which 
# does not intersect W

# W has 24 orbits on triples of roots with same CartanMat as W.
# N acts simply transitively on these orbits
# W has 4 orbits of triples of reflections with same braidrels as W.
# N acts via its quotient A4 on these.

# all distinct reflections of W
refs:=Set(List(Reflections(W),x->Position(Reflections(W),x)));

triples:=List(Combinations(refs,3),s->Reflections(W){s});
triples:=Filtered(triples,x->Size(Group(x,()))=54);
orbits:=Orbits(W,List(triples,Set),OnSets);

# for a phi of order 3 find a Ruth-quad which stays in same W-orbit and
# return an element of W.phi which stabilizes quadruple
test:=function(x)local e,reps,sets,i,scalar,v,m;
  # quads of roots which have the same CartanMat and give representatives of 
  # W-orbits of quads of reflections satisfying Ruth relations
  sets:=[[1,2,3,44],[21,3,1,32],[3,11,2,36],[22,3,2,16]];
  reps:=List(sets,x->Reflections(W){x});
  Print("x^",OrderPerm(x),"=1");
  for i in [1..4] do
    e:=RepresentativeOperation(W,Set(sets[i]),OnSets(Set(sets[i]),x),OnSets);
    if e<>false then 
      x:=x*e^-1;
      scalar:=AsRootOfUnity(ProportionalityCoefficient(W.roots[sets[i][3]^x],
        W.roots[sets[i][3]]));
      v:=scalar+[0..2]/3;
      m:=Minimum(List(v,Denominator));
      m:=PositionProperty(v,n->Denominator(n)=m);
      scalar:=E(Denominator(v[m]))^Numerator(v[m]);
      z:=PermMatX(W,E(3)*IdentityMat(3));
      x:=z^(m-1)*x;
      scalar:=ProportionalityCoefficient(W.roots[sets[i][3]^x],
        W.roots[sets[i][3]]);
      Print(" i=",i," scalar=",scalar," induces ",
       PermListList(W.roots{OnTuples(sets[i],x)},W.roots{sets[i]}*scalar),"\n");
      return x;
    fi;
  od;
  Print(" false\n");
  return false;
end;

# for a phi of order 4 find a quadruple which stays in same W-orbit and
# return an element of W.phi which stabilizes quadruple
test4:=function(x)local op,reps,sets,i,scalar,v,m,res;
  # quads of roots which have the same CartanMat and give representatives of 
  # W-orbits of appropriate quads of reflections 
  sets:=[[1,50,3,12],[3,52,2,23],[1,16,3,43],[2,37,3,15],[50,3,52,38],[1,23,3,22]];
       # [[1,2,3,12],[3,1,2,20],[1,15,3,9],[2,6,3,15],[2,3,1,15],[1,20,3,5]];
  Print("x^",OrderPerm(x),"=1\n");
  res:=[];
  for i in sets do
    v:=W.rootInclusion{i};
    y:=Set(v);
    op:=RepresentativeOperation(W,y,OnSets(y,x),OnSets);
    if op<>false then
      Print("i=",Position(sets,i)," op=",
       PermListList(v,OnTuples(v,x*op^-1)),"\n");fi;
    if op<>false and PermListList(v,OnTuples(v,x*op^-1))=(1,2,3,4) then
      Print("gone here\n");
      x:=x*op^-1;
      scalar:=List([1..4],j->ProportionalityCoefficient(W.roots[i[j]^x],
        W.roots[i[1+(j mod 4)]]));
      if false in scalar then Print("***rev ");
      scalar:=List([1..4],j->ProportionalityCoefficient(W.roots[i[j]^x],
        W.roots[i[1+((j-2) mod 4)]]));
      fi;
      if Length(Set(scalar))>1 then Error();fi;
      scalar:=AsRootOfUnity(scalar[1]);
      v:=scalar+[0..2]/3;
      m:=Minimum(List(v,Denominator));
      m:=PositionProperty(v,n->Denominator(n)=m);
      scalar:=E(Denominator(v[m]))^Numerator(v[m]);
      z:=PermMatX(W,E(3)*IdentityMat(3));
      x:=z^(m-1)*x;
      scalar:=List([1..4],j->ProportionalityCoefficient(W.roots[i[j]^x],
        W.roots[i[1+(j mod 4)]]));
      if false in scalar then Print("***rev ");
      scalar:=List([1..4],j->ProportionalityCoefficient(W.roots[i[j]^x],
        W.roots[i[1+((j-2) mod 4)]]));
      fi;
      Print(" i=",i," scalar=",scalar[1], " induces ",
       PermListList(W.roots{OnTuples(i,x)},W.roots{i}*scalar[1]),"\n");
      Add(res,x);
    fi;
  od;
  if Length(res)>0 then return res[1];fi;
  Print(" false\n");
  return false;
end;

# find all tuples of roots of W with CartanMat c
findrootscartan:=function(W,c)local l,i;
  l:=[[]];
  for i in [1..Length(c)] do
    l:=Concatenation(List(l,x->List([1..Length(W.roots)],y->Concatenation(x,[y]))));
    l:=Filtered(l,x->CartanMat(W,x)=c{[1..i]}{[1..i]});
    Print(i,"-tuples:",Length(l),"\n");
  od;
  return l;
end;

sets:=[[1,50,3,12],[3,52,2,23],[1,16,3,43],[2,37,3,15],[50,3,52,38],[1,23,3,22]];
rot:=function(o,x)
  return Filtered([1..Length(o)], function(i)local op,y;
    y:=o[i][1];
    op:=RepresentativeOperation(W,Set(y),OnSets(Set(y),x),OnSets);
    if op=false then return false;
    else 
      Print("i=",i,"=>",PermListList(y,OnTuples(y,x*op^-1)),"\n");
      return PermListList(y,OnTuples(y,x*op^-1))=(1,2,3,4);
    fi;
    end);
end;
#act1:=function(t,x)
# return Set(OnSets(
