#############################################################################
##
#A algebra.gap              ALGEBRA package             C'edric Bonnaf'e 2005
##
#Y  Copyright (C) 2005 -   C'edric Bonnaf'e
##
##  This file contains functions to work with finite-dimensional algebras
##  
###   Code and documentation reviewed by Jean Michel 4/2010
#############################################################################

###################################
## Utility functions              #
###################################

# Digits(n[,base])  digits of n in base <base> (default 10)
Digits:=function(arg) local n,base,d;
  n:=arg[1];
  if Length(arg)=1 then base:=10; else base:=arg[2]; fi;
  d:=[];
  while n<>0 do Add(d,n mod base); n:=QuoInt(n,base); od;
  return Reversed(d);
end;

# ByDigits(l[,base]) integer which has digits l in base <base> (default 10)
ByDigits:=function(arg)local d,res,base;
  if Length(arg)=1 then base:=10; else base:=arg[2]; fi;
  res:=0;
  for d in arg[1] do res:=base*res+d;od;
  return res;
end;

# Part of n product of factors of p (of elements of p if p is a list)
PiPart:=function(n,p)local l;
  if IsInt(p) then 
    if p=0 then return n;fi;
    if IsPrime(p) then l:=1;
      while n mod p=0 do n:=n/p;l:=l*p;od;
      return l;
    else p:=Set(FactorsInt(p));
    fi;
  fi;
  return Product(p,i->PiPart(n,i));
end;

# SignedCompositions(n) . . . . .  the set of signed compositions of n
# that   is,   tuples   of   non-zero   integers   [i_1,...i_r]  such  that
# |i_1|+...+|i_r|=n.
SignedCompositions:=n->Set(Concatenation(List(OrderedPartitions(n),
    l->Cartesian(List(l,i->[i,-i])))));

# SignedPartitions(n) . . . . .  the set of signed partitions of n, that  is,
# tuples  of  integers  [i_1,...,i_r,j_1,...,j_s]$  such that i_k>0, j_k<0,
# |i_1|+...+|i_r|+|j_1|+...+|j_s|=n, i_1>=...>=i_r and |j_1|>=...>=|j_s|.
SignedPartitions:=n->Set(Concatenation(List([0..n],i->
    List(Cartesian(Partitions(i),-Partitions(n-i)),Concatenation))));

# reduction mod p of the p-integral cyclotomic x
CyclotomicModP:=function(x,p) local r,n,zeta,np,pp;
  if IsList(x) then return List(x,y->CyclotomicModP(y,p));fi;
  if IsInt(x) then return x*Z(p)^0;fi;
  x:=COEFFSCYC(x); n:=Length(x); # n:=NofCyc(x);x:=CoeffsCyc(x,n); 
  pp:=PiPart(n,p);np:=n/pp;r:=OrderMod(p,np);
  zeta:=Z(p^r)^((p^r-1)/np*GcdRepresentation(pp,p^r-1)[1]); #n-th root of unity
  if zeta^n<>Z(p)^0 then Error();fi;
  return Sum([1..n],i->zeta^(i-1)*Numerator(x[i])/Denominator(x[i]));
  #Stupidly, product of finite field element and rational is not defined
end;

# p-blocks of finite group G
PBlocks:=function(G,p)local T,l;
  if not IsBound(G.pBlocks) then G.pBlocks:=rec();fi;
  if not IsBound(G.pBlocks.(p)) then
    T:=CharTable(G); l:=Length(T.irreducibles);
    T:=List(T.irreducibles,chi->List([1..l],
	j->CyclotomicModP(T.classes[j]*chi[j]/chi[1],p)));
    G.pBlocks.(p):=Set(CollectBy([1..l],T));
  fi;
  return G.pBlocks.(p);
end;

# p-part of element g of G
PiComponent:=function(G,g,p)local n,np;
  n:=Order(G,g); np:=PiPart(n,p);
  return g^(GcdRepresentation(np,n/np)[2]*n/np);
end;

# collect indices of classes of same pi-component
PiSections:=function(G,p)local pi;
  pi:=List(ConjugacyClasses(G),c->PiComponent(G,Representative(c),p));
  return CollectBy([1..Length(pi)],List(pi,x->PositionClass(G,x)));
end;

# collect indices of classes of same pi'-component
PiPrimeSections:=function(G,p)local pi;
  pi:=List(List(ConjugacyClasses(G),Representative),x->x/PiComponent(G,x,p));
  return CollectBy([1..Length(pi)],List(pi,x->PositionClass(G,x)));
end;

PRank:=function(P,p)local res,pval;
  pval:=n->Number(FactorsInt(n),x->x=p);
  P:=SylowSubgroup(P,p);
  if IsElementaryAbelian(P) then return pval(Size(P));fi;
  #elif IsAbelian(P) then res:=Size(P)/Size(FrattiniSubgroup(P));
  #  return Length(FactorsInt(res));
  res:=List(ConjugacyClassesSubgroups(P),Representative);
  res:=Filtered(res,IsElementaryAbelian);
  return pval(Maximum(List(res,Size)));
end;

################################################################
## TablePrint is a variation on CharTable printing suitable   ##
## for our algebras                                           ##
################################################################

TablePrint:=function(table)local i,j,clen,cl,rl,nbcols,cwidths,t,x;
  Print("\n");
  clen:=List(table.columns,Length); cl:=Maximum(clen); nbcols:=Length(clen);
  rl:=Maximum(List(table.rows,Length));
  t:=List(table.matrix,r->List(r,function(x)
    if x=0*x then return ".";else return Format(x);fi;end));
  cwidths:=[];
  for i in [1..nbcols] do 
    cwidths[i]:=Maximum(List([1..Length(table.rows)],j->Length(t[j][i])));
    if IsList(table.columns[i]) and not IsString(table.columns[i]) and
       Length(table.columns[i])>0
    then
      cwidths[i]:=Maximum(cwidths[i],Maximum(List(table.columns[i],Length)));
    fi;
  od;
  for j in [1..cl] do 
    Print(String("",rl+1));
    for i in [1..nbcols] do 
      if j<=cl-clen[i] then x:="";else x:=table.columns[i][j-cl+clen[i]];fi;
      Print(String(x,2+cwidths[i]));
    od;
    Print("\n");
  od;
  Print("\n");
  for i in [1..Length(table.rows)] do 
    Print(String(table.rows[i],rl)," ");
    for j in [1..nbcols] do Print(String(t[i][j],2+cwidths[j]));od;
    Print("  \n");
  od;
end;

####################################
## Operations on algebra elements ##
####################################

AlgebraEltOps:=OperationsRecord("AlgebraEltOps");

AlgebraElement:=function(arg)local r;
  r:=rec(operations:=AlgebraEltOps,domain:=arg[1],coefficients:=arg[2]);
  if Length(arg)=3 and arg[3]<>arg[1] then r.info:=arg[3];fi;
  return r;
end;

IsAlgebraElement:=r->IsRec(r) and IsBound(r.operations) and
  r.operations=AlgebraEltOps;

AlgebraEltOps.String:=function(r)local res,i,A,f,s,info;
  A:=r.domain; 
  if IsBound(r.info) then info:=r.info;else info:=A;fi;
  if IsBound(info.string) then return info.string(r);fi;
  res:=r.coefficients;
  if Length(res)=0 then return "0";fi;
  f:=function(coef) local res;
    if coef[1]=-A.field.one and coef[1]<>A.field.one then res:="-";
    elif coef[1]<>A.field.one then res:=SPrint(coef[1],"*");
    else res:="";
    fi;
    return SPrint(res,info.basisname,"(",info.parameters[coef[2]],")");
  end;
  s:=f(res[1]);
  for i in [2..Length(res)] do 
    if res[i][1]<0 or (res[i][1]=-A.field.one and A.field.char<>2) 
    then PrintToString(s,"-",f([-res[i][1],res[i][2]]));
    else PrintToString(s,"+",f(res[i]));
    fi;
  od;
  return s;
end;

AlgebraEltOps.Print:=function(r)Print(String(r));end;

AlgebraEltOps.\=:=function(x,y)
  if not IsRec(x) then return false;fi;
  if not IsRec(y) then return false;fi;
  if x.domain.identification<>y.domain.identification then return false;fi;
  if x.coefficients<>y.coefficients then return false;fi;
  if IsBound(x.info) then 
   if not IsBound(y.info) or x.info<>y.info then return false;fi;
  fi;
  return true;
end;

AlgebraEltOps.Normalize:=function(x)
  x.coefficients:=List(CollectBy(x.coefficients,x->x[2]),
    m->[Sum(m,x->x[1]),m[1][2]]);
  x.coefficients:=Filtered(x.coefficients,x->x[1]<>0*x[1]);
  return x;
end;

AlgebraEltOps.\+:=function(x,y)local sum;
  if x.domain.identification <> y.domain.identification then 
    Error("Elements must lie in the same algebra");
  fi;
  if IsBound(x.info) and (not IsBound(y.info) or x.info<>y.info) then
    Error("Elements must belong to the same basis");
  fi;
  sum:=ShallowCopy(x);
  sum.coefficients:=Concatenation(x.coefficients,y.coefficients);
  return AlgebraEltOps.Normalize(sum);
end;

AlgebraEltOps.\-:=function(x,y)return x+(-1)*y;end;

AlgebraEltOps.\*:=function(x,y) local A,res,i,j,m;
  if IsList(x) then return List(x,a->a*y);
  elif IsList(y) then return List(y,a->x*a);
  elif not IsRec(x) then 
    res:=ShallowCopy(y);
    res.coefficients:=List(y.coefficients,i->[x*i[1],i[2]]);
  elif not IsRec(y) then
    res:=ShallowCopy(x);
    res.coefficients:=List(x.coefficients,i->[y*i[1],i[2]]); A:=x.domain;
  else   
    A:=x.domain;
    if A.identification <> y.domain.identification then 
      Error("Elements must lie in the same algebra");
    fi; 
    if IsBound(x.info) and (not IsBound(y.info) or x.info<>y.info) then
      Error("Elements must belong to the same basis");
    fi;
    res:=ShallowCopy(x);res.coefficients:=[];
    if IsBound(x.info) then m:=x.info.multiplication;
    else m:=A.multiplication;fi;
    for i in x.coefficients do for j in y.coefficients do 
      Append(res.coefficients,List(m(i[2],j[2]),k->[i[1]*j[1]*k[1],k[2]]));
    od; od;
  fi;
  return AlgebraEltOps.Normalize(res);
end;    

AlgebraEltOps.\^:=function(x,n)local res,A;
  A:=x.domain;
  if n<0 then
    x:=A.EltToVector(x*A.basis);
    if DeterminantMat(x)=0 then Error("x is not invertible");fi;
    x:=A.VectorToElt(A.EltToVector(A.one)*x^-1);
    n:=-n;
  fi;
  if n=0 then return A.one;fi;
  if n=1 then return x;fi;
  res:=false; 
  while n>0 do
    if n mod 2 <> 0 then 
      if res=false then res:=x;else res:=res*x;fi;
    fi;
    x:=x*x;
    n:=QuoInt(n,2);
  od;
  return res;
end;

AlgebraEltOps.Degree:=function(r)local l;
  if not IsBound(r.degree) then 
    l:=CharacterDegrees(r.domain.group);
    r.degree:=Sum(r.coefficients, i-> i[1]*l[i[2]]);
  fi;
  return r.degree;
end;

AlgebraEltOps.Representative:=function(r)local Q;
  Q:=r.domain.field;
  return Sum(r.coefficients, i-> i[1]*Indeterminate(Q)^(i[2]-1));
end;

AlgebraEltOps.Coefficients:=function(r)local A,res,c;
  A:=r.domain; res:=[1..A.dimension]*0;
  for c in r.coefficients do res[c[2]]:=c[1];od;
  return res;
end;

#AlgebraEltOps.in:=function(x,A) 
#  return x.domain.identification = A.identification;
#end;

############################
## Operations on algebras ##
############################

FDAlgebraOps:=OperationsRecord("FDAlgebraOps");

# StructureConstants(A[,info][,structureconstants])
FDAlgebraOps.StructureConstants:=function(arg)local A,d,info;
  A:=arg[1];arg:=arg{[2..Length(arg)]};
  if Length(arg)>0 and IsString(arg[1]) then
    info:=A.(arg[1]);arg:=arg{[2..Length(arg)]};
  else info:=A;fi;
  if not IsBound(info.structureconstants) then 
    if Length(arg)>0 then info.structureconstants:=arg[1];
    else d:=A.dimension;
      info.structureconstants:=List([1..d],i->List([1..d],
           j->info.multiplication(i,j)));
    fi;
    info.multiplication:=function(i,j)return info.structureconstants[i][j];end;
  fi;
  return info.structureconstants;
end;

FDAlgebraOps.underlyingspace:=function(A)local d; d:=A.dimension;
  A.basis:=List([1..d],i->AlgebraElement(A,[[A.field.one,i]]));
  A.vectorspace:=(A.field)^d;
  A.EltToVector:=function(x) local res;
    if IsList(x) then return List(x,A.EltToVector);fi;
    if Length(x.coefficients)=0 then return Zero(A.vectorspace);fi;
    return Sum(x.coefficients,i->i[1]*A.vectorspace.basis.vectors[i[2]]);
  end;
  A.VectorToElt:=function(vec)
    if IsMat(vec) or (vec=[] and A.dimension<>0) then 
      return List(vec,A.VectorToElt);fi;
    return AlgebraElement(A,
      Filtered(List([1..Length(vec)],i->[vec[i],i]),x->x[1]<>0*x[1]));
  end;
# shrink list to independent vectors
  A.IndependentEntries:=list->A.VectorToElt(CanonicalBasis(
    Subspace(A.vectorspace,A.EltToVector(list))).vectors);
# repeatedly multiply on the left  list by gens
  A.SaturateLeft:=function(gens,list)local dim;
    list:=A.IndependentEntries(list);
    repeat dim:=Length(list); 
      list:=A.IndependentEntries(Concatenation(List(gens,i->i*list)));
    until Length(list)=dim;
    return list;
  end;
# repeatedly multiply on the right  list by gens
  A.SaturateRight:=function(gens,list)local dim;
    list:=A.IndependentEntries(list);
    repeat dim:=Length(list); 
      list:=A.IndependentEntries(Concatenation(List(gens,i->list*i)));
    until Length(list)=dim;
    return list;
  end;
end;

FDAlgebraOps.\in:=function(x,R) 
  if not IsAlgebraElement(x) then return false;
  else return x.domain.identification=R.identification;
  fi;
end;

FDAlgebraOps.MinimalPolynomial:=function(R,x)local res,V,y,vects,v;
  V:=R.vectorspace; vects:=R.EltToVector([R.one]);
  y:=x;v:=R.EltToVector(y);
  while not (v in Subspace(V,vects)) do
    Add(vects,v); y:=x*y;v:=R.EltToVector(y);
  od;
  res:=-Coefficients(Basis(Subspace(V,vects),vects),v);
  Add(res,R.field.one);
  return Polynomial(R.field,res);
end;

# Basis(A,name[,rec])
FDAlgebraOps.Basis:=function(arg)local A,l,m,name,info;
  A:=arg[1];
  if Length(arg)=1 then return A.basis;fi;
  name:=arg[2];
  if Length(arg)=2 then
    if not IsBound(A.(name)) then Error(A," has no ",name," basis");
    else return A.(name).basis;
    fi;
  fi;
  if IsList(arg[3]) then A.(name):=rec(value:=arg[3]);
  else A.(name):=arg[3];
  fi;
  info:=A.(name); 
  if not IsBound(info.basisname) then info.basisname:=name;fi;
  if not IsBound(info.parameters) then info.parameters:=[1..A.dimension];fi;
  l:=info.value;
  if Length(l)<>A.dimension then Error(l," is not a basis\n");fi;
  info.m:=List(l,function(r)local v,c;v:=[1..A.dimension]*0;
    for c in r.coefficients do v[c[2]]:=c[1];od;
    return v;end);
  if DeterminantMat(info.m)=0 then Error(l," is not a basis\n");fi;
  info.inv:=info.m^-1;
  info.basis:=List([1..A.dimension],i->AlgebraElement(A,[[1,i]],info));
  info.ord:=function(r)return A.basis*info.m[r];end;
  info.alt:=function(r)
    if Length(r.coefficients)=0 then return AlgebraElement(A,[],info);
    else return Sum(r.coefficients,x->x[1]*info.basis*info.inv[x[2]]);
    fi;
  end;
  info.one:=info.alt(A.one);
  info.multiplication:=function(x,y)
    return info.alt(info.ord(x)*info.ord(y)).coefficients;
  end;
  return info.basis;
end;

#########################################################################
## AlgebraHomomorphismByLinearity(A,B [,l[,nocheck]]) morphism A -> B  ##
## that sends A.basis to the list l of elements of B (default B.basis) ##
#########################################################################
AlgebraHomomorphismByLinearity:=function(arg) local A,B,l,f,d;
  A:=arg[1];B:=arg[2];
  if Length(arg)>=3 then l:=arg[3];else l:=B.basis;fi;
  d:=Dimension(A);
  if Length(l)<>d then Error("list of images must have length dim(A)");fi;
  f:=function(element) 
    if element.coefficients=[] then return B.zero;
    else return Sum(element.coefficients,i->i[1]*l[i[2]]);
    fi;
  end;
  if Length(arg)<4 then
     if ForAny([1..d],i->ForAny([1..d],j->
       f(A.basis[i])*f(A.basis[j])<>f(A.basis[i]*A.basis[j])))
     then Error("This is not a morphism of algebras");
     fi;
  fi;
  return f;
end; 

##########################################################
## The function SubAlgebra(A,l) computes the subalgebra ##
## of A generated by the list l                         ##
##########################################################
##
## SUB:=SubAlgebra(A,l) is an algebra
## It is endowed with SUB.injection : SUB -> A
##
SubAlgebra:=function(A,l) local subspace,SUB,base;
  SUB:=rec(operations:=OperationsRecord("SubAlgebraOps",FDAlgebraOps),
    field:=A.field,
    parent:=A,
    identification:=["subalgebra",A.identification,l]);
  SUB.operations.Print:=function(r)Print("SubAlgebra(",A,",",l,")");end;
  SUB.zero:=AlgebraElement(SUB,[]);
  SUB.basisinclusion:=A.SaturateLeft(Concatenation([A.one],l),
    Concatenation([A.one],l));
  SUB.dimension:=Length(SUB.basisinclusion);
  FDAlgebraOps.underlyingspace(SUB);
  subspace:=Subspace(A.vectorspace,A.EltToVector(SUB.basisinclusion));
  base:=CanonicalBasis(subspace);
  SUB.multiplication:=function(i,j)local ab;
    ab:=A.EltToVector(SUB.basisinclusion[i]*SUB.basisinclusion[j]);
    ab:=Coefficients(base,ab);
    ab:=List([1..Length(ab)], i-> [ab[i],i]);
    return Filtered(ab, i-> i[1] <> 0);
  end;
  SUB.injection:=AlgebraHomomorphismByLinearity(SUB,A,SUB.basisinclusion);
  SUB.belongs:=v->A.EltToVector(v) in subspace;
  SUB.restriction:=v->SUB.VectorToElt(Coefficients(base,A.EltToVector(v)));
  SUB.one:=SUB.restriction(A.one);
  SUB.string:=r->String(SUB.injection(r));
  SUB.generators:=List(l, SUB.restriction);
  return SUB;
end;

FDAlgebraOps.Centre:=function(A)local Q,SUB,comm,space,subspace,base,basis,l;
  Q:=A.field;
  SUB:=rec(operations:=OperationsRecord("FDSubAlgebraOps",FDAlgebraOps),
    field:=Q,
    parent:=A,
    basisname:="SUB",
    identification:=["centre",A.identification]);
  SUB.zero:=AlgebraElement(SUB,[]);
  SUB.one:=AlgebraElement(SUB,[[Q.one,1]]);
  space:=A.vectorspace;
  if IsBound(A.generators) then l:=A.generators; else l:=A.basis;fi;
  comm:=List(l,i->NullspaceMat(A.EltToVector(i*A.basis-A.basis*i)));
  subspace:=Intersection(List(comm,i->Subspace(space,i)));
  basis:=CanonicalBasis(subspace);
  base:=A.VectorToElt(basis.vectors);
  SUB.dimension:=Length(base);
  FDAlgebraOps.underlyingspace(SUB);
  SUB.multiplication:=function(i,j)local ab;
    ab:=Coefficients(basis,A.EltToVector(base[i]*base[j]));
    ab:=List([1..Length(ab)], i-> [ab[i],i]);
    return Filtered(ab, i-> i[1] <> 0);
  end;
  SUB.injection:=AlgebraHomomorphismByLinearity(SUB,A,base);
  SUB.belongs:=v->A.EltToVector(v) in subspace;
  SUB.restriction:=v->SUB.VectorToElt(Coefficients(basis,A.EltToVector(v)));
  SUB.one:=SUB.restriction(A.one);
  SUB.string:=r->String(SUB.injection(r));
  SUB.centre:=SUB;
  SUB.operations.Print:=function(r) Print("Centre(",A,")");end;
  A.centre:=SUB;
  return SUB;
end;

#######################################################################
## IsAbelian(A) returns true if A is commutative and false otherwise ##
#######################################################################
FDAlgebraOps.IsAbelian:=function(A)local l;
  if not IsBound(A.isAbelian) then
    if IsBound(A.generators) then l:=A.generators;else l:=A.basis;fi;
    A.isAbelian:=ForAll(l,x->ForAll(l,y->x*y=y*x));
  fi;
  return A.isAbelian;
end;

###########################################################################
## IsAssociative(A) returns true if A is associative and false otherwise ##
###########################################################################
FDAlgebraOps.IsAssociative:=function(A)local b;
  if not IsBound(A.isAssociative) then
    b:=A.basis;
    A.isAssociative:=ForAll(b,x->ForAll(b,y->ForAll(b,z->(x*y)*z=x*(y*z))));
  fi;
  return A.isAssociative;
end;

###################################################################
## The function CentralizerAlgebra(A,l) computes the centralizer ##
## of the elements of l in A                                     ##
###################################################################
##
## SUB:=CentralizerAlgebra(A,l) is an algebra
## It is endowed with SUB.injection : SUB -> A
##
###################################################################
CentralizerAlgebra:=function(A,l)local SUB,comm,subspace,base,basis;
  SUB:=rec(operations:=FDAlgebraOps,
    field:=A.field,
    parent:=A,
    basisname:="SUB",
    identification:=["centralizer",A.identification,l]);
  SUB.zero:=AlgebraElement(SUB,[]);
  comm:=List(l,i->NullspaceMat(A.EltToVector(i*A.basis-A.basis*i)));
  subspace:=Intersection(List(comm,i->Subspace(A.vectorspace,i)));
  basis:=CanonicalBasis(subspace);
  base:=A.VectorToElt(basis.vectors);
  SUB.dimension:=Length(base);
  SUB.operations.Print:=function(r) Print("Centralizer(",A,",",l,")");end;
  FDAlgebraOps.underlyingspace(SUB);
  SUB.multiplication:=function(i,j)local ab;
    ab:=Coefficients(basis,A.EltToVector(base[i]*base[j]));
    ab:=List([1..Length(ab)], i-> [ab[i],i]);
    return Filtered(ab, i-> i[1] <> 0);
  end;
  SUB.injection:=AlgebraHomomorphismByLinearity(SUB,A,base);
  SUB.belongs:=v->A.EltToVector(v) in subspace;
  SUB.restriction:=v->SUB.VectorToElt(Coefficients(basis,A.EltToVector(v)));
  SUB.one:=SUB.restriction(A.one);
  SUB.string:=r->String(SUB.injection(r));
  return SUB;
end;

IsCentralElement:=function(A,x)local l;
  if not (x in A) then Error("<x> must lie in <A>");fi;
  if IsBound(A.generators) then l:=A.generators; else l:=A.basis; fi;
  return ForAll(l,i->x*i=i*x);
end;

#########################################################################
##                           Ideals
## I.parent:=A
## I.basis is a K-basis of I
## I.dimension is the dimension of I
## I.(left or right)traces is the list of the traces of the elements of 
##   A.basis in their action on I
##########################################################################

#######################################################################
## LeftIdeal(A,vect) returns the left ideal I of A generated by vect ##
#######################################################################
LeftIdeal:=function(A,vect) local ideal,base;
  if IsAlgebraElement(vect) then vect:=[vect];fi;
  base:=A.IndependentEntries(vect);
  if IsBound(A.generators) then
       base:=A.SaturateLeft(Concatenation([A.one],A.generators),base);
  else base:=A.IndependentEntries(Concatenation(List(A.basis,i->i*vect)));
  fi;
  ideal:=rec(operations:=rec(Print:=
      function(r)Print("LeftIdeal(",A,",",vect,")");end),
    parent:=A,
    generators:=vect,
    basis:=base,
    base:=CanonicalBasis(Subspace(A.vectorspace,A.EltToVector(base))));
  ideal.dimension:=Length(ideal.basis);
  if Dimension(ideal)=0 then ideal.operations.LeftTrace:=x->0*A.field.one;
  else ideal.operations.LeftTrace:=x->TraceMat(List(ideal.base.vectors,v->
        Coefficients(ideal.base,A.EltToVector(x*A.VectorToElt(v)))));
  fi;
  return ideal;
end;

LeftTraces:=function(A,I) 
  if not IsBound(I.operations.LeftTrace) then 
    Error(I," must be a left ideal");
  fi;
  if not IsBound(I.lefttraces) then 
    I.lefttraces:=List(A.basis, i-> I.operations.LeftTrace(i));
  fi;
  return I.lefttraces;
end;

RightTraces:=function(A,I) 
  if not IsBound(I.righttraces) then 
    I.righttraces:=List(A.basis, i-> I.operations.RightTrace(i));
  fi;
  if not IsBound(I.operations.RightTrace) then 
    Error(I," must be a right ideal");
  fi;
  return I.righttraces;
end;

LeftIndecomposableProjectives:=function(A)
  if not IsBound(A.leftprojectives) then 
    A.leftprojectives:=List(Idempotents(A),i->LeftIdeal(A,[i]));
  fi;
  return A.leftprojectives;
end;

#########################################################################
## RightIdeal(A,vect) returns the right ideal I of A generated by vect ##
#########################################################################
RightIdeal:=function(A,vect) local ideal,base;
  if IsAlgebraElement(vect) then vect:=[vect];fi;
  base:=A.IndependentEntries(vect);
  if IsBound(A.generators) then
       base:=A.SaturateRight(Concatenation([A.one],A.generators),base);
  else base:=A.IndependentEntries(Concatenation(List(A.basis,i->vect*i)));
  fi;
  ideal:=rec(operations:=rec(Print:=function(r) 
     Print("RightIdeal(",A,",",vect,")");end),
    parent:=A,
    generators:=vect,
    basis:=base,
    base:=CanonicalBasis(Subspace(A.vectorspace,A.EltToVector(base))));
  ideal.dimension:=Length(ideal.basis);
  if Dimension(ideal)=0 then ideal.operations.RightTrace:=x->0*A.field.one;
  else ideal.operations.RightTrace:=x->TraceMat(List(ideal.base.vectors,v->
        Coefficients(ideal.base,A.EltToVector(A.VectorToElt(v)*x))));
  fi;
  return ideal;
end;

####################################################################
## TwoSidedIdeal(A,vect) two-sided ideal I of A generated by vect ##
####################################################################
TwoSidedIdeal:=function(A,vect)local ideal,base;
  if IsAlgebraElement(vect) then vect:=[vect];fi;
  base:=A.IndependentEntries(vect);
  if IsBound(A.generators) then
    base:=A.SaturateLeft(Concatenation([A.one],A.generators),base);
    base:=A.SaturateRight(Concatenation([A.one],A.generators),base);
  else 
    base:=A.IndependentEntries(Concatenation(List(A.basis,i->i*vect)));
    base:=A.IndependentEntries(Concatenation(List(A.basis,i->vect*i)));
  fi;
  ideal:=rec(operations:=rec(Print:=
      function(r)Print("TwoSidedIdeal(",A,",",vect,")");end),
    parent:=A,
    generators:=vect,
    basis:=base,
    base:=CanonicalBasis(Subspace(A.vectorspace,A.EltToVector(base))));
  ideal.dimension:=Length(ideal.basis);
  if Dimension(ideal)=0 then 
    ideal.operations.LeftTrace:=x->0*A.field.one;
    ideal.operations.RightTrace:=x->0*A.field.one;
  else 
    ideal.operations.LeftTrace:=x->TraceMat(List(ideal.base.vectors,v->
        Coefficients(ideal.base,A.EltToVector(x*A.VectorToElt(v)))));
    ideal.operations.RightTrace:=x->TraceMat(List(ideal.base.vectors,v->
        Coefficients(ideal.base,A.EltToVector(A.VectorToElt(v)*x))));
  fi;
  return ideal;
end;

VectorSpaceProjection:=function(quotient, vector) local res;
  res:=rec(operations:=SpaceCosetRowSpaceOps,isDomain:=true,isSpaceCoset:=true);
  res.representative:=vector;
  res.factorDen:=quotient.factorDen;
  res.representative:=CanonicalRepresentative(res);
  return res;
end;

QuotientAlgebra:=function(A,g) local I,i,j,k,res,VA,VI,B,base,basis,space;
  if IsAlgebraElement(g) or IsList(g) then I:=TwoSidedIdeal(A,g);
  else I:=g; fi;
  VA:=A.vectorspace;
  VI:=Subspace(VA,A.EltToVector(I.basis));
  B:=rec(operations:=OperationsRecord("QuotientFDAlgebraOps",FDAlgebraOps),
   field:=A.field,
   parent:=A,
   basisname:="QUOTIENT",
   type:="Quotient algebra",
   dimension:=Dimension(VA)-Dimension(VI),
   identification:=["Quotient algebra",A.identification,I.generators]);
  B.operations.Print:=function(r) Print("QuotientAlgebra(",A,",",I,")");end;
  B.zero:=AlgebraElement(B,[]);
  #SUB.one:=AlgebraElement(B,[[Q.one,1]]);
  space:=VA/VI;
  basis:=Basis(space);
  B.parameters:=[1..B.dimension];
  FDAlgebraOps.underlyingspace(B);
  B.multiplication:=function(i,j) local a,b,ab;
    a:=A.VectorToElt(CanonicalRepresentative(
      VectorSpaceProjection(space,basis.vectors[i])));
    b:=A.VectorToElt(CanonicalRepresentative(
      VectorSpaceProjection(space,basis.vectors[j])));
    ab:=VectorSpaceProjection(VA/VI,A.EltToVector(a*b));
    ab:=space.semiEchelonBasis.operations.Coefficients(basis,ab);
    ab:=List([1..Length(ab)], i-> [ab[i],i]);
    return Filtered(ab,i->i[1]<>0*i[1]);
  end;
  base:=List(A.basis,i->VectorSpaceProjection(space,A.EltToVector(i)));
  base:=List(base,i->space.semiEchelonBasis.operations.Coefficients(basis,i));
  base:=List(base,B.VectorToElt);
  B.projection:=AlgebraHomomorphismByLinearity(A,B,base);
  B.one:=B.projection(A.one);
  B.representative:=function(r) local i,rep;
    if r.coefficients=[] then return A.zero;fi;
    rep:=Sum(r.coefficients, i-> 
      i[1]*VectorSpaceProjection(space,basis.vectors[i[2]]));
    return A.VectorToElt(Representative(rep));
  end;
 #B.string:=r->SPrint("Class(",B.representative(r),")");
  B.operations.Idempotents:=R->Filtered(List(Idempotents(B.parent), 
    B.projection),i->i<>B.zero);
  if IsBound(A.generators) then 
    if Length(A.generators) < B.dimension then 
      B.generators:=List(A.generators, i-> B.projection(i));
    fi;
  fi;
  return B;
end;

#####################################################
##
## RadicalPower(A,n) computes Radical(A)^n
## This is again a two-sided ideal
## 
RadicalPower:=function(A,n) local res,subspace,can,ideal,base;
  if n=0 then return A;
  elif n=1 then return Radical(A);
  elif IsBound(A.radicalpowers) and Length(A.radicalpowers)>=n then 
    return A.radicalpowers[n];
  fi;
  res:=Concatenation(List(Radical(A).basis,i->i*RadicalPower(A,n-1).basis));
  base:=CanonicalBasis(Subspace(A.vectorspace,A.EltToVector(res)));
  ideal:=rec(operations:=rec(Print:=
    function(r)Print("TwoSidedIdeal(",A,",",r.generators,")");end),
    parent:=A,
    generators:=A.VectorToElt(base.vectors),
    basis:=A.VectorToElt(base.vectors));
  ideal.dimension:=Length(ideal.basis);
  if Dimension(ideal)=0 then 
    ideal.operations.LeftTrace:=x->0*A.field.one;
    ideal.operations.RightTrace:=x->0*A.field.one;
  else 
    ideal.operations.LeftTrace:=x->TraceMat(List(base.vectors,v->
	Coefficients(base,A.EltToVector(x*A.VectorToElt(v)))));
    ideal.operations.RightTrace:=x->TraceMat(List(base.vectors,v->
	Coefficients(base,A.EltToVector(A.VectorToElt(v)*x))));
  fi;
  A.radicalpowers[n]:=ideal;
  return ideal;
end;  

#####################################################
##
## LoewyLength(A) computes the Loewy length of A, 
## that is, the smallest natural number d such 
## that RadicalPower(A,d)=0
##
#####################################################
LoewyLength:=function(A) local d,res,gen;
  if not IsBound(A.loewylength) then 
    if Length(Radical(A).generators)=1 then 
      gen:=Radical(A).generators[1]; res:=gen;
      d:=1;while res<>A.zero do d:=d+1; res:=res*gen;od;
    else 
      d:=1;while RadicalPower(A,d).dimension>0 do d:=d+1; od;
    fi;
    A.loewylength:=d;
  fi;
  return A.loewylength;
end; 

#####################################################
## CharacterDecomposition(A,chi)  decomposition of ##
## the character chi of A as a sum of irreducibles ##
#####################################################
CharacterDecomposition:=function(A,character)local irr;
  irr:=CharTable(A).basistraces;
  irr:=Basis(Subspace((A.field)^Length(character), irr),irr);
  return Coefficients(irr, character);
end;

CentralIdempotents:=function(A) local idem,perm,mat;
  if not IsBound(A.centralidempotents) then 
    if IsBound(A.operations.CentralIdempotents) then 
      A.centralidempotents:=A.operations.CentralIdempotents(A);
    else 
      idem:=Idempotents(A);
      mat:=List(idem,i->List(LeftIndecomposableProjectives(A),
        j->Dimension(Subspace(A.vectorspace,A.EltToVector(i*j.basis)))));
      A.blocks:=DecomposedMat(mat);
      A.centralidempotents:=List(A.blocks,i->Sum(idem{i}));
      perm:=Concatenation(A.blocks);
      A.cartanmatrix:=mat{perm}{perm};
    fi;
  fi;
  return A.centralidempotents;
end;

FDAlgebraOps.CartanMatrix:=function(A)local rows;
  if not IsBound(A.cartanmatrix) then 
    CentralIdempotents(A);# computes A.blocks and A.cartanmatrix
  fi;
  if IsBound(A.operations.CharTable) and IsBound(CharTable(A).rows)
  then rows:=CharTable(A).rows; 
  else rows:=List([1..Length(Idempotents(A))],String);
  fi;
  rows:=rows{Concatenation(A.blocks)};
  return rec(operations:=rec(Print:=TablePrint), blocks:=A.blocks,
    rows:=rows, columns:=List(rows,i->[i]),
    matrix:=A.cartanmatrix, field:=A.field);
end;

######################################################################
## PolynomialQuotientAlgebra(P) returns the algebra A=K[T]/(P(T)) ##
######################################################################
##
## The function A.class sends a polynomial to its image in A
## An element of A is printed as "Class(polynomial)"
##
######################################################################

PolynomialQuotientAlgebra:=function(P) local Q,A,q;
  Q:=P.baseRing;q:=Indeterminate(Q);
  if not IsBound(q.name)then q.name:="q";fi;
  A:=rec(operations:=OperationsRecord("PQAlgebraOps",FDAlgebraOps),
   field:=Q,
   dimension:=Degree(P),
   isAbelian:=true,
   basisname:="Class",
   identification:=["PolynomialQuotientAlgebra",P]);
  A.operations.Print:=function(A)Print(Q,"[",q.name,"] / (",P,")");end; 
  FDAlgebraOps.underlyingspace(A);  
  A.zero:=AlgebraElement(A,[]);
  A.one:=AlgebraElement(A,[[Q.one,1]]);
  A.string:=r->SPrint("Class(",Sum(r.coefficients, i-> i[1]*q^(i[2]-1)),")");
  A.class:=function(pol) local res; res:=EuclideanRemainder(pol,P);
    return AlgebraElement(A,List([1..Length(res.coefficients)], 
      i-> [res.coefficients[i],res.valuation+i]))+A.zero;
  end;
  A.multiplication:=function(i,j)
    return A.class(EuclideanRemainder(q^(i+j-2),P)).coefficients;
  end;
  FDAlgebraOps.StructureConstants(A);
  A.generators:=[A.basis[2]];
  A.centre:=A;
  A.operations.Radical:=function(B) local res,der;
    if Q.char=0 then 
      der:=Derivative(P); res:=Gcd(P/Gcd(der,P),der);
      A.radical:=TwoSidedIdeal(A,[A.class(res)]);
      A.radicalpowers:=[A.radical];
      return A.radical;
    fi;
    res:=Product(List(Filtered(Collected(Factors(P)),i->i[2]>1),i->i[1]));
    A.radical:=TwoSidedIdeal(A,[A.class(res)]);
    A.radicalpowers:=[A.radical];
    return A.radical;
  end;
  A.operations.Idempotents:=function(A)local i,result,a,b,bezout,res;
    result:=[];
    for i in Collected(Factors(P)) do 
      a:=i[1]^i[2]; b:=P/a;
      bezout:=GcdRepresentation(a,b);
      res:=bezout[1]*a+bezout[2]*b;
      res:=res.coefficients[1];
      res:=bezout[2]*b/res;
      Add(result,A.class(res));
    od;
    return result;
  end;
  A.operations.CentralIdempotents:=Idempotents;
  A.isGroup:=true;
  A.operations.CharacterDegrees:=function(A)
    A.characterDegrees:=List([1..Length(Set(Factors(P)))], i-> 1);
    return A.characterDegrees;
  end;
  A.operations.CharTable:=function(A) local table,i,j,k,factors,mat,deg,coefs;
    factors:=Set(Factors(P));
    table:=rec(domain:=A,field:=A.field,
     operations:=rec(Print:=TablePrint),
     characterDegrees:=List(factors, Degree),
     columns:=List(factors,i-> [String(i)]),
     rows:=List([1..A.dimension], i-> String(q^(i-1))));
    table.matrix:=NullMat(A.dimension,A.dimension);
    for k in [1..Length(factors)] do 
      coefs:=Concatenation([1..factors[k].valuation]*A.field.zero,
         factors[k].coefficients);
      deg:=Degree(factors[k]);
      mat:=List([1..deg],i->[1..deg]*A.field.zero);
      for i in [1..deg-1] do mat[i+1][i]:=1;od;
      for i in [1..deg] do mat[i][deg]:=-coefs[i];od;
      for j in [1..A.dimension] do 
        table.matrix[j][k]:=TraceMat(mat^(j-1));
      od;
    od;
    table.irreducibles:=table.matrix;
    table.basistraces:=table.matrix;
    return table;
  end;
  A.operations.CartanMatrix:=function(A)local factors;
    factors:=Collected(Factors(P));
    return rec(field:=Rationals,domain:=A,
      operations:=rec(Print:=TablePrint),
      rows:=List(factors,i->String(i[1])),
      columns:=List(factors,i->[String(i[1])]),
      matrix:=DiagonalMat(List(factors,x->x[2])));
  end;
  return A;
end;

####################################################################
## GroupAlgebra(G,Q) renvoie l'algebre de groupe A=Q[G] (si Q est 
## omis, on prend pour Q le corps des rationnels)
##
## En plus des champs existants, A contiendra : 
##   A.group = G
##   A. structureconstants = constantes de structures
## De plus :
##   A.parameters = [1..Size(G)] 
##   A.basisname  = "e" (par defaut)
## 
## Les elements seront representes sous la 
## forme a1 * e(i1) + a2* e(i2) + ... ou ak est un element 
## du corps et e(ik) designe l'element numerote ik dans la 
## liste des elements du groupe
##
## La fonction GroupAlgebra munit G du champ 
##   G.law 
## codant la table de multiplication de G
##
## La fonction Augmentation code le morphisme d'augmentation
##################################################################

GroupAlgebraOps:=OperationsRecord("GroupAlgebraOps",FDAlgebraOps);

GroupAlgebraOps.CharacterDegrees:=function(A) 
  if A.field.char=0 then 
    A.characterDegrees:=CharacterDegrees(A.group);
  else Error("Cannot compute CharacterDegrees in positive characteristic");
  fi;
  return A.characterdegrees;
end;

GroupAlgebraOps.CharTable:=function(A)local conj;
  if A.field.char>0 then 
    Error("Cannot compute CharTable in positive characteristic");
  fi;
  A.charTable:=CharTable(A.group);
  conj:=List(Elements(A.group),i->PositionClass(A.group,i));
  A.charTable.basistraces:=List(A.charTable.irreducibles,chi->chi{conj});
  A.charTable.Print:=function(table)Display(table);end;
  return A.charTable; 
end;

GroupAlgebraOps.Radical:=function(A) 
  if A.field.char<>0 then 
    Error("Cannot compute CentralIdempotents in positive characteristic");
  fi;
  return TwoSidedIdeal(A,[]);
end;

GroupAlgebraOps.Print:=function(A) 
  Print("GroupAlgebra(",A.group,",",A.field,")");
end;

GroupAlgebra:=function(arg) local G,Q,A,e,d;
  G:=arg[1];
  if Length(arg)=1 then Q:=Rationals; else Q:=arg[2]; fi;
  e:=Elements(G); d:=Length(e);
  A:=rec(operations:=Copy(GroupAlgebraOps),
    field:=Q,
    group:=G,
    type:="Group algebra",
    dimension:=d,
    parameters:=[1..d],
    basisname:="e");
  A.identification:=[A.type,A.group,A.field];
  A.zero:=AlgebraElement(A,[]);
  A.one:=AlgebraElement(A,[[Q.one,Position(e,G.identity)]]);
  FDAlgebraOps.underlyingspace(A);
  A.generators:=List(G.generators,i->A.basis[Position(e,i)]);
  if d < 1201 then 
    if not IsBound(G.law) then 
      G.law:=List([1..d],i->OnTuples([1..d],PermListList(e,e[i]*e)));
    fi;
    FDAlgebraOps.StructureConstants(A,
      List([1..d],i->List([1..d],j->[[Q.one,G.law[i][j]]])));
  else 
    A.multiplication:=function(i,j) return [[Q.one,Position(e,e[i]*e[j])]];end;
  fi;
  if Q.char=0 then 
    A.radical:=TwoSidedIdeal(A,[]);
    A.radicalpowers:=[A.radical];
  fi;
  A.isGroup:=true;
  A.embedding:=function(g)return A.basis[Position(e,g)];end;
  return A;
end;

GroupAlgebraOps.CentralIdempotents:=function(A)local conj,e,id,pid;
  if A.field.char=0 then 
    e:=Elements(A.group);pid:=Position(e,A.group.identity);
    conj:=List(Elements(A.group),i->PositionClass(A.group,i^-1));
    return List(CharTable(A).irreducibles,
      chi->chi[conj[pid]]*chi{conj}*A.basis*Size(A.group)^-1);
  fi; 
  id:=CentralIdempotents(GroupAlgebra(A.group));
  pid:=List(PBlocks(A.group,A.field.char),i->Sum(id{i}));
  return List(pid,i->AlgebraElement(A,List(i.coefficients, 
    j->[CyclotomicModP(j[1],A.field.char),j[2]]))+A.zero);
end;

Augmentation:=function(r)return Sum(r.coefficients,i->i[1]);end;

## GrothendieckRing(G,Q) renvoie l'anneau de Grothendieck A=Q Irr(G) 
##(si Q est omis, on prend pour Q le corps des rationnels)
##
## En plus des champs existants, A contiendra : 
##   A.group = G
## De plus :
##   A.parameters = [1..Length(CharacterDegrees(G))] 
##   A.basisname  = "CHI" (par defaut)
## 
## Les elements seront representes sous la 
## forme a1 * CHI(i1) + a2 * CHI(i2) + ... ou ak est un element 
## du corps et CHI(ik) designe le caractere numerote ik dans la 
## table de caracteres du groupe
##
## La fonction GrothendieckRing munit G du champs 
##   G.tensorproducts
## codant la decomposition de produits tensoriels de 
## caracteres irreductibles et du champs 
##   G.grothendieckring = A
##
## La fonction Degre envoie un element de l'anneau de 
## Grothendieck sur son degre (virtuel a  priori bien sur)

GrothendieckRingOps:=OperationsRecord("GrothendieckRingOps",FDAlgebraOps);

GrothendieckRingOps.CharacterDegrees:=function(A)local l,i;
  if A.field.char=0 then 
    A.characterDegrees:=List([1..A.dimension], i-> 1);
  else 
    l:=ConjugacyClasses(A.group);
    l:=Filtered(l, i-> Gcd(Order(A.group,Representative(i)),A.field.char)=1);
    A.characterDegrees:=List([1..Length(l)], i-> 1);
  fi;
  return A.characterDegrees;
end;

GrothendieckRingOps.CharTable:=function(A) local orb,res,rows;
  orb:=ConjugacyClasses(A.group);
  res:=rec(
    inverse:=List(orb,i->PositionClass(A.group,Representative(i)^-1)),
    field:=A.field,
    operations:=rec(Print:=TablePrint),
    columns:=List([1..A.dimension],i->[SPrint("X.",i)]),
    domain:=A,
    irreducibles:=TransposedMat(CharTable(A.group).irreducibles),
    centralizers:=CharTable(A.group).centralizers);
  if A.field.char=0 then rows:=[1..A.dimension];
  else rows:=Filtered([1..Length(orb)], i-> 
     Gcd(Order(A.group,Representative(orb[i])),A.field.char)=1);
    res.irreducibles:=List(rows,
          i->CyclotomicModP(res.irreducibles[i],A.field.char));
    res.centralizers:=res.centralizers{rows};
  fi;
  res.parameters:=rows;
  res.rows:=List(rows,i->SPrint("MU.",i));
  res.basistraces:=res.irreducibles;
  res.matrix:=res.irreducibles;
  return res;
end;

GrothendieckRingOps.Print:=function(A)
  Print("GrothendieckRing(",A.group,",",A.field,")");
end;

GrothendieckRingOps.Radical:=function(A)local p,e,mat,res;
  if A.field.char>0 then 
    if (Size(A.group) mod A.field.char) <> 0 then 
      A.radical:=TwoSidedIdeal(A,[]);
      A.radicalpowers:=[A.radical];
    else 
      p:=A.field.char;
      mat:=A.EltToVector(List(A.basis,i->i^p));
      e:=1;
      while p^e < A.dimension do 
        e:=e+1;
        mat:=mat^p;
      od;
      res:=NullspaceMat(mat);
      if Length(res)=0 then res:=[A.zero];else res:=A.VectorToElt(res);fi;
      A.radical:=LeftIdeal(A,res);
      A.radicalpowers:=[A.radical];
    fi;
  fi;
  return A.radical;
end;

GrothendieckRingOps.CartanMatrix:=function(A)local res;
  res:=rec(operations:=rec(Print:=TablePrint),field:=Rationals);
  if A.field.char=0 then 
    res.matrix:=IdentityMat(A.dimension);
    res.rows:=List([1..A.dimension],i->SPrint("MU.",i));
    res.columns:=List([1..A.dimension],i->[SPrint("IND.",i)]);
  else 
    res.matrix:=DiagonalMat(List(PiPrimeSections(A.group,A.field.char),Length));
    res.rows:=List(CharTable(A).parameters,i->SPrint("MU.",i));
    res.columns:=List(CharTable(A).parameters,i->[SPrint("IND.",i)]);
  fi;
  return res;
end;

GrothendieckRing:=function(arg)local G,Q,A,d,T,irr;
  G:=arg[1];
  if Length(arg)=1 then Q:=Rationals; else Q:=arg[2]; fi;
  if IsBound(G.tensorproducts) then 
    d:=Length(G.tensorproducts);
  else
    T:=CharTable(G);
    irr:=T.irreducibles;
    d:=Length(irr);
  fi;
  A:=rec(operations:=Copy(GrothendieckRingOps),
    group:=G,
    field:=Q,
    type:="Grothendieck ring",
    dimension:=d,
    parameters:=[1..d],
    basisname:="X",
    isAbelian:=true);
  A.zero:=AlgebraElement(A,[]);
  A.identification:=[A.type,A.group,A.field];
  FDAlgebraOps.underlyingspace(A);
  if not IsBound(G.tensorproducts) then 
    G.tensorproducts:=
       List([1..d],i->MatScalarProducts(T,Tensored([irr[i]],irr),irr));
  fi;
  A.one:=AlgebraElement(A,[[Q.one,Position(G.tensorproducts,IdentityMat(d))]]);
  FDAlgebraOps.StructureConstants(A,List([1..d], i-> List([1..d], j-> 
      Filtered(List([1..d], k->[G.tensorproducts[i][j][k],k]),l->l[1]<>0))));
  #CharacterDegrees(G);
  A.isGroup:=true;
  if Q.char=0 then 
    A.radical:=TwoSidedIdeal(A,[]);
    A.radicalpowers:=[A.radical];
  fi;
  A.centre:=A;
  return A;
end;

GrothendieckRingOps.Idempotents:=function(A) local T,e;
  if A.field.char=0 then 
    T:=CharTable(A);
    return List([1..A.dimension], i-> 
      (T.irreducibles[T.inverse[i]]/T.centralizers[i])*A.basis);
  fi;
  e:=Idempotents(GrothendieckRing(A.group));
  e:=List(PiPrimeSections(A.group,A.field.char),i->Sum(e{i}));
  return List(e,i->AlgebraElement(A,List(i.coefficients,
    j->[CyclotomicModP(j[1],A.field.char),j[2]]))+A.zero);
end;
      
GrothendieckRingOps.CentralIdempotents:=Idempotents;

GrothendieckRingOps.PrincipalIdempotent:=function(A)local G,e,c,res;
  G:=A.group;
  e:=Idempotents(GrothendieckRing(G));
  c:=PositionClass(G,G.identity);
  res:=Sum(e{First(PiPrimeSections(G,A.field.char),i->c in i)});
  return AlgebraElement(A,List(res.coefficients,i-> 
    [A.field.one*Numerator(i[1])/Denominator(i[1]),i[2]]));
end;

GrothendieckRingOps.PrincipalLoewyLength:=function(A)local e;
  e:=GrothendieckRingOps.PrincipalIdempotent(A);
  LoewyLength(A);
  e:=List(A.radicalpowers,i->LeftIdeal(A,i.basis*e));
  A.principalloewyseries:=Filtered(e,i->i.dimension>0);
  return Length(A.principalloewyseries)+1;
end;

RestrictionHomomorphism:=function(A,B)local k;
  k:=InductionTable(B.group,A.group).scalar;
  return AlgebraHomomorphismByLinearity(A,B,k*B.basis);
end;

AdamsOperation:=function(A,n)local r,G;
  G:=A.group;
  r:=List(ConjugacyClasses(G),x->PositionClass(G,Representative(x)^n));
  r:=List(G.charTable.irreducibles,chi->chi{r});
  r:=MatScalarProducts(G.charTable,G.charTable.irreducibles,r);
  return AlgebraHomomorphismByLinearity(A,A,r*A.basis);
end;

## GroupAlgebraCentre(G,Q) renvoie le centre de l'algebre de groupe Q[G] 
##(si Q est omis, on prend pour Q le corps des rationnels)
GroupAlgebraCentreOps:=OperationsRecord("GroupAlgebraCentreOps",FDAlgebraOps);

GroupAlgebraCentreOps.Print:=function(A) 
  Print("GroupAlgebraCentre(",A.group,",",A.field,")");
end;

GroupAlgebraCentreOps.Radical:=function(A) local p,e,mat,res;
  if A.field.char > 0 then 
    if (Size(A.group) mod A.field.char) <> 0 then 
      A.radical:=TwoSidedIdeal(A,[]);
      A.radicalpowers:=[A.radical];
    else 
      p:=A.field.char;
      mat:=A.EltToVector(List(A.basis,i->i^p));
      e:=1;
      while p^e<A.dimension do 
        e:=e+1;
        mat:=mat^p;
      od;
      res:=NullspaceMat(mat);
      if Length(res)=0 then res:=[A.zero];else res:=A.VectorToElt(res);fi;
      A.radical:=LeftIdeal(A,res);
      A.radicalpowers:=[A.radical];
    fi;
  fi;
  return A.radical;
end;

GroupAlgebraCentre:=function(arg) local G,Q,A,d,T,irr,prod,i,j,invs,pid;
  G:=arg[1];
  if Length(arg)=1 then Q:=Rationals; else Q:=arg[2]; fi;
  T:=CharTable(G);d:=Length(T.irreducibles);
  A:=rec(operations:=Copy(GroupAlgebraCentreOps),
    group:=G,
    field:=Q,
    basisname:="Cl",
    type:="Group algebra centre",
    dimension:=d,
    parameters:=[1..d]);
  if IsBound(T.classnames) then A.parameters:=T.classnames;fi;
  A.identification:=[A.type,A.group,A.field];
  A.zero:=AlgebraElement(A,[]);
  FDAlgebraOps.underlyingspace(A);
  pid:=PositionClass(G,G.identity);
  A.one:=AlgebraElement(A,[[Q.one,pid]]);
  if not IsBound(G.groupalgebracentre) then 
    invs:=List(ConjugacyClasses(G),x->PositionClass(G,Representative(x)^-1));
    prod:=List([1..d],i->[]);
    irr:=T.irreducibles;
    for i in [1..d] do 
      for j in [1..i] do 
        prod[i][j]:=T.classes[i]*T.classes[j]*
	  TransposedMat(irr){invs}*List(irr,l->l[i]*l[j]/l[pid])/Size(G);
	prod[j][i]:=prod[i][j];
      od;
    od;
    G.groupalgebracentre:=prod;
  fi;
  FDAlgebraOps.StructureConstants(A,List([1..d], i-> List([1..d], j-> 
     Filtered(List([1..d], k->[G.groupalgebracentre[i][j][k],k]),l->l[1]<>0))));
#  A.isGroup:=true;
  if Q.char=0 then 
    A.radical:=TwoSidedIdeal(A,[]);
    A.radicalpowers:=[A.radical];
  fi;
  A.isAbelian:=true;
  A.centre:=A;
  return A;
end;

#################################################################################
## La fonction SolomonAlgebra(W,Q) construit l'algebre de Solomon 
## d'un groupe de Coxeter fini W sur le corps Q (s'il est omis, 
## la fonction prend pour Q le corps des rationnels). Le resultat 
## est une algebre A munie en plus des champs suivants (ici, r est 
## le rang semi-simple de W) :
##   A.group  = W
##   A.xbase  = fonction qui envoie la partie I de [1..r] sur x_I
##   A.ybase  = fonction qui envoie la partie I de [1..r] sur y_I
##   A.radical = radical de A : 
##        A.radical.base = base de A.radical (x_I-x_J)_{I et J conjugues}
##        A.radical.dimension = dimension de A.radical
##   A.radicalpowers = [A.radical]
##   A.injection = injection dans l'algebre de groupe (disponible 
##                 si |W| <= 2000 : jusqu'a D5 donc...)
##
## De plus
##   A.base       = base x_I
##   A.parameters = parties de [1..r] (comme nombres : ex. 123 pour [1,2,3])
##   A.basename   = "X" (par defaut)
##   
## La fonction SolomonAlgebra munit W du record W.solomon contenant 
## les champs suivants :
##   W.solomon.mackey = constantes de structure (coeff de x_I * x_J sur x_K)
##   W.solomon.conjugacy = classes de conjugaison de parties de [1..r] 
##                         reperees par leur position dans A.parameters
##   W.solomon.domain = A
##   W.solomon.subsets = parties de [1..r]


SolomonAlgebraOps:=Copy(FDAlgebraOps);

SolomonAlgebraOps.Print:=function(A) 
  if A.type="Solomon algebra" then 
    Print("SolomonAlgebra(",A.group,",",A.field,")");
  else 
    Print("GeneralizedSolomonAlgebra(",A.group,",",A.field,")");
  fi;
end;

SolomonAlgebraOps.CharacterDegrees:=function(A) 
  if A.field.char=0 then 
    if A.type="Solomon algebra" then 
      A.characterDegrees:=List(A.group.solomon.conjugacy, i-> 1);
    elif A.type="Generalized Solomon algebra" then 
      A.characterDegrees:=List(A.group.generalizedsolomon.conjugacy, i-> 1);
    else
      Error("<A> must be a Solomon algebra");
    fi;
  else Error("Cannot compute character degrees of Solomon algebra in positive characteristic");
  fi;
  return A.characterDegrees;
end;

SolomonAlgebraOps.CharTable:=function(A)local W,conj,orb,cox,res,irr,inc;
  W:=A.group;
  conj:=List(W.solomon.conjugacy,i->i[1]);
  orb:=List(conj,i->ReflectionSubgroup(W,W.solomon.subsets[i]));
  cox:=List(orb,function(g)
    if g.generators=[] then return ();else return Product(g.generators);fi;
    end);
  res:=rec(field:=A.field,
    domain:=A,
    irreducibles:=List(conj,i->List(W.solomon.mackey{conj},function(j) 
      if IsBound(j[i][i]) then return j[i][i]; else return 0; fi;end)),
    rows:=List(cox,i->IntListToString(CoxeterWord(W,i))),
    columns:=List(cox,i->List(CoxeterWord(W,i),String)),
    operations:=rec(Print:=TablePrint));
  res.columns[Length(res.columns)]:=["0"];
  res.basistraces:=List([1..Length(res.irreducibles)],function(ii)local r;
    r:=Concatenation(List([1..Length(W.solomon.conjugacy)],
      i->List(W.solomon.conjugacy[i], j-> [j,res.irreducibles[ii][i]])));
    Sort(r); return List(r,i->i[2]);end);
  res.matrix:=res.irreducibles;
  if A.field.char=0 then return res;fi;
  irr:=CyclotomicModP(res.irreducibles,A.field.char);
  inc:=CollectBy([1..Length(irr)],irr);Sort(inc);
  res.rows:=List(inc, i-> res.rows[i[1]]);
  res.irreducibles:=List(inc, i-> irr[i[1]]);
  res.matrix:=res.irreducibles;
  res.basistraces:=CyclotomicModP(res.basistraces,A.field.char);
  res.parameters:=List(inc, i-> A.parameters[W.solomon.conjugacy[i[1]][1]]);
  res.invertiblematrix:=List(inc,i-> List(inc,j->irr[i[1]][j[1]]));
  res.incidence:=inc;
  return res;
end;

SolomonAlgebraOps.injection:=function(A) local B,W,X;
  W:=A.group;
  B:=GroupAlgebra(A.group,A.field);
  if A.type="Solomon algebra" then 
    X:=List(W.solomon.subsets, i-> 
      Sum(ReducedRightCosetRepresentatives(W, ReflectionSubgroup(W,i)), 
        i-> B.embedding(i^-1)));
  else
    X:=List(W.generalizedsolomon.subgroups, i-> 
      Sum(ReducedRightCosetRepresentatives(W, i),i-> B.embedding(i^-1)));
  fi;
  A.injection:=AlgebraHomomorphismByLinearity(A,B,X,"no check");
  return A.injection;
end;

# SolomonAlgebra(W[,K]) K is the field
SolomonAlgebra:=function(arg) local W,Q,S,A,r,d,n,l,I,inclusion,f,i,j,k,
  c,II,JJ,conj,orb,XX,B,subspace,base,trad,Isorted;
  W:=arg[1];
  if Length(arg)=1 then Q:=Rationals; else Q:=arg[2]; fi;
  A:=rec(operations:=Copy(SolomonAlgebraOps),
   group:=W,
   field:=Q,
   type:="Solomon algebra");
  A.identification:=[A.type,A.group,A.field];
  A.zero:=AlgebraElement(A,[]);
  A.one:=AlgebraElement(A,[[Q.one,1]]);  
  n:=Size(W);
  r:=W.semisimpleRank;
  S:=W.reflections{[1..r]};
  d:=2^r;
  A.dimension:=d;
  FDAlgebraOps.underlyingspace(A);
  if r <=7 then A.generators:=List([2..r+1],i->A.basis[i]);fi;
  I:=Concatenation(List([r,r-1..0],i->Combinations([1..r],i)));
  inclusion:=List([1..d], i-> List(Combinations(I[i]),j->Position(I,j)));
  if r < 10 then A.parameters:=List(I,IntListToString);
  else A.parameters:=I;
  fi;
  ###################################
  if not IsBound(W.solomon) then 
    if n>2000 then
      InfoAlgebra("# Computing structure constants...\n");
      InfoAlgebra("# Nb. elements to examine: ",n,".\c");
    fi;
    c:=n+1; 
    Isorted:=Set(I);trad:=List(Isorted,j->Position(I,j));
    W.solomon:=rec(mackey:=List([1..d],i->List([1..d],j->[])));
    ForEachElement(W,function(w)local d,g,i,j,l,k,II,JJ,p;
      g:=Position(I,Difference([1..r],LeftDescentSet(W,w)));
      d:=Position(I,Difference([1..r],RightDescentSet(W,w)));
      l:=List([1..r],i->Position(S,S[i]^w));
      for i in inclusion[g] do II:=Set(l{I[i]});
	for j in inclusion[d] do
	  JJ:=ShallowCopy(II);
	  IntersectSet(JJ,I[j]);
	  k:=trad[Position(Isorted,JJ)];p:=W.solomon.mackey[i][j];
	  if IsBound(p[k]) then p[k]:=p[k]+1; else p[k]:=1; fi;
	od;
      od;
      c:=c-1; if c mod 2000=0 then InfoAlgebra(c,".\c");fi;
    end);
    if n>2000 then InfoAlgebra("...done!\n");fi;
  fi;
  ###################################
  W.solomon.subsets:=I;
  if n>10000 then InfoAlgebra("# Computing orbits of parabolic subgroups...\c");fi;
  orb:=[1..d];
  if not IsBound(W.solomon.conjugacy) then 
    conj:=[];
    while Length(orb)>0 do 
      i:=orb[1];
      j:=Filtered(orb,k->IsBound(W.solomon.mackey[i][k][k]) and 
		   Length(I[i])=Length(I[k]));
      Add(conj,j);
      orb:=Difference(orb,j);
    od;
    W.solomon.conjugacy:=conj;
  fi;
  if n>10000 then InfoAlgebra("done!\n");fi;
  ###################################
  A.xbasis:=function(arg) local a;
    if Length(arg)=1 then 
      if IsList(arg[1]) then a:=Set(arg[1]);
      else a:=Set(Digits(arg[1]));
      fi;
    else a:=Set(arg);
    fi;
    return A.basis[Position(W.solomon.subsets,a)];
  end;
  if A.field.char=0 then 
    A.xprimebasis:=function(arg) local a,i;
      if Length(arg)=1 then 
        if IsList(arg[1]) then a:=Set(arg[1]);
        else a:=Set(Digits(arg[1]));
        fi;
      else a:=Set(arg);
      fi;
      return Sum(Combinations(a), i-> (-1/2)^(Length(a)-Length(i))*A.xbasis(i));
    end;
  fi;
  A.ybasis:=function(arg) local a,res,i;
    if Length(arg)=1 then 
      if IsList(arg[1]) then a:=Set(arg[1]);
      else a:=Set(Digits(arg[1]));
      fi;
    else a:=Set(arg);
    fi;
    res:=Copy(A.zero);
    for i in Filtered([1..Length(I)],i->IsSubset(I[i],a)) do 
      res:=res+(-1)^(Length(I[i])-Length(a))*A.basis[i];
    od;
    return res;
  end;
  A.basisname:="X";
  ###################################
  A.multiplication:=function(i,j) local res,k;
    res:=[];
    for k in [1..d] do 
      if IsBound(W.solomon.mackey[i][j][k]) then 
        Add(res, [W.solomon.mackey[i][j][k]*Q.one, k]);
      fi;
    od;
    return Filtered(res, i-> i[1] <> Q.zero);
  end;
  FDAlgebraOps.StructureConstants(A);
  ###################################
  if Q.char=0 then 
    if n>10000 then InfoAlgebra("# Computing radical...\c");fi;
    k:=Concatenation(List(W.solomon.conjugacy,i->
      List(Difference(i,[i[1]]),j->A.basis[i[1]]-A.basis[j])));
    A.radical:=rec(parent:=A,
      generators:=k,
      basis:=k);
    A.radical.dimension:=Length(A.radical.basis);
    A.radical.operations:=rec(Print:=function(r) 
      Print("TwoSidedIdeal(",A,",",k,")");end);
    if Dimension(A.radical)=0 then 
      A.radical.operations.LeftTrace:=x->0*A.field.one;
      A.radical.operations.RightTrace:=x->0*A.field.one;
    else 
      subspace:=Subspace(A.vectorspace,A.EltToVector(k));
      base:=CanonicalBasis(subspace);
      A.radical.operations.LeftTrace:=x->TraceMat(List(base.vectors,v->
	  Coefficients(base,A.EltToVector(x*A.VectorToElt(v)))));
      A.radical.operations.RightTrace:=x->TraceMat(List(base.vectors,v->
	  Coefficients(base,A.EltToVector(A.VectorToElt(v)*x))));
    fi;
    A.radicalpowers:=[A.radical];
    if n>10000 then InfoAlgebra("done!\n");fi;
  fi;
  ###################################
  A.isGroup:=true;
  if Q.char=0 then 
    A.characterDegrees:=List([1..Length(W.solomon.conjugacy)], i-> 1);
  fi;
  return A;
end;

PermutationCharacterByInductiontable:=function(W,S)local t,i;
  t:=InductionTable(S,W).scalar;
  i:=PositionId(CharTable(S));
  return TransposedMat(t)[i]*CharTable(W).irreducibles;
end;

SolomonHomomorphism:=function(x) local A,i,W,T,irr,mat,WI,p;
  A:=x.domain; W:=A.group; T:=CharTable(W); irr:=T.irreducibles;
  mat:=List(irr, i-> 0);
  for i in x.coefficients do 
    if A.type="Solomon algebra" then 
      WI:=ReflectionSubgroup(W,W.solomon.subsets[i[2]]);
    elif A.type="Generalized Solomon algebra" then 
      WI:=W.generalizedsolomon.subgroups[i[2]];
    fi;
#   p:=PermutationCharacter(W,WI);
    p:=PermutationCharacterByInductiontable(W,WI); # faster for large groups
    mat:=mat+i[1]*MatScalarProducts(T,irr,[p])[1];
  od;
  return GrothendieckRing(W,A.field).VectorToElt(mat);
end;

SolomonAlgebraOps.xprimePrint:=function(r) local res,i,A,f,xprime;
  A:=r.domain;
  A.xprimebasisname:="X'";
  res:=A.EltToVector(r);
  xprime:=A.EltToVector(List(A.group.solomon.subsets,A.xprimebasis));
  res:=A.VectorToElt(res*xprime^-1);
  res:=res.coefficients;
  if IsBound(A.xprimeprint) then A.xprimeprint(r);
  else 
    f:=function(coef) 
    if coef[1]=A.field.one then 
      Print(A.xprimebasisname,"(",A.parameters[coef[2]],")");
    elif coef[1]=-A.field.one then 
      Print(" - ",A.xprimebasisname,"(",A.parameters[coef[2]],")");
    else Print(coef[1],"*",A.xprimebasisname,"(",A.parameters[coef[2]],")");
    fi;
    end;
  if Length(res)=0 then Print("0*",A.xprimebasisname,"(",
     A.parameters[A.one.coefficients[1][2]],")");
  else 
    f(res[1]);
    for i in [2..Length(res)] do 
      if res[i][1] < 0 or (res[i][1]=-A.field.one and A.field.char <> 2) 
         then Print(" - ");f([-res[i][1],res[i][2]]);
      else Print(" + ");f(res[i]);
      fi;
    od;
  fi;
  fi;
end;

ProjectionMatrix:=function(quotient) local res;
  quotient.zero:=quotient.operations.Zero(quotient);
  quotient.generators:=quotient.operations.Generators(quotient);
  Basis(quotient);
  res:=List(Basis(quotient.factorNum).vectors, 
    i-> quotient.semiEchelonBasis.operations.Coefficients(Basis(quotient),
    VectorSpaceProjection(quotient,i)));
  return res;
end;

SolomonAlgebraOps.Radical:=function(A)local B,space,quotient,mat;
  if A.field.char=0 then return A.radical;fi;
  B:=GrothendieckRing(A.group,GF(A.field.char));
  space:=B.vectorspace;
  quotient:=space/Subspace(space,B.EltToVector(Radical(B).basis));
  mat:=List(A.basis,i->B.EltToVector(SolomonHomomorphism(i)));
  mat:=mat*ProjectionMatrix(quotient);
  A.radical:=TwoSidedIdeal(A,A.VectorToElt(NullspaceMat(mat)));
  A.radicalpowers:=[A.radical];
  return A.radical;
end;

SolomonAlgebraOps.Idempotents:=function(A) local W,M,I,e,r,f,n,v,pol,i;
  W:=A.group;
  if A.field.char=0 then 
    M:=CharTable(A).irreducibles;
    if A.type="Solomon algebra" then 
      I:=List(W.solomon.conjugacy, i-> W.solomon.subsets[i[1]]);
    elif A.type="Generalized Solomon algebra" then 
      I:=List(W.generalizedsolomon.conjugacy, i-> 
        W.generalizedsolomon.signedcompositions[i[1]]);
    fi;
  else 
    #if A.type="Solomon algebra" then 
      M:=CharTable(A).invertiblematrix;
      I:=CharTable(A).parameters;
    #fi;
  fi;
  e:=TransposedMat(M^-1)*List(I,A.xbasis);
  r:=Length(e);
  f:=List([1..r],i->Sum(e{[i..r]}));
  n:=LoewyLength(A);
  v:=Indeterminate(A.field); v.name:="v";
  pol:=Sum([0..n],j->Binomial(2*n,j)*A.field.one*v^(2*n-j)*(1-v)^j);
  InfoAlgebra("# Computations to do: ",r);
  f[1]:=A.one;
  for i in [2..r] do 
    f[i]:=Value(pol,f[i-1]*f[i]*f[i-1]);
    InfoAlgebra(".",r-i,"\c");
  od;
  InfoAlgebra("\n");
  e[r]:=f[r]; for i in [1..r-1] do e[i]:=f[i]-f[i+1];od;
  return e;
end;

SolomonAlgebraOps.CartanMatrix:=function(A)local t;
  t:=FDAlgebraOps.CartanMatrix(A);
  t.columns:=List(t.rows,List);
  return t;
end;

## GeneralizedSolomonAlgebra(n[,Q])    construit   l'algebre   de   Solomon
## generalisee  de CoxeterGroup("B",n)  sur le  corps Q  (s'il est omis, la
## fonction  prend pour  Q le  corps des  rationnels). Le  resultat est une
## algebre  A  munie  en  plus  des  champs  suivants  (ici,  r est le rang
## semi-simple de W) :
##   A.group  = W
##   A.xbase  = fonction qui envoie une composition signee C de r 
##              sur l'element x_C
##   A.ybase  = fonction qui envoie une composition signee C de r 
##              sur l'element y_C
##   A.radical = radical de A : 
##        A.radical.base = base de A.radical (x_C-x_D)_{C et D conjugues}
##        A.radical.dimension = dimension de A.radical
##   A.radicalpowers = [A.radical]
##   A.injection = injection dans l'algebre de groupe (disponible 
##                 si |W| <= 400 : jusqu'a B4 donc...)
##
## De plus
##   A.base        = base des x_C
##   A.parameters  = compositions signees de r (sous forme de chaines de caracteres)
##   A.basename    = "X" (par defaut)
##   
## La fonction GeneralizedSolomonAlgebra munit W du record 
## W.generalizedsolomon contenant les champs suivants :
##   W.generalizedsolomon.mackey = constantes de structure de la multiplication
##                      x_C * x_D
##   W.generalizedsolomon.conjugacy = classes de conjugaison de compositions 
##                       signees de r reperees par leur position dans A.parameters
##   W.generalizedsolomon.domain = A
##   W.generalizedsolomon.subgroups = sous-groupes de W de la forme W_C 
##                                    ou C est une composition signee
##   W.generalizedsolomon.inclusion = inclusion des y_D dans les x_C
##   W.generalizedsolomon.signedcompositions = composition signees de r

GeneralizedSolomonAlgebraOps:=OperationsRecord("GeneralizedSolomonAlgebraOps",
  SolomonAlgebraOps);

GeneralizedSolomonAlgebraOps.CharTable:=function(A)
  local W,conj,orb,cox,l,res,TW,rows,irr,inc;
  W:=A.group;
  conj:=List(W.generalizedsolomon.conjugacy, i-> i[1]);
  orb:=List(conj, i-> W.generalizedsolomon.subgroups[i]);
  rows:=List(conj, i-> W.generalizedsolomon.signedcompositions[i]);
  #cox2:=List(conj, i-> W.generalizedsolomon.signedcompositions[i]);
  cox:=List(orb,function(S)
    if S.generators=[] then return (); else return Product(S.generators); fi;
  end);
  l:=Length(conj);
  TW:=CharTable(W);
  conj:=List(cox,i->PositionClass(W,i));
  res:=rec(irreducibles:=List([1..l], i-> PermutationCharacter(W,orb[i])),
    rows:=List(rows,Join),
    columns:=List(rows, i-> List(i, String)),
    operations:=rec(Print:=TablePrint),
    domain:=A,
    field:=A.field);
  res.irreducibles:=TransposedMat(res.irreducibles){conj};
  res.matrix:=res.irreducibles;
  #res.orb:=orb;
  res.basistraces:=List([1..Length(res.irreducibles)] , ii-> 
     List([1..Length(W.generalizedsolomon.conjugacy)], i-> 
     List(W.generalizedsolomon.conjugacy[i], j-> res.irreducibles[ii][i])));
  res.basistraces:=List(res.basistraces, Concatenation);
  if A.field.char=0 then return res;fi;
  irr:=CyclotomicModP(res.irreducibles,A.field.char);
  inc:=CollectBy([1..Length(irr)],irr);Sort(inc);
  res.rows:=List(inc, i-> res.rows[i[1]]);
  res.irreducibles:=List(inc, i-> irr[i[1]]);
  res.matrix:=res.irreducibles;
  res.basistraces:=CyclotomicModP(res.basistraces,A.field.char);
  res.parameters:=List(inc, i->W.generalizedsolomon.signedcompositions[
    W.generalizedsolomon.conjugacy[i[1]][1]]);
  res.invertiblematrix:=List(inc,i->List(inc,j->irr[i[1]][j[1]]));
  res.incidence:=inc;
  return res;
end;

GeneralizedSolomonAlgebraOps.CartanMatrix:=function(A)local t;
  t:=FDAlgebraOps.CartanMatrix(A);
  t.columns:=List(t.rows,Split);
  return t;
end;
  
GeneralizedSolomonAlgebra:=function(arg) local W,Q,A,B,r,d,i,j,k,S,T,ST,I,SC,
  montee,eta,AC,etas,montees,inclusion,n,f,mm,mminv,mackey,gmackey,w,x,ii,jj,
  kk,bips,conj,II,XX,inclusioninv,c,SCsimple,base,subspace;
  W:=CoxeterGroup("B",arg[1]);
  if Length(arg)=1 then Q:=Rationals; else Q:=arg[2]; fi;
  r:=W.semisimpleRank;
  d:=2*3^(r-1);
  A:=rec(operations:=Copy(GeneralizedSolomonAlgebraOps),
    group:=W,
    field:=Q,
    dimension:=d,
    type:="Generalized Solomon algebra");
  A.identification:=[A.type,A.group,A.field];
  A.one:=AlgebraElement(A,[[Q.one,1]]);  
  A.zero:=AlgebraElement(A,[]);
  FDAlgebraOps.underlyingspace(A);
  if not IsBound(W.generalizedsolomon) then W.generalizedsolomon:=rec();fi;
  S:=W.reflections{[2..r]};
  T:=[W.reflections[1]];
  for i in [1..r-1] do T[i+1]:=S[i]*T[i]*S[i];od;
  ST:=Concatenation(S,T);
  I:=Reversed(SignedCompositions(r));
  SC:=function(c) local i,k,l,res;
    l:=0;
    res:=[];
    for i in [1..Length(c)] do 
      k:=AbsInt(c[i]);Append(res,S{[l+1..l+k-1]});
      if c[i]>= 0 then Append(res,T{[l+1..l+k]});fi;
      l:=l+k;
    od;
    return Set(res);
  end;
  SCsimple:=function(c) local i,k,l,res;
    l:=0;
    res:=[];
    for i in [1..Length(c)] do 
      k:=AbsInt(c[i]);Append(res, S{[l+1..l+k-1]});
      if c[i]>=0 then Append(T{[l+1..l+1]});fi;
      l:=l+k;
    od;
    return Set(res);
  end;
  Sort(I,function(i,j) 
    return [Size(Subgroup(W,SC(List(i,AbsInt)))),Size(Subgroup(W,SC(i)))]>
           [Size(Subgroup(W,SC(List(j,AbsInt)))),Size(Subgroup(W,SC(j)))];
  end);
  montee:=function(w)local l; l:=CoxeterLength(W,w);
    return Set(Filtered(ST,i->CoxeterLength(W,w*i)>l));end;
  eta:=function(c) local i,k,l,res,WI,WJ;
    res:=();
    for i in Reversed([1..Length(c)]) do 
      l:=Sum([1..i],j->AbsInt(c[j]));
      k:=Sum([1..i-1],j->AbsInt(c[j]));
      if c[i]>0 then 
        WI:=ReflectionSubgroup(W,[2..l]);
        WJ:=ReflectionSubgroup(W,Difference([2..l],[k+1]));
      else 
        WI:=ReflectionSubgroup(W,[1..l]);
        WJ:=ReflectionSubgroup(W,Difference([1..l],[k+1]));
      fi;
      res:=res*LongestCoxeterElement(WI)*LongestCoxeterElement(WJ);
    od;
    return res;
  end;
  AC:=function(c) return montee(eta(c));end;
  etas:=List(I,eta);
  montees:=List(etas,montee);
  inclusion:=List(I, i-> List(I, function(j)
    if IsSubset(AC(i),SC(j)) then return 1;else return 0;fi;end));
  W.generalizedsolomon.subgroups:=List(I, i-> 
     ReflectionSubgroup(W, List(SC(i), j-> Position(Reflections(W),j))));
  inclusioninv:=inclusion^-1;
  inclusion:=TransposedMat(inclusion);
  n:=Size(W);
  ###################################
  InfoAlgebra("# Computing structure constants...\n");
  if not IsBound(W.generalizedsolomon.mackey) then 
    InfoAlgebra("# Elements to examine: ",Size(W),".\c");
    mackey:=[];
    for i in [1..d] do mackey[i]:=List([1..d],j->[]);
    od;
    f:=function(w) local i,j,k,p;
      i:=montee(w);
      i:=Position(montees,i);
      for k in [1..d] do 
        j:=montee(w^-1*etas[k]);
        j:=Position(montees,j);
        if IsBound(mackey[i][j]) then 
	  p:=mackey[i][j];
          if IsBound(p[k]) then p[k]:=p[k]+1; else p[k]:=1; fi;
        else mackey[i][j]:=[];mackey[i][j][k]:=1;
        fi;
      od;
      c:=c+1;
      if IsInt((n-c)/200) then InfoAlgebra(n-c,".\c");fi;
    end;
    c:=0;
    ForEachElement(W,f);
    gmackey:=List([1..d], i-> List([1..d], j-> []));
    mm:=List(inclusion, i-> Filtered([1..d], j-> i[j] <> 0));
    mminv:=List(inclusioninv,i->Filtered(List([1..d],j->[j,i[j]]),x->x[2]<>0));
    InfoAlgebra("# Computations to do: ",d,".\c");
    for i in [1..d] do 
      for j in [1..d] do
        for k in [1..d] do 
          x:=0;
          for ii in mm[i] do 
	    for jj in mm[j] do 
	      for kk in mminv[k] do 
		if IsBound(mackey[ii][jj][kk[1]]) then 
		  x:=x + mackey[ii][jj][kk[1]]*kk[2];
		fi;
	      od;
	    od;
          od;
          if x <> 0 then gmackey[i][j][k]:=x;fi;
        od;
      od;
      InfoAlgebra(d-i,".\c");
    od;
    InfoAlgebra(" done!\n");
    W.generalizedsolomon.mackey:=gmackey;
  fi;
  ###################################
  InfoAlgebra("# Computing orbits of parabolic subgroups...\c");
  bips:=SignedPartitions(r);
  if not IsBound(W.generalizedsolomon.conjugacy) then 
    II:=Copy(I);
    for i in [1..d] do Sort(II[i]);od;
    for i in [1..Length(bips)] do Sort(bips[i]);od;
    conj:=List(bips, i-> Filtered([1..d], j-> II[j]=i));
    Sort(conj,function(i,j) return i[1] < j[1];end);
    W.generalizedsolomon.conjugacy:=conj;
  fi;
  InfoAlgebra("done!\n");
  ###################################
  A.xbasis:=function(arg) local c;
    c:=arg; if IsList(c[1]) then c:=c[1];fi;
    if Sum(c, AbsInt) <> r then return "ERROR";fi;
    return A.basis[Position(I,c)];
  end;
  if A.field.char=0 then 
    A.xprimebasis:=function(arg) local c,aa;
    c:=arg;
    if IsList(c[1]) then c:=c[1];fi;
    if Sum(c, AbsInt) <> r then return "ERROR";fi;
    aa:=List(SCsimple(c), i-> Position(Reflections(W),i));
    return Sum(Combinations(aa), i-> (-1/2)^(Length(aa)-Length(i))*
     A.basis[Position(W.generalizedsolomon.subgroups,ReflectionSubgroup(W,i))]);
    end;
  fi;
  A.ybasis:=function(arg) local i,j,res,c;
    c:=arg;
    if IsList(c[1]) then c:=c[1];fi;
    j:=Position(I,c);
    res:=Sum([1..d], i-> inclusioninv[i][j]*A.basis[i]);
    return res;
  end;
  ###################################
  A.multiplication:=function(i,j) local res,k;
    res:=[];
    for k in [1..d] do 
      if IsBound(W.generalizedsolomon.mackey[i][j][k]) then 
        Add(res,[W.generalizedsolomon.mackey[i][j][k],k]);
      fi;
    od;
    return res;
  end;
  FDAlgebraOps.StructureConstants(A);
  ###################################
  if Q.char=0 then 
    InfoAlgebra("# Computing radical...\c");
    A.radical:=rec(parent:=A);
    A.radical.basis:=Concatenation(List(W.generalizedsolomon.conjugacy,i->
      List(Difference(i,[i[1]]),j->A.basis[i[1]]-A.basis[j])));
    A.radical.dimension:=Length(A.radical.basis);
    A.radical.generators:=A.radical.basis;
    A.radical.operations:=rec(Print:=function(r) 
      Print("TwoSidedIdeal(",A,",",A.radical.basis,")");end);
    if Dimension(A.radical)=0 then 
      A.radical.operations.LeftTrace:=x->0*A.field.one;
      A.radical.operations.RightTrace:=x->0*A.field.one;
    else 
    subspace:=Subspace(A.vectorspace,A.EltToVector(A.radical.basis));
    base:=CanonicalBasis(subspace);
    A.radical.operations.LeftTrace:=x->TraceMat(List(base.vectors,v->
        Coefficients(base,A.EltToVector(x*A.VectorToElt(v)))));
    A.radical.operations.RightTrace:=x->TraceMat(List(base.vectors,v->
        Coefficients(base,A.EltToVector(A.VectorToElt(v)*x))));
    fi;
    A.radicalpowers:=[A.radical];
    InfoAlgebra(" done!\n");
  fi;
  ####################################
  A.parameters:=List(I,Join);
  W.generalizedsolomon.signedcompositions:=I;
  A.basisname:="X";
  W.generalizedsolomon.inclusion:=inclusion;
  A.isGroup:=true;
  if Q.char=0 then 
    A.characterDegrees:=List([1..Length(W.generalizedsolomon.conjugacy)],i->1);
  fi;
  return A;
end;

QuaternionAlgebra:=function(arg) local a,b,Q,A;
  a:=arg[1];
  b:=arg[2];
  if Length(arg) = 3 then Q:=arg[3]; else Q:=Rationals;fi;
  A:=rec(operations:=Copy(FDAlgebraOps),
    dimension:=4,
    field:=Q,
    type:="Quaternion algebra",
    identification:=["Quaternion algebra",a,b,Q],
    basisname:="X",
    parameters:=[1..4]);
  FDAlgebraOps.underlyingspace(A);
  A.one:=AlgebraElement(A,[[Q.one,1]]);
  A.zero:=AlgebraElement(A,[]);
  # 1,i,j,k are the 4 basis elements
  FDAlgebraOps.StructureConstants(A,
  [ [[[Q.one,1]],[[Q.one,2]],[[Q.one,3]],[[Q.one,4]]],
    [[[Q.one,2]], [[a,1]],[[Q.one,4]],[[a,3]]],
    [[[Q.one,3]], [[-Q.one,4]],[[b,1]],[[-b,2]]],
    [[[Q.one,4]], [[-a,3]],[[b,2]],[[-a*b,1]]] ]);
  A.involution:=function(r) local res;
    res:=A.EltToVector(r);
    return A.VectorToElt([res[1],-res[2],-res[3],-res[4]]);
  end;
  A.norm:=r->r*A.involution(r);
  return A;
end;

#########################################

ZeroHeckeAlgebraOps:=OperationsRecord("ZeroHeckeAlgebraOps",FDAlgebraOps);

ZeroHeckeAlgebraOps.Print:=function(A)
  Print("ZeroHeckeAlgebra(",A.group,")");end;

# ZeroHeckeAlgebra(W[,K]) Hecke algebra of W with parameters 0,1 over K
ZeroHeckeAlgebra:=function(arg) local W,Q,A,e,d,c,ILD;
  W:=arg[1];
  if Length(arg)=1 then Q:=Rationals; else Q:=arg[2]; fi;
  e:=CoxeterElements(W);# we use they are stored by increasing CoxeterLength
  d:=Length(e);
  A:=rec(operations:=Copy(ZeroHeckeAlgebraOps),
    group:=W,
    field:=Q,
    type:="Zero Hecke algebra",
    mots:=CoxeterWords(W), # they are in same order as CoxeterElements
    dimension:=d,
    basisname:="T",
    subsets:=Combinations([1..W.semisimpleRank]));
  A.identification:=[A.type,A.group,A.field];
  A.parameters:=List(A.mots,IntListToString);
  A.one:=AlgebraElement(A,[[Q.one,1]]);
  A.zero:=AlgebraElement(A,[]);
  FDAlgebraOps.underlyingspace(A);
  A.generators:=List([1..W.semisimpleRank], i-> A.basis[1+i]);
  ILD:=W.operations.IsLeftDescending;
  A.multiplication:=function(i,j) local res,k;
    if i=1 then return [[Q.one,j]];fi;
    if j=1 then return [[Q.one,i]];fi;
    res:=e[j];
    for k in Reversed(A.mots[i]) do 
      if not ILD(W,res,k) then res:=W.reflections[k]*res;fi;
    od;
    return [[Q.one,Position(e,res)]];
  end;
  A.tbasis:=function(arg) local mot;
    if Length(arg)>1 then mot:=arg;
    elif IsList(arg[1]) then mot:=arg[1];
    elif IsInt(arg[1]) then mot:=Digits(arg[1]);
    else mot:=CoxeterWord(W,arg[1]);
    fi;
    if mot=[] then return A.basis[Position(A.mots,[])];
    else return Product(mot,i->A.basis[Position(A.mots,[i])]);
    fi;
  end;
  c:=List(CollectBy(List([1..d],i->[i,A.mots[i]]),x->Set(x[2])),
     x->List(x,x->x[1]));
  A.radical:=List(c,i->List([2..Length(i)],j->A.basis[i[1]]-A.basis[i[j]]));
  A.radical:=TwoSidedIdeal(A,Concatenation(A.radical));
  A.radicalpowers:=[A.radical];
  A.isGroup:=true;
  A.embedding:=function(g) return A.basis[Position(e,g)];end;
  return A;
end;
