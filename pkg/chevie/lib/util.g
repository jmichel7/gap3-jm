#############################################################################
##
#A  util.g              CHEVIE library          Frank Luebeck, Jean Michel
##
#Y  Copyright (C) 1992 - 2010  Lehrstuhl D fuer Mathematik, RWTH Aachen,
#Y  and   University  Paris VII.
##
##  This  file  contains  routines  of  general  utility  which  have  been
##  developed  during  the  writing  of  the Chevie package. These routines
##  should  really put into some  other part of the  GAP library, they have
##  nothing to do with Coxeter groups.
##

############################################################################
# Replace(s,s1,r1,s2,r2,....)
#  replace in string (or more generally list) s all occurences of s1 by r1
#  then of s2 by r2, etc... return the resulting string or list
#
Replace:=function(arg)local pat,rep,s,l,i;
  s:=arg[1];arg:=arg{[2..Length(arg)]};
  while Length(arg)>0 do
    pat:=arg[1];rep:=arg[2];arg:=arg{[3..Length(arg)]};
    l:=Length(pat);i:=1;
    while i<=Length(s)-l+1 do
      if s{i+[0..l-1]}=pat then 
	s:=Concatenation(s{[1..i-1]},rep,s{[i+l..Length(s)]});
	i:=i+Length(rep);
      else i:=i+1;
      fi;
    od;
  od;
  return s;
end;

#############################################################################
##
#F  Rotation(<l>, <i>) . . . . . Rotate l by i
##  
Rotation:=function(s,i)local e;e:=Length(s);i:=i mod e;
  return Concatenation(s{[i+1..e]},s{[1..i]});
end;

#############################################################################
##
#F  Rotations(<l>) . . . . . All rotations of list l
##  
Rotations:=s->List([0..Length(s)-1],i->Rotation(s,i));

#############################################################################
##
#F  Drop(<l>, <i>) . . . . . Drop i-th element from list l or
#                            elements in positions given by list i
##  
Drop:=function(l,i)
  if IsInt(i) then return l{Concatenation([1..i-1],[i+1..Length(l)])};fi;
  return l{Difference([1..Length(l)],i)};
end;

#############################################################################
##
#F  DifferenceMultiSet( <l>, <s> ) . . . . . difference of multisets
##  
DifferenceMultiSet:=function(l,s)local e,p;
  for e in s do 
    p:=Position(l,e);
    if p=false then Error(s," is not a sub-multi-set of ",l,"\n");fi;
    l:=Drop(l,p);
  od;
  return l;
end;

#############################################################################
##
#F  SymmetricDifference( <M>, <N> ) . . . . . symmetric difference of two sets
##  
SymmetricDifference:=function(x,y)
  return Difference(Union(x,y),IntersectionSet(x,y));
end;

#############################################################################
##
#F  InverseListsMap( <l> ) . . . . . . . . . . . . . . . . . . . see below
##  
##  <l>  should  be a  list   of  lists  of  positive  integers such  that
##  Concatenation(<l>) contains no element  twice. The function returns  a
##  list r, such that r[i]=[k,l] if <l>[k][l]=i.
##  
InverseListsMap:=function(l) local   res,  i; res:=[];
  for i in [1..Length(l)] do
    res{l[i]}:=List([1..Length(l[i])],j->[i,j]);
  od;
  return res;
end;

#############################################################################
##
#F  Inherit(<rec1>,<rec2>[,<fields>]) . . . . . Copy to rec1 fields of rec2
##    if <fields> specified, copies only those fields.
##
Inherit:=function(arg)local ff,f;
  if Length(arg)=3 then ff:=arg[3];else ff:=RecFields(arg[2]);fi;
  for f in ff do arg[1].(f):=arg[2].(f);od;
  return arg[1];
end;

# CollectCoefficients(rec)
# Normalizes rec.elm and rec.coeff
# where the elms are sorted and the coeffs of duplicated elms have been added.
#  could be almost (apart from zero elimination) implemented as
#  CollectCoefficients:=function(t)
#    t.elm:=Set(t.elm);t.coeff:=List(CollectBy(t.coeff,t.elm),Sum);
#  end;
CollectCoefficients:=function(t)local i,j;
  if t.elm=[] then return;fi;
  if Length(t.elm)=1 then 
    if t.coeff[1]=0*t.coeff[1] then t.elm:=[];t.coeff:=[];fi;
    return;
  fi;
  SortParallel(t.elm,t.coeff);j:=1;
  for i in [2..Length(t.elm)] do
    if t.elm[i]=t.elm[j] then t.coeff[j]:=t.coeff[j]+t.coeff[i];
    else if t.coeff[j]<>0*t.coeff[j] then j:=j+1;fi;
      t.coeff[j]:=t.coeff[i];t.elm[j]:=t.elm[i];
    fi;  
  od;
  if t.coeff[j]=0*t.coeff[j]  then j:=j-1;fi;
  if j<>Length(t.elm) then t.elm:=t.elm{[1..j]};t.coeff:=t.coeff{[1..j]};fi;
  t.elm:=Set(t.elm); # to speed up Position
end;

#############################################################################
##
#F  CartesianAt([l],i) returns Cartesian(List(l,j->[1..j]))[i]
#                      or false if i is not in [1..Product(l)]
##
CartesianAt:=function(l,pos)local res,i,j;
  res:=[];pos:=pos-1;
  for j in [Length(l),Length(l)-1.. 1] do 
    i:=pos mod l[j];
    pos:=(pos-i)/l[j];
    res[j]:=i;
  od;
  if pos<>0 then return false;else return res+1;fi;
end;

#############################################################################
##
# PositionCartesian([l1,..,ln],[i1,..,in]) returns
#  Position(Cartesian([1..l1],..,[1..ln]),[i1,..,in])
##
PositionCartesian:=function(l,ind)local res,prod,i;
  res:=1;prod:=1;
  for i in [Length(l),Length(l)-1..1] do
    if ind[i]=false then return false;fi;
    res:=res+(ind[i]-1)*prod;
    prod:=prod*l[i];
  od;
  return res;
end;

#############################################################################
##
#F  TwoTree( <m> ) . . . . . . . . . 
# Let m be a matrix whose zeroes are in symmetric positions with respect
# to the diagonal.  Let G be the graph where  vertices i,j are connected
# if m[i][j]<>0. If G is a line this  routine returns it. Else if G is a
# 3-branch tree this routine returns returns [c,b1,b2,b3] where c is the
# central node and b1,b2,b3 are  the branches oriented starting from the
# central node, listed by increasing length. Otherwise it returns 'false'
#
TwoTree:=function(m)local branch,r,line,l,p,tmp;
  branch:=function(x)
    repeat x:=PositionProperty([1..r],i->m[x][i]<>0 and not i in line);
     if x<>false then Add(line,x);fi;
    until x=false;
  end;
  r:=Length(m);line:=[1];branch(1);l:=Length(line);branch(1);
  line:=line{Concatenation([Length(line),Length(line)-1..l+1],[1..l])};
  l:=Length(line);
  if ForAny([1..l],i->ForAny([1..i-2],j->m[line[j]][line[i]]<>0)) then 
    return false;fi;
  if l=r then return line;fi;
  p:=Difference([1..r],line);
  p:=PositionProperty(line,x->ForAny(m[x]{p},i->i<>0));
  branch(line[p]);
  if Length(line)<>r then return false;fi;
  tmp:=line[p];line:=[line{[p-1,p-2..1]},line{[p+1..l]},line{[l+1..r]}];
  SortBy(line,Length);
  return Concatenation([tmp],line);
end;

#############################################################################
##
#F  Applyword(<word>,<gens>) . . . . . Evaluate word <W> on generators <gens>
##
##  Use the convention that a negative entry in word means an inverse.
ApplyWord:=function(w,gens)
  if Length(w)=0 then return gens[1]^0;
  else return Product(w,function(x)
    if x>0 then return gens[x];else return gens[-x]^-1;fi;end);
  fi;
end;

#############################################################################
##
#F  CharRepresentationWords(<rep>,<elts>) . . . . character of representation 
#F  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  on words
##  
##  <rep> is a list of matrices representing images of generators,
##  <elts> a list of words in the generators.  Returns the character on <elts>.
##  Tries to be efficient by not computing twice the same product of matrices.
## 
EvScheme:=function(elts)local l,res,i,j,p,n,w,reps,repinv;
  l:=[1..Length(elts)]; elts:=ShallowCopy(elts);SortParallel(elts,l);
  reps:=[];
  for i in [1..Length(elts)] do
    res:=[];w:=elts[i];n:=Length(w)-1;
    while Length(w)>0 do
      if n<=1 then Add(res,[w[1]]);w:=w{[2..Length(w)]};n:=Length(w);
      else p:=PositionSorted(elts,w{[1..n]});
        if p<=Length(elts) and elts[p]=w{[1..n]} then 
          Add(res,l[p]); w:=w{[n+1..Length(w)]};n:=Length(w);
        else n:=n-1;
        fi;
      fi;
    od;
    reps[l[i]]:=res;
  od;
  return reps;
end;

# EvalWords(gens,elts[,ff])
# evaluate elts on gens. If ff given do ff on each result. Otherwise
# return list of results.
EvalWords:=function(arg)local gens,elts,ff,l,res,i,f,v,reps,repinv,uses;
  gens:=arg[1];elts:=arg[2];
  l:=EvScheme(elts); InfoChevie2(Length(l),":\c");
  if IsBound(arg[3]) then uses:=List(elts,x->0);
    for v in l do for i in v do if not IsList(i)then uses[i]:=uses[i]+1;fi;od;od;
  fi;
  reps:=[];repinv:=[];
  f:=function(i)local p,res;
    if IsBound(reps[i]) then return reps[i];fi;
    res:=gens[1]^0;
    for p in l[i] do
      if IsList(p) then p:=p[1];
        if p>0 then res:=res*gens[p];
        else p:=-p;
          if not IsBound(repinv[p]) then repinv[p]:=gens[p]^-1;fi;
          res:=res*repinv[p];
        fi;
      else res:=res*f(p);
        if IsBound(arg[3]) then
           uses[p]:=uses[p]-1;if uses[p]=0 then Unbind(reps[p]);fi;
	fi;
      fi;
    od; 
    if not IsBound(arg[3]) or uses[i]>0 then reps[i]:=res;fi;
    return res;
  end;
  if not IsBound(arg[3]) then return List([1..Length(l)],f);fi;
  res:=[];
  for i in [1..Length(l)] do 
#   InfoChevie2(Position(l,i),".\c");
#   if i mod 10=0 then InfoChevie2(".\c");fi;
    res[i]:=arg[3](f(i));
#   Print(Number([1..Length(reps)],i->IsBound(reps[i])),".\c");
  od;
  InfoChevie2(Number([1..Length(l)],i->IsBound(reps[i])),"\n");
# if res<>CharRepresentationWordsO(gens,elts) then Error();fi;
  return res;
end;

CharRepresentationWords:=function(rep,elts)local F;
  if IsRec(rep) then F:=rep.F;rep:=rep.gens;fi;
  if Length(rep)=0 then 
    if IsBound(F) then return List(elts,x->TraceMat(F)); 
    else return List(elts,x->1);
    fi;
  fi;
  return EvalWords(rep,elts,function(x)
    if IsBound(F) then return TraceMat(x*F); 
    else return TraceMat(x); 
    fi;
  end);
end;

############################################################################
##
#F  InductionTable( <u>, <g> ) . . . . . . . . table of induced characters 
##  
##  'InductionTable'  returns a record describing  the decomposition of the
##  irreducible characters of the subgroup <u> induced to the group <g>.
##  
##  The result can be displayed using 'Display'.
##  
##  This function also works for Spets (Reflection Cosets)
##  
InductionTableOps:=OperationsRecord("InductionTableOps");

InductionTableOps.String:=t->SPrint(t.what,"Table(",t.u,", ",t.g,")");

InductionTableOps.Print:=function(t)Print(String(t));end;

InductionTableOps.Format:=function(t,option)
  option.rowLabels:=t.gNames(t,option);option.columnLabels:=t.uNames(t,option);
  return SPrint(t.head(t,option),"\n",FormatTable(List(t.scalar,v->List(v,
   function(a)if a=0 then return ".";else return Format(a,option);fi;end)),
       option));
end;

InductionTableOps.Display:=function(t,opt)
  opt.screenColumns:=SizeScreen()[1];
  Print(Format(t,opt));
end;

InductionTable:=function(u,g) local tu,tg,res;
  tu:=CharTable(u); tg:=CharTable(g);
  StoreFusion(tu,tg,FusionConjugacyClasses(u,g));
  res:=rec(scalar:=MatScalarProducts(tu,tu.irreducibles,
                                  Restricted(tg,tu,tg.irreducibles)));
  res.u:=u; res.g:=g; res.what:="Induction";
  res.head:=function(t,option)local n;
    n:=List([t.u,t.g],function(g)
      if IsBound(g.operations.ReflectionName) 
      then if IsBound(option.TeX) 
           then return SPrint("$",ReflectionName(g,option),"$");
           else return ReflectionName(g,option);fi;
      else return Format(g,option);
      fi;
    end);
    if IsBound(option.TeX) 
    then return SPrint(t.what," from ",n[1]," to ",n[2],"\n");
    else return SPrint(t.what," from ",n[1]," to ",n[2]);
    fi;
  end;
  res.uNames:=function(res,option)return CharNames(res.u,option);end;
  res.gNames:=function(res,option)return CharNames(res.g,option);end;
  res.operations:=InductionTableOps;
  return res;
end;

############################################################################
##
#F  AbelianGenerators( <l>) . . . . . . . minimal genators of abelian group
##
##  l is the list of elements of an abelian group G; this is a
##  (non-optimal) routine to return a minimal set of generators for G
##
AbelianGenerators:=function(l)local res,o,ord;
  ord:=function(x)local i,p,id;
    i:=1;p:=x;id:=x^0;while p<>id do p:=p*x;i:=i+1;od;
    return i;
  end;
  res:=[];
  while l<>[] do
    o:=List(l,ord);
    Add(res,l[Position(o,Maximum(o))]);
    l:=Difference(l,Elements(ApplyFunc(Group,res)));
  od;
  return res;
end;

# decompose tensor product arg{[2..Length(arg)]} (given as indices in chars)
DecomposeTensor:=function(arg)local t,ct,c;
  ct:=CharTable(arg[1]);c:=ct.irreducibles;
  t:=List(TransposedMat(c{arg{[2..Length(arg)]}}),Product);
  return MatScalarProducts(ct,c,[t])[1]; 
end;

# orbit of list of functions gens on item v
FOrbit:=function(gens,v)local res,n,new,o;
  res:=[]; new:=[v];
# Print("#I orbit length:");
  repeat
    n:=Set(Concatenation(List(gens,g->List(new,g))));
    o:=Length(res); UniteSet(res,new);
#   Print(Length(res)," \c");
    new:=Difference(n,res);
# until Length(new)=0;
  until Length(res)=o;
  return res;
end;

# orbits of list of functions gens on list v
FOrbits:=function(gens,l)local res,o;
  res:=[];
  while Length(l)>0 do
    o:=FOrbit(gens,l[1]);
    Add(res,o);
    l:=Difference(l,o);
#   Print(" #remain: ",Length(l),"\n");
  od;
  return res;
end;

# gens is a list of functions which can operate on element e
# returns minimal word w such that Composition(gens{w}) applied to e 
# satifies cond
MinimalWordProperty:=function(e,gens,cond)
  local bag,new,elements,cayleyGraph,nbLength,p,res;
    if cond(e) then return [];fi;
    elements:=[e]; nbLength:=[1]; cayleyGraph:=[[]]; bag:=Set(elements);
  InfoChevie("#I ");
  while true do
    new:=List([1+Sum(nbLength{[1..Length(nbLength)-1]})..Length(elements)],
      h->List([1..Length(gens)],function(g)local n;
       n:=gens[g](elements[h]);
       return [n,[h,g]];end));
    new:=List(CollectBy(Concatenation(new),x->x[1]),x->x[1]);
    new:=Filtered(new,x->not x[1] in bag);
    Append(cayleyGraph,List(new,x->x[2]));
    new:=List(new,x->x[1]);
    p:=PositionProperty(new,cond);
    if p<>false then InfoChevie("\n");
      res:=[]; p:=Length(elements)+p;
      while p<>1 do Add(res,cayleyGraph[p][2]);p:=cayleyGraph[p][1];od;
      return res;
    fi;
    if ForAll(new,x->x in elements) then Error("no solution");fi;
    Append(elements,new);
    UniteSet(bag,new);
#   if Length(new)>10 then 
       InfoChevie(Length(new)," \c");
#   fi;
    Add(nbLength,Length(new));
  od;
end;

#############################################################################
##
#F  PointsAndRepresentativesOrbits( <g>[, <m>] ) . . . . . . orbits of points 
#F  under permutation group <g>
##  
##  <g>  must  be  a permutation  group.   'PointAndRepresentativesOrbits'
##  returns a  list [orb,rep].  Here  orb is  a list  of the <g>-orbits on
##  [1..LargestMovedPoint(<g>)] or,  if  given, [1..<m>].  rep[i][j] is an
##  element of <g>, such that orb[i][1]^rep[i][j] = orb[i][j].
##  
##  ??? better to give list of points instead of <m> ???
PointsAndRepresentativesOrbits:=function ( arg )
  local  G, orbs, orb, max, g, gs, new, gen, Ggen, p, pnt, fst, img;
  G:=arg[1];

  if Length(arg)>1 then max:=arg[2];
  else max:=PermGroupOps.LargestMovedPoint(G);
  fi;

  Ggen:=G.generators;
  new := BlistList( [ 1 .. max ], [ 1 .. max ] );
  orbs := [  ]; gs:=[];
  fst := 1;
  while fst <> false  do
    orb := [ fst ]; g:=[()]; new[fst] := false;
    p:=1;
    while p<=Length(orb)  do
      for gen  in Ggen  do
	img := orb[p] ^ gen;
	if new[img]  then
	  Add( orb, img ); Add(g,g[p]*gen); new[img] := false;
	fi;
      od;
      p:=p+1;
    od;
    Add( orbs, orb ); Add( gs, g );
    fst := Position( new, true, fst );
  od;
  return [orbs,gs];
end;

#############################################################################
# A Dictionary data type.
# Use is:
#   d:=Dictionary() : make a new dictionary.
#   d.Get(k) : get value in d for key k.
#   d.Insert(k,v) : set in d value of key k to be v; returns v.
Dictionary:=function()local res;
  res:=rec(keys:=[],values:=[]);
  res.Get:=function(k)local p; p:=PositionSet(res.keys,k);
   if p=false then return p;else return res.values[p];fi;end;
  res.Insert:=function(k,v)local p;p:=PositionSorted(res.keys,k);
    if p<=Length(res.keys) and res.keys[p]=k then res.values[p]:=v;
    else res.keys:=Concatenation(res.keys{[1..p-1]},[k],
               res.keys{[p..Length(res.keys)]});
         res.values:=Concatenation(res.values{[1..p-1]},[v],
               res.values{[p..Length(res.values)]});
    fi;
    return v;
  end;
  res.operations:=rec(Print:=function(d)
    Print("Dictionary with ",Length(d.keys)," entries");end);
  return res;
end;

############################################################################
GetRootWarned:=[];

############################################################################
# GetRoot(x, [n, [msg]])
#  return an n-th (default 2) root of x when possible, else signals an error.
#  if msg is present and InfoChevie=Print warns about choice of root made
# after printing msg
#
#  For now, it is possible to find an nth-root of:
#  (1) a monomial of the form a.x^b when b is divisible by n and
#        we know how to find an n-th  root of a.
#  (2) a root of unity
#  (3) a rational, when n=2 or x is a perfect power
#  (4) a product of (2) by (3)
#
GetRoot:=function(arg)local x,n,msg,i,j,k,ret,d,n1,error;
  x:=arg[1];
  if Length(arg)=1 or not IsRat(arg[2]) then n:=2;arg:=arg{[2..Length(arg)]};
  else n:=arg[2];arg:=arg{[3..Length(arg)]};fi;
  if n=1 then return x;fi;
  msg:=ApplyFunc(SPrint,arg);
  if msg="" and IsBound(CHEVIE.checkroots) then msg:="checkroots";fi;
  error:=msg<>"no";
  if msg<>"" and error then
    ret:=function(res) 
      if not x in GetRootWarned then
	InfoChevie("#warning: ",msg,": ",Format(res)," chosen as ",Ordinal(n),
	  " root of ",x,"\n");
	AddSet(GetRootWarned,x);
      fi;
      return res;
    end;
  else ret:=x->x;
  fi;
  if x=1 then return ret(1);
  elif IsPolynomial(x) and x.valuation mod n=0 and Length(x.coefficients)=1 then
    return ret(Polynomial(x.baseRing,[GetRoot(x.coefficients[1],n,"no")],x.valuation/n));
  elif IsInt(x) and x<>-1 then # treat -1 as other roots of unity
    if x>0 or n mod 2=1 then i:=RootInt(x,n); if i^n=x then return ret(i);fi;
    else i:=RootInt(-x,n); if i^n=-x then return i*E(2*n);fi;
    fi;
    if n=2 then return ret(ER(x));fi;
  elif IsRat(x) and x<>-1 then
    return GetRoot(Numerator(x),n)/GetRoot(Denominator(x),n);
  elif IsCyc(x) then
    i:=AsRootOfUnity(x);
    if i<>false then
      d:=Denominator(i);j:=1;n1:=n;
      repeat k:=Gcd(n1,d);n1:=n1/k;j:=j*k; until k=1;
      return ret(E(j*d)^(Numerator(i)*GcdRepresentation(n1,d)[1]));
    fi;
    if IsInt(x^2) and n mod 2=1 then
      i:=RootInt(x^2,n); if i^n=x^2 then return ret(ER(i));fi;
    fi;
    i:=CoeffsCyc(x,NofCyc(x));
    j:=Gcd(List(i,Numerator))/Lcm(List(i,Denominator));
    if AsRootOfUnity(x/j)<>false then
      return ret(GetRoot(x/j,n,"no")*GetRoot(j,n,"no"));
    fi;
    i:=Quadratic(x*Denominator(j));
    if i<>false and n=2 then
      k:=List([1,-1],f->(i.a+f*ER(i.a^2-i.root*i.b^2))/(2*i.root*i.d));
      k:=Filtered(List(Filtered(k,IsRat),GetRoot),x->x<>false);
      k:=List(k,x->i.b/(2*x*i.d)+ER(i.root)*x);
      if Length(k)>0 then return k[1]/ER(Denominator(j));fi;
    fi;
  elif IsRec(x) and IsBound(x.operations) and IsBound(x.operations.GetRoot) then
    return ret(x.operations.GetRoot(x,n));
  fi;
  if error then Error(msg,": unable to compute ",n,"-th root of ",x,"\n");fi;
  return false;
end;

############################################################################
# EvalPolRoot(pol, x, n, p)
#  evaluate polynomial pol at p*x^(1/n)
#  
#  The point of this routine is to avoid unnecessary root extractions
#  during evaluation (e.g., if pol has no terms of odd degree and n=2,
#  then no root extraction is necessary).
#
EvalPolRoot:=function(pol,x,n,p)local P,l,r; 
  if IsPolynomial(pol) then
    if Length(pol.coefficients)=0 then return 0;fi;
    P:=Concatenation([1..pol.valuation mod n]*0,pol.coefficients);
    P:=List([0..n-1],i->Value(Polynomial(pol.baseRing,
           P{Filtered([1..Length(P)],j->(j-1) mod n=i)},
	   (pol.valuation-pol.valuation mod n)/n),x*p^n));
    pol:=Polynomial(pol.baseRing,P{[1..First([n,n-1..1],i->P[i]<>0*P[i])]},0);
    if Length(pol.coefficients)=0 then return 0;fi;
    l:=pol.valuation-1+Filtered([1..Length(pol.coefficients)],
	     i->pol.coefficients[i]<>0*pol.coefficients[i]);
    Add(l,n); r:=Gcd(l);
    pol:=Polynomial(pol.baseRing,
           pol.coefficients{[1,1+r..Length(pol.coefficients)]},
           pol.valuation/r);
    n:=n/r;
    return Value(pol,GetRoot(x,n)*p^r);
  elif IsRec(pol) and IsBound(pol.operations) 
       and IsBound(pol.operations.EvalPolRoot)
  then return pol.operations.EvalPolRoot(pol,x,n,p);
  fi;
end;

###########################################################################
#  This returns the smallest i such that 
#  stst... (i terms)= tsts... (i terms)
# (loops indefinitely if no such i exists)
BraidRelation:=function(s,t)local p1,p2,i;
  p1:=s;p2:=t;i:=1;s:=[s,t];
  while p1<>p2 do
    p1:=p1*s[1+i mod 2];p2:=p2*s[2-i mod 2];
    i:=i+1;
  od;
  return i;
end;

# CheckRelation(gens,rel[,f])
# Check  that  the  homogeneous  relation  rel[1]=rel[2] holds for gens (in
# particular  that no  left factor  homogeneous relation  holds). If given,
# call f in case of failure.
CheckRelation:=function(arg)local gens,rel,f,r,l,i,p;
  gens:=arg[1];rel:=arg[2];if Length(arg)=3 then f:=arg[3];else f:=Ignore;fi;
  p:=r->SPrint(IntListToString(r[1]),"=",IntListToString(r[2]));
  l:=gens[rel[1][1]];r:=gens[rel[2][1]];
  i:=1;
  while i<Length(rel[1]) do
    if l=r then 
      f(" relation ",p(rel)," already holds as ",p(List(rel,x->x{[1..i]})));
      return false;
    fi;
    i:=i+1;
    l:=l*gens[rel[1][i]];r:=r*gens[rel[2][i]];
  od;
  if l=r then return true;fi;
  f(" relation ",p(rel),"failed");return false;
end;

#############################################################################
# The       next      function      is       a      replacement      for
# Indeterminate(ApplyFunc(DefaultRing,V))   since   GAP   3.3   function
# DefaultRing  fails  for  DefaultRing(X(Cyclotomics),E(3))  instead  of
# returning X(Cyclotomics)

ChevieIndeterminate:=function(V)local domain,myin,elm;
  myin:=function(a,b)
    if IsRec(b) and IsBound(b.domain) then
      if a in b.domain then return true;
      else while IsBound(b.baseRing) do 
             b:=b.baseRing; if a in b then return true;fi;od;
      fi;
    fi;
    return false;
  end;
# if ForAll( V, IsInt )  then return X(Integers); elif 
  if ForAll( V, IsRat )  then return X(Rationals);
  elif ForAll( V, IsCyc )  then return X(Cyclotomics);
  else
    for elm  in V  do
      if ForAll( V, l-> myin(l,elm)) then return X(elm.domain); fi;
    od;
#   Error( "sorry, the elements of <V> lie in no common ring domain" );
    return Mvp("foo");
  fi;
end;

# get a word in W.generators for element x of W
GetWord:=function(W,x)
  return TietzeWordAbstractWord(Factorization(W,x),W.abstractGenerators);
end;

# Position of deterlinant character in chartable
PositionDet:=function(W)local ci;
  ci:=ChevieCharInfo(W);
  if not IsBound(ci.positionDet) then return false;
  else return ci.positionDet;
  fi;
end;

# permutation of Chartable of W effected by \chi->\chi\otimes\det
DetPerm:=function(W)local ct,sgn;
  ct:=CharTable(W).irreducibles;
  sgn:=ct[PositionDet(W)];
  return List(ct,x->Position(ct,Zip(x,sgn,function(a,b)return a*b;end)));
end;

##########################################################################
# Timing utility functions

# Returns time since last call
Dtime:=function()local oldtime;
  if not IsBound(CHEVIE.time) then CHEVIE.time:=Runtime();fi;
  oldtime:=CHEVIE.time; CHEVIE.time:=Runtime();
  return CHEVIE.time-oldtime;
end;

# Nicely formats result of Dtime
Stime:=function()return StringTime(Dtime());end;

Elapsed:=t->StringTime(Runtime()-t);
