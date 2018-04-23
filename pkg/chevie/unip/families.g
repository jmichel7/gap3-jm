#############################################################################
##
#A  families.g      Families of Unipotent Characters         (C)  Jean Michel
##
#Y  Copyright 1996-2017,  Univ. Paris VII
##  
##  This file contains basic functions for families of unipotent characters of
##  Spetses, and data for the most common families.
##
##  A family has at least the following fields:
##  .fourierMat  The Fourier matrix S
##  .size        Length(.fourierMat]
##  .eigenvalues The list T of Frobenius eigenvalues
##  .name        Some name describing the family (uniquely determined by S,T)
##                  if not given is ?size
##  .explanation A longer description of the family, if relevant
##                  if not given is .name
##  .charLabels  A description of the .size elements of the family
##                  if not given is List([1..size],String)
##
##  If it is in a Spets then it must have the additional fields:
##  .charNumbers Which unipotent characters are concerned
##  .special     Position of the special character
##  .coSpecial   Position of the cospecial character. If not given is .special
##  .qEigen      If .eigenvalue has fractional power of q, give it.
##
##  The following fields may be bound:
##  .signs       Signs applied to the (standard for this family) fourierMat
##  .deltas      Signs applied *on the right* to .fourierMat (in E7,E8,twisted)
##
##  Drinfeld doubles of a group G have the following additional fields:
##  .group       The group G
##  .xy          List of representatives of orbits of commuting pairs in G
##  .mellinLabels Labels for .xy
##  .mellin      Base change (x,y)->(x,\chi) fourierMat^(mellin^-1) is perm
##  .lusztig     Bound (and true) if Lusztig's variation
#############################################################################
#
#   Operations for families
#
FamilyOps:=OperationsRecord("FamilyOps");

IsFamily:=f->IsRec(f) and IsBound(f.operations) and f.operations=FamilyOps;

FamilyOps.String:=function(f)local res;
  res:=SPrint("Family(",FormatGAP(TeXStrip(f.name)));
  if IsBound(f.charNumbers) then
    PrintToString(res,",",FormatGAP(f.charNumbers));
  fi;
  return SPrint(res,")");
end;

FamilyOps.Print:=function(f)Print(String(f));end;

FamilyOps.Format:=function(f,opt)local t,collab,rowlab,res;
  opt:=ShallowCopy(opt);
  if IsBound(f.charLabels) then opt.rowLabels:=f.charLabels;
    if not IsBound(opt.TeX) then opt.rowLabels:=List(opt.rowLabels,TeXStrip);fi;
  else  opt.rowLabels:=List([1..Size(f)],x->Format(x,opt));
  fi;
  t:=[List(f.eigenvalues,x->Format(x,opt))];
  if IsBound(opt.TeX)then collab:=["\\Omega"];else collab:=["eigen"];fi;
  if IsBound(f.signs) then Add(t,f.signs);Add(collab,"signs"); fi;
  Append(t,List(f.fourierMat,x->List(x,y->Format(y,opt))));
  if Maximum(List(opt.rowLabels,Length))<=4 then Append(collab,opt.rowLabels);
  else Append(collab,List(opt.rowLabels,x->" "));
  fi;
  if IsBound(opt.TeX) then opt.rowsLabel:="\\hbox{label}";
  else opt.rowsLabel:="label";fi;
  if IsBound(f.name) then 
    if not IsBound(opt.TeX) then res:=SPrint(TeXStrip(f.name),"\n");
    else res:=SPrint("$",f.name,"$\n");fi;
  else res:="";fi;
  return Concatenation(res,
    FormatTable(TransposedMat(t),Inherit(opt,rec(columnLabels:=collab))));
end;

FamilyOps.Display:=function(x,opt)Print(Format(x,opt));end;
FamilyOps.Fourier:=f->f.fourierMat;
FamilyOps.Eigenvalues:=f->f.eigenvalues;
FamilyOps.Frobenius:=f->f.eigenvalues;

FamilyOps.\*:=function(arg)local res,a;
  arg:=List(arg,Family);
  res:=rec(
    charLabels:=List(Cartesian(List(arg,x->x.charLabels)),
      v->Join(v,"\\otimes ")),
    fourierMat:=ApplyFunc(KroneckerProduct,List(arg,x->x.fourierMat)),
    eigenvalues:=List(Cartesian(List(arg,x->x.eigenvalues)),Product),
    name:=Join(List(arg,f->f.name),"\\otimes "),
    explanation:=SPrint("Tensor(",Join(List(arg,f->f.explanation)),")"),
    size:=Product(arg,Size));
  if Length(arg)=1 then Inherit(res,arg[1]);fi;
  if ForAll(arg,f->IsBound(f.charNumbers)) then 
    res.charNumbers:=Cartesian(List(arg,x->x.charNumbers));
  fi;
  if ForAll(arg,f->IsBound(f.charLabels)) then 
    res.charLabels:=List(Cartesian(List(arg,x->x.charLabels)),Join);
  fi;
  if ForAll(arg,f->IsBound(f.special)) then 
    res.special:=PositionCartesian(List(arg,Size),List(arg,f->f.special));
    res.cospecial:=PositionCartesian(List(arg,Size),
      List(arg,function(f)if IsBound(f.cospecial) then return f.cospecial;
                                          else return f.special;fi;end));
    if res.cospecial=res.special then Unbind(res.cospecial);fi;
  fi;
  if ForAll(arg,f->IsBound(f.perm) or Size(f)=1) then 
    res.perm:=PermListList(Cartesian(List(arg,x->[1..Size(x)])),
        Cartesian(List(arg,function(x)if IsBound(x.perm) then return
          Permuted([1..Size(x)],x.perm);else return [1];fi;end)));
  fi;
  if ForAll(arg,f->IsBound(f.lusztig) or Size(f)=1) then 
    res.lusztig:=true;
  fi;
  if ForAny(arg,f->IsBound(f.qEigen)) then
    res.qEigen:=List(Cartesian(List(arg,function(f) 
      if IsBound(f.qEigen) then return f.qEigen;else return f.eigenvalues*0;fi;
      end)),Sum);
  fi;
  res.operations:=FamilyOps;
  return res;
end;

# permutes a family f by permutation or GaloisCyc
FamilyOps.\^:=function(f,p)local n;
  f:=ShallowCopy(f);
  if IsPerm(p) then
    for n in ["x","chi","charNumbers","eigenvalues","unpdeg","fakdeg",
      "mellinLabels","charLabels"] do
      if IsBound(f.(n)) then f.(n):=Permuted(f.(n),p);fi;
    od;
    for n in ["fourierMat","mellin"] do
      if IsBound(f.(n)) then f.(n):=OnMatrices(f.(n),p);fi;
    od;
    for n in ["perm","special"] do
      if IsBound(f.(n)) then f.(n):=f.(n)^p;fi;
    od;
    f.explanation:=SPrint("Permuted(",p,",",f.explanation,")");
  elif IsInt(p) then
    f.fourierMat:=List(f.fourierMat,x->List(x,y->GaloisCyc(y,p)));
    f.eigenvalues:=List(f.eigenvalues,x->GaloisCyc(x,p));
    if IsBound(f.sh) then f.sh:=List(f.sh,x->GaloisCyc(x,p));fi;
    if p=-1 then f.name:=Concatenation("\\overline ",f.name);
      f.explanation:=SPrint("ComplexConjugate(",f.explanation,")");
    else f.name:=SPrint("Gal(",p,",",f.name,")");
      f.explanation:=SPrint("GaloisCyc(",p,",",f.explanation,")");
    fi;
  fi;
  return f;
end;

OnFamily:=FamilyOps.\^; # for compatibility

FamilyOps.ComplexConjugate:=f->f^-1;

# Family(family[,charNumbers[,opt]]) where opt is a record containing
# -  signs (for complex groups)
# -  which character is special, cospecial
# -  or changing something (including the name)
Family:=function(arg)local f,ind; f:=arg[1];
  if IsString(f) then ind:=f;f:=ShallowCopy(CHEVIE.families.(f));
    if not IsBound(f.name) then f.name:=ind;fi;
  else f:=ShallowCopy(f);
  fi;
  f.operations:=FamilyOps; f.size:=Length(f.eigenvalues);
  ind:=[1..f.size];
  if not IsBound(f.charLabels) then f.charLabels:=List(ind,String);fi;
  if not IsBound(f.name) then f.name:=SPrint("?",f.size);fi;
  if not IsBound(f.explanation) then f.explanation:=f.name;fi;
  if not IsBound(f.special) then f.special:=1;fi;
  if Length(arg)=1 then return f;fi;
  f.charNumbers:=arg[2];
  if Length(arg)=2 then return f;fi;
  Inherit(f,arg[3]);
  if IsBound(f.signs) then
    f.fourierMat:=List(ind,i->f.signs[i]*List(ind,
       j->f.signs[j]*f.fourierMat[i][j]));
    if IsBound(f.perm) and -1 in Flat(f.fourierMat^2) then Unbind(f.perm);fi;
  fi;
  return f;
end;

###############################################################################
#
#    Most common families
#
# family for singletons
CHEVIE.families.C1:=Family(rec(group:="C1", name:="C_1",
  explanation:="trivial",
  charLabels:=[""], fourierMat:=[[1]], eigenvalues:=[1],
  mellin:=[[1]],mellinLabels:=[""]));
# variant occuring e.g. in unitary groups
CHEVIE.families.("C'1"):=rec(group:="C1", name:="C'_1",
  explanation:="-trivial",
  charLabels:=[""],
  fourierMat:=[[-1]],
  eigenvalues:=[-1],
  sh:=[1]);
# the Drinfeld double of Z/2
CHEVIE.families.C2:=rec(group:="C2", name:="C_2",
  explanation:="DrinfeldDouble(Z/2)",
  charLabels:=["(1,1)", "(g_2,1)", "(1,\\varepsilon)", "(g_2,\\varepsilon)"],
  fourierMat:=1/2*[[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],
  eigenvalues:=[1,1,1,-1],
  perm:=(),
  mellin:=[[1,1,0,0],[1,-1,0,0],[0,0,1,1],[0,0,1,-1]],
  mellinLabels:=["(1,1)","(1,g2)","(g2,1)","(g2,g2)"]);
# variation occuring in E7 and E8
CHEVIE.families.("C'2"):=rec(group:="C2",name:="C'_2",
  explanation:="TwistedDrinfeldDouble(Z/2)",
  charLabels:=["(1,1)",  "(1,\\varepsilon)", "(g_2,1)","(g_2,\\varepsilon)"],
  fourierMat:=1/2*[[1,1,-1,-1],[1,1,1,1],[-1,1,1,-1],[-1,1,-1,1]],
  eigenvalues:=[1,1,E(4),-E(4)],
  qEigen:=[0,0,1/2,1/2],
  perm:=(3,4),
  lusztig:=true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  cospecial:=2);
# same but non-lusztig version. Occurs in  H3 and H4 and G24,26,33,34,27
CHEVIE.families.("C'\"2"):=rec(group:="C2", name:="C'''_2",
  explanation:="TwistedDrinfeldDouble(Z/2)'",
  charLabels:=["(1,1)", "(1,\\varepsilon)", "(g_2,1)", "(g_2,\\varepsilon)"],
  fourierMat:=1/2*[[1,1,-1,-1],[1,1,1,1],[-1,1,-1,1],[-1,1,1,-1]],
  eigenvalues:=[1,1,E(4),-E(4)],
  qEigen:=[0,0,1/2,1/2],
  perm:=(3,4),
  cospecial:=2);

# n-th exterior power of the character ring of the cyclic group
# of order n, according to Michael Cuntz.
CHEVIE.families.ExtPowCyclic:=function(e,n)local g;
  g:=rec();
# S and T matrices for the cyclic group of order e
  if e mod 2 = 0 then
    g.eigenvalues:=E(24)^(e-1)*List([0..e-1],i->E(2*e)^(i*i+e*i));
  else
    g.eigenvalues:=E(24)^(e-1)*List([0..e-1],i->E(e)^QuoInt((i*i+e*i),2));
  fi;
  g.eigenvalues:=DiagonalOfMat(ExteriorPower(DiagonalMat(g.eigenvalues),n));
  g.fourierMat:=ExteriorPower(
    List([0..e-1],i->List([0..e-1],j->E(e)^(i*j)))/ER(e),n);
  g.charSymbols:=Combinations([0..e-1],n);
  g.charLabels:=List(g.charSymbols,
    s->Join(List(s,x->Format(E(e)^x,rec(TeX:=true))),"\\!\\wedge\\!"));
  if n>1 then g.name:=SPrint("R(\\BZ/",e,")^{\\wedge ",n,"}");
    g.explanation:=SPrint(Ordinal(n)," exterior power of char. ring of Z/",e);
  else g.name:=SPrint("R(\\BZ/",e,")");
    g.explanation:=SPrint("character ring of Z/",e);
  fi;
  g.eigenvalues:=g.eigenvalues/g.eigenvalues[1];
  g.special:=1;
  g.operations:=FamilyOps;
  return g;
end;

# The families for cyclic groups (big family in Z_n is X(n(n-1)/2))
# It is  the exterior square of the character ring for Z/pZ,
# However for a different normalisation than the above: the Fourier
# matrix is the opposite, and the eigenvalues are different if e<>3
CHEVIE.families.X:=function(p)local ss;
  ss:=Combinations([0..p-1],2);
  return rec(
    name:=SPrint("R_{\\BZ/",p,"}^{\\wedge 2}"),
    explanation:=SPrint("DoubleTaft(",p,")"),
    charSymbols:=ss,
    charLabels:=List(ss,s->SPrint(Format(E(p)^s[1],rec(TeX:=true)),
      "\\!\\wedge\\!",Format(E(p)^s[2],rec(TeX:=true)))),
    eigenvalues:=List(ss,s->E(p)^Product(s)),
    fourierMat:=List(ss,i->List(ss,j->(E(p)^(i*Reversed(j))-E(p)^(i*j))/p)),
    special:=1,cospecial:=p-1,operations:=FamilyOps);
end;

# Subring of family f determined by ind and the scaling factor scal
# ind is function(f,i) returns true if i in subring
SubFamily:=function(f,ind,scal,label)local res;
  ind:=Filtered([1..Length(f.eigenvalues)],i->ind(f,i));
  res:=rec(fourierMat:=f.fourierMat{ind}{ind}*scal,
   eigenvalues:=f.eigenvalues{ind},
   charLabels:=f.charLabels{ind},
   operations:=FamilyOps);
  res.name:=SPrint(f.name,"_{[",label,"]}");
  if IsBound(f.charSymbols) then res.charSymbols:=f.charSymbols{ind};fi;
  if IsBound(f.group) then res.group:=f.group;fi;
  if f.special in ind then res.special:=Position(ind,f.special);fi;
  return res;
end;

SubFamilyij:=function(f,i,j,scal)local g;
  g:=SubFamily(f,
  function(f,k)return Sum(f.charSymbols[k]) mod j=i;end,scal,Join([i,j]));
  g.explanation:=SPrint("subfamily(sum(charsymbols)mod ",j,"=",i," of ",
    f.explanation,")");
  return g;
end;

CHEVIE.families.X5:=SubFamilyij(CHEVIE.families.X(6),1,3,1-E(3));
CHEVIE.families.X5.cospecial:=5;

#  family made from QZ(3) by Gunter 22/9/99:
CHEVIE.families.Z9:=CHEVIE.families.ExtPowCyclic(9,1);
if CHEVIE.families.Z9.eigenvalues<>List([0..8],i->E(9)^(5*i^2))then Error();fi;
Inherit(CHEVIE.families.Z9, rec(
  perm:=(2,9)(3,8)(4,7)(5,6),
  qEigen:=[0,2/3,1/3,0,2/3,1/3,0,2/3,1/3]));

CHEVIE.families.Z4:=CHEVIE.families.ExtPowCyclic(4,1);
Inherit(CHEVIE.families.Z4, rec(
  fourierMat:=CHEVIE.families.Z4.fourierMat/E(4),
  eigenvalues:=CHEVIE.families.Z4.eigenvalues/CHEVIE.families.Z4.eigenvalues[2],
  special:=2,
  qEigen:=[1/2,0,1/2,0]));

CHEVIE.families.Z2:=CHEVIE.families.ExtPowCyclic(2,1);
CHEVIE.families.Z2.eigenvalues:=CHEVIE.families.Z2.eigenvalues*E(8)^3;
# the commonest families for complex groups... (also appearing as conjugates
# or in tensor products)

# Drinfeld double of Z/nZ
CHEVIE.families.QZ:=function(n)local g;
  g:=ComplexReflectionGroup(n,1,1);
  g.name:=SPrint("\\BZ/",n);
  CharTable(g).classnames:=CharNames(g,rec(TeX:=true));
  return DrinfeldDouble(g);
end;
  
# the Drinfeld double of S3, Lusztig's version
CHEVIE.families.S3:=rec(group:="S3", name:="D(S_3)",
  explanation:="Drinfeld double of S3, Lusztig's version",
  charLabels:=[ "(1,1)", "(g_2,1)", "(g_3,1)", "(1,\\rho)", "(1,\\varepsilon)",
		"(g_2,\\varepsilon)", "(g_3,\\zeta_3)", "(g_3,\\zeta_3^2)"],
  fourierMat:=[[1, 3, 2, 2,1, 3, 2, 2],[3, 3, 0, 0,-3,-3, 0, 0],
		[2, 0, 4,-2,2, 0,-2,-2],[2, 0,-2, 4, 2, 0,-2,-2],
		[1,-3, 2, 2,1,-3, 2, 2],[3,-3, 0, 0,-3, 3, 0, 0],
		[2, 0,-2,-2,2, 0, 4,-2],[2, 0,-2,-2, 2, 0,-2, 4]]/6,
  eigenvalues:=[1,1,1,1,1,-1,E(3),E(3)^2],
  perm:=(7,8),
  lusztig:=true, # does not satisfy (ST)^3=1 but (SPT)^3=1
  mellin:=[[1,0,0,2,1,0,0,0],[0,1,0,0,0,1,0,0],[0,0,1,0,0,0,1,1],[1,0,0,-1,1,0,
   0,0],[1,0,0,0,-1,0,0,0],[0,1,0,0,0,-1,0,0],[0,0,1,0,0,0,E(3),E(3)^2],
   [0,0,1,0,0,0,E(3)^2,E(3)]],
  mellinLabels:=["(1,1)","(g2,1)","(g3,1)","(1,g3)","(1,g2)","(g2,g2)",
    "(g3,g3)","(g3,g3^2)"]);

# for types B and D
# make list of families from list sym of symbols [S,T] with |S|>=|T|.
FamiliesClassical:=function(sym) local t,f,res,c,Z1,l;
# As in Lusztig Book chap. 4 we define
# Z2=Intersection(S,T), Z1=SymmetricDifference(S,T) 
# M0=elements in positions in Z1 different(mod 2) from |Z1|
# Msharp=SymmetricDifference(T-Z2,M0)
  t:=List(sym,function(ST)local f,D;ST:=FullSymbol(ST);
    f:=rec(Z1:=ApplyFunc(SymmetricDifference,ST));D:=Length(f.Z1) mod 2;
    f.Msharp:=SymmetricDifference(Difference(f.Z1,ST[2]),
				     f.Z1{[1+D,3+D..Length(f.Z1)-1]});
    if D=1 and Length(f.Msharp) mod 2<>0 then 
      f.Msharp:=Difference(f.Z1,f.Msharp);
    fi;
    f.content:=Concatenation(ST);Sort(f.content);return f;
  end);

  res:=[];
  for l in CollectBy([1..Length(t)],i->t[i].content) do
    f:=rec(content:=t[l[1]].content,charNumbers:=l);
    f.Msharp:=List(t{l},x->x.Msharp);
    if Length(l)=2 then # separate (S,S) gathered together in type D
      Add(res,rec(content:=f.content,charNumbers:=[l[2]],Msharp:=[f.Msharp[2]]));
      f:=rec(content:=f.content,charNumbers:=[l[1]],Msharp:=[f.Msharp[1]]);
    fi;
    Add(res,f);
  od;
  for f in res do
    Z1:=Filtered(f.content,x->Number(f.content,y->y=x)=1);
    f.fourierMat:=2^(-QuoInt(Length(Z1)-1,2))*
            List(f.Msharp,x->List(f.Msharp,y->(-1)^Length(Intersection(x,y))));
    f.eigenvalues:=List(f.charNumbers,x->(-1)^QuoInt(DefectSymbol(sym[x])+1,4));
# Now convert Msharp to a pair (x, \chi). Lusztig's recipe: express in the
# basis e_i=z_i+z_{i+1}, identify e1,e3,.. to elts and e2,e4,.. to chars.
# But this is not an isom to (Z/2Z)^r. Needed e.g. to replace elts with
# basis e1,e1+e3,e1+e3+e5... to get a dual basis to e2,e4,...
# In case Dn we work modulo the vector e1+e3+e5+...
    if Length(f.eigenvalues)=1 then f.charLabels:=[""];f.special:=1;
    else f.charLabels:=List(f.Msharp,function(M)local v,D,v1,v2,s;
      v:=List(Z1,z->Number(M,y->y>=z)mod 2);
      D:=Length(v);
      v1:=v{[2,4..D-(D mod 2)]};
      v2:=v{[3,5..D-1+(D mod 2)]};
      if D mod 2=1 then Add(v1,0);fi;
      # v1, v2 is coordinates in (e1,e3,e5,..) and in (e2,e4,..) basis
      v1:=List([1..Length(v2)],i->Sum(v1{[i,i+1]})mod 2);
      # coordinates in e1, e1+e3, e1+e3+e5, ...
      s:="+-";
      return ConcatenationString(s{v2+1},",",s{v1+1});end);
      f.special:=PositionProperty(f.charLabels,x->ForAll(x,y->y in "+,"));
    fi;
 #  if Length(f.eigenvalues)>1 and not 
 #    f.eigenvalues=List(f.charLabels,function(s)local s1,s2;
 #    s1:=Position(s,',');s2:=s{[s1+1..Length(s)]};s1:=s{[1..s1-1]};
 #    return (-1)^Number([1..Length(s1)],i->s1[i]='-' and s2[i]='-');end)
 #  then Error();
 #  fi;
    f.name:=IntListToString(f.content);
    f.explanation:="classical family";
    f.perm:=();
    f.size:=Length(f.charNumbers);
    f.operations:=FamilyOps;
  od;
  return res;
end;

# The big family in dihedral groups. For e=5 occurs in H3, H4
CHEVIE.families.Dihedral:=function(e)local f,c,nc,e1; e1:=QuoInt(e,2);
# the cuspidal chars are S(k,l) where 0<k<l<e-k
  nc:=Concatenation(List([1..e1-1],k->List([k+1..e-k-1],l->[k,l])));
  if e mod 2=0 then
# the principal series chars in f are:[S(0,l) with 0<l<e1]+[S(0,e1)',S(0,e1)'']
    nc:=Concatenation([[0,e1,1],[0,e1,-1]],List([1..e1-1],l->[0,l]),nc);
  else
# The principal series chars in f are:[S(0,l) with 0<l<e1+1]
    nc:=Concatenation(List([1..e1],l->[0,l]),nc);
  fi;
  c:=a->E(e)^a+E(e)^(-a);
  f:=rec(
    eigenvalues:=List(nc,s->E(e)^-Product(s{[1,2]})),
    size:=Length(nc),
    parameters:=nc,
    charLabels:=List(nc,FormatGAP),
    name:=IntListToString(Concatenation([1..e-2]*0,[1,1])),
    explanation:=SPrint("Dihedral(",e,") family"),
    operations:=FamilyOps);
  if e mod 2=0 then
    f.fourierMat:=List(nc,i->List(nc,function(j)
      if Length(i)=2 then 
	if Length(j)=2 then return (c(j*[i[2],-i[1]])-c(j*[-i[1],i[2]]))/e;
	else return  ((-1)^i[1]-(-1)^i[2])/e;
	fi;
      elif Length(i)=3 then 
	if Length(j)=2 then return ((-1)^j[1]-(-1)^j[2])/e;
	elif i=j then return (1-(-1)^e1+e)/2/e;
	else return (1-(-1)^e1-e)/2/e;
	fi;
      fi;end));
    f.special:=3;
    f.lusztig:=true;
  else
# The associated symbol to S(0,l) is s_i=[0] for i\ne 0,l and s_0=s_l=[1].
    f.fourierMat:=List(nc,i->List(nc,j->
# (-1)^Number([i[1],j[1]],x->x=0)*  This sign is in
# [Malle, "Unipotente Grade", Beispiel 6.29]
                                          (c(i*Reversed(j))-c(i*j))/e));
    f.special:=1;
  fi;
  c:=Filtered([1..Length(nc)],function(i)local p;
    p:=Position(nc,[nc[i][1],e-nc[i][2]]);return p<>false and p>i;end);
  f.perm:=Product(c,function(i)local p;
    p:=Position(nc,[nc[i][1],e-nc[i][2]]);return (i,p);end);
  return f;
end;

# size of non abelian fourier matrix of group g
NrDrinfeldDouble:=g->Sum(ConjugacyClasses(g),c->
    Length(ConjugacyClasses(Centralizer(g,Representative(c)))));

# compute the Drinfeld double of group g
# The result is a record with fields:
#   .xchi: labels of lines of the fourier matrix by
#          pairs [x,chi]: no. class of g , no character of C_g(x)
#   .charLabels: labels for xchi
#   .group is the group 
#   .special is the position of the special line (here 1 where (x,chi)=(1,1) is)
#   .eigenvalues are the eigenvalues chi(x)/chi(1)
#   .fourierMat is the Fourier matrix
#   .mellin is such that ApplyFunc(DiagonalMat,M) is the Mellin transform mat
#   then fourierMat^ApplyFunc(DiagonalMat,M)  is a permutation
#   .mellinLabels: labels for pairs [x,y] of commuting elements
# DrinfeldDouble(G,rec(lusztig:=true)) computes the matrix T
# DrinfeldDouble(G) computes the matrix S=\Delta T (see conv.tex)
DrinfeldDouble:=function(arg)local g,res,p,opt,r,lu; 
  g:=arg[1];if Length(arg)=1 then opt:=rec();else opt:=arg[2];fi;
  res:=rec(group:=g);
  lu:=IsBound(opt.lusztig) and opt.lusztig;if lu then res.lusztig:=lu;fi;
  res.classinfo:=Zip(ConjugacyClasses(g),
                     ClassNamesCharTable(CharTable(g)),
    function(c,n)local r,t;
    r:=rec(elt:=Representative(c),name:=n);
    if r.elt=g.identity then r.name:="1";fi;
    r.centralizer:=Centralizer(g,r.elt);
    r.centelms:=List(ConjugacyClasses(r.centralizer),Representative);
    t:=CharTable(r.centralizer);
    r.charNames:=CharNames(r.centralizer,rec(TeX:=true));
    r.names:=ClassNamesCharTable(t);
    r.names[Position(r.centelms,g.identity)]:="1";
    r.chars:=t.irreducibles;
    r.charNames[PositionProperty(r.chars,y->y=y*0+1)]:="1";
    r.centralizers:=t.centralizers;
    return r;end);
  res.charLabels:=Concatenation(List(res.classinfo,r->List(r.charNames,
    y->SPrint("(",r.name,",",y,")"))));
  if IsAbelian(g) then for r in res.classinfo do
    r.names:=List(r.centelms,x->First(res.classinfo,s->s.elt=x).name);
  od;fi;
  res.size:=Length(res.charLabels);
  res.eigenvalues:=Concatenation(List(res.classinfo,
    function(r)local ct;ct:=TransposedMat(r.chars);
      return Zip(ct[PositionClass(r.centralizer,r.elt)],
                 ct[PositionClass(r.centralizer,g.identity)],
       function(x,y)return x/y;end);end));
  if lu then res.name:="L";res.explanation:="Lusztig's";
  else res.name:="";res.explanation:="";fi;
  PrintToString(res.name,"D(",g,")");
  PrintToString(res.explanation,"DrinfeldDouble(",g,")");
  IsString(res.explanation);
  res.operations:=FamilyOps;
  res.mellin:=ApplyFunc(DiagonalMat,List(res.classinfo,r->
     ComplexConjugate(List(r.chars,
       x->Zip(x,r.centralizers,function(x,y)return x/y;end)))^-1));
  res.mellinLabels:=Concatenation(List(res.classinfo,
    x->List(x.names,y->SPrint("(",x.name,",",y,")"))));
  res.xy:=Concatenation(List(res.classinfo,r->List(r.centelms,y->[r.elt,y])));
  p:=Concatenation(List(res.classinfo,r->List(r.centelms,
    function(y)local p,r1;r1:=res.classinfo[PositionClass(g,y^-1)];
      return Position(res.xy,[r1.elt,r1.centelms[PositionClass(r1.centralizer,
        r.elt^RepresentativeOperation(g,y^-1,r1.elt))]]);
      end))); # Fourier permutation of the Mellin basis
  Unbind(res.classinfo);
  res.fourierMat:=IdentityMat(res.size){p}^res.mellin;
  if lu then 
    res.perm:=PermListList(ComplexConjugate(TransposedMat(res.mellin)),
      TransposedMat(res.mellin));
    res.fourierMat:=Permuted(res.fourierMat,res.perm);
  fi;
  res.special:=Position(res.charLabels,"(1,1)");
  return res;
end;

FusionAlgebra:=function(arg) local S,params,A,d,special,s,opt;
  if IsMat(arg[1]) then S:=arg[1];else S:=arg[1].fourierMat;fi;
  if Length(arg)=2 then opt:=arg[2];fi;
  if Length(arg)=1 then params:=[1..Length(S)]; else params:=arg[2].params; fi;
  if IsMat(arg[1]) then 
    if IsBound(opt) then special:=opt.special;else special:=1;fi;
  else special:=arg[1].special;fi;
  A:=rec(fourier:=S, field:=Rationals,
   operations:=OperationsRecord("FusionAlgebraOps",FDAlgebraOps),
   type:="Fusion algebra",parameters:=params,basisname:="T");
  A.identification:=[A.type,A.field,S];
  A.zero:=AlgebraElement(A,[]);
  A.one:=AlgebraElement(A,[[Rationals.one,special]]);
  d:=Length(S);
  A.dimension:=d;
  A.operations.underlyingspace(A);
  A.involution:=SignedPermListList(S,ComplexConjugate(S));
  A.Involution:=function(r)local b; b:=SignPermuted(A.basis,A.involution);
    return Sum(r.coefficients,x->x[1]*b[x[2]]);
  end;
  A.operations.CharTable:=function(A)local tbl;
    tbl:=rec(domain:=A,field:=A.field,operations:=rec(Print:=TablePrint));
    if IsBound(A.idempotents) then
      tbl.irreducibles:=List(A.idempotents,e->List(A.basis,b->
        ProportionalityCoefficient(Coefficients(b*e),Coefficients(e))));
      tbl.basistraces:=tbl.irreducibles;
      tbl.matrix:=tbl.irreducibles;
    fi;
    tbl.characterDegrees:=[1..A.dimension]*0+1;
    tbl.columns:=List(A.parameters,String);
    tbl.rows:=List(A.parameters,String);
    return tbl;
  end;
  s:=List(TransposedMat(ComplexConjugate(S)),x->x/ComplexConjugate(x[special]));
  A.structureconstants:=List([1..d],i->List([1..i],j->
    Filtered(TransposedMat([List([1..d],k->S[i][k]*S[j][k])*s,[1..d]]),
     x->x[1]<>0)));
  A.multiplication:=function(i,j)
    if i>=j then return A.structureconstants[i][j];
            else return A.structureconstants[j][i];fi;end;
  A.operations.Print:=function(A)Print(A.type," dim.",A.dimension);end;
  A.idempotents:=SignPermuted(ComplexConjugate(DiagonalMat(S[special]))*
     TransposedMat(S)*A.basis,A.involution);
  d:=Zip(CharTable(A).irreducibles,TransposedMat(S),ProportionalityCoefficient);
  if false in d then return A;fi;
  A.cDim:=d[special]^2;
  A.qDim:=List(d,x->d[special]/x);
  return A;
end;

# Michael Cuntz's formula
MQF:=function(g,N)local f,L,Y,Y1,W,S;
  f:=DrinfeldDouble(g);
  L:=[1..Length(f.xchi)];
  Y:=Filtered(L,function(i)local r;r:=f.classinfo[f.xchi[i][1]];
    return ForAll([1..r.ncl],j->r.chars[f.xchi[i][2]][j]=1 or
      not r.centelms[j] in N);end);
  Y1:=Filtered(L,function(i)local e;e:=f.classinfo[f.xchi[i][1]].elt;
    return e=g.identity or not e in N;end);
  W:=List(L,i->ER(Size(Intersection(N,Centralizer(g,f.classinfo[f.xchi[i][1]].elt))))*List(Y,function(j)local p,r;
    r:=f.classinfo[f.xchi[i][1]]; 
    if r.elt in N then
      if i=j then return 1; else return 0;fi;
    else
      if i=j then return -(r.ncl-2)/(r.ncl-1);
      elif r.elt=f.classinfo[f.xchi[j][1]].elt then return 1/(r.ncl-1);
      else return 0;
      fi;
    fi;end));
  S:=List(L,i->List(L,function(j)local p,q,r,s;
   p:=f.xchi[i];q:=f.xchi[j];
   r:=f.classinfo[p[1]]; s:=f.classinfo[q[1]];
   return Sum(Filtered(Elements(g),x->Comm(r.elt,s.elt^(x^-1))=g.identity),
     x->r.chars[p[2]][PositionClass(r.centralizer,s.elt^(x^-1))]*
        s.chars[q[2]][PositionClass(s.centralizer,r.elt^x)]/Size(r.centralizer)
	/Size(s.centralizer));end));
  if S=f.fourierMat then Print("michael's formula for fourierMat agrees\n");fi;
  return f.fourierMat{Y1}{L}*W;
end;
