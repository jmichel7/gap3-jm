##########################################################################
#    Contribution to the Chevie package
#    Meinolf Geck -Jean Michel   May 1995
#    
#    This file contains functions implementing the algorithm described
#    in our paper ``On ``good'' elements in the conjugacy classes of
#     finite Coxeter groups and their eigenvalues on the irreducible
#     representations of Iwahori-Hecke algebras.''
#       Proc.\ London Math.\ Soc. (to appear).
#
#    At the end of the file, an example shows the instructions used
#    to get the character table for type E8
##########################################################################
# The following function transforms a polynomial to double list
# [degrees of non 0 coeffs, coeffs]
poltocouple:=function(p)local ind;
  if IsMvp(p) then
    return [List(p.elm,function(e)
      if Length(e.coeff)=0 then return 0; else return e.coeff[1];fi;end),
      p.coeff];
  fi;
  ind:=Filtered([1..Length(p.coefficients)],i->p.coefficients[i]<>0);
  return [p.valuation+ind-1,p.coefficients{ind}];
end;

# The following function keeps only those eigenvalues which can contribute
# to $\chi(T_w)$, i.e. those $v^a$ such that $a/2$ is integral if $\chi$
# is rational or such that $a$ is integral otherwise. The second argument
# is a list of pairs |[i,j]| such that $\chi_i$ is Galois-conjugate to
# $\chi_j$. |eig| contains for each $\chi$ just
# the list $n_{\chi,1},\ldots,n_{\chi,r}$ such that the eigenvalues are
# $v^{2 n_{\chi,1}},\ldots,v^{2n_{\chi,r}}$

reduce:=function(eig,gal)
  return List([1..Length(eig)],i->Filtered(eig[i],
	      function(x) if i in Concatenation(gal)
	                  then return IsInt(2*x);
			  else return IsInt(x);fi;end));
end;

# Now is a series of routines giving relations 1--7. The result of all these
# routines is always a list of vectors $[v_1,\ldots,v_r]$ such that each $v_i$
# is of same length $1+A$ where |A=Sum(eig,Length)|
# (the number of variables $a_{\chi,i}$) and represents a relation
#  $v_{A+1} + \sum_{j=1}^{j=A} a_j v_i=0$ (where $a_j$ is the $j$th variable
#  $a_{\chi,i}$).

# The next function returns the relations given by the value
# $\chi_(T_w)_{v\mapsto 1}$. The first argument is the index of the
# class of $w$ in the character table of $W$.
#
relsvalue:=function(index,W,eig)local i,v,res,model,ti;
  InfoChevie("using relations coming from character values\n");
  ti:=CharTable(W);
  res:=[];model:=List(eig,x->List(x,y->0));
  for i in [1..Length(eig)] do
    v:=Copy(model);v[i]:=v[i]+1;
    if Length(eig[i])>0 then
       v:=Concatenation(v);
       Add(v,-ti.irreducibles[i][index]);
       Add(res,v);
    fi;
  od;
  return res;
end;

# The next function returns the relations coming from Curtis-Alvis duality.
# The first argument is the length of $w$.

relsdual:=function(len,W,eig)local dual,model,res,i,c,v,ti,j;
  InfoChevie("using relations coming from Curtis-Alvis duality\n");
  ti:=CharTable(W);
  res:=[]; model:=List(eig,x->List(x,y->0));
  dual:=List(ti.irreducibles,x->Position(ti.irreducibles,
    List([1..Length(x)],i->x[i]*(-1)^Length(ti.classtext[i]))));
  InfoChevie("dual",dual,"\n");
  for i in [1..Length(ti.classtext)] do
   for c in [1..Length(eig[i])] do
     v:=Copy(model);
     v[i][c]:=-1;
     j:=Position(eig[dual[i]],len-eig[i][c]);
     if j<>false then v[dual[i]][j]:=v[dual[i]][j]+(-1)^len;fi;
# the addition is to handle correctly a self-dual coefficient
# the 'if' is because integrality condition may have been stronger for dual
     v:=Concatenation(v);Add(v,0);Add(res,v);
   od;
  od;
  return res;
end;

# The next function takes an argument |gal| as |reduce| does and returns
# relations expressing that $\chi_i$ and $\chi_j$ are Galois-conjugate

relsgal:=function(eig,gal)local model,c,v,res,x;
  InfoChevie("using relations coming from Galois-conjugation\n");
  res:=[]; model:=List(eig,x->List(x,y->0));
  for x in gal do
    for c in [1..Length(eig[x[1]])] do
      v:=Copy(model);
      v[x[1]][c]:=-1; v[x[2]][c]:=(-1)^(2*eig[x[1]][c]);
      v:=Concatenation(v);
      Add(v,0);
      Add(res,v);
    od;
  od;
  return res;
end;

# The next function (used in the one just after) returns the values on $T_w$
# of the exterior powers of the reflection character of the Hecke
# algebra of $W$ as polynomials in $v$.
# The first argument is the |CoxeterWord| representing $w$, the second
# argument is the Hecke algebra of $W$.
#
refchars:=function(w,H)local l,p,T,q,r;
  l:=HeckeReflectionRepresentation(H);
  if Length(w)=0 then l:=l[1]^0;else l:=Product(l{w});fi;
  q:=H.parameter[1][1];r:=Group(H).semisimpleRank;
  if IsMvp(q) then p:=Coefficients(DeterminantMat(l^0*Mvp("foo")-l),"foo");
  else p:=DeterminantMat(l^0*X(q.domain)-l).coefficients;
  fi;
  return List([1..r+1],i->p[i]*q^(Length(w)*(i-r))*(-1)^(r+1-i));
end;

# The next function returns the relations expressing the value on $T_w$ of
# the exterior powers of the reflection character.
# The first argument is the |CoxeterWord| of $w$, and the second is the
# Hecke algebra.
#
relsrefl:=function(w,H,eig)local i,j,m,c,v,res,refs,model,pos,refind,ti;
  InfoChevie("using values of the exterior powers of reflection character\n");
  ti:=CharTable(H);
  refs:=refchars(w,H);
  refind:=Reversed(ChevieCharInfo(Group(H)).extRefl);
  res:=[]; model:=List(eig,x->List(x,y->0));
  for i in [1..Length(refind)] do
    m:=poltocouple(refs[i]);
    if IsBound(H.rootParameter[1]) then m[1]:=m[1]/2;fi;
    v:=Copy(model);
    c:=IdentityMat(Length(eig[refind[i]]));
    for j in [1..Length(eig[refind[i]])] do
     v[refind[i]]:=c[j];
     pos:=Position(m[1],eig[refind[i]][j]);
     if pos=false then
       Add(res,Concatenation(Concatenation(v),[0]));
     else
       Add(res,Concatenation(Concatenation(v),[-m[2][pos]]));
       m[2][pos]:=0;
     fi;
    od;
    if Number(m[2],x->x<>0)>0 then 
 InfoChevie(" refs[",i,"]:",refs[i]," eig[",refind[i],"]:",eig[refind[i]],"\n");
     Error("value of refl. chars");fi;
  od;
  return res;
end;
    
# The next function returns relations expressing
# $\varphi(T_w)_{v\mapsto \exp(\pi/d)}=0$ for virtual characters orthogonal
# to $\Phi_d$-projectives (the kernel of the $\Phi_e$-decomposition matrix)
#  The third argument is a list of vectors $v_i$ such that $v_i$ is the
# list of coefficients of some $\varphi$ on the irreducibles.

relsd:=function(d,eig,null)local v,i,res,n,model;
  InfoChevie("using relations from the Phi_",d,"-decomposition matrix\n");
  res:=[];
  model:=List(eig,k->List(2*k,x->E(2*d)^x));
  for i in null do
     v:=Concatenation(List([1..Length(model)],k->model[k]*i[k]));
     Add(v,0);
     n:=NofCyc(v);
     v:=TransposedMat(List(v,x->CoeffsCyc(x,n)));
     v:=Filtered(v,x->Number(x,y->y<>0)>0);
     Append(res,v);
  od;
  return res;
end;

# The next function returns the relations coming from 
# $\sum_\chi \chi(T_w)d_\chi=0$. The second argument is the Hecke algebra.
#
relsdegg:=function(eig,H)local offs,v,gendeg,res,i,j,k,degs,factor;
  InfoChevie("using relations coming from generic degrees\n");
  offs:=List(eig,Length);
  offs:=Concatenation([0],List([1..Length(offs)],i->Sum(offs{[1..i]})));
  v:=SchurElements(H);
  gendeg:=List(v,x->v[1]/x);
  res:=List([1..1+4*Group(H).N],x->[1..offs[Length(offs)]]*0);
  for i in [1..Length(eig)] do
   degs:=poltocouple(gendeg[i]);
   if IsBound(H.rootParameter[1]) then degs[1]:=degs[1]/2;fi;
   for j in [1..Length(eig[i])] do
     res{2*(eig[i][j]+degs[1])}[j+offs[i]]:=degs[2];
   od;
  od;
  res:=Filtered(res,x->Number(x,y->y<>0)>0);
  for v in res do Add(v,0);od;
  return res;
end;

# To save space we apply each new set of relations as soon as possible.
# We represent our current knowledge as a record |M| with 3 fields:
# \begin{itemize}
#  \item |M.known| is the indices of the variables $a_j$ that we know already.
#  \item |M.values| is the values of the variables $a_j$ that we know already.
#  \item |M.relations| is a set of vectors representing the relations on
#     the remaining variables.
# \end{itemize}
# The next function takes a bunch of new relations and modifies |M| accordingly.
# It uses the routine  |TriangulizeMat| to find all completely known basis
# vectors resulting from |M.relations|, and then supresses the corresponding
# columns from |M.relations|, adding entries instead to |M.known| and
# |M.values|.
# It returns |true| iff at the end of the process all values are known.

apply:=function(newrels,M)local compl,i,newval,newind,NotZero;
  if Length(newrels)=0 then return false;fi;
  NotZero:=x->x<>0*x;
  compl:=Filtered([1..Length(newrels[1])-1],i->not i in M.known);
  InfoChevie(Length(newrels)," new relations\c");
  if Length(M.known)>0 then
    newrels:=List(newrels,
                x->Concatenation(x{compl},[x[Length(x)]+x{M.known}*M.values]));
  fi;
  newrels:=BaseMat(newrels);InfoChevie("(",Length(newrels)," independent) \c");
  Append(M.relations,newrels);
  M.relations:=BaseMat(M.relations);
  if Length(M.relations)>0 and 
     PositionProperty(M.relations[Length(M.relations)],y->NotZero(y))=
              Length(M.relations[1])
  then Error("contradictory relations");fi;
  i:=List(M.relations,x->Number(x{[1..Length(x)-1]},y->NotZero(y))=1);
  newval:=ListBlist(M.relations,i);
  M.relations:=ListBlist(M.relations,List(i,x->not x));
  newind:=List(newval,x->PositionProperty(x,y->NotZero(y)));
  M.relations:=List(M.relations,x->Concatenation(
	x{Filtered([1..Length(x)-1],y->not y in newind)},[x[Length(x)]]));
  Append(M.known,compl{newind});
  Append(M.values,List(newval,x->-x[Length(x)]));
  SortParallel(M.known,M.values);
  InfoChevie("known:",Length(M.known),
              " unknown:",M.total-Length(M.known)-Length(M.relations));
  if Length(M.relations)>0 then 
    InfoChevie("=",Length(M.relations[1])-1,"-",Length(M.relations));
  fi;
  Print("\n");
  return Length(M.relations)=0;
end;

##########################################################################
##
#F RelationsDecMats( <W> ) . . . . . .  compute relations derived from
## . . . . . . . . . . . . . . . . . . Phi_e-modular decomposition numbers
##
## 'RelationsDecMats' returns a list (with holes) which contains at the
## e-th position a set of vectors which yield a linear  relation among
## the Phi_e-modular reductions of the irreducible characters of Hecke(W,q).
## e runs over all possible values for which Phi_e has multiplicity at least
## two in the Poincare polynomials.
##
## Example: gap> RelationsDecMats(CoxeterGroup("E",6)); 
##
## The relations are obtained by induction of all possible linear
## independent relations from the maximal (proper) parabolic subalgebras.
##
RelationsDecMats:=function(W)
  local field,HW,r,V,W,ls,ind,e,null,rels,v,b,l,i,ct,ct1,q;
  
  field:=function(e,min) local i,q;
    i:=1;
    while true do
      q:=1+2*e*i;
      if q>min and IsPrime(q) then return [i,q]; fi;
      i:=i+1;
    od;
  end; 

  ind:=[];ct:=[];q:=X(Rationals);HW:=Hecke(W,q^2,q);r:=W.semisimpleRank;
  for i in [1..r] do
    V:=ReflectionSubgroup(W,Difference([1..r],[i]));
    InfoChevie("subalgebra of type ",ReflectionName(V),":\n");
    InfoChevie("Computing induce/restrict matrices ... \c");
    ind[i]:=TransposedMat(InductionTable(V,W).scalar);
    InfoChevie("\nComputing CharTable ... \c");
    ct[i]:=CharTable(HeckeSubAlgebra(HW,Difference([1..r],[i]))).irreducibles;
    InfoChevie("\n");
  od;
  null:=[];
  ls:=Union(List(ReflectionDegrees(W),DivisorsInt));
  for e in ls{[2..Length(ls)]} do
    v:=field(e,2*Maximum(Factors(Size(W))));
    InfoChevie("e = ",e,"; using finite field of order ",v[2]);
    rels:=[];
    for i in [1..r] do
      ct1:=Value(ct[i],E(2*e));
      for b in NullspaceMat(Value(ct[i],Z(v[2])^v[1])) do
        b:=List(b,function(x)x:=Int(x);if x>(v[2]+1)/2 then x:=x-v[2];fi;
	  return x;end);
        if b*ct1=0*b then Add(rels,b*ind[i]);
        else InfoChevie("wrong relation \n");
        fi;
      od;
    od;
    rels:=BaseMat(rels);
    InfoChevie(" got ",Length(rels)," relations\n");
    if rels<>[] then null[e]:=Copy(rels);fi;
  od;
  return null;
end;

# Now we put everything together to get a routine which returns all
# character values on $T_w$. The arguments are:
# \begin{itemize}
#  \item The index of the class of $w$ in the character table of $W$.
#  \item The Hecke algebra of $W$.
#  \item A list |gal| as explained in |reduce|.
# \end{itemize}
#
# The result is the list of character values on $T_w$ as polynomials in $v$.
#
getchar:=function(ind,H)local eig,W,M,values,e,ti,w,d;
  InfoChevie("trying element ",ind," ...\n");

# The next function returns the polynomial described by |eig| and a record
# |M| such that |M.known| is everybody and |M.relations| is empty.
# Actually we accept that |M.know| is not everybody and then give the
# arbitrary value 999 to unknown entries, which is quite distinctive
# and allows us at a glance to see what is known or not in
# a partially known character value.
#
  values:=function(M,eig)local vals,offs;
    vals:=999+[1..Sum(eig,Length)]*0;
    vals{M.known}:=M.values;
    offs:=List(eig,Length);
    offs:=Concatenation([0],List([1..Length(offs)],i->Sum(offs{[1..i]})));
    vals:=List([1..Length(eig)],i->vals{[offs[i]+1..offs[i+1]]});
    return 
      List([1..Length(vals)],function(i)
	if Length(vals[i])=0 then return 0;
	else return vals[i]*List(eig[i],function(n)
              if IsInt(n) then return H.parameter[1][1]^n;
                          else return H.rootParameter[1]^(2*n);
              fi;end);
	fi;end);
  end;

  W:=Group(H);
  w:=ChevieClassInfo(W).classtext[ind];
  d:=OrderPerm(EltWord(W,w));
  eig:=List(HeckeCharValuesGood(H,w),x->poltocouple(x)[1]/(2*d));
  eig:=reduce(eig,H.galoisperm);
  M:=rec(known:=[],values:=[],relations:=[],total:=Sum(eig,Length));
  if apply(relsvalue(ind,W,eig),M) or
     apply(relsgal(eig,H.galoisperm),M) or
     apply(relsrefl(w,H,eig),M) then return values(M,eig);fi;
  if not IsBound(H.nulldecmat) then H.nulldecmat:=RelationsDecMats(W);fi;
  for e in Filtered([1..Length(H.nulldecmat)],j->IsBound(H.nulldecmat[j])) do
    if apply(relsd(e,eig,H.nulldecmat[e]),M) then return values(M,eig);fi;
  od;
  if apply(relsdegg(eig,H),M) or
     apply(relsdual(Length(w),W,eig),M) then return values(M,eig);fi;
  InfoChevie("not enough relations\n");return values(M,eig);
end;

# return indices of the cuspidal classes in the list of classes of W 
CuspidalClasses:=function(W)local cl;
  cl:=ChevieClassInfo(W).classtext;
  return Filtered([1..Length(cl)],i->Set(cl[i])=[1..W.semisimpleRank]);
end;

# The following lines use the above routines to write all the character
# values on cuspidal classes of the Hecke algebra of $E_8$.

v:=X(Rationals); v.name:="v"; # define $v$
W:=CoxeterGroup("E",8); H:=Hecke(W,v^2,v); # define the hecke algebra.

# Define the galois-relations
H.galoisperm:=[[62,105],[63,106]]; # for E8
#H.galoisperm:=[]; # in general
#H.galoisperm:=[[59,60]]; # for E7

# And go!
l:=List(CuspidalClasses(W),i->getchar(i,H));
