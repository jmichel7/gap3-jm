#############################################################################
##
#A  matrix.g              CHEVIE library                          Jean Michel
##
#Y  Copyright (C) 2002 - 2010  University Paris VII.
##
##  This  file contains routines which augment the capabilities of GAP3 for
##  manipulating   matrices,  developed  while  writing  of  Chevie.  These
##  routines  should be in some other part  of the GAP library, as they are
##  not specific to Chevie.
##
############################################################################
##
#F  DecomposedMat( <M> ) . . . . . . . . . Find if the square matrix M
#F  admits a block decomposition.
##  
##  Define  a graph G with vertices [1..Length(M)] and with an edge between
##  i  and j if either M[i][j] or M[j][i] is not zero or false.
##  DecomposedMat  returns a list of lists I such that I[1],I[2], etc.. are
##  the  vertices  in  each  connected  component  of  G.  In  other words,
##  M{I[1]}{I[1]},M{I[2]}{I[2]},etc... are blocks of M.
## 
DecomposedMat:=function(M)local l, i, j, k, cc, cj, nz; 
  # cc[i]: in which component is i, initialized to different components
  l:=Length(M); cc:=[1..l]; nz:=x->x<>0*x and x<>false;
  for i in [1..l] do for j in [i+1..l] do
    # if new relation i~j then merge components:
    if (nz(M[i][j]) or nz(M[j][i])) and cc[i]<>cc[j] then cj:=cc[j];
      for k in [1..l] do if cc[k]=cj then cc[k]:=cc[i]; fi; od;
    fi;
  od; od;
  return Set(CollectBy([1..l],cc));
end;

############################################################################
#F BlocksMat( <M> ) . . . . . . . . . .  Find if the rectangular matrix M
#F admits a block decomposition.
##  
##  Define a bipartite graph G with vertices [1..Length(M)],
##  [1..Length(M[1])] and with an edge between i and j if either M[i][j] is
##  not  zero or false. BlocksMat  returns a list of  pairs of lists I such
##  that  [I[1][1],I[1][2]],  etc..  are  the  vertices  in  each connected
##  component of G. In other words, M{I[1][1]}{I[1][2]},
##  M{I[2][1]}{I[2][2]}, etc... are blocks of M.
## 
BlocksMat:=function(M)local  l, c, comps, p, q; comps:=[];
  for l in [1..Length(M)] do for c in [1..Length(M[1])] do
    if M[l][c]<>0*M[l][c] and M[l][c]<>false then 
      p:=PositionProperty(comps,x->l in x[1]);
      q:=PositionProperty(comps,x->c in x[2]);
      if p=false then
        if q=false then Add(comps,[[l],[c]]);
        else AddSet(comps[q][1],l);
	fi;
      elif q=false then AddSet(comps[p][2],c);
      elif p=q then AddSet(comps[p][1],l);AddSet(comps[p][2],c);
      else UniteSet(comps[p][1],comps[q][1]); UniteSet(comps[p][2],comps[q][2]);
        comps:=Drop(comps,q);
      fi;
    fi;
  od;od;
  Sort(comps);return comps;
end;

#############################################################################
##
#F  IsNormalizing( <lst>, <mat> ) . . . . . . . true if matrix <mat> lets
#F  set <lst> of vectors invariant
##  
IsNormalizing:=function(l,M)
  return Set(l*M)=Set(l);
end;

############################################################################
# EigenvaluesMat: returns eigenvalues of a square matrix over the cyclotomics
# which are roots of unity or 0. Other eigenvalues are not returned.
#
EigenvaluesMat:=function(mat)local p,res,x;
  p:=CycPol(CharacteristicPolynomial(mat));
  res:=[1..p.valuation]*0;
  for x in p.vcyc do 
    Append(res,[1..x[2]]*0+E(Denominator(x[1]))^Numerator(x[1]));
  od;
  if IsPolynomial(p.coeff) and Degree(p.coeff)=1 then
    Add(res,-p.coeff.coefficients[1]/p.coeff.coefficients[2]);
  fi;
  return res;
end;

############################################################################
#
# ExteriorPower(A,m) . . returns the m-th exterior power of square matrix m
# 
ExteriorPower:=function(A,m) local basis;
  basis:=Combinations([1..Length(A)],m); 
  return List(basis,i->List(basis,j->DeterminantMat(A{i}{j})));
end;

############################################################################
#
# SymmetricPower(A,m) . . returns the m-th symmetric power of square matrix m
# 
SymmetricPower:=function(A,m) local basis,f;
  f:=j->Product(List(Collected(j),x->x[2]),Factorial);
  basis:=UnorderedTuples([1..Length(A)],m); 
  return List(basis,i->List(basis,j->Permanent(A{i}{j})/f(i)));
end;

############################################################################
#
# SchurFunctor(M,l) . . returns the Schur functor of the matrix M
#  corresponding to partition l (e.g. if l=[n] returns the n-th symmetric
# power and if l=[1,1,1] returns the 3rd exterior power).
# 
# The current algorithm (from Littlewood) is rather inefficient 
#
SchurFunctor:=function(A,la)local n,r,m,S,M,f,rep,basis;
  n:=Sum(la);
  S:=CoxeterGroupSymmetricGroup(n);
# r:=Representations(Hecke(CoxeterGroup("A",n-1)),Position(Partitions(n),la));
  r:=Representations(S,Position(Partitions(n),la));
  rep:=function(x)x:=CoxeterWord(S,x);
    if Length(x)=0 then return r[1]^0;
    else return Product(r{x});fi;end;
  f:=j->Product(List(Collected(j),x->x[2]),Factorial);
  basis:=UnorderedTuples([1..Length(A)],n); 
  M:=Sum(Elements(S),x->KroneckerProduct(rep(x),
     List(basis,function(i)
       i:=Permuted(i,x);
       return List(basis,function(j)local p,k;p:=1;
         for k in [1..n] do p:=p*A[i[k]][j[k]];od;return p;end)/f(i);end)));
# Print(Length(M),"=>");
  M:=Filtered(M,x->x<>x*0);
  M:=TransposedMat(M);
  M:=Filtered(M,x->x<>x*0);
  m:=Set(List(M,x->Position(M,x)));
  m:=List(m,x->Filtered([1..Length(M)],i->M[i]=M[x]));
  M:=M{List(m,x->x[1])};M:=TransposedMat(M);
  M:=List(m,x->Sum(M{x}));
# Print(Length(M),"\n");
  return M;
end;

if false then
# semistandard tableaux filled with elements of [1..n]
sst:=function(S,n)local rim,res,s,t,w,p,i,u,min;
  rim:=P->Filtered([1..Length(P)],j->0<>P[j] and (j=Length(P) or  P[j+1]<P[j]));
  w:=Sum(S);if w=0 then return [List(S,x->[])];fi;
  res:=[];
  for p in rim(S) do
    s:=Copy(S);s[p]:=s[p]-1;s:=sst(s,n);
    for t in s do 
      min:=Flat(t);
      if Length(min)=0 then min:=1;else min:=Maximum(min);fi;
      if p>1 and min=t[p-1][Length(t[p])+1] then min:=min+1;fi;
      for i in [min..n] do
        u:=Copy(t);Add(u[p],i);Add(res,u);
      od;
    od;
  od;
  return Set(res);
end;

# the following is false but one hopes for a similar formula to be true
#SchurFunctor2:=function(m,p)local t;
#  t:=sst(p,Length(m));
#  return List(t,i->List(t,j->
#    Product([1..Maximum(p)],function(c)local I,J;
#      I:=List(Filtered(i,x->Length(x)>=c),y->y[c]);
#      J:=List(Filtered(j,x->Length(x)>=c),y->y[c]);
#      return DeterminantMat(m{I}{J});end)));
#end;
fi;

#############################################################################
##
#F  IndependentLines( <M>) . . . . . . 
# Returns the smallest (for lexicographic order) subset I of [1..Length(M)]
# Such that the rank of M{I} is equal to the rank of M.
##
IndependentLines:=function(M)
  if M=[] then return [];else
    return LinearIndependentColumns(TransposedMat(M));
  fi;end;

#############################################################################
##
#F  ProportionalityCoefficient( <v>, <w>) . . . . . . 'v/w'
#F  v and w are vectors. Returns scalar l such that v=l*w if such
##  exists, false otherwise.
##
ProportionalityCoefficient:=function(v,w)local i,found,coeff;
  found:=false;coeff:=0;
  for i in [1..Length(w)] do
    if not found then
      if w[i]=0 then 
        if v[i]<>0 then return false;fi;
      else found:=true;coeff:=v[i]/w[i];
      fi;
    elif v[i]<>coeff*w[i] then return false;
    fi;
  od;
  return coeff;
end;

#############################################################################
##
#F  RepresentativeDiagonalConjugations(M,N) . . diagonal matrix D such that
#   square matrices M and N are related by N=M^D or false if none exist
#
#   One has N[i][j]=M[i][j]*d[j]/d[i]
#
RepresentativeDiagonalConjugation:=function(M,N)local d,n,i,j,c;
  d:=M[1]*0;d[1]:=1;n:=Length(M);
  for i in [1..n] do
    for j in [i+1..n] do
      if M[i][j]<>0 then
	if N[i][j]=0 then return false;fi;
	if d[i]<>0 then 
	  c:=d[i]*N[i][j]/M[i][j];
	  if d[j]<>0 then if c<>d[j] then return false;fi;
	  else d[j]:=c;
	  fi;
	fi;
      fi;
    od;
  od;
  if 0 in d then return false;fi;
  if N<>M^DiagonalMat(d) then return false;fi;
  return d;
end;
    
#############################################################################
# Transporter(l1, l2) . . . . . . find matrix transporter from l1 to l2
#
# l1 and l2 should be lists of same length of square matrices of the same size.
# Transporter returns a basis of the vector space of matrices A
# such that for any i we have  A*l1[i] = l2[i]*A
# (the basis is empty if this vector space is 0)
Transporter:=function(l1,l2)local n,M,m,at,i,j,k,v;
  n:=Length(l1[1]);
  at:=function(i,j)return j+(i-1)*n;end;
  M:=[];
  for i in [1..n] do for j in [1..n] do for m in [1..Length(l1)] do
    v:=[1..n^2]*0;
    for k in [1..n] do
      v[at(i,k)]:=v[at(i,k)]+l1[m][k][j];
      v[at(k,j)]:=v[at(k,j)]-l2[m][i][k];
    od;
    Add(M,v);
  od;od;od;
  TriangulizeMat(M);
  M:=Filtered(M,x->x<>0*x);
  v:=NullspaceMat(TransposedMat(M));
  if Length(v)=0 then return false;fi;
  return List(v,w->List([1..n],i->w{[1..n]+(i-1)*n}));
end;

#############################################################################
##
#F  OnMatrices(M,p) . . . Simultaneaous action on rows and columns of a 
#                        permutation p on the square matrix M
#
OnMatrices:=function(M,p)return Permuted(List(M,y->Permuted(y,p)),p);end;

#############################################################################
##
# MatStab([g,]M[,extra]) . . . stabilizer of square matrix M for OnMatrices
#
# <g> if given should be a subgroup of Symmetricgroup(Length(M)) [default].
# returns stabilizer in g of M and of vector <extra> if given.
#
MatStab:=function(arg)local stab,M,g,r,l,I,p,s,k,n,i,j,e,blocks,extra;
  if IsGroup(arg[1]) then g:=arg[1];arg:=arg{[2..Length(arg)]};fi;
  if Length(arg)=2 then extra:=arg[2];fi;
  M:=arg[1];k:=Length(M);
  blocks:=I->CollectBy(I,function(i)local inv;
    inv:=List(I,i->[Collected(M[i]{I}),Collected(M{I}[i]),M[i][i]]);
    if IsBound(extra) then Add(inv,extra[i]);fi;
    return inv;end);
  stab:=function(I)local ind,g,p,iM;
    ind:=blocks(I);
    if Length(ind)>1 then return Concatenation(List(ind,J->stab(J)));
    elif Length(I)>1 then
      if Length(I)>7 then InfoChevie("#I Large Block:",I,"\n");fi;
      g:=MatStab(CoxeterGroupSymmetricGroup(Length(I)),M{I}{I});
      p:=MappingPermListList([1..Length(I)],I);
      return [rec(gens:=OnTuples(g.generators,p),ind:=I)];
    else return [];
    fi;
  end;
  if IsBound(g) then
    for r in blocks([1..k]) do g:=Stabilizer(g,Set(r),OnSets);od;
    s:=Concatenation(List([1..k],i->List([1..k],j->[M[i][j],k*(i-1)+j])));
    Sort(s);
    l:=[];j:=0;
    for i in [1..Length(s)] do
      if i=1 or s[i][1]<>s[i-1][1] then j:=j+1;l[j]:=[];fi;
      Add(l[j],s[i][2]);
    od;
    n:=Cartesian([1..k],[1..k]);
    e:=Group(List(g.generators,p->PermListList(n,List(n,x->OnTuples(x,p)))),());
    for s in l do e:=Stabilizer(e,s,OnSets);od;
    return Group(List(e.generators,p->PermList(List([1..k],i->n[i^p][2]))),());
  fi;
  g:=Group(());I:=[];
  for r in stab([1..k]) do
    Append(I,r.ind);p:=MappingPermListList(I,[1..Length(I)]);
    g:=Group(Concatenation(g.generators,OnTuples(r.gens,p)),());
    g:=MatStab(g,M{I}{I});
  od;
  return Group(List(g.generators,x->x^(p^-1)),());
end;

##########################################################################
##
#F DistHelpedRepresentativeOperation(G,x,y,opr,dist)   .   .   Heuristic
## version of 'RepresentativeOperation' when a distance is known between
## x and y.
#
# In the case where x and y live in a space on which we have a distance,
# we  can try to  get from x  to y in  G by multiplying by the generator
# which  brings us closer.  Coupled with the  trick of applying a random
# perturbation  if we fall in  a hole which is  not our goal, this often
# works surprisingly well. 'dist(x,y)' should return rational numbers.
#
DistHelpedRepresentativeOperation:=function(G,x,y,opr,dist)local p,d,prev,cv,x1;
  InfoChevie("#I  group:",Size(G)," too big - trying random walk\n");
  # best generator in G towards x=y
  cv:=function(x)local minimum,mp,i,nn;
    mp:=();
    minimum:=dist(x,y);
    for i in G.generators do
      nn:=dist(opr(x,i),y);
      if nn<minimum then minimum:=nn; mp:=i;fi;
    od;
    if mp<>() then InfoChevie(Position(G.generators,mp),"->",minimum," \c");fi;
    return mp;
  end;

  p:=();
  InfoChevie(dist(x,y)," \c");
  while true do
    x1:=opr(x,p);
    prev:=dist(x1,y);
    if prev=0 then InfoChevie("\n"); return p;fi;
    p:=p*cv(x1);
    d:=dist(opr(x,p),y);
    if d=prev then 
      # Print(Format(opr(x,p)-y),"\n");
      InfoChevie("\n#I stalled -- restarting at a random element of G\n");
      p:=p*Random(G);
    fi;
  od;
end;

#############################################################################
##
# PermMatMat(M, N[, m ,n]) . . . . p such that OnMatrices(M,p)=N
# If  in  addition  the  vectors  m  and  n  are  given,  p  should satisfy
# Permuted(m,p)=n.
#
# Efficient version of 
#   RepresentativeOperation(SymmetricGroup(Length(M)),M,N,OnMatrices).
#
# Test CartanMat("D",12) with p=( 1, 5, 2, 8,12, 4, 7)( 3, 9,11, 6)
PermMatMat:=function(arg)local ind,l,I,J,r,p,e,g,s,h,sg,M,N,m,n;
  M:=arg[1];N:=arg[2];
  if Length(arg)>2 then m:=arg[3];n:=arg[4];
  elif M=N then return ();fi;
  if Length(M)<>Length(N) then 
    InfoChevie("# matrices do not have same dimensions");return false;
  fi;
  if Length(M[1])<>Length(M) or Length(N[1])<>Length(N) then 
    Error("matrices are not square");return false;
  fi;
  sg:=n->Group(Concatenation(List([1..n-1],i->List([i+1..n],j->(i,j)))),());
  ind:=function(I,J)local iM,iN,p;
    iM:=List(I,i->[Collected(M[i]{I}),Collected(M{I}[i]),M[i][i]]);
    iN:=List(J,i->[Collected(N[i]{J}),Collected(N{J}[i]),N[i][i]]);
    if IsBound(m) then 
      iM:=Zip(iM,m{I},function(x,y)Add(x,y);return x;end);
      iN:=Zip(iN,n{J},function(x,y)Add(x,y);return x;end);
    fi;
    if Set(iM)<>Set(iN) then return false;fi;
    iM:=CollectBy(I,iM); iN:=CollectBy(J,iN);
    if List(iM,Length)<>List(iN,Length) then return false;fi;
    p:=Zip(iM,iN,function(I,J)local g;
      if Length(I)>7 then InfoChevie("#I  large block:",Length(I),"\n");
        if Length(iM)=1 then
          p:=DistHelpedRepresentativeOperation(sg(Length(I)),M{I}{I},N{J}{J},
	   OnMatrices,function(M,N)return Sum(M-N,x->Number(x,y->y<>0));end);
	elif IsBound(m) then p:=PermMatMat(M{I}{I},N{J}{J},m{I},n{J});
	else p:=PermMatMat(M{I}{I},N{J}{J});
	fi;
      else p:=RepresentativeOperation(sg(Length(I)),M{I}{I},N{J}{J},OnMatrices);
      fi;
      if p=false then return false;fi;
      I:=Permuted(I,p);
      g:=MatStab(M{I}{I});
      p:=MappingPermListList([1..Length(I)],I);
      g:=Group(OnTuples(g.generators,p),());
      return [I,J,g];end);
 #  Print("p=",p,"\n");
    if false in p then return false; else return p;fi;
  end;
  l:=ind([1..Length(M)],[1..Length(N)]);
  if l=false then return false;fi;
  I:=[];J:=[];g:=Group(());
  for r in l do
    Append(I,r[1]);Append(J,r[2]);
    s:=Length(r[1]);
    g:=Group(Concatenation(g.generators,r[3].generators),());
    p:=MappingPermListList(I,[1..Length(I)]);
    h:=Group(OnTuples(g.generators,p),());
    if M{I}{I}<>N{J}{J} then
      InfoChevie("# I=",Length(I)," stab=",Size(g),"\n");
      e:=RepresentativeOperation(h,M{I}{I},N{J}{J},OnMatrices);
      if e=false then return false;
      else I:=Permuted(I,e);
      fi;
    fi;
    h:=Stabilizer(h,M{I}{I},OnMatrices);
    g:=Group(OnTuples(h.generators,p^-1),());
  od;
  return MappingPermListList(I,J);
end;

PermMatMat2:=function(m1,m2)local mm,rg,rg1,iv,inv,i,g,p,perm,best,s,dist;
  dist:=function(m1,m2)return Sum([1..Length(m1)],
       i->Number([1..Length(m1[1])],j->m1[i][j]<>m2[i][j]));
  end;
  if Length(m1)<>Length(m2) then 
    Error("matrices do not have same dimensions");return false;
  fi;
  if Length(m1)=0 then return [(),()];fi;
  if Length(m1[1])<>Length(m1) or Length(m2[1])<>Length(m2) then 
    Error("matrices are not square");return false;
  fi;
  iv:=function(gg,x)return List(gg,g->Collected(x{g}));end;
  mm:=[m1,m2];
  InfoChevie(dist(mm[1],mm[2]),"\c");
  perm:=[(),()];
  rg1:=[[1..Length(m1)]];
  repeat
    rg:=rg1;rg1:=[]; inv:=[];
    for g in rg do
      for i in [1,2] do
        inv[i]:=List(mm[i]{g},x->iv(rg,x));
	p:=MappingPermListList(Flat(CollectBy(g,inv[i])),g);
	perm[i]:=perm[i]*p; mm[i]:=OnMatrices(mm[i],p);
        Sort(inv[i]);
      od;
      if inv[1]<>inv[2] then return false;fi;
      Append(rg1,CollectBy(g,inv[1]));
    od;
    InfoChevie("=>",dist(mm[1],mm[2]),"\c");
  until rg=rg1;
  best:=function(l)local d,G,m,e;
    d:=dist(mm[1],mm[2]);
    G:=Elements(Group(List([1..Length(l)-1],i->(l[i],l[i+1])),()));
    m:=Minimum(List(G,x->dist(OnMatrices(mm[1],x),mm[2])));
    if m<d then 
      e:=First(G,x->dist(OnMatrices(mm[1],x),mm[2])=m);
      mm[1]:=OnMatrices(mm[1],e);
      perm[1]:=perm[1]*e;
      InfoChevie("->",m,"\c");
      return true;
    fi;
    return false;
  end;
  repeat s:=false; for g in rg do s:=s or best(g); od; until not s;
  InfoChevie("\n");
  if mm[1]<>mm[2] then Error("PermMatMat failed");fi;
  p:=perm[1]/perm[2];
  if OnMatrices(m1,p)<>m2 then Error("theory");fi;
  return p;
end;

############################################################################
# PermutedByCols(m) . . Operation of a permutation on the columns of a  matrix
#
PermutedByCols:=function(m,p)return List(m,x->Permuted(x,p));end;

############################################################################
# RepresentativeRowColPermutation(m1,m2) . . . . . . . . . . . . whether
#  matrix m1 is conjugate to matrix m2 by row/col permutations
#
#  m1 and m2 should be rectangular matrices of the same dimensions.
#  The function returns a pair of permutations [p1,p2] such that
#  PermutedByCols(Permuted(m1,p[1]),p[2])=Permuted(PermutedByCols(m1,p2),p1)=m2
#  if such permutations exist, and false otherwise.
#
RepresentativeRowColPermutation:=function(m1,m2)
  local mm,rg,rg1,cg,cg1,inv,i,g,p,nr,nc,rperm,cperm,best,s,dist;
  nr:=Length(m1);nc:=Length(m1[1]);
  if nr<>Length(m2) then Error("not same dimensions");fi;
  if nr=0 then return [(),()];fi;
  if nc<>Length(m2[1]) then Error("not same dimensions");fi;
  dist:=function(arg)
    if Length(arg)=2 then 
        return Sum([1..nr],i->Number([1..nc],j->arg[1][i][j]<>arg[2][i][j]));
    elif arg[3]=Permuted then 
         return Sum(arg[4],i->Number([1..nc],j->arg[1][i][j]<>arg[2][i][j]));
    else return Sum(arg[4],j->Number([1..nr],i->arg[1][i][j]<>arg[2][i][j]));
    fi;
  end;
  mm:=[m1,m2];
  InfoChevie("# ",dist(m1,m2),"\c");
  rperm:=[(),()];cperm:=[(),()];
  rg1:=[[1..nr]];cg1:=[[1..nc]];
  repeat
    rg:=rg1;rg1:=[]; cg:=cg1;cg1:=[]; inv:=[];
    for g in rg do
      for i in [1,2] do
        inv[i]:=List(mm[i]{g},x->List(cg,g->Collected(x{g})));
	p:=MappingPermListList(Flat(CollectBy(g,inv[i])),g);
	rperm[i]:=rperm[i]*p; mm[i]:=Permuted(mm[i],p);
        Sort(inv[i]);
      od;
      if inv[1]<>inv[2] then return false;fi;
      Append(rg1,CollectBy(g,inv[1]));
    od;
    for g in cg do
      for i in [1,2] do
        inv[i]:=List(TransposedMat(mm[i]){g},x->List(rg,g->Collected(x{g})));
	p:=MappingPermListList(Flat(CollectBy(g,inv[i])),g);
	cperm[i]:=cperm[i]*p; mm[i]:=PermutedByCols(mm[i],p);
        Sort(inv[i]);
      od;
      if inv[1]<>inv[2] then return false;fi;
      Append(cg1,CollectBy(g,inv[1]));
    od;
    InfoChevie("=>",dist(mm[1],mm[2]),"\c");
  until rg=rg1 and cg=cg1;
  best:=function(l,opr)local d,G,m,e;
    if Length(l)=1 then return false;fi;
    d:=dist(mm[1],mm[2],opr,l);
    G:=Elements(Group(List([1..Length(l)-1],i->(l[i],l[i+1])),()));
    for e in G do
      m:=dist(opr(mm[1],e),mm[2],opr,l);
      if m<d then
	InfoChevie(m-d,"\c");
	if opr=Permuted then rperm[1]:=rperm[1]*e;else cperm[1]:=cperm[1]*e;fi;
	mm[1]:=opr(mm[1],e);return true;
      fi;
    od;
    return false;
  end;
  repeat s:=false;
    for g in rg do s:=s or best(g,Permuted); od;
    for g in cg do s:=s or best(g,PermutedByCols); od;
  until not s;
  InfoChevie("\n");
  if mm[1]<>mm[2] then Error("RepresentativeRowColOperation failed");fi;
  p:=[rperm[1]/rperm[2],cperm[1]/cperm[2]];
# if PermutedByCols(Permuted(m1,p[1]),p[2])<>m2 then Error("theory");fi;
  return p;
end;

#############################################################################
##
# CoFactors(M) . . . .  Cofactors of the square matrix M.
# 
# This returns a matrix such that CoFactors(M)*M=DeterminantMat(M)*M^0
#
CoFactors:=function(m)local v;
  if Length(m)=1 then return [[1]];fi;
  v:=[1..Length(m)];
  return TransposedMat(List(v,i->List(v,j->(-1)^(i+j)*
    DeterminantMat(m{Filtered(v,x->x<>i)}{Filtered(v,x->x<>j)}))));
end;

#############################################################################
##
# BigCellDecomposition(M [, b]) . . . .  Decompose in the big Bruhat cell
# M should be a square matrix such that the principal minors which are
# union of blocks should be non-zero.
# The  function decomposes  M as  a product  P*L*tP where  P is lower block
# unitriangular  and tp  upper block-unitriangular  (with identity diagonal
# blocks)  and L block-diagonal according to the  block structure b; b is a
# list  of lists of  union [1..Length(M)]. If  not given, the trivial block
# structure [[1],..,[Length(M)]] is assumed.
# If  M is  symmetric then  tP=TransposedMat(P) and  the result is the pair
# [tP,L]. else the result is [P,L,tP]
#
BigCellDecomposition:=function(arg)local M,b,P,tP,L,j,i,c,cb,db;
  M:=arg[1];
  if Length(arg)=1 then b:=List([1..Length(M)],i->[i]);else b:=arg[2];fi;
  L:=M^0;P:=M^0;
  if M=TransposedMat(M) then
    for j in [1..Length(b)] do
      L{b[j]}{b[j]}:=M{b[j]}{b[j]}-Sum([1..j-1],k->P{b[j]}{b[k]}*L{b[k]}{b[k]}*
	TransposedMat(P{b[j]}{b[k]}));
      cb:=CoFactors(L{b[j]}{b[j]});
      db:=DeterminantMat(L{b[j]}{b[j]});
      for i in [j+1..Length(b)] do
	P{b[i]}{b[j]}:=(M{b[i]}{b[j]}-Sum([1..j-1],k->P{b[i]}{b[k]}*
	 L{b[k]}{b[k]}*TransposedMat(P{b[j]}{b[k]})))*cb/db;
      od;
    od;
    return [TransposedMat(P),L];
  fi;
  tP:=M^0;
  for j in [1..Length(b)] do
    L{b[j]}{b[j]}:=M{b[j]}{b[j]}-Sum([1..j-1],k->P{b[j]}{b[k]}*L{b[k]}{b[k]}*
      tP{b[k]}{b[j]});
    cb:=CoFactors(L{b[j]}{b[j]});
    db:=DeterminantMat(L{b[j]}{b[j]});
    for i in [j+1..Length(b)] do
      P{b[i]}{b[j]}:=(M{b[i]}{b[j]}-Sum([1..j-1],k->P{b[i]}{b[k]}*
       L{b[k]}{b[k]}*tP{b[k]}{b[j]}))*cb/db;
      tP{b[j]}{b[i]}:=cb*(M{b[j]}{b[i]}-Sum([1..j-1],k->P{b[j]}{b[k]}*
       L{b[k]}{b[k]}*tP{b[k]}{b[i]}))/db;
    od;
  od;
  return [P,L,tP];
end;
