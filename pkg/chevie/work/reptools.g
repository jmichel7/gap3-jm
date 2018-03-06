# Tools to build reps of W

# describe a positive integer vector as a list of indices with multiplicity
desc:=function(l)local r,i,j;r:=[];
  for i in [1..Length(l)] do for j in [1..l[i]] do Add(r,i);od;od;
  return r;
end;

# describe schur functors of Representations(W,j) when their dim is <maxdim
schurrep:=function(W,j,maxdim)local res,i,c,t,l,para;res:=[];
  t:=CharTable(W);
  for i in [2..5] do
    c:=CharTable(CoxeterGroupSymmetricGroup(i));
    para:=List(c.irredinfo,x->x.charparam[1]);
    l:=Symmetrisations(t,t.irreducibles{[j]},c);
    l:=TransposedMat([para,l]);
    l:=Filtered(l,x->x[2][1]<=maxdim);
    l:=List(l,x->[x[1],desc(DecomposeChar(W,x[2]))]);
    Append(res,l);
  od;
  return res;
end;

# find which reps are schur functors<=5 of others or tensor of 2 others
# # descreps(W[,indices of reps desired])
descreps:=function(arg)local W,t,i,j,l,k,v,dims,maxdim,v,res,c,reps;
  res:=[];W:=arg[1];t:=CharTable(W);
  dims:=List(t.irreducibles,x->x[1]);maxdim:=Maximum(dims);
  if Length(arg)=1 then reps:=[1..NrConjugacyClasses(W)];else reps:=arg[2];fi;
  for j in reps do
    l:=Filtered(schurrep(W,j,maxdim),x->Length(x[2])=1);
    for v in l do p:=v[2][1];
      if dims[p]>dims[j] then
	Print("Schur_",IntListToString(v[1]),"(",j,")=",p," \t\c");
	Add(res,[p,j,v[1]]);
      fi;
    od;
  od;
  l:=Filtered([1..NrConjugacyClasses(W)],i->dims[i]>1 and 
    dims[i]<=maxdim/Minimum(Filtered(dims,x->x<>1)));
  for j in Filtered(Combinations(l,2),p->Product(dims{p})<=maxdim) do
    v:=DecomposeTensor(W,j[1],j[2]);
    if Sum(v)=1 then  Print(j[1]," tensor ",j[2],"=",Position(v,1)," \t\c");
      Add(res,[Position(v,1),j[1],j[2]]);
    fi;
  od;
  Sort(res);return res;
end;

# the rest of the file builds a matrix for a rep which appears with multiplicity
# 1 in a tensor product
TensRep:=function(W,a,b) 
  if IsInt(a) then a:=Representations(W,a);fi;
  if IsInt(b) then b:=Representations(W,b);fi;
  return Zip(a,b,KroneckerProduct);
end;

wordsclasses:=function(W)local cl,decode,i,res,w,e;
  for w in List(ConjugacyClasses(W),Representative) do
    w:=MinimalWord(W,w);
  od;
# assumes MinimalWord(W,longest elt) has been done; computes cl
  decode:=function(i)local w;w:=[];
    while i<>1 do Add(w,W.lastmult[i][2]); i:=W.lastmult[i][1];od;
    return Reversed(w);
  end;
  res:=List([1..NrConjugacyClasses(W)],x->[]);
  for i in [1..Size(W)] do
    w:=decode(i);e:=EltWord(W,w);
    Add(res[PositionClass(W,e)],w);
    if i mod 100=1 then Print(".\c");fi;
  od;
  return res;
end;

# cl:=wordsclasses(W)
# computes sc for Repa tensor Repb
ClassSums2:=function(W,a,b,cl)local sc,e,i,elt,cnt;
  elt:=function(rep,w)local res,i;
    res:=rep[1]^0;
    for i in w do if i>0 then res:=res*rep[i];else res:=res*rep[-i]^-1;fi;od;
    return res;
  end;
  if IsInt(a) then a:=Representations(W,a);fi;
  if IsInt(b) then b:=Representations(W,b);fi;
  sc:=[];
  for i in [1..Length(cl)] do
    Print("\nStarting class ",i," Length ",Length(cl[i]),":");
    sc[i]:=IdentityMat(Length(a[1])*Length(b[1]))*0;
    cnt:=1;
    for e in cl[i] do
      sc[i]:=sc[i]+KroneckerProduct(elt(a,e),elt(b,e));
      cnt:=cnt+1;
      if (cnt mod 100)=0 then Print(".\c");fi;
    od;
  od;
  return sc;
end;

# computes sc for ExteriorPower(Repa,n)
ClassSums3:=function(W,a,n,cl)local sc,e,i,elt,cnt;
  elt:=function(rep,w)local res,i;
    res:=rep[1]^0;
    for i in w do if i>0 then res:=res*rep[i];else res:=res*rep[-i]^-1;fi;od;
    return res;
  end;
  if IsInt(a) then a:=Representations(W,a);fi;
  sc:=[];
  for i in [1..Length(cl)] do
    Print("\nStarting class ",i," Length ",Length(cl[i]),":");
    sc[i]:=IdentityMat(NrCombinations([1..Length(a[1])],n))*0;
    cnt:=1;
    for e in cl[i] do
      sc[i]:=sc[i]+ExteriorPower(elt(a,e),n);
      cnt:=cnt+1;
      if (cnt mod 100)=0 then Print(".\c");fi;
    od;
  od;
  return sc;
end;

# get rep with the following after sc:=ClassSums2
rep:=function(W,i,a,b)local M,M1,ct;
  ct:=CharTable(W).irreducibles;
  Print(Stime()," starting sums\n");
  M:=Sum([1..NrConjugacyClasses(W)],j->
   ComplexConjugate(ct[i][j])*sc[j])*ct[i][1]/Size(W);
  Print(Stime()," starting nullspace\n");
  M:=Concatenation(NullspaceMat(M-M^0),NullspaceMat(M));
  Print(Stime()," starting inverse\n");
  M1:=M^-1;
  Print(Stime()," starting conjugations\n");
  M:=List(TensRep(W,a,b),x->M*x*M1);
  Print(Stime(),"\n");
  return List(M,x->x{[1..ct[i][1]]}{[1..ct[i][1]]});end;

compact:=r->List(r,x->List(x,function(y)local l;
 l:=Filtered([1..Length(y)],i->y[i]<>0);return Concatenation(l,y{l});end));
