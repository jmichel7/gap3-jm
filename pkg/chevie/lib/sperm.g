#############################################################################
##
#A  sperm.g          CHEVIE library                               Jean Michel
##
#Y  Copyright (C) 2017 - 2017   University  Paris VII.
##
##  This  file  contains  routines  to work with signed permutations.
##  A signed permutation is a compact way to work with permutations
##  of [-n,..,-1,1,..n]. Such a permutation is represented by a
##  list of length n of the integers in [1..n] with signs.
##

## Apply signed permutation sp to list l
SignPermuted:=function(l,sp)
  return List(sp,i->l[AbsInt(i)]*SignInt(i));
end;

# SignedPerm to matrix m such that m*l=SignPermuted(sp,l)
SignedPermutationMat:=function(sp)local n,res,i;
  n:=Length(sp);
  res:=NullMat(n);
  for i in [1..n] do res[i][AbsInt(sp[i])]:=SignInt(sp[i]);od;
  return res;
end;

# Find if exists signed perm which permutes list a to list b
SignedPermListList:=function(a,b)local a1,b1,p,i,res;
  a1:=List(a,x->[x,-x]);for p in a1 do Sort(p);od;
  b1:=List(b,x->[x,-x]);for p in b1 do Sort(p);od;
  p:=PermListList(a1,b1);
  if p=false then return false;fi;
  res:=Permuted([1..Length(a)],p);
  for i in [1..Length(a)] do
    if b[i^p]<>a[i] then res[i]:=-res[i];fi;
  od;
  return res;
end;
