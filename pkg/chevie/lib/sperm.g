#############################################################################
##
#A  sperm.g          CHEVIE library                               Jean Michel
##
#Y  Copyright (C) 2017 - 2018   University  Paris VII.
##
##  This  file contains routines to work with signed permutations. A signed
##  permutation is a compact way to work with permutations of [-n..-1,1..n]
##  which preserve the pairs [-i,i]. Such a permutation is represented by a
##  list  l  of  length  n  of  the  integers  in  [1..n] with signs, which
##  represents  the permutation sending i to l[i],  or by an element of the
##  hyperoctaedral group, acting on [n+1..2*n] instead of [-n..-1].
##

## Apply signed permutation sp to list l
SignPermuted:=function(l,sp)local res,i;
  if IsPerm(sp) then 
    return Permuted(Concatenation(l,Reversed(-l)),sp){[1..Length(l)]};
  else res:=ShallowCopy(l);
    for i in [1..Length(l)] do 
      res[AbsInt(sp[i])]:=l[i]*SignInt(sp[i]);
    od;
    return res;
  fi;
end;

# SignedPerm to matrix m such that m*l=SignPermuted(sp,l)
SignedPermutationMat:=function(arg)local n,res,i,sp,d,f;
  sp:=arg[1];
  if IsList(sp) then
    n:=Length(sp);
    res:=NullMat(n);
    for i in [1..n] do res[i][AbsInt(sp[i])]:=SignInt(sp[i]);od;
  else
    if Length(arg)=2 then n:=arg[2];
    else n:=QuoInt(LargestMovedPointPerm(sp)+1,2);fi;
    return SignedPermutationMat(SignedPerm(sp,n));
  fi;
  return res;
end;

# Transforms perm or matrix or perm+signs to a list signed perm,
# and a list signed perm to a perm
SignedPerm:=function(arg)local ls,n,i,nz,res;
  ls:=arg[1];
  if IsMat(ls) then
    n:=Length(ls);
    res:=[];
    for i in [1..n] do
      nz:=Filtered([1..n],x->ls[i][x]<>0);
      if Length(nz)<>1 then return false;fi;
      nz:=nz[1];
      if ls[i][nz]=1 then res[i]:=nz;
      elif ls[i][nz]=-1 then res[i]:=-nz;
      else return false;
      fi;
    od;
    return res;
  elif IsList(ls) then
    n:=Length(ls);
    ls:=Concatenation(
      List(ls,function(i)if i<0 then return i+2*n+1;else return i;fi;end),
      List(Reversed(ls),
        function(i)if i<0 then return -i;else return -i+2*n+1;fi;end));
    return PermList(ls);
  else
    if Length(arg)=2 then n:=arg[2];
      if IsList(n) then return Permuted(Zip([1..Length(n)],n,
          function(x,y)return x*y;end),ls^-1);fi;
    else n:=QuoInt(LargestMovedPointPerm(ls)+1,2);fi;
    return SignPermuted([1..n],ls^-1);
  fi;
end;

CyclesSignedPerm:=function(arg)local l,ps,n;
  ps:=arg[1];
  if IsList(ps) then return CyclesSignedPerm(SignedPerm(ps),Length(ps));fi;
  if Length(arg)=2 then n:=arg[2];
  else n:=QuoInt(LargestMovedPointPerm(ps)+1,2);
  fi;
  l:=Filtered(Cycles(ps),x->ForAny(x,x->x<=n));
  l:=List(l,x->List(x,function(y)if y>n then return y-2*n-1;
                                        else return y;fi;end));
  l:=l{Filtered([1..Length(l)],i->not ForAny(l{[1..i-1]},s->-l[i][1] in s))};
  return l;
# return List(l,function(s)
#   if -s[1] in s then return s{[1..QuoInt(Length(s),2)+1]};
#   else return s;fi;end);
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
