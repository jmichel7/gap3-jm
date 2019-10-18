#############################################################################
##
#A  sperm.g          CHEVIE library                               Jean Michel
##
#Y  Copyright (C) 2017 - 2019   University  Paris VII.
##
##  This  file contains routines to work with signed permutations. A signed
##  permutation is a compact way to work with permutations of [-n..-1,1..n]
##  which preserve the pairs [-i,i]. Such a permutation is represented by a
##  list  l  of  length  n  of  the  integers  in  [1..n] with signs, which
##  represents the permutation sending i to l[i]. It may be also convenient
##  in  computations  to  work  with  elements of the hyperoctaedral group,
##  acting  on [2,4..2*n]  instead of  [-1..-n] and  [1,3,2*n-1] instead of
##  [1..n].
##

SignedPermOps:=OperationsRecord("SignedPermOps");

SignedPermOps.String:=function(p)local l;
  l:=Cycles(p);
  if Length(l)=0 then return "()";
  else return Join(List(Cycles(p),c->SPrint("(",Join(c),")")),"");
  fi;
end;

SignedPermOps.Print:=function(x)Print(String(x));end;

SignedPermOps.\^:=function(p,i)
  if IsInt(i) then return SignedPerm(SignedPermOps.HO(p)^i,Length(p.l));
  elif IsInt(p) then return i.l[p]; 
  elif IsPerm(i) then return p^SignedPerm(ListPerm(i));
  else return i^-1*p*i; # i is a SignedPerm
  fi;
end;

SignedPermOps.\*:=function(p,q)
  if IsList(p) then return List(p,x->x*q);fi;
  if IsList(q) then return List(q,x->p*x);fi;
  return SignedPerm(SignedPermOps.HO(p)*SignedPermOps.HO(q),
    Maximum(Length(p.l),Length(q.l)));
end;

SignedPermOps.Comm:=function(p,q)return p^-1*q^-1*p*q;end;

SignedPermOps.\=:=function(a,b)local u;
  if not IsSignedPerm(b) then return false;fi;
  if Length(a.l)>Length(b.l) then 
    u:=[Length(b.l)+1..Length(a.l)];
    return a.l{[1..Length(b.l)]}=b.l and a.l{u}=u;
  else
    u:=[Length(a.l)+1..Length(b.l)];
    return b.l{[1..Length(a.l)]}=a.l and b.l{u}=u;
  fi;
end;

SignedPermOps.\<:=function(a,b)return a.l<b.l;end;

# Apply SignedPerm sp to list l
SignedPermOps.Permuted:=function(l,sp)local res,i;
  res:=ShallowCopy(l);
  for i in [1..Length(sp.l)] do 
    if sp.l[i]>0 then res[sp.l[i]]:=l[i];
    else res[-sp.l[i]]:=-l[i];
    fi;
  od;
  return res;
end;

# SignedPerm to matrix m such that m*l=Permuted(l,sp)
SignedPermOps.PermutationMat:=function(sp)local res,i,n;
  n:=Length(sp.l);
  res:=NullMat(n);
  for i in [1..n] do res[i][AbsInt(sp.l[i])]:=SignInt(sp.l[i]);od;
  return res;
end;

# argument SignedPerm or hyperoctaedral perm
SignedPermOps.Cycles:=function(p)local l,ps,n;
  ps:=SignedPermOps.HO(p);
  n:=Length(p.l);
  l:=Filtered(Cycles(ps),x->ForAny(x,x->x mod 2=1));
  l:=List(l,x->List(x,function(y)if y mod 2=1 then return QuoInt(y+1,2);
                                 else return -QuoInt(y,2);fi;end));
  l:=l{Filtered([1..Length(l)],i->not ForAny(l{[1..i-1]},s->-l[i][1] in s))};
  return l;
end;

# hyperoctaedral perm from SignedPerm
SignedPermOps.HO:=function(p)local n,ls,i,j;
  n:=Length(p.l);
  ls:=[1..2*n]*0;
  for j in [1..n] do
    i:=p.l[j];
    if i>0 then ls[2*j-1]:=2*i-1;ls[2*j]:=2*i;
    else        ls[2*j-1]:=-2*i;ls[2*j]:=-2*i-1;
    fi;
  od;
  return PermList(ls);
end;

IsSignedPerm:=x->IsRec(x) and IsBound(x.operations) 
    and x.operations=SignedPermOps;

#underlying Perm of a SignedPerm
Perm:=p->PermList(List(p.l,AbsInt));

#underlying Signs of a SignedPerm
Signs:=p->Permuted(List(p.l,SignInt),Perm(p));

# We have the properties p=SignedPerm(Perm(p),Signs(p)) and
# if N=OnMatrices(M,p) then
#    M=OnMatrices(N,Perm(p))^DiagonalMat(Signs(p)))

# Transforms perm or matrix or perm+signs to a list signed perm,
# and SignedPerm to an hyperoctaedral perm
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
    return SignedPerm(res);
  elif IsList(ls) then
    return rec(l:=ls,operations:=SignedPermOps,
    domain:=rec(operations:=rec(
      \in:=function(a,b)return IsSignedPerm(a);end,
      Group:=function(D,gens,id) return rec(isDomain:=true,isGroup:=true,
	identity:=id,generators:=gens,operations:=GroupOps,isFinite:=true);
      end)));
  else # hyperoctaedral perm, length
# Apply hyperoctaedral permutation sp to list l
    if Length(arg)=2 then n:=arg[2];
      if IsList(n) then  # perm, signs
        return SignedPerm(Permuted(Zip([1..Length(n)],n,
          function(x,y)return x*y;end),ls^-1));fi;
    elif ls=() then n:=1;
    else n:=QuoInt(LargestMovedPointPerm(ls)+1,2);fi;
    return SignedPerm(List(OnTuples([1,3..2*n-1],ls),function(i)
      if i mod 2=0 then return -i/2;else return QuoInt(i+1,2);fi;end));
  fi;
end;

# Find if exists signed perm which permutes list a to list b
SignedPermListList:=function(a,b)local p,i,res;
  p:=PermListList(List(a,x->Set([x,-x])),List(b,x->Set([x,-x])));
  if p=false then return false;fi;
  res:=Permuted([1..Length(a)],p);
  for i in [1..Length(a)] do
    if b[i^(p^-1)]<>a[i] then res[i]:=-res[i];fi;
  od;
  return SignedPerm(res);
end;

# duplicate lines and cols of M so HOgroup operates
SignedPermOps.dup:=function(M)local res,i,j;
  res:=List([1..2*Length(M)],i->[1..2*Length(M)]*0);
  for i in [1..Length(M)] do for j in [1..Length(M)] do
    res{[2*i-1,2*i]}{[2*j-1,2*j]}:=[[M[i][j],-M[i][j]],[-M[i][j],M[i][j]]];
  od;od;
  return res;
end;

# SignedMatStab(M [,extra]) find permutations with signs stabilizing M
# (and such that the associated permutations in addition stabilizes extra
#  -- which could be for instance a list of eigenvalues)
SignedMatStab:=function(arg)local blocks,stab,g,r,gens,I,p,ss,M,extra,k,gr;
  M:=arg[1];k:=Length(M);
  ss:=x->Set([x,-x]);
  if M<>TransposedMat(M) then Error("M should be symmetric");fi;
  if Length(arg)>1 then extra:=arg[2];else extra:=List(M,x->1);fi;
  blocks:=CollectBy([1..k],function(i)local inv;
    inv:=[Collected(List(M[i],ss)),M[i][i]];
    if IsBound(extra) then Add(inv,extra[i]);fi;
    return inv;end);
  g:=Group(PermList([]));I:=[];
  for r in blocks do
    if Length(r)>5 then InfoChevie("#I Large Block:",r,"\n");fi;
    gr:=MatStab(CoxeterGroupHyperoctaedralGroup(Length(r)),SignedPermOps.dup(M{r}{r}));
    p:=MappingPermListList([1..Length(r)],r);
    gens:=List(gr.generators,x->SignedPerm(x,Length(r))^p);
    g:=ApplyFunc(Group,Concatenation(g.generators,gens));
    Append(I,r); 
    p:=MappingPermListList(I,[1..Length(I)]);
    g:=ApplyFunc(Group,List(g.generators,x->x^p));
    g:=Group(List(g.generators,SignedPermOps.HO),());
    g:=MatStab(g,SignedPermOps.dup(M{I}{I}));
    g:=Group(List(g.generators,x->SignedPerm(x,k)),SignedPerm([]));
    g:=ApplyFunc(Group,List(g.generators,x->x^(p^-1)));
  od;
  return g;
end;

# SignedPermMatMat(M, N[, extra1, extra2]) find p such that PsOnMatrices(M,p)=N
# [and such that Permuted(extra1,p)=extra2]
SignedPermMatMat:=function(arg)
  local ind,l,I,J,r,p,e,g,n,h,trans,tr,q,ss,M,N,extra1,extra2,PsOnMatrices;
  PsOnMatrices:=function(M,p)return OnMatrices(M,SignedPerm(p,Length(M)));end;
  M:=arg[1];N:=arg[2];
  ss:=x->Set([x,-x]);
  if M<>TransposedMat(M) then Error("M should be symmetric");fi;
  if N<>TransposedMat(N) then Error("N should be symmetric");fi;
  if Length(arg)=2 then extra1:=List(M,x->1);extra2:=List(N,x->1);fi;
  ind:=function(I,J)local iM,iN,p,n;
    iM:=List(I,function(i)local inv;
      inv:=[Collected(List(M[i]{I},ss)),M[i][i]];
      if IsBound(extra1) then Add(inv,extra1[i]);fi;
      return inv;end);
    iN:=List(J,function(i)local inv;
      inv:=[Collected(List(N[i]{J},ss)),N[i][i]];
      if IsBound(extra2) then Add(inv,extra2[i]);fi;
      return inv;end);
    if Collected(iM)<>Collected(iN) then 
       InfoChevie("content differs");return false;fi;
    iM:=CollectBy(I,iM); iN:=CollectBy(J,iN);
    if Length(iM)=1 then
      if Length(I)>6 then InfoChevie("large block:",Length(I),"\n");
        p:=DistHelpedRepresentativeOperation(Group(Reflections(CoxeterGroupHyperoctaedralGroup(Length(I))),()),
           M{I}{I},N{J}{J},PsOnMatrices,
	   function(M,N)return Sum(M-N,x->Number(x,y->y<>0));end);
      else p:=RepresentativeOperation(CoxeterGroupHyperoctaedralGroup(Length(I)),
         M{I}{I},N{J}{J},PsOnMatrices);
      fi;
      if p=false then InfoChevie("could not match block");return false;fi;
      return [[I,J,SignedPerm(p,Length(I))]];
    else p:=Zip(iM,iN,ind);
      if false in p then return false;
      else return Concatenation(p);
      fi;
    fi;
  end;
  l:=ind([1..Length(M)],[1..Length(N)]);
  if l=false then return false;fi;
  I:=[];J:=[];g:=Group(SignedPerm([]));tr:=SignedPerm([]);
  for r in l do
#   Print("r=",r,"\n");
    n:=Length(r[1]);
    q:=MappingPermListList([1..Length(I)],[1..Length(I)]);
    p:=MappingPermListList([1..n],[1..n]+Length(I));
    Append(I,r[1]);Append(J,r[2]);
#   Print("#I=",Length(I),"\c");
    if Comm(r[3]^p,tr^q)<>SignedPerm([]) then Error("noncomm");fi;
    tr:=tr^q*r[3]^p;
    h:=OnTuples(SignedMatStab(M{r[1]}{r[1]}).generators,p);
    g:=Group(Concatenation(OnTuples(g.generators,q),h),SignedPerm([]));
#   Print(" #g=",Size(g),"\c");
    e:=RepresentativeOperation(g,M{I}{I},
      OnMatrices(N{J}{J},tr^-1),OnMatrices);
    if e=false then return false;
    else if e^-1*e^tr<>SignedPerm([]) then 
            Print("*** tr does not commute to e\n");
         fi;
         tr:=e*tr;
    fi;
    g:=MatStab(Group(List(g.generators,SignedPermOps.HO),()),SignedPermOps.dup(M{I}{I}));
#   Print(" #stab=",Size(g),"\n");
    g:=Group(List(g.generators,x->SignedPerm(x,Length(I))),SignedPerm([]));
  od;
  # transporter of a ps from [1..Length(I)] to I
  trans:=I->SignedPerm(ListPerm(MappingPermListList([1..Length(I)],I)));
  return trans(I)^-1*tr*trans(J);
end;
