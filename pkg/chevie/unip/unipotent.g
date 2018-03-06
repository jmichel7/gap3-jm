###########################################################################
##
#A  unipotent.g           CHEVIE library      Olivier Dudas and Jean Michel
##
#Y  Copyright (C) 2010 - 2015  Universite Paris VII
##
##  This file contains some programs for working with unipotent radicals of
##  Borel subgroups and their elements.
###########################################################################

# UnipotentGroup(W) represents the unipotent radical U of a Borel of the
# algebraic group attached to W.
#
# The reference for determining the structure of U is
#
#  R.W. Carter, "Simple groups of Lie type",  Wiley 1972, section 4.2.
#
# A  Chevalley basis of the Lie algebra of  U is a basis e_r indexed by the
# positive  roots  such  that  [e_r,e_s]=N_{r,s}  e_{r+s}  for some integer
# constants N_{r,s}.
#
# Let  < be an order on roots induced by a total order on the vector space,
# for  instance the lexicographic order on roots.  A pair (r,s) of roots is
# special if 0<r<s and r+s is a root.
#
# Constants C_{r,s,i,j} are defined [Carter 5.2.3] by
#
# u_s(u)u_r(t)=u_r(t)u_s(u)\prod_{i,j>0} u_{ir+js}(C_{r,s,i,j}(-t)^iu^j)
#
# Where ir+js runs over the roots of this shape.
#
# U is represented by  a record with the following data in its fields:
#
# .weylGroup           the underlying Weyl group
# .specialPairs        contains the triples  of indices of  roots (r,s,r+s) 
#                      where (r,s) is special, ordered by (r+s,r), followed
#                      by the triples (s,r,r+s) for the same list.
# .chevalleyConstants  the constants N_{r,s} for .specialPairs
# .order               an order on the positive roots
# .commutatorConstants stores the C_{r,s,i,j} by storing for each special pair
#                      (r,s) the list of quadruples [i,j,ir+js,C_{r,s,i,j}].

UnipotentElementOps:=OperationsRecord("UnipotentElementsOps");

UnipotentElement:=function(U,l)
  return rec(list:=l,group:=U,operations:=UnipotentElementOps);
end;

IsUnipotentElement:=u->IsRec(u) and IsBound(u.operations) and
                      u.operations=UnipotentElementOps;

UnipotentElementOps.Format:=function(u,opt)
  if u.list=[] then return "()";fi;
  return Join(List(u.list,function(p)local r;
    if IsBound(opt.root) then 
      r:=IntListToString(Parent(u.group.weylGroup).roots[p[1]]);
    else r:=String(p[1]);
    fi;
    if IsBound(opt.TeX) then return SPrint("u_{",r,"}(",p[2],")");
    else return SPrint("u",r,"(",p[2],")");
    fi;end)," * ");
end;

UnipotentElementOps.String:=u->Format(u,rec());

UnipotentElementOps.Display:=function(u,opt)
  Cut(Format(u,opt),rec(after:=" ",before:="*"));
end;

UnipotentElementOps.Print:=function(u)Print(String(u));end;

UnipotentElementOps.\*:=function(u,v)
  return u.group.Element(Concatenation(u.list,v.list));
end;

UnipotentElementOps.\/:=function(u,v)return u*v^-1;end;

UnipotentElementOps.\^:=function(u,n)local p,s,W,U;
  U:=u.group;W:=Parent(U.weylGroup);
  if IsInt(n) then 
    if n=0 then return UnipotentElement(U,[]);
    elif n=-1 then return U.Element(Reversed(List(u.list,p->[p[1],-p[2]])));
    elif n<0 then return (u^-1)^(-n);
    else p:=false;
      while n>0 do
       if n mod 2 <> 0 then if p=false then p:=u; else p:=p*u; fi; fi;
       u:=u*u; n:=QuoInt(n,2);
      od;
      return p;
    fi;
  elif IsUnipotentElement(n) then return n^-1*u*n; 
  elif IsSemisimpleElement(n) then 
     return UnipotentElement(U,List(u.list,p->[p[1],n^W.roots[p[1]]*p[2]]));
  else # assume n is an element of W
    # each root of u must remain positive by n
    if n=() then return u;fi;
    s:=FirstLeftDescending(W,n);
    p:=Filtered(u.list,x->x[1]^n>W.parentN);
    if p<>[] then
      Error(u," should have no coefficient on root ",p[1][1],"\n");fi;
    return U.Element(List(u.list,i->[i[1]^n,U.Eta(s,i[1])*i[2]]))
                     ^(W.reflections[s]*n);
  fi;
end;

UnipotentGroupOps:=OperationsRecord("UnipotentGroupOps");

# RootOrder(W)
# returns  the ordering  of the  roots of  W according  to the  linear form
# determined  by the weights Reversed(1+1/2,1+1/4,...)  on the simple roots
# (it is checked that this induces a total order on the roots).
# These weights are chosen so that the resulting order is almost always
# the default order on the roots.
RootOrder:=function(W)local r,i,weights;
  r:=W.roots{[1..W.N]};
  i:=Length(r[1]);weights:=1+List([i,(i-1)..1],j->1/2^j);
  if Length(Set(r*weights))<>Length(r) then
    Error("the linear form ",weights," takes the same value on two roots");
  fi;
  i:=SortingPerm(r*weights);
  if i<>() then InfoChevie("#I Need ",i," to sort the roots of ",W,"\n");fi;
  return Permuted(W.rootInclusion{[1..W.N]},i);
end;

UnipotentGroup:=function(W)local U,r,s,pos,order,i,r1,s1,l,M,PW;
  U:=rec(weylGroup:=W, operations:=UnipotentGroupOps);
  if W.roots{W.N+[1..W.N]}<>-W.roots{[1..W.N]} then Error();fi;
# We need an order on the roots coming from the ambient vector space.
  order:=RootOrder(W);
  order:=W.rootInclusion{[1..W.N]};
  PW:=Parent(W);
# compute special pairs
  U.specialPairs:=[];
  for s in order do
    for r in order{[1..Position(order,s)-1]} do
      pos:=Position(PW.roots,PW.roots[r]+PW.roots[s]);
      if pos<>false then Add(U.specialPairs,[r,s,pos]);fi;
    od;
  od;
  SortBy(U.specialPairs,x->[Position(order,x[3]),Position(order,x[1])]);
  U.ns:=Length(U.specialPairs);
  Append(U.specialPairs,List(U.specialPairs,x->x{[2,1,3]}));
  l:=function(r) # Length of root r
    if r=false then return 1;
    else return PW.rootLengths[PW.orbitRepresentative[r]];fi;
  end;
# Compute N_{a,b} for non-necessarily positive roots (see Carter 4.2.1 (ii))
  U.N:=function(a,b)local c,ra,rb,pos;
    if a<0 then ra:=-PW.roots[-a];else ra:=PW.roots[a];fi;
    if b<0 then rb:=-PW.roots[-b];else rb:=PW.roots[b];fi;
    c:=Position(PW.roots,ra+rb);
    if c=false then return 0;elif c>PW.N then c:=-(c-PW.N);fi;
    if a>0 then
      if b>0 then return U.chevalleyConstants[Position(U.specialPairs,[a,b,c])];
      elif c<0 then return  U.N(-c,a)*l(-c)/l(-b);
      else          return -U.N(-b,c)*l(c)/l(a);
      fi;
    elif b<0 then return -U.N(-a,-b);
    elif c<0 then return  U.N(b,-c)*l(-c)/l(-a);
    else          return -U.N(c,-a)*l(c)/l(b);
    fi;
  end;
# Compute N_{r,s} for each special pair.
# see formula in proof of Carter 4.2.2
  U.chevalleyConstants:=[];
  for i in [1..U.ns] do
    if i=1 or U.specialPairs[i-1][3]<>U.specialPairs[i][3] then # extraspecial
      r:=U.specialPairs[i][1]; s:=U.specialPairs[i][2];
      U.chevalleyConstants[i]:=Number([0..3],
                               j->PW.roots[s]-PW.roots[r]*j in PW.roots); 
      U.chevalleyConstants[i+U.ns]:=-U.chevalleyConstants[i];
    else # special of sum r+s
      r1:=U.specialPairs[i][1]; s1:=U.specialPairs[i][2]; 
      U.chevalleyConstants[i]:=l(U.specialPairs[i][3])/U.N(r,s)*
       (U.N(s,-r1)*U.N(r,-s1)/l(Position(PW.roots,PW.roots[s]-PW.roots[r1]))
       +U.N(-r1,r)*U.N(s,-s1)/l(Position(PW.roots,PW.roots[r]-PW.roots[r1])));
      U.chevalleyConstants[i+U.ns]:=-U.chevalleyConstants[i];
    fi;
  od;
  M:=function(r,s,i)local d,j,m; # M_{r,s,i} of [Carter, bottom of page 61]
    m:=1;
    for j in [1..i] do
      d:=Position(PW.roots,PW.roots[r]+PW.roots[s]);
      if d=false or d>PW.N then return 0;fi;
      m:=m*U.chevalleyConstants[Position(U.specialPairs,[r,s,d])]/j;
      s:=d;
    od;
    return m;  
  end;
  U.commutatorConstants:=List(U.specialPairs,function(p)local L,C,c;
    L:=[];
    for c in [[1,1],[2,1],[1,2],[3,1],[1,3],[3,2],[2,3]] do
    # possible (i,j) such that there may exist a root is+jr.
    # see [Carter, top of page 76] for the formulas
      if c[2]=1 then C:=M(p[1],p[2],c[1]);
      elif c[1]=1 then C:=(-1)^c[2]*M(p[2],p[1],c[2]);
      elif c{[1,2]}=[3,2] then C:=M(p[3],p[1],2)/3;
      elif c{[1,2]}=[2,3] then C:=-2*M(p[3],p[2],2)/3;
      else C:=0; 
      fi;
      if C<>0 then Add(L,[c[1],c[2],
        Position(PW.roots,c[1]*PW.roots[p[1]]+c[2]*PW.roots[p[2]]),C]);fi;
    od;
    return L;end);
  U.order:=W.rootInclusion{[1..W.N]}; # default ordering of roots
# Computes the constants \eta_{r,s} defined in [Carter, 6.4.2 and 6.4.3]
# allowing to conjugate a unipotent element by a reflection.
  U.Eta:=function(a,b)local p,q,L,i,eta,N;
    L:=List([-4..4],j->Position(PW.roots,PW.roots[a]*j+PW.roots[b]));
    L:=Filtered([-4..4],j->L[5+j]<>false and L[5+j]<=PW.N);
    p:=-L[1];q:=L[Length(L)];
    eta:=(-1)^p;
    for i in [0..Maximum(p-1,q-1)] do
      N:=U.N(a,Position(PW.roots,PW.roots[a]*(i-p)+PW.roots[b]));
      if i<p or i<q then eta:=eta*SignInt(N);fi;
    od;
    return eta;
  end;
# U.CanonicalForm(list[,order])
# returns list in canonical order U.order or the argument order if given.
  U.CanonicalForm:=function(arg)local i,l,order,c,res;
    l:=arg[1];
    if Length(arg)=2 then order:=arg[2];else order:=U.order;fi;
    i:=1;
    while i<=Length(l) do 
      if l[i][2]=0*l[i][2] then
	l:=Concatenation(l{[1..i-1]},l{[i+1..Length(l)]});
	if i>1 then i:=i-1;fi;
      elif i<Length(l) and l[i][1]=l[i+1][1] then 
	l:=Concatenation(l{[1..i-1]},[[l[i][1],l[i][2]+l[i+1][2]]],
	 l{[i+2..Length(l)]});
      elif i<Length(l) and Position(order,l[i][1])>Position(order,l[i+1][1]) 
      then 
  # Transform u_s(u) u_r(t) by Chevalley relation
  # 
  # u_s(u)u_r(t)=u_r(t)u_s(u)\prod_{i,j>0} u_{ir+js}(C_{r,s,i,j}(-t)^iu^j)
  #
  # Here l[i]=[s,u] et l[i+1]=[r,t]. 
	res:=Concatenation(l{[1..i-1]},l{[i+1,i]});
	c:=Position(PW.roots,PW.roots[l[i+1][1]]+PW.roots[l[i][1]]);
	if c<>false then c:=Position(U.specialPairs,[l[i+1][1],l[i][1],c]);fi;
	if c<>false then c:=List(U.commutatorConstants[c],
	    k->[k[3],k[4]*(-l[i+1][2])^k[1]*l[i][2]^k[2]]);
	  Append(res,Filtered(c,x->x[2]<>0*x[2]));
	fi;
	Append(res,l{[i+2..Length(l)]});
	l:=res;
	if i>1 then i:=i-1;fi;
      else i:=i+1;
      fi;
    od;
    return l;
  end;
  U.Element:=function(arg)local i,l,order,c,res;
    if IsList(arg[1]) then return UnipotentElement(U,U.CanonicalForm(arg[1]));
    elif Length(arg)=1 then return UnipotentElement(U,[[arg[1],1]]);
    else return UnipotentElement(U,U.CanonicalForm(List([1..Length(arg)/2],
      i->arg{2*i+[-1,0]})));
    fi;
  end;
  return U;
end;

UnipotentGroupOps.Print:=function(U)
  Print("UnipotentGroup(",U.weylGroup,")");end;

# Projection of u on U/D(U)
UnipotentAbelianPart:=u->UnipotentElement(u.group,
  Filtered(u.list,x->x[1] in u.group.weylGroup.rootInclusion
    {u.group.weylGroup.generatingReflections}));

# Decompose u on subgroups U\cap (U^-)^w and U\cap U^w
UnipotentDecompose:=function(w,u)local order,U,N;
  U:=u.group;N:=U.weylGroup.parentN;
  order:=Concatenation(Filtered(U.order,i->i^w>N),Filtered(U.order,i->i^w<=N)); 
  u:=U.CanonicalForm(u.list,order);
  return [UnipotentElement(U,Filtered(u,i->i[1]^w>N)),
          UnipotentElement(U,Filtered(u,i->i[1]^w<=N))];
end;
