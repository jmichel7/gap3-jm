#############################################################################
##
#A  po.g              CHEVIE library          Jean Michel
##
#Y  Copyright (C) 2003 - 2010  University  Paris VII.
##
##  This  file  contains  useful routines to deal with posets, equivalence
##  relations, set partitions.
##
##  Posets are records with at least one of the two fields:
#   .incidence, a boolean matrix of entry [i][j] true iff i<=j
#   .hasse a list representing the Hasse diagram: the ith entry is the list
#   of  indices of  elements which  are immediate  successors to i (covers)
#   (i.e. the list of j such that i<j and there is no k such that i<k<j)
##
# The arguments should be partitions p1,..,pn of a set.
# Returns the finest partition refined by all argument partitions
# It is also the 'or' of the equivalence relations.
LcmPartitions:=function(arg)local lcm2,res,x;
  lcm2:=function(a,b)local res,p; res:=[];
  for p in b do Add(res,Union(Filtered(a,x->Intersection(x,p)<>[])));od;
  b:=[];
  for p in res do Add(b,Union(Filtered(res,x->Intersection(x,p)<>[])));od;
  return Set(b);
  end;
  res:=arg[1];
  for x in arg{[2..Length(arg)]} do res:=lcm2(res,x);od;
  return res;
end;

# The arguments should be partitions p1,..,pn of a set.
# Returns the coarsest partition which refines by all argument partitions
# It is also the 'and' of the equivalence relations.
GcdPartitions:=function(arg)local gcd2,res,x;
  gcd2:=function(a,b)local res; 
  res:=List(a,x->List(b,y->Intersection(x,y)));
  return Set(Filtered(Concatenation(res),x->x<>[]));
  end;
  res:=arg[1];
  for x in arg{[2..Length(arg)]} do res:=gcd2(res,x);od;
  return res;
end;

# Given Hasse diagram of a poset returns a compatible linear order
LinearExtension:=function(P)local n,v,Q,res,x,ord;
  ord:=Hasse(P);n:=List(ord,x->0);
  for v in ord do for x in v do n[x]:=n[x]+1;od;od;
  Q:=Filtered([1..Length(n)],x->n[x]=0);
  res:=[];
  while Length(Q)>0 do
    Add(res,Q[1]);
    for x in ord[Q[1]] do 
      n[x]:=n[x]-1;
      if n[x]=0 then Add(Q,x);fi;
    od;
    Q:=Q{[2..Length(Q)]};
  od;
  if Sum(n)>0 then Error("cycle");fi;
  return res;
end;

# Hasse diagram of poset p (list of covers)
Hasse:=function(p)local ind,m;
  if not IsBound(p.hasse) then
    if false then #old algorithm
      ind:=[1..Length(p.incidence)];
      p.hasse:=List(ind,function(i)local gt;gt:=ListBlist(ind,p.incidence[i]);
	return Filtered(gt,j->SizeBlist(p.incidence{gt}[j])=2);end);
    else # 6 times faster for a matrix 1500x1500
      m:=List(p.incidence,x->List(x,
        function(y)if y then return 1;else return 0;fi;end));
      p.hasse:=List(m*m,x->Filtered([1..Length(x)],y->x[y]=2));
    fi;
  fi;
  return p.hasse;
end;

# Incidence matrix of poset p
Incidence:=function(p)local n,i,x;
  if not IsBound(p.incidence) then
    n:=LinearExtension(p);
    p.incidence:=List([1..Length(n)],x->List([1..Length(n)],y->y=x));
    for i in [Length(n)-1,Length(n)-2..1] do
      for x in p.hasse[n[i]] do
	p.incidence[n[i]]:=UnionBlist(p.incidence[n[i]],p.incidence[x]);
      od;
    od;
  fi;
  return p.incidence;
end;
   
# return a set of chains covering the Hasse diagram h
Chains:=function(P)local ch,i,j,p,h;ch:=[];
  h:=Hasse(P);
  for i in LinearExtension(P) do
    for j in h[i] do
      p:=PositionProperty(ch,c->i=c[Length(c)]);
      if p<>false then Add(ch[p],j);
      else Add(ch,[i,j]);
      fi;
    od;
  od;
  return ch;
end;

# Given  a poset  p returns  the associated  equivalence relation (elements
# with same predecessors and successors are equivalent).
# If computed from incidence matrix this partition is in linear order
# compatible with p.
Partition:=function(p)local l,ind,I;
  if IsBound(p.hasse) then
    l:=Reversed(p);
    l:=CollectBy([1..Length(p.hasse)],i->[l.hasse[i],p.hasse[i]]);
  else I:=Incidence(p);ind:=[1..Length(I)];
    l:=List(ind,i->[List(ind,j->j<>i and I[i][j]),
		    List(ind,j->j<>i and I[j][i])]);
    l:=List(Set(l),x->Filtered(ind,i->l[i]=x));
  fi;
  return l;
end;

PosetOps:=OperationsRecord("PosetOps");

PosetOps.Format:=function(x,opt)local s,labels,p,sep;
  p:=Partition(x);s:=Hasse(x);
  s:=Poset(List(p,x->Set(List(s[x[1]],y->PositionProperty(p,z->y in z)))));
  labels:=List(p,y->Join(List(y,function(n)
    if IsBound(x.label) then return x.label(x,n,opt);else return String(n);fi;
  end)));
  if IsBound(opt.symbol) then sep:=opt.symbol;
  elif IsBound(opt.TeX) then sep:="{<}";
  else sep:="<";fi;
  s:=List(Chains(s),x->Join(labels{x},sep));
  if IsBound(opt.TeX) then 
    return SPrint("\\noindent",
      Concatenation(List(s,x->SPrint("$",x,"$\\hfill\\break\n"))));
  else return Concatenation(List(s,x->SPrint(x,"\n")));
  fi;
end;

PosetOps.Display:=function(x,opt)Print(Format(x,opt));end;

PosetOps.String:=x->SPrint("Poset with ",Size(x)," elements");

PosetOps.Print:=function(x)Print(String(x));end;

# Restricted(poset,indices) restriction of poset to indices
PosetOps.Restricted:=function(a)local p,ind,res;
  p:=a[1];ind:=a[2]; 
  res:=ShallowCopy(p);
  if Length(ind)=Size(p) and Set(ind)=[1..Size(p)] then # just renumbering
    if IsBound(res.hasse) then
      res.hasse:=List(res.hasse{ind},x->List(x,y->Position(ind,y)));
    fi;
    if IsBound(res.incidence) then
       res.incidence:=res.incidence{ind}{ind};
    fi;
  else Incidence(p);
    res.incidence:=List([1..Length(ind)],i->List([1..Length(ind)],function(j)
      if i<>j and ind[i]=ind[j] then return false;
      else return p.incidence[ind[i]][ind[j]];
      fi;end));
    Unbind(res.hasse);
  fi;
  res.size:=Length(ind);
  if IsBound(p.label) then res.indices:=ind;
    res.label:=function(x,n,opt)return p.label(x,x.indices[n],opt);end;
  fi;
  return res;
end;

# opposite poset
PosetOps.Reversed:=function(p)local res,i,j;
  res:=ShallowCopy(p);
  if IsBound(p.incidence) then res.incidence:=TransposedMat(p.incidence);fi;
  if IsBound(p.hasse) then
    res.hasse:=List(p.hasse,x->[]);
    for i in [1..Length(p.hasse)] do 
      for j in p.hasse[i] do Add(res.hasse[j],i);od;
    od;
  fi;
  return res;
end;

# creates poset from m either Hasse diagram or incidence matrix
Poset:=function(arg)local res,m; m:=arg[1];
  if IsList(m) then
    if IsBound(m[1]) and IsList(m[1]) and IsBound(m[1][1]) and IsBool(m[1][1])
    then return rec(incidence:=m,size:=Length(m),operations:=PosetOps);
    else return rec(hasse:=m,size:=Length(m),operations:=PosetOps);
    fi;
  elif IsRec(m) and IsBound(m.operations) and IsBound(m.operations.Poset) then
    return ApplyFunc(m.operations.Poset,arg);
  else Error(m," has no method for Poset\n");
  fi;
end;

# checks if elements of Poset have joins (upper bounds)
IsJoinLattice:=function(P)local subl,i,j,l,k,ll,ord;
  ord:=Incidence(P);
  subl:=[];
  for i in [1..Length(ord)] do
    for j in [1..i-1] do
      if not ord[j][i] or ord[i][j] then
        l:=IntersectionBlist(ord[i],ord[j]); # elts >= both i and j
	if not l in subl then
          if not ForAny(ListBlist(ord,l),y->l=IntersectionBlist(l,y)) then
	    for k in ListBlist([1..Length(ord)],l) do 
	      ll:=ShallowCopy(ord[k]);ll[k]:=false;SubtractBlist(l,ll);
	    od;
	    InfoChevie("# ",i," and ",j," have bounds",
	      ListBlist([1..Length(ord)],l),"\n");
	    return false;fi;
          AddSet(subl,l);
	fi;
      fi;
    od;
  od;
  return true;
end;

# checks if elements of Poset have meets (lower bounds)
IsMeetLattice:=function(P)local subl,i,j,l,k,ll,ord;
  ord:=TransposedMat(Incidence(P));
  subl:=[];
  for i in [1..Length(ord)] do
    for j in [1..i-1] do
      if not ord[j][i] or ord[i][j] then
        l:=IntersectionBlist(ord[i],ord[j]); # elts >= both i and j
	if not l in subl then
          if not ForAny(ListBlist(ord,l),y->l=IntersectionBlist(l,y)) then
	    for k in ListBlist([1..Length(ord)],l) do 
	      ll:=ShallowCopy(ord[k]);ll[k]:=false;SubtractBlist(l,ll);
	    od;
	    InfoChevie("# ",i," and ",j," have bounds",
	      ListBlist([1..Length(ord)],l),"\n");
	    return false;fi;
          AddSet(subl,l);
	fi;
      fi;
    od;
  od;
  return true;
end;
