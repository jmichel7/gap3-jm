# sparse.g, (c) Jean Michel dec. 2003
#
# A few functions to work with 'sparse' matrices, which are reporesented as
# list v of records  where v[i].j is bound iff m[i][j]<>0
#
SparseOps:=OperationsRecord("SparseOps");

SparseOps.Print:=function(m)
  Print("<sparse ",Length(m.v),"x",Length(m.v)," matrix with ",
    Sum(m.v,x->Length(RecFields(x)))," bound entries");
end;

SparseOps.Transpose:=function(m)local res,i,g,number;
  number:=function(s)local res,c;
    res:=0;for c in s do res:=10*res+Position("0123456789",c)-1;od;
    return res;
  end;
  res:=rec(v:=List(m.v,x->rec()),operations:=SparseOps);
  for i in [1..Length(m.v)] do 
    for g in RecFields(m.v[i]) do
      res.v[number(g)].(i):=m.v[i].(g);
    od;
  od;
  return res;
end;

WGraphToSparseRepresentation:=function(rk,gr,v)local l,x,y,i,j,n,S,V,mu;
  V:=[];
  for S in gr[1] do 
    if  IsInt(S) then Append(V,List([1..S],i->V[Length(V)]));
    else Add(V,S);
    fi;
  od;
  n:=Length(V);
  S:=List([1..rk],i->rec(v:=[],operations:=SparseOps));
  for i in [1..rk] do for j in [1..n] do S[i].v[j]:=rec((j):=v^2);od;od;
  for j in [1..n] do for i in V[j] do S[i].v[j].(j):=-v^0;od;od;
  for i in gr[2] do 
    if IsList(i[1]) then mu:=i[1];else mu:=[i[1],i[1]];fi;
    for l in i[2] do 
      x:=l[1];
      for y in l{[2..Length(l)]} do
        for j in Difference(V[y],V[x]) do S[j].v[y].(x):=mu[2]*v;od;
        for j in Difference(V[x],V[y]) do S[j].v[x].(y):=mu[1]*v;od;
      od;
    od;
  od;
  return S;
end;

# Trace of a 'sparse matrix'
SparseOps.Trace:=function(r)local i,res;
  res:=0;
  for i in [1..Length(r.v)] do
    if IsBound(r.v[i].(i)) then res:=res+r.v[i].(i);fi;
  od;
  return res;
end;

SparseOps.\+:=function(r,s)local res,i,f,sum;
  res:=Copy(r);
  for i in [1..Length(r.v)] do 
    for f in RecFields(s.v[i]) do
      if IsBound(r.v[i].(f)) then sum:=r.v[i].(f)+s.v[i].(f);
        if sum<>0 then res.v[i].(f):=sum;else Unbind(res.v[i].(f));fi;
      else res.v[i].(f):=s.v[i].(f);
      fi;
    od;
  od;
  return res;
end;

SparseOps.\-:=function(r,s)return r+(-1)*s;end;

SparseOps.\*:=function(r,s)local i,j,res,f,rf,sf,sum,ff;
  if not IsRec(r) or not IsBound(r.operations) or r.operations<>SparseOps
  then return rec(v:=List(s.v,function(R)
    res:=rec();for f in RecFields(R) do res.(f):=r*R.(f);od;return res;end),
    operations:=SparseOps);
  fi;
  res:=rec(v:=List(r.v,x->rec()),operations:=SparseOps);
  s:=SparseOps.Transpose(s);
  rf:=List(r.v,x->Set(RecFields(x)));
  sf:=List(s.v,x->Set(RecFields(x)));
  for i in [1..Length(r.v)] do for j in [1..Length(r.v)] do 
      f:=Copy(rf[i]);IntersectSet(f,sf[j]);
      if f<>[] then 
        sum:=0;for ff in f do sum:=sum+r.v[i].(ff)*s.v[j].(ff);od;
        if sum<>0 then res.v[i].(j):=sum;fi;
      fi;
  od; od;
  return res;
end;

SparseOps.\^:=function(m,i)local i,y;
  if i=0 then 
    return rec(v:=List([1..Length(m.v)],i->rec((i):=1)),operations:=SparseOps);
  elif i>0 then 
    while i>0 do
     if i mod 2 <> 0 then if IsBound(y) then y:=y*m;else y:=m;fi;fi;
     if i>=2 then m:=m*m;fi;
     i:=QuoInt(i,2);
    od;
    return y;
  fi;
end;

SparseOps.\=:=function(r,s)
  return ForAll([1..Length(r.v)],i->RecFields(r.v[i])=RecFields(s.v[i])
    and ForAll(RecFields(r.v[i]),f->r.v[i].(f)=s.v[i].(f)));end;
