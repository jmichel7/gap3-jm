###############################################################################
##
##  grape.g  (Version 2.31)       GRAPE Library               Leonard Soicher
##
##
##  Copyright 1992-1997 Leonard Soicher, School of Mathematical Sciences, 
##                      Queen Mary and Westfield College, London E1 4NS, U.K.
##
#
# This file contains the GRAPE 2.31 library.  
#
# Main changes since GRAPE 2.2:
#
# (1) The new  CompleteSubgraphsOfGivenSize,  which allows for searching in
# a vertex-weighted graph for cliques with a given vertex-weight sum.
#
# (2) The new function  PartialLinearSpaces,  which classifies partial linear
# spaces with given point graph and parameters  s,t.  
#
# (3) The new function VertexTransitiveDRGs, which determines the 
# distance-regular generalized orbital graphs for a given transitive 
# permutation group.
#
# (4) New functions  CayleyGraph  and  SwitchedGraph.
#
# (5) A one-vertex graph is now considered to be bipartite, 
# with bicomponents = [[],[1]] (to be consistent with considering a 
# zero-vertex graph to be bipartite, with bicomponents = [[],[]]). 
#
# (6) A bug fixed in the function UnderlyingGraph. That bug had
# the effect that if the returned graph had loops, then it 
# might have had its isSimple component erroneously set to true.
#
# The development of GRAPE was partially supported by the European Union
# HCM grant in "Computational Group Theory".
#
 
GRAPE_RANDOM := false; # Determines if random methods are to be used
		       # in  GRAPE  functions. 
		       # The default is that random methods are not 
		       # used (GRAPE_RANDOM=false).
		       # If random methods are used (GRAPE_RANDOM=true), 
		       # the correctness of the graph theoretical results 
		       # do *not* depend on the correctness of the randomly 
		       # computed objects. The use of random methods only 
		       # influences the time taken for certain graph 
		       # theoretical functions.

GRAPE_NRANGENS := 12;  # The number of random generators taken for a subgroup
		       # when  GRAPE_RANDOM=true.

MakeStabChainGivenLimit := function(G,limit)
#
# Makes a stabilizer chain for the permutation group  G,  
# given an upper bound  limit  on the size of  G.  
# 
local H;
if not IsPermGroup(G) or not IsInt(limit) then
   Error("usage: MakeStabChainGivenLimit( <PermGroup>, <Int> )");
fi;
if IsBound(G.stabChain) then
   return;
fi;
H:=Group(G.generators,());
StabChain(H,rec(random:=800,limit:=limit));
if Size(H) > limit then
   Error("upper bound <limit> on size of <G> incorrect");
elif Size(H)=limit then 
   # Size(G)=limit
   StabChain(G,rec(size:=limit));
else
   # Size(G) < limit  or  random method failed.
   StabChain(G);
fi;
end;

OrbitRepresentatives := function(arg)
#
# Let  G=arg[1],  L=arg[2],  op=arg[3]  (default: OnPoints).  Then this 
# function returns a list of representatives of  Orbits(G,L,op)  
# (one representative per orbit),  although these representatives should 
# not be regarded as canonical in any way.
#
local G,L,op,i,j,len,f;
G:=arg[1];
L:=arg[2];
if IsBound(arg[3]) then 
   op:=arg[3];
else
   op:=OnPoints;
fi;
if not (IsGroup(G) and IsList(L) and IsFunc(op)) then 
   Error("usage: OrbitRepresentatives( <Group>, <List> [, <Func> ] )");
fi;
len:=Length(L);
f:=BlistList([1..len],[1..len]);
for i in [1..len] do
   if f[i] then
      for j in [i+1..len] do 
	 if f[j] and RepresentativeOperation(G,L[i],L[j],op)<>false then
	    f[j]:=false;
	 fi;
      od;
   fi;
od;
return ListBlist(L,f);
end;

OrbitNumbers := function(G,n)
#
# Returns the orbits of  G  on  [1..n]  in the form of a record 
# containing the orbit representatives and a length  n  list 
# orbitNumbers,  such that  orbitNumbers[j]=i  means that point 
# j  is in the orbit of the i-th representative.
#
local i,j,orbnum,reps,im,norb,g,orb;
if not IsPermGroup(G) or not IsInt(n) then 
   Error("usage: OrbitNumbers( <PermGroup>, <Int> )");
fi;
orbnum:=[];
for i in [1..n] do
   orbnum[i]:=0;
od;
reps:=[];
norb:=0;
for i in [1..n] do
   if orbnum[i]=0 then      # new orbit
      Add(reps,i);
      norb:=norb+1;
      orbnum[i]:=norb;
      orb:=[i];
      for j in orb do 
	 for g in G.generators do
	    im:=j^g;
	    if orbnum[im]=0 then 
	       orbnum[im]:=norb;
	       Add(orb,im); 
	    fi; 
	 od; 
      od; 
   fi;
od;
return rec(representatives:=reps,orbitNumbers:=orbnum);
end;

NumbersToSets := function(vec)
#
# Returns a list of sets as described by the numbers in vec, i.e.
# i is in the j-th set iff vec[i]=j>0.
#
local list,i,j;
if not IsList(vec) then 
   Error("usage: NumbersToSets( <List> )");
fi;
if Length(vec)=0 then
   return [];
fi;
list:=[];
for i in [1..Maximum(vec)] do
   list[i]:=[];
od;
for i in [1..Length(vec)] do
   j:=vec[i];
   if j>0 then 
      Add(list[j],i);
   fi;
od;
for i in [1..Length(list)] do
   IsSet(list[i]);
od;
return list;
end;

MaximumMovedPoint := function(perms)
#
# Returns the maximum largest point moved by an element of  perms.
# 0  is returned iff  perms  is empty or consists only of identity 
# permutations.
#
local m,t,g;
if not IsList(perms) then 
   Error("usage: MaximumMovedPoint( <List> )");
fi;
m:=0;
for g in perms do
   if g <> () then
      t:=LargestMovedPointPerm(g);
      if t > m then 
	 m:=t;
      fi;
   fi;
od;
return m;
end;

IntransitiveGroupGenerators := function(arg)
local conjperm,i,newgens,gens1,gens2,max1,max2;
gens1:=arg[1];
gens2:=arg[2];
if IsBound(arg[3]) then
   max1:=arg[3];
else
   max1:=MaximumMovedPoint(gens1);
fi;
if IsBound(arg[4]) then
   max2:=arg[4];
else
   max2:=MaximumMovedPoint(gens2);
fi;
if not (IsList(gens1) and IsList(gens2) and IsInt(max1) and IsInt(max2)) then 
   Error(
   "usage: IntransitiveGroupGenerators( <List>, <List> [,<Int> [,<Int> ]] )");
fi;
if Length(gens1)<>Length(gens2) then
   Error("Length(<gens1>) <> Length(<gens2>)");
fi;
conjperm:=PermList(Concatenation(List([1..max2],x->x+max1),[1..max1]));
newgens:=[];
for i in [1..Length(gens1)] do
   newgens[i]:=gens1[i]*(gens2[i]^conjperm);
od;
return newgens;
end;

ProbablyStabilizer := function(G,pt)
#
# Returns a subgroup of  Stabilizer(G,pt),  which is very often 
# the full stabilizer. In fact, if GRAPE_RANDOM=false, then it
# is guaranteed to be the full stabilizer.
#
local sch,orb,x,y,im,k,gens,g,stabgens,i;
if not IsPermGroup(G) or not IsInt(pt) then
   Error("usage: ProbablyStabilizer( <PermGroup>, <Int> )");
fi;
if not GRAPE_RANDOM then
   return Stabilizer(G,pt);
fi;
k:=Length(G.generators);
if k=0 then
   return Group(());
fi;
# Make a Schreier vector of permutations for the orbit  pt^G.
sch:=[];
orb:=[pt];
sch[pt]:=();
for x in orb do
   for g in G.generators do
      im:=x^g;
      if not IsBound(sch[im]) then
	 sch[im]:=g;
	 Add(orb,im);
      fi;
   od;
od; 
# Now make a randomish generating sequence  gens  for  G.
gens:=Copy(G.generators);
if k > 1 then
   for i in [k+1..GRAPE_NRANGENS] do
      gens[i]:=gens[RandomList([1..k])];
   od;
   for i in [1..Maximum(k,GRAPE_NRANGENS)] do
      gens[i]:=Product(gens);
   od;
fi;
# Now make a list  stabgens  of random elements of the stabilizer of  pt.
x:=();
stabgens:=[];
for i in [1..GRAPE_NRANGENS] do
   x:=x*RandomList(gens);
   im:=pt^x;
   while im<>pt do
      x:=x/sch[im];
      im:=im/sch[im];
   od;
   if x<>() then
      Add(stabgens,x);
   fi;
od;
return Group(stabgens,());
end;

ProbablyStabilizerOrbitNumbers := function(G,pt,n)
#
# Returns the "orbit numbers" record for a subgroup of  Stabilizer(G,pt), 
# in its action on  [1..n].  
# This subgroup is very often the full stabilizer, and in fact, 
# if  GRAPE_RANDOM=false,  then it is guaranteed to be the full stabilizer. 
#
if not IsPermGroup(G) or not IsInt(pt) or not IsInt(n) then
   Error(
   "usage: ProbablyStabilizerOrbitNumbers( <PermGroup>, <Int>, <Int>  )");
fi;
return OrbitNumbers(ProbablyStabilizer(G,pt),n);
end;

RepWord := function(gens,sch,r)
#
# Given a sequence  gens  of group generators, and a  (word type)
# schreier vector  sch  made using  gens,  this function returns a 
# record containing the orbit representative for  r  (wrt  sch),  and
# a word in  gens  taking this representative to  r. 
# (We assume  sch  includes the orbit of  r.)
#
local word,w;
word:=[]; 
w:=sch[r];
while w > 0 do
   Add(word,w); 
   r:=r/gens[w]; 
   w:=sch[r];
od;
return rec(word:=Reversed(word),representative:=r);
end;
   
NullGraph := function(arg)
#
# Returns null graph with  n  vertices and group  G=arg[1].
# If  arg[2]  is bound then  n=arg[2],  otherwise  n  is the maximum 
# largest moved point of the generators of  G.
# The  names,  autGroup,  and  canonicalLabelling  components of the 
# returned null graph are left unbound; however, the  isSimple  
# component is set (to true).
#
local G,n,gamma,nadj,sch,orb,i,j,k,im,gens;
G:=arg[1];
if not IsPermGroup(G) or (IsBound(arg[2]) and not IsInt(arg[2])) then
   Error("usage: NullGraph( <PermGroup>, [, <Int> ] )");
fi;
n:=MaximumMovedPoint(G.generators); 
if IsBound(arg[2]) then
   if arg[2] < n  then
      Error("<arg[2]> too small");
   fi;
   n:=arg[2];
fi;
gamma:=rec(isGraph:=true,order:=n,group:=G,schreierVector:=[],
	   adjacencies:=[],representatives:=[],isSimple:=true);
if gamma.order=0 then
   return gamma;
fi;
#
# Make  gamma.schreierVector  and  gamma.adjacencies.
#
sch:=gamma.schreierVector; 
gens:=gamma.group.generators; 
nadj:=0;
for i in [1..n] do 
   sch[i]:=0; 
od;
for i in [1..n] do
   if sch[i]=0 then      # new orbit
      Add(gamma.representatives,i);
      nadj:=nadj+1;
      sch[i]:=-nadj;     # tells where to find the adjacency set.
      gamma.adjacencies[nadj]:=Set([]);
      orb:=[i];
      for j in orb do 
	 for k in [1..Length(gens)] do
	    im:=j^gens[k];
	    if sch[im]=0 then 
	       sch[im]:=k; 
	       Add(orb,im); 
	    fi; 
	 od; 
      od; 
   fi;
od;
return gamma;
end;

CompleteGraph := function(arg)
#
# Returns a complete graph with  n  vertices and group  G=arg[1].
# If  arg[2]  is bound then  n=arg[2],  otherwise  n  is the maximum 
# largest moved point of the generators of the permutation group  G.
# If the boolean argument  arg[3]  is bound and has 
# value true then the complete graph will have all possible loops, 
# otherwise it will have no loops (the default).
#
# The  names,  autGroup,  and  canonicalLabelling  components of the 
# returned complete graph are left unbound; however, the  isSimple  
# component is set.
#
local G,n,gamma,i,mustloops;
G:=arg[1];
if not IsPermGroup(G) or (IsBound(arg[2]) and not IsInt(arg[2]))
		      or (IsBound(arg[3]) and not IsBool(arg[3])) then
   Error("usage: CompleteGraph( <PermGroup>, [, <Int> [, <Bool> ]] )");
fi;
n:=MaximumMovedPoint(G.generators); 
if IsBound(arg[2]) then
   if arg[2] < n  then
      Error("<arg[2]> too small");
   fi;
   n:=arg[2];
fi;
if IsBound(arg[3]) then
   mustloops:=arg[3];
else 
   mustloops:=false;
fi;
gamma:=NullGraph(G,n);
if gamma.order=0 then
   return gamma;
fi;
if mustloops then 
   gamma.isSimple:=false;
fi;
for i in [1..Length(gamma.adjacencies)] do
   gamma.adjacencies[i]:=Set([1..n]);
   if not mustloops then
      RemoveSet(gamma.adjacencies[i],gamma.representatives[i]);
   fi;
od;
return gamma;
end;

Graph := function(arg)
#
# First suppose that  arg[5]  is unbound or has value  false.
# Then  L=arg[2]  is a list of elements of a set  S  on which 
# G=arg[1]  acts with action  act=arg[3].  Also  rel=arg[4]  is a boolean
# function defining a  G-invariant relation on  S  (so that 
# for  g in G,  rel(x,y)  iff  rel(act(x,g),act(y,g)) ). 
# Then function  Graph  returns the graph  gamma  with vertex names
# Concatenation(Orbits(G,L,act)),  and  x  is joined to  y
# in  gamma  iff  rel(VertexName(gamma,x),VertexName(gamma,y)).
#
# If  arg[5]  has value  true  then it is assumed that  L=arg[2] 
# is invariant under  G=arg[1]  with action  act=arg[3]. Then
# the function  Graph  behaves as above, except that  gamma.names
# becomes a copy of  L.
#
local G,L,act,rel,invt,gamma,i,vertices,reps,H,orb,x,y;
G:=arg[1];
L:=arg[2];
act:=arg[3];
rel:=arg[4];
if IsBound(arg[5]) then
   invt:=arg[5];
else
   invt:=false;
fi;
if not (IsGroup(G) and IsList(L) and IsFunc(act) and IsFunc(rel) 
	and IsBool(invt)) then
   Error("usage: Graph( <Group>, <List>, <Func>, <Func> [, <Bool> ] )");
fi;
if invt then
   vertices:=Copy(L);
else
   vertices:=Concatenation(Orbits(G,L,act));
fi;
gamma:=NullGraph(Operation(G,vertices,act),Length(vertices));
if not GRAPE_RANDOM then
   if (IsBound(G.size) and G.size<>"infinity") or 
      (IsPermGroup(G) and IsBound(G.stabChain)) or
      (IsBound(G.operations) and IsBound(G.operations.name) 
       and G.operations.name="SymmetricPermGroupOps") then
      #
      # Size(G) is an integer and should be easy to compute.
      #
      MakeStabChainGivenLimit(gamma.group,Size(G));
   fi;
fi;
Unbind(gamma.isSimple);
gamma.names:=vertices;
reps:=gamma.representatives;
for i in [1..Length(reps)] do
   H:=ProbablyStabilizer(gamma.group,reps[i]);
   x:=gamma.names[reps[i]];
   if Length(H.generators)=0 then  #  H  is trivial.
      gamma.adjacencies[i]:=Filtered([1..gamma.order],i->rel(x,gamma.names[i]));
      IsSet(gamma.adjacencies);
   else
      for orb in Orbits(H,[1..gamma.order]) do
	 y:=gamma.names[orb[1]];
	 if rel(x,y) then
	    UniteSet(gamma.adjacencies[i],orb);
	 fi;
      od;
   fi;
od;
return gamma;
end;

JohnsonGraph := function(n,e)
#
# Returns the Johnson graph, whose vertices are the e-subsets
# of {1,...,n},  with x joined to y iff  Intersection(x,y)
# has size  e-1.
#
local rel,J;
if not IsInt(n) or not IsInt(e) then 
   Error("usage: JohnsonGraph( <Int>, <Int> )");
fi;
if e<0 or n<e then
   Error("must have 0 <= <e> <= <n>");
fi;
rel := function(x,y)
   return Length(Intersection(x,y))=e-1; 
end;
J:=Graph(SymmetricGroup(n),Combinations([1..n],e),OnSets,rel,true);
J.isSimple:=true;
return J;
end;

IsGraph := function(obj)
#
# Returns  true  iff  obj  is a graph.
#
return IsRec(obj) and IsBound(obj.isGraph) and obj.isGraph;
end;

OrderGraph := function(gamma)
#
# returns the order of  gamma.
#
if not IsGraph(gamma) then
   Error("usage: OrderGraph( <Graph> )");
fi;
return gamma.order;
end;

Vertices := function(gamma)
#
# Returns the vertex set of graph  gamma.
#
if not IsGraph(gamma) then
   Error("usage: Vertices( <Graph> )");
fi;
return [1..gamma.order];
end;

IsVertex := function(gamma,v)
#
# Returns  true  iff  v  is vertex of  gamma.
#
if not IsGraph(gamma) then
   Error("usage: IsVertex( <Graph>, <obj> )");
fi;
return IsInt(v) and v >= 1 and v <= gamma.order;
end;

AssignVertexNames := function(gamma,names)
#
# Assign vertex names for  gamma,  so that the (external) name of 
# vertex  i  becomes  names[i].
#
if not IsGraph(gamma) or not IsList(names) then
   Error("usage: AssignVertexNames( <Graph>, <List> )");
fi;
gamma.names:=Copy(names);
end;

VertexName := function(gamma,v)
#
# Returns (a copy of) the (external) name of the vertex  v  of  gamma.
#
if not (IsGraph(gamma) and IsInt(v)) then
   Error("usage: VertexName( <Graph>, <Int> )");
fi;
if IsBound(gamma.names) then 
   return Copy(gamma.names[v]);
else 
   return v;
fi;
end;

VertexDegree := function(gamma,v)
#
# Returns the vertex (out)degree of vertex  v  in the graph  gamma.
#
local rw,sch;
if not IsGraph(gamma) or not IsInt(v) then
   Error("usage: VertexDegree( <Graph>, <Int> )");
fi;
if v<1 or v>gamma.order then
   Error("<v> is not a vertex of <gamma>");
fi;
sch:=gamma.schreierVector;
rw:=RepWord(gamma.group.generators,sch,v);
return Length(gamma.adjacencies[-sch[rw.representative]]); 
end;

VertexDegrees := function(gamma)
#
# Returns the set of vertex (out)degrees for the graph  gamma.
#
local adj,degs;
if not IsGraph(gamma) then
   Error("usage: VertexDegrees( <Graph> )");
fi;
degs:=Set([]);
for adj in gamma.adjacencies do
   AddSet(degs,Length(adj));
od;
return degs;
end;

IsVertexPairEdge := function(gamma,x,y)
#
# Assuming that  x,y  are vertices of  gamma,  returns true
# iff  [x,y]  is an edge of  gamma.
#
local w,sch,gens;
sch:=gamma.schreierVector;
gens:=gamma.group.generators;
w:=sch[x];
while w > 0 do
   x:=x/gens[w];
   y:=y/gens[w];
   w:=sch[x];
od;
return y in gamma.adjacencies[-w];
end;

IsEdge := function(gamma,e)
#
# Returns  true  iff  e  is an edge of  gamma.
#
if not IsGraph(gamma) then 
   Error("usage: IsEdge( <Graph>, <obj> )");
fi;
if not IsList(e) or Length(e)<>2 or not IsVertex(gamma,e[1])
		 or not IsVertex(gamma,e[2]) then
   return false;
fi;
return IsVertexPairEdge(gamma,e[1],e[2]);
end;

Adjacency := function(gamma,v)
#
# Returns (a copy of) the set of vertices of  gamma  adjacent to vertex  v.
#
local w,adj,rw,gens,sch;
sch:=gamma.schreierVector;
if sch[v] < 0 then 
   return Copy(gamma.adjacencies[-sch[v]]);
fi;
gens:=gamma.group.generators;
rw:=RepWord(gens,sch,v);
adj:=gamma.adjacencies[-sch[rw.representative]]; 
for w in rw.word do 
   adj:=OnTuples(adj,gens[w]); 
od;
return Set(adj);
end;

IsSimpleGraph := function(gamma)
#
# Returns  true  iff graph  gamma  is simple (i.e. has no loops and 
# if [x,y] is an edge then so is [y,x]).  Also sets the isSimple 
# field of  gamma  if this field was not already bound.
#
local adj,i,x,H,orb;
if not IsGraph(gamma) then 
   Error("usage: IsSimpleGraph( <Graph> )");
fi;
if IsBound(gamma.isSimple) then
   return gamma.isSimple;
fi;
for i in [1..Length(gamma.adjacencies)] do
   adj:=gamma.adjacencies[i];
   x:=gamma.representatives[i];
   if x in adj then    # a loop exists
      gamma.isSimple:=false;
      return false;
   fi;
   H:=ProbablyStabilizer(gamma.group,x);
   for orb in Orbits(H,adj) do
      if not IsVertexPairEdge(gamma,orb[1],x) then
	 gamma.isSimple:=false;
	 return false;
      fi;
   od;
od;
gamma.isSimple:=true;
return true;
end;

DirectedEdges := function(gamma)
#
# Returns the set of directed (ordered) edges of  gamma.
#
local i,j,edges;
if not IsGraph(gamma) then 
   Error("usage: DirectedEdges( <Graph> )");
fi;
edges:=[];
for i in [1..gamma.order] do
   for j in Adjacency(gamma,i) do
      Add(edges,[i,j]);
   od;
od;
IsSet(edges);  # edges is a set.
return edges;
end;

UndirectedEdges := function(gamma)
#
# Returns the set of undirected edges of  gamma,  which must be 
# a simple graph.
#
local i,j,edges;
if not IsGraph(gamma) then 
   Error("usage: UndirectedEdges( <Graph> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> must be a simple graph");
fi;
edges:=[];
for i in [1..gamma.order-1] do
   for j in Adjacency(gamma,i) do
      if i<j then 
	 Add(edges,[i,j]);
      fi;
   od;
od;
IsSet(edges);  # edges is a set.
return edges;
end;

AddEdgeOrbit := function(arg) 
#
# Let  gamma=arg[1]  and  e=arg[2].
# If  arg[3]  is bound then it is assumed to be  Stabilizer(gamma.group,e[1]).
# This procedure adds edge orbit  e^gamma.group  to the edge set of  gamma.
#
local w,word,sch,gens,gamma,e,x,y,orb,u,v;
gamma:=arg[1]; 
e:=arg[2];
if not IsGraph(gamma) or not IsList(e) 
    or (IsBound(arg[3]) and not IsPermGroup(arg[3])) then
   Error("usage: AddEdgeOrbit( <Graph>, <List>, [, <PermGroup> ] )");
fi;
if Length(e)<>2 or not IsVertex(gamma,e[1]) or not IsVertex(gamma,e[2]) then
   Error("invalid <e>");
fi;
sch:=gamma.schreierVector;
gens:=gamma.group.generators;
x:=e[1];
y:=e[2];
w:=sch[x];
word:=[];
while w > 0 do
   Add(word,w);
   x:=x/gens[w];
   y:=y/gens[w];
   w:=sch[x];
od;
if not(y in gamma.adjacencies[-sch[x]]) then
   #  e  is not an edge of  gamma
   if not IsBound(arg[3]) then
      orb:=Orbit(Stabilizer(gamma.group,x),y);
   else
      if ForAny(arg[3].generators,x->e[1]^x<>e[1]) then
	 Error("<arg[3]>  not equal to  Stabilizer(<gamma.group>,<e[1]>)");
      fi;
      orb:=[];
      for u in Orbit(arg[3],e[2]) do
	 v:=u;
	 for w in word do
	    v:=v/gens[w];
	 od;
	 Add(orb,v);
      od;
   fi;
   UniteSet(gamma.adjacencies[-sch[x]],orb);
   if e[1]=e[2] then
      gamma.isSimple:=false;
   elif IsBound(gamma.isSimple) and gamma.isSimple then
      if not IsVertexPairEdge(gamma,e[2],e[1]) then 
	 gamma.isSimple:=false; 
      fi;
   else 
      Unbind(gamma.isSimple);
   fi;
   Unbind(gamma.autGroup);
   Unbind(gamma.canonicalLabelling);
fi;   
end;

RemoveEdgeOrbit := function(arg) 
#
# Let  gamma=arg[1]  and  e=arg[2].
# If  arg[3]  is bound then it is assumed to be  Stabilizer(gamma.group,e[1]).
# This procedure removes the edge orbit  e^gamma.group  from the edge set 
# of  gamma, if this orbit exists, and otherwise does nothing.
#
local w,word,sch,gens,gamma,e,x,y,orb,u,v;
gamma:=arg[1]; 
e:=arg[2];
if not IsGraph(gamma) or not IsList(e) 
    or (IsBound(arg[3]) and not IsPermGroup(arg[3])) then
   Error("usage: RemoveEdgeOrbit( <Graph>, <List>, [, <PermGroup> ] )");
fi;
if Length(e)<>2 or not IsVertex(gamma,e[1]) or not IsVertex(gamma,e[2]) then
   Error("invalid <e>");
fi;
sch:=gamma.schreierVector;
gens:=gamma.group.generators;
x:=e[1];
y:=e[2];
w:=sch[x];
word:=[];
while w > 0 do
   Add(word,w);
   x:=x/gens[w];
   y:=y/gens[w];
   w:=sch[x];
od;
if y in gamma.adjacencies[-sch[x]] then
   #  e  is an edge of  gamma
   if not IsBound(arg[3]) then
      orb:=Orbit(Stabilizer(gamma.group,x),y);
   else
      if ForAny(arg[3].generators,x->e[1]^x<>e[1]) then
	 Error("<arg[3]>  not equal to  Stabilizer(<gamma.group>,<e[1]>)");
      fi;
      orb:=[];
      for u in Orbit(arg[3],e[2]) do
	 v:=u;
	 for w in word do
	    v:=v/gens[w];
	 od;
	 Add(orb,v);
      od;
   fi;
   SubtractSet(gamma.adjacencies[-sch[x]],orb);
   if IsBound(gamma.isSimple) and gamma.isSimple then
      if IsVertexPairEdge(gamma,e[2],e[1]) then 
	 gamma.isSimple:=false; 
      fi;
   else 
      Unbind(gamma.isSimple);
   fi;
   Unbind(gamma.autGroup);
   Unbind(gamma.canonicalLabelling);
fi;   
end;

EdgeOrbitsGraph := function(arg)
#
# Let  G=arg[1],  E=arg[2].
# Returns the (directed) graph with vertex set {1,...,n} and edge set 
# the union over  e in E  of  e^G,  where  n=arg[3]  if  arg[3]  is bound,
# and  n=MaximumMovedPoint(G.generators)  otherwise.
# (E can consist of just a singleton edge.)
#
local G,E,n,gamma,e;
G:=arg[1];
E:=arg[2];
if IsInt(E[1]) then   # assume  E  consists of a single edge.
   E:=[E];
fi;
if IsBound(arg[3]) then 
   n:=arg[3];
else
   n:=MaximumMovedPoint(G.generators);
fi;
if not IsPermGroup(G) or not IsList(E) or not IsInt(n) then
   Error("usage: EdgeOrbitsGraph( <PermGroup>, <List> [, <Int> ] )");
fi;
gamma:=NullGraph(G,n);
for e in E do 
   AddEdgeOrbit(gamma,e); 
od;
return gamma;
end;

VertexColouring := function(gamma)
#
# Returns a proper vertex-colouring for the graph  gamma,  which must be
# simple.
#
# Here a (proper) vertex-colouring is a list  C  of natural numbers,
# of length  gamma.order,  such that  C[i]<>C[j]  whenever 
# [i,j]  is an edge of  gamma.
#
# At present a greedy algorithm is used.  
#
local i,j,g,c,C,orb,a,adj,adjs,adjcolours,maxcolour,im,gens;
if not IsGraph(gamma) then 
   Error("usage: VertexColouring( <Graph> )");
fi;
if gamma.order=0 then
   return [];
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
C:=List([1..gamma.order],x->0);
maxcolour:=0;
gens:=gamma.group.generators;
for i in [1..Length(gamma.representatives)] do
   orb:=[gamma.representatives[i]];
   adjs:=[];
   adj:=gamma.adjacencies[i];
   adjs[orb[1]]:=adj;
   # colour vertex  orb[1]
   adjcolours:=BlistList([1..maxcolour+1],[]);
   for a in adj do
      if C[a]>0 then
	 adjcolours[C[a]]:=true;
      fi;
   od;
   c:=1;
   while adjcolours[c] do
      c:=c+1;
   od;
   C[orb[1]]:=c;
   if c>maxcolour then
      maxcolour:=c;
   fi;
   for j in orb do 
      for g in gens do
	 im:=j^g;
	 if C[im]=0 then 
	    Add(orb,im);
	    adj:=OnTuples(adjs[j],g);
	    adjs[im]:=adj;
	    # colour vertex  im
	    adjcolours:=BlistList([1..maxcolour+1],[]);
	    for a in adj do
	       if C[a]>0 then
		  adjcolours[C[a]]:=true;
	       fi;
	    od;
	    c:=1;
	    while adjcolours[c] do
	       c:=c+1;
	    od;
	    C[im]:=c;
	    if c>maxcolour then
	       maxcolour:=c;
	    fi;
	 fi; 
      od; 
      Unbind(adjs[j]);   
   od; 
od;
return C;
end;

CollapsedAdjacencyMat := function(arg)
#
# Returns the collapsed adjacency matrix  A  for  gamma=arg[2]  wrt  
# group  G=arg[1],  assuming  G <= Aut(gamma). 
# The rows and columns of  A  are indexed by the orbits 
# orbs[1],...,orbs[n], say, of  G  on the vertices of  
# gamma, and the entry  A[i][j]  of  A  is defined as follows:
#    Let  reps[i]  be a representative of the  i-th  G-orbit  orbs[i].
#    Then  A[i][j] equals the number of neighbours (in  gamma)
#    of  reps[i]  in  orbs[j]. 
# Note that this definition does not depend on the choice of 
# representative  reps[i].
#
# *** New for Grape 2.3: In the special case where this function 
# is given just one argument, then we must have  gamma=arg[1], 
# we must have  gamma.group  transitive on the vertices of  gamma,
# and then the returned collapsed adjacency matrix for  gamma  is 
# w.r.t. the stabilizer in  gamma.group  of  1.  Additionally 
# [1]=orbs[1].  This feature is to 
# conform with the definition of collapsed adjacency matrix in 
# Praeger and Soicher, "Low Rank Representations and Graphs for 
# Sporadic Groups", CUP, Cambridge, 1997.  (In GRAPE we allow a collapsed
# adjacency matrix to be more general, as we can collapse w.r.t. to
# an arbitrary subgroup of  Aut(gamma),  and  gamma  need not 
# even be vertex-transitive.)  
#
local G,gamma,orbs,i,j,n,A,orbnum,reps;
if Length(arg)=1 then
   gamma:=arg[1];
   if not IsGraph(gamma) then
      Error("usage: CollapsedAdjacencyMat( [<PermGroup>,] <Graph>)");
   fi;
   if gamma.order=0 then 
      return []; 
   fi; 
   if not IsTransitive(gamma.group,[1..gamma.order]) then
      Error(
       "<gamma.group> not transitive on vertices of single argument <gamma>"); 
   fi;
   G := Stabilizer(gamma.group,1);
else   
   G := arg[1];
   gamma := arg[2];
   if not IsPermGroup(G) or not IsGraph(gamma) then
      Error("usage: CollapsedAdjacencyMat( [<PermGroup>,] <Graph> )");
   fi;
   if gamma.order=0 then
      return [];
   fi;
fi;
orbs:=OrbitNumbers(G,gamma.order);
orbnum:=orbs.orbitNumbers;
reps:=orbs.representatives;
n:=Length(reps);
A:=NullMat(n,n);
for i in [1..n] do 
   for j in Adjacency(gamma,reps[i]) do
      A[i][orbnum[j]]:=A[i][orbnum[j]]+1;
   od;
od;
return A;
end;

OrbitalGraphColadjMats := function(arg)
#
# This function returns a sequence of collapsed adjacency 
# matrices for the the orbital graphs of the (assumed) transitive  
# G=arg[1].  The matrices are collapsed w.r.t.  Stabilizer(G,1),  so
# that these are collapsed adjacency matrices in the sense of 
# Praeger and Soicher, "Low Rank Representations and Graphs for 
# Sporadic Groups", CUP, Cambridge, 1997. 
# The matrices are collapsed w.r.t. a fixed ordering of the G-suborbits,
# with the trivial suborbit  [1]  coming first.
#
# If  arg[2]  is bound, then it is assumed to be  Stabilizer(G,1).
#
# An alternative name for this function is  
# OrbitalGraphIntersectionMatrices, 
# for backwards compatibility with GRAPE 2.2.
# Note that the sequence of matrices returned by this function
# are closely related to, but not necessarily exactly the 
# same as, the intersection matrices for the
# association scheme whose classes are the non-trivial
# G-orbitals, as defined in Brouwer, Cohen and Neumaier,
# "Distance-regular Graphs", Springer, Berlin/New York, 1989.
#
local G,H,orbs,deg,i,j,k,n,intmats,A,orbnum,reps,gamma;
G:=arg[1];
if not IsPermGroup(G) or (IsBound(arg[2]) and not IsPermGroup(arg[2])) then
   Error("usage: OrbitalGraphColadjMats( <PermGroup> [, <PermGroup> ] )");
fi;
if IsBound(arg[2]) then
   H:=arg[2];
   if ForAny(H.generators,x->1^x<>1) then
      Error("<H> does not fix the point 1");
   fi;
else
   H:=Stabilizer(G,1);
fi;
deg:=Maximum(MaximumMovedPoint(G.generators),1);
if not IsTransitive(G,[1..deg]) then
   Error("<G> not transitive");
fi;
gamma:=NullGraph(G,deg);
orbs:=OrbitNumbers(H,gamma.order);
orbnum:=orbs.orbitNumbers;
reps:=orbs.representatives;
if reps[1]<>1 then # this cannot happen!
   Error("internal error");
fi;
n:=Length(reps);
intmats:=[];
for i in [1..n] do
   AddEdgeOrbit(gamma,[1,reps[i]],H);
   A:=NullMat(n,n);
   for j in [1..n] do 
      for k in Adjacency(gamma,reps[j]) do
	 A[j][orbnum[k]]:=A[j][orbnum[k]]+1;
      od;
   od;
   intmats[i]:=A;
   if i < n then
      RemoveEdgeOrbit(gamma,[1,reps[i]],H);
   fi;
od;
return intmats;
end;

OrbitalGraphIntersectionMatrices := OrbitalGraphColadjMats;

LocalInfo := function(arg)
#
# Calculates  "local info"  for  gamma=arg[1]  from point of view of vertex  
# set (or list or singleton vertex)  V=arg[2].
#
# Returns record containing the  "layer numbers"  for gamma w.r.t.
# V,  as well as the the local diameter and local girth of gamma w.r.t.  V.
# ( layerNumbers[i]=j>0 if vertex  i  is in layer[j]  (i.e. at distance
# j-1  from  V),  layerNumbers[i]=0  if vertex  i
# is not joined by a path to some element of  V.)
# Also, if a local parameter  ci[V], ai[V], or bi[V]  exists then
# this information is recorded in  localParameters[i+1][1], 
# localParameters[i+1][2], or localParameters[i+1][3], 
# respectively (otherwise a -1 is recorded). 
#
# *** If  gamma  is not simple then local girth and the local parameters
# may not be what you think. The local girth has no real meaning if 
# |V| > 1.
#
# *** But note: If arg[3] is bound and arg[3] > 0
# then the procedure stops after  layers[arg[3]]  has been determined.
# If  arg[4]  is bound (a set or list or singleton vertex), then the
# procedure stops when the first layer containing a vertex in  arg[4]  is 
# complete.
# Also, if  arg[4]  is bound then the local info record contains 
# a distance field whose value is  min{ d(v,w)) | v in V, w in arg[4] }.
#
# If  arg[5]  is bound then it is assumed to be a subgroup of Aut(gamma)
# stabilising  V  setwise.
#
local gamma,V,layers,localDiameter,localGirth,localParameters,i,j,x,y,next,
      nprev,nhere,nnext,sum,orbs,orbnum,laynum,lnum,
      stoplayer,stopvertices,distance,loc,reps,layerNumbers;
gamma:=arg[1];
V:=arg[2];
if IsInt(V) then 
   V:=[V];
fi;
if not IsGraph(gamma) or not IsList(V) then 
   Error("usage: LocalInfo( <Graph>, <Int> or <List>, ... )");
fi;
if not IsSet(V) then
   V:=Set(V);
fi;
if V=Set([]) or not IsSubset([1..gamma.order],V) then 
   Error("<V> must be non-empty set of vertices of <gamma>");
fi;
if IsBound(arg[3]) then 
   stoplayer:=arg[3]; 
   if not IsInt(stoplayer) or stoplayer < 0 then
      Error("<stoplayer> must be integer >= 0");
   fi;
else 
   stoplayer:=0; 
fi;
if IsBound(arg[4]) then 
   stopvertices:=arg[4]; 
   if IsInt(stopvertices) then
      stopvertices:=[stopvertices];
   fi;
   if not IsSet(stopvertices) then
      stopvertices:=Set(stopvertices);
   fi;
   if not IsSubset([1..gamma.order],stopvertices) 
      then
	 Error("<stopvertices> must be a set of vertices of <gamma>");
   fi;
else 
   stopvertices:=Set([]); 
fi;
if IsBound(arg[5]) then 
   if not IsPermGroup(arg[5]) then 
      Error("<arg[5]> must be a permutation group (<= Stab(<V>)");
   fi;
   orbs:=OrbitNumbers(arg[5],gamma.order);
else
   if Length(V)=1 then
      if IsBound(gamma.autGroup) then
	 orbs:=ProbablyStabilizerOrbitNumbers(gamma.autGroup,V[1],gamma.order);
      else
	 orbs:=ProbablyStabilizerOrbitNumbers(gamma.group,V[1],gamma.order);
      fi;
   else
      orbs:=rec(representatives:=[1..gamma.order],
		orbitNumbers:=[1..gamma.order]); 
   fi;
fi;
orbnum:=orbs.orbitNumbers;
reps:=orbs.representatives;
laynum:=[];
for i in [1..Length(reps)] do 
   laynum[i]:=0;
od;
localGirth:=-1; 
distance:=-1;
localParameters:=[]; 
next:=Set([]);
for i in V do
   AddSet(next,orbnum[i]);
od;
sum:=Length(next);
for i in next do 
   laynum[i]:=1;
od;
layers:=[]; 
layers[1]:=next; 
i:=1; 
if Length(Intersection(V,stopvertices)) > 0 then
   stoplayer:=1; 
   distance:=0;
fi;
while stoplayer<>i and Length(next)>0 do
   next:=Set([]);
   for x in layers[i] do 
      nprev:=0; 
      nhere:=0; 
      nnext:=0;
      for y in Adjacency(gamma,reps[x]) do
	 lnum:=laynum[orbnum[y]];
	 if i>1 and lnum=i-1 then 
	    nprev:=nprev+1;
	 elif lnum=i then 
	    nhere:=nhere+1;
	 elif lnum=i+1 then
	    nnext:=nnext+1;
	 elif lnum=0 then
	    AddSet(next,orbnum[y]); 
	    nnext:=nnext+1;
	    laynum[orbnum[y]]:=i+1;
	 fi;
      od;
      if (localGirth=-1 or localGirth=2*i-1) and nprev>1 then 
	 localGirth:=2*(i-1); 
      fi;
      if localGirth=-1 and nhere>0 then 
	 localGirth:=2*i-1; 
      fi;
      if not IsBound(localParameters[i]) then 
	 localParameters[i]:=[nprev,nhere,nnext];
      else
	 if nprev<>localParameters[i][1] then 
	    localParameters[i][1]:=-1; 
	 fi;
	 if nhere<>localParameters[i][2] then 
	    localParameters[i][2]:=-1; 
	 fi;
	 if nnext<>localParameters[i][3] then 
	    localParameters[i][3]:=-1; 
	 fi;
      fi;
   od;
   if Length(next)>0 then 
      i:=i+1; 
      layers[i]:=next; 
      for j in stopvertices do 
	 if laynum[orbnum[j]]=i then 
	    stoplayer:=i; 
	    distance:=i-1;
	 fi;
      od;
      sum:=sum+Length(next); 
   fi;
od;
if sum=Length(reps) then 
   localDiameter:=Length(layers)-1; 
else 
   localDiameter:=-1; 
fi;
# now change  orbnum  to give the layer numbers instead of orbit numbers.
layerNumbers:=orbnum;
for i in [1..gamma.order] do
   layerNumbers[i]:=laynum[orbnum[i]];
od;
loc:=rec(layerNumbers:=layerNumbers,localDiameter:=localDiameter,
	 localGirth:=localGirth,localParameters:=localParameters);
if Length(stopvertices) > 0 then 
   loc.distance:=distance;
fi;
return loc;
end;

LocalInfoMat := function(A,rows)
#
# Calculates local info on a graph using a collapsed adjacency matrix 
#  A  for that graph.
# This local info is from the point of view of the set of vertices
# represented by the set  rows  of row indices of  A.
# The elements of  layers[i]  will be the row indices representing 
# the vertices of the i-th layer.  
# No  distance  field will be calculated.
#
# *** If  A  is not the collapsed adjacency matrix for a simple graph 
# then  localGirth  and localParameters may not be what you think.
# If  rows  does not represent a single vertex then  localGirth  has 
# no real meaning. 
#
local layers,localDiameter,localGirth,localParameters,i,j,x,y,next,
      nprev,nhere,nnext,sum,laynum,lnum,n;
if IsInt(rows) then 
   rows:=[rows];
fi;
if not IsMat(A) or not IsList(rows) then 
   Error("usage: LocalInfoMat( <Mat>, <Int> or <List> )");
fi;
if not IsSet(rows) then
   rows:=Set(rows);
fi;
n:=Length(A);
if rows=Set([]) or not IsSubset([1..n],rows) then 
   Error("<rows> must be non-empty set of row indices");
fi;
laynum:=[];
for i in [1..n] do
   laynum[i]:=0;
od;
localGirth:=-1; 
localParameters:=[]; 
next:=Set([]);
for i in rows do
   AddSet(next,i);
od;
for i in next do 
   laynum[i]:=1;
od;
layers:=[]; 
layers[1]:=next; 
i:=1; 
sum:=Length(rows);
while Length(next)>0 do
   next:=Set([]);
   for x in layers[i] do 
      nprev:=0; 
      nhere:=0; 
      nnext:=0;
      for y in [1..n] do
	 j:=A[x][y];
	 if j>0 then
	    lnum:=laynum[y];
	    if i>1 and lnum=i-1 then 
	       nprev:=nprev+j;
	    elif lnum=i then 
	       nhere:=nhere+j;
	    elif lnum=i+1 then
	       nnext:=nnext+j;
	    elif lnum=0 then
	       AddSet(next,y); 
	       nnext:=nnext+j;
	       laynum[y]:=i+1;
	    fi;
	 fi;
      od;
      if (localGirth=-1 or localGirth=2*i-1) and nprev>1 then 
	 localGirth:=2*(i-1); 
      fi;
      if localGirth=-1 and nhere>0 then 
	 localGirth:=2*i-1; 
      fi;
      if not IsBound(localParameters[i]) then 
	 localParameters[i]:=[nprev,nhere,nnext];
      else
	 if nprev<>localParameters[i][1] then 
	    localParameters[i][1]:=-1; 
	 fi;
	 if nhere<>localParameters[i][2] then 
	    localParameters[i][2]:=-1; 
	 fi;
	 if nnext<>localParameters[i][3] then 
	    localParameters[i][3]:=-1; 
	 fi;
      fi;
   od;
   if Length(next)>0 then 
      i:=i+1; 
      layers[i]:=next; 
      sum:=sum+Length(next); 
   fi;
od;
if sum=n then 
   localDiameter:=Length(layers)-1; 
else 
   localDiameter:=-1; 
fi;
return rec(layerNumbers:=laynum,localDiameter:=localDiameter,
	   localGirth:=localGirth,localParameters:=localParameters);
end;

InducedSubgraph := function(arg) 
#
# Returns the subgraph of  gamma=arg[1]  induced on the vertex list  V=arg[2].
# If  arg[3]  is unbound, then the trivial group is the group associated 
# with the returned induced subgraph. 
# If  arg[3]  is bound, this function assumes that  G=arg[3]   fixes  
# V  setwise, and is a group of automorphisms of the induced subgraph 
# when restriced to  V.  In this case, the image of  G  acting on  V  is 
# the group associated with the returned induced subgraph.
#
local gamma,V,G,indu,i,j,W,VV,X,x,gens;
gamma:=arg[1];
V:=arg[2];
if not IsGraph(gamma) or not IsList(V) 
    or (IsBound(arg[3]) and not IsPermGroup(arg[3])) then
   Error("usage: InducedSubgraph( <Graph>, <List>, [, <PermGroup> ] )");
fi;
if not IsSubset([1..gamma.order],V) then
   Error("<V> must be a list of vertices of <gamma>");
fi;
if IsBound(arg[3]) then 
   G:=arg[3];
else
   G:=Group(());
fi;
W:=[];
for i in [1..Length(V)] do 
   W[V[i]]:=i; 
od;
gens:=[]; 
for i in [1..Length(G.generators)] do 
   gens[i]:=[];
   for j in V do 
      gens[i][W[j]]:=W[j^G.generators[i]]; 
   od;
   gens[i]:=PermList(gens[i]);
od;
indu:=NullGraph(Group(gens,()),Length(V));
if IsBound(gamma.isSimple) and gamma.isSimple then 
   indu.isSimple:=true;
else 
   Unbind(indu.isSimple); 
fi;
VV:=Set(V);
for i in [1..Length(indu.representatives)] do
   X:=Intersection(VV,Adjacency(gamma,V[indu.representatives[i]]));
   for x in X do 
      AddSet(indu.adjacencies[i],W[x]); 
   od;
od;
indu.names:=[];
if IsBound(gamma.names) then 
   for i in V do 
      indu.names[W[i]]:=gamma.names[i]; 
   od;
else 
   for i in V do 
      indu.names[W[i]]:=i; 
   od;
fi;
return indu;
end;

Distance := function(arg)
#
# Let  gamma=arg[1],  X=arg[2],  Y=arg[3].
# Returns the distance  d(X,Y)  in the graph  gamma, where  X,Y 
# are singleton vertices or lists of vertices.
# (Returns  -1  if no (directed) path joins  X  to  Y  in  gamma.)
# If  arg[4]  is bound, then it is assumed to be a subgroup
# of  Aut(gamma)  stabilizing  X  setwise.
#
local gamma,X,Y;
gamma:=arg[1];
X:=arg[2];
if IsInt(X) then 
   X:=[X];
fi;
Y:=arg[3];
if IsInt(Y) then
   Y:=[Y];
fi;
if not (IsGraph(gamma) and IsList(X) and IsList(Y)) then 
   Error(ConcatenationString("usage: Distance( <Graph>, <Int> or <List>, ",
			     "<Int> or <List> [, <PermGroup> ] )"));
fi;
if IsBound(arg[4]) then 
   return LocalInfo(gamma,X,0,Y,arg[4]).distance;
else
   return LocalInfo(gamma,X,0,Y).distance;
fi;
end;

Diameter := function(gamma)
#
# Returns the diameter of  gamma. 
# A diameter of  -1  means that gamma is not (strongly) connected.  
#
local r,d,loc;
if not IsGraph(gamma) then 
   Error("usage: Diameter( <Graph> )");
fi;
if gamma.order=0 then
   Error("<gamma> has no vertices");
fi;
d:=-1;
for r in gamma.representatives do 
   loc:=LocalInfo(gamma,r);
   if loc.localDiameter=-1 then
      return -1; 
   fi;
   if loc.localDiameter > d then
      d:=loc.localDiameter;
   fi;
od;
return d;
end;

Girth := function(gamma)
#
# Returns the girth of  gamma,  which must be a simple graph. 
# A girth of  -1  means that gamma is a forest.  
#
local r,g,locgirth,maxlook,adj;
if not IsGraph(gamma) then 
   Error("usage: Girth( <Graph> )");
fi;
if gamma.order=0 then
   return -1;
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
adj:=gamma.adjacencies[1];
if adj<>[] and Intersection(adj,Adjacency(gamma,adj[1]))<>[] then
   return 3;
fi;
g:=-1;
maxlook:=0;
for r in gamma.representatives do 
   locgirth:=LocalInfo(gamma,r,maxlook).localGirth;
   if locgirth=3 then 
      return 3;
   fi;
   if locgirth<>-1 then
      if g=-1 or locgirth<g then
	  g:=locgirth;
	  maxlook:=Int((g+1)/2);
      fi;
   fi;
od;
return g;
end;

IsRegularGraph := function(gamma)
#
# Returns  true  iff the graph  gamma  is (out)regular.
#
local deg,i;
if not IsGraph(gamma) then 
   Error("usage: IsRegularGraph( <Graph> )");
fi;
if gamma.order=0 then 
   return true;
fi;
deg:=Length(gamma.adjacencies[1]);
for i in [2..Length(gamma.adjacencies)] do
   if deg <> Length(gamma.adjacencies[i]) then 
      return false;
   fi;
od;
return true;
end;

IsNullGraph := function(gamma)
#
# Returns  true  iff the graph  gamma  has no edges.
#
local i;
if not IsGraph(gamma) then 
   Error("usage: IsNullGraph( <Graph> )");
fi;
for i in [1..Length(gamma.adjacencies)] do
   if Length(gamma.adjacencies[i])<>0 then 
      return false;
   fi;
od;
return true;
end;

IsCompleteGraph := function(arg)
#
# Returns  true  iff the graph  gamma=arg[1]  is a complete graph.
# arg[2] is true iff all loops must exist for  gamma  to be considered
# a complete graph (default: false), otherwise loops are ignored (except to 
# possibly set  gamma.isSimple). 
#
local deg,i,notnecsimple,gamma,mustloops;
gamma:=arg[1];
if IsBound(arg[2]) then 
   mustloops := arg[2];
else
   mustloops := false;
fi;
if not IsGraph(gamma) or not IsBool(mustloops) then 
   Error("usage: IsCompleteGraph( <Graph> [, <Bool> ] )");
fi;
notnecsimple := not IsBound(gamma.isSimple) or not gamma.isSimple;
for i in [1..Length(gamma.adjacencies)] do
   deg := Length(gamma.adjacencies[i]);
   if deg < gamma.order-1 then 
      return false;
   fi;
   if deg=gamma.order-1 then
      if mustloops then
	 return false;
      fi;
      if notnecsimple and 
       (gamma.representatives[i] in gamma.adjacencies[i]) then
	 gamma.isSimple := false;
	 return false;
      fi;
   fi;
od;
return true;
end;

IsLoopy := function(gamma)
#
# Returns  true  iff graph  gamma  has a loop.
#
local i;
if not IsGraph(gamma) then 
   Error("usage: IsLoopy( <Graph> )");
fi;
for i in [1..Length(gamma.adjacencies)] do
   if gamma.representatives[i] in gamma.adjacencies[i] then
      gamma.isSimple := false;
      return true;
   fi;
od;
return false;
end;

IsConnectedGraph := function(gamma)
#
# Returns true iff  gamma  is (strongly) connected.
#
if not IsGraph(gamma) then 
   Error("usage: IsConnectedGraph( <Graph> )");
fi;
if gamma.order=0 then
   return true;
fi;
if IsSimpleGraph(gamma) then
   return LocalInfo(gamma,1).localDiameter > -1;
else
   return Diameter(gamma) > -1;
fi;
end;

ConnectedComponent := function(gamma,v)
#
# Returns the set of all vertices in  gamma  which can be reached by 
# a path starting at vertex  v.  The graph  gamma  must be simple.
#
local comp,j,laynum;
if not IsGraph(gamma) or not IsInt(v) then 
   Error("usage: ConnectedComponent( <Graph>, <Int> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
laynum:=LocalInfo(gamma,v).layerNumbers;
comp:=[];
for j in [1..gamma.order] do
   if laynum[j] > 0 then 
      Add(comp,j);
   fi; 
od;
IsSet(comp);
return comp;
end;

ConnectedComponents := function(gamma)
#
# Returns a list of the vertex sets of the connected components
# of  gamma,  which must be a simple graph.
#
local comp,used,i,j,x,cmp,laynum;
if not IsGraph(gamma) then 
   Error("usage: ConnectedComponents( <Graph> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
comp:=[]; 
used:=BlistList([1..gamma.order],[]);
for i in [1..gamma.order] do 
   if not used[i] then   # new component
      cmp:=[];
      laynum:=LocalInfo(gamma,i).layerNumbers;
      for j in [1..gamma.order] do
	 if laynum[j] > 0 then 
	    Add(cmp,j);
	 fi;
      od;
      IsSet(cmp);
      for x in Orbit(gamma.group,cmp,OnSets) do
	 Add(comp,x);
	 for j in x do 
	    used[j]:=true;
	 od;
      od;
   fi; 
od;
return comp;
end;

ComponentLocalInfos := function(gamma)
#
# Returns a sequence of localinfos for the connected components of  
# gamma  (w.r.t. some vertex in each component).
# The graph  gamma  must be simple.
#
local comp,used,i,j,k,laynum;
if not IsGraph(gamma) then 
   Error("usage: ComponentLocalInfos( <Graph> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
comp:=[]; 
used:=BlistList([1..gamma.order],[]);
k:=0;
for i in [1..gamma.order] do 
   if not used[i] then   # new component
      k:=k+1; 
      comp[k]:=LocalInfo(gamma,i);
      laynum:=comp[k].layerNumbers;
      for j in [1..gamma.order] do
	 if laynum[j] > 0 then
	    used[j]:=true;
	 fi;
      od;
   fi; 
od;
return comp;
end;

Bicomponents := function(gamma)
#
# If  gamma  is bipartite, returns a length 2 list of
# bicomponents, or parts, of  gamma,  else returns the empty list.
# *** This function is for simple  gamma  only.
#
# Note: if gamma.order=0 this function returns [[],[]], and if 
# gamma.order=1 this function returns [[],[1]] (unlike GRAPE 2.2
# which returned [], which was inconsistent with considering 
# a zero vertex graph to be bipartite).
#
local bicomps,i,lnum,loc,locs;
if not IsGraph(gamma) then 
   Error("usage: Bicomponents( <Graph> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
bicomps:=[Set([]),Set([])]; 
if gamma.order=0 then 
   return bicomps;
fi;
if IsNullGraph(gamma) then 
   return [Set([1..gamma.order-1]),Set([gamma.order])];
fi;
locs:=ComponentLocalInfos(gamma);
for loc in locs do
   for i in [2..Length(loc.localParameters)] do
      if loc.localParameters[i][2]<>0 then
	 return [];
      fi;
   od;
   for i in [1..Length(loc.layerNumbers)] do 
      lnum:=loc.layerNumbers[i];
      if lnum>0 then
	 if lnum mod 2 = 1 then
	    AddSet(bicomps[1],i);
	 else
	    AddSet(bicomps[2],i);
	 fi;
      fi; 
   od;
od;
return bicomps;
end;

IsBipartite := function(gamma)
#
# Returns  true  iff  gamma  is bipartite. 
# *** This function is only for simple  gamma.
#
# Note: Now the one vertex graph is considered to be bipartite 
# (as well as the zero vertex graph). This is a change from the inconsistent 
# GRAPE 2.2 view that a zero vertex graph is bipartite, but not a one 
# vertex graph.
#
if not IsGraph(gamma) then 
   Error("usage: IsBipartite( <Graph> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
return Length(Bicomponents(gamma))=2;
end;

Layers := function(arg)
#
# Returns the list of vertex layers of  gamma=arg[1],  
# starting from  V=arg[2],  which may be a vertex list or singleton vertex. 
# Layers[i]  is the set of vertices at distance  i-1  from  V.
# If  arg[3]  is bound then it is assumed to be a subgroup 
# of  Aut(gamma)  stabilizing  V  setwise.
#
local gamma,V;
gamma:=arg[1];
V:=arg[2];
if IsInt(V) then
   V:=[V];
fi;
if not (IsGraph(gamma) and IsList(V)) then 
   Error("usage: Layers( <Graph>, <Int> or <List>, [, <PermGroup>] )");
fi;
if IsBound(arg[3]) then
   return NumbersToSets(LocalInfo(gamma,V,0,[],arg[3]).layerNumbers);
else
   return NumbersToSets(LocalInfo(gamma,V).layerNumbers);
fi;
end;

LocalParameters := function(arg)
#
# Returns the local parameters of simple, connected  gamma=arg[1],  
# w.r.t to vertex list (or singleton vertex)  V=arg[2].
# The nonexistence of a local parameter is denoted by  -1.
# If  arg[3]  is bound then it is assumed to be a subgroup 
# of  Aut(gamma)  stabilizing  V  setwise.
#
local gamma,V,loc;
gamma:=arg[1];
V:=arg[2];
if IsInt(V) then
   V:=[V];
fi;
if not IsGraph(gamma) or not IsList(V) then 
   Error("usage: LocalParameters( <Graph>, <Int> or <List>, [, <PermGroup>] )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
if Length(V)>1 and not IsConnectedGraph(gamma) then
   Error("<gamma> not a connected graph");
fi;
if IsBound(arg[3]) then
   loc:=LocalInfo(gamma,V,0,[],arg[3]);
else
   loc:=LocalInfo(gamma,V);
fi;
if loc.localDiameter=-1 then
   Error("<gamma> not a connected graph");
fi;
return loc.localParameters;
end;

GlobalParameters := function(gamma)
#
# Determines the global parameters of connected, simple graph  gamma.
# The nonexistence of a global parameter is denoted by  -1.
#
local i,j,k,reps,pars,lp,loc;
if not IsGraph(gamma) then 
   Error("usage: GlobalParameters( <Graph> )");
fi;
if gamma.order=0 then
   return [];
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
reps:=gamma.representatives;
loc:=LocalInfo(gamma,reps[1]);
if loc.localDiameter=-1 then
   Error("<gamma> not a connected graph");
fi;
pars:=loc.localParameters;
for i in [2..Length(reps)] do
   lp:=LocalInfo(gamma,reps[i]).localParameters;
   for j in [1..Maximum(Length(lp),Length(pars))] do
      if not IsBound(lp[j]) or not IsBound(pars[j]) then
	 pars[j]:=[-1,-1,-1];
      else
	 for k in [1..3] do
	    if pars[j][k]<>lp[j][k] then
	       pars[j][k]:=-1;
	    fi;
	 od;
      fi;
   od;
od;
return pars;
end;

IsDistanceRegular := function(gamma)
#
# Returns  true  iff  gamma  is distance-regular 
# (a graph must be simple to be distance-regular).
#
local i,reps,pars,lp,loc,d;
if not IsGraph(gamma) then 
   Error("usage: IsDistanceRegular( <Graph> )");
fi;
if gamma.order=0 then
   return true;
fi;
if not IsSimpleGraph(gamma) then
   return false;
fi;
reps:=gamma.representatives;
loc:=LocalInfo(gamma,reps[1]);
pars:=loc.localParameters;
d:=loc.localDiameter;
if d=-1 then  # gamma not connected
   return false;
fi;
if -1 in Flat(pars) then 
   return false;
fi;
for i in [2..Length(reps)] do
   loc:=LocalInfo(gamma,reps[i]);
   if loc.localDiameter<>d then 
      return false;
   fi;
   if pars <> loc.localParameters then
      return false;
   fi;
od;
return true;
end;

DistanceSet := function(arg)
#
# Let  gamma=arg[1],  distances=arg[2],  V=arg[3].
# Returns the set of vertices  w  of  gamma,  such that  d(V,w)  is in
# distances (a list or singleton distance). 
# If  arg[4]  is bound, then it is assumed to be a subgroup
# of  Aut(gamma)  stabilizing  V  setwise.
#
local gamma,distances,V,maxlayer,distset,laynum,x,i;
gamma:=arg[1];
distances:=arg[2];
V:=arg[3];
if IsInt(distances) then   # assume  distances  consists of a single distance.
   distances:=[distances];
fi;
if not (IsGraph(gamma) and IsList(distances) and (IsList(V) or IsInt(V))) then
   Error(ConcatenationString("usage: DistanceSet( <Graph>, <Int> or <List>, ",
			     "<Int> or <List> [, <PermGroup> ] )"));
fi;
if not IsSet(distances) then 
   distances:=Set(distances);
fi;
distset:=[];
if Length(distances)=0 then
   return distset;
fi;
maxlayer:=Maximum(distances)+1;
if IsBound(arg[4]) then
   laynum:=LocalInfo(gamma,V,maxlayer,[],arg[4]).layerNumbers;
else
   laynum:=LocalInfo(gamma,V,maxlayer).layerNumbers;
fi;
for i in [1..gamma.order] do
   if laynum[i]-1 in distances then
      Add(distset,i);
   fi;
od;
IsSet(distset);
return distset;
end;

DistanceSetInduced := function(arg)
#
# Let  gamma=arg[1],  distances=arg[2],  V=arg[3].
# Returns the graph induced on the set of vertices  w  of  gamma,  
# such that  d(V,w)  is in distances (a list or singleton distance). 
# If  arg[4]  is bound, then it is assumed to be a subgroup
# of  Aut(gamma)  stabilizing  V  setwise.
#
local gamma,distances,V,distset,H;
gamma:=arg[1];
distances:=arg[2];
V:=arg[3];
if IsInt(distances) then   # assume  distances  consists of a single distance.
   distances:=[distances];
fi;
if IsInt(V) then
   V:=[V];
fi;
if IsBound(arg[4]) then
   H:=arg[4];
elif Length(V)=1 then
   H:=ProbablyStabilizer(gamma.group,V[1]);
else
   H:=Group(());
fi;
if not (IsGraph(gamma) and IsList(distances) and IsList(V) and IsPermGroup(H))
  then
   Error(ConcatenationString("usage: DistanceSetInduced( <Graph>, ",
	 "<Int> or <List>, <Int> or <List> [, <PermGroup> ] )"));
fi;
distset:=DistanceSet(gamma,distances,V,H);
return InducedSubgraph(gamma,distset,H);
end;

DistanceGraph := function(gamma,distances)
#
# Returns graph  delta  with the same vertex set, names,  and group as 
# gamma,  and  (x,y)  is an edge of  delta  iff  d(x,y)  (in gamma)
# is in  distances. 
#
local r,delta,d,i;
if IsInt(distances) then
   distances:=[distances];
fi;
if not IsGraph(gamma) or not IsList(distances) then 
   Error("usage: DistanceGraph( <Graph>, <Int> or <List> )");
fi;
delta:=rec(isGraph:=true,order:=gamma.order,group:=gamma.group,
	   schreierVector:=gamma.schreierVector,adjacencies:=[],
	   representatives:=gamma.representatives);
if IsBound(gamma.names) then
   delta.names:=Copy(gamma.names);
fi;
for i in [1..Length(delta.representatives)] do 
   delta.adjacencies[i]:=DistanceSet(gamma,distances,delta.representatives[i]);
od;
if not (0 in distances) and IsBound(gamma.isSimple) and gamma.isSimple then
   delta.isSimple:=true;
fi;
return delta;
end;

ComplementGraph := function(arg)
#
# Returns the complement of the graph  gamma=arg[1]. 
# arg[2] is true iff loops/nonloops are to be complemented (default:false).
#
local gamma,comploops,i,delta,notnecsimple;
gamma:=arg[1];
if IsBound(arg[2]) then
   comploops:=arg[2];
else
   comploops:=false;
fi;
if not IsGraph(gamma) or not IsBool(comploops) then 
   Error("usage: ComplementGraph( <Graph> [, <Bool> ] )");
fi;
notnecsimple:=not IsBound(gamma.isSimple) or not gamma.isSimple;
delta:=rec(isGraph:=true,order:=gamma.order,group:=gamma.group,
	   schreierVector:=gamma.schreierVector,adjacencies:=[],
	   representatives:=gamma.representatives);
if IsBound(gamma.names) then
   delta.names:=Copy(gamma.names);
fi;
if IsBound(gamma.autGroup) then
   delta.autGroup:=gamma.autGroup;
fi;
if IsBound(gamma.isSimple) then
   if gamma.isSimple and not comploops then
      delta.isSimple:=true;
   fi;
fi;
for i in [1..Length(delta.representatives)] do 
   delta.adjacencies[i]:=Difference([1..gamma.order],gamma.adjacencies[i]);
   if not comploops then
      RemoveSet(delta.adjacencies[i],delta.representatives[i]);
      if notnecsimple and (gamma.representatives[i] in gamma.adjacencies[i])
       then
	 AddSet(delta.adjacencies[i],delta.representatives[i]);
      fi;
   fi;
od;
return delta;
end;

PointGraph := function(arg)
#
# Assuming that  gamma=arg[1]  is simple, connected, and bipartite, 
# this function returns the connected component containing  
# v=arg[2]  of the distance-2  graph of  gamma=arg[1]  
# (default:  arg[2]=1,  unless  gamma has zero
# vertices, in which case a zero vertex graph is returned). 
# Thus, if  gamma  is the incidence graph of a (connected) geometry, and 
# v  represents a point, then the point graph of the geometry is returned.
#
local gamma,delta,bicomps,comp,v,gens,hgens,i,g,j,outer;
gamma:=arg[1];
if IsBound(arg[2]) then 
   v:=arg[2];
else
   v:=1;
fi;
if not IsGraph(gamma) or not IsInt(v) then
   Error("usage: PointGraph( <Graph> [, <Int> ])");
fi;
if gamma.order=0 then
   return Copy(gamma);
fi;
bicomps:=Bicomponents(gamma);
if Length(bicomps)=0 or not IsSimpleGraph(gamma) 
		     or not IsConnectedGraph(gamma) then
   Error("<gamma> not  simple,connected,bipartite");
fi;
if v in bicomps[1] then 
   comp:=bicomps[1];
else
   comp:=bicomps[2];
fi;
delta:=DistanceGraph(gamma,2);
# construct Schreier generators for the subgroup of  gamma.group 
# fixing  comp.
gens:=gamma.group.generators;
hgens:=Set([]);
for i in [1..Length(gens)] do
   g:=gens[i];
   if v^g in comp then
      AddSet(hgens,g);
      if IsBound(outer) then
	 AddSet(hgens,outer*g/outer);
      fi;
   else    # g is an "outer" element
      if IsBound(outer) then
	 AddSet(hgens,g/outer);
	 AddSet(hgens,outer*g);
      else 
	 outer:=g;
	 for j in [1..i-1] do
	    AddSet(hgens,outer*gens[j]/outer);
	 od;
	 g:=g^2;
	 if g <> () then
	    AddSet(hgens,g);
	 fi;
      fi;
   fi;
od;
return InducedSubgraph(delta,comp,Group(hgens,()));
end;

EdgeGraph := function(gamma)
#
# Returns the edge graph, also called the line graph, of the 
# (assumed) simple graph  gamma.
# This edge graph  delta  has the unordered edges of  gamma  as 
# vertices, and  e  is joined to  f  in delta precisely when 
# e<>f,  and  e,f  have a common vertex in  gamma.
#
local delta,i,j,k,edgeset,adj,r,e,f;
if not IsGraph(gamma) then 
   Error("usage: EdgeGraph( <Graph> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
edgeset:=UndirectedEdges(gamma);
delta:=NullGraph(Operation(gamma.group,edgeset,OnSets),Length(edgeset));
if not IsBound(gamma.names) then
   delta.names:=edgeset;
else
   delta.names:=[];
   for i in [1..Length(edgeset)] do
      delta.names[i]:=
	 Set([gamma.names[edgeset[i][1]],gamma.names[edgeset[i][2]]]);
   od;
fi;
delta.isSimple:=true;
for i in [1..Length(delta.representatives)] do
   r:=delta.representatives[i];
   e:=edgeset[r];
   adj:=delta.adjacencies[i];
   for k in [1,2] do
      for j in Adjacency(gamma,e[k]) do
	 f:=Set([e[k],j]);
	 if e<>f then
	    AddSet(adj,Position(edgeset,f));
	 fi;
      od;
   od;
od;
return delta;
end;

QuotientGraph := function(gamma,R)
#
# Returns the quotient graph  delta  of  gamma  defined by the 
# smallest  gamma.group  invariant equivalence relation  S 
# (on the vertices of  gamma)  containing the relation  R  (given
# as a list of ordered pairs of vertices of  gamma).
# The vertices of this quotient  delta  are the equivalence 
# classes of  S,  and  [X,Y]  is an edge of  delta  iff  
# [x,y]  is an edge of  gamma  for some  x in X,  y in Y.
#
local root,Q,F,V,W,i,j,r,q,x,y,names,gens,delta,g,h,m,pos;

root := function(x)
#
# Returns the root of the tree containing  x  in the forest represented 
# by  F,  and compresses the path travelled in this tree to find the root.
# F[x]=-x  if  x  is a root, else  F[x]  is the parent of  x.
#
local t,u;
t:=F[x];
while t>0 do 
   t:=F[t];
od;
# compress path
u:=F[x];
while u<>t do
   F[x]:=-t;
   x:=u;
   u:=F[x];
od;  
return -t;
end;

if not IsGraph(gamma) or not IsList(R) then
   Error("usage: QuotientGraph( <Graph>, <List> )");
fi;
if gamma.order<=1 or Length(R)=0 then
   delta:=Copy(gamma);
   delta.names:=List([1..gamma.order],i->[VertexName(gamma,i)]);
   return delta;
fi;
if IsInt(R[1]) then   # assume  R  consists of a single pair.
   R:=[R];
fi;
F:=[];
for i in [1..gamma.order] do
   F[i]:=-i;
od;
Q:=[];
for r in R do
   x:=root(r[1]);
   y:=root(r[2]);
   if x<>y then
      if x>y then
	 Add(Q,x);
	 F[x]:=y;
      else
	 Add(Q,y);
	 F[y]:=x;
      fi;
   fi;
od;
for q in Q do 
   for g in gamma.group.generators do
      x:=root(F[q]^g);
      y:=root(q^g);
      if x<>y then
	 if x>y then
	    Add(Q,x);
	    F[x]:=y;
	 else
	    Add(Q,y);
	    F[y]:=x;
	 fi;
      fi;
   od;
od;
for i in Reversed([1..gamma.order]) do
   if F[i] < 0 then
      F[i]:=-F[i];
   else
      F[i]:=root(F[i]);
   fi;
od;
V:=Set(F);
W:=[];
names:=[];
m:=Length(V);
for i in [1..m] do
   W[V[i]]:=i;
   names[i]:=Set([]);
od;
for i in [1..gamma.order] do
   AddSet(names[W[F[i]]],VertexName(gamma,i));
od; 
gens:=[];
for g in gamma.group.generators do
   h:=[];
   for i in [1..m] do
      h[i]:=W[F[V[i]^g]];
   od;
   Add(gens,PermList(h));
od;
delta:=NullGraph(Group(gens,()),m);
delta.names:=names;
IsSet(delta.representatives);
for i in [1..gamma.order] do
   pos:=Position(delta.representatives,W[F[i]]);
   if pos<>false then
      for j in Adjacency(gamma,i) do
	 AddSet(delta.adjacencies[pos],W[F[j]]);
      od;
   fi;
od;
if IsLoopy(delta) then
   delta.isSimple:=false;
elif IsBound(gamma.isSimple) and gamma.isSimple then
   delta.isSimple:=true;
else
   Unbind(delta.isSimple);
fi;
return delta;
end;
 
BipartiteDouble := function(gamma)
#
# Returns the bipartite double of  gamma,  as defined in BCN.
#
local gens,g,delta,n,i,adj;
if not IsGraph(gamma) then 
   Error("usage: BipartiteDouble( <Graph> )");
fi;
if gamma.order=0 then 
   return Copy(gamma);
fi;
n:=gamma.order;
gens:=IntransitiveGroupGenerators
	 (gamma.group.generators,gamma.group.generators,n,n);
g:=[];
for i in [1..n] do
   g[i]:=i+n;
   g[i+n]:=i;
od;
Add(gens,PermList(g));
delta:=NullGraph(Group(gens,()),2*n);
for i in [1..Length(delta.adjacencies)] do
   adj:=Adjacency(gamma,delta.representatives[i]);
   if Length(adj)=0 then
      delta.adjacencies[i]:=[];
   else
      delta.adjacencies[i]:=adj+n;
   fi;
od;
if not IsBound(gamma.isSimple) or not gamma.isSimple then
   Unbind(delta.isSimple);
fi;
delta.names:=[];
for i in [1..n] do
   delta.names[i]:=[VertexName(gamma,i),"+"];
od;
for i in [n+1..2*n] do
   delta.names[i]:=[VertexName(gamma,i-n),"-"];
od;
return delta;
end;
   
UnderlyingGraph := function(gamma)
#
# Returns the underlying graph  delta  of  gamma.
# This graph has the same vertex set as  gamma,  and 
# has an edge  [x,y]  precisely when  gamma  has an 
# edge  [x,y]  or  [y,x].
# This function also sets the  isSimple  fields of 
# gamma  (via IsSimpleGraph)  and  delta.
#
local delta,adj,i,x,orb,H;
if not IsGraph(gamma) then 
   Error("usage: UnderlyingGraph( <Graph> )");
fi;
delta:=Copy(gamma);
if IsSimpleGraph(gamma) then
   delta.isSimple:=true;
   return delta;
fi;
for i in [1..Length(delta.adjacencies)] do
   adj:=delta.adjacencies[i];
   x:=delta.representatives[i];
   H:=ProbablyStabilizer(delta.group,x);
   for orb in Orbits(H,adj) do
      if not IsVertexPairEdge(delta,orb[1],x) then
	 AddEdgeOrbit(delta,[orb[1],x]);
      fi;
   od;
od;
delta.isSimple := not IsLoopy(delta); 
return delta;
end;

NewGroupGraph := function(G,gamma)
#
# Returns a copy of  delta  of  gamma, except that  delta.group=G.
#
local delta,i;
if not IsPermGroup(G) or not IsGraph(gamma) then 
   Error("usage: NewGroupGraph( <PermGroup>, <Graph> )");
fi;
delta:=NullGraph(G,gamma.order);
if IsBound(gamma.isSimple) then
   delta.isSimple:=gamma.isSimple;
else
   Unbind(delta.isSimple);
fi;
if IsBound(gamma.autGroup) then
   delta.autGroup:=gamma.autGroup;
fi;
if IsBound(gamma.canonicalLabelling) then
   delta.canonicalLabelling:=gamma.canonicalLabelling;
fi;
if IsBound(gamma.names) then
   delta.names:=Copy(gamma.names);
fi;
for i in [1..Length(delta.representatives)] do
   delta.adjacencies[i]:=Adjacency(gamma,delta.representatives[i]);
od;
return delta;
end;

GeodesicsGraph := function(gamma,x,y)
#
# Returns the graph induced on the set of geodesics between 
# vertices  x  and  y,  but not including  x  or  y. 
# *** This function is only for simple  gamma.
#
local i,n,locx,geoset,H,laynumx,laynumy,w,g,h,rwx,rwy,gens,sch,orb,pt,im;
if not IsGraph(gamma) or not IsInt(x) or not IsInt(y) then 
   Error("usage: GeodesicsGraph( <Graph>, <Int>, <Int> )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
locx:=LocalInfo(gamma,x,0,y);
if locx.distance=-1 then
   Error("<x> not joined to <y>");
fi;
laynumx:=locx.layerNumbers;
laynumy:=LocalInfo(gamma,y,0,x).layerNumbers;
geoset:=[];
n:=locx.distance;
for i in [1..gamma.order] do 
   if laynumx[i]>1 and laynumy[i]>1 and laynumx[i]+laynumy[i]=n+2 then
      Add(geoset,i);
   fi;
od;
H:=ProbablyStabilizer(gamma.group,y);
gens:=gamma.group.generators;
rwx:=RepWord(gens,gamma.schreierVector,x);
rwy:=RepWord(gens,gamma.schreierVector,y);
g:=();
if rwx.representative=rwy.representative then
   for w in Reversed(rwx.word) do 
      g:=g*gens[w]^-1;
   od;
   for w in rwy.word do 
      g:=g*gens[w];
   od;
   pt:=y^g;
   sch:=[];
   orb:=[pt];
   sch[pt]:=();
   for i in orb do
      for h in H.generators do
	 im:=i^h;
	 if not IsBound(sch[im]) then
	    sch[im]:=h;
	    Add(orb,im);
	 fi;
      od;
   od; 
   if IsBound(sch[x]) then
      i:=x;
      h:=();
      while i<>pt do
	 h:=h/sch[i];
	 i:=x^h;
      od;
      g:=g*h^-1;
   fi;
fi;
H:=ProbablyStabilizer(H,x);
if x^g=y and y^g=x then
   Add(H.generators,g);
   H:=Group(H.generators,());
fi;
return InducedSubgraph(gamma,geoset,H);
end;

IndependentSet := function(arg)
#
# Returns a (hopefully large) independent set (coclique) of  gamma=arg[1].
# The returned independent set will contain the (assumed) independent set  
# arg[2]  (default [])  and not contain any element of arg[3]
# (default [], in which case the returned independent set is maximal).
# An error is signalled if arg[2] and arg[3] have non-trivial intersection.
# A "greedy" algorithm is used, and the graph  gamma  must be simple.
#
local gamma,is,forbidden,i,j,k,poss,adj,mindeg,minvert,degs;
gamma:=arg[1];
if not IsBound(arg[2]) then
   is:=Set([]);
else
   is:=Set(Copy(arg[2]));
fi;
if not IsBound(arg[3]) then
   forbidden:=Set([]);
else
   forbidden:=Set(Copy(arg[3]));
fi;
if not IsGraph(gamma) or not IsSet(is) or not IsSet(forbidden) then 
   Error("usage: IndependentSet( <Graph> [, <List> [, <List> ]] )");
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
if Length(Intersection(is,forbidden))>0 then
   Error("<is> and <forbidden> have non-trivial intersection");
fi;
if gamma.order=0 then 
   return Set([]);
fi;
poss:=Difference([1..gamma.order],forbidden); 
SubtractSet(poss,is);
for i in is do
   SubtractSet(poss,Adjacency(gamma,i));
od;
#  is  contains the independent set so far.
#  poss  contains the possible new elements of the independent set.
degs:=[];
for i in poss do
   degs[i]:=Length(Intersection(Adjacency(gamma,i),poss));
od;
while poss <> [] do
   minvert:=poss[1];
   mindeg:=degs[poss[1]];
   for i in [2..Length(poss)] do
      k:=degs[poss[i]];
      if k < mindeg then
	 mindeg:=k;
	 minvert:=poss[i];
      fi;
   od;
   AddSet(is,minvert);
   RemoveSet(poss,minvert);
   adj:=Intersection(Adjacency(gamma,minvert),poss);
   SubtractSet(poss,adj); 
   for i in adj do
      for j in Intersection(poss,Adjacency(gamma,i)) do
	 degs[j]:=degs[j]-1;
      od;
   od;
od;
return is;
end;

CollapsedIndependentOrbitsGraph := function(arg)
#
# Given a subgroup  G=arg[1]  of the automorphism group of
# the graph  gamma=arg[2],  this function returns a graph  delta 
# defined as follows.  The vertices of  delta  are 
# those G-orbits of V(gamma) that are independent sets,
# and  x  is *not* joined to  y  in  delta  iff  x union y  is an
# independent set in  gamma.
# If  arg[3]  is bound then it is assumed to be a subgroup of 
# Aut(gamma) preserving the set of orbits of  G  on the vertices
# of  gamma  (for example, the normalizer of  G  in gamma.group). 
#
local G,gamma,N,orb,orbs,i,j,L,X,rel;
G:=arg[1];
gamma:=arg[2];
if IsBound(arg[3]) then
   N:=arg[3];
else
   N:=G;
fi;
if not IsPermGroup(G) or not IsGraph(gamma) or not IsPermGroup(N) then
   Error(ConcatenationString("usage: CollapsedIndependentOrbitsGraph( ",
	 "<PermGroup>, <Graph> [, <PermGroup>] )"));
fi;
orbs:=Orbits(G,[1..gamma.order]);
i:=1;
L:=[];
rel:=[];
for orb in orbs do
  if Length(Intersection(orb,Adjacency(gamma,orb[1])))=0 then
      Add(L,orb);
      for j in [i+1..i+Length(orb)-1] do
	 Add(rel,[i,j]);
      od;
      i:=i+Length(orb);
   fi;
od;
X:=InducedSubgraph(gamma,Flat(L),N);
return QuotientGraph(X,rel);
end;

CollapsedCompleteOrbitsGraph := function(arg)
#
# Given a subgroup  G=arg[1]  of the automorphism group of
# the graph  gamma=arg[2] (assumed simple),  this function returns
# a graph  delta defined as follows.  The vertices of  delta  are 
# those G-orbits of V(gamma) that are complete subgraphs,
# and  x  is joined to  y  in  delta  iff  x union y  is a
# complete subgraph of  gamma.
# If  arg[3]  is bound then it is assumed to be a subgroup of 
# Aut(gamma) preserving the set of orbits of  G  on the vertices
# of  gamma  (for example, the normalizer of  G  in gamma.group). 
#
local G,gamma,N;
G:=arg[1];
gamma:=arg[2];
if IsBound(arg[3]) then
   N:=arg[3];
else
   N:=G;
fi;
if not IsPermGroup(G) or not IsGraph(gamma) or not IsPermGroup(N) then
   Error(ConcatenationString("usage: CollapsedCompleteOrbitsGraph( ",
	 "<PermGroup>, <Graph> [, <PermGroup>] )"));
fi;
if not IsSimpleGraph(gamma) then
   Error("<gamma> not a simple graph");
fi;
return 
   ComplementGraph(CollapsedIndependentOrbitsGraph(G,ComplementGraph(gamma),N));
end;

CompleteSubgraphs := function(arg)
#
# Let  gamma=arg[1]  (must be simple),  k=arg[2].
#
# This function returns a set  K  of complete subgraphs of  gamma. 
# A complete subgraph is represented by its vertex set.
# If  k>=0  then the elements of  K  have size  k,  else these 
# elements represent maximal complete subgraphs of  gamma.
# The default for  k=arg[2]  is -1  (i.e. maximal complete subgraphs).
# 
# arg[3] is used to control how many complete subgraphs are returned.
# If  arg[3]  is true (the default), then  K  will contain (perhaps 
# properly) a set of  gamma.group  orbit-representatives
# of the size  k  (if k>=0)  or maximal  (if k<0)  complete 
# subgraphs of  gamma.
# If  arg[3]  is false then  K  will contain at most one element. In this
# case, if  k<0  then  K  will contain just one maximal complete subgraph,
# and if  k>=0  then  K  will contain a complete subgraph of size  k 
# iff such a subgraph is contained in  gamma.
#
local gamma,k,allsubs,CompleteSubgraphsSearch;

CompleteSubgraphsSearch := function(gamma,k,forbidden)
#
# Function called by  CompleteSubgraphs  to do all the work.
# The variable  allsubs  is global.
# forbidden  is a set of vertex-names which are not 
# allowed to be in any returned complete subgraph.
# This function assumes that on the initial call, gammma.names is 
# bound, and is equal to [1..gamma.order].
#
local n,i,j,a,delta,adj,rep,ans,ans1,names;
n:=gamma.order;
if k>n then 
   return [];
fi;
if k=0 or n=0 then
   return [[]];
fi;
ans:=[];
names:=gamma.names;
for i in [1..Length(gamma.representatives)] do
   adj:=gamma.adjacencies[i];
   rep:=gamma.representatives[i];
   if not (names[rep] in forbidden) and Length(adj)>=k-1 then
      delta:=InducedSubgraph(gamma,adj,ProbablyStabilizer(gamma.group,rep));
      ans1:=CompleteSubgraphsSearch(delta,k-1,
				    Intersection(forbidden,delta.names));
      for a in ans1 do
	 AddSet(a,names[rep]);
	 AddSet(ans,a);
      od;
      if not allsubs and Length(ans)>0 then
	 return ans;
      fi;
   fi;
   if i<Length(gamma.representatives) then 
      for j in Orbit(gamma.group,rep) do
	 AddSet(forbidden,names[j]);
      od;
   fi;
od;
return ans;
end;

#
# begin  CompleteSubgraphs
#
gamma:=ShallowCopy(arg[1]);
gamma.names:=[1..gamma.order];
if IsBound(arg[2]) then
   k:=arg[2];
else
   k:=-1;
fi;
if IsBound(arg[3]) then
   allsubs:=arg[3];
else
   allsubs:=true;
fi;
if not IsGraph(gamma) or not IsInt(k) or not IsBool(allsubs) then
   Error("usage: CompleteSubgraphs( <Graph>, [<Int> [,<Bool> ]] )");
fi;
if not IsSimpleGraph(gamma) then 
   Error("<gamma> not a simple graph");
fi;
return CompleteSubgraphsSearch(gamma,k,[]);
end;

NextSizeCompleteSubgraphs := function(gamma,K)
#
# Given a set  K  of (vertex sets of) complete subgraphs of the
# (assumed simple) graph  gamma,  this function returns a set  L
# of (vertex sets of) complete subgraphs of  gamma  with the
# following properties:
#   (1) Each element of  L  contains an element of  K  together
#       with an additional vertex of  gamma.
#   (2) For each non-maximal  y  in  K,  there is in  L  a gamma.group
#       representative of each  Size(y)+1  complete subgraph of  gamma
#       containing  y.
#
#  *** Note: This function uses set stabilizers in  gamma.group.
# 
local L,y,int,orbs,orb;
if not IsGraph(gamma) or not IsSet(K) then
   Error("usage: NextSizeCompleteSubgraphs( <Graph>, <Set> )");
fi;
if not IsSimpleGraph(gamma) then 
   Error("<gamma> not a simple graph");
fi;
L:=[];
for y in K do
   if y=[] then
      int:=[1..gamma.order];
   else
      int:=Intersection(List(y,x->Adjacency(gamma,x)));
   fi;
   if int <> [] then
      orbs:=Orbits(Stabilizer(gamma.group,y,OnSets),int);
      for orb in orbs do 
	 AddSet(L,Union(y,[orb[1]]));
      od;
   fi;
od;
return L;
end;

GraphImage := function(gamma,perm)
#
# Returns the image  delta  of the graph  gamma,  under the permutation  perm,
# which must be a permutation of  [1..gamma.order].
#
local delta,i,perminv;
if not IsGraph(gamma) or not IsPerm(perm) then 
   Error("usage: GraphImage( <Graph>, <Perm> )");
fi;
if MaximumMovedPoint([perm]) > gamma.order then
   Error("<perm> must be a permutation of [1..<gamma.order>]");
fi;
delta:=NullGraph(Group(List(gamma.group.generators,x->x^perm),()),gamma.order);
if IsBound(gamma.group.stabChain) then
   StabChain(delta.group,rec(size:=Size(gamma.group)));
fi;
if IsBound(gamma.isSimple) then
   delta.isSimple:=gamma.isSimple;
else
   Unbind(delta.isSimple);
fi;
if IsBound(gamma.autGroup) then
   delta.autGroup:=Group(List(gamma.autGroup.generators,x->x^perm),());
   if IsBound(gamma.autGroup.stabChain) then
      StabChain(delta.autGroup,rec(size:=Size(gamma.autGroup)));
   fi;
fi;
if IsBound(gamma.canonicalLabelling) then
   delta.canonicalLabelling:=gamma.canonicalLabelling*perm;
fi;
perminv:=perm^-1;
if IsBound(gamma.names) then
   delta.names:=List([1..delta.order],i->Copy(gamma.names[i^perminv]));
fi;
for i in [1..Length(delta.representatives)] do
   delta.adjacencies[i]:=
      OnSets(Adjacency(gamma,delta.representatives[i]^perminv),perm);
od;
return delta;
end;

CompleteSubgraphsOfGivenSize := function(arg)
#
# Let  gamma=arg[1]  be a simple graph, and  k=arg[2]  a non-negative
# integer.
#
# This function returns a set  K  of complete subgraphs of  gamma
# of given size  k  (if such subgraphs exist, else the empty set 
# is returned). A complete subgraph is represented by its vertex set.
# 
# The (optional) boolean parameter  arg[3]  is used
# to control how many complete subgraphs are returned.
# If  arg[3]  is true (the default), then  K  will contain 
# (perhaps properly) a set of  gamma.group  orbit-representatives
# of the size  k  complete subgraphs of  gamma.
# If  arg[3]  is false then  K  will contain at most one element,
# and will contain one element iff  gamma  contains a complete 
# subgraph of size  k.
#
# If  arg[4]  is bound and has value  true,  then it is assumed
# that all complete subgraphs of size  k  of  gamma  are maximal.
# (This assumption can help the function run faster.)
# Otherwise no such assumption is made.
#
# If  arg[5]  is bound then it contains  true,  false,  or a rational
# number  colnum.  If not bound or false then this parameter has no effect.
# Otherwise the following takes place.  
# If arg[5]=true, then first a reasonable value of colnum is calculated. 
# We then try (properly) vertex-colouring  gamma  when  k>n/colnum,  
# in the hope of finding that gamma can be so-coloured with fewer 
# than  k  colours, and so has no complete subgraph of size  k. 
#
# If the parameter  arg[6]  is bound, then it is assumed to be a 
# gamma.group invariant list (where the action permutes the list indices) 
# of length  gamma.order  of positive integer weights
# corresponding to the vertices, in which case, a complete subgraph 
# of size  k  means a complete subgraph whose vertex-weights sum to  k.
#
local gamma,k,allsubs,allmaxes,colnum,weights,weighted,perm,col,mm,i,cw,
      CompleteSubgraphsOfGivenSizeSearch;

CompleteSubgraphsOfGivenSizeSearch := function(gamma,k)
#
# Function called by  CompleteSubgraphsOfGivenSize  to do all the work.
# The variables  allmaxes, allsubs, colnum, weights, and weighted are 
# global. 
#
local n,i,j,delta,adj,rep,a,ans,ans1,names,G,orb,ll,mm,forbidden,
      A,nactive,wt,bigwt,cw,sortperm,CompleteSubgraphsOfGivenSizeSearch1;

CompleteSubgraphsOfGivenSizeSearch1 := function(active,k)
#
# This function does the work of  CompleteSubgraphsOfGivenSizeSearch  
# when the group associated with the graph is trivial.  We assume the 
# graph is represented by the rows  active  of the  n by n  adjacency 
# matrix  A.
# The variables n, A, names, allmaxes, allsubs, colnum, weights, and 
# weighted are global.
# The set  active  may be changed by this function.
#
local mask,a,b,c,col,flag,i,j,ans,ans1,ll,mm,nactive,wt,bigwt,cw,sortperm,rep;
if weighted then
   nactive:=Sum(List(active,a->weights[names[a]]));
else
   nactive:=Length(active);
fi;
if k > nactive then 
   return [];
fi;
if k = 0 then
   return [[]];
fi;
mask := BlistList([1..n],active);
repeat
   flag := true;
   mm := -1;
   bigwt := weights[names[active[1]]]; 
   for i in active do
      b:=IntersectionBlist(mask,A[i]);
      if weighted then
	 ll:=Sum(List(ListBlist(names,b),a->weights[a]));
      else
	 ll:=SizeBlist(b);
      fi;
      wt := weights[names[i]];
      if ll < k-wt or wt > k then  # remove vertex i
	 flag:=false;
	 nactive:=nactive-wt;
	 if nactive < k then
	    return [];
	 fi;
	 mask[i]:=false;
	 RemoveSet(active,i);
      fi;   
      if flag and ll > mm and wt = bigwt then
	 mm := ll;
	 j := i;
      fi;
   od;
until flag;
if k = nactive then  # we are down to a complete graph of size k
   return [Set(List(active,a->names[a]))];
fi;
if j = active[1] then 
   sortperm := ();
else 
   sortperm := (active[1],j);
fi;
if k > nactive/colnum then  # (properly) vertex-colour the graph
   col:=[];
   mm:=0;  # max. colour used so far
   for i in [1..Length(active)] do  # col[i] := colour of active[i]  
      c:=BlistList([1..mm+1],[]);
      b:=IntersectionBlist(mask,A[active[i]]);
      for j in [1..i-1] do
	 if b[active[j]] then
	    c[col[j]]:=true;
	 fi;
      od;
      j:=1;
      while c[j] do
	 j:=j+1;
      od;
      col[i]:=j;
      if j>mm then
	 mm:=j;
      fi;
   od;
   if weighted then
      cw := List([1..mm],i->0); 
      #  cw[j]  will contain the maximum weight of a vertex having colour j.
      for i in [1..Length(active)] do
	 if weights[names[active[i]]] > cw[col[i]] then
	    cw[col[i]] := weights[names[active[i]]];
	 fi;
      od; 
      mm:=Sum(cw); 
   fi;
   if mm < k then  # there is no complete subgraph of size k.
      return [];
   fi;
fi;
mm:=IntersectionBlist(mask,A[active[1]^sortperm]);
bigwt:=weights[names[active[1]^sortperm]];
ans:=[];
for i in active do
   rep := i^sortperm;
   if not (mm[rep] and (allmaxes or 
		       (not allsubs and weights[names[rep]]=bigwt))) then
      #
      # We must search for complete subgraphs containing rep.
      # Otherwise, we know that all the subgraphs we need to find containing
      # rep either contain v=active[1^sortperm] or contain a vertex not adjacent
      # to v, and we will find these subgraphs without searching here.
      # (This is essentially Theorem 8.2 of Reingold, Nievergelt and Deo,
      # "Combinatorial Algorithms: Theory and Practice", Prentice-Hall,
      # 1977.)
      # 
      b:=IntersectionBlist(mask,A[rep]);
      ans1:=CompleteSubgraphsOfGivenSizeSearch1(ListBlist([1..n],b),
						k-weights[names[rep]]);
      for a in ans1 do
	 AddSet(a,names[rep]);
	 AddSet(ans,a);
      od;
      if Length(ans)>0 and not allsubs then
	 return ans;
      fi;
      mask[rep]:=false;
   fi;
od;
return ans;
end;

#
# begin  CompleteSubgraphsOfGivenSizeSearch
#
G:=gamma.group;
n:=gamma.order;
names:=gamma.names;
if Length(G.generators)=0 or 
   (G.generators[1]=() and MaximumMovedPoint(G.generators)=0) then  
   # Special case:  G is trivial.
   # Make  A,  the adjacency matrix of gamma.
   A:=List([1..n],i->BlistList([1..n],gamma.adjacencies[i]));
   return CompleteSubgraphsOfGivenSizeSearch1([1..n],k);
fi;
if weighted then
   nactive:=Sum(List(names,a->weights[a]));
else
   nactive:=n;
fi;
if k > nactive then 
   return [];
fi;
if k = 0 then
   return [[]];
fi;
forbidden:=[];
mm := -1;
bigwt := weights[names[gamma.representatives[1]]];
for i in [1..Length(gamma.representatives)] do
   rep:=gamma.representatives[i];
   adj:=gamma.adjacencies[i];
   if weighted then
      ll:=Sum(List(adj,a->weights[names[a]]));
   else
      ll:=Length(adj);
   fi;
   wt := weights[names[rep]];
   if ll < k-wt or wt > k then  # delete vertex-orbit containing rep
      UniteSet(forbidden,Orbit(G,rep));
   fi; 
   if forbidden = [] and ll > mm and wt = bigwt then
      j := i;
      mm := ll;
   fi;
od;
if Length(forbidden)=n then
   return [];
fi;
if forbidden<>[] then
   delta:=InducedSubgraph(gamma,Difference([1..n],forbidden),G);
   return CompleteSubgraphsOfGivenSizeSearch(delta,k);
fi;
#  forbidden=[].
if k=nactive then  # we are down to a complete graph of size k
   return [Set(names)];
fi;
if j=1 then
   sortperm := ();
else 
   sortperm := (1,j);
fi;
if k > nactive/colnum then
   col:=VertexColouring(gamma);
   mm:=Length(Set(col));
   if weighted then
      # cw[j]  will contain the maximum weight of a vertex having colour j.
      cw := List([1..mm],i->0);
      for i in [1..gamma.order] do
	 if weights[names[i]] > cw[col[i]] then
	    cw[col[i]] := weights[names[i]];
	 fi;
      od; 
      mm:=Sum(cw);
   fi;
   if mm < k then  # no complete subgraph of size k.
      return [];
   fi;
fi;
mm:=gamma.adjacencies[1^sortperm];
bigwt:=weights[names[gamma.representatives[1^sortperm]]];
ans:=[];
for i in [1..Length(gamma.representatives)] do
   rep:=gamma.representatives[i^sortperm];
   orb:=Set(Orbit(G,rep));
   if not ((allmaxes or (not allsubs and weights[names[rep]]=bigwt))
	    and IsSubset(mm,orb)) then
      #
      # We must search for complete subgraphs containing rep.
      # Otherwise, we know that all the subgraphs we need to find containing
      # rep either contain v=gamma.representatives[1^sortperm] 
      # or contain a vertex not adjacent
      # to v, and we will find these subgraphs without searching here.
      # (This uses a generalization of Theorem 8.2 of Reingold, Nievergelt 
      # and Deo, "Combinatorial Algorithms: Theory and Practice", 
      # Prentice-Hall, 1977.)
      # 
      adj:=gamma.adjacencies[i^sortperm];
      delta:=InducedSubgraph(gamma,Difference(adj,forbidden),
			     ProbablyStabilizer(G,rep));
      ans1:=CompleteSubgraphsOfGivenSizeSearch(delta,k-weights[names[rep]]);
      for a in ans1 do
	 AddSet(a,names[rep]);
	 AddSet(ans,a);
      od;
      if Length(ans)>0 and not allsubs then
	 return ans;
      fi;
      if i<Length(gamma.representatives) then
	 UniteSet(forbidden,orb);
      fi;
   fi;
od;
return ans;
end;

#
# begin  CompleteSubgraphsOfGivenSize
#
gamma:=ShallowCopy(arg[1]);
gamma.names:=[1..gamma.order];
k:=arg[2];
if k<0 then 
   Error("<k> must be non-negative");
fi;
if IsBound(arg[3]) then
   allsubs:=arg[3];
else
   allsubs:=true;
fi;
allmaxes := IsBound(arg[4]) and arg[4]=true;
if not (IsGraph(gamma) and IsInt(k) and IsBool(allsubs) and IsBool(allmaxes)) 
   then Error(ConcatenationString("usage: CompleteSubgraphsOfGivenSize( ",
	"<Graph>, <Int> [,<Bool> [,<Bool> [, ...]]] )"));
fi;
if not IsSimpleGraph(gamma) then 
   Error("<gamma> not a simple graph");
fi;
if k=0 then 
   return [[]];
fi;
if gamma.order=0 then
   return [];
fi;
if IsBound(arg[6]) then
   weights := arg[6];
   if not IsList(weights) or Length(weights)<>gamma.order then
      Error("<weights> not a list of length <gamma.order>");
   fi;
   if ForAny(weights,w->not IsInt(w) or w<=0) then
      Error("all weights must be positive integers");
   fi;
   if ForAny(gamma.group.generators,g->
	     ForAny([1..Length(weights)],i->weights[i^g]<>weights[i])) then
      Error("<weights> not <gamma.group>-invariant");
   fi;
   if ForAll([1..gamma.order-1],i->weights[i]>=weights[i+1]) then
      perm := ();
   else
      perm := Sortex(Copy(weights))
	      *PermList(Reversed(List([1..gamma.order],x->x))); 
	      # note: function List applyed to range to avoid bug (in Sun cc ?).
      gamma := GraphImage(gamma,perm);
   fi;
else
   weights := List([1..gamma.order],i->1);
   perm := ();
fi;
# Now the vertices of  gamma  are ordered in non-increasing order of weight.
weighted := weights[gamma.names[1]] > 1;
if IsBound(arg[5]) and arg[5]<>false then
   if arg[5]=true then  # calculate a reasonable value for  colnum
      col:=VertexColouring(gamma);
      mm:=Length(Set(col));
      if weighted then
	 # cw[j]  will contain the maximum weight of a vertex having colour j.
	 cw := List([1..mm],i->0);
	 for i in [1..gamma.order] do
	    if weights[gamma.names[i]] > cw[col[i]] then
	       cw[col[i]] := weights[gamma.names[i]];
	    fi;
	 od; 
	 mm:=Sum(cw);
      fi;
      colnum:=gamma.order/mm;
      if mm < k then  # no complete subgraphs of size k.
	 return [];
      fi;
   else
      colnum:=arg[5];
   fi;
else
   colnum:=1;
fi;
return CompleteSubgraphsOfGivenSizeSearch(gamma,k);
end;

CayleyGraph := function(arg)
#
# Given a group  G=arg[1]  and a list  gens=arg[2]  of 
# generators for  G,  this function constructs a Cayley graph 
# for  G  w.r.t.  the generators  gens.  The generating list  
# arg[2]  is optional, and if omitted, then we take  gens:=G.generators.  
# The boolean argument  arg[3]  is also optional, and if true (the default)
# then the returned graph is undirected (as if  gens  was closed 
# under inversion whether or not it is). 
#
# The Cayley graph  caygraph  which is returned is defined as follows:
# the vertices (actually the vertex names) of  caygraph  are the elements
# of  G;  if  arg[3]=true  (the default) then vertices  x,y  are 
# joined by an edge iff there is a  g  in  gens with  y=g*x  
# or  y=g^-1*x;  if  arg[3]=false  then vertices  x,y  are 
# joined by an edge iff there is a  g  in  gens with  y=g*x.  
# 
# *Note* It is not checked whether  G = <gens>.  However, even if  G  
# is not generated by  gens,  the function still works as described 
# above (as long as  gens  is contained in  G), but returns a 
# "Cayley graph" which is not connected.
# 
local G,gens,elms,undirected,caygraph;
G:=arg[1];
if IsBound(arg[2]) then 
   gens:=arg[2];
else
   gens:=G.generators;
fi;
if IsBound(arg[3]) then
   undirected:=arg[3];
else
   undirected:=true;
fi;
if not(IsGroup(G) and IsList(gens) and IsBool(undirected)) then
   Error("usage: CayleyGraph( <Group> [, <List> [, <Bool> ]] )");
fi;
elms:=Elements(G);
if IsPermGroup(G) then
   StabChain(G,rec(size:=Length(elms)));
fi;
if not IsSet(gens) then
   gens:=Set(gens);
fi;
caygraph := Graph(G,elms,OnRight,
		  function(x,y) return y*x^-1 in gens; end,true);
#
# Note that  caygraph.group  comes from the right regular action of
# G  as a group of automorphisms of the Cayley graph constructed.  
#
StabChain(caygraph.group,rec(size:=caygraph.order));
if undirected then
  caygraph:=UnderlyingGraph(caygraph);
fi;
return caygraph;
end;

SwitchedGraph := function(arg)
#
# Returns the switched graph  delta  of graph  gamma=arg[1],  
# w.r.t to vertex list (or singleton vertex)  V=arg[2].
#
# The returned graph  delta  has vertex set the same as  gamma. 
# If vertices  x,y  of  delta  are both in  V  or both not in
# V,  then  [x,y]  is an edge of  delta  iff  [x,y]  is an edge
# of  gamma;  otherwise  [x,y]  is an edge of delta  iff  [x,y]
# is not an edge of  gamma.
# 
# If  arg[3]  is bound then it is assumed to be a subgroup 
# of  Aut(gamma)  stabilizing  V  setwise.
#
local gamma,delta,n,V,W,H,A,i;
gamma:=arg[1];
V:=arg[2];
if IsInt(V) then
   V:=[V];
fi;
if not IsGraph(gamma) or not IsList(V) then 
   Error("usage: SwitchedGraph( <Graph>, <Int> or <List>, [,<PermGroup>] )");
fi;
n:=gamma.order;
V:=Set(V);
if not IsSubset([1..n],V) then 
   Error("<V> must be a subset of [1..<n>]");
fi;
if Length(V) > n/2 then
   V:=Difference([1..n],V);
fi;
if V=[] then 
   return Copy(gamma);
fi;
if IsBound(arg[3]) then 
   H:=arg[3];
elif Length(V)=1 then
   H:=ProbablyStabilizer(gamma.group,V[1]);
else
   H:=Group(());
fi;
if not IsPermGroup(H) then
   Error("usage: SwitchedGraph( <Graph>, <Int> or <List>, [, <PermGroup>] )");
fi;
delta:=NullGraph(H,n);
if IsBound(gamma.isSimple) then
   delta.isSimple:=gamma.isSimple;
else
   Unbind(delta.isSimple);
fi;
if IsBound(gamma.names) then
   delta.names:=Copy(gamma.names);
fi;
W:=Difference([1..n],V);
for i in [1..Length(delta.representatives)] do
   A:=Adjacency(gamma,delta.representatives[i]);
   if delta.representatives[i] in V then
      delta.adjacencies[i]:=Union(Intersection(A,V),Difference(W,A));
   else
      delta.adjacencies[i]:=Union(Intersection(A,W),Difference(V,A));
   fi;
od;
return delta;
end;

VertexTransitiveDRGs := function(gpin)
#
# If the input to this function is a permutation group  G=gpin, 
# then it must be transitive, and we first set 
# coladjmats:=OrbitalGraphColadjMats(G).  
# 
# Otherwise, we take  coladjmats:=gpin,  which must be a list of collapsed 
# adjacency matrices for the orbital digraphs of a transitive permutation
# group  G  (on a set V say), collapsed w.r.t. a point stabilizer
# (such as the list of matrices produced by  OrbitalGraphColadjMats  
# or  EnumColadj).  
#
# In either case, this function returns a record (called  result), 
# which gives information on  G. 
# The most important component of this record is the list  
# orbitalCombinations,  whose elements give the combinations of 
# the (indices of) the G-orbitals whose union gives the edge-set 
# of a distance-regular graph with vertex-set  V. 
# The component  intersectionArrays  gives the corresponding 
# intersection arrays. The component  degree  is the degree of
# the permutation group  G,  rank  is its (permutation) rank, and  
# isPrimitive  is true if  G  is primitive, and false otherwise.
# It is assumed that the orbital/suborbit indexing used is the same 
# as that for the rows (and columns) of each of the matrices and 
# also for the indexing of the matrices themselves, with the trivial 
# suborbit first, so that, in particular,  coladjmat[1]  must be an 
# identity matrix.  
#
# The techniques used in this function are described in:
# Praeger and Soicher, "Low Rank Representations and Graphs for 
# Sporadic Groups", CUP, Cambridge, 1997. 
#
# *Warning* This function checks all subsets of [2..result.rank], so 
# the rank of  G  must not be large!
#
local coladjmats,i,rank,comb,loc,sum,degree,prim,result;
if not IsList(gpin) and not IsPermGroup(gpin) then 
   Error("usage: VertexTransitiveDRGs( <List> or <PermGroup> )"); 
fi;
if IsPermGroup(gpin) then 
   # remark: OrbitalGraphColadjMats will check if  gpin  transitive,
   #         so we do not do this here.
   coladjmats := OrbitalGraphColadjMats(gpin);
else 
   coladjmats := gpin;
fi; 
if coladjmats=[] or not IsMat(coladjmats[1]) 
      or not IsInt(coladjmats[1][1][1]) then 
   Error("<coladjmats> must be a non-empty list of integer matrices");
fi;
prim:=true;
rank:=Length(coladjmats);
for i in [2..rank] do 
   if LocalInfoMat(coladjmats[i],1).localDiameter=(-1) then
      prim:=false;
   fi;
od;
degree:=Sum(Sum(List(coladjmats,a->a[1])));
result:=rec(degree:=degree, rank:=rank, isPrimitive:=prim, 
	    orbitalCombinations:=[], intersectionArrays:=[]);
for comb in Combinations([2..rank]) do
   if comb<>[] then
      sum:=Sum(Sublist(coladjmats,comb)); 
      loc:=LocalInfoMat(sum,1);
      if loc.localDiameter <> -1 and loc.localParameters[2][1]=1  
	    and not (-1 in Flat(loc.localParameters)) then
	 Add(result.orbitalCombinations,comb);
	 Add(result.intersectionArrays,loc.localParameters);
      fi;   
   fi;   
od;
return result;
end;

#######################################################################
#
# Next comes the part of  GRAPE  depending on B.D.McKay's  nauty 
# system.
#
# define some global variables so as not to get warning messages
GRAPE_dr_sgens:=0; 
GRAPE_dr_base:=0; 
GRAPE_dr_canon:=0;

AutGroupGraph := function(arg) 
#
# Returns automorphism group of (directed) graph arg[1], using B.McKay's 
# dreadnaut, nauty  programs.
# If arg[2] exists then it is a vertex-colouring (not necessarily proper)
# for the graph, and the subgroup of Aut(graph) preserving this colouring 
# is returned instead of the full automorphism group.
# (Here a vertex-colouring is a list of sets, forming a partion of the 
# vertices. The set for the last colour may be omitted, and the "sets"
# may in fact be lists.)
#
local gamma,gp,col,ftmp,fdre,fg,i; 
gamma:=arg[1];
if not IsGraph(gamma) or (IsBound(arg[2]) and not IsList(arg[2])) then
   Error("usage: AutGroupGraph( <Graph> [, <List> ] )");
fi;
if gamma.order=0 then 
   return Group(());
fi;
if IsBound(arg[2]) then 
   col:=List(arg[2],Set); 
else 
   col:=[]; 
fi;
if IsBound(gamma.autGroup) and (col=[] or col=[Vertices(gamma)]) then
   return gamma.autGroup;
fi;
ftmp:=TmpName(); Exec(Concatenation("touch ",ftmp));
fdre:=TmpName(); Exec(Concatenation("touch ",fdre));
fg:=TmpName(); Exec(Concatenation("touch ",fg));
PrintTo(ftmp,gamma.order,"[\n"); 
for i in [1..gamma.order-1] do 
   AppendTo(ftmp,Adjacency(gamma,i),",\n");
od;
AppendTo(ftmp,Adjacency(gamma,gamma.order),"]\n",col,"\n");
PrintTo( fdre, "d\n" );
ExecPkg( "grape", "bin/gap3todr", Concatenation("<",ftmp,">>",fdre), "." );
AppendTo( fdre, "> ", ftmp, " xq\n" );
ExecPkg( "grape", "bin/dreadnaut", Concatenation("<",fdre), "." );
ExecPkg( "grape", "bin/drtogap3", Concatenation("<",ftmp,">",fg), "." );
Read(fg);
gp:=Group(GRAPE_dr_sgens,());
MakeStabChainStrongGenerators(gp,GRAPE_dr_base,GRAPE_dr_sgens);
if col=[] or col=[Vertices(gamma)] then
   gamma.autGroup:=gp;
fi;
Exec(ConcatenationString("rm -f ",ftmp," ",fdre," ",fg));
return gp;
end;

SetAutGroupCanonicalLabelling := function(gamma) 
#
# Sets the  autGroup  and  canonicalLabelling  fields of  gamma,
# if these fields are not bound.
#
local ftmp1,ftmp2,fdre,fg,adj,i;
if not IsGraph(gamma) then
   Error("usage: SetAutGroupCanonicalLabelling( <Graph> )");
fi;
if IsBound(gamma.autGroup) and IsBound(gamma.canonicalLabelling) then
   return(); 
fi;
if gamma.order=0 then 
   gamma.autGroup:=Group(());
   gamma.canonicalLabelling:=();
   return();
fi;
ftmp1:=TmpName(); Exec(Concatenation("touch ",ftmp1));
ftmp2:=TmpName(); Exec(Concatenation("touch ",ftmp2));
fdre:=TmpName(); Exec(Concatenation("touch ",fdre));
fg:=TmpName(); Exec(Concatenation("touch ",fg));
PrintTo(ftmp1,gamma.order,"[\n"); 
for i in [1..gamma.order-1] do 
   AppendTo(ftmp1,Adjacency(gamma,i),",\n");
od;
AppendTo(ftmp1,Adjacency(gamma,gamma.order),"]\n",[],"\n");
PrintTo(ftmp2,gamma.order,"\n"); 
PrintTo( fdre, "d\n" );
ExecPkg( "grape", "bin/gap3todr", Concatenation("<",ftmp1,">>",fdre), "." );
AppendTo( fdre, "> ", ftmp1, " cx\n>> ", ftmp2, " bq\n" );
ExecPkg( "grape", "bin/dreadnaut", Concatenation("<",fdre), "." );
ExecPkg( "grape", "bin/drtogap3", Concatenation("<",ftmp1," >",fg), "." );
ExecPkg( "grape", "bin/drcanon3", Concatenation("<",ftmp2," >>",fg), "." );
Read(fg);
Exec(ConcatenationString("rm -f ",ftmp1," ",ftmp2," ",fdre," ",fg));
if not IsBound(gamma.autGroup) then 
   gamma.autGroup:=Group(GRAPE_dr_sgens,());
   MakeStabChainStrongGenerators(gamma.autGroup,GRAPE_dr_base,GRAPE_dr_sgens);
fi;
if not IsBound(gamma.canonicalLabelling) then 
   gamma.canonicalLabelling:=GRAPE_dr_canon; 
fi;
end;

IsIsomorphicGraph := function(gamma1,gamma2)
#
# Returns  true  iff  gamma1  and  gamma2  are isomorphic.           
#
local g,i,j,adj1,adj2,x,aut1,aut2,reps1;
if not IsGraph(gamma1) or not IsGraph(gamma2) then
   Error("usage: IsIsomorphicGraph( <Graph>, <Graph> )");
fi;
if gamma1.order <> gamma2.order 
  or VertexDegrees(gamma1) <> VertexDegrees(gamma2) then 
   return false; 
fi;
SetAutGroupCanonicalLabelling(gamma1);
SetAutGroupCanonicalLabelling(gamma2);
aut1:=gamma1.autGroup;
aut2:=gamma2.autGroup;
if Size(aut1) <> Size(aut2) then 
   return false; 
fi;
x:=gamma1.canonicalLabelling^-1*gamma2.canonicalLabelling;
for g in aut1.generators do 
   if not g^x in aut2 then
      return false;
   fi;
od; 
reps1:=OrbitNumbers(aut1,gamma1.order).representatives;
for i in reps1 do 
   adj1:=Adjacency(gamma1,i);
   adj2:=Adjacency(gamma2,i^x);
   if Length(adj1)<>Length(adj2) then 
      return false;
   fi;
   for j in adj1 do
      if not (j^x in adj2) then 
	 return false; 
      fi;
   od;
od;
return true;
end;

GraphIsomorphism := function(gamma1,gamma2)
#
# Returns an isomorphism from  gamma1  to  gamma2,  if  gamma1  and
# gamma2  are isomorphic,  else returns  false.
#
if not IsGraph(gamma1) or not IsGraph(gamma2) then
   Error("usage: GraphIsomorphism( <Graph>, <Graph> )");
fi;
if not IsIsomorphicGraph(gamma1,gamma2) then
   return false;
else
   SetAutGroupCanonicalLabelling(gamma1);
   SetAutGroupCanonicalLabelling(gamma2);
   return gamma1.canonicalLabelling^-1*gamma2.canonicalLabelling;
fi;
end;

PartialLinearSpaces := function(arg)
#
# Let  s  and  t  be positive integers.  Then a *partial linear space*  
# (P,L),  with *parameters*  s,t,  consists of a set  P  of *points*, 
# together with a set  L  of (s+1)-subsets of  P  called *lines*, 
# such that every point is in exactly  t+1  lines, and 
# every pair of (distinct) points is contained in at most one line.
# The *point graph* of a partial linear space  S  having point-set
# P  is the graph with vertex-set  P  and having  (p,q)  an edge iff 
# p<>q  and  p,q  lie on a common line of  S. Two partial linear 
# spaces  (P,L)  and  (P',L')  (with parameters  s,t)  are said 
# to be *isomorphic* if there is a bijection  P-->P'  which induces
# a bijection  L-->L'.
#
# This function returns a list of representatives of distinct isomorphism 
# classes of partial linear spaces with (simple) point graph  ptgraph=arg[1],  
# and parameters  s=arg[2],t=arg[3].  The default is that representatives
# for all isomorphism classes are returned.  
# 
# The integer argument  nspaces=arg[4]  is optional, and has 
# default value  -1,  which means that representatives for all
# isomorphism classes are returned.  If  nspaces>=0  then exactly  nspaces
# representatives are returned if there are at least  nspaces  isomorphism
# classes, otherwise representatives for all isomorphism classes are returned.
#
# In the output of this function, a partial linear space  S  is represented by 
# its incidence graph  delta.  The naming (and ordering) of point-vertices in 
# delta  corresponds to the numbering  1,...,ptgraph.order  of the
# vertices of  ptgraph.  The naming of the line-vertices is by (s+1)-subsets
# of point-vertex names.  The group  delta.group  associated with 
# the incidence graph  delta  is the automorphism group of  S
# (acting on points and lines, and preserving both sets).
#
# If  arg[5]  is bound then it controls the printlevel  (default 0).
# Permitted values for  arg[5]  are 0,1,2.  
# 
# If  arg[6]  is bound then it is assumed to be a list (without repeats)
# of the (s+1)-cliques of  ptgraph.  If known, this can help the function
# to run faster. 
# 
local ptgraph,aut,X,printlevel,I,K,s,t,deg,search,cliques,nlines,
      ans,lines,pts,i,j,k,adj,nspaces;
ptgraph:=arg[1];
s:=arg[2];
t:=arg[3];
if IsBound(arg[4]) then
   nspaces:=arg[4];
else
   nspaces:=-1;
fi;
if IsBound(arg[5]) then
   printlevel:=arg[5];
else
   printlevel:=0;
fi;
if not (IsGraph(ptgraph) and IsInt(s) and IsInt(t) and 
	IsInt(nspaces) and IsInt(printlevel)) 
    or (IsBound(arg[6]) and not IsList(arg[6])) then
   Error(ConcatenationString("usage: PartialLinearSpaces(",
	 "<Graph>, <Int>, <Int>, [<Int>, [<Int>, [<List>]]])"));
fi;
if not printlevel in [0,1,2] then
   Error("<printlevel> must be 0, 1, or 2");
fi;
if s<1 or t<1 then
   Error("<s> and <t> must be positive integers");
fi;
if not IsSimpleGraph(ptgraph) then
   Error("<ptgraph> must be a simple graph");
fi;
nlines:=(t+1)*ptgraph.order/(s+1);      # number of lines
if not IsInt(nlines) or nspaces=0 then  # no partial linear spaces
   return [];
fi;
if IsBound(arg[6]) then
   cliques:=arg[6];
   if ForAny(cliques,x->not IsSet(x) or Length(x)<>s+1) then
      Error("<arg[6]> incorrect");
   fi;
   if Length(cliques) < nlines then   
      return [];
   fi;
fi;  
deg:=s*(t+1);
if ForAny(ptgraph.adjacencies,x->Length(x)<>deg) then
   return [];
fi;
aut:=AutGroupGraph(ptgraph);
ptgraph:=NewGroupGraph(aut,ptgraph);
if not IsBound(cliques) then
   if IsCompleteGraph(ptgraph) then
      cliques:=Combinations([1..ptgraph.order],s+1); 
   else
      K:=CompleteSubgraphsOfGivenSize(ptgraph,s+1,true,false,
	       ptgraph.order/Length(Set(VertexColouring(ptgraph))) );
      cliques:=Concatenation(Orbits(ptgraph.group,K,OnSets));
   fi;
   if Length(cliques) < nlines then 
      return [];
   fi;
fi;
if nlines=t*(s+1)+1 then  # line graph is complete graph
   adj:=function(x,y) return x<>y and Length(Intersection(x,y))<>1; end;
else
   adj:=function(x,y) return x<>y and Length(Intersection(x,y))>1; end;
fi;
X:=Graph(aut,cliques,OnSets,adj,true);
X.isSimple:=true;
Unbind(X.names);
MakeStabChainGivenLimit(X.group,Size(aut));
#
# The set  S  of independent sets of size  nlines  in  X  is in 1-to-1 
# correspondence with (the line sets of) the partial linear spaces
# with point graph  ptgraph,  and parameters  s,t.
# 
# Moreover, the orbits of  X.group  (induced by  Aut(ptgraph))  on  S  
# are in  1-to-1  correspondence with the isomorphism classes 
# of the partial linear spaces with point graph  ptgraph,  
# and parameters  s,t.
# 
# We shall classify  S  modulo  X.group,  in order to classify the 
# required partial linear spaces up to isomorphism
# (except that we stop if  nspaces>0  isomorphism classes are required
# and we find that number of them).
#
if Size(X.group) < Size(aut) then
   #
   # non-faithful action of  aut  on (s+1)-cliques, which is impossible if 
   # a subset of these cliques form the line set of a partial linear space
   # with parameters  s,t > 1.
   #
   return [];
fi;
if printlevel > 0 then
   Print("X.order=",X.order," VertexDegrees(X)=",VertexDegrees(X),"\n");
fi;
I:=IndependentSet(ptgraph);
if printlevel > 0 then
   Print("Length(I)=",Length(I),"\n");
fi;
I:=Concatenation(I,Difference([1..ptgraph.order],I));
#
# We shall build the possible partial linear spaces by determining
# the possible line-sets through  I[1],I[2],I[3],...  (in order)
# and backtracking when necessary.  
#
# It appears to be a good strategy to start  I  with the vertices 
# of a maximal independent set of  ptgraph.
#
ans:=[]; 
#
# Representatives of new isomorphism classes of the required partial
# linear spaces are put in  ans  as and when they are found.
#

search := function ( i, sofar, live, H )
#
# This is the function for the backtrack search.
#
# Given in  sofar  the vertices of  X  indexing the (s+1)-cliques 
# forming all the lines through the points  I[1],...,I[i-1],  this function
# determines representatives for the new isomorphism classes 
# (not in  ans)  of the required partial linear spaces  S,
# such that the line-set of  S  contains all the cliques indexed by elements 
# of  sofar,  but no clique not indexed by an element in the union of  
# sofar  and  live.  (The clique indexed by  v  is simply  cliques[v].)
# 
# This function assumes that  sofar  and  live  are disjoint sets, and
# that  H <= X.group  stabilizes each of  sofar  and  live  (setwise).
# It is also assumed that  X.names  is unbound, or  X.names=[1..X.order].
#
# Note: On entry to this function it is assumed that  nspaces=-1
# or  nspaces > 0.  If  nspaces > 0  then it is assumed that (on entry)
# nspaces  isomorphism classes have not yet been found.  
# If  nspaces > 0  then this function terminates if  nspaces  
# isomorphism classes are found.  
# 
local  j, L, K, ind, k, forbid, nlinesreq;
if printlevel > 1 then
   Print("\ni=",i," Size(H)=",Size(H));
fi;
if Length(sofar)=nlines then
   #
   # partial linear space found.
   # check if its isomorphism class is new.
   #
   j := 0;
   repeat
      j := j+1;
      if j > Length(ans) then
	 if printlevel > 1 then
	    Print("\n",List(sofar,x->cliques[x]),"\n");
	 fi;     
	 Add(ans,sofar);
	 return;
      fi;
      if RepresentativeOperation(X.group,sofar,ans[j],OnSets)<>false then
	 return;
      fi;
   until false;
fi;
nlinesreq:=(t+1)-Length(Filtered(sofar,x->I[i] in cliques[x]));
#
#  nlinesreq  is the number of new lines through  I[i]  that 
# must be found.
#
if nlinesreq=0 then
   search(i+1,sofar,live,H);
   return;
fi;
L:=Filtered( live, x->I[i] in cliques[x]);
if printlevel > 1 then    
   Print(" nlinesreq=",nlinesreq,"  Length(L)=",Length(L));
fi;    
if Length(L) < nlinesreq then
   return;
fi;
H := Stabilizer( H, L, OnSets );
ind := ComplementGraph(InducedSubgraph( X, L, H ));
K := CompleteSubgraphsOfGivenSize( ind, nlinesreq, true, true, 
     ind.order / Length(Set(VertexColouring(ind))) );
# 
#  K  contains the sets of possible additional lines through  I[i].
# 
if K = [] then
   return;
fi;
if Length(K) > 1 and Size(ind.group) > 1 then
   #
   # We perform some partial isomorph rejection.
   #
   K := OrbitRepresentatives( ind.group, K, OnSets );
fi;
if printlevel > 1 then    
   Print("  Length(K)=",Length(K));
fi;    
for k in K  do
   L := List( k, x->ind.names[x] );
   forbid := Union( L, Union(List(L,x->Adjacency(X,x))) );
   search( i+1, Union( sofar, L), Difference(live,forbid),
	   Stabilizer(H,L,OnSets) ); 
   if nspaces>=0 and Length(ans)=nspaces then
      return;
   fi;
od;
end;

search(1,[],[1..X.order],X.group);   
for i in [1..Length(ans)] do
   #
   # Determine the incidence graph of the partial linear space
   # whose lines are  List(ans[i],x->cliques[x]),  and store the 
   # result in  ans[i].
   #
   pts:=List([1..ptgraph.order],x->[]);
   for j in ans[i] do
      for k in cliques[j] do
	 AddSet(pts[k],j);
      od;
   od;
   aut:=Operation(Stabilizer(X.group,ans[i],OnSets),pts,OnSets);
   lines:=Set( List(ans[i],x->cliques[x]) );
   ans[i]:=Graph(aut,
		 Concatenation([1..ptgraph.order],lines),
		 function(x,g)
		    if IsInt(x) then
		       return x^g;
		    else
		       return OnSets(x,g);
		    fi;
		 end,
		 function(x,y)
		    if IsInt(x) then
		       return IsSet(y) and x in y;
		    else 
		       return IsInt(y) and y in x;
		    fi;
		 end,
		 true);
   ans[i].isSimple:=true;
od;
return ans;
end;

#######################################################################
#
# Now comes the part of  GRAPE  depending on the  enum3  and  enum3ca  
# standalones. 
#

# define some global variables so as not to get warning messages
GRAPE_coladjmatseq:=0;
GRAPE_tc_permgens:=0; 

Enum := function(arg) 
#
# Function to run the  enum  coset enumerator. 
# arg[1] is the file containing coset enumerator input (default: stdin).
# arg[2] is true iff a list of group generators should be returned 
#        by the function (default: true).
# arg[3] is true iff various information is to be printed  (default: true). 
# arg[4] is the file to keep a record of the enumeration performed 
#        (default: record not kept). 
#
# *Warning* The use of this function causes various files with names 
# begining with  GRAPE_  to be created and deleted in the current 
# directory. If you run another program using Enum or EnumColadj 
# (or enum3 or enum3ca) then run it in a different directory to 
# avoid file clashes. 
#
# The enumerator input in the ASCII file arg[1] should have the format 
# described, below. 
# 
# The input has five main sections which should be separated by full stops
# as follows:
# 
#    all generators.generators not known to be involutions.
#    subgroup generators.coxeter relations.other relations.
# 
# Here are two presentations for Alt(5) in this form (the first
# enumerates over the identity, the second over <a,b>):
#
#   ab.b...b3,(ab)5.  
#   abc..a,b.ab5bc3.(abc)5.
#
# Throughout, blanks and linefeeds are ignored, commas and semicolons are
# equivalent, and upper and lower case letters are distinguished. In
# addition, in sections 1, 2, and 4, commas and semicolons are ignored.
#
# The generators must come from the set {A,B,...,Z,a,b,...,z}.
# Generators to be used are given in the first section.
#
# The generators x for which you do NOT want the assumption that 1=xx
# should be put in the second section.
#
# The third section contains words for the generators of the subgroup
# over which you want to enumerate. These words should be separated by
# commas (or semicolons).
#
# We give the syntax for a word ({ } means zero or more times):
#     word ::= null | term word 
#     term ::= '1' | factor power 
#     factor ::= generator | '(' word ')' | commutator 
#     commutator ::= '[' word ',' word {',' word} ']' 
#     power ::= null | unsigned integer | '-' | '-' unsigned integer 
#
# Multiplication is denoted by juxtaposition, a factor is inverted by
# following it with '-', and a factor is raised to the n-th power when
# followed by the integer n.  Commutators are left normed so that
# [a,b,c,...] means [[a,b],c,...].  A term consisting of just '1' is the
# empty term, and so a word consisting of just '1' is the empty word.
# Additions to the above syntax are that left parenthesis '(' and left
# square bracket '[' are equivalent in words, as are ')' and ']'.
#
# The fourth section contains the Coxeter relations. If you do not want
# any of the implications which follow then leave this section empty. If
# X and Y are generators, and n is an unsigned integer then 'XYn' in this
# section denotes the relator (XY)n.  If  A,B,C,...,V,W,X,Y,Z,...  are
# generators, and  k,m,n,...  are unsigned integers then
# 'ABC ... VWXkYmZn ...' is a shorthand for 
# 'AB3 BC3 ... VW3 WXk WYm WZn ...'.  
# Such a shorthand expression should always end with an (unsigned)
# integer.  Each unordered pair of generators should be specified in this
# section at most once. If a pair X,Y is not specified then the relator
# (XY)2 is assumed. Specifying 'XY0' means that no relator is entered or
# assumed for the X,Y pair.  Note: no assumptions or implications about
# the orders of the generators are made in this section.  Here is an
# example presentation for the Weyl group of E6 (to be enumerated over
# the Weyl group of D5):
#   abcdef..a,b,c,d,f.abcde3 cf3..  
#
# The fifth section is for relations which are not of Coxeter type.
# They should be separated by commas (or semicolons), and each should be
# a word by itself (a relator) or an expression of the form
# word1=word2=...=wordk . This is translated to be 
# (word1)-1 word2,...,(word1)-1 wordk, so it is usually most efficient
# if word1 is the shortest word amongst word1,word2,...,wordk.  Here is
# an example presentation for Alt(5)wr2:
#  wabcd.w.w,a,b,c.abcd3 wa0b0c0.w3,[w,a],[w,b],a=(wc)3.
#
local infile,gensfile,printfile;
if not IsBound(arg[1]) then 
   infile:=""; 
else 
   infile:=arg[1]; 
fi;
if infile<>"" then 
   infile:=ConcatenationString("< ",infile); 
fi;
if not IsBound(arg[2]) then 
   arg[2]:=true; 
fi;
if arg[2] then 
   gensfile:="GRAPE_tc_permgens.g"; 
else 
   gensfile:="''"; 
fi;
if not IsBound(arg[3]) then 
   arg[3]:=true; 
fi;
if arg[3] then 
   printfile:=""; 
else 
   printfile:="> /dev/null"; 
fi;
if not IsBound(arg[4]) then 
   arg[4]:="''"; 
fi;
ExecPkg( "grape", "bin/enum3",
	Concatenation( arg[4], " ", gensfile, " ", infile, " ", printfile ),
	"." );
if arg[2] then 
   Read(gensfile); 
   Exec(ConcatenationString("rm ",gensfile)); 
   return GRAPE_tc_permgens;
fi;
end;

EnumColadj := function(arg) 
#
# Function to run the  enum3ca  programs, and return a sequence
# of intersection matrices for the orbital graphs for
# the (transitive) permutation group produced by the enumerator.
# (Note: before GRAPE 2.1 the intersection matrix for the trivial 
# orbital graph was not included.)
#
# arg[1] is the file containing the coset enumerator input (default: stdin).
# See the documentation for function  Enum  for the format of this input.
#
# *Warning* The use of this function causes various files with names 
# begining with  GRAPE_  to be created and deleted in the current 
# directory. If you run another program using Enum or EnumColadj 
# (or enum3 or enum3ca) then run it in a different directory to 
# avoid file clashes. 
#
local infile,n;
if not IsBound(arg[1]) then 
   infile:=""; 
else 
   infile:=arg[1]; 
fi;
if infile<>"" then 
   infile:=ConcatenationString("< ",infile); 
fi;
ExecPkg( "grape", "bin/enum3ca", Concatenation("> /dev/null ",infile), "." );
Read("GRAPE_coladj.g"); 
Exec("rm GRAPE_coladj.g GRAPE_coladj.tex");
if GRAPE_coladjmatseq=[] then
   n:=1;
else
   n:=Length(GRAPE_coladjmatseq[1]);
fi;
return Concatenation([IdentityMat(n)],GRAPE_coladjmatseq);
end;
