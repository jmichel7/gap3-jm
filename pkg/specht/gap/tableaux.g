#########################################################################
##                                                                     ##
##   SPECHT This file contains functions for generating                ##
##     and manipulating (semistandard) tableaux                        ##
##                                                                     ##
##   Andrew Mathas                                                     ##
##                                                                     ##
#########################################################################

#D Functions for generating and manipulating (semistandard) tableaux

#C New functions for version 3.1:                                        
#C  o ActOnMultiTableau                                                  
#C  o ActOnTableau                                                       
#C  o AreLadderEquivalentTableaux                                        
#C  o BumpingTableau                                                     
#C  o BumpTableau                                                        
#C  o BumpTableaux                                                       
#C  o ConjugateTableau                                                   
#C  o CoxeterWordTableau                                                 
#C  o DescentSetTableau                                                  
#C  o DominatesTableau                                                   
#C  o FirstSemistandardTableauType                                       
#C  o FirstTableau                                                       
#C  o InitialSuperTableau                                                
#C  o InitialTableau                                                     
#C  o IsComposition                                                      
#C  o IsPartition                                                        
#C  o IsRowSemistandardTableau                                           
#C  o IsRowStandardTableau                                               
#C  o IsSemistandardTableau                                              
#C  o IsStandardTableau                                                  
#C  o LadderClassTableau                                                 
#C  o LadderSequencePartition                                            
#C  o LastSemistandard                                                   
#C  o LastSemistandardTableauType                                        
#C  o LastTableau                                                        
#C  o MappingPermTableaux                                                
#C  o ResidueSequence                                                    
#C  o ReverseRobinsonSchenstedCorrespondent                              
#C  o RobinsonSchenstedCorrespondent                                     
#C  o Subtableau                                                         
#C  o SuperSemistandardFromStandard                                      
#C  o SuperSemistandardtableaux                                          
#C  o SymmPermTableau                                                    
#C  o TypedSemistandardTableaux                                          
#C  o TypedSuperSemistandardTableaux                                     
#C  o TypeTableau                                                        
#C  o UnbumpTableaux                                                     
#C  o VeryTightTableau                                                   

## January 2004
##   o hacked SemistandardTableaux so that it can also return the
##     set of (m|n)-semistandard tableaux of type (mu|nu). Added
##     various routines for semistandard tableaux

## December 1998
##   o small conditional bug in SemistandardTableaux() fixed.

## July 1997
##   o fixed a bug (reported by Schmuel Zelikson) in 
##     SemistandardTableaux(); added Type and Shape
##     functions for tableaux and allowed the type of
##     a tableau to be a composition which has parts 
##     which are zero.

## March 1996

#C the dual tableaux of t = [ [t_1, t_2, ...], [t_k, ...], ... ]
ConjugateTableau:=function(t)
  local  d, r, c;
  if t=[[]] then return t; fi;
  d:=List([1..Length(t[1])], r->[]);
  for r in [1..Length(t)] do
    for c in [1..Length(t[r])] do
      d[c][r]:=t[r][c];
    od;
  od;
  return d;
end;   ## ConjugateTableau

#U ConjugateMultiTableau(t)......................................
#M Return the multi-tableau which is conjugate to <t>.           
ConjugateMultiTableau:=function(t)
  return List(Reversed(t), ConjugateTableau);
end;

#C return true if <mu> is a partition
IsPartition:=function(mu)
  return mu=[] or (IsList(mu) and ForAll(mu, IsInt) 
         and ForAll([1..Length(mu)-1], i->mu[i]>=mu[i+1])
	 and mu[Length(mu)]>=0);
end;
IsMultiPartition:=function(mu) local m;
  return IsList(mu) and ForAll(mu,m->IsPartition(m)); 
end;

#C return true if <mu> is a composition
IsComposition:=function(mu)
  return mu=[] or (IsList(mu) and ForAll(mu, IsInt));
end;

# Return <true> if the tableau <t> is row standard
IsRowStandardTableau:=function(t)
  return ForAll(t, IsSet);
end;

# Return <true> if the tableau <t> is row semistandard; that is, the
# entries of t weakly increase along rows and strictly increase down
# columns.
IsRowSemistandardTableau:=function(t) local row, i;
  return ForAll(t, row->ForAll([1..Length(row)-1],i->row[i+1]>=row[i]));
end;

# return true if <t> is super row standard; that is, ...
#IsSuperRowSemistandard:=function(m,t)
#  Error("I don't know what this means yet...\n");
#end;

# Return <true> if the tableau <t> is standard
IsStandardTableau:=function(t)
  return ForAll(t, IsSet) and ForAll(ConjugateTableau(t), IsSet);
end;

# Return <true> if the tableau <t> is a standard multi-tableau
IsMultiStandardTableau:=function(tab) local t;
  return IsList(tab) and ForAll(tab, t->IsStandardTableau(t));
end;


# Return <true> if the tableau <t> is row semistandard; that is, the
# entries of t weakly increase along rows and strictly increase down
# columns.
IsSemistandardTableau:=function(t)
  return ForAll(t, row->ForAll([1..Length(row)-1], i->row[i+1]>=row[i]))
     and ForAll(ConjugateTableau(t),IsSet);
end;

#C Returns the list of semi-standard nu-tableaux of content mu.
##   Usage: SemistandardTableaux(nu, mu) or SemistandardTableaux(nu).
##      or  SemistandardTableaux(m,n,lambda,mu,nu)
## The last usage returns the set of (m|n)-semistandard
## lambda-tableaux of type (mu|nu).
## If mu is omitted the list of all semistandard nu-tableaux is return;
## otherwise only those semistandard nu-tableaux of type mu is returned.
##   Nodes are placed recursively via FillTableau in such a way so as
## so avoid repeats.
SemistandardTableaux:=function(arg) 
  local FillTableau, ss, i, j, nu, mu, N, M, t, nss;
  ## FillTableau adds upto <n>nodes of weight <i> into the tableau <t> on 
  ## or below its <row>-th row (similar to the Littlewood-Richardson rule).
  FillTableau:=function(t, i, n, row)
    local row, max, nn, nodes, nt, r;
    if row>M+Length(mu) then return;
    elif n=0 then            # nodes from i-th row of mu have all been placed
      if i=M+Length(mu) then # t is completely full
	      Add(ss, t); 
        return;
      else
        while i<M+Length(mu) and n=0 do
          i:=i + 1;          # start next row of mu
          n:=mu[i-M];        # mu[i] nodes to go into FillTableau
        od;
        row:=1;              # starting from the first row
      fi;
    fi;
    for r in [row..Minimum(Length(nu),Length(t))] do
      if Length(t[r]) < nu[r] then
        if r = 1 then max:=nu[1]-Length(t[1]);
        elif Length(t[r-1]) > Length(t[r]) and t[r-1][Length(t[r]+1)]<i and Length(t[r-1])>=n then
          max:=Position(t[r-1], i);
          if max=false then max:=Length(t[r-1])-Length(t[r]); #no i in t[r-1]
          else max:=max-1-Length(t[r]);
          fi;
        else max:=0;
        fi;
        max:=Minimum(max, n, nu[r]-Length(t[r]));
        nodes:=[];
        for nn in [1..max] do
          Add(nodes, i);
          nt:=Copy(t);
          Append(nt[r], nodes);
          FillTableau(nt, i, n - nn, r + 1);
        od;
      fi;
    od; 
    r:=Length(t);
    if r < Length(nu)  and n <= nu[r+1] and n <= Length(t[r]) 
    and t[r][n] < i  then
      Add(t, List([1..n], nn->i));
      if i < M+Length(mu) then FillTableau(t, i+1, mu[i+1-M], 1);
      else Add(ss, t); 
      fi;
    fi;
  end;

  M:=0; # only non-zero when contructing super semistandard tableau.
  if Length(arg)=5 and IsInt(arg[1]) and IsInt(arg[2]) 
  and ForAll(arg{[3..5]}, IsList) then ## supersemistandard
    M:=arg[1]; N:=arg[2]; 
    nu:=Copy(arg[3]); mu:=Copy(arg[4]); 
    if Length(mu)>M or Length(arg[5])>N or Sum(nu)<>Sum(mu)+Sum(arg[5]) 
    then Error("SuperSemistandard(<m>,<n>,<tau>,<mu>,<nu>)");
    fi;
    # necessary conditions for semistandard <nu>-tableau of type
    # (mu|tau) to exist - may also be sufficient ???
    Append(mu, [1..Sum(arg[5])]*0+1);
    N:=Copy(arg[5]); Append(N, [1..Sum(arg[4])]*0+1);
    if not ( Dominates(nu,mu) and Dominates(ConjugatePartition(N),nu) ) then
      return [];
    fi;
    mu:=Copy(arg[4]);
    if Length(mu)=0 then ss:=[ [[]] ];
    elif Length(mu)=1 then ss:=[ [[1..mu[1]]*0+1] ];
    else
      ss:=[]; M:=0;
      FillTableau([ [1..mu[1]]*0+1 ], 2, mu[2], 1);
    fi;
    nu:=ConjugatePartition(nu); mu:=arg[5]; M:=arg[1];
    if Length(mu)=0 then return Set( ss );
    else
      if ss=[[[]]] then
	ss:=[];
	nss:=[ ConjugateTableau([[1..mu[1]]*0+M+1]) ]; 
	# we've just placed the first row of mu
	M:=M+1;
	mu:=mu{[2..Length(mu)]};
      else nss:=ss; ss:=[]; 
      fi;
      if mu=[] then return nss;
      else
        for t in nss do
          FillTableau(ConjugateTableau(t), M+1, mu[1], 1);
        od;
      fi;
      return Set(List(ss, ConjugateTableau));
    fi;
  elif Length(arg)=2 and IsList(arg[1]) and IsList(arg[2]) then
    nu:=Copy(arg[1]); mu:=Copy(arg[2]);
    if Sum(nu) <> Sum(mu) then
      Error("<nu> and <mu> must be partitions of the same integer.\n");
    fi;

    ## no semi-standard nu-tableau with content mu
    if not Dominates(nu, mu) then    
      if mu<>nu then return [];
      else return [ List([1..Length(mu)], i->[1..mu[i]]*0+i) ];
      fi;
    elif Length(mu)=1 then # then mu=nu=[n], say
       return [ [[1..mu[1]]*0+1] ];
    fi;

    ss:=[];              ## will hold the semi-standard tableaux
    FillTableau([ [1..mu[1]]*0+1 ], 2, mu[2], 1);

    return Set(ss); 
  elif Length(arg)=1 and IsList(arg[1]) then
    nu:=Flat(arg);
    ss:=[];
    for mu in Partitions(Sum(nu)) do
      if mu = nu then  # add semistandard mu-tableau of type mu.
        Append(ss,[ List([1..Length(mu)], i->[1..mu[i]]*0+i ) ]);
      elif Dominates(nu,mu) then FillTableau([[1..mu[1]]*0+1], 2, mu[2], 1);
      fi;
    od;
    return Set(ss);
  else 
    Error("Usage: SemistandardTableaux(<nu>,<mu>)\n",
          "    or SemistandardTableaux(<nu>)");
  fi;
end;   ## SemistandardTableaux

# return the set of (m|n)-semistandard lambda-tableaux of type (mu|nu)
SuperSemistandardtableaux:=function(m,n,lambda,mu,nu)
  return SemistandardTableaux(m,n,lambda,mu,nu);
end;

# return the set of semistandard tableaux of type mu
TypedSemistandardTableaux:=function(mu) 
  local tabs, lambda;
  tabs:=[];
  for lambda in Partitions(Sum(mu)) do
    if Dominates(lambda,mu) then
      Append(tabs, SemistandardTableaux(lambda,mu));
    fi;
  od;
  return tabs;
end;

# return the set of (m|n)-semistandard tableaux of type (mu|nu)
TypedSuperSemistandardTableaux:=function(m,n,mu,nu) 
  local tabs, lambda;
  tabs:=[];
  for lambda in Partitions(Sum(mu)+Sum(nu)) do
    Append(tabs, SemistandardTableaux(m,n,lambda,mu,nu));
  od;
  return tabs;
end;

#M Add [i..n] to row <row> or lower of t so as to obtain a 
#M standard mu-tableau, where mu is a partition of n.       
SPECHT.FillStandardTableau:=function(t,i,n,row,mu) local tabs, s, r;
  tabs:=[];
  r:=row;
  while r<=Length(mu) do
    if Length(t[r])<mu[r] 
      and (r=1 or Length(t[r])<Length(t[r-1]) ) then
      s:=Copy(t);
      Add(s[r],i);
      if i=n then Add(tabs,s);
      else Append(tabs, SPECHT.FillStandardTableau(s,i+1,n,1,mu));
      fi;
    fi;
    r:=r+1;
  od;
  return tabs;
end;

#C Return a list of the standard tableau of shape  mu)
StandardTableaux:=function(arg) local mu, n, emptyTab;
  mu:=Flat(arg);
  if not IsPartition(mu) then
    Error("usage, StandardTableaux(<partition>");
  fi;
  n:=Sum(Flat(mu));
  emptyTab:=List(mu,r->[]);
  return SPECHT.FillStandardTableau(emptyTab,1,n,1,mu);
end;

# Return <true> if the tableau <t> is standard
IsStandardTableau:=function(t)
  return ForAll(t, IsSet) and ForAll(ConjugateTableau(t), IsSet);
end;

#U StandardTableauxSemistandardTableau(T)..................................
#M Return the set of standard tableaux which correspond to the semistandard
#M tableaux T.                                                             
StandardTableauxSemistandardTableau:=function(T)
  local n, m, lengths, tabs, row, len, newtabs, s, c, rows, t, i, type;
  n:=Length(Flat(T));  # number of boxes in T
  m:=Maximum(Flat(T)); # biggest number in T
  lengths:=List([1..m], i->[1..Length(T)]*0);
  for row in [1..Length(T)] do
    for c in Collected(T[row]) do
      lengths[c[1]][row]:=c[2];
    od;
  od;
  tabs:=[List([1..Length(T)], row->[])];
  m:=1;
  row:=1;
  type:=1;
  while m<=n do
    len:=Sum(lengths[type]);
    if len>0 then
      newtabs:=[];
      for rows in Multicombinations([m..m+len-1],lengths[type]{[1..row]}) do
        for t in tabs do
          s:=Copy(t);
          for i in [1..row] do
            Append(s[i],rows[i]);
          od;
          Add(newtabs,s);
        od;
      od;
      m:=m+len;
      type:=type+1;
      if row<Length(T) then row:=row+1; fi;
      tabs:=newtabs;
    fi;
  od;
  return tabs;
end;

#M Add [i..n] to row <row> or lower of t[comp] so as to obtain a 
#M standard mu-tableau, where mu is a multipartition of n.       
SPECHT.FillMultiStandardTableau:=function(t,i,n,comp,row,mu)
  local tabs, s, r;
  tabs:=[];
  r:=row;
  while r<=Length(mu[comp]) do
    if Length(t[comp][r])<mu[comp][r] 
      and (r=1 or Length(t[comp][r])<Length(t[comp][r-1]) ) then
      s:=Copy(t);
      Add(s[comp][r],i);
      if i=n then Add(tabs,s);
      else Append(tabs, SPECHT.FillMultiStandardTableau(s,i+1,n,1,1,mu));
      fi;
    fi;
    r:=r+1;
  od;
  if comp<Length(mu) then
    Append(tabs, SPECHT.FillMultiStandardTableau(t,i,n,comp+1,1,mu));
  fi;
  return tabs;
end;

#C return the type of the tableau <tab>; note the slight complication
## because the composition which is the type of <tab> may contain zeros.
TypeTableau:=function(tab)
  tab:=Flat(tab);
  Append(tab, [1..Maximum(tab)]);
  tab:=Collected(tab);
  return tab{[1..Length(tab)]}[2]-1;
end;

#C return the shape of the tableau <tab>
ShapeTableau:=function(tab)
  local shape, l;
  shape:=List(tab, Length);
  l:=Length(shape);
  while shape[l]=0 do 
    Unbind(shape[l]); 
    l:=l-1;
  od;
  return shape;
end;

##############################################################

# return the tableau obtained from <t> by acting on it with the
# PERMUTATION <w> - note this is an honest permutation and not a
# chevie permutation.
ActOnTableau:=function(t,w) local row,i;
  return List(t, row->List(row, i->i^w));
end;

# same as previous but for an arbitrarily deeply nested list
ActOnMultiTableau:=function(t,w) local row;
  if IsList(t) then return List(t, row->ActOnMultiTableau(row,w));
  else return t^w;
  fi;
end;

# return the initial (standard) tableau t^mu of shape mu
InitialTableau:=function(arg) local mu,tab,s, m;
  mu:=Flat(arg);
  tab:=[];
  s:=0;
  for m in mu do
    Add(tab, [s+1..s+m]);
    s:=s+m;
  od;
  return tab;
end;

# return the initial super tableau of type (mu|nu)
InitialSuperTableau:=function(mu,nu)
  local t, i;
  t:=ConjugateTableau(List([1..Length(mu)],i->[1..mu[i]]*0+i));
  for i in [1..Length(nu)] do
    Append(t[1], [1..nu[i]]*0+Sum(mu)+1);
  od;
  return ConjugateTableau(t);
end;

# return the tableau T^mu or T^(mu|nu); that is the unique tableau T of
# type mu (resp. (mu|nu)), such that shape(S)\gedom shape(T) whenever
# S is a tableau of type mu (resp. (mu|nu)).
FirstSemistandardTableauType:=function(arg) local m, mu, nu, t, i;
  if Length(arg)=3 then m:=arg[1]; mu:=arg[2]; nu:=arg[3]; 
  elif Length(arg)=1 then mu:=arg[1]; nu:=[]; m:=Sum(mu);
  fi;
  if not ( Length(arg) in [1,3] and IsComposition(mu) and IsInt(m)
    and IsComposition(nu) ) then
    Error("usage: FirstSemistandardTableauType(mu)\n",
          "          or  FirstSemistandardTableauType(m,mu,nu) "); 
  fi;
  if mu=[] then t:=[[]];
  else t:=List([1..Length(mu)], row->[1..mu[row]]*0+row);
  fi;
  if nu=[] then return t; 
  else
    t:=ConjugateTableau(t);
    nu:=ConjugatePartition(nu);
    for i in [1..Length(nu)] do
      Append(t[1], [1..nu[i]]*0+m+i);
    od;
    return ConjugateTableau(t);
  fi;
end;

#U FirstStandardTableau(mu).................................................
#M Return the first row/standard tableaux of shape mu.                      
FirstStandardTableau:=function(mu) local d, tab, row;
  d:=0;
  tab:=[];
  for row in mu do
    Add(tab, [1..row]+d);
    d:=d+row;
  od;
  return tab;
end;

# return the tableau T^mu or T^(mu|nu); that is the unique tableau T of
# type mu (resp. (mu|nu)), such that shape(T)\gedom shape(S) whenever
# S is a tableau of type mu (resp. (mu|nu)).
LastSemistandardTableauType:=function(arg) 
  local mu, nu, t, m, i, j;
  mu:=arg[1];
  if Length(arg)=2 then nu:=arg[2]; else nu:=[]; fi;
  if not IsComposition(mu) or Length(arg)>2 
    or (Length(arg)=2 and not IsComposition(nu)) then 
    Error("usage: LastSemistandardTableauType(mu [,nu])"); 
  fi;
  t:=[[]];
  for i in [1..Length(mu)] do
    Append(t[1], [1..mu[i]]*0+i);
  od;
  if nu=[] then return t; 
  else
    m:=Sum(mu);
    for i in [1..Length(nu)] do
      for j in [1..nu[i]] do
	if j>Length(t) then t[j]:=[m+i];
	else Add(t[j], m+i);
	fi;
      od;
    od;
    return t;
  fi;
end;

# Given a semistandard tableaux <S> return the first semistandard
# tableau s of "type" <S>.
FirstTableau:=function(S) local type, n, row, a, i;
  type:=Collected(Flat(S));
  S:=Copy(S); n:=-1;
  for a in [1..Length(type)] do
    row:=1;
    while type[a][2]>0 do
      i:=Position(S[row], type[a][1]);
      if IsInt(i) then # a match
	S[row][i]:=n; 
	n:=n-1;
	type[a][2]:=type[a][2]-1;  # one less number to find
      else row:=row+1;
      fi;
    od;
  od;
  return -S;
end;

# Given a semistandard tableaux <S> return the last semistandard
# tableau s of "type" <S>.
LastTableau:=function(S)
  return ConjugateTableau(FirstTableau(ConjugateTableau(S)));
end;

# return the last semistandard mu-tableau of type nu
LastSemistandard:=function(mu,nu)
  local tab, Ki, row, dk, i;
  if not ( Sum(mu)=Sum(nu) and Dominates(mu,nu) ) then return false; fi;
  tab:=[ [1..nu[1]]*0+1 ]; # first row
  for i in [2..Length(nu)] do
    Ki:=nu[i];
    if Length(tab)=Length(mu) then row:=Length(tab);
    else
      Add(tab, [1..Minimum(Ki,mu[Length(tab)+1])]*0+i);
      Ki:=Ki-Length(tab[Length(tab)]);
      row:=Length(tab)-1;
    fi;
    while Ki>0 do
      if Length(tab[row])<mu[row] then
	if row=1 then
	  if Length(tab[1])+Ki>mu[1] then
	    Error("Does not fit!!!");
	  else
	    Append(tab[row], [1..Ki]*0+i);
	    Ki:=0;
	  fi;
	else
	  dk:=Minimum(Ki,Length(tab[row-1])-Length(tab[row]));
	  Append(tab[row], [1..dk]*0+i);
	  Ki:=Ki-dk;
	fi;
      fi;
      row:=row-1;
    od;
  od;
  return tab;
end;

# given a (standard) tableau <t> and a composition <mu> return the
# (row) semistandard tableau mu(t)
SemistandardFromStandard:=function(mu, t)
  local tmu, row, c;
  tmu:=Flat(FirstSemistandardTableauType(mu));
  return List(t, row->List(row, c->tmu[c]));
end;

# given a (standard) tableau <t> and a composition <mu> return the
# (row) semistandard tableau mu(t)
SuperSemistandardFromStandard:=function(m,mu,nu,t)
  local tmu, row, c;
  tmu:=Flat(FirstSemistandardTableauType(m,mu,nu));
  return List(t, row->List(row, c->tmu[c]));
end;

## Given two standard tableau <s> and <t> - which are assumed to be of
## the same shape - return the permuation w such that <t>=<s>w.
MappingPermTableaux:=function(s,t)
  return MappingPermListList(Flat(s),Flat(t));
end;

# return distinguished coset rep corresponding to the tableau t
SymmPermTableau:=function(t) 
  return MappingPermTableaux(InitialTableau(ShapeTableau(t)),t);
end;

# return distinguished coset rep corresponding to the tableau t
CoxeterWordTableau:=function(t) 
  return CoxeterWordSymmPerm(
                MappingPermTableaux(InitialTableau(ShapeTableau(t)),t));
end;

# Insert <p> into the tableau <p> via bumping and record the position
# <q> of the new node in <Q>. This functions changes <P> and <Q>.
BumpTableaux:=function(P,Q, p, q)
  local row, j, k;
  row:=1;
  while row>0 do
    if row>Length(P) then Add(P, [p]); Add(Q,[q]); return; # new last row
    elif p>P[row][Length(P[row])] then 
      Add(P[row], p); Add(Q[row],q);
      return; # finished
    else      # put p into position j and bump k into the next row
      j:=First([1..Length(P[row])], k->p<P[row][k]);
      k:=P[row][j];
      P[row][j]:=p;
      p:=k;
      row:=row+1;
    fi;
  od;
end;

# Given a permutation <w> in Sym(<n>) return the pair of tableaux
# which correspond to <w> under the Robinson-Schensted correspondence.
RobinsonSchenstedCorrespondent:=function(n,w) local i, P, Q;
  P:=[]; Q:=[];
  for i in [1..n] do
    BumpTableaux(P,Q,i^w,i);
  od;
  return rec(P:=P, Q:=Q);
end;

# Remove the largest number from the recording tableau <Q> and unbump,
# and return, the corresponding number from <P>. Both <P> and <Q> are
# changed by this function.
UnbumpTableaux:=function(P,Q) local p,q,row,i,newp;
  # first find <q> in <Q> - assume no repeats
  q:=Maximum(Flat(Q));
  row:=1;
  repeat
    i:=Position(Q[row],q);
    if not IsInt(i) then row:=row+1; fi;
  until IsInt(i);
  if i=1 then Unbind(Q[row]); else Unbind(Q[row][i]); fi;
  # now unbump from <P>
  p:=P[row][i];
  if i=1 then Unbind(P[row]); else Unbind(P[row][i]); fi;
  row:=row-1;
  while row>0 do
    i:=First([Length(P[row]),Length(P[row])-1..1], i->p>P[row][i]);
    newp:=P[row][i];
    P[row][i]:=p;
    p:=newp;
    row:=row-1;
  od;
  return p;
end;

# Return the permuation corresponding to (P,Q) under the RS
# correspondence
ReverseRobinsonSchenstedCorrespondent:=function(arg) local n,w,P,Q;
  if Length(arg)=1 and IsRec(arg[1]) then
    P:=arg[1].P;
    Q:=arg[1].Q;
  elif Length(arg)=2 and ForAll(arg, IsList) then
    P:=arg[1];
    Q:=arg[2];
  fi;
  n:=Length(Flat(P));
  w:=[];
  while P<>[] do
    Add(w, UnbumpTableaux(P,Q));
  od;
  return PermListList([1..n], Reversed(w));
end;

# insert <i> into the tableau <tab> via bumping
BumpTableau:=function(tab, i)
  local row, j, k;
  row:=1;
  tab:=Copy(tab);
  while row>0 do
    if row>Length(tab) then Add(tab, [i]); row:=0; # new last row
    elif i>tab[row][Length(tab[row])] then 
      Add(tab[row], i);
      row:=0; # finished
    else      # put i into position j and bump k into the next row
      j:=First([1..Length(tab[row])], k->i<tab[row][k]);
      k:=tab[row][j];
      tab[row][j]:=i;
      i:=k;
      row:=row+1;
    fi;
  od;
  return tab;
end;

# return the tableau obtained by bumping the sequence
# i1=arg[1],...ik=arg[k] into <tab>=arg[1]
BumpingTableau:=function(arg)
  local tab, i;
  tab:=Copy(arg[1]);
  for i in arg{[2..Length(arg)]} do
    tab:=BumpTableau(tab, i);
  od;
  return tab;
end;

# return the decent set of the tableau <tab>
DescentSetTableau:=function(tab)
  local n, pos, row, i;
  n:=Sum(tab, Length);
  pos:=[1..n]*0;
  for row in [1..Length(tab)] do
    for i in tab[row] do
      pos[i]:=row;
    od;
  od;
  return Filtered([1..n-1], row->pos[row]<pos[row+1]);
end;

##############################################################

#C If <tab> is a tableau its *residue sequence* is (r_1,r_2,...,r_n)
#C This function returns the residue sequence of a tableau <t>; that 
## is the sequence (r_1,r_2,...,r_n) where r_i is the residue of i
## in <tab>.
## ResidueSequence(e|H, t)
ResidueSequence:=function(e, tab) local res, r, c;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not IsInt(e) then Error("usage, ResidueSequence(e|H, t)"); fi;

  res:=[];
  for r in [1..Length(tab)] do
    for c in [1..Length(tab[r])] do
      res[ tab[r][c] ] := (c-r) mod e;
    od;
  od;
  return res;
end;

# Return a list [l_1,l_2,...] where l_i is the number of nodes on
# ladder i on mu. Notice that because the ladders have slope e-1
# the node (r,c) lies on the same ladder as the node (r+(e-1)(c-1),1); 
# we say that this is ladder r+(e-1)(c-1).
LadderSequencePartition:=function(arg)
  local e, mu, ladder, r, c;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, LadderSequencePartition(<H>,<mu>)");
  fi;

  mu:=Flat(arg{[2..Length(arg)]});
  r:=Maximum(List([1..Length(mu)], r->r+(e-1)*(mu[r]-1)));
  ladder:=[1..r]*0;
  for r in [1..Length(mu)] do
    for c in [1..mu[r]] do
      ladder[r+(e-1)*(c-1)]:=ladder[r+(e-1)*(c-1)]+1;
    od;
  od;
  return ladder;
end;

# Sort the entries of <tab> according to the ladder that they lie on
LadderClassTableau:=function(H,tab)
  local e, mu, r, ladders, c;
  if IsRec(H) and IsBound(H.e) then e:=H.e;
  else Error("usage, LadderClassTableau(<H>,<tab>)");
  fi;
  mu:=ShapeTableau(tab);
  r:=Maximum(List([1..Length(mu)], r->r+(e-1)*(mu[r]-1)));
  ladders:=List([1..r], c->[]);
  for r in [1..Length(tab)] do
    for c in [1..Length(tab[r])] do
      AddSet(ladders[r+(e-1)*(c-1)], tab[r][c]);
    od;
  od;
  return ladders;
end;

# determine whether of not the two tableau <s> and <t> are ladder
# equivalent
AreLadderEquivalentTableaux:=function(H,s,t)
  return LadderClassTableau(H,s)=LadderClassTableau(H,t);
end;

####################################################################

# The following functions are used for printing tableaux            
VeryTightStringList:=function(list) local str,a;
  str:="";
  for a in list do
    PrintToString(str, a);
  od;
  return str;
end;

VeryTightTableau:=function(tab) local str,row;
  if tab=[] then return ""; fi;
  str:=VeryTightStringList(tab[1]);
  for row in [2..Length(tab)] do
    PrintToString(str, "/",VeryTightStringList(tab[row]));
  od;
  return str;
end;

VeryTightMultiTableau:=function(tab) local str, t;
  str:=VeryTightTableau(tab[1]);
  for t in [2..Length(tab)] do
    PrintToString(str,":",VeryTightTableau(tab[t]));
  od;
  return str;
end;

#C returns a `tight' string for the tableau <t>
StringTableau:=function(t) local tab, p;
  tab:="[[";
  for p in t do
    PrintToString(tab, TightStringList(p), "],[");
  od;
  tab:=tab{[1..Length(tab)-2]};
  PrintToString( tab, "]" );
  return tab;
end;

# return a string containing TeX codes for printing a tableau
TeXTableau:=function(tab)
  local pre, stab, row;
  pre:="\\tab("; stab:="";
  for row in tab do
    stab:=ConcatenationString(stab, pre, 
	   ApplyFunc(ConcatenationString, List(row,String)));
    pre:=",";
  od;
  return ConcatenationString(stab, ")");
end;

#U Subtableau(tab, k)
#M Return the subtableau of the tableau <tab> which has all M entries 
#M less than <k>. We assume that <tab> is row semistandard; that is,
#M in each row the entries are in weakly increasing order.
Subtableau:=function(tab,k) local row, c;
  return List(tab, row->Filtered(row, c->c<k));
end;

#U DominatesTableau(s,t)
#M Return true if the tableau <s> dominates <t>; that is,
#M    shape(s,k) \gedom shape(t, k) for k=1,2,...,n,
#M where shape(s,k) is the shape of the subteableau of s containing as
#M entries the numbers 1..k. Here <s> and <t> should be (at least) row 
#M semistandard.
DominatesTableau:=function(s,t)
  local n, k;
  n:=Maximum(Flat(s));
  return ForAll([2..n], k->Dominates(ShapeTableau(Subtableau(s,k)),
                                     ShapeTableau(Subtableau(t,k))) );
end;
 
#U MultiStandardTableaux(mu)..........................................
#M Return the list of standard tableaux of shape mu where mu is a 
#M multi-partition.
MultiStandardTableaux:=function(mu)
  local tabs, n, nu, r, l, t;
  if ForAll(mu,c->c=[]) then return [ Copy(mu) ]; fi;
  tabs:=[];
  n:=Sum(Flat(mu));
  for r in Reversed([1..Length(mu)]) do
    # we treat the last removable node differently as doing
    if mu[r]<>[] then
      l:=Length(mu[r]);
      nu:=Copy(mu);
      if mu[r][l]>1 then
        nu[r][l]:=nu[r][l]-1;
        for t in MultiStandardTableaux(nu) do
          Add(t[r][l],n);
          Add(tabs,t);
        od;
      else
        Unbind(nu[r][l]);
        for t in MultiStandardTableaux(nu) do
          Add(t[r],[n]);
          Add(tabs,t);
        od;
      fi;
      while l>1 do
        l:=l-1;
        if mu[r][l]>mu[r][l+1] then
          nu:=Copy(mu);
          nu[r][l]:=nu[r][l]-1;
          for t in MultiStandardTableaux(nu) do
            Add(t[r][l],n);
            Add(tabs,t);
          od;
        fi;
      od;
    fi;
  od;
  return tabs;
end;
