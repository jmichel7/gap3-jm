#######################################################################
##  SPECHT - symmcomb.g : Combinatorial functions on partitions.     ##
##     coming from theory of crystal graphs and Littelmann paths     ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Andrew Mathas                                                 ##
##                                                                   ##
#######################################################################

#U GoodNodes(H|e,mu)......................................................
#U GoodNodes(H|e,mu,I)....................................................
#M Returns the numbers of the rows which end in one of Kleshchev's "good  
#M nodes" (see [LLT] or one of Kleshchev's papers for a description).     
#M Basically, reading from the top down, count +1 for a *removable* node  
#M of residue r and -1 for an *addable* node of residue r. The last       
#M removable r-node with all of these tallies strictly positive is the    
#M (unique) good node of residue r - should it exist.                     
#M If <I> is supplied the number of the row containing the unique good    
#M node of residue I is return.                                           
GoodNodes:=function(arg) local e, mu, I, goodnodes, res, i, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)=3 and IsList(arg[2]) and IsInt(arg[3]) then
    mu:=arg[2]; I:=arg[3];
  else mu:=Flat(arg{[2..Length(arg)]});
  fi;
  if not IsBound(e) or (IsBound(I) and (I<0 or I>=e) ) then
    Error("usage: GoodNodes(<e|H>, mu [, I])\n");
  fi;
  goodnodes:=List([1..e],i->false); ## will hold the good nodes
  res:=List([1..e], i->0);          ## tally of #removable-#addable r-nodes
  for i in [1..Length(mu)] do
    r:=(mu[i]-i) mod e;
    r:=r+1; 
    if i=Length(mu) or mu[i]>mu[i+1] then  ## removable r-node
      if res[r]=0 then goodnodes[r]:=i; 
      else res[r]:=res[r]+1;
      fi;
    fi;
    if r=e then r:=1; else r:=r+1; fi;
    if i=1 or mu[i]<mu[i-1] then           ## addable r-node
      res[r]:=res[r]-1;
    fi;
  od;
  if IsBound(I) then return goodnodes[I+1];
  else return goodnodes;
  fi;
end;

#M Given an e-regular partition mu this function returns the corresponding
## good node sequence (= path is Kleshchev's e-good partition lattice).
##   usage: GoodNodeSequence(e|H, mu);
GoodNodeSequence:=function(arg) local e, mu, goodnodeseq,  row, res, r;
  if IsInt(arg[1])  then e := arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e := arg[1].e;
  else Error("usage, GoodNodeSequence(<H>,<mu>) or ",
              "GoodNodeSequence(<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]}); 
  if not IsERegular(e,mu) then
    Error("GoodNodeSequence(<e>,<mu>): <mu> must be <e>-regular\n");
  fi;
  goodnodeseq:=[];
  while mu<>[] do
    row:=1;
    while row<Length(mu) and mu[row]=mu[row+1] do
      row:=row+1;      ## there is a good node with the same residue as
    od;                ## the first removable node
    r:=(mu[row]-row) mod e;
    res:=0;
    repeat
      if r=(mu[row]-row) mod e and (row=Length(mu) or mu[row]>mu[row+1])
      then
         if res=0 then 
           if mu[row]=1 then Unbind(mu[row]);
           else mu[row]:=mu[row]-1;
           fi;
           Add(goodnodeseq, r);
         else res:=res+1; 
         fi;
       elif r=(mu[row]+1-row) mod e and mu[row]<mu[row-1] then ## addable
         res:=res-1;
       fi;
       row:=row+1;
    until row>Length(mu);
  od;
  return goodnodeseq{[Length(goodnodeseq),Length(goodnodeseq)-1..1]};
end;

#U VeryGoodNodes(e,mu)....................................................
#M Return the set of i-good nodes for <mu> which have no addable i-nodes  
#M above them.                                                            
VeryGoodNodes:=function(arg) local e, mu, goods, adds, vgs, i;
  if IsInt(arg[1]) then e := arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e := arg[1].e;
  else Error("usage, GoodNodeSequence(<H>|<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  goods:=GoodNodes(e,mu);
  adds:=AddableINodes(e,mu);
  vgs:=[];
  for i in [1..e] do
    if IsInt(goods[i]) and (adds[i]=[] or goods[i]<Minimum(adds[i])) 
    then vgs[i]:=goods[i];
    else vgs[i]:=false;
    fi;
  od;
  return vgs;
end;

#U LadderNodeSequences(e,mu)..............................................
#M Return a list ol all of the ldder node sequences from the empty        
#M partition to <mu>.                                                     
LadderNodeSequences:=function(e,mu) local lseqs, ladders, nu, s, l;
  if mu=[] then return [[]]; fi;
  lseqs:=[];
  ladders:=VeryGoodNodes(e,mu);
  for s in [1..e] do
    if IsInt(ladders[s]) then
      nu:=Copy(mu);
      if nu[ladders[s]]=1 then Unbind(nu[ladders[s]]);
      else nu[ladders[s]]:=nu[ladders[s]]-1;
      fi;
      for l in LadderNodeSequences(e,nu) do
	Add(l,s-1);
	Add(lseqs,l);
      od;
    fi;
  od;
  return lseqs;
end;


#U TopGoodNodeSequence(H,mu)..............................................
#U TopGoodNodeSequence(e,mu)..............................................
#M Return a good node sequence for <mu> where at each stage we pick the   
#M highest good node.                                                     
TopGoodNodeSequence:=function(arg) local e, mu, t, gs, r;
  if IsInt(arg[1]) then e := arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e := arg[1].e;
  else Error("usage, GoodNodeSequence(<H>|<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  if mu=[] then return []; fi;
  t:=Minimum(GoodNodes(e,mu));
  if not IsInt(t) then return t; fi;
  r:=(mu[t]-t) mod e;
  if mu[t]=1 then Unbind(mu[t]); else mu[t]:=mu[t]-1;fi;
  gs:=TopGoodNodeSequence(e,mu);
  Add(gs,r);
  return gs;
end;

#M Returns the list of all good node sequences for the partition <mu>
GoodNodeSequences:=function(arg) local e, mu, r, gnss, nu, s, res;
  if IsInt(arg[1]) then e := arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e := arg[1].e;
  else Error("usage, GoodNodeSequence(<H>|<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  if not IsERegular(e,mu) then
    Error("GoodNodeSequence(<e>,<mu>): <mu> must be <e>-regular\n");
  fi;
  if mu=[1] then gnss:=[ [0] ]; 
  else
    gnss:=[];
    for r in GoodNodes(e,mu) do
      if r<>false then
        nu:=Copy(mu);
        nu[r]:=nu[r]-1;
        if nu[r]=0 then Unbind(nu[r]); fi;
        res:=(mu[r]-r) mod e;
        for s in GoodNodeSequences(e,nu) do
          Add(s,res); 
          Add(gnss, s);
        od;
      fi;
    od;
  fi;
  return gnss;
end;

#M Given a good node sequence this function returns the corresponding
## partition, or false if the sequence is not a good node sequence.
##   usage: GoodNodeSequence(H|e, mu)
PartitionGoodNodeSequence:=function(arg) local e, gns, mu, r, i, res, row;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, PartitionGoodNodeSequence(<H>|<e>, <gns>)");
  fi;
  gns:=Flat(arg{[2..Length(arg)]});
  mu:=[];
  for r in gns do
    row:=0;
    res:=0;
    for i in [1..Length(mu)] do
      if r=(mu[i]-i) mod e and (i=Length(mu) or mu[i]>mu[i+1]) and res<0 
      then res:=res+1;
      elif r=(mu[i]+1-i) mod e and (i=1 or mu[i]<mu[i-1]) then
        if res=0 then row:=i; fi;
        res:=res-1;
      fi;
    od;
    if res=0 and r=(-Length(mu))mod e then mu[Length(mu)+1]:=1;
    elif row>0 then mu[row]:=mu[row]+1;
    else return false;  ## bad sequence
    fi;
  od;
  return mu;
end;

#U TableauGoodNodeSequence(e,gns).........................................
#U TableauGoodNodeSequence(H,gns).........................................
#M Given a good node sequence this function returns the corresponding     
#M tableau, or false if the sequence is not a good node sequence.         
TableauGoodNodeSequence:=function(arg) local e, gns, tab, row, res, r, i;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, TableauGoodNodeSequence(<H>|<e>, <gns>)");
  fi;
  gns:=Flat(arg{[2..Length(arg)]});
  tab:=[];
  for r in [1..Length(gns)] do
    row:=0;
    res:=0;
    for i in [1..Length(tab)] do
      if gns[r]=(Length(tab[i])-i) mod e 
	and (i=Length(tab) or Length(tab[i])>Length(tab[i+1])) and res<0 
      then res:=res+1;
      elif gns[r]=(Length(tab[i])+1-i) mod e 
	and (i=1 or Length(tab[i])<Length(tab[i-1])) then
        if res=0 then row:=i; fi;
        res:=res-1;
      fi;
    od;
    if res=0 and gns[r]=(-Length(tab))mod e then tab[Length(tab)+1]:=[r];
    elif row>0 then Add(tab[row],r);
    else return false;  ## bad sequence
    fi;
  od;
  return tab;
end;

#U GoodNodeLatticePath(H|e,mu)............................................
#M Returns one path in the good partition lattice from the empty partition
#M to the partition <mu>.                                                 
GoodNodeLatticePath:=function(arg) local e, i, gns;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, GoodNodeLatticePath(<H>|<e>, <mu>)");
  fi;
  gns:=GoodNodeSequence(e,Flat(arg{[2..Length(arg)]}));
  return List([1..Length(gns)],i->PartitionGoodNodeSequence(e,gns{[1..i]}));
end;

#U GoodNodeLatticePaths(H|e,mu)...........................................
#M Returns the list of all paths in the good partition lattice from the   
#M empty partition to the partition <mu>.                                 
GoodNodeLatticePaths:=function(arg) local e, gns, g, i;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, GoodNodeLatticePaths(<H>|<e>, <mu>)");
  fi;
  gns:=GoodNodeSequences(e,Flat(arg{[2..Length(arg)]}));
  return List(gns, g->List([1..Length(g)],
              i->PartitionGoodNodeSequence(e,g{[1..i]})));
end;

#U LatticePathGoodNodeSequence(<H>|<e>,<gns>).............................
#M Returns the partitions in the path in the e-good partition lattice     
#M which is given by the good node sequence <gns>.                        
LatticePathGoodNodeSequence:=function(arg) local e, gns, i;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, LatticePathGoodNodeSequence(<H>|<e>, <gns>)");
  fi;
  gns:=Flat(arg{[2..Length(arg)]});
  gns:=List([1..Length(gns)],i->PartitionGoodNodeSequence(e,gns{[1..i]}));
  if false in gns then return gns{[1..Position(gns,false)]};
  else return gns;
  fi;
end;

# return the list of nodes which are "strongly normal" in the sense
# that they are (i) removable and (ii) no addable nodes of the same
# residue occur above them
StronglyNormalNodes:=function(e,mu)
  local strong, addables, res, row;
  strong:=List([1..e], i->[]); # will contain the strongly normal nodes
  addables:=[1..e]*0;
  for row in [1..Length(mu)] do
    res:=(mu[row]-row) mod e+1;
    if addables[res]=0 and (row=Length(mu) or mu[row]>mu[row+1]) then 
      Add(strong[res],row);
    fi;
    if row=1 or mu[row]<mu[row-1] then addables[res mod e+1]:=1;
    fi;
  od;
  return strong;
end;

#U NormalNodes(H,mu)......................................................
#U NormalNodes(H,mu,i)....................................................
#M Returns the numbers of the rows which end in one of Kleshchev's "normal
#M nodes" (see [LLT] or one of Kleshchev's papers for a description).     
NormalNodes:=function(arg) local e, mu, I, normalnodes, res, i, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)=3 and IsList(arg[2]) and IsInt(arg[3]) then
    mu:=arg[2]; I:=arg[3] mod e;
  else mu:=Flat(arg{[2..Length(arg)]});
  fi;
  if not (IsBound(e) and IsPartition(mu)) then
    Error("usage: NormalNodes(<e|H>, mu [, I])\n");
  fi;
  normalnodes:=List([1..e],i->[]);  ## will hold the normal nodes
  res:=List([1..e], i->0);          ## tally of #removable-#addable r-nodes
  for i in [1..Length(mu)] do
    r:=(mu[i]-i) mod e;
    r:=r+1; 
    if i=Length(mu) or mu[i]>mu[i+1] then  ## removable r-node
      if res[r]=0 then Add(normalnodes[r],i); 
      else res[r]:=res[r]+1;
      fi;
    fi;
    if r=e then r:=1; else r:=r+1; fi;
    if i=1 or mu[i]<mu[i-1] then           ## addable r-node
      res[r]:=res[r]-1;
    fi;
  od;
  if IsBound(I) then return normalnodes[I+1];
  else return normalnodes;
  fi;
end;

#U NormalNodesAbacusRunners(r,runners)
#M Return the list of normal beads on runner <r> of the list <runners>
#M of abacus runners 
NormalNodesAbacusRunners:=function(r, runners)
  local rminus, count, normals, max, a;
  if r=1 then rminus:=Length(runners);
  else rminus:=r-1;
  fi;
  count:=0;
  normals:=[];
  max:=Maximum(Flat([runners[r], runners[rminus]]));
  for a in [max,max-1..0] do
    if a in runners[r] and not (a in runners[rminus]) then # removable
      if count>=0 then Add(normals, a); fi;
      if count<0 then count:=count+1; fi;
    elif not a in runners[r] and a in runners[rminus] then #addable
      count:=count-1;
    fi;
  od;
  return normals;
end;

#U NormalAndConormalNodes(e,mu)
#M Return a list [rec_1,...,rec_e] or records, where each record has two  
#M components: a list normals of rows ending in normal nodes, a list      
#M conormals of rows ending in conormal nodes.                            
NormalAndConormalNodes:=function(e,mu) local cns, l, r;
  cns:=List([1..e],i->rec(normals:=[],conormals:=[]));
  l:=Length(mu);
  r:=-l mod e+1;  # addable node at bottom of mu.
  Add(cns[r].conormals,l);
  while l>0 do
    r:=(mu[l]-l) mod e +1;
    if l=Length(mu) or mu[l]>mu[l+1] then # removable
      Add(cns[r].normals,l);
    fi;
    if r=e then r:=1;
    else r:=r+1;
    fi;
    if l=1 or mu[l]<mu[l-1] then           # removable
      if cns[r].normals<>[] then Unbind(cns[r].normals[Length(cns[r].normals)]);
      else Add(cns[r].conormals,l);
      fi;
    fi;
    l:=l-1;
  od;
  return cns;
end;


#U ConormalNodes(H,mu)....................................................
#U ConormalNodes(H,mu,i)..................................................
#M Returns the numbers of the rows which end in a conormal node. That is, 
#M an addable node which has at least as many addable nodes as removable  
ConormalNodes:=function(arg) local e, mu, I, conormalnodes, res, i, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)=3 and IsList(arg[2]) and IsInt(arg[3]) then
    mu:=arg[2]; I:=arg[3] mod e;
  else mu:=Flat(arg{[2..Length(arg)]});
  fi;
  if not (IsBound(e) and IsPartition(mu)) then
    Error("usage: ConormalNodes(<e|H>, mu [, I])\n");
  fi;
  conormalnodes:=List([1..e],i->[]);## will hold the conormal nodes
  res:=List([1..e], i->0);          ## tally of #addable-#removable r-nodes
  AddSet(conormalnodes[(-Length(mu)) mod e+1],Length(mu)+1);
  # last addable node is coconormal
  for i in Reversed([1..Length(mu)]) do
    r:=((mu[i]-i) mod e)+1;         ## the residue of the node
    if i=Length(mu) or (i<Length(mu) and mu[i]>mu[i+1]) then  ## removable r-node
      res[r]:=res[r]-1;             ## increment acount as removable node
    fi;
    if r=e then r:=1; else r:=r+1; fi;     ## the residue of the `added' node
    if i=1 or mu[i]<mu[i-1] then           ## addable r-node
      if res[r]=0 then AddSet(conormalnodes[r], i);
      else res[r]:=res[r]+1;
      fi;
    fi;
    #Print("i = ", i, ", r = ", r,", res = ", res, "\n");
  od;
  if IsBound(I) then return conormalnodes[I+1];
  else return conormalnodes;
  fi;
end;

## usage: RemoveNormalNodes(H|e, mu, i)
## returnsthe partition obtained from <mu> by removing all the normal
## nodes of residue <i>.
RemoveNormalNodes:=function(e, mu, I) local res,i,r;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not ( IsInt(e) and IsList(mu) and IsInt(I) and I<e and I>=0) then 
    Error("usage, RemoveNormalNodes(e|H, mu, I)"); 
  fi;
  mu:=Copy(mu);               ## we are going to change this so...
  res:=0;                     ## tally of #removable-#addable I-nodes
  for i in [1..Length(mu)] do
    r:=(mu[i]-i) mod e;
    if r=I and (i=Length(mu) or mu[i]>mu[i+1]) then  ## removable I-node
      if res=0 then mu[i]:=mu[i]-1;                  ## normal I-node
      else res:=res+1;
      fi;
    fi;
    if r=e then r:=1; else r:=r+1; fi;
    if r=I and (i=1 or mu[i]<mu[i-1]) then           ## addable I-node
      res:=res-1;
    fi;
  od;
  return mu;
end;

#U UpPartition(e, mu)
#M Let beta be a set of beta numbers for the partition <mu>. Define       
#M           U = {x in beta | x - e \not\in beta}.                        
#M Then <mu> is an e-core iff U=\emptyset, in which case set Up(mu)=mu.   
#M If U\ne\emptyset then let u=Max(U) and define                          
#M    V ={ v>p | v-e \in beta, v \notin beta and v\not\equiv u mod e}.    
#M Put v=Min(V) and define Up(mu) to be the partition with beta set       
#M beta-{u}\cup{v}. This, and the four functions which follow. are the    
#M key combinatorial definitions in the paper of Ariki, Kreiman and       
#M Tsuchioka which gives a non-recursive description of the Kleshchev     
#M multipartitions.                                                       
UpPartition:=function(e,mu)
  local beta, U, u, v, b;
  if not IsInt(e) then e:=e.e; fi;
  beta:=BetaSet(mu);
  U:=Filtered(beta, b->b-e>=0 and not b-e in beta);
  if U=[] then return mu; fi;
  u:=Maximum(U);
  v:=u;
  repeat
    v:=v+1;
  until Mod(v-u,e)>0 and v-e in beta and not v in beta;
  RemoveSet(beta,u);
  AddSet(beta,v);
  return PartitionBetaSet(beta);
end;

#U RoofPartition(e,mu)....................................................
#M Return the UpPartition-stable partition   which is obtained from <mu>  
#M by applying UpPartition() enough times.                                
RoofPartition:=function(e,mu) local nu, upnu;
  if not IsInt(e) then e:=e.e; fi;
  upnu:=Copy(mu);
  repeat
    nu:=upnu;
    upnu:=UpPartition(e,upnu);
  until upnu=nu;
  return nu;
end;

#U DownPartition(e, mu)
#M Let beta be a set of beta numbers for the partition <mu>. Define       
#M           U = {x in beta | x - e \not\in beta}.                        
#M Then <mu> is an e-core iff U=\emptyset, in which case set Down(mu)=mu. 
#M If U\ne\emptyset then let u=Min(U)-e and define                        
#M           W ={ w>u | w \in beta, w+e \notin U }.                       
#M Put v=Min(V) and define Down(mu) to be the partition with beta set     
#M beta-{u}\cup{v}
DownPartition:=function(e,mu)
  local beta, U, u, v, b;
  if not IsInt(e) then e:=e.e; fi;
  beta:=BetaSet(mu);
  U:=Filtered(beta, b->b-e>=0 and not b-e in beta);
  if U=[] then return mu; fi;
  u:=Minimum(U)-e;
  v:=u;
  repeat
    v:=v+1;
  until (v in beta and not v+e in beta) or v=u+e;
  if v>u+e then v:=v+e; fi;
  RemoveSet(beta,v);
  AddSet(beta,u);
  return PartitionBetaSet(beta);
end;

#U BasePartition(e,mu)....................................................
#M Return the DownPartition-stable partition   which is obtained from <mu>
#M by applying DownPartition() enough times.                              
BasePartition:=function(e,mu) local nu, downnu;
  if not IsInt(e) then e:=e.e; fi;
  downnu:=Copy(mu);
  repeat
    nu:=downnu;
    downnu:=DownPartition(e,downnu);
  until downnu=nu;
  return nu;
end;

#U TauPartition(e,m,mu)...................................................
#M Return tau_{e,m}(mu), as defined by Ariki, Kreiman and Tsucioka. That  
#M is, if beta={\mu_i-i+m|i\ge0} is the set of m-beta numbers then define 
#M M_i(mu)=max{b in beta|b\equiv i mod e}. Then tau_{e,m}(mu) is the      
#M partition with beta numbers beta\cup{M_{i_1}(mu)+e,...,M_{i_m}(mu)+e}, 
#M where we order the M_i(mu) so that M_{i_1}(mu)>...>M_{i_e}(mu).        
TauPartition:=function(e,m,mu)
  local beta, M, i, b;
  beta:=List([1..Length(mu)], i->mu[i]-i); # start with "charge" 0
  M:=[];
  for i in [1..e]+Length(mu) do
    M[(-i)mod e+1]:=-i;
  od;
  for b in beta do
    i:=b mod e;
    if M[i+1]<b then M[i+1]:=b; fi;
  od;
  M:=Reversed(Set(M));
  Append(beta, M{[1..m]}+e);
  beta:=Reversed(Set(beta));
  return List([1..Length(beta)], i->beta[i]+i-m); # extract as charge m
end;

#U CrystalWeightPartition(e,mu)...........................................
#M Given a positive integer <e> and an e-regular partition <mu> return the
#M "crystal weight" of mu. This is the sequence [w_0,...,w_{e-1}] where   
#M w_i is the number of i-nodes in the diagram of <mu>
CrystalWeightPartition:=function(e,mu) local wt, alpha, i, row, col;
  # we set wt = Lambda_0
  wt:=[1..e]*0; wt[1]:=1;
  alpha:=List([1..e], i->[1..e]*0);
  if e=2 then alpha:=[[2,-2],[-2,2]];
  else
    for i in [1..e] do
      alpha[i][i]:=2;
      if i>1 then alpha[i][i-1]:=-1; else alpha[i][e]:=-1; fi;
      if i<e then alpha[i][i+1]:=-1; else alpha[i][1]:=-1; fi;
    od;
  fi;
  for row in [1..Length(mu)] do
    for col in [1..mu[row]] do
      i:=(col-row)mod e;
      wt:=wt-alpha[i+1];
    od;
  od;
  return wt;
end;

#U WeylGroupActionPartition(e,w,mu).......................................
#M Given <e>, a word <w> in {0,1,...,e-1} and an e-regular partition <mu> 
#M return w*mu, where the action of s_i on mu is given by                 
#M         s_i mu = { \tilde f_i^{ wt(mu)(h_i) mu, if wt(mu)(h_i)\ge0     
#M                  { \tilde e_i^{-wt(mu)(h_i) mu, if wt(mu)(h_i)\le0     
WeylGroupActionPartition:=function(e,w,mu)
  local wmu, wt, normals, conormals, i, x;
  if not IsERegular(e,mu) then return false; fi;
  if IsInt(w) then w:=[w]; fi;  # allow w=i, rather than a word
  wmu:=Copy(mu);
  for i in w do
    wt:=CrystalWeightPartition(e,wmu);
    if wt[i+1]>0 then
      conormals:=ConormalNodes(e,wmu)[i+1];
      for x in conormals{[1..wt[i+1]]} do
        if IsBound(wmu[x]) then wmu[x]:=wmu[x]+1;
        else wmu[x]:=1;
        fi;
      od;
    else
      normals:=NormalNodes(e,wmu)[i+1];
      for x in normals{[1..-wt[i+1]]} do
        if wmu[x]>1 then wmu[x]:=wmu[x]-1;
        else Unbind(wmu[x]);
        fi;
      od;
    fi;
  od;
  return wmu;
end;

#C We now set up a record for computing the crystal amplitude maps. These 
#C are the maps which embed the crystal B(Lambda_i) into B(h\Lambda_i). If
#C mu is an e-regular partition and [i_1,...,i_n) is a good node sequence 
#C for 
CrystalAmplitudeOps:=OperationsRecord("Crystal amplitude");
CrystalAmplitudeOps.Print:=function(self) local i;
  Print(" (",LabelPartition(self.cores[1]),")^",self.rat[1]);
  for i in [2..self.h] do
    Print(" x (",LabelPartition(self.cores[i]),")^",self.rat[i]);
  od;
  Print("\n");
end;

#U AmplitutdeMap(e,h,mu)
AmplitutdeMap:=function(e,h,mu) local i;
  if IsECore(e,mu) then 
    return rec(e:=e, h:=h, rat:=[1..h]*0+1/h, cores:=List([1..h],i->mu),
               operations:=CrystalAmplitudeOps);
  fi;
end;
