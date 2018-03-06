#######################################################################
##  SPECHT - symmcomb.g : Combinatorial functions on partitions.     ##
##                                                                   ##
##     This file contains most of the combinatorial functions used   ##
##     by Specht. Most are standard operations on Young diagrams     ##
##     or partitions.                                                ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Andrew Mathas                                                 ##
##                                                                   ##
#######################################################################

#D Combinatorial functions coming from symmetric groups and Hecke algebras

## 2.4: October 1997:
##  - added funtions MullineuxSymbol, PartitionMullineuxSymbol,
##    NormalNodes (plus undocumented friends), BetaSet, PartitionBetaSet.

## 2.2: June 1996:
##  - mainly change of function names to make it more compatible with
##    GAP conventions.

## 2.1: April 1996:
##  - added functions for finding paths in the good node partition
##    lattice.

## 2.0: March 1996 : symmcomb.g file created, by breaking specht.g
##  - added functions for finding Kleshchev's "good nodes" and implemented
##    his (and James') algorithm for the Mullineux map.

## 1.0: December 1995: initial release.

######################################################################

#U Lexicographic(mu,nu)                                                   
#M Return true if <mu> appears before <nu> in the lexicographic ordering  
#M on the set of partitions.                                              
Lexicographic:=function(lambda,mu) return lambda=mu or lambda>mu; end;

#U LengthLexicographic(mu,nu)                                             
#M Return true if <mu> appears before <nu> in the total order on the      
#M set of partitions where first we compare the lengths of <mu> and <nu>  
#M and trhen we compare then lexicographically. By default this is used   
#M to sort the rows of DecompositionMatrix().                             
LengthLexicographic:=function(mu,nu)
  if Length(mu)=Length(nu) then
    return mu=nu or mu>nu;
  else return Length(mu)<Length(nu);
  fi;
end;

#U LengthReverseLexicographic(mu,nu)                                      
#M The total ordering of the set of partitions given by first             
#M comparing lengths and then reverse lexicographically.                  
LengthReverseLexicographic:=function(mu,nu)
  if Length(mu)=Length(nu) then
    return mu=nu or Reversed(mu)<Reversed(nu);
  else return Length(mu)<Length(nu);
  fi;
end;

#U ReverseDominance(nu,mu)                                                
#M The total ordering of the set of partitions given by reverse           
#M dominance; that is, nu is greater than mu if                           
#M     sum_{i>k} nu_i \ge sum_{i>k mu_i,   for all k.                     
#M Surprisingly, this is a total order.                                   
ReverseDominance:=function(nu,mu) local i, Mu, Nu;
  if mu=nu then return true;
  elif Length(nu)=Length(mu) then
    i:=Length(mu);
    Mu:=0; Nu:=0;
    while i > 0 do
      Mu:=Mu + mu[i]; Nu:=Nu + nu[i];
      if Nu < Mu then return true;
      elif Nu > Mu then return false;
      fi;
      i:=i - 1;
    od;
  else return Length(nu)<Length(mu);
  fi;
end; #ReverseDominance

#U Dominates(mu,nu)
#M The partial ordering of the set of partitions given by 
#M dominance; that is, nu is greater than mu if
#M     sum_{i>=1}^k nu_i \ge sum_{i=1}^k mu_i,   for all k.
Dominates:=function(mu, nu) local i, m, n, lnu, lmu;
  if nu=mu then return true;
  elif nu>mu then return false;
  fi;
  m:=0; n:=0; # partial sums
  i:=1; lnu:=Length(nu); lmu:=Length(mu);
  while i <=lnu and i<=lmu do
    m:=m + mu[i]; n:=n + nu[i];
    if n > m then return false; fi;
    i:=i + 1;
  od;
  return true;
end;  # Dominates

#U ConjugatePartition(mu)
#U ConjugatePartition(mu_1,mu_2,...)
#M Returns the conjugate which is conjugate to the composition <mu>;
#M that is, the partition whose i-th entry is the number of boxes in
#M the i-th column of the diagram of <mu>.
ConjugatePartition:=function(arg) 
  local comp, M, con, i, j;
  comp:=Flat(arg);
  if comp=[] then return []; fi;
  M:=Maximum(comp);
  con:=[1..M]*0;
  for i in comp do
    for j in [1..i] do
      con[j]:=con[j]+1;
    od;
  od;
  return con;
end; # ConjugatePartition

#U ConjugateMultiPartition(mu)...........................................
#M Return the conjugate multipartition to <mu>.                          
ConjugateMultiPartition:=function(mu)
  return List(Reversed(mu),ConjugatePartition);
end;

#U ReverseConjugateLexicographic(mu,nu)
#M The partial ordering of the set of partitions given by 
#M mu>nu if nu'>mu', where ' denotes the conjugate partition.
ReverseConjugateLexicographic:=function(mu,nu)
  return Lexicographic(ConjugatePartition(nu),ConjugatePartition(mu));
end;

#U LengthLexicographicallyByQuotient(e,mu,nu).............................
#M The total on partitions given by ordering their <e>-quotients length   
#M lexicographically. For compatibiltiy with the other orderings we cannot
#M just set H.Ordering=LengthLexiographicallyByQuotient as this function  
#M has too many arguments. You can get around this by defining a suitable 
#M wrapper yourself or passing this order function to DecompositionMatrix 
#M to set it as the ordering (it will make the wrapped function for you). 
#M This function is most useful when mu and nu are in the same block.     
LengthLexicographicallyByQuotient:=function(e,mu,nu)
  local core, w, beads, qmu, qmus, qnu, qnus, q;
  # We want to make sure that the largest partition in the block has
  # quotient ((0)...,(0),(w)), so we are going to compute everything
  # by hand - slowly, unfortunately - and construct the required abacuses. 
  # We assume that mu and nu are in the same block.
  if mu=nu then return true; fi;
  core:=ECore(e,mu);
  w:=(Sum(mu)-Sum(core))/e; # the e-weight => need >= len(core)+w beads
  beads:=Length(core)+w*e;
  qmu:=List(EAbacusRunners(e,mu,beads), q->PartitionBetaSet(Set(q)));
  qmus:=List(qmu, Sum);
  qnu:=List(EAbacusRunners(e,nu,beads), q->PartitionBetaSet(Set(q)));
  qnus:=List(qnu, Sum);
  q:=1;
  if qmus[1]=0 and qnus[1]>0 then return true;
  elif qmus[1]>0 and qnus[1]=0 then return false;
  elif qmus[1]=0 then  # so qnus[1]>0 too
    if qmus=qnus then
      q:=e;
      while q>1 and qnu[q]=qmu[q] do q:=q-1; od;
      return LengthLexicographic(qmu[q],qnu[q]);
    else return Reversed(qmus)>Reversed(qnus);
    fi;
  else  # qmu[1]=qnu[1]=[]
    if qmus=qnus then
      q:=First([1..e], q->qmu[q]<>qnu[q]);
      return LengthLexicographic(qnu[q],qmu[q]);
    else return Reversed(qnus)<Reversed(qmus);
    fi;
  fi;
end;

#U QuotientOrder(qmu,qnu).................................................
#M Return true if the quotient qmu beats the quotient qnu in the ordering 
#M |qmu^(1)|+...+qmu^(i)+\sum_{j=1}^k qmu^(i+1)_j                         
#M             \ge |qnu^(1)|+...+qnu^(i)+\sum_{j=1}^k qnu^(i+1)_j,        
#M for all i and k. This is a partial order.                              
QuotientOrder:=function(qmu,qnu) local i, muSum, nuSum, j;
  muSum:=0;
  nuSum:=0;
  i:=1;
  while i<=Minimum(Length(qmu),Length(qnu)) do
    j:=1;
    while j<=Minimum(Length(qmu[i]),Length(qnu[i])) do
      muSum:=muSum+qmu[i][j]; 
      nuSum:=nuSum+qnu[i][j]; 
      if nuSum>muSum then return false; fi;
      j:=j+1;
    od;
    if j<Length(qmu[i]) then 
      muSum:=muSum+Sum(qmu[i]{[j..Length(qmu[i])]});
    fi;
    if j<Length(qnu[i]) then 
      nuSum:=nuSum+Sum(qnu[i]{[j..Length(qnu[i])]});
      if nuSum>muSum then return false; fi;
    fi;
    i:=i+1;
  od;
  return true;
end;

#U LittlewoodRichardsonRule(alpha,beta)...................................
#M Return the list of partitions, with multiplicity, obtained by applying 
#M the Littlewood-Richardson rule to the partitions <alpha> and <beta>.   
#M That is, return the set of partitions <nu> such that                   
#M     Ind_{ Sym(m) x Sym(n) }^{ Sym(m+n) }( S(alpha) x S(beta) )         
#M                \bigoplus_{nu} S(nu)                                    
#M where |alpha|=m, |beta|=n and S(mu) is the Specht module indexed by    
#M <mu>.                                                                  
#C The algorithm has (at least), one obvious improvement in that it should
#C collect like terms using something like H.operations.Collect after     
#C wrapping on each row of beta.                                          
LittlewoodRichardsonRule:=function(alpha, beta)
  local lrr, newlrr, x, i, j, row, Place, max, newbies;

  # place k nodes in row r>1 and above; max is the maximum number
  # of new nodes which may be added on this row and below (so
  # this is dependent upon p.new=the number of nodes added to a
  # given row from the previous row of beta).
  Place:=function(p, k, r, max) local newp, np, m, i, M;

    if r > Length(p.lam) then  # top of the partition
      Add(p.lam, k); Add(p.new, k); return [ p ]; else
      if r > 1 and p.lam[r]=p.lam[r-1] then
        max:=max + p.new[r];
        p.new[r]:=0;
        return Place(p, k, r+1, max);
      else
        if r=1 then            # m number of nodes that can be new
          m:=Minimum(k, max);  # to row r
        else m:=Minimum(p.lam[r-1]-p.lam[r], k, max);
        fi;
        if m >=0 and k-m <=p.lam[r] then  # something may fit
          newp:=[];
          for i in [0..m] do  # i nodes on row r
            if k-i <=p.lam[r] then      # remaining nodes will fit on top
              M:=max - i + p.new[r];    # uncovered nodes in previous rows
              np:=Copy(p);
              if k-i=0 and m > 0 then   # all gone
                np.lam[r]:=np.lam[r] + i;
                np.new[r]:=i;
                Add(newp, np);
              else                      # more nodes can still be placed
                np.new[r]:=i;
                for np in Place(np, k-i, r+1, M) do
                  np.lam[r]:=np.lam[r] + i;
                  Add(newp, np);
                od;
              fi;
            fi;
          od;
          return newp;
        fi;
        return [];
      fi;
    fi;
  end;  # end of Place; LRR internal

  if alpha=[] or alpha=[0] then return [ beta ];
  elif beta=[] or beta=[0] then return [ alpha ];
  elif Length(beta)*Sum(beta) > Length(alpha)*Sum(alpha) then
    return LittlewoodRichardsonRule(beta, alpha); 
  else
    lrr:=Place(rec(lam:=Copy(alpha),       # partition
                 new:=List(alpha, i->0)),  # new nodes added from this row
                 beta[1], 1, beta[1]);
    for i in [2..Length(beta)] do
      newlrr:=[];
      for x in lrr do
        row:=1;
        while x.new[row]=0 do row:=row + 1; od;
        max:=x.new[row];
        x.new[row]:=0;
        Append(newlrr, Place(x, beta[i], row+1, max));
      od;
      lrr:=newlrr;
    od;
    return List(lrr, x->x.lam);
  fi;
end;  # LittlewoodRichardsonRule

#U InverseLittlewoodRichardsonRule(lambda)................................
#U InverseLittlewoodRichardsonRule(lambda,mu).............................
#M The first form returns all pairs (sigma,tau) with multiplicity         
#M a^\lambda_{sigma,tau}=LR coefficient. The second form returns only     
#M those pairs of the form (mu,tau), again with multiplicity.             
InverseLittlewoodRichardsonRule:=function(arg)
  local initialise,fill,alpha,start,n,invlr,npp,new,mu,row,newp,r,max,l,p,x;

  initialise:=function(p, r) local M, np, newp, i;
    if r=1 then newp:=[ ]; M:=alpha[1];
    else newp:=[ Copy(p) ]; M:=Minimum(alpha[r], p[r-1]);
    fi;
    for i in [1..M] do
      np:=Copy(p);
      np[r]:=i;
      if r < Length(alpha) then Append(newp, initialise(np, r+1));
      else Add(newp, np);
      fi;
    od;
    return newp;
  end;

  fill:=function(p, row, r, max) local m, M, np, newp, i, x;
    newp:=[];
    m:=Minimum(Minimum(p.total[r-1],alpha[r])-p.total[r], max);
    if row > 1 then m:=Minimum(m, p.mu[row-1]-p.mu[row]); fi;
    max:=max + p.new[r];
    for i in [0..m] do
      np:=Copy(p);
      np.new[r]:=i;
      np.mu[row]:=np.mu[row] + i;
      if r=Length(alpha) then np.total[r]:=np.total[r] + i; Add(newp, np);
      else
        for x in fill(np, row, r+1, max-i) do
           x.total[r]:=x.total[r] + i; Add(newp, x);
        od;
      fi;
    od;
    return newp;
  end;

  if Length(arg)=2 and ForAll(arg,IsList) then 
    alpha:=arg[1];
    if not Dominates(alpha, arg[2]) or Length(alpha)<Length(arg[2]) then return [];
    elif alpha=arg[2] then return [ [alpha,[]] ];
    fi;
    start:=[ arg[2] ];
    invlr:=[ ];
  elif Length(arg)<>1 then
    Error("Usage: InverseLittlewoodRichardsonRule(lambda)\n",
          "    or InverseLittlewoodRichardsonRule(lambda,mu),");
  elif arg[1]=[] then return [ [[],[]] ];
  else
    alpha:=Flat(arg);
    start:=initialise([],1);
    invlr:=[ [ [], alpha ] ];
  fi;
  n:=Sum(alpha);
  for l in start do
    npp:=[rec(total:=Copy(l), new:=List(alpha, r -> 0), mu:=[])];
    for r in [Length(npp[1].total)+1..Length(alpha)] do
      npp[1].total[r]:=0;
    od;
    row:=1;
    while npp<>[] do
      newp:=[];
      for p in npp do
        if row > 1 then r:=row - 1;
        else r:=1;
        fi;
        max:=0;
        while r < Length(p.total) and p.total[r]=alpha[r] do
          max:=max + p.new[r]; p.new[r]:=0; r:=r + 1;
        od;
        p.mu[row]:=alpha[r] - p.total[r];
        if row=1 or p.mu[row] <=max then
          if row=1 then max:=p.total[1];
          else max:=max + p.new[r] - p.mu[row];
          fi;
          p.new[r]:=p.mu[row];
          if r < Length(alpha) then
            for x in fill(p, row, r+1, max) do
              x.total[r]:=x.total[r] + p.mu[row];
              if Sum(x.total)=n then Add(invlr, [l, x.mu]); 
              else Add(newp, x);
              fi;
            od;
          else Add(invlr, [l, p.mu]);
          fi;
        fi;
      od;
      row:=row + 1;
      npp:=newp;
    od;
  od;
  if not ( Length(arg)=2 and ForAll(arg,IsList) ) then
    invlr[Length(invlr)]:=[alpha,[]];   ## rough hack...
  fi;
  return invlr;
end;

#U LittlewoodRichardsonCoefficient(mu,nu,lambda)..........................
#M Return the Littlewood-Richardson coefficient a^lamdba_{mu,nu}.         
LittlewoodRichardsonCoefficient:=function(mu,nu,lambda) local x;
  if Sum(lambda)<>Sum(mu)+Sum(nu) then return 0;
  else return Length(Filtered(InverseLittlewoodRichardsonRule(lambda,mu),
                       x->x=[mu,nu]));
  fi;
end;

#U SpechtDimension(mu)....................................................
#U SpechtDimension(x).....................................................
#M Returns the dimension of the Specht module S(<mu>), or the linear      
#M combination of Specht modules <x>.                                     
SpechtDimension:=function(arg) local Dim,y;
  Dim:=function(mu) local mud, i,j,d;
    mud:=ConjugatePartition(mu);
    d:=Factorial(Sum(mu));
    for i in [1..Length(mu)] do
      for j in [1..mu[i]] do
        d:=d/(mu[i] + mud[j] - i - j + 1);
      od;
    od;
    return d;
  end;

  if Length(arg)=1 and IsSpecht(arg[1]) then
    return Sum([1..Length(arg[1].coeffs)],
               y->arg[1].coeffs[y]*Dim(arg[1].parts[y]));
  else return Dim(Flat(arg));
  fi;
end; #SpechtDimension

#U BetaNumbers(mu)........................................................
#U BetaNumbers(mu,l)......................................................
#M Returns a list of beta numbers for the partition mu. In the second form
#M the list of beta numbers of length <l> is returned.                    
BetaNumbers:=function(arg) local mu, l, beta;
  if Length(arg)=1 then mu:=arg[1]; l:=Length(mu);
  elif Length(arg)=2 and IsList(arg[1]) then 
    mu:=arg[1]; 
    l:=arg[2];
  else 
    mu:=Flat(arg);
    l:=Length(mu);
  fi;
  if not IsPartition(mu) or l<Length(mu) then 
    Error("usage, BetaNumbers(mu [,l])");
  fi;
  if mu=[] then beta:=[0];
  else beta:=mu + [Length(mu)-1, Length(mu)-2..0];
  fi;
  if Length(beta)<l then 
    l:=l-Length(beta);
    beta:=beta+l;
    Append(beta, [l-1,l-2..0]);
  fi;
  return beta;
end;

## returns a set of the beta numbers for the partition mu
## JM: conflicts slightly with BetaSet in standard library ctsymmet.g
#BetaSet:=function(mu)
#  if mu=[] then return [0];
#  else return Reversed(mu) + [0..Length(mu)-1];
#  fi;
#end;

## given a beta set return the corresponding partition
PartitionBetaSet:=function(beta) local i;
  if beta=[] or beta[Length(beta)]=Length(beta)-1 then return []; fi;
  if beta[1]=0 then 
    i:=First([1..Length(beta)], i->beta[i]>i-1);
    beta:=beta-i+1;
    beta:=beta{[i..Length(beta)]};
  fi;
  beta:=beta-[0..Length(beta)-1];
  return Reversed(beta);
end;

#U EAbacusRunners(e, mu)
#U EAbacusRunners(e, mu, length])
#M Return the <e>-abacus runners for the partition <mu>. That is, a
#M list of length <e> where each component of the list gives the beta
#M numbers for the beads which appear on the corresponding runners of the
#M <e>-abacus of <mu>. In the second form, <length> specifies the total 
#M number of beads on the abacus.
EAbacusRunners:=function(arg) local e, mu, beta, extraBeads, b, aba, i;
  e:=arg[1];
  mu:=arg[2];
  ## first we find a set of beta numbers for mu; we want an e-multiple
  ## of (strictly) decreasing beta numbers for mu.
  beta:=BetaNumbers(mu);
  if Length(arg)=3 then extraBeads:=Maximum(0,arg[3]-Length(beta));
  else extraBeads:=-Length(beta) mod e;  ## -> an e-multiple of beads
  fi;
  aba:=List([1..e], i->[]);
  if extraBeads>0 then                ## add extraBeads beads to beta
    beta:=beta+extraBeads;
    Append(beta,[extraBeads-1,extraBeads-2..0]);
  fi;
  for i in beta do
    Add(aba[ (i mod e)+1 ], Int(i/e) );
  od;
  return aba;
end;

#U FullEAbacusRunners(e, mu)
#M Given a partition <mu> and an integer <e> return the corresponding
#M <e>-abacus runners for <mu>. Here <mu> can be either a partition or
#M a record (eg. an AbacusSymbopl) containing an EAbacusRunner operation.
FullEAbacusRunners:=function(e, mu)
  local beta, core, w, l, runners, r;
  if IsRec(mu) and IsBound(mu.operations) and
       IsBound(mu.operations.FullEAbacusRunners) 
  then return mu.operations.FullEAbacusRunners(e, mu);
  elif not IsPartition(mu) then
    Error("I don't know how to construct an abacus for <mu>");
  fi;
  # mu is a partition
  beta:=BetaNumbers(mu);
  core:=ECore(e,mu);
  w:=(Sum(mu)-Sum(core))/e;  # the e-weight
  l:=Length(core)+w*e;
  if Length(beta)<l then
    beta:=beta+l-Length(beta);
    Append(beta,[l-Length(beta)-1,l-Length(beta)-2..0]);
  fi;
  runners:=List([1..e], r->[]);
  for r in beta do
    Add(runners[ r mod e + 1], Int(r/e) );
  od;
  return List(runners, Set);
end;

#U PartitionAbacusRunners(abacus)
#M Given an e-abacus return the corresponding partition.
PartitionAbacusRunners:=function(abacus) local e, r;
  e:=Length(abacus);
  return PartitionBetaSet(Set(Flat(List([0..e-1], r->r+e*abacus[r+1]))));
end; 

#U ECore(e|H,mu)
#U ECore(e|H,mu_1,mu_2,...)
#M Return the <e>-core of the partition <mu>.
ECore:=function(arg) local core, beta, i, j, e;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, ECore(<H>,<mu>), ECore(<e>,<mu>)");
  fi;

  if e=0 then return Flat(arg{[2..Length(arg)]}); fi;  
  beta:=List(EAbacusRunners(e,Flat(arg{[2..Length(arg)]})), i->Length(i));
  beta:=beta - Minimum(beta);  # remove all fully occupied rows
  if Maximum(beta)=0 then return [];
  else
    ## at present beta contains the number of beads on each runner of
    ## the abacus. next we get the beta numbers for all of the beads
    core:=[];
    for i in [1..e] do
      Append(core, List([1..beta[i]], j->e*(j-1)+i-1));
    od;
    Sort(core);
    if core[1]=0 then ## remove the beads which don't affect the beta numbers
      if core[Length(core)]=Length(core)-1 then return []; fi; #empty
      i:=First([1..Length(core)], i->core[i]<>i-1);
      core:=core{[i..Length(core)]}-i+1;
    fi;
    
    ## finally, we unravel the beta numbers of our core
    core:=List([1..Length(core)],i->core[i]-i+1);
    return core{[Length(core),Length(core)-1..1]};
  fi;
end;

#U IsECore(e|H,mu)
#U IsECore(e|H,mu)1,mu_2,...)
#M Returns true if <mu> is an <e>-core, and false othewise.
IsECore:=function(arg) local e, r;
  if IsInt(arg[1]) then e:=arg[1];
  elif IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, IsECore(<H>,<mu>), EWeight(<e>,<mu>)");
  fi;

  return ForAll(EAbacusRunners(e,Flat(arg{[2..Length(arg)]})), 
                          r->r=[] or Length(r)=r[1]+1);
end;

#U EWeight(e|H,mu)
#U EWeight(e,mu_1,mu_2,...)
#M Returns the <e>-weight of the partition <mu>.
EWeight:=function(arg) local e, r, wt;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EWeight(<H>,<mu>), EWeight(<e>,<mu>)");
  fi;

  # this is slightly more efficient than computing (Sum(mu)-Sum(ECore(e,mu))/e
  if e=0 then return 0;
  else return Sum(EAbacusRunners(e,Flat(arg{[2..Length(arg)]})),
                  r->Sum(r)-Length(r)*(Length(r)-1)/2);
  fi;
end;

#U EQuotient(e|H,mu)
#U EQuotient(e|H,mu_1,mu_2,...)
#M Return the <e>-quotient of the partition <mu> partition. 
EQuotient:=function(arg) local e, q, mu, d, i, j, qj, x;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EQuotient(<H>,<mu>), EQuotient(<e>,<mu>)");
  fi;

  if e=0 then return []; fi;
  mu:=Flat(arg{[2..Length(arg)]});
  d:=ConjugatePartition(mu);
  q:=List([1..e], j->[]);
  # we use the star of mu to compute the quotient
  for i in [1..Length(mu)] do
    x:=0;
    qj:=(mu[i]-i) mod e;
    for j in [1..mu[i]] do
      if (j-d[j]-1) mod e=qj then x:=x + 1; fi;
    od;
    if x<>0 then Add(q[qj+1], x); fi;
  od;
  return q;
end;  # EQuotient

#U CanonicalQuotient(e|H,mu)..............................................
#U CanonicalQuotient(e|H,mu_1,mu_2,...)...................................
#M Return the canonical <e>-quotient of the partition <mu> partition where
#M the components of the quotient are ordered via `Richards ordering' of  
#M the abacus; that is, according to the positions of the last bead on    
#M each runner.                                                           
CanonicalQuotient:=function(arg) local e, mu, abacus, order;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EQuotient(<H>,<mu>), EQuotient(<e>,<mu>)");
  fi;
  if e=0 then return []; fi;
  mu:=Flat(arg{[2..Length(arg)]});
  abacus:=List(EAbacusRunners(e, mu),Set);
  order:=[1..e];
  SortParallel(List([1..e],i->i-1+e*Length(abacus[i])),order);
  return List(order, i->PartitionBetaSet(abacus[i]));
end;  # EQuotient

#U CombineCanonicalQuotientECore(quot,core)...............................
#M Return the partition which has canonical quotient <quot> and e-core    
#M <core.                                                                 
CombineCanonicalQuotientECore:=function(quot,core)
  local e, coreabacus, order, quotbetas, diff, i;
  e:=Length(quot);
  # get an abacus for <core> with enough beads to fit <quot>
  coreabacus:=EAbacusRunners(e,core,Length(core)+e*Sum(Flat(quot)));
  # get beta numbers for the quotient
  quotbetas:=List(quot, BetaNumbers);
  # find out the ordering of quot
  order:=[1..e];
  SortParallel(List([1..e],i->i-1+e*Length(coreabacus[i])),order);
  # now put <quot> onto coreabacus
  for i in [1..e] do
    if Length(quotbetas[i])<Length(coreabacus[order[i]]) then
      # pad out the beta numbers in each component to the right length
      diff:=Length(coreabacus[order[i]])-Length(quotbetas[i]);
      quotbetas[i]:=quotbetas[i]+diff;
      Append(quotbetas[i], [diff-1,diff-2..0]);
    fi;
    coreabacus[order[i]]:=quotbetas[i];
  od;
  return PartitionAbacusRunners(coreabacus);
end;

#U LabelCanonicalQuotient(e,mu)...........................................
#M Return a string label using the canonical quotient for the partition   
#M <mu>.                                                                  
LabelCanonicalQuotient:=function(e,mu) local quot, label, i;
  if not IsInt(e) then e:=e.e; fi;
  quot:=CanonicalQuotient(e,mu);
  label:="";
  for i in [1..e] do
    if quot[i]<>[] then
      PrintToString(label,RealLabelPartition(quot[i]));
    fi;
    if i<e then PrintToString(label,"|"); fi;
  od;
  return label;
end;

#U TeXCanonicalQuotient(e,mu)...........................................
#M Return a TeXable string for the canonical quotient for the partition 
#M <mu>. This string has &s separating the components.                  
TeXCanonicalQuotient:=function(e,mu) local quot, label, i;
  if not IsInt(e) then e:=e.e; fi;
  quot:=CanonicalQuotient(e,mu);
  label:="";
  for i in [1..e] do
    if quot[i]<>[] then
      PrintToString(label,RealLabelPartition(quot[i]));
    fi;
    if i<e then PrintToString(label,"&"); fi;
  od;
  return label;
end;

#U TeXQuotient(quot)...................................................
#M Return a TeXable string for a quotient <quot>. This string has |s   
#M separating the components.                                          
TeXQuotient:=function(quot) local sep, label, mu;
  sep:="";
  label:="\\<";
  for mu in quot do
    PrintToString(label,sep);
    if mu<>[] then
      PrintToString(label,RealLabelPartition(mu));
    fi;
    sep:="|";
  od;
  PrintToString(label,"\\>");
  return label;
end;

#U CanonicalQuotientDominance(mu,nu)....................................
#M This ordering only really makes sense on a given block, however, we  
#M extend this to all partitions by using dominance on the cores. If    
#M <mu> and <nu> are in the same block then we return true if           
#M    \sum_{s=1}^{t-1} |mu^s|+\sum_{j=1}^i mu^s_j                       
#M           \ge \sum_{s=1}^{t-1} |nu^s|+\sum_{j=1}^i nu^s_j            
#M for 1<=t<=e and all i\ge1. Otherwise we return false.                
CanonicalQuotientDominance:=function(e,mu,nu)
  local cmu, cnu, m, n, s, row;
  if not IsInt(e) then e:=e.e; fi; ## assume OK...
  if mu=nu then return true; fi;
  cmu:=ECore(e,mu);
  cnu:=ECore(e,nu);
  if cmu<>cnu then return Dominates(cmu,cnu); fi;
  cmu:=CanonicalQuotient(e,mu);
  cnu:=CanonicalQuotient(e,nu);
  m:=0; n:=0; # cumulative sums
  s:=1;
  while s<=Length(cmu) do
    row:=1;
    while row<=Length(cmu[s]) and row<=Length(cnu[s]) do
      m:=m+cmu[s][row];
      n:=n+cnu[s][row];
      if n>m then return false; fi;
      row:=row+1;
    od;
    if row<Length(cmu[s]) then m:=m+Sum(cmu[s]{[row..Length(cmu[s])]});
    elif row<Length(cnu[s]) then 
      n:=n+Sum(cnu[s]{[row..Length(cnu[s])]});
      if n>m then return false; fi;
  fi;
    s:=s+1;
  od;
  return true;
end;

#U CanonicalQuotientLengthLexicographic(e,mu,nu)
#M A total order that looks like it is comptabible with the partial
#M order above.
CanonicalQuotientLengthLexicographic:=function(e,mu,nu)
  local cmu, cnu, s;
  if not IsInt(e) then e:=e.e; fi; ## assume OK...
  if mu=nu then return true; fi;
  cmu:=ECore(e,mu);
  cnu:=ECore(e,nu);
  if cmu<>cnu then return Dominates(cmu,cnu); fi;
  cmu:=CanonicalQuotient(e,mu);
  cnu:=CanonicalQuotient(e,nu);
  s:=1;
  while s<=Length(cmu) do
    if cmu[s]=cnu[s] then s:=s+1;
    elif Length(cmu[s])<Length(cnu[s]) then return true;
    elif Length(cmu[s])=Length(cnu[s]) and cmu[s]>cnu[s] then return true;
    else return false;
    fi;
  od;
  Error("We should never end up here");
end;

EAbacusOperations:=OperationsRecord("AbacusOps");
EAbacusOperations.Print:=function(abacus)
  local m, j, i;
  if abacus.mu=[0] then
    for j in [1..abacus.e] do Print("  ."); od;
    Print("\n\n");
  else
    m:=Maximum(Flat(abacus.runners)) + 1;
    for i in [0..m] do
      for j in [1..abacus.e] do
        if  i in abacus.runners[j] then Print("  0");
        else Print("  .");
        fi;
      od;
      Print("\n");
    od;
  fi;
end;
EAbacusOperations.TeXps:=function(abacus)
  local str, m, j, i;
  str:=SPrint("\\abacus{\n  ");
  for j in [1..abacus.e] do PrintToString(str," &"); od; 
  PrintToString(str,"\\\\\n");
  if abacus.mu=[0] then
    PrintToString(str,"  ");
    for j in [1..abacus.e] do PrintToString(str,"-&"); od;
    PrintToString(str,"\\\\\n");
    m:=0;
  else
    m:=Maximum(Flat(abacus.runners)) + 1;
    for i in [0..m] do
      PrintToString(str,"  ");
      for j in [1..abacus.e] do
        if  i in abacus.runners[j] then PrintToString(str," &");
        else PrintToString(str,"-&");
        fi;
      od;
      PrintToString(str,"\\\\\n");
    od;
  fi;
  PrintToString(str,"  "); 
  for j in [1..abacus.e] do PrintToString(str,"\\,&");od;
  PrintToString(str,"\\\\\n");
  for j in [1..abacus.e] do 
    PrintToString(str,"\\ncline{1,",j,"}{",m+3,",",j,"}\n");
  od;
  PrintToString(str,"}\n");
  return str;
end;
# using the tikz abacus macro
EAbacusOperations.TeX:=function(abacus) local beta, i, rows;
  beta:=[];
  for i in [1..abacus.e] do
    Append(beta, abacus.runners[i]*abacus.e+i-1);
  od;
  if beta=[] then
    return SPrint("\\abacus{",abacus.e,"}{",2,"}{",TightStringList([0..abacus.e-1]),"}");
  else
    rows:=Maximum(Flat(abacus.runners))+1;
    return SPrint("\\abacus{",abacus.e,"}{",rows,"}{",TightStringList(Set(beta)),"}");
  fi;
end;

#U EAbacus(e|H,mu)
#U EAbacus(e|H,mu,beads)
#U EAbacus(e|H,mu_1,mu_2,...)
#M Prints the <e>-abacus for the partition <mu>. The number of beads
#M on the abacus is the smallest possible integer multiple of <e>, except
#M if the syntax of the second call is used when the specifed number of
#M beads is used.
EAbacus:=function(arg) local e, mu, beads, core, w, operations;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, EAbacus(<H>,<mu>), EAbacus(<e>,<mu>)");
  fi;
  if Length(arg)=3 and IsPartition(arg[2]) then
    mu:=arg[2];
    beads:=arg[3];
  else
    mu:=Flat(arg{[2..Length(arg)]});
    core:=ECore(e,mu);
    w:=(Sum(mu)-Sum(core))/e;
    beads:=Length(core)+w*e;
  fi;
  if mu=[] then mu:=[0]; fi;
  return rec(e:=e, runners:=EAbacusRunners(e,mu,beads), 
             mu:=mu,
             operations:=EAbacusOperations);
end;

#U FullEAbacus(H,mu)
#M Print the e-abacus for <mu> which has at least <w> beads on 
#M each runner and exactly w beads on the first runner. Here w is
#M the e-weight of <mu>, <mu> is a partition and <H> is either
#M a positive integer or a Specht record.
FullEAbacus:=function(arg)
  local e, mu, beta, w, l, B, b;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, FullAbacus(<H>,<mu>), FullAbacus(<e>,<mu>)");
  fi;
  mu:=Flat( arg{[2..Length(arg)]} );
  beta:=BetaNumbers(mu);
  w:=EWeight(e,mu);
  l:=Length(ECore(e,mu))+w*e;
  if Length(beta)<l then
    beta:=beta+l-Length(beta);
    Append(beta,[l-Length(beta)-1,l-Length(beta)-2..0]);
  fi;
  # the beta number at the end of the last row
  B:=e*Int(beta[1]/e)+e-1;
  for b in [0..B] do
    if b mod e =0 then Print("\n    "); fi;
    if b in beta then Print(" O"); else Print(" ."); fi;
  od;
  Print("\n\n");
end;

#U FullEAbacuses(H,mu,nu)
#M Prints the full e-abcuses for <mu> and <nu> side by side.
FullEAbacuses:=function(e,mu,nu)
  local mubeta, w, nubeta, rows, r, b;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not IsInt(e) then
    Error("usage, FullEAbacuss(<e>,<mu>,<nu>)");
  fi;
  mubeta:=BetaNumbers(mu);
  w:=EWeight(e,mu);
  mubeta:=mubeta+w*e;
  Append(mubeta,[w*e-1,w*e-2..0]);
  ## we assume that nu and mu are in the same block
  nubeta:=BetaNumbers(nu);
  b:=Length(mubeta)-Length(nubeta);
  nubeta:=nubeta+b;
  Append(nubeta,[b-1,b-2..0]);
  # the maximum number of rows-1 in the two abacuses
  rows:=Maximum(Int(mubeta[1]/e), Int(nubeta[1]/e));
  for r in [0..rows] do
    Print("\n    ");
    for b in [e*r..e*(r+1)-1] do
      if b in mubeta then Print(" O"); else Print(" ."); fi;
    od;
    Print("    ");
    for b in [e*r..e*(r+1)-1] do
      if b in nubeta then Print(" O"); else Print(" ."); fi;
    od;
  od;
  Print("\n\n");
end;  

#U CombineEQuotientECore(quot,core);  
#M Returns the unique partition which has <e>-quotient <quot> and
#M <e>-core <core>.
CombineEQuotientECore:=function(q, c) local e, aba, m, beta, i, j;
  e:=Length(q);
  aba:=EAbacusRunners(e,c); # abacus with an e-multiple of runners to which
                            # we need to add m beads to fit the quotient
  m:=Maximum(List([1..e], i->Length(q[i])-Length(aba[i])));
  m:=Maximum(m, 0);
  beta:=[];
  for i in [1..e] do
    if q[i]<>[] then
      q[i]:=q[i] + Length(aba[i]) + m;
      for j in [1..Length(q[i])] do Add(beta, (q[i][j]-j)*e + i - 1); od;
    fi;
    for j in [1..Length(aba[i])+m-Length(q[i])] do
      Add(beta, (j-1)*e + i - 1);
    od;
  od;
  Sort(beta);
  if Length(beta)>0 and beta[1]=0 then  ## remove irrelevant beta numbers; see ECore()
    if beta[Length(beta)]=Length(beta)-1 then return []; fi;
    m:=First([1..Length(beta)],i->beta[i]<>i-1);
    beta:=beta{[m..Length(beta)]}-m+1;
  fi;
  beta:=List([1..Length(beta)], i->beta[i]-i+1);
  return beta{[Length(beta),Length(beta)-1..1]};
end;  # CombineEQuotientECore

#U IsERegular(e|H,mu)                                                     
#U IsERegular(e|H,mu_1,mu_2,...)                                          
#M Returns true if <mu> is an <e>-regular partition.                      
IsERegular:=function(arg) local mu, e, i;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsRec(arg[1]) and IsBound(arg[1].e) then  e:=arg[1].e;
  else Error("usage, IsERegular(<H>,<mu>) or IsERegular(<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  if e=0 then return false;
  else ## assume that mu is ordered
    e:=e-1;
    return ForAll([1..Length(mu)-e], i->mu[i]<>mu[i+e]);
  fi;
end;

#U ERegularPartitions(e,n)                                                
#M Returns a list of the <e>-regular partitions of <n>.                   
ERegularPartitions:=function(e,n) local p,i;
  if IsRec(e) and IsBound(e.e) then e:=e.e; 
  elif not IsInt(e) then
    Error("usage, ERegularPartitions(<H>,<mu>) or ",
          "ERegularPartitions(<e>,<mu>)");
  fi;
  if n<2 then return [ [n] ];
  elif e=0 then return Partitions(n);
  fi;
  e:=e-1;
  return Filtered(Partitions(n),
                  p->ForAll([1..Length(p)-e], i->p[i]<>p[i+e]) );
end;

#U IsERestrictedPartition(e|H,mu)                                         
#U IsERestrictedPartition(e|H,mu_1,mu_2,...)                              
#M Returns true if <mu> is an <e>-restricted partition.                   
IsERestrictedPartition:=function(arg) local mu, e, i;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsRec(arg[1]) and IsBound(arg[1].e) then  e:=arg[1].e;
  else Error("usage, IsERestrictedPartition(<e|H>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  if e=0 then return false;
  elif mu=[] then return true;
  else ## assume that mu is ordered
    return mu[Length(mu)]<e and ForAll([1..Length(mu)-1], i->mu[i]-mu[i+1]<e);
  fi;
end;

#U ERestrictedPartitions(e,n)                                             
#M Returns a list of the <e>-regular partitions of <n>.                   
ERestrictedPartitions:=function(e,n) local mu,i;
  if IsRec(e) and IsBound(e.e) then e:=e.e; 
  elif not IsInt(e) then
    Error("usage, ERestrictedPartitions(<e|H>,<mu>)");
  fi;
  if n=0 then return [[]];
  else return Filtered(Partitions(n), mu->mu[Length(mu)]<e 
          and ForAll([1..Length(mu)-1], i->mu[i]-mu[i+1]<e) );
  fi;
end;

#U PartitionsOfGivenWeight(e,n,w)                                         
#M Return the list of the partitions of <n> which have <e>-weight <w>.    
PartitionsOfGivenWeight:=function(e, n, w) local mu;
  return Filtered(Partitions(n), mu->EWeight(e,mu)=w);
end;

#U ERegularPartitionsOfGivenWeight(e,n,w)                                 
#M Return the list of the partitions of <n> which have <e>-weight <w>.    
ERegularPartitionsOfGivenWeight:=function(e, n, w) local mu;
  return Filtered(ERegularPartitions(e,n), mu->EWeight(e,mu)=w);
end;

#U EResidueDiagram(e,mu) 
#U EResidueDiagram(x) 
#M The second form prints the residue diagrams of the e-regular partitions 
#M in x
EResidueDiagram:=function(arg) local e, rs, r, x, PrintEResidueDiagram;
  PrintEResidueDiagram:=function(e,mu) local i, j, res;
    if e=0 then res:=r->r;
    else res:=r ->r mod e;
    fi;
    if mu=[] then Print("\n");
    else
      for i in [1..Length(mu)] do
        for j in [1..mu[i]] do 
          Print(String(res(j-i),4)); 
        od;
        Print("\n");
      od;
    fi;
  end;

  if IsSpecht(arg[1]) or (Length(arg)=2 and IsSpecht(arg[2])) then
    if IsSpecht(arg[1]) then x:=arg[1]; else x:=arg[2]; fi;
    rs:=ListERegulars(x);
    if rs=[] or IsInt(rs[1]) then PrintEResidueDiagram(x.H.e, rs);
    else
      for r in rs do
        if r[1]<>1 then Print(r[1],"*"); fi;
        Print(r[2],"\n");
        PrintEResidueDiagram(x.H.e, r[2]);
      od;
      if Length(rs) > 1 then
        Print("# There are ", Length(rs), " ", x.H.e, 
                "-regular partitions.\n");
      fi;
    fi;
  else 
    if IsInt(arg[1]) then e:=arg[1];
    elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
    else Error("usage, EResidueDiagram(<H>,<mu>), EResidueDiagram(<e>,<mu>)");
    fi;

    PrintEResidueDiagram(e,Flat(arg{[2..Length(arg)]}));
  fi;
end; # EResidueDiagram

#U EResidueSequenceTableau(H,t)
#M return the residue sequence of the tableau t; that is, the sequence
#M ( res_t(1), res_t(2), ..., ,res_t(n) )
EResidueSequenceTableau:=function(H,t)
  local res, r, c;
  res:=[];
  for r in [1..Length(t)] do
    for c in [1..Length(t[r])] do
      res[ t[r][c] ] := (c-r) mod H.e;
    od;
  od;
  return res;
end;

#U EResidueSequencesPartition(H,mu)
#M return the set of all residues sequences of the standard mu-tableau.
EResidueSequencesPartition:=function(H,mu)
  local res, t;
  res:=[];
  for t in StandardTableaux(mu) do
    AddSet(res, EResidueSequenceTableau(H,t));
  od;
  return res;
end;

#M Returns the partion obtained from mu by pushing nodes to the top
## of their e-ladders (see [JK]; there the notation is mu^R).
ERegularizationPartition:=function(arg) local e, mu, ladder, r, c, C, k;
  if IsInt(arg[1]) then e:=arg[1]; 
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e; 
  else Error("usage, ERegularizationPartition(<H>,<mu>),\n",
             "  or   ERegularizationPartition(<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  ladder:=List(mu,r->List([1..r],c->0));
  for r in [2..Length(mu)] do  
    for c in [1..mu[r]] do
      k:=r-(e-1)*Int(r/(e-1));
      if k<1 then k:=k+e-1; fi;
      while k<r do
        C:=c+(r-k)/(e-1);
        if IsBound(ladder[k][C]) then k:=k+e-1;
        else 
          ladder[k][C]:=0; 
          Unbind(ladder[r][c]);
          k:=r;
        fi;
      od;
    od;
  od;
  return List(Filtered(ladder,r->Length(r)>0), r->Length(r));
end;  ## ERegularizationPartition

#P Print the valuation of the hook lengths in the diagram of <mu>
#M Gordon James calls this the p-power diagram in the Sym(n) case.
EPowerDiagram:=function(arg) local valuation, mu, mud, i, j;
  if IsRec(arg[1]) and IsBound(arg[1].valuation) then
    valuation:=arg[1].valuation;
  else Error("usage, EPowerDiagram(<H>,<mu>), EPowerDiagram(<e>,<mu>)");
  fi;
  mu:=Flat(arg{[2..Length(arg)]});
  mud:=ConjugatePartition(mu);
  for i in [1..Length(mu)] do
    for j in [1..mu[i]] do
        Print("  ", valuation(mu[i]+mud[j]-i-j+1) );
    od;
    Print("\n");
  od;
end;

#P hook length diagram
HookLengthDiagram:=function(arg) local mu, mud, i, j;
  mu:=Flat(arg);
  mud:=ConjugatePartition(mu);
  for i in [1..Length(mu)] do
    for j in [1..mu[i]] do
      Print(String(mu[i]+mud[j]-i-j+1, 4));
    od;
    Print("\n");
  od;
end; #HookLengthDiagram

# return an array containing the hook lengths in [mu]
HookLengthsPartition:=function(arg) local mu, mud, hooks, i, j;
  mu:=Flat(arg);
  mud:=ConjugatePartition(mu);
  hooks:=[];
  for i in [1..Length(mu)] do
    hooks[i]:=[];
    for j in [1..mu[i]] do
      hooks[i][j]:=mu[i]+mud[j]-i-j+1;
    od;
  od;
  return hooks;
end; #HookLengthDiagram

# return the list of rows which end in an addable node
AddableNodes:=function(arg) local mu, addables, i;
  mu:=Flat(arg);
  if not ForAll(mu, IsInt) then Error("usage: AddableNodes(<mu>)");fi;
  addables:=[1];
  for i in [2..Length(mu)] do
    if mu[i-1]>mu[i] then AddSet(addables, i); fi;
  od;
  AddSet(addables, Length(mu)+1);
  return addables;
end;

# return the list of rows which end in an addable i-node
AddableINodes:=function(e,mu) local addables, i;
  if IsRec(e) then e:=e.e; fi;
  addables:=List([1..e], i->[]);
  if mu=[] then addables[1]:=[1];
  else Add(addables[(mu[1]+1-1)mod e+1], 1);
  fi;
  for i in [2..Length(mu)] do
    if mu[i-1]>mu[i] then AddSet(addables[(mu[i]+1-i)mod e+1], i); fi;
  od;
  AddSet(addables[(-Length(mu)) mod e+1], Length(mu)+1);
  return addables;
end;

# return the list of rows which end in a removable node
RemovableNodes:=function(arg) local mu, removables, i;
  mu:=Flat(arg);
  if not ForAll(mu, IsInt) then Error("usage: RemovableNodes(<mu>)");fi;
  removables:=[];
  for i in [1..Length(mu)-1] do
    if mu[i]>mu[i+1] then Add(removables, i); fi;
  od;
  Add(removables, Length(mu));
  return removables;
end;

# return the list of rows which end in a removable i-node
RemovableINodes:=function(e,mu) local mu, removables, i;
  if IsRec(e) then e:=e.e; fi;
  removables:=List([1..e], i->[]);
  for i in [1..Length(mu)-1] do
    if mu[i]>mu[i+1] then Add(removables[(mu[i]-i) mod e+1], i); fi;
  od;
  Add(removables[(mu[Length(mu)]-Length(mu)) mod e+1], Length(mu));
  return removables;
end;

# remove all nodes of residue i
AddAllINodes:=function(e,mu,i) local r;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not IsInt(e) or (i<0 or i>e) then
    Error("usage: AddAllINodes(e,mu,i)");
  fi;
  mu:=Copy(mu);
  for r in [1..Length(mu)] do
    if (r=1 or mu[r]<mu[r-1]) and (mu[r]+1-r) mod e=i then 
      mu[r]:=mu[r]+1; 
    fi;
  od;
  if (-Length(mu)) mod e =i then Add(mu,1); fi;
  return mu;
end;

# remove all nodes of residue i
RemoveAllINodes:=function(e,mu,i) local r;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not IsInt(e) or (i<0 or i>e) then
    Error("usage: RemoveAllINodes(e,mu,i)");
  fi;
  mu:=Copy(mu);
  for r in [1..Length(mu)] do
    if (r=Length(mu) or mu[r]>mu[r+1]) and (mu[r]-r) mod e =i then 
      mu[r]:=mu[r]-1; 
      if mu[r]=0 then Unbind(mu[r]); fi;
    fi;
  od;
  return mu;
end;

#U MullineuxSymbol(<H>|<e>, <mu>).........................................
#M returns the Mullineux symbol of the <e>-regular partition <mu>. The    
#M algorithm is basically to shuffle the first column hooks lengths; this 
#M is a reformulation of Mullineux's approach.                            
#M E.g. if e=3 and mu=[4,3,2] then we do the following:                   
#M    betanums =  [6, 4, 2, 0]                                            
#M             -> [4, 3, 1, 0] :6->4, 4->3 ( we want 2 but can only       
#M                                           remove 1 more node as e=3 )  
#M             -> [3, 2, 1, 0].                                           
#M To get the Mullineux symbols we record the number of beads removed at  
#M each stage and also the number of signiciant numbers in the previous   
#M beta number (i.e. the numebr of rows); here we get                     
#M                  5, 3, 1                                               
#M                  3, 2, 1                                               
MullineuxSymbol:=function(arg) 
  local e, mu, betaset, newbetaset, tally, difference,i,ms;

  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  else Error("usage, MullineuxSymbol(<H>|<e>, <mu>)");
  fi;

  mu:=arg{[2..Length(arg)]};
  if mu=[] or mu=[0] then return [ [0],[0] ];
  elif IsList(mu[1]) then mu:=mu[1]; 
  fi;
  betaset:=BetaSet(mu);
  ms:=[ [],[] ];
  while betaset<>[] do
    newbetaset:=Copy(betaset);
    RemoveSet(newbetaset, newbetaset[Length(newbetaset)]);
    AddSet(newbetaset,0);
    difference:=betaset-newbetaset;
    tally:=0;
    Add(ms[1], 0);
    Add(ms[2], Length(betaset));
    for i in [Length(betaset),Length(betaset)-1..1] do
      tally:=tally+difference[i];
      if tally>=e then
        newbetaset[i]:=newbetaset[i]+tally-e;
        ms[1][Length(ms[1])]:=ms[1][Length(ms[1])]+e;
        tally:=0;
      fi;
    od;
    ms[1][Length(ms[1])]:=ms[1][Length(ms[1])]+tally;
    betaset:=newbetaset;
    if not IsSet(betaset) then return false; fi; ## can happen?
    if betaset[1]=0 then
      if betaset[Length(betaset)]=Length(betaset)-1 then 
        betaset:=[];
      else
        i:=First([1..Length(betaset)], i->betaset[i]<>i-1);
        betaset:=betaset{[i..Length(betaset)]}-i+1;
      fi;
    fi;
  od;
  return ms;
end;

#M given a Mullineux Symbol <ms> and an integer <e>, return the corresponding 
#M <e>-regular partition.
PartitionMullineuxSymbol:=function(arg) local e, ms, betaset, i,tally,betaN;
  if IsInt(arg[1])  then e:=arg[1];
  elif IsRec(arg[1]) and IsBound(arg[1].e) then e:=arg[1].e;
  fi;
  if Length(arg)<>2 or not IsBound(e) then
    Error("usage, MullineuxSymbol(<H>|<e>, <mu>)");
  fi;
  ms:=Copy(arg[2]);

  betaset:=[0..ms[2][1]-1];
  ms[2]:=ms[2][1]-ms[2]+1;  # significant numbers in betaset
  i:=Length(ms[1]);
  while i>0 do
    tally:=0;
    betaN:=ms[2][i];
    repeat
      if tally=0 then
        tally:=ms[1][i] mod e;
        if tally=0 then tally:=e; fi;
        ms[1][i]:=ms[1][i]-tally;
      fi;
      if betaN=Length(betaset) then 
        betaset[betaN]:=betaset[betaN]+tally;
        tally:=0;
      else
        if betaset[betaN+1]-betaset[betaN]>tally then
          betaset[betaN]:=betaset[betaN]+tally;
          tally:=0;
        else
          tally:=tally-betaset[betaN+1]+betaset[betaN];
          betaset[betaN]:=betaset[betaN+1];
        fi;
      fi;
      betaN:=betaN+1;
    until (tally=0 and ms[1][i]=0) or betaN>Length(betaset);
    if tally>0 or ms[1][i]>0 then return false; fi;
    i:=i-1;
  od; ## while
  return PartitionBetaSet(betaset);
end;

#M removes the rim hook from mu which corresponding to the 
## (row,cols)-th hook.
RemoveRimHook:=function(arg) local mu, row, col, mud, r, c, x, nx;
  mu:=Copy(arg[1]);
  row:=arg[2];
  col:=arg[3];
  if Length(arg)=3 then mud:=ConjugatePartition(mu);
  else mud:=arg[4];
  fi;
  r:=mud[col];
  x:=col;
  while r >=row do
    nx:=mu[r];
    if x=1 then Unbind(mu[r]);
    else mu[r]:=x - 1;
    fi;
    x:=nx;
    r:=r - 1;
  od;
  return mu;
end;

#U AddRimHook(nu, row, h)                                                 
#M Given a partition <nu>, a row n number <row> and an integer <h> return 
#M the the pair [mu, l] where mu is the partition obtained by adding a    
#M rim hook with foot in row <row>, of length of length <h>, and l is     
#M the leg length of the added hook. If the resulting diagram is not the  
#M diagram of a partition then false is returned.                         
AddRimHook:=function(nu, row, h) local r;
  nu:=Copy(nu);
  r:=row;
  if r=Length(nu) + 1 then nu[r]:=0;
  elif r > Length(nu) then h:=0;
  fi;
  while r > 1 and h > 0 do
    h:=h-nu[r-1]+nu[r]-1;
    if h > 0 then 
      nu[r]:=nu[r-1]+1; 
      r:=r-1;
    elif h < 0 then 
      nu[r]:=h+nu[r-1]+1;
    fi;
  od;
  if h > 0 then nu[1]:=nu[1] + h; r:=1;
  elif h=0 then return false;
  fi;
  return [nu, row - r];
end;

#U RemoveNodeFromRow(mu, row)                                            
#M Remove the last node from row <row> of the partition <mu>.            
#M No check is made as to whether or not this node is removable.         
RemoveNodeFromRow:=function(mu, row) 
  local nu;
  nu:=Copy(mu);
  nu[row]:=nu[row]-1;
  if nu[row]=0 then Unbind(nu[row]); fi;
  return nu;
end;

## return a list of the compositions of k into m parts
## (allowing zero parts)
PartitionsOfLength:=function(k,m) local mu;
  return Filtered(Partitions(k), mu->Length(mu)<=m);
end;

## return a list of the compositions of k into m parts
## (allowing zero parts)
CompositionsOfLength:=function(k,m) local l, mcomps, comps,c;
  if m=1 then return [[k]]; fi;  
  comps:=[];
  for l in [0..k] do
    mcomps:=CompositionsOfLength(l,m-1);
    for c in mcomps do
      Add(c, k-l);
    od;
    Append(comps, mcomps);
  od;
  return comps;
end;

# return the cycle of the permutation <w> as a partition of <n>
CycleTypePermutation:=function(n, w) local type;
  type:=List(Cycles(w, [1..n]), Length);
  Sort(type);
  return Reversed(type);
end;

#U CycleTypeTypeBWord(word)..............................................
#M Given a "word" for a Coxeter element of type B return its "cycle type"
#M which is a bipartition.                                               
CycleTypeTypeBWord:=function(n,word) local type, i, j;
  type:=[[],[]];
  i:=1;
  while i<=Length(word) do
    j:=i;
    if word[i]<0 then 
      while j<Length(word) and word[j+1]=j-i-word[i] do
        j:=j+1;
      od;
      Add(type[1],j-i+1);
    else 
      while j<Length(word) and word[j+1]=j-i+word[i]+1 do
        j:=j+1;
      od;
      Add(type[2],j-i+2);
    fi;
    i:=j+1;
  od;
  Sort(type[1]); type[1]:=Reversed(type[1]);
  Sort(type[2]); type[2]:=Reversed(type[2]);
  j:=n-Sum(Flat(type));
  if j>0 then Append(type[2], [1..j]*0+1); fi;
  return type;
end;

#U PermutationCycleType(mu)...............................................
#M Given a composition <mu>  return a permutation which ash <mu> as its   
#M cycle type.                                                            
PermutationCycleType:=function(arg)
  local mu, s, w, i, j;
  mu:=Flat(arg);
  if not IsComposition(mu) then Error("usage: PermutationCycleType(mu"); fi;
  s:=0; w:=();
  for i in [1..Length(mu)] do
    if mu[i]>1 then w:=w*Product([1..mu[i]-1]+s, j->(j,j+1)); fi;
    s:=s+mu[i];
  od;
  return w;
end;

# Return the set of generators for the Young subgroup Sym(mu) given
# a composition <mu>
GeneratorsYoungSubgroupComposition:=function(mu) local J,s,m;
  J:=[];
  s:=0;
  for m in mu do
    Append(J, [s+1..s+m-1]);
    s:=s+m;
  od;
  return J;
end;

# return the quantum/Guassian integer [n]_q
QuantumInteger:=function(q,n) local i;
  if n=0 then return 0*q^0;
  elif n>0 then return Sum([0..n-1], i->q^i);
  else return -q^n*QuantumInteger(q,-n);
  fi;
end;

# return quantum/Guassian factorial [n]_q!
QuantumFactorial:=function(q,n) local i;
  return Product([2..n], i->QuantumInteger(q,i));
end;

# return the index polynomial [n]!/[mu]!
IndexPolynomialComposition:=function(q,mu) local i;
  return QuantumFactorial(q,Sum(mu))/Product(mu,i->QuantumFactorial(q,i));
end;

#M Calculates the determinant of the Gram matrix of the Specht module
## (defined over the integers) 
IntegralGramDeterminant:=function(arg)
    local  mu, mud, g, hooklen, c, row, r, s, v;
    mu := Flat( arg{[ 1 .. Length( arg ) ]} );
    Sort( mu );
    mu := mu{[ Length( mu ), Length( mu ) - 1 .. 1 ]};
    mud := ConjugatePartition( mu );
    hooklen := [  ];
    for r  in [ 1 .. Length( mu ) ]  do
      hooklen[r] := [  ];
      for c  in [ 1 .. mu[r] ]  do
        hooklen[r][c] := mu[r] + mud[c] - r - c + 1;
      od;
    od;
    g := 1;
    for c  in [ 1 .. mu[1] ]  do
      for row  in [ 1 .. mud[1] ]  do
        for r  in [ row + 1 .. mud[1] ]  do
          if mu[row] >= c and mu[r] >= c  then
            v := hooklen[row][c] / hooklen[r][c];
            if v <> 1  then
              s := AddRimHook(RemoveRimHook(mu,r,c,mud), row, hooklen[r][c]);
              if s <> false  then
                g:=g*v^((-1)^(s[2]+mud[c]-r)*SpechtDimension(s[1]));
              fi;
            fi;
          fi;
        od;
      od;
    od;
    return g;
end;

## #F Calculates the Specht modules in sum_{i>0}S^lambda(i) using the
## ## q-analogue of Schaper's theorem. 
## ## Uses H.valuation.
## ##   Usage:  Schaper(H,mu);
## GramDeterminant:=function(arg)
##   local H, mu, mud, top, bot, hooklen, c, row, r, s, v;
##   
##   if arg=[] or not ( IsRec(arg[1]) and IsBound(arg[1].IsSpecht) ) then
##     Error("usage, GramDeterminant(<H>,<mu>)");
##   fi;
##   H:=arg[1];
##   mu:=Flat(arg{[2..Length(arg)]});
## 
##   Sort(mu); mu:=mu{[Length(mu),Length(mu)-1..1]};
##   mud:=ConjugatePartition(mu);
##   hooklen:=[];
##   for r in [1..Length(mu)] do
##     hooklen[r]:=[];
##     for c in [1..mu[r]] do
##       hooklen[r][c]:=mu[r] + mud[c] - r - c + 1;
##     od;
##   od;
## 
##   top:=[1..Sum(mu)]*0;
##   bot:=Copy(top);
##   for c in [1..mu[1]] do
##     for row in [1..mud[1]] do
##       for r in [row+1..mud[1]] do
##         if mu[row] >=c and mu[r] >=c then
##           v:=hooklen[row][c] - hooklen[r][c];
##           if v<>0 then
##             s:=AddRimHook(RemoveRimHook(mu,r,c,mud),row,hooklen[r][c]);
##             if s<>false then
##               s[1]:=SpechtDimension(s[1]);
##               if (s[2]+mud[c]-r) mod 2 = 0 then
##                 top[hooklen[row][c]]:=top[hooklen[row][c]]+s[1];
##                 bot[hooklen[r][c]]:=bot[hooklen[r][c]]+s[1];
##               else 
##                 bot[hooklen[row][c]]:=bot[hooklen[row][c]]+s[1];
##                 top[hooklen[r][c]]:=top[hooklen[r][c]]+s[1];
##               fi;
##             fi;
##           fi;
##         fi;
##       od;
##     od;
##   od;
##   return schaper;
## end;  # Schaper()

#U MullineuxMap(e|H|d, mu)                                                
#U MullineuxMap(x)                                                        
#M This function returns the image of <mu> under the Mullineux map using  
#M the Kleshcehev(-James) algorihm, or the supplied decomposition matrix. 
#M Alternatively, given a "module" x it works out the image of x under    
#M Mullineux.                                                             
MullineuxMap:=function(arg) local e, mu, x, v, module; 
  if Length(arg)=1 and IsSpecht(arg[1]) then   ## MullineuxMap(x)
    x:=arg[1]; 
    if x=false or not IsERegular(x.H.e,x.parts[Length(x.parts)]) then   
      Print("# The Mullineux map is defined only for e-regular partitions\n");
      return false;
    fi;
    module:=x.module;
    if x=false or x=0*x then return false; fi;
    if x.module="S" then
      return x.H.operations.New("S",x.coeffs,List(x.parts,ConjugatePartition));
    elif x.module="Sq" then
      v:=x.H.info.Indeterminate;
      return x.H.operations.New("Sq",
            List([1..Length(x.coeffs)],
                mu->Value(v^-EWeight(x.H.e,x.parts[mu])*x.coeffs[mu],v^-1)),
            List(x.parts,ConjugatePartition) );
    elif Length(x.module)=1 then  # either P or D
      return x.H.operations.New(x.module,x.coeffs[mu],
                     List(x.parts,mu->MullineuxMap(x.H.e,mu)));
    else                          # either Pq or Dq
      v:=x.H.info.Indeterminate;
      return x.H.operations.New(x.module,
            List([1..Length(x.coeffs)],
                mu->Value(v^-EWeight(x.H.e,x.parts[mu])*x.coeffs[mu],v^-1)),
                     List(x.parts,mu->MullineuxMap(x.H.e,mu)));
    fi;
  elif Length(arg)>1 then
    e:=arg[1];
    mu:=Flat(arg{[2..Length(arg)]});
    if IsDecompositionMatrix(e) then          ## MullineuxMap(d,mu)
      x:=e.H.P(e,mu);
      if x=false or x=0*x then
        Print("MullineuxMap(<d>,<mu>), P(<d>,<mu>) not known\n");
        return false;
      else return ConjugatePartition(x.parts[1]);
      fi;
    fi;
    if IsRec(e) then
      if IsBound(e.e) then e:=e.e; elif IsBound(e.H) then e:=e.H.e; fi;
    fi; 
    if IsInt(e) then
      if not IsERegular(e,mu) then                     ## q-Schur algebra
        Error("# The Mullineux map is defined only for e-regular ",
              "partitions\n");
      else return PartitionGoodNodeSequence(e,
                    List(GoodNodeSequence(e,mu),x->-x mod e));
      fi;
    fi;
  fi;
  Error("usage: MullineuxMap(e|H|d, mu) or MullineuxMap(x)\n");
end;

#U Schaper(H,mu)                                                          
#M Calculates the Specht modules in sum_{i>0}S^lambda(i) using the        
#M q-analogue of Schaper's theorem. Uses H.valuation.                     
Schaper:=function(arg)
  local H, mu, mud, schaper, hooklen, c, row, r, s, v;
  if arg=[] or not ( IsRec(arg[1]) and IsBound(arg[1].valuation) ) then
    Error("usage, Schaper(<H>,<mu>)");
  fi;
  H:=arg[1];
  mu:=Flat(arg{[2..Length(arg)]});
  Sort(mu); mu:=mu{[Length(mu),Length(mu)-1..1]};
  mud:=ConjugatePartition(mu);
  hooklen:=[];
  for r in [1..Length(mu)] do
    hooklen[r]:=[];
    for c in [1..mu[r]] do
      hooklen[r][c]:=mu[r] + mud[c] - r - c + 1;
    od;
  od;
  schaper:=H.operations.New("S",[0],[[]]);
  for c in [1..mu[1]] do
    for row in [1..mud[1]] do
      for r in [row+1..mud[1]] do
        if mu[row] >=c and mu[r] >=c then
          v:=H.valuation(hooklen[row][c]) 
                - H.valuation(hooklen[r][c]);
          if v<>0 then
            s:=AddRimHook(RemoveRimHook(mu,r,c,mud),row,hooklen[r][c]);
            if s<>false then
              schaper:=schaper+H.operations.New("S",[(-1)^(s[2]+mud[c]-r)*v],[s[1]]);
            fi;
          fi;
        fi;
      od;
    od;
  od;
  return schaper;
end;  # Schaper()

#U SchaperMatrix(d)                                                       
#M Returns the matrix of upper bounds on the entries in the decomposition 
#M matrix <d> given by the Jantzen sum formula.                           
SchaperMatrix:=function(d) local r, C, c, coeff, sh, shmat;
  shmat:=d.H.operations.NewDecompositionMatrix(d.rows,d.cols,true);
  shmat.operations:=Copy(shmat.operations);
  Unbind(shmat.operations.Induce);
  shmat.d:=List(shmat.cols, c->rec(parts:=[],coeffs:=[]));
  C:=Length(d.cols)+1; ## this keeps track of which column we're up to
  for r in [Length(d.rows),Length(d.rows)-1..1] do
    if d.rows[r] in d.cols then C:=C-1; fi;
    sh:=Schaper(d.H,d.rows[r]);
    for c in [C..Length(d.cols)] do
      coeff:=InnerProduct(sh,d.P(d,d.cols[c]));
      if coeff<>false and coeff<>0*coeff then 
        Add(shmat.d[c].parts,r);
        Add(shmat.d[c].coeffs,coeff);
      fi;
    od;
  od;
  sh:=[];
  for c in [1..Length(d.d)] do
    Add(shmat.d[c].parts, Position(shmat.rows,shmat.cols[c]));
    Add(shmat.d[c].coeffs,1);
  od;
  shmat.matname:="Schaper matrix";
  return shmat;
end;

#U JantzenCoefficients(e,p,partitions)...................................
#M Given a set of partitions <partitions> return the matrix (J_{mu,nu}), 
#M where mu and nu run over <partitions> and J_{mu,nu} is the            
#M corresponding Jantzen coefficient.                                    
JantzenCoefficients:=function(e,p,partitions)
  local J, S, Jmu, pnu, mu, nu;
  J:=List(partitions,mu->[1..Length(partitions)]*0);
  if p=0 then S:=Schur(e);
  else S:=Schur(e,p);
  fi;
  LabelPartition:=mu->LabelByEQuotient(e,mu);
  for mu in [1..Length(partitions)] do
    Jmu:=Schaper(S,partitions[mu]);
    for nu in [1..Length(Jmu.parts)] do
      pnu:=Position(partitions,Jmu.parts[nu]);
      if IsInt(pnu) then J[mu][pnu]:=Jmu.coeffs[nu]; fi;
    od;
  od;
  return DecompositionMatrixMatrix(S,J,partitions,partitions);
end;

#U PPowerGramDeterminant(p,mu)...........................................
#M Returns the power of <p> wihch divided the determinant of the Gram    
#M matrix of the Specht module S<mu>.                                    
PPowerGramDeterminant:=function(p, mu)
    local  mu, mud, gp, hooklen, c, row, r, s, v;
    mud := ConjugatePartition( mu );
    hooklen := [  ];
    for r  in [ 1 .. Length( mu ) ]  do
      hooklen[r] := [  ];
      for c  in [ 1 .. mu[r] ]  do
        hooklen[r][c] := mu[r] + mud[c] - r - c + 1;
      od;
    od;
    gp := 0;
    for c  in [ 1 .. mu[1] ]  do
      for row  in [ 1 .. mud[1] ]  do
        for r  in [ row + 1 .. mud[1] ]  do
          if mu[row] >= c and mu[r] >= c  then
            v := PPower(p, hooklen[row][c]) - PPower(p, hooklen[r][c]);
            if v <> 0  then
              s := AddRimHook(RemoveRimHook(mu,r,c,mud), row, hooklen[r][c]);
              if s <> false  then
                gp:=gp+v*((-1)^(s[2]+mud[c]-r)*SpechtDimension(s[1]));
              fi;
            fi;
          fi;
        od;
      od;
    od;
    return gp;
end;

#U FrobeniusSymbolPartition(mu)                                           
#U FrobeniusSymbolPartition(mu_1,mu_2,...)                                
#M Returns the Frobenius symbol of the partition <mu>.                    
FrobeniusSymbolPartition:=function(arg)
  local mu, muc, a, b, i;
  mu:=Flat(arg);
  if not ForAll(mu,IsInt) then Error("usage: FrobeniusSymbolPartition(mu)"); fi;
  muc:=ConjugatePartition(mu);
  a:=[];b:=[];
  i:=1;
  while i<=Length(mu) and i<=Length(muc) do
    if IsBound(mu[i]) and mu[i]>=i then Add(a,mu[i]-i+1);fi;
    if IsBound(muc[i]) and muc[i]>=i then Add(b,muc[i]-i+1);fi;
    i:=i+1;
  od;
  return [a,b];
end;

#U IsCarterPartition(H,mu)                                                
#M Returns <true> if the partition <mu> satisfies the Carter criterion    
IsCarterPartition:=function(H,mu)
  local hooks, row, i,j;
  # compute the H-valuation of the hook lengths in [mu']
  hooks:=List(HookLengthsPartition(ConjugatePartition(mu)), 
	   row->List(row, H.valuation));
  return ForAll([1..Length(hooks)], 
           i->ForAll([2..Length(hooks[i])], j->hooks[i][j]=hooks[i][1]));
end;

#U HookPartitions(n)                                                      
#M Return a list of the hook partitions of n.                             
HookPartitions:=function(n) local l,p,hooks;
  hooks:=[];
  for l in [1..n-1] do
    p:=[0..l]*0+1; 
    p[1]:=n-l;
    Add(hooks, p);
  od;
  return hooks;
end;

#U ScopesClassPartition(H,mu)                                             
#M Return a list containing the partitions which are in the same Scopes   
#M class as the partition <mu>.                                           
ScopesClassPartition:=function(H,mu)
  local HigherInScopesClass, LowerInScopesClass, scopes;
  # return the list of all bigger partitions in the same Scopes class
  # as mu
  HigherInScopesClass:=function(H,mu)
    local iaddables, iremovables, scopes, nu, i, r;
    iaddables:=AddableINodes(H,mu); 
    iremovables:=RemovableINodes(H,mu); 
    scopes:=[];
    for i in [1..H.e] do
      if iaddables[i]<>[] and iremovables[i]=[] then
	nu:=Copy(mu);
	for r in iaddables[i] do
	  if r>Length(nu) then Add(nu,1);
	  else nu[r]:=nu[r]+1;
	  fi;
	od;
	Add(scopes, nu);
	Append(scopes, HigherInScopesClass(H,nu));
      fi;
    od;
    return scopes;
  end;
  # return the list of all smaller partitions in the same Scopes class
  # as mu
  LowerInScopesClass:=function(H,mu)
    local iaddables, iremovables, scopes, nu, i, r;
    iaddables:=AddableINodes(H,mu); 
    iremovables:=RemovableINodes(H,mu); 
    scopes:=[];
    for i in [1..H.e] do
      if iremovables[i]<>[] and iaddables[i]=[] then
	nu:=Copy(mu);
	for r in iremovables[i] do
	  nu[r]:=nu[r]-1;
	  if nu[r]=0 then Unbind(nu[r]);fi;
	od;
	Add(scopes, nu);
	Append(scopes, LowerInScopesClass(H,nu));
      fi;
    od;
    return scopes;
  end;
  if IsInRouquierBlock(H,mu) then
    Print("This is a Rouquier block\n");
    return [mu];
  fi;
  scopes:=[mu];
  Append(scopes, HigherInScopesClass(H,mu));
  Append(scopes, LowerInScopesClass(H,mu));
  return Set(scopes);
end;

##U ScopesReducedPartition(H|e,mu).........................................
##M Given a partition <mu> return the minimal partition obtained by        
##M applying Scopes moves to the <e>-abacus of <mu>; that is, swapping     
##M adjacent runners i+1 and i whenever runner i+1 has w or more beads     
##M than runner i.                                                         
#ScopesReducedPartition:=function(arg)
#  local e, mu, w, abacus, j, ai, i;
#  e:=arg[1];
#  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
#  mu:=Flat(arg{[2..Length(arg)]});
#  if not (IsInt(e) and IsPartition(mu)) then
#    Error("usage: ScopesScopesPartition(H|e,mju)");
#  fi;
#  w:=EWeight(e,mu);
#  abacus:=EAbacusRunners(e,mu, Length(mu));
#  for i in [2..e] do
#    j:=i;
#    while j>1 and Length(abacus[j-1])<=Length(abacus[i])-w do
#      j:=j-1;
#    od;
#    if j<i then
#      ai:=abacus[i];
#      abacus[i]:=abacus[j];
#      abacus[j]:=ai;
#      Print("i = ", i, ", j = ", j, ", abacus = ", abacus,"\n");
#      Print("  --->  ", PartitionAbacusRunners(abacus),"\n");
#    fi;
#  od;
#  return PartitionAbacusRunners(abacus);
#end;

#U SwapRunnersDown(e,mu).................................................
#M Given partition <mu> return the list of *smaller* partitions which are
#M obtained by swapping adjacent runners in the abacus of <mu>.          
SwapRunnersDown:=function(e,mu) local swaps, aba, res, abai, nu, i;
  swaps:=[];
  aba:=FullEAbacusRunners(e,mu);
  # find the residue of each runner by finding the highest beta number
  res:=e;
  while res>1 and Length(aba[res])<Length(aba[res-1]) do
    res:=res-1;
  od;
  # the residue of runner i is now i-res+mu[1]-1.
  for i in [1..e-1] do
    if Length(aba[i+1])>Length(aba[i]) then
      abai:=Copy(aba);
      abai[i]:=aba[i+1];
      abai[i+1]:=aba[i];
      nu:=PartitionAbacusRunners(abai);
      AddSet(swaps, [nu,i, (i-res+mu[1]) mod e]); 
    fi;
  od;
  return swaps;
end;

#U PartitionsInSameBlock(H,mu)                                            
#M Return a SET of the partitions which are in the same block as <mu>.    
PartitionsInSameBlock:=function(e,mu) local core, q;
  if IsRec(e) then e:=e.e; fi;
  core:=ECore(e,mu);
  return Set( List(PartitionTuples((Sum(mu)-Sum(core))/e, e), 
                     q->CombineEQuotientECore(q,core)) );
end;

#U CoxeterWordSymmPerm(x)                                                 
#M Given a permutation <x> in Sym(n) return the corresponding CoxeterWord 
#M (we assume that <x> does indeed belong to W=Sym(n) (certainly it       
#M belongs M in some symmetric group).                                    
CoxeterWordSymmPerm:=function(x) local w, i;
  w:=[];
  while x<>() do
    i:=1; while i^x<(i+1)^x do i:=i+1; od;
    Add(w,i);
    x:=(i,i+1)*x;
  od;
  return w;

end;

#U SymmWordSymmPerm(n,w)                                                  
#M Given a permutaiton <w> in Sym(<n>) return the word [w_1,...,w_n] for  
#M <w> which is given by i^w=w_i.                                         
SymmWordSymmPerm:=function(n,w) local i;
  return List([1..n], i->i^w);
end;

#U SymmPermSymmWord(word)                                                 
#M Given a list <word>=[w_1,...,w_n] return the corresponding             
#M permutation w which is given by i^w=w_i.                               
SymmPermSymmWord:=function(word)
  return MappingPermListList([1..Length(word)], word);
end;

##########################################################################

#C Formal characters of (Specht) modules, in the sense of Grojnowski      

#M The following internal record constructs a "module" in with which we   
#M can compute sums of formal JM characters.                              
SPECHT.JMCharacter:=ModuleRecord("JMFormalCharacters",
	                   ch->SPrint("e(",TightStringList(ch),")"),
			   Integers
);

#U JMCharacterTableau(H,tab)..............................................
#M The JM-character of a tableaux module is the linear combination of     
#M terms of the form e(i_1..i_n), where m_t L_j = i_j m_t + higher terms. 
JMCharacterTableau:=function(H,tab)
  local ch, r, c;
  ch:=[];
  for r in [1..Length(tab)] do
    for c in [1..Length(tab[r])] do
      ch[tab[r][c]]:=(c-r) mod H.e;
    od;
  od;
  return MakeModuleElement(SPECHT.JMCharacter, ch);
end;

#U JMCharacterSpechtModule(H,mu)..........................................
#M Return the formal character of the Specht module S(mu).                
JMCharacterSpechtModule:=function(H,mu) local chs,t;
   return Sum(StandardTableaux(mu),t->JMCharacterTableau(H,t));
end;

#U JMcharacter(x)..........................................................
#M Compute the JM character of a linear combination of SPECHT modules,
#M but first converting them to a linear combination of Specht modules.
JMCharacter:=function(a) local H, b, i;
  H:=a.H;
  b:=H.S(a);
  return Sum([1..Length(b.parts)], 
               i->b.coeffs[i]*JMCharacterSpechtModule(H,b.parts[i]));
end;

#U ParityPartition(e,mu)..................................................
#M Return the e-parity of the partition mu; that is, the sum of the leg   
#M lengths of mu modulo 2.                                                
ParityPartition:=function(e,mu) local runners, r, b, newRunners, lleg, s;
  if IsRec(e) and IsBound(e.e) then e:=e.e; fi;
  if not (IsInt(e) and IsPartition(mu)) then
    Error("usage, ParityPartition(e,mu)");
  fi;
  if IsECore(e,mu) then return 0; fi;
  runners:=EAbacusRunners(e,mu);
  r:=First([1..e], r->runners[r]<>[] and runners[r][1]>Length(runners[r])-1);
  b:=First([1..Length(runners[r])], b->not runners[r][b]-1 in runners[r]);
  newRunners:=Copy(runners);
  newRunners[r][b]:=newRunners[r][b]-1;
  lleg:=Length(Filtered([1..r-1],s->runners[r][b] in runners[s]))
        +Length(Filtered([r+1..e],s->runners[r][b]-1 in runners[s]));
  return Mod(ParityPartition(e,PartitionAbacusRunners(newRunners))+lleg,2);
end;

#U IsContainedInPartitions(mu,nu).........................................
#M Return true if the partition <mu> is contained in the partition <nu>,  
#M or rather [mu]\subseteq [nu].                                          
IsContainedInPartitions:=function(mu,nu) local i;
  if Length(mu)>Length(nu) then return false;
  else return ForAll([1..Length(mu)],i->mu[i]<=nu[i]);
  fi;
end;

DonkinNormalFormPartitionOps:=OperationsRecord("Partition normal form");
DonkinNormalFormPartitionOps.String:=function(self)
  return SPrint(self.s-1,"*(",TightStringList([self.l-1,self.l-2..0]),
              ") + ",self.d,"*(1^",self.l,") + ",self.s,
              "*(",TightStringList(self.chi),")");
end;
DonkinNormalFormPartitionOps.Print:=function(self)Print(String(self));end;
DonkinNormalFormPartitionOps.TeX:=String;

#U DonkinNormalFormPartition:=funtion(e,p,l,mu)...........................
#M Given a partition <mu> or length at least <l> return the normal form of
#M <mu>. That is, write the partition in the form                         
#M     mu = (s(mu)-1) rho_l + d(mu) omega_l + s chi(mu)                   
#M where s(mu)=max{s in {1,e,ep,ep^2,...} | mu_i-mu_{i+1}\equiv -1\mod s  
#M                                                for 1<= i<l            }
#M   d(mu) = mu[l] mod sm and chi(mu) is the unique partition satisfying  
#M equation (*).                                                          
DonkinNormalFormPartition:=function(e,p,l,mu) local nf, s, m, i;
  if Length(mu)>l or (p>0 and Gcd(e,p)>1) then 
    Error("DonkinNormalFormPartition(e,p,l,mu): len(mu) mnust be at least l!");
  fi;
  nf:=rec(e:=e, p:=p, l:=l, operations:=DonkinNormalFormPartitionOps);
  if Length(mu)<l-1 then 
    nf.s:=1; 
    nf.d:=0;
    nf.chi:=Copy(mu);
  else
    nf.s:=1; m:=e;
    while m>0 and ForAll([1..l-1], i->(mu[i]-mu[i+1]+1) mod(nf.s*m)=0) do
      nf.s:=nf.s*m;
      m:=p;
    od;
    if not IsBound(mu[l]) then nf.d:=0;
    else nf.d:=mu[l] mod nf.s;
    fi;
    i:=1;
    nf.chi:=[];
    for i in [1..l] do
      m:=(nf.s-1)*(l-i)+nf.d;
      if m<mu[i] then nf.chi[i]:=(mu[i]-m)/nf.s; fi;
    od;
  fi;
  return nf;
end;
    
