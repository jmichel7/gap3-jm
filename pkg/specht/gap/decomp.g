#########################################################################
##  SPECHT - decomp.g : routines for working with SPECHT decomposition ##
##           matrices                                                  ##
##                                                                     ##
##     These programs, and the enclosed libraries, are distributed     ##
##     under the usual licensing agreements and conditions of GAP.     ##
##                                                                     ##
##     Andrew Mathas                                                   ##
##                                                                     ##
#########################################################################

#D Functions for manipulating decomposition matrices.

#C These functions were previously contained in specht.g. The idea was   
#C to move all functions associated with decomposition matrices into the 
#C one place, but due to the cumbersome nature of the Specht() records   
#C this was not quite possible. New functions for version 3.1 include:   
#C  o BlockDecompositionMatrix                                           
#C  o DecompositionNumberMultiplicityOne                                 
#C  o DecompositionNumberRouquierBlock                                   
#C  o DominatingDecompositionMatrix                                      
#C  o EmptyDecompositionMatrix                                           
#C  o FormalCharactersDecompositionMatrix                                
#C  o IsInRouquierBlock                                                  
#C  o LinkageClassesDecompositionMatrix                                  
#C  o MultiplicityOneDecomposition                                       
#C  o RestrictedBlockDecompositionMatrix                                 
#C  o RouquierBlock                                                      
#C  o RowLimitedDecompositionMatrix                                      
#C  o SubmatrixDecompositionMatrix                                       
#C  o WeightDecompositionMatrix                                          
#C  o WeightTwoDecompositionMatrix                                       

#C yup; this, and its associates, could do with shorter names...
IsDecompositionMatrix:=function(d)
  return IsRec(d) and IsBound(d.IsDecompositionMatrix);
end;

## DecompositionMatrices
## Decomposition matrices 'd' in Specht are represented as records with the
## following components:
##   
##   d : a list, indexed by d.cols, where each entry is a record
##       corresponding to a column of 'd'; this record has components
##       two * sets* coeffs and parts, where parts is the index of the 
##       corresponding partition in d.rows.
##   rows : the *set* of the partitions which make up the rows or 'd'.
##   cols : the *set* of the partitions which make up the rows or 'd'.
##   inverse : a list of records containing the inverse of 'd'. These
##          records are computed only as needed.
##   dimensions : a list of the dimenions of the simle modules; again
##          comuted only as needed.
##   IsDecompositionMatrix : false if 'd' is a crystallized decomposition
##          matrix, and true otherwise.
##   H :a pointer back to the corresponding algebra
##   operations :
##     = : equality.
##     Print, TeX, Matrix : printing, TeX, and a GAP Matrix.
##     AddIndecomposable : for adding a PIM to 'd'.
##     Store : for updating Specht's internal record of 'd'.
##     S, P, D: for accessing the entries of 'd' and using 'd' to
##              convert between the various types of 'H'--modules.
##              These are actually records, each containing three
##              functions S(), P(), and D(); so X.Y() tells 'd' how
##              to write an X-module as a linear comination of Y-modules.
##     Invert : calculates D(mu) using 'd'.
##     IsNewIndecomposable : the heart of the 'IsNewIndecomposable'
##              function.
##     Induce : for inducing decomposition matrices (non--crystallized).
##   P : a short-hand for d.H.P('d',<mu>).

#U MultiplictyOneSimpleDecomposition(H,mu)
#  Assuming that <x> is multiplicity free rewrite it as a sum of
#  simple modules by using Janzten-Schaper. These assumptions are met
#  when x=H.S(a,b) or when x=H.S(mu) where mu is a partition of weight
#  2 and e>2.
MultiplicityOneDecomposition:=function(H, mu)
  local simples, schaper, nu;
  if not H.IsSpecht or IsERegular(H,mu) then simples:=H.D(mu);
  else simples:=0*H.D();
  fi;
  schaper:=Schaper(H,mu);
  if schaper<>0*schaper then
    for nu in [1..Length(schaper.parts)] do
      simples:=simples
	  +schaper.coeffs[nu]*MultiplicityOneDecomposition(H,schaper.parts[nu]);
    od;
  fi;
  simples.coeffs:=simples.coeffs*0+1; # change any n's to 1 
  return simples;
end;

# compute d_{mu,nu} using MultiplicityOneDecomposition().
DecompositionNumberMultiplicityOne:=function(H,mu,nu)
  return Coefficient(MultiplicityOneDecomposition(H,mu),nu);
end;

#U IsInRouquierBlock(e,mu)
#M return true if mu belongs to a Rouquier block. That is, for any
#M abacus configuration of <mu>
#M e-weight w
IsInRouquierBlock:=function(e,mu)
  return Set(Flat(PyramidPartition(e,mu).a))=[0];
  #w:=EWeight(e,mu);
  #abacuslens:=List(FullEAbacusRunners(e,mu), Length);
  #Sort(abacuslens);
  #return ForAll([1..e-1],r->AbsInt(abacuslens[r]-abacuslens[r+1])>=w-1);
end;

# compute the decomposition number [S(mu):D(nu)] when both partitions
# below to a Rouquier block.
DecompositionNumberRouquierBlock:=function(H,mu,nu)
  local RouquierMultiplicity, core, w, muaba, nuaba, b, quotients, i, d;
  # recursively compute 
  #   sum prod a^{q^i}{alpha,beta^i}a^{q^{i+1}}_{beta^i,beta^{i+1}}....
  RouquierMultiplicity:=function(quotients, a, i) local b;
    if i=Length(quotients) then 
      return Length(InverseLittlewoodRichardsonRule(quotients[i],a));
    elif Mod(i,2)=0 then
      return Sum(Collected(InverseLittlewoodRichardsonRule(quotients[i],a)), 
              b->b[2]*RouquierMultiplicity(quotients, ConjugatePartition(b[1][2]),i+1));
    else 
      return Sum(Collected(InverseLittlewoodRichardsonRule(quotients[i],a)), 
              b->b[2]*RouquierMultiplicity(quotients, b[1][2],i+1));
    fi;
  end;
  if ECore(H,mu)<>ECore(H,nu) then return 0; fi;
  if not IsInRouquierBlock(H,mu)  then return false; fi;
  core:=ECore(H,mu);                      # the e-core of mu and nu
  if core=mu or mu=nu then return 1; 
  elif not Dominates(nu,mu) then return 0;
  fi;
  if H.IsSpecht and not IsERegular(H,nu) then return false; fi;
  w:=(Sum(mu)-Sum(core))/H.e;             # the e-weight of mu and nu
  if H.p>0 and w>=H.e then return false; fi;
  # we want an abacus with enough beads for every partition in the block.
  muaba:=EAbacusRunners(H.e,mu,Length(core)+w*H.e);  
  nuaba:=EAbacusRunners(H.e,nu,Length(core)+w*H.e);  
  # Reorder the runners using Scopes equivalence (i.e. order them in
  # in terms of increasing number of beads). This is legal as, by 
  # definition, the bead gaps are at least the weight for a Rouquier core.
  SortParallel(nuaba, muaba, function(a,b) return Length(a)<Length(b);end);
  # now find the quotients
  quotients:=[];
  for i in [1..H.e] do
    quotients[2*i-1]:=PartitionBetaSet(Set(nuaba[i]));
    quotients[2*i]:=PartitionBetaSet(Set(muaba[i]));
  od;
  Add(quotients,[0]);
  # finally, we compute the decomposition number, which is the sum of
  # all (alpha^i,beta^i), i=1..e, of the Littlewood-Richardson
  # coefficients c^{mu^i}_{alpha^i,beta^i}*c^{nu^i}_{beta^i,alpha^{i+1}}
  return RouquierMultiplicity(quotients,quotients[1],2);
end;

#C CrystalizedDecompositionNumberRouquierBlock(H,mu,nu)...................
#C CrystalizedDecompositionNumberRouquierBlock:=function(H,mu,nu)

#U RouquierBlock(H,z,w)
#M Return the decomposition matrix of the Rouquier block for <H>
#M which has weight <w> and core rho(<z>)=RouquierCore({H},{Z}
RouquierBlock:=function(H,z,w) local rhoabacus, rho, i;
  rhoabacus:=[];
  for i in [1..H.e] do
    rhoabacus[i]:=[1..(i-1)*(z-1)]-1;
  od;
  rho:=PartitionAbacusRunners(rhoabacus);
  if rho=[] then rho:=[w*H.e];
  else rho[1]:=rho[1]+w*H.e;
  fi;
  if H.p=0 or w<H.p then
    return BlockDecompositionMatrix(H,rho,DecompositionNumberRouquierBlock);
  else return BlockDecompositionMatrix(H,rho);
  fi;
end;


#U DecompositionNumber(H,mu,nu)                                          
#U DecompositionNumber(d,mu,nu)                                          
#M Returns the decomposition number d_{mu,nu}; row and column removal are
#M used if the projective P(nu) is not already known. If unable to       
#M calculate the decomposition number we return false. Note that         
#M H.IsSpecht is false if we are looking at decomposition matrices of a  
#M q-Schur algebra and true for a Hecke algebra.                         
DecompositionNumber:=function(H,mu,nu) 
  local Pnu, w, RowAndColumnRemoval, a;
  SpechtInfo("# DecompositionNumber: mu=", mu, ", nu=", nu,"\n");
  if IsDecompositionMatrix(H) then
    Pnu:=H.P(H,nu);
    if Pnu<>false then 
      SpechtInfo("# DecompositionNumber: known PIM - , mu=", mu, ", nu=", nu,"\n"); 
      return Coefficient(Pnu,mu); 
    fi;
    H:=H.H;
  elif mu=nu then 
    SpechtInfo("# DecompositionNumber: equality - , mu=", mu, ", nu=", nu,"\n"); 
    return 1;
  elif not Dominates(nu,mu) then 
    SpechtInfo("# DecompositionNumber: dominance - , mu=", mu, ", nu=", nu,"\n"); 
    SPECHT.Info:="Dominance";
    return 0;
  elif IsRec(H) and IsBound(H.IsSpecht) then
    Pnu:=H.operations.P.S(H.operations.New("P",[1],[nu]),true);  
    if Pnu<>false then 
      SpechtInfo("# DecompositionNumber: known PIM - , mu=", mu, ", nu=", nu,"\n"); 
      return Coefficient(Pnu,mu); 
    fi;
  else Error("usage, DecompositionMatrix(<d> or <H>, <mu>, <nu>)");
  fi;
  if H.IsSpecht and not IsERegular(H.e, nu) then
    return false;
    #Error("DecompositionNumber(H,mu,nu), <nu> is not ",H.e,"-regular");
  fi;

  # check that mu and nu are in the same block
  if ECore(H.e,mu)<>ECore(H.e,nu) then return 0; fi;

  w:=EWeight(H,nu);
  if (H.e=2 and w=2) then
    SpechtInfo("# Weight two: mu=", mu, ", nu=", nu,"\n"); 
    return DecompositionNumberWeightTwoPartitions(H, mu,nu);
  elif (H.e>3 and w=3) or Length(mu)=2 then
    SpechtInfo("# Multiplicity free Specht module: mu=", mu, ", nu=", nu,"\n"); 
    SPECHT.Info:="Multiplicity free decomposition number";
    return Coefficient(MultiplicityOneDecomposition(H,mu),nu);
  fi;

  # if mu and nu belong to a Rouquier block then we just compute 
  if w<H.p and IsInRouquierBlock(H.e,mu) then
    SpechtInfo("# RouquierBlock: mu = ", mu, ", nu=", nu,"\n"); 
    return DecompositionNumberRouquierBlock(H,mu,nu);
  fi;

  ## Next we try row and column removal (James, Theorem 6.18)
  ## (here fn is either the identity or conjugation).
  RowAndColumnRemoval:=function(fn) local m,n,i,d1,d2;
    ## H, mu, and nu as above
    mu:=fn(Copy(mu)); nu:=fn(Copy(nu));
    SpechtInfo("# Row and column removal: mu=", mu, ", nu=", nu,"\n"); 
    m:=0; n:=0; i:=1; 
    while i<Length(nu) and i<Length(mu) do
      m:=m+mu[i]; n:=n+nu[i];
      if m=n then 
        d2:=DecompositionNumber(H, fn(mu{[i+1..Length(mu)]}),
                   fn(nu{[i+1..Length(nu)]}));
        if d2=0 then 
          SPECHT.Info:="Row and Column removal";
	  return d2;
        elif IsInt(d2) then 
          d1:=DecompositionNumber(H, fn(mu{[1..i]}),fn(nu{[1..i]}));
          if IsInt(d1) then 
            SPECHT.Info:="Row and Column removal";
	    return d1*d2; 
	  fi;
        fi;

      fi;
      i:=i+1;
    od;
    return false;
  end;
  Pnu:=RowAndColumnRemoval(a->a);
  if Pnu=false then 
    Pnu:=RowAndColumnRemoval(ConjugatePartition);
    if IsInt(Pnu) then 
      SpechtInfo("# DecompositionNumber: column removal - , mu=", mu, ", nu=", nu,"\n"); 
    fi;
  else SpechtInfo("# DecompositionNumber: row removal - , mu=", mu, ", nu=", nu,"\n"); 
  fi;
  return Pnu;
end;

#C Returns a list of those e-regular partitions mu such that Px-P(mu)
#M has positive coefficients (ie. those partitions mu such that P(mu)
#M could potentially split off Px). Simple minded, but useful.
Obstructions:=function(d,Px) local obs, mu, Pmu, possibles;
  obs:=[];
  if d.H.IsSpecht then 
    possibles:=Filtered(Px.parts, mu->IsERegular(Px.H.e, mu));
  else possibles:=Px.parts;
  fi;
  for mu in possibles do
    if mu<>Px.parts[Length(Px.parts)] then
      Pmu:=d.P(d,mu);
      if Pmu=false or PositiveCoefficients(Px-Pmu) then Add(obs,mu); fi;
    fi;
  od;
  return obs{[Length(obs),Length(obs)-1..1]};
end;

## Interface to d.operations.IsNewDecompositionMatrix. Returns true
## if <Px> contains an indecomposable not listed in <d> and false
## otherwise. Note that the value of <Px> may well be changed by
## this function. If the argument <mu> is used then we assume
## that all of the decomposition numbers down given by <Px> down to 
## <mu> are correct. Note also that if d is the decomposition matrix
## for H(Sym_{r+1}) then the decomposition matrix for H(Sym_r) is passed
## to IsNewDecompositionMatrix.
##   Usage: IsNewIndecomposable(<d>,<Px> [,<mu>]);
## If <mu> is not supplied then we set mu:=true; this
## turns on the message printing in IsNewIndecomposable().
IsNewIndecomposable:=function(arg) local d, oldd, mu;
  if Length(arg)<2 or
  not (IsDecompositionMatrix(arg[1]) and IsSpecht(arg[2]) ) then
    Error("usage, IsNewIndecomposable(<d>,<Px> [,<mu>]");
  fi;
   
  d:=arg[1];
  if Length(arg)=2 then mu:=true;
  else mu:=arg{[3..Length(arg)]};
  fi;
  oldd:=d.H.operations.FindDecompositionMatrix(Sum(d.rows[1])-1);
  return d.operations.IsNewIndecomposable(d,arg[2],oldd,mu);
end;

#P A front end to d.operations.AddIndecomposable. This funciton adds <Px>
## into the decomposition matrix <d> and checks that it is compatible with
## its image under the Mullineux map, if this is already in <d>, and
## inserts it if it is not.
AddIndecomposable:=function(d, Px)
  if IsSpecht(Px) and Px.module="S" and IsDecompositionMatrix(d) then
    if Position(d.cols, Px.parts[Length(Px.parts)])=false then 
      Print("# The projective P(",TightStringList(Px.parts[Length(Px.parts)]),
            ") is not listed in <D>\n");
    else d.operations.AddIndecomposable(d,Px,true);
    fi;
  else Error("usage: AddIndecomposable(<d>,<Px>)\n");
  fi;
end; # AddIndecomposable

#P Removes the columns for <Px> in <d>
RemoveIndecomposable:=function(arg) local d, r, c;
  d:=arg[1];
  c:=Position(d.cols, Flat(arg{[2..Length(arg)]}));
  if c=false then 
    Print("RemoveIndecomposable(<d>,<mu>), <mu> is not listed in <d>\n");
  else Unbind(d.d[c]);
  fi;
end;

#U MissingIndecomposables(d)                                             
#U MissingIndecomposables(d,silent)                                      
#M In the first form print the list of missing PIMs from <d>. In the
#M second form return the list of missing partitions.
MissingIndecomposables:=function(arg) local d, c, missing;
  d:=arg[1];
  missing:=Filtered([1..Length(d.cols)], c->not IsBound(d.d[c]) );
  if missing<>[] and Length(arg)=1 then
    Print("The following projectives are missing from <d>:\n  ");
    for c in Reversed(missing) do
      Print("  ", d.cols[c]); 
    od;
    Print("\n");
  fi;
  if Length(arg)>1 then return d.cols{missing}; fi;
end; 

## When no ordering is supplied then rows are ordered first by length and
## then lexicographically. The rows and columns may also be explicitly
## assigned.
## Usage:
##   DecompositionMatrix(H, n [,ordering]);
##   DecompositionMatrix(H, <file>) ** force Specht() to read <file>
SpechtOps.DecompositionMatrix:=function(arg) local H, d, Px, c, n;
  if arg=[] or not (IsRec(arg[1]) and IsBound(arg[1].IsSpecht)) 
  or Length(arg)=1 or Length(arg)>3 then
    Error("usage, DecompositionMatrix(<H>, <n>|<file> [,Ordering])");
  fi;
  H:=arg[1];
  
  if IsString(arg[2]) then
    if IsBound(arg[3]) and IsFunc(arg[3]) then H.Ordering:=arg[3]; fi;
    d:=H.operations.ReadDecompositionMatrix(H,arg[2],false);
    if d<>false and not IsBound(d.matname) then ## override and copy
      d.operations.Store(d);
      MissingIndecomposables(d);
    fi;
  elif not IsInt(arg[2]) then
    Error("usage, DecompositionMatrix(<H>, <n>|<file> [,Ordering])");
  else
    n:=arg[2];
    if IsBound(arg[3]) and IsFunc(arg[3]) then 
      if arg[3]<>LengthLexicographicallyByQuotient then H.Ordering:=arg[3]; 
      else H.Ordering:=function(mu,nu) 
	     return LengthLexicographicallyByQuotient(H.e,mu,nu);
           end;
      fi;
    fi;
    d:=H.operations.FindDecompositionMatrix(n);
    
    if d=false then 
      if H.p>0 and n>2*H.e then  ## no point even trying
        Print("# This decomposition matrix is not known; use ",
              "CalculateDecompositionMatrix()\n# or ",
              "InducedDecompositionMatrix() to calculate with this matrix.",
              "\n");
        return d;
      fi;
      if H.IsSpecht then c:=ERegularPartitions(H.e,n);
      else c:=Partitions(n);
      fi;
      d:=H.operations.NewDecompositionMatrix(Partitions(n),c,true);
    fi;
    if ForAny([1..Length(d.cols)],c->not IsBound(d.d[c])) then
      for c in [1..Length(d.cols)] do
        if not IsBound(d.d[c]) then
          Px:=H.operations.P.S(H.operations.New("P",[1],[d.cols[c]]),true);
          if Px<>false then d.operations.AddIndecomposable(d,Px,false); 
          else Print("# Projective indecomposable P(", 
                     TightStringList(d.cols[c]),") not known.\n");
          fi;
        fi;
      od;
      d.operations.Store(d);
    fi;
  fi;
  if d<>false then   ## can't risk corrupting the internal matrix lists
    d:=ShallowCopy(d);    
  fi;
  return d;
end;  ## DecompositionMatrix

#C Return the submatrix of <d> with entries in the specified rows 
## and columns
SubmatrixDecompositionMatrix:=function(arg)
  local d, rows, cols, name, b, p, cpos, rowpos, r, p, c; 
  d:=arg[1];
  if Length(arg)<3 or Length(arg)>4 or not IsDecompositionMatrix(d) then
    Error("usage: SubmatrixDecompositionMatix(<d>, <rows>, <cols>), or\n",
      "       SubmatrixDecompositionMatix(<d>, <rows>, <cols>, <name>)\n");
  fi;
  rows:=arg[2]; cols:=arg[3];
  b:=d.H.operations.NewDecompositionMatrix(rows, cols,
        d.IsDecompositionMatrix);
  rowpos:=List(d.rows, p->Position(b.rows,p));
  for c in [1..Length(b.cols)] do
    cpos:=Position(d.cols, b.cols[c]);
    if IsBound(d.d[cpos]) then
      b.d[c]:=rec(coeffs:=[], parts:=[]);
      for r in [1..Length(d.d[cpos].parts)] do
        if rowpos[d.d[cpos].parts[r]] <> false then
          Add(b.d[c].coeffs, d.d[cpos].coeffs[r]);
	  Add(b.d[c].parts, rowpos[d.d[cpos].parts[r]]);
	fi;
      od;
    fi;
  od;
  if Length(arg)=4 then
    if IsBound(d.matname) then b.matname:=Concatenation(d.matname,",",arg[4]);
    else b.matname:=arg[4];
    fi;
  fi;
  return b;
end;

#U DominatedDecompositionMatrix(d,mu)                                    
#U DominatedDecompositionMatrix(d,mu1,mu2,...)                           
#M Return the submatrix of <d> whose rows and columns are indexed by     
#M the partitions which are dominated by the partition <mu>.             
DominatedDecompositionMatrix:=function(arg) local d, mu, p;
  d:=arg[1]; mu:=Flat(arg{[2..Length(arg)]});
  if not (IsDecompositionMatrix(d) and IsList(mu) ) then
    Error("usage: DominatedDecompositionMatrix(<d>, <mu>)"); 
  fi;
  return SubmatrixDecompositionMatrix(d,
           Filtered(d.rows, p->Dominates(mu,p)),
           Filtered(d.cols, p->Dominates(mu,p)),
	   Concatenation("Dominates:", LabelPartition(mu))
	 );
end;

#U EmptyDecompositionMatrix(H,rows,cols)..................................
#M Returns an empty <H>-decomposition matrix with the specified <rows> and
#M <cols>.                                                                
EmptyDecompositionMatrix:=function(H,rows,cols)
  return H.operations.NewDecompositionMatrix(rows,cols,true);
end;

#U EmptyCrystalizedDecompositionMatrix(H,rows,cols).......................
#M Returns an empty crystalized  <H>-decomposition matrix with the        
#M specified <rows> and <cols>.                                           
EmptyCrystalizedDecompositionMatrix:=function(H,rows,cols)
  return H.operations.NewDecompositionMatrix(rows,cols,false);
end;

# return the submatrix of the decomposition matrix with  rows and
# columns indexed by the partitions of <n> of weight 2. Here the
# function <DecompNo> could be either DecompositionNumber or
# DecompositionNumberWeightTwo, which is useful for comparison
# purposes...of course, both functions should give the same martix.
WeightTwoDecompositionMatrix:=function(H,n, DecompNo)
  local rows, cols, dmat, decnums, d, mu, nu;
  rows:=PartitionsOfGivenWeight(H.e, n,2);
  cols:=Filtered(rows, mu->IsERegular(H,mu));
  dmat:=EmptyDecompositionMatrix(H,rows,cols);
  for mu in [1..Length(cols)] do
    decnums:=rec(coeffs:=[], parts:=[]);
    for nu in [1..Length(rows)] do
      d:=DecompNo(H, rows[nu], cols[mu]);
      if d<>0 then # add into dmat
        Add(decnums.coeffs, d);
        Add(decnums.parts, nu);
      fi;
    od;
    Add(dmat.d, decnums);
  od;
  return dmat;
end;

#C Return the submatrix of <d> with entries indexed by partitions of
## weight <w>.
WeightDecompositionMatrix:=function(d, w) local d, mu;
  if not (IsDecompositionMatrix(d) and IsInt(w) ) then
    Error("usage: WeightDecompositionMatrix(<d>, <w>)"); 
  fi;
  if Filtered(d.cols, mu->EWeight(d.H,mu)=w)<>[] then
    return SubmatrixDecompositionMatrix(d,
           Filtered(d.rows, mu->EWeight(d.H,mu)=w),
           Filtered(d.cols, mu->EWeight(d.H,mu)=w), 
	   Concatenation("weight=", TightStringList(w))
	 );
  else return false;
  fi;
end;

#U BlockDecompostionMatrix(d,<mu>)........................................
#U BlockDecompositionMatrix(H,<mu>).......................................
#U BlockDecompositionMatrix(H,<mu>,<func>)................................
#M The first form returns the submtrix of <d> indexed by the partitions   
#M which have the same core as <mu>; the second form returns the          
#M decomposition matrix with rows and columns indexed by the partitions   
#M which are in the same block as <mu> and where the decomposition numbers
#M are computed using DecompositionNumber(), or the function <func> if    
#M supplied.                                                              
BlockDecompositionMatrix:=function(arg) 
  local d, mu, H, DecompNum, rows, cols, b, decnums, nu, p;
  if IsDecompositionMatrix(arg[1]) then
    d:=arg[1];
    mu:=ECore(d.H, Flat(arg{[2..Length(arg)]}) );
    if not (IsDecompositionMatrix(d) and IsList(mu) ) then
      Error("usage: BlockDecompositionMatrix(<d>, <mu>)"); 
    fi;
    return SubmatrixDecompositionMatrix(d,
             Filtered(d.rows, p->ECore(d.H,p)=mu),
             Filtered(d.cols, p->ECore(d.H,p)=mu), 
	     Concatenation("block=", TightStringList(mu))
	   );
  elif IsRec(arg[1]) and IsBound(arg[1].IsSpecht) and IsList(arg[2]) 
    and ( Length(arg)>1 and ForAll(arg[2],IsInt) )
    and ( Length(arg)=2 or  IsFunc(arg[3]) ) then
    H:=arg[1]; 
    mu:=arg[2]; 
    if Length(arg)=2 then DecompNum:=DecompositionNumber;
    else DecompNum:=arg[3];
    fi;
    rows:=PartitionsInSameBlock(H,mu);
    if H.IsSpecht then
      cols:=Filtered(rows,nu->IsERegular(H,nu));
    else cols:=rows;
    fi;
    b:=EmptyDecompositionMatrix(H,rows,cols);
    for mu in [1..Length(cols)] do
      decnums:=rec(coeffs:=[], parts:=[]);
      for nu in [1..Length(rows)] do
        d:=DecompNum(H, rows[nu], cols[mu]);
        if d<>0 then # add into b
          Add(decnums.coeffs, d);
          Add(decnums.parts, nu);
        fi;
      od;
      Add(b.d, decnums);
    od;
    return b;
  fi;
end;

#U BlockCrystalizedDecompositionMatrix(H,mu)..............................
#M Return the crystalized decomposition matrix for the <H>-block which    
#M contains the partitions <mu>.                                          
BlockCrystalizedDecompositionMatrix:=function(H,mu) local rows,cols,b,a,n,r,nu,c,v;
  rows:=PartitionsInSameBlock(H,mu);
  if H.IsSpecht then cols:=Filtered(rows,b->IsERegular(H,b));
  else cols:=rows;
  fi;
  b:=EmptyCrystalizedDecompositionMatrix(H,rows,cols);
  # I experiemented using the Aq() function but it does't seem to be faster...
  if false and H.IsSpecht then ## regular columns so we'll try using Aq() for speed
    v:=H.info.Indeterminate;
    for c in [1..Length(b.cols)] do
      a:=H.Aq(b.cols[c]);
      n:=1;
      while n<Length(a.parts) do
        if a.coeffs[Length(a.parts)-n].valuation>0 then n:=n+1;
        else
          r := Copy( a.coeffs[Length(a.parts)-n] );
          nu:=a.parts[Length(a.parts)-n];
          if Length(r.coefficients) < 1-r.coefficients then
            Append(r.coefficients,List([1..Length(r.coefficients)-1-r.valuation],i->0));
          fi;
          r.coefficients := r.coefficients{[1..1-r.valuation]};
          Append(r.coefficients, Reversed(r.coefficients{[1..-r.valuation]}));
          a:=a-r*b.P(b,nu);
          if nu in a.parts then n:= n+1; fi;
        fi;
      od;
      r := List(a.coeffs, s->s <> 0 * v ^ 0);
      if false in r then
        a.coeffs := ListBlist( a.coeffs, r );
        a.parts := ListBlist( a.parts, r );
      fi;
      b.operations.AddIndecomposable(b,a,false);
    od;
  else
    for c in [1..Length(b.cols)] do
      b.operations.AddIndecomposable(b,H.Pq(b.cols[c]),false);
    od;
  fi;
  return b;
end;

#U MultiplicityFreeBlockDecompositionMatrix(H,mu).........................
#M Compute a decomposition matrix using the Janzten sum formula assuming  
#M that the multiplicities are always zero or one.                        
MultiplicityFreeBlockDecompositionMatrix:=function(H,mu)
  local rows, b, spechts, np, specht, schaper, nu, p;
  rows:=Set(PartitionsInSameBlock(H,mu));
  if not H.IsSpecht then b:=EmptyDecompositionMatrix(H,rows,rows);
  else b:=EmptyDecompositionMatrix(H,rows,Filtered(rows,nu->IsERegular(H,nu)));
  fi;
  spechts:=[];
  for nu in Reversed([1..Length(b.rows)]) do
    np:=Position(b.cols, b.rows[nu]);
    if IsInt(np) then 
      specht:=H.D(b.rows[nu]); 
      b.d[np]:=rec(parts:=[],coeffs:=[]);
    else specht:=0*H.D();
    fi;
    schaper:=Schaper(H,b.rows[nu]);
    if schaper<>0*schaper then
      for p in [1..Length(schaper.parts)] do
        specht:=specht+schaper.coeffs[p]*spechts[Position(b.rows,schaper.parts[p])];
      od;
      specht.coeffs:=specht.coeffs*0+1; # change any n's to 1 
    fi;
    spechts[nu]:=specht;
    for p in specht.parts do
      np:=Position(b.cols,p); 
      AddSet(b.d[np].parts, nu);
      Add(b.d[np].coeffs,1);
    od;
  od;
  return b;
end;
    

#U LinkageClassesDecompositionMatrix(d)...................................
#M Return the list of submatrices corresponding to linkage classes in the 
#M decompostion matrix <d>. If d is the full decomposition matrix of      
#M S(n,n) or H(n) then this is determined by the e-cores; however, in     
#M general there may be more blocks.                                      
LinkageClassesDecompositionMatrix:=function(d) 
  local blocks, newblocks, b, rows, cols, c, r;
  # start by putting each partition into a singleton block and then
  # merge these "blocks" whenever they are linked by a column of d
  blocks:=[1..Length(d.rows)];
  for c in [1..Length(d.cols)] do
    newblocks:=Set(blocks{d.d[c].parts});
    if Length(newblocks)>1 then
      b:=Minimum(newblocks);
      # merge all of these blocks into one
      for r in [1..Length(blocks)] do
        if blocks[r] in newblocks then blocks[r]:=b;
        fi;
      od;
    fi;
  od;
  # now break up the decomposition matrix according to <blocks>
  newblocks:=[];
  for b in Set(blocks) do
    rows:=d.rows{Filtered([1..Length(d.rows)], r->blocks[r]=b)};
    cols:=Filtered(d.cols, c->c in rows);
    Add(newblocks, SubmatrixDecompositionMatrix(d,rows,cols));
  od;
  # finally, sort the blocks according to their core
  for b in newblocks do
    b.core:=ECore(d.H.e,b.rows[1]);
  od;
  Sort(newblocks, function(a,b) return Dominates(a.core,b.core);end);
  return newblocks;
end;

#C Return the submatrix of <d> indexed by the partitions of length
## less than or equal to <len>
RowLimitedDecompositionMatrix:=function(d,len) local p;
  return SubmatrixDecompositionMatrix(d,
           Filtered(d.rows, p->Length(p)<=len),
           Filtered(d.cols, p->Length(p)<=len),
	   Concatenation("length<=", String(len))
	 );
end;

#C Tries to calulcate the decomposition matrix d_{H,n} from scratch.
## At present will return only those columns indexed by the partitions
## of e-weight less than 2.
CalculateDecompositionMatrix:=function(H,n) local d, c, Px;
  if H.IsSpecht then c:=ERegularPartitions(H.e,n);
  else c:=Partitions(n);
  fi;
  d:=H.operations.NewDecompositionMatrix(Partitions(n),c,true);
  for c in [1..Length(d.cols)] do
    if not IsBound(d.d[c]) then
      Px:=H.operations.P.S(H.operations.New("P",[1],[d.cols[c]]),true);
      if Px<>false then d.operations.AddIndecomposable(d,Px,false);
       else Print("# Projective indecomposable P(",
                  TightStringList(d.cols[c]),") not known.\n");
      fi;
    fi;
  od;
  return d;
end;

#U CrystalizedDecompositionMatrix(H,n)
#U CrystalizedDecompositionMatrix(H,n, ordering)
#M Returns a crystallized decomposition matrix
CrystalizedDecompositionMatrix:=function(arg) local H, d, Px, c, n;
  if arg=[] or not (IsRec(arg[1]) and IsBound(arg[1].IsSpecht)) then
    Error("usage, CrystalizedDecompositionMatrix(<H>,<n>)");
  elif arg[1].p<>0 then
    Error("Crystal decomosition matrices are defined only when p=0\n");
  fi;
  H:=arg[1];

  if IsInt(arg[2]) then
    n:=arg[2];
    if IsBound(arg[3]) and IsFunc(arg[3]) then H.Ordering:=arg[3];
    fi;
  else Error("usage, CrystalizedDecompositionMatrix(<H>,<n>)");
  fi;

  d:=H.operations.ReadDecompositionMatrix(H,n,true);
  if d<>false then d:=ShallowCopy(d);
  elif H.IsSpecht then
    d:=H.operations.NewDecompositionMatrix(
              Partitions(n),ERegularPartitions(H.e,n),false);
  else
    c:=Partitions(n);
    d:=H.operations.NewDecompositionMatrix(c,c,false);
  fi;
  for c in [1..Length(d.cols)] do
    if not IsBound(d.d[c]) then
      d.operations.AddIndecomposable(d,H.Pq(d.cols[c]),false);
    fi;
  od;
  return d;
end; 

#U CartanDeterminant(H, n)................................................
#U CartanDeterminant(d)...................................................
#M In the first form, this returns the determinant of the Cartan matrix   
#M for <H>_<n>; in the second form, it returns the Cartan determinant of  
#M the decomposition matrix <d>.                                          
CartanDeterminant:=function(arg) local d;
  if Length(arg)=1 and IsDecompositionMatrix(arg[1]) then
    d:=MatrixDecompositionMatrix( arg[1] );
  elif IsRec(arg[1]) and IsBound(arg[1].operations) and
    arg[1].operations.name="SpechtOps" then
    d:=MatrixDecompositionMatrix(DecompositionMatrix(arg[1], arg[2]));
  fi;
  return DeterminantMat(TransposedMat(d)*d);
end;
  
## Front end to the Induced function inside d.operations. The reason it is
## done this way is more historical than good sense.
#U InducedDecompositionMatrix(d)
#M Computes, as much as possible, the decomposition matrix obtained by
#M  inducing the PIMs described by the columns of <d>.
InducedDecompositionMatrix:=function(d)
  return d.operations.Induced(d);
end;

#C Returns the inverse of (the e-regular part of) d. We invert the matrix
## 'by hand' because the matrix routines can't handle polynomial entries.
## This should be much faster than it is???
#U InvertDecompositionMatrix(d)
#M Returns the inverse of the decomposition matrix <d>. 
InvertDecompositionMatrix:=function(d) local inverse, c, r;
  inverse:=d.H.operations.NewDecompositionMatrix(d.cols,d.cols,
                                             d.IsDecompositionMatrix);
  inverse.operations:=Copy(inverse.operations);
  Unbind(inverse.operations.Induce);

  ## for some reason I can't put this inside the second loop (deleting
  ## the first because d.inverse is not updated this way around...).
  for c in [1..Length(inverse.cols)] do
    d.operations.Invert(d,d.cols[c]);
  od;
  for c in [1..Length(inverse.cols)] do
    if IsBound(d.inverse[c]) then
      inverse.d[c]:=rec(parts:=[], coeffs:=[]);
      for r in [1..c] do
        if IsBound(d.inverse[r]) and c in d.inverse[r].parts then 
          Add(inverse.d[c].parts,r);
          Add(inverse.d[c].coeffs,
              d.inverse[r].coeffs[Position(d.inverse[r].parts,c)]);
        fi;
      od;
      if inverse.d[c]=rec(parts:=[], coeffs:=[]) then Unbind(inverse.d[c]); fi;
    fi;
  od;
  inverse.matname:="Inverse matrix";
  return inverse;
end;
  
## Saves a full decomposition matrix; actually, only the d, rows, and cols
## records components are saved and the rest calculated when read back in.
## The decomposition matrices are saved in the following format:
##   A_Specht_Decomposition_Matrix:=rec(
##   d:=[[r1,...,rk,d1,...dk],[...],...[]],rows:=[..],cols:=[...]);
## where r1,...,rk are the rows in the first column with corresponding
## decomposition numbers d1,...,dk (if di is a polynomial then it is saved
## as a list [di.valuation,<sequence of di.coffcients]; in particular we
## don't save the polynomial name).
#U SaveDecompositionMatrix(<d>) 
#U SaveDecompositionMatrix(<d>,<filename>)
#M Saves the decomposition matrix <d> in a compact format which can be
#M reader back in by SPECHT; see ReadDecompositionMatrix().
SaveDecompositionMatrix:=function(arg)
  local d,TightList,n,file,SaveDm,size, r, c,str;

  if not Length(arg) in [1,2] or not IsDecompositionMatrix(arg[1])
  or (Length(arg)=2 and not IsString(arg[2])) then
    Error("usage: SaveDecompositionMatrix(<d>)\n",
          "   or  SaveDecompositionMatrix(<d>,<filename>)\n");
  fi;
  d:=arg[1];
  n:=Sum(d.rows[1]);
  if Length(arg)=2 then file:=arg[2];
  elif IsBound(d.matname) then
    file:=Concatenation(d.H.HeckeRing,".",d.matname{[1]},String(n));
  elif d.IsDecompositionMatrix then
    file:=Concatenation(d.H.HeckeRing,".",String(n));
  else  ## crystallized decomposition matrix
    file:=Concatenation("e", String(d.H.e), "crys.", String(n));
  fi;

  size:=SizeScreen();    ## SizeScreen(0 shouldn't affect PrintTo()
  SizeScreen([80,40]);  ## but it does; this is our protection.
  
  TightList:=function(list) local l, str;
    str:="[";
    for l in list{[1..Length(list)]} do 
      if IsList(l) then 
        Print(str); 
        TightList(l);
      else Print(str,l); 
      fi;
      str:=",";
    od;
    Print("]"); 
  end;

  if d=false then Error("SaveDecompositionMatrix(<d>), d=false!!!\n");
  elif Length(arg)=1 and d.H.HeckeRing="unknown" then
    Print("SaveDecompositionMatrix(d): \n     the base ring of the Hecke ",
          "algebra is unknown.\n     You must set <d>.H.HeckeRing in ",
          "order to save <d>.\n");
  fi;

  SaveDm:=function()
    Print("## This is a GAP library file generated by \n## SPECHT ", 
          d.H.info.version, "\n\n## This file contains ");
    if IsBound(d.matname) then
      Print("a(n) ", d.matname, " for n = ", Sum(d.rows[1]),"\n");
    else
      if not d.IsDecompositionMatrix then Print("the crystallized "); fi;
      Print("the decomposition matrix\n## of the ");
      if d.H.IsSpecht then
        if d.H.e<>d.H.p then Print("Hecke algebra of "); 
        else Print("symmetric group ");
        fi;
      else Print("q-Schur algebra of "); 
      fi;
      Print("Sym(",n,") over a field\n## ");
      if d.H.p=0 then Print("of characteristic 0 with ");
      elif d.H.p=d.H.e then Print("of characteristic ",d.H.p,".\n\n");
      else Print("with HeckeRing = ", d.H.HeckeRing, ", and ");
      fi;
      if d.H.p<>d.H.e then Print("e=", d.H.e, ".\n\n");fi;
    fi;

    Print("A_Specht_Decomposition_Matrix:=rec(\nd:=[");
    str:="[";
    for c in [1..Length(d.cols)] do
      if not IsBound(d.d[c]) then Print(str,"]");
      else
        for r in d.d[c].coeffs do
          if IsPolynomial(r) then
          Print(str,"[",r.valuation,",",TightStringList(r.coefficients),"]");
          else Print(str,r); 
          fi;
          str:=",";
        od;
        for r in d.d[c].parts do
          Print(str,r);
        od;
        Print("]"); 
        str:=",[";
      fi;
    od;
    Print("],rows:="); TightList(d.rows);
    Print(",cols:="); TightList(d.cols);
    if not d.IsDecompositionMatrix then
      Print(",crystal:=true");
    fi;
    if IsBound(d.matname) then Print(",matname:=\"",d.matname,"\""); fi;
    Print(");\n");
  end;

  ## the actual saving of d
  SpechtInfo("#I* ", ReadIndent, "SaveDecompositionMatrix( \"",
            file, "\")\n");
  PrintTo(file,SaveDm());
  
  ## now we put d into DecompositionMatrices
  if not IsBound(d.matname) then d.operations.Store(d); fi;

  SizeScreen(size); # restore screen.
end; # SaveDecompositionMatrix()

#U AdjustmentMatrix(dp,d) 
#M Returns the 'adjustment matrix' A forthe pair of decomposition
#M matrices <d> and <dp>; thus, <dp>=<d>*A.
AdjustmentMatrix:=function(dp,d) local ad, c, x;
  if d.rows<>dp.rows then return false; fi;

  ad:=dp.H.operations.NewDecompositionMatrix(d.cols,dp.cols,true);
  ad.operations:=Copy(ad.operations);
  Unbind(ad.operations.Induce);
  ad.matname:="Adjustment matrix";
  c:=1;
  while ad<>false and c<=Length(dp.cols) do
    if dp.cols[c] in d.cols then
      x:=dp.P(dp, dp.cols[c]);
      x.H:=d.H;
      x:=d.H.P(d,x);
      if x=false then ad:=false;
      else d.operations.AddIndecomposable(ad,x,false);
      fi;
    fi;
    c:=c+1;
  od;
  return ad;
end;

#U MatrixDecompositionMatrix(d)
#M Returns a GAP matrix for the decomposition matrix <d>. Note that 
#M the rows and columns and <d> are ordered according to H.info.Ordering.
MatrixDecompositionMatrix:=function(d) local r,c, rows, cols, m;
  rows:=Copy(d.rows);
  if d.H.Ordering<>Lexicographic then
    Sort(rows,d.H.Ordering);
    rows:=List(rows,r->Position(d.rows,r));
  else rows:=[Length(rows),Length(rows)-1..1];
  fi;
  cols:=Copy(d.cols);
  if d.H.Ordering<>Lexicographic then
    Sort(cols,d.H.Ordering);
    cols:=List(cols,r->Position(d.cols,r));
  else cols:=[Length(cols),Length(cols)-1..1];
  fi;
  m:=[];
  for r in [1..Length(rows)] do
    m[r]:=[];
    for c in [1..Length(cols)] do
      if IsBound(d.d[cols[c]]) and rows[r] in d.d[cols[c]].parts then
        m[r][c]:=d.d[cols[c]].coeffs[Position(d.d[cols[c]].parts,rows[r])];
      else m[r][c]:=0;
      fi;
    od;
  od;
  return m;
end;

#U DecompositionMatrixMatrix(H,m,n)
#U DecompositionMatrixMatrix(H,m,rows,cols)
#M Given a matrix <m>, and a Specht record <H>, return the corresponding SPECHT 
#M decomposition matrix. In the second form the partitions which index the rows 
#M and columns of <m> are passed explicitly as lists; in the first
#M form these are inferred from the size of <m> and the integer <n>.
DecompositionMatrixMatrix:=function(arg) 
  local H, m, n, rows, cols, d, c, r;
  H:=arg[1];
  m:=arg[2];
  if not (IsBound(H.IsSpecht) and IsMatrix(m) 
  and ( (Length(arg)=3 and IsInt(arg[3])) 
  or (Length(arg)=4 and IsList(arg[3]) and IsList(arg[4])) )) then
    Error("usage: DecompositionMatrixMatrix(H,m, [n|row,cols])");
  fi;
  if Length(arg)=3 then
    n:=arg[3];
    rows:=Partitions(n);
    cols:=ERegularPartitions(H,n);
    if Length(rows)<>Length(m) then rows:=cols; fi;
    if Length(cols)<>Length(m[1]) then cols:=rows; fi;
  else
    rows:=arg[3];
    cols:=arg[4];
  fi;
  if Length(rows)<>Length(m) or Length(cols)<>Length(m[1]) then
    Error("usage: DecompositionMatrixMatrix(H,m, [n|row,cols])");
  fi;
  # test to see if matrix is crystalized
  d:=H.operations.NewDecompositionMatrix(Set(rows), Set(cols), 
            ForAll(m, r->ForAll(r, c->IsInt(c))) );

  rows:=List(rows, r->Position(d.rows,r));
  cols:=List(cols, r->Position(d.cols,r));
  for c in [1..Length(cols)] do
     d.d[c]:=rec(parts:=[], coeffs:=[]);
     for r in [1..Length(rows)] do
       if IsBound(m[rows[r]][cols[c]])     ## maybe a polynomial
       and m[rows[r]][cols[c]]<>0*m[rows[r]][cols[c]] then 
         Add(d.d[c].parts, r);
         Add(d.d[c].coeffs, m[rows[r]][cols[c]]); 
       fi;
     od;
     if d.d[c].parts=[] then Unbind(d.d[c]); fi;
  od;
  return d;
end;

#U FormalCharactersDecompositionMatrix(d)
#M return the matrix giving the formal characters for all of the *simple* 
#M modules as determined by the decomposition matrix <d>. This only
#M makes sense if <d> is the decomposition matrix of a
#M quasi-hereditary quotient of a Schur algebra.
## Note that we need to conjugate the partitions because the Weyl
## module S.W(mu) should be labelled S.W(mu').
FormalCharactersDecompositionMatrix:=function(d)
  local S, chs, Fmu, mu, nu, l;
  S:=d.H;
  chs:=List(d.rows, mu->[1..Length(d.cols)]*0);
  for mu in [1..Length(d.rows)] do
    Fmu:=S.W(S.F(d,d.rows[mu]));
    for nu in [mu..Length(d.rows)] do
      chs[mu][nu]:=Sum([1..Length(Fmu.parts)], l->Fmu.coeffs[l]
          *Length(SemistandardTableaux(ConjugatePartition(Fmu.parts[l]),
                                       ConjugatePartition(d.cols[nu]))));
    od;
  od;
  return DecompositionMatrixMatrix(S,chs,d.rows,d.cols);
end;


#U SSRestrictedModule(x,k1,r1,k2,r2...)                                 
#M This does what I thought that SRestrictedModule did; namely, it      
#M applies the restriction sequence r1^k1 r2^k2... to the module <x>.   
SSRestrictedModule:=function(arg)
  local m, r, res;
  m:=arg[1];
  if Length(arg)=2 and IsList(arg[2]) then res:=arg[2];
  elif IsOddInt(Length(arg)) then res:=arg{[2..Length(arg)]};
  else Error("usage:, SSRestrictedModule(x,k1,r1,k2,r2...)");
  fi;
  r:=1;
  while r<Length(res) do
    m:=SRestrictedModule(m,res[r], res[r+1]);
    r:=r+2;
  od;
  return m;
end;

#U RestrictedBlockDecompositionMatrix(B, i1 [,i2,...,ik])                 
#U RestrictedBlockDecompositionMatrix(B, b, i1 [,i2,...,ik])              
#M Compute and return as much as possible of the decomposition matrix     
#M of some block obtained by restricting the columns of the decomposition 
#M matrix <b> according to the residue sequence <i1>, ..., <ik>.          
#M In the second form, a (partially complete) form the of block that is   
#M being restricted to can also be given, in which case any new columns   
#M are added to this matrix, which is then returned.                      
RestrictedBlockDecompositionMatrix:=function(arg)
  local B, res, resP, resB, mu, rows, misses, stillalive, DN, missus;
  B:=arg[1];
  DN:=SPECHT.DecompositionNumberOK;
  SPECHT.DecompositionNumberOK:=false;
  if IsDecompositionMatrix(arg[2]) then
    resB:=arg[2];
    res:=arg{[3..Length(arg)]};
  else res:=arg{[2..Length(arg)]};
  fi;
  if Length(res)=1 then res:=[1,res[1]]; fi;
  if not IsDecompositionMatrix(B) or res=[] or not ForAll(res, IsInt) then
    Error("usage, RestrictDecompositionMatrix(B, i1 [,i2,...,ik]),");
  fi;
  resP:=SSRestrictedModule(B.P(B,B.cols[1]), res);
  if resP=0*resP then 
    return EmptyDecompositionMatrix(B.H, [],[]);
  fi;
  if not IsBound(resB) then
    rows:=PartitionsInSameBlock(B.H,resP.parts[1]);
    if B.H.IsSpecht then
      resB:=EmptyDecompositionMatrix(B.H, rows, 
	                           Filtered(rows, mu->IsERegular(B.H,mu)));
    else 
      resB:=EmptyDecompositionMatrix(B.H,rows, rows);
    fi;
    resB.IsDecompositionMatrix:=B.IsDecompositionMatrix;
  fi;
  misses:=[];   # these will be the projectives that we cannot
                # decompose on the first attempt
  if IsNewIndecomposable(resB,resP) then AddIndecomposable(resB,resP); 
  else Add(misses, resP);
  fi;
  for mu in [2..Length(B.cols)] do
    if IsBound(B.d[mu]) then
      resP:=SSRestrictedModule(B.P(B,B.cols[mu]), res);
      Print("resP = ", resP,"\n");
      if IsNewIndecomposable(resB,resP) then AddIndecomposable(resB,resP);
      else Add(misses, resP); 
      fi;
    fi;
    #Print(".\c");
  od;
  # Finally, we look again at those projectives which we could not
  # decompose; perhaps we can decompose some of them now.
  stillalive:=true;
  while stillalive and misses<>[] and MissingIndecomposables(resB,true)<>[] do
    SpechtInfo("# ...looping through missed projectives...\n");
    missus:=[]; # the new misses
    stillalive:=false;  # this will become true only if we find a new PIM
    for resP in misses do
      SpechtInfo("resP = ",resP,"\n");
      if IsNewIndecomposable(resB,resP) then 
        AddIndecomposable(resB,resP);
	stillalive:=true;
	SpechtInfo("new one!!\n");
      else Add(missus, resP); 
      fi;
    od;

    misses:=missus;
    #Print(".\c");
  od;
  Print("\n");
  SPECHT.DecompositionNumberOK:=DN;
  return resB;
end;

#C The following functions are for changing the labelling of the rows    
#C of decomposition matrices.                                            

#U NormalLabels(H)                                                       
#M Return to the default labelling of the decomposition matrix.          
NormalLabels:=function(H)
  LabelPartition:=RealLabelPartition;
  H.Ordering:=Lexicographic;
end;

#U QuotientLabels(H)                                                      
#M Change labelling of decomposition matrices so that uses the e-quotient.
QuotientLabels:=function(H)
   LabelPartition:=function(mu) return LabelByEQuotient(H.e, mu); end;
   H.Ordering:=function(mu,nu)
     return LengthLexicographicallyByQuotient(H.e,mu,nu);
   end;
end;

#U RichardsLabels(H)                                                      
#M Change labelling of decomposition matrices so that uses Richard's      
#M symbols.                                                               
RichardsLabels:=function(H)
   LabelPartition:=function(mu) return RichardsLabelPartition(H.e, mu); end;
   H.Ordering:=function(mu,nu) return RichardsOrdering(H.e,mu,nu); end;
end;

#U CanonicalLabelling(H)
CanonicalLabelling:=function(H)
   LabelPartition:=function(mu) 
     return LabelCanonicalQuotient(H.e, mu); 
   end;
   #H.Ordering:=function(mu,nu) 
   #  return CanonicalQuotientLengthLexicographic(H.e,mu,nu); 
   #end;
end;

#U IsColumnMinimal(d).....................................................
#M Returns true if the decomposition matrix <d> has the property that no  
#M can be subtracted from a previous one without making some of the       
#M entries negative.                                                      
IsColumnMinimal:=function(d)
  local P, cols, mu, nu;
  if not IsDecompositionMatrix(d) then Error("usage: IsColumnMinimal(d)"); fi;
  P:=List(d.cols, mu->d.P(d,mu));
  # have to avoid the first "trivial" column for Schur algebras
  if Length(d.rows)=Length(d.cols) then cols:=2; else cols:=1; fi;
  for mu in [cols..Length(d.cols)-1] do
    for nu in [mu+1..Length(d.cols)] do
      if PositiveCoefficients(P[nu]-P[mu]) then return false; fi;
    od;
  od;
  return true;
end;

#U MaximumsDecompositionMatrix(b).........................................
#M Return the maximum entries in a decomposition matrix along with the    
#M partitions which index them.                                           
MaximumsDecompositionMatrix:=function(b)
  local M, maxes, rows, mu, i;
  M:=Maximum(Flat(MatrixDecompositionMatrix(b)));
  maxes:=[M];
  for mu in [1..Length(b.d)] do
    if M in b.d[mu].coeffs then
      rows:=Filtered([1..Length(b.d[mu].coeffs)],i->b.d[mu].coeffs[i]=M);
      Add(maxes, [b.cols[mu], b.rows{rows}]);
    fi;
  od;
  return maxes;
end;

#U FindAllCrystalizedDecompositionNumbers(w)                              
#M Generate a list of all crystalized decomposition numbers for a given   
#M weight w using Fayers' result from his "generalized runner removal"    
#M paper.                                                                 
FindAllCrystalizedDecompositionNumbers:=function(w)
  local H, v, currentPis, nextPis, c, dqs, level, qMat, found, pdqs, mq, p, q;
  H:=Specht(2*w);
  v:=H.info.Indeterminate;
  currentPis:=[];
  # We start with the principle block of weight w when e=2w and then move
  # up the adjacant pyramids. This is more efficient than first generating
  # all pyramids, especially as we throw out the conjugate pyramids along
  # the way.
  # Initially did this using sets for next and current Pis but for some
  # reason this lead to duplicate pyramids being stored in these lists.
  nextPis:=[ PyramidPartition(H,[w*H.e]) ];
  c:=0;
  dqs:=[0];                          # to stop 0*v^0 being added !
  while nextPis<>[] do
    currentPis:=nextPis;
    nextPis:=[];
    level:=Sum(Flat(currentPis[1].a));
    Print(level,"--------------------------------------\n");
    for p in [1..Length(currentPis)] do
      qMat:=CrystalizedDecompositionMatrixPyramid(H,currentPis[p]);
      found:=0;
      pdqs:=Set(Flat(MatrixDecompositionMatrix(qMat)));
      for q in pdqs do
        if not q in dqs then
          AddSet(dqs,q);
          found:=found+1;
          mq:=v^w*Value(q*v^0,v^-1); # poly from m(p)
          if not mq in dqs then
            AddSet(dqs,mq);
            found:=found+1;
          fi;
          if found=0 then Print(currentPis[p],"\n"); fi;
        fi;
      od;
      c:=c+1;
      Print("   ",String(p,3),"/",String(Length(currentPis),-3)," (",
            String(SPrint(level,")"),-4)," :  ",pdqs, "\n                 :  ",dqs);
      if found>0 then Print(" - ", found, " new!\n");
      else Print("\n");
      fi;
      for q in SmallerAdjacentPyramids(currentPis[p]) do
        if not (q in nextPis or ConjugatePyramid(q) in nextPis) then 
          Add(nextPis,q); 
        fi;
      od;
    od;
  od;
  return dqs;
end;

