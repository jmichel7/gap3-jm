#M#####################################################################
#M  SPECHT - specht.g : the kernel of SPECHT                         ##
#M                                                                   ##
#M     A GAP package for calculating the decomposition numbers of    ##
#M     Hecke algebras of type A (over fields of characteristic       ##
#M     zero). The functions provided are primarily combinatorial in  ##
#M     nature. Many of the combinatorial tools from the (modular)    ##
#M     representation theory of the symmetric groups appear in the   ##
#M     package.                                                      ##
#M                                                                   ##
#M     These programs, and the enclosed libraries, are distributed   ##
#M     under the usual licensing agreements and conditions of GAP.   ##
#M                                                                   ##
#M     Andrew Mathas                                                 ##
#M                                                                   ##
#M#####################################################################

#D The heart of the SPECHT package: defines the Specht record structures

#M New changes
#M  - removed cached uses of CrystalMatrices[] as this was having the
#M    unfortunate side effect of forcing Specht to compute the set of all 
#M    partitions of n when it computed Pq(mu), for mu a partition of n.
#M    In turn this significantly slowed down the calculation of the
#M    crystalized decompositon matrices of large blocks indexed by
#M    pyramids when doing calculations with a fixed weight.
#M  - implemented a version of the LTL algorithm for the q-Schur
#M    algebras in characteristic zero. What we actually do is that we
#M    use my "equating decomposition numbers for different primes"
#M    paper with Gordon James to compute singular singular columns
#M    using a larger regular partition for an (e+1)th root of unity.
#M    As a result we now cache crystal decomposition numbers (this has
#M    not been extensively tested so there could still be some bugs here).
#M  - added <H>.PartitionLength to limit LLT calculations to
#M    partitions of a given length - untested.
#M  - amalgamated various version of files lying around; new
#M    functions available include:
#M      DominatedDecompositionMatrix()
#M      SubmatrixDecompositionMatrix()
#M      BlockDecompositionMatrix()
#M      RowLimitedDecompositionMatrix()
#M      IntegralGramDeterminant()
#M    These functions have now been moved to decomp.g.
#M  - changed PrintModule() to use H.Ordering()
#M  - added CartanMatrix() and CartanDeterminant()
#M  - added ProjectiveDimension()

#M Change log
#M 2.4:
#M  - fixed more bugs in H.valuation; returned incorrect answers before
#M    when e=0 or e=p (symmetric group case).
#M  - fixed bug in Dq(), via Sasha Kleshchev.

#M 2.3:
#M  - fixed bug in H.valuation  reported by Johannes Lipp
#M  - fixed bug in Sq() reported by Johannes Lipp.

#M 2.2: June 1996: various changes requested by the referee.
#M  - mainly changing function names.
#M  - DecompositionMatrix changed so that it no longer attempts to
#M    calculate decomposition matrices in the finite field case; added
#M    CalculateDecompositionMatrix() to do this.
#M  - replaced Matrix() function with MatrixDecompositionMatrix() and
#M    added the function DecompositionMatrixMatrix().

#M 2.1: April 1996:
#M  - Added a filename argument to SaveDecompositionMatrix and made it
#M    save the non-decomposition matrices under (more) sensible names; the 
#M    later is done using the existence of a record component d.matname.
#M    Need to do something here about reading such matrices in (this can be
#M    done using DecompositionMatrix)..
#M  - Changed ReadDecompositionmatrix so that it automatically reads the
#M    format of version 1.0 files (and deleted ReadOldDecompositionMatrix).
#M    Also fixed a bug here which was causing confusion between Hecke algebra
#M    and Schur algebra matrices.
#M  - Renamed FindDecompositionMatrix as KnownDecompositionMatrix because it
#M    doesn't try to calculate a crytalized decomposition matrix; the new
#M    FindDecompositionMatrix will return the decomposition matrix if at all
#M    possible. In particular, this fixes a 'bug' in SimpleDimension.
#M  - Rewrote AdjustmentMatrix so that it actually works.
#M  - Changed P()->S() module conversions so that it no longer requires
#M    the projective module to have positive coefficients (and fixed bug in
#M    matrix ops).

#M 2.0: March 1996: 
#M   - LLT algorithm implemented for calculating crystal basis of
#M     the Fock space and hence by specialization the decomposition
#M     matrices for all Hecke algebras over fields of characteristic
#M     zero; most of the work is done by the internal function Pq(). 
#M     This required a new set of 'module' types ("Pq", "Sq", and "Dq"),
#M     with correspondining operation sets H.operations.X. In particular,
#M     these include q-induction and q-restriction, effectively allowing
#M     computations with $U_q(\hat sl_e)$ on the Fock space.
#M   - crystallized decomposition matrices added and decomposition
#M     matrix type overhauled. d.d[c] is now a list of records with
#M     partitions replaced by references to d.rows etc. Changed format
#M     of decomposition matrix library files so as to handle polynomial
#M     entries in the matrices..
#M   - printing of Specht records changed allowing more compact and
#M     readable notation. eg. S(1,1,1,1)->S(1^4). (see SpechtPrintFn).
#M   - reversed order of parts and coeffs records in H.S(), d.d etc;
#M     these lists are now sets which improves use of Position().
#M   - reorganised the Specht function and the record H() which it returns; 
#M     in particular added H.info.
#M   - extended SInducedModule() and SRestrict to allow multiple inducing and
#M     restricting (all residues).

#M 1.0: December 1995: initial release.

#M####################################################################

#M The infamous Specht identifier (for the modules H.S etc, not the
#M records returned by Specht); not completely robust but...
IsSpecht:=function(x)
  return IsRec(x) and IsBound(x.parts);
end;

#M####################################################################

#M Here is a description of the structure of the main records used 
#M in SPECHT

#M 1. Specht()
#M Specht() is the main function in specht.g, it returns a record 'H'
#M which represents the family of Hecke algebras, or Schur algebras, 
#M corresponding to some fixed field <R> and parameter <q>. 'H' has 
#M the components:
#M   IsSpecht      : is either 'true' or 'false' depending on whether
#M                   'H' is a Hecke algebra or Schur algebra resp.
#M   S(), P(), D() : these three functions return records which 
#M                   represent Specht modules, PIMs, and simple
#M                   'H' modules repectively. These functions exist
#M                   only when 'H' is a Hecke algebra record.
#M   W(), P(), F() : these are the corresponding functions for Weyl
#M                   modules, PIMs, and simple modules of Schur algbras.
#M   info          : this is a record with components
#M                     version: SPECHT version number,
#M                     Library: path to SPECHT library files
#M                     SpechtDirectory: addition directory searched by
#M                            SPECHT (defaults to current directory)
#M                     Indeterminate: the indedeterminate used by SPECHT
#M                            in the LLT algorithm when H.p=0.
#M   operations    : apart from the obvious things like the Print()
#M                   function for 'H' this record also contains the
#M                   operation records for the modules S(), P() etc.
#M                   as well as functions for manipulating these records
#M                   and accessing decomposition matrix files. The most
#M                   most important of these are:
#M                     S, P, D, Pq, Sq, Dq : operations records for modules
#M                     New : creation function for modules. Internally
#M                       modules are created via
#M                         H.operations.new(module,coeffs,parts)
#M                       where module is one of "S", "P", "D", "Sq", "Pq",
#M                       or "Dq" ("S" and "D" are used even for Schur
#M                       algebras), coeffs is an integer or a *set* of
#M                       integers and parts is a partition or a *set* of
#M                       partitions. In any programs the use of New is 
#M                       better than H.S(mu), for example, because the 
#M                       function names are different for Hecke and Schur 
#M                       algebras. Note that coeffs and parts must be 
#M                       ordered reverse lexicographically (ie. they are 
#M                       *sets*).
#M                     Collect: like New() except that  coeffs and parts
#M                       need not be sets (and may contain repeats).
#M                     NewDecompositionMatrix : creates a decomposition
#M                       matrix.
#M                     ReadDecompositionMatrix : reads, and returns, a 
#M                       decomposition matrix file.
#M                     KnownDecompositionMatrix : returns a decomposition
#M                       matrix of a given size; will either extract this
#M                       matrix from Specht's internal lists or call
#M                       ReadDecompositionMatrix(), or try to calculate 
#M                       the decomposition matrix (without using the
#M                       crystalized decomposition matrices).
#M                     FindDecompositionMatrix : like KnownDM except that
#M                       it will calculate the crystalized DM if needed.
#M   Ordering      : a function for ordering partitions; controls how
#M                   decomposition matrices for H are printed.
#M   e             : order of <q> in <R>
#M   p             : characteristic of <R>
#M   valuation     : the valuation map of [JM2]; used primarily by
#M                   the q-Schaper theorem.
#M   HeckeRing     : bookkeeping string used primarily in searching for
#M                   library files.
#M   Pq(), Sq()    : Functions for computing elements of the Fock space
#M                   when H.p=0 (used in LLT algorithm). Note that there is
#M                   no Dq; also unlike their counter parts S(), P(), and
#M                   D() they accept only partitions as arguments.
#M
#M 2. The module functions S(), P() and D() (and Schur equivalents)
#M These functions return record 'x' which represents some 'H'--module.
#M 'x' is a record with the following components:
#M   H      : a pointer back to the corresponding algebra
#M   module : one of "S", "P", "D", "Sq", "Pq", or "Dq", (not "W", or "F").
#M   coeffs : a *set* of coefficients
#M   parts  : the corresponding *set* of partitions
#M   operations :
#M       + - * / : for algebric manipulations
#M       Print : calls PrintModule
#M       Coefficient : returns the coefficient of a given partition
#M       PositiveCoefficients : true if all coefficients are non-negative
#M       IntegralCoefficients : true if all coefficients are integral
#M       InnerProduct : computes the 'Kronecker' inner product
#M       Induce, Restrict, SInduce, SRestrict : induction and restriction
#M                 functions taylored to 'x'. These functions convert 'x'
#M                 to a linear combination of Specht modules, induce and
#M                 then convert back to the type of 'x' (if possible).
#M                 Quantized versions are applied as appropriate.
#M       S, P, D : functions for rewriting 'x' into the specified type
#M                 (so, for example, S('x') rewrites 'x' as a linear
#M                 combination of Specht modules).
#M 'x'.operations is a pointer to 'H'.operations.('x'.module).

#M Calculates the dimensions of the projective modules in d
#M Usage:  ProjectiveDimension(d)   -> prints all simple dimensions
#M         ProjectiveDimension(H,n) -> prints all again
#M  ProjectiveDimension(H,mu) or ProjectiveDimension(d,mu) -> dim P(mu)
ProjectiveDimension:=function(arg) 
  local d, mu, r, c, x, collabel, M, cols;

  if IsDecompositionMatrix(arg[1]) then 
    d:=arg[1]; mu:=Flat(arg{[2..Length(arg)]});
  elif IsRec(arg[1]) and IsBound(arg[1].IsSpecht) and Length(arg)>1 then
    mu:=Flat(arg{[2..Length(arg)]});
    d:=arg[1].operations.FindDecompositionMatrix(Sum(mu));
    if d=false then 
      SpechtInfo("# ProjectiveDimension(H,n), the decomposition matrix of H_n is ",
            "not known.\n");
      return false;
    fi;
    if Length(mu)=1 and not IsList(arg[2]) then mu:=[]; fi;
  else
    Error("usage, ProjectiveDimension(<d>), ProjectiveDimension(<H>,<n>), ",
          "or ProjectiveDimension(<d>|<H>,<mu>)");
  fi;

  if not d.H.IsSpecht then 
    Error("# ProjectiveDimension() not implemented for Schur algebras\n");
  elif mu=[] then 
    cols:=Copy(d.cols);
    if d.H.Ordering=Lexicographic then
      cols:=cols{[Length(cols),Length(cols)-1..1]};
    else Sort(cols, d.H.Ordering);
    fi;
    cols:=List(cols, c->Position(d.cols,c));
    collabel:=List([1..Length(cols)], c->LabelPartition(d.cols[cols[c]]));
    M:=Maximum(List(collabel, Length))+1;
 
    for c in [1..Length(cols)] do
      Print(String(collabel[c],-M),": ");
      if IsBound(d.pdimensions[cols[c]]) then 
        Print(d.pdimensions[cols[c]],"\n");
      else
        x:=d.P(d,d.cols[cols[c]]);
        if x=false then Print("not known\n");
        else
          d.pdimensions[cols[c]]:=Sum([1..Length(x.parts)],
                             r->x.coeffs[r]*SpechtDimension(x.parts[r]));
          Print(d.pdimensions[cols[c]],"\n");
        fi;
      fi;
    od;
  else
    c:=Position(d.cols,mu); 
    if c=false then 
      SpechtInfo("# ProjectiveDimension(<d>,<mu>), <mu> is not in <d>.cols\n");
      return false;
    else
      if not IsBound(d.pdimensions[c]) then 
        x:=d.H.P(d,d.cols[c]);
        if x=false then return false;
        else d.pdimensions[c]:=Sum([1..Length(x.parts)],
                            r->x.coeffs[r]*SpechtDimension(x.parts[r]));
        fi;
      fi;
      return d.pdimensions[c];
    fi;
  fi;
end;

#M Calculates the dimensions of the simple modules in d
#M Usage:  SimpleDimension(d)   -> prints all simple dimensions
#M         SimpleDimension(H,n) -> prints all again
#M         SimpleDimension(H,mu) or SimpleDimension(d,mu) -> dim D(mu)
SimpleDimension:=function(arg) local d, mu, r, c, x, collabel, M, cols;
  if IsDecompositionMatrix(arg[1]) then 
    d:=arg[1]; mu:=Flat(arg{[2..Length(arg)]});
  elif IsRec(arg[1]) and IsBound(arg[1].IsSpecht) and Length(arg)>1 then
    mu:=Flat(arg{[2..Length(arg)]});
    d:=arg[1].operations.FindDecompositionMatrix(Sum(mu));
    if d=false then 
      SpechtInfo("# SimpleDimension(H,n), the decomposition matrix of H_n is ",
            "not known.\n");
      return false;
    fi;
    if Length(mu)=1 and not IsList(arg[2]) then mu:=[]; fi;
  else
    Error("usage, SimpleDimension(<d>), SimpleDimension(<H>,<n>), ",
          "or SimpleDimension(<d>|<H>,<mu>)");
  fi;

  if not d.H.IsSpecht then 
    Error("# SimpleDimension() not implemented for Schur algebras\n");
  elif mu=[] then 
    cols:=Copy(d.cols);
    if d.H.Ordering=Lexicographic then
      cols:=cols{[Length(cols),Length(cols)-1..1]};
    else Sort(cols, d.H.Ordering);
    fi;
    cols:=List(cols, c->Position(d.cols,c));
    collabel:=List([1..Length(cols)], c->LabelPartition(d.cols[cols[c]]));
    M:=Maximum(List(collabel, Length))+1;
 
    for c in [1..Length(cols)] do
      Print(String(collabel[c],-M),": ");
      if IsBound(d.dimensions[cols[c]]) then 
        Print(d.dimensions[cols[c]],"\n");
      else
        x:=d.H.D(d,d.cols[cols[c]]);
        if x=false then Print("not known\n");
        else
          d.dimensions[cols[c]]:=Sum([1..Length(x.parts)],
                             r->x.coeffs[r]*SpechtDimension(x.parts[r]));
          Print(d.dimensions[cols[c]],"\n");
        fi;
      fi;
    od;
  else
    c:=Position(d.cols,mu); 
    if c=false then 
      SpechtInfo("# SimpleDimension(<d>,<mu>), <mu> is not in <d>.cols\n");
      return false;
    else
      if not IsBound(d.dimensions[c]) then 
        x:=d.H.D(d,d.cols[c]);
        if x=false then return false;
        else d.dimensions[c]:=Sum([1..Length(x.parts)],
                            r->x.coeffs[r]*SpechtDimension(x.parts[r]));
        fi;
      fi;
      return d.dimensions[c];
    fi;
  fi;
end;

#P returns a list of the e-regular partitions occurring in x
ListERegulars:=function(x) local p;
  if x.H.e=0 then return x.parts;
  elif x=0*x then return [];
  else return List(Filtered([Length(x.parts),Length(x.parts)-1..1],
         p->IsERegular(x.H.e,x.parts[p])),p->[x.coeffs[p], x.parts[p]]);
  fi;
end; # ListERegulars

#P Print the e-regular partitions in x if IsSpecht(x); on the other hand,
#M if IsDecompositionMatrix(x) then return the e-regular part of the
#M decomposition matrix.
ERegulars:=function(x) local y, regs, len, r;
  if IsDecompositionMatrix(x) then
    regs:=rec(operations:=x.operations);
    for y in RecFields(x) do 
      if not y in ["d","rows","labels"] then regs.(y):=x.(y); fi;
    od;
    regs.d:=[]; 
    for y in [1..Length(x.cols)] do
      if IsBound(x.d[y]) then
        regs.d[y]:=rec(parts:=[], coeffs:=[]);
        for r in [1..Length(x.d[y].parts)] do
          len:=Position(x.cols,x.rows[x.d[y].parts[r]]);
          if len<>false then 
            Add(regs.d[y].parts,len);
            Add(regs.d[y].coeffs,x.d[y].coeffs[r]);
          fi;
        od;
      fi;
    od;
    regs.rows:=regs.cols;
    return regs;
  else
    len:=0;
    regs:=ListERegulars(x);
    if regs=[] or IsInt(regs[1]) then Print(regs, "\n");
    else
      for y in regs do
        if (len + 5 + 4*Length(y[2])) > 75 then len:=0; Print("\n"); fi;
        if y[1]<>1 then Print(y[1], "*"); len:=len + 3; fi;
        Print(y[2], "  ");
        len:=len + 5 + 4*Length(y[2]);
      od;
      Print("\n");
    fi;
  fi;
end; # eRegular

#M Returns true if S(mu)=D(mu) - note that this implies that mu is e-regular
#M (if mu is not e-regular, false is returned).     -- see [JM2]
#M IsSimle(H,mu)
#M   ** uses H.valuation
IsSimpleModule:=function(arg) local H,mu, mud, simple, r, c, v;
  if arg=[] or not ( IsRec(arg[1]) and IsBound(arg[1].valuation) ) then
    Error("usage, IsSimpleModule(<H>,<mu>)");
  fi;

  H:=arg[1]; mu:=Flat(arg{[2..Length(arg)]});
  if not IsERegular(H.e,mu) then return false; 
  elif mu=[] then return true; fi;

  mud:=ConjugatePartition(mu);
  simple:=true; c:=1;
  while simple and c <=mu[1] do
    v:=H.valuation(mu[1]+mud[c]-c);
    simple:=ForAll([2..mud[c]], r->v=H.valuation(mu[r]+mud[c]-c-r+1));
    c:=c+1;
  od;
  return simple;
end; #IsSimpleModule
  
#M Split an element up into compontents which have the same core.
#M Usage: SplitECores(x) - returns as list of all block components
#M        SplitECores(x,lambda) - returns a list with (i) core lambda,
#M (ii) the same core as lambda, or (iii) the same core as the first
#M element in lambda if IsSpecht(lambda).
SplitECores:=function(arg) local cores, c, cpos, y, cmp;
  if arg=[] or arg[1]=false then return []; 
  elif not IsSpecht(arg[1]) then
    Error("usage, SplitECores(<x>), SplitECores(<x>,<mu>), or ",
              "SplitECores(<x>,<y>)");
  elif arg[1]=0*arg[1] then return [];
  fi;

  if Length(arg)=1 then
    cores:=[]; cmp:=[];
    for y in [1..Length(arg[1].parts)] do
      c:=ECore(arg[1].H.e, arg[1].parts[y]);
      cpos:=Position(cores, c);
      if cpos=false then 
        Add(cores, c); 
        cpos:=Length(cores);
        cmp[cpos]:=[[],[]];
      fi;
      Add(cmp[cpos][1], arg[1].coeffs[y]);
      Add(cmp[cpos][2], arg[1].parts[y]);
    od;
    for y in [1..Length(cmp)] do
      cmp[y]:=arg[1].H.operations.New(arg[1].module,cmp[y][1],cmp[y][2]);
    od;
  else
    if Length(arg)=2 and IsSpecht(arg[2]) then 
      c:=ECore(arg[2].H.e, arg[2].parts[Length(arg[1].parts)]);
    else c:=ECore(arg[1].H.e, Flat(arg{[2..Length(arg)]}));
    fi;
    cmp:=[ [],[] ];
    for y in [1..Length(arg[1].parts)] do
      if ECore(arg[1].H.e, arg[1].parts[y])=c then 
        Add(cmp[1], arg[1].coeffs[y]); 
        Add(cmp[2], arg[1].parts[y]); 
      fi;
    od;
    cmp:=arg[1].H.operations.New(arg[1].module, cmp[1], cmp[2]);
  fi;
  return cmp;
end; #SplitECores

#M Limit the length of the partitions taht Specht computes with. This
#M is useful because Specht computes decomposition matrices by filling in
#M the columns.
LimitPartitionLength:=function(H,l)
  H.PartitionLength:=l;
end;

SpechtOps:=OperationsRecord("SpechtOps");

#M#########################################################################

#M Specht() is the main function in the package, although in truth it is
#M little more than a wrapper for the funcions S(), P(), and D().
#M Originally, I had these as external functions, but decided that it
#M was better to tie these functions to e=H.e as strongly as possible.
Specht:=function(arg) 
  local H, a, i, 
        Handler, AddModules, MultiplyModules, PrintModule, 
        InducedModule, sinduced, SInducedModule, qsinduced, qSInducedModule,
        RestrictedModule, srestricted, SRestrictedModule, qrestricted,
        qSRestrictedModule, InnerProduct, Coefficient, PositiveCoefficients, 
        IntegralCoeffucuents, PositiveCrystalCoefficients, 
        IntegralCrystalCoefficients, PrintDecompositionMatrix, 
        CrystalMatrixOps, DecompositionMatrixOps, v, Pq, Sq, Dq, HSC, BUG;
  
  ## This function is described above. It interprets the function calls
  ## H.X(*) where X = S(), P() or D().
  Handler:=function(module, arg) local usage,mu,d,z,n;

    usage:=function()
      Error("usage: ",module,"(<mu1,mu2,...>), ", module,"(<x>), ",
             module,"(<d>,<mu1,mu2,...>), or ",module,"(<d>,<x>).");
    end;

    if arg=[] then
      return H.operations.New(module,[1],[[]]);
    elif IsDecompositionMatrix(arg[1]) then
      d:=arg[1];
      if Length(arg)=1 then usage();
      elif Length(arg)=2 and IsSpecht(arg[2]) then
        return d.operations.(arg[2].module{[1]}).(module)(d,arg[2]);
      elif Length(arg)=2 and arg[2]=false then return arg[2];
      else mu:=Flat(arg{[2..Length(arg)]});
      fi;
    elif Length(arg)=1 and IsSpecht(arg[1]) then
      return arg[1].operations.(module)(arg[1],false);
    else mu:=Flat(arg);
    fi;

    if not ForAll(mu,z->IsInt(z)) then usage(); fi;
    z:=Copy(mu);
    Sort(mu, function(a,b) return a>b;end); # non-increasing
    if mu<>z then 
      Print("## ",module,"(mu), warning <mu> is not a partition.\n");
    fi;
    if Length(mu)>0 and mu[Length(mu)]<0 then
      Error("## ", module,"(mu): <mu> contains negative parts.\n");
    fi;
    z:=Position(mu,0);
    if z<>false then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu

    if not IsERegular(H.e,mu) and (H.IsSpecht or Length(module)=2)
    and module{[1]}<>"S" then
      Error(module,"(mu): <mu>=[",TightStringList(mu),
              "] must be ", H.e,"-regular\n\n");
    fi;

    if IsBound(d) then
      if module{[1]}="S" then
        return d.operations.S.D(d,d.H.operations.New("S",[1],[mu]));
      else return d.operations.(module{[1]}).S(d,
                     H.operations.New(module,[1],[mu]));
     fi;
    else return H.operations.New(module,[1],[mu]);
    fi;
  end; # Handler

  ## Specht() returns the record H. This record will also contain an
  ## operations record which in turn will contain all of the functions
  ## for operations on the various modules S(), P() and D().
  H:=rec(IsSpecht:=true,

    S:=function(arg) return Handler("S",arg); end,
    P:=function(arg) return Handler("P",arg); end,
    D:=function(arg) return Handler("D",arg); end,
    
    ## bits and pieces about H
    info:=rec(version:=SPECHT.Version),

    ## ordering used when printing decomposition matricwes
    Ordering:=Lexicographic,

    operations:=SpechtOps

  );  ## H:=rec(...)
  Inherit(H.operations,rec(
      \=:=function(H,K) return H.IsSpecht=K.IsSpecht and H.e=K.e and H.p=K.p; end,

      ## for ordering the rows of the decomposition matrices
      ## (as it is common to all decomposition matrices it lives here)
    
      Print:=function(x)
        if H.IsSpecht then Print("Specht("); else Print("Schur("); fi;
        Print("e=", x.e, ", ");
        if x.p<>0   then Print("p=", x.p, ", "); fi;
        if H.IsSpecht then Print("S(), P(), D()");
        else Print("W(), P(), F()");
        fi;
        if IsBound(H.Pq) then Print(", Pq()"); fi;
        if H.p<>H.e and H.p<>0 then 
          Print(", HeckeRing=\"", H.HeckeRing, "\")");
        else Print(")");
        fi;
      end,

      ## The following two functions are used by P(), and elsewhere.
      ##   generate the hook (k,1^n-k)  - as a list - where k=arg
      ##   actually not quite a hook since if arg is a list (n,k1,k2,...)
      ##   this returns (k1,k2,...,1^(n-sum k_i))
      Hook:=function(arg) local n, k, K, i;
        n:=arg[1];
        K:=arg{[2..Length(arg)]};
        k:=Sum(K);
        if k < n then Append(K, List([1..(n-k)], i->1));
        elif k > n then Error("hook, partition ", k, " bigger than ",n, "\n");
        fi;
        return K;
      end,
    
      DoubleHook:=function(n,x,y,a) local s, i;
        s:=[x];
        if y<>0 then Add(s,y); fi;
        if a<>0 then Append(s, List([1..a],i->2)); fi;
        i:=Sum(s);
        if i < n then
          Append(s, List([1..n-i], i->1));
          return s;
        elif i=n then return s;
        else return [];
        fi;
      end,
    
      ## Returns p(n) - p(n-1,1) + p(n-2,1^2) - ... + (-1)^(n-1)*p(1^n).
      ## So, S(mu)*Omega(n) is the linear combination of the S(nu)'s where
      ## nu is obtained by wrapping an n-hook onto mu and attaching the
      ## sign of the leg length.
      Omega:=function(module,n) local x;
        return H.operations.New(module,List([1..n],x->(-1)^(x)),
                         List([1..n],x->H.operations.Hook(n,x)));
      end,
    
      New:=function(m, c, p)
        return rec(H:=H,module:=m,coeffs:=c,parts:=p,
                   operations:=H.operations.(m));
      end,
    
      ## Takes two lists, one containing coefficients and the other the
      ## corresponding partitions, and orders them lexicogrphcailly collecting
      ## like terms on the way. We use a variation on quicksort which is
      ## induced by the lexicographic order (if parts contains partitions of
      ## different integers this can lead to an error - which we don't catch).
      Collect:=function(module, coeffs, parts)
        local newx, i, Place, Unplace, places;
    
        ## inserting parts[i] into places. if parts[i]=[p1,p2,...] then
        ## we insert it into places at places[p1][p2][...] stopping
        ## at the fist empty position (say places[p1], or places[p1][p2]
        ## etc.). Here we are trying to position parts[i] at
        ## places[p1]...[pd]. Actually, we just insert i rather than
        ## parts[i] and if we find that parts[i]=parts[j] for some j then
        ## we set coeffs[i]:=coeffs[i]+coeffs[j] and don't place j.
        ## Surprisingly, these routines are significantly faster than
        ## usin Set() to sort and collect the terms.
        Place:=function(i, places, d) local t;
          if IsInt(places) then
            t:=places;
            if parts[t]=parts[i] then coeffs[t]:=coeffs[t]+coeffs[i];
            else
              places:=[];
              places[parts[t][d]]:=t;
              if parts[i][d]<>parts[t][d] then places[parts[i][d]]:=i;
              else places:=Place(i, places, d);
              fi;
            fi;
          elif places=[] or not IsBound(places[parts[i][d]]) then  
            # must be a list
            places[parts[i][d]]:=i;
          else places[parts[i][d]]:=Place(i, places[parts[i][d]], d+1);
          fi;
          return places;
        end;
        Unplace:=function(places) local p, newp, np;
          newp:=[[],[]];
          for p in places do
            if IsInt(p) then if coeffs[p]<>0*coeffs[p] then
              Add(newp[1], coeffs[p]);
              Add(newp[2], Copy(parts[p])); fi;
            else
              np:=Unplace(p);
              Append(newp[1], np[1]);
              Append(newp[2], np[2]);
            fi;
          od;
          return newp;
        end;
        if IsBound(H.PartitionLength) then # strip out long partitions
	        places:=Filtered([1..Length(parts)],
                            i->Length(parts[i])<=H.PartitionLength);
	        parts:=parts{places};
	        coeffs:=coeffs{places};
        fi;
        if parts=[] then return H.operations.New(module,[0],[[]]);
        else 
          places:=[];
          for i in [1..Length(parts)] do places:=Place(i, places, 1); od;
          newx:=Unplace(places);
          if newx=[[],[]] then return H.operations.New(module,[0],[[]]);
          else return H.operations.New(module,newx[1],newx[2]);
          fi;
        fi;
     end  ## H.operations.Collect
    
    )); ## H.operations(); although we are not quite finished with it yet:
  ## Next, we run through the arguments to Specht(), and add the
  ## appropriate things to H (this is more flexible than is 
  ## explained in the manual).
  if Length(arg)=1 and IsList(arg[1]) then arg:=arg[1]; fi; ## see Schur()
  for a in arg do
    if IsInt(a) then
      if IsBound(H.e) then  ## set H.p
        if IsPrime(a) then H.p:=a;
        else Error("Specht(<e>,<p>), <p> must be a prime number");
        fi;
        if H.e=H.p then 
          H.HeckeRing:=Concatenation("p",String(H.p),"sym");
          ## return the exponent of the maximum power of p dividing x
          H.valuation:=function(x) local i;
            i:=0;
            while x mod H.p=0 do
              i:=i+1;
              x:=x/H.p;
            od;
            return i;
          end;
        elif not IsBound(H.valuation) then
          ## return the exponent of the maximum power of p that
          ## divides e^x-1.
          H.valuation:=function(x) local i;
            if x mod H.e<>0 then return 0;
            else
              i:=0;
	      x:=x/H.e;
              while x mod H.p=0 do
                i:=i+1;
                x:=x/H.p;
              od;
              return 1+H.p^i;
            fi;
          end;
        fi;
      else                 ## Set H.e
        H.e:=a;
        H.info.Library:=ConcatenationString(SPECHT.Library,"e",String(H.e),"/");
      fi;
    elif IsFunc(a) then
      H.valuation:=a;
      if not IsBound(H.HeckeRing) then H.HeckeRing:="unknown"; fi;
      if not IsBound(H.e) then
        H.e:=1; 
        while H.valuation(H.e)=0 do H.e:=H.e+1; od;
      fi;
    elif IsString(a) then H.HeckeRing:=a;
    else Error("usage: Specht(<e> [ [,<p>] [,<valuation>] [,<ring>] ])");
    fi;
  od;
  if not IsBound(H.e) then
    Error("usage: Specht(<e> [ [,<p>] [,<valuation>] [,<ring>] ])");
  fi;
  if not IsBound(SPECHT.Cache[H.e]) then SPECHT.Cache[H.e]:=[]; fi;
  if not IsBound(H.p) then H.p:=0; fi;
  if not IsBound(H.valuation) then
    if H.e=0 then H.valuation:=i->i;
    else
      H.valuation:=function(x)
        if x mod H.e = 0 then return 1;
        else return 0;
        fi;
      end;
    fi;
  fi;
  if not IsBound(H.HeckeRing) then
    H.HeckeRing:=Concatenation("e",String(H.e), "p",String(H.p));
  fi;

  ## The next five functions are for restricting and inducing Specht
  ## modules. They all assume that their arguments are indeed Specht
  ## modules; conversations are done in H.operations.X.Y() as necessary.

  ## r-induction: on Specht modules:
  InducedModule:=function(a, e, r) local ind, x, i, j, np;
    ind:=[[],[]];
    for x in [1..Length(a.parts)] do
      for i in [1..Length(a.parts[x])] do
        if (a.parts[x][i] + 1 - i) mod e=r then
          if i=1 or a.parts[x][i-1] > a.parts[x][i] then
            np:=Copy(a.parts[x]);
            np[i]:=np[i]+1;
            Add(ind[1], a.coeffs[x]);
            Add(ind[2], np);
          fi;
        fi;
      od;
      if ( -Length(a.parts[x]) mod e)=r then
        np:=Copy(a.parts[x]); Add(np, 1);
        Add(ind[1],a.coeffs[x]);
        Add(ind[2],np);
      fi;
    od;
    if ind=[ [],[] ] then return H.operations.New("S",[0],[[]]);
    else return H.operations.Collect("S", ind[1], ind[2]);
    fi;
  end;

  # add n nodes of residue r to the partition y from the i-th row down
  sinduced:=function(y, n, e, r, i) local ny, j, z;
    ny:=[];
    for j in [i..Length(y)-n+1] do
      if r=(y[j] - j + 1) mod e then 
        if j=1 or y[j] < y[j-1] then
          z:=Copy(y);
          z[j]:=z[j] + 1; # only one node of residue r can be added
          if n=1 then Add(ny, z);   # no more nodes to add
          else Append(ny, sinduced(z, n-1, e, r, j+1));
          fi;
        fi;
      fi;
    od;
    return ny;
  end;

  ## String-induction: add s r's from each partition in x (ignoring
  ## multiplicities). Does both standard and q-induction.

  ## We look at the size of x.module to decide whether we want to use
  ## ordinary indcution or q-induction (in the Fock space). We could
  ## write H.operations.X.SInduce to as to make this choice for us, or
  ## do q-induction always, setting v=1 afterwards, but this seems the
  ## better choice.
  SInducedModule:=function(x, e, s, r) local coeffs, parts, y, z;
    if s=0 then return x.H.operations.New(x.module,[1],[[]]); fi;
    coeffs:=[]; parts:=[];
    for y in [1..Length(x.parts)] do
      Append(parts,sinduced(x.parts[y], s, e, r, 1));
      Append(coeffs,List([1..Length(parts)-Length(coeffs)],r->x.coeffs[y]));
      if r=( -Length(x.parts[y]) mod e) then # add a node to the bottom
        z:=Copy(x.parts[y]);
        Add(z,1);
        if s > 1 then                        # need to add some more nodes
          Append(parts,sinduced(z, s-1, e, r, 1));
          Append(coeffs,List([1..Length(parts)-Length(coeffs)],
                             r->x.coeffs[y]));
        else Add(coeffs, x.coeffs[y]); Add(parts, z);
        fi;
      fi;
    od;

    if coeffs=[] then return H.operations.New(x.module,[0],[[]]);
    else return H.operations.Collect(x.module, coeffs, parts);
    fi;
  end;  # SInduce

  ## r-restriction
  RestrictedModule:=function(a, e, r) local ind, x, i, j, np;
    ind:=[[],[]];
    for x in [1..Length(a.parts)] do
      for i in [1..Length(a.parts[x])] do
        if (a.parts[x][i] - i) mod e=r then
          np:=Copy(a.parts[x]);
          if i=Length(a.parts[x]) or np[i] > np[i+1] then
            np[i]:=np[i] - 1;
            if np[i]=0 then Unbind(np[i]); fi;
            Add(ind[1], a.coeffs[x]); Add(ind[2], np);
          fi;
        fi;
      od;
    od;
    if ind=[ [],[] ] then return H.operations.New("S",[0],[[]]);
    else return H.operations.Collect("S", ind[1], ind[2]);
    fi;
  end;

  ## remove n nodes from y from the ith row down
  srestricted:=function(y, n, e, r, i) local ny, j, z;
    ny:=[];
    for j in [i..Length(y)-n+1] do
      if r=(y[j] - j) mod e then
        if j=Length(y) or y[j] > y[j+1] then
          z:=Copy(y);
          z[j]:=z[j] - 1;
          if z[j]=0 then   # n must be 1
            Unbind(z[j]);
            Add(ny, z);
          elif n=1 then Add(ny, z); # no mode nodes to remove
          else Append(ny, srestricted(z, n-1, e, r, j+1));
          fi;
        fi;
      fi;
    od;
    return ny;
  end;

  ## string-restriction: remove m r's from each partition in x
  SRestrictedModule:=function(x,e,s,r) local e, coeffs, parts, y, i;
    coeffs:=[]; parts:=[];
    for y in [1..Length(x.parts)] do
      Append(parts, srestricted(x.parts[y], s, e, r, 1));
      Append(coeffs,List([1..Length(parts)-Length(coeffs)],i->x.coeffs[y]));
    od;
    if parts=[] then return x.H.operations.New("S",[0],[[]]);
    else return x.H.operations.Collect("S", coeffs, parts);
    fi;
  end;  # SRestrict

  InnerProduct:=function(a,b) local pr, x, y;
    if a=0*a or b=0*b then return 0;
    elif a.module<>b.module then
      a:=a.operations.S(a,true);
      b:=b.operations.S(b,true);
    fi;

    pr:=0; x:=1; y:=1;  # use the fact that a.parts and b.parts are ordered
    while x <=Length(a.parts) and y <=Length(b.parts) do
      if a.parts[x]=b.parts[y] then
        pr:=pr + a.coeffs[x]*b.coeffs[y];
        x:=x + 1; y:=y + 1;
      elif a.parts[x]<b.parts[y] then x:=x + 1;
      else y:=y + 1;
      fi;
    od;
    return pr;
  end;  # InnerProduct

  #F Returns the Coefficient of p in x
  Coefficient:=function(arg) local x,p,pos;
    x:=arg[1];p:=arg{[2..Length(arg)]};
    pos:=Position(x.parts, p);
    if pos=false then return 0;
    else return x.coeffs[pos];
    fi;
  end;
    
  #F Returns true if all coefficients are non-negative
  PositiveCoefficients:=function(x) local c;
    return ForAll(x.coeffs, c->c>=0);
  end;

  #F Returns true if all coefficients are positive
  IntegralCoefficients:=function(x) local c;
    return ForAll(x.coeffs, c->IsInt(c));
  end;

  ## addition for partitions. Note that it must keep the records
  ## in lexigraphical order and make sure that no partitions appear
  ## with coefficient Zero.
  AddModules:=function(a,b) local i, j, ab, x;
    if a=false or b=false then return false;
    elif a=0*a then return b; 
    elif b=0*b then return a; 
    elif a.H<>b.H then 
      Error("modules belong to different Grothendieck rings");
    fi;
    if a.module<>b.module then # only convert to Specht modules if different
      if Length(a.module)<>Length(b.module) then
        Error("AddModule(<a>,<b>): can only add modules of same type.");
      fi;
      a:=a.operations.S(a,false);  
      b:=b.operations.S(b,false); 
      if a=false or b=false then return false;fi;
    fi;
    ## profiling shows _convincingly_ that the method used below to add
    ## a and b is faster than using SortParallel or H.operations.Collect.
    ab:=[[],[]];
    i:=1; j:=1;
    while i <=Length(a.parts) and j <=Length(b.parts) do
      if a.parts[i]=b.parts[j] then
        x:=a.coeffs[i]+b.coeffs[j];
        if x<>0*x then 
          Add(ab[1],x); 
          Add(ab[2], a.parts[i]); 
        fi;
        i:=i+1; j:=j+1;
      elif a.parts[i] < b.parts[j] then
        if a.coeffs[i]<>0*a.coeffs[i] then 
          Add(ab[1], a.coeffs[i]);
          Add(ab[2], a.parts[i]); 
        fi;
        i:=i+1;
      else
        if b.coeffs[j]<>0*b.coeffs[j] then 
          Add(ab[1], b.coeffs[j]);
          Add(ab[2], b.parts[j]); 
        fi;
        j:=j+1;
      fi;
    od;
    if i <=Length(a.parts) then 
      Append(ab[1], a.coeffs{[i..Length(a.coeffs)]});
      Append(ab[2], a.parts{[i..Length(a.parts)]});
    elif j <=Length(b.parts) then 
      Append(ab[1], b.coeffs{[j..Length(b.coeffs)]});
      Append(ab[2], b.parts{[j..Length(b.parts)]});
    fi;
    if ab=[[],[]] then ab:=[ [0],[[]] ]; fi;
    return H.operations.New(a.module, ab[1], ab[2]);
  end;  # AddModules
    
  MultiplyModules:=function(a,b) local x, y, ab, abcoeff, xy, z;
    if a=false or b=false then return false;
    elif a=0 or b=0 then return H.operations.New(b.module,[0],[[]]);
    elif not ( IsRec(a) and IsBound(a.H) ) then
      return H.operations.New(b.module, a*b.coeffs, b.parts);
    elif not ( IsRec(b) and IsBound(b.H) ) then
      return H.operations.New(a.module, b*a.coeffs, a.parts);
    elif a.coeffs=[0] or b.coeffs=[0] then 
      return H.operations.New(b.module,[0],[[]]);
    elif a.H<>b.H then 
      Error("modules belong to different Grothendieck rings");
    fi;
    a:=a.operations.S(a,false);
    ab:=[[],[]];
    for x in [1..Length(a.parts)] do
      for y in [1..Length(b.parts)] do
        abcoeff:=a.coeffs[x]*b.coeffs[y];
        if abcoeff<>0*abcoeff then
          z:=LittlewoodRichardsonRule(a.parts[x], b.parts[y]);
          Append(ab[1], List(z, xy->abcoeff));
          Append(ab[2], z);
        fi;
      od;
    od;
    if ab=[] then return H.operations.New(a.module,[0],[[]]);
    else return H.operations.Collect(b.module, ab[1], ab[2]);
    fi;
  end;  # MultiplyModules
    
  ## The Print function for the "module" elements.
  PrintModule:=function(module, a) 
    local r, x, len, n, star, parts, coeffs;

    if IsBound(a.name) then Print(a.name,"\n");
    else
      coeffs:=Copy(a.coeffs);
      parts:=Copy(a.parts);
      
      # if H.PartitionLength is bound then only print those
      # partitions which are short enough.
      if IsBound(H.PartitionLength) then
        n:=Filtered([1..Length(parts)], 
                        r->Length(parts[r])<=H.PartitionLength);
        coeffs:=coeffs{n};
        parts:=parts{n};
        
        if parts=[] then Print("0"); return; fi;
      fi;

      if a.H.Ordering<>Lexicographic then
        SortParallel(parts, coeffs, a.H.Ordering);
      fi;

      if Length(parts) > 1 then
        n:=Sum(parts[1]);
        if not ForAll(parts, x->Sum(x)=n) then
          Error(module,"(x), <x> contains mixed partitions\n\n");
        fi;
      fi;
      if SpechtPrintFn=LabelPartition then star:=""; ## pretty printing
      else star:="*";
      fi;

      if Length(parts)=0 or coeffs=[0] then
        coeffs:=[ 0 ]; parts:=[ [] ];
        Print("0",star,module,"()");
      else
        for x in [Length(parts),Length(parts)-1..1] do
          if IsPolynomial(coeffs[x]) then
            if Length(coeffs[x].coefficients)=1 then
              if coeffs[x].valuation=0 then
                if coeffs[x].coefficients=[-1] then Print(" - ");
                elif coeffs[x].coefficients[1]<0 then 
                  Print(coeffs[x].coefficients[1],star);
                else
                  if x<Length(parts) then Print(" + "); fi;
                  if coeffs[x].coefficients<>[1] then
                    Print(coeffs[x].coefficients[1],star);
                  fi;
                fi;
              else
                if x<Length(parts) and coeffs[x].coefficients[1]>0 then 
                  Print(" + ");
                fi;
                Print(coeffs[x],star);
              fi;
            elif coeffs[x].coefficients[Length(coeffs[x].coefficients)]<0 
            then Print(" - (",-coeffs[x],")", star);
            else
              if x<Length(parts) then Print(" + "); fi;
              Print("(",coeffs[x],")", star);
            fi;
          else
            if coeffs[x]=-1 then Print(" - ");
            elif coeffs[x]<0 then Print(coeffs[x],star);
            else
              if x<Length(parts) then Print(" + "); fi;
              if coeffs[x]<>1 then Print(coeffs[x],star); fi;
            fi;
          fi;
          Print(module,"(",SpechtPrintFn(parts[x]),")");
        od;
      fi;
    fi;
  end;  # PrintModule
    
  H.operations.S:=rec(
    \+:=function(a,b) return AddModules(a,b); end,
    \-:=function(a,b) 
      if a=false or b=false then return false;
      else 
        b:=Copy(b); b.coeffs:=-b.coeffs;
        return AddModules(a,b); 
      fi;
    end,
    \*:=function(a,b) return MultiplyModules(a,b); end,
    \/:=function(b,n) local x;
      if n=0 then Error("can't divide by 0!\n");
      else return b.H.operations.New(b.module, b.coeffs/n, b.parts);
      fi;
    end,
    
    Print:=function(x)
      if x.H.IsSpecht then PrintModule("S",x);
      else PrintModule("W",x);
      fi;
    end,
    
    Coefficient:=Coefficient,
    PositiveCoefficients:=PositiveCoefficients,
    IntegralCoefficients:=IntegralCoefficients,
    InnerProduct:=InnerProduct,
    
    ## Induction and restriction; for S()
    InducedModule:=function(x, list) local r;
      if x=false or x=0*x then return x;
      elif list=[] then return InducedModule(x,1,0);
      elif H.e=0 then
        Error("Induce, r-induction is not defined when e=0.");
      elif ForAny(list,r-> r>=H.e or r<0) then
        Error("Induce, r-induction is defined only when 0<=r<e.\n");
      else 
        for r in list do
          x:=InducedModule(x,H.e,r);
        od;
        return x;
      fi;
    end,

    RestrictedModule:=function(x, list) local r;
      if x=false or x=0*x then return x;
      elif list=[] then return RestrictedModule(x,1,0);
      elif H.e=0 then
        Error("Restrict, r-restriction is not defined when e=0.");
     elif ForAny(list,r-> r>=H.e or r<0) then
        Error("Restrict, r-restriction is defined only when 0<=r<e.\n");
      else 
        for r in list do
          x:=RestrictedModule(x,H.e,r);
        od;
        return x;
      fi;
    end,

    SInducedModule:=function(x, list) local r;
      if x=false or x=0*x then return x;
      elif Length(list)=1 then
        list:=list[1];
        if list=0 then return H.operations.New("Sq",[1],[[]]); fi;
        while list > 0 do
          x:=SInducedModule(x,1,1,0);
          list:=list-1;
        od;
        return x;
      elif H.e=0 then
        Error("SInduce, r-induction is not defined when e=0.");
      elif list[2]>H.e or list[2]<0 then
        Error("SInduce, r-induction is defined only when 0<=r<e.\n");
      else return SInducedModule(x, H.e, list[1], list[2]);
      fi;
    end,

    SRestrictedModule:=function(x, list) local r;
      if x=false or x=0*x then return x;
      elif Length(list)=1 then
        list:=list[1];
        if list=0 then return H.operations.New("Sq",[1],[[]]); fi;
        while list > 0 do
          x:=SRestrictedModule(x,1,1,0);
          list:=list-1;
        od;
        return x;
      elif H.e=0 then
        Error("SRestrict, r-restriction is not defined when e=0.");
      elif list[2]>H.e or list[2]<0 then
        Error("SRestrict, r-restriction is defined only when 0<=r<e.\n");
      else return SRestrictedModule(x, H.e, list[1], list[2]);
      fi;
    end,
      
    ## Finally the conversion functions S(), P() and D(). All take
    ## a linear combination of Specht modules and return corresponding
    ## linear combinations of Specht, indecomposables, and simples resp.
    ## If they have a problem they return false and print an error
    ## message unless silent=true.
    
    S:=function(x, silent) return x; end,             # S() -> S()
    
    ## Here I only allow for linear combinations of projectives which 
    ## have non-negative coefficients; the reason for this is that I
    ## can't see how to do it otherwise. The problem is that in the
    ## Grothendieck ring there are many ways to write a given linear
    ## combination of Specht modules (or PIMs).
    P:=function(x, silent) local proj;                # S() -> P()
      if x=false or x=0*x then return x; 
      elif x.parts=[[]] then return H.operations.New("P",x.coeffs[1],[]);
      fi;
    
      proj:=H.operations.New("P",[0],[[]]);
      while x<>false and x<>0*x and 
      ( not H.IsSpecht or IsERegular(H.e,x.parts[Length(x.parts)]) ) do
        proj:=proj+H.operations.New("P",[x.coeffs[Length(x.parts)]], 
                                        [x.parts[Length(x.parts)]]);
        x:=x+H.operations.P.S(
                  H.operations.New("P",[-x.coeffs[Length(x.parts)]],
                                       [x.parts[Length(x.parts)]]),true);
      od;
      if x=false or x<>0*x then 
        if not silent then 
          Print("# P(<x>), unable to rewrite <x> as a sum of projectives\n");
        fi;
      else return proj;
      fi;
      return false;
    end, 
    
    D:=function(x,silent) local y, d, simples, r, c;  # S() -> D()
      if x=false or x=0*x then return x;
      elif x.parts=[[]] then return H.operations.New("D",[x.coeffs[1]],[[]]);
      fi;
      
      d:=H.operations.KnownDecompositionMatrix(Sum(x.parts[1]));
      if d<>false then 
        y:=d.operations.S.D(d,x); 
        if y<>false then return y; fi;
      fi;
    
      ## since that didn't work, we use the LLT algorithm when IsBound(H.Pq)
      if IsBound(H.Pq) then
        return Sum([1..Length(x.parts)], 
                   r->x.coeffs[r]*Specialized(Sq(x.parts[r]))); 
      fi;

      # next, see if we can calculate the answer.
      d:=Concatenation(H.HeckeRing,"D");
      # finally, we can hope that only partitions of e-weight<2 appear in x
      r:=1; simples:=H.operations.New("D",[0],[[]]); 
      while simples<>false and r <= Length(x.parts) do
        if IsSimpleModule(x.H, x.parts[r]) then
          simples:=simples+H.operations.New("D",[x.coeffs[r]], [x.parts[r]]);
        elif IsERegular(x.H.e,x.parts[r]) and EWeight(x.H.e,x.parts[r])=1 
        then
          y:=H.operations.New("S",[1],[ECore(x.H.e,x.parts[r])])
                             * H.operations.Omega("S",x.H.e);
          c:=Position(y.parts,x.parts[r]); ## >1 since not IsSimpleModule
          simples:=simples
                    +H.operations.New("D",[1,1]*x.coeffs[r],[y.parts[c],y.parts[c-1]]);
        elif IsBound(x.operations.(d)) then 
          simples:=simples+x.coeffs[r]*x.operations.(d)(x.parts[r]);
        else simples:=false;
        fi;
        r:=r+1;
      od;
      if simples=false and not silent then
        Print("# D(<x>), unable to rewrite <x> as a sum of simples\n");
        return false;
      else return simples;
      fi;
    end
  );   ## H.operations.S
    
  H.operations.P:=rec(
    \+:=function(a,b) return AddModules(a,b); end,
    \-:=function(a,b) 
      if a=false or b=false then return false;
      else 
        b:=Copy(b); b.coeffs:=-b.coeffs;
        return AddModules(a,b);
      fi;
    end,
    \*:=function(a,b) local x, nx;
      if not( IsRec(a) and IsBound(a.module) ) then 
        return H.operations.New(b.module,a*b.coeffs,b.parts);
      elif a.module=b.module then 
        x:=MultiplyModules(a.operations.S(a,false),b.operations.S(b,false));
        nx:=x.operations.P(x,true);
        if nx<>false then return nx; else return x; fi;
      else   
        return MultiplyModules(a.operations.S(a,false),
                               b.operations.S(b,false));
      fi;
    end,
    \/:=H.operations.S.\/,

    Print:=function(x) PrintModule("P",x); end,

    Coefficient:=Coefficient,
    PositiveCoefficients:=PositiveCoefficients,
    IntegralCoefficients:=IntegralCoefficients,
    InnerProduct:=InnerProduct,
  
    InducedModule:=function(x,list) local nx;
      x:=H.operations.S.InducedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.P(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    SInducedModule:=function(x,list) local nx;
      x:=H.operations.S.SInducedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.P(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    RestrictedModule:=function(x,list) local nx;
      x:=H.operations.S.RestrictedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.P(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    SRestrictedModule:=function(x,list) local nx;
      x:=H.operations.S.SRestrictedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.P(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    ## Next the functions S(), P(), and D(). The P->S functions are
    ## quite involved.
   
    #F Writes x, which is a sum of indecomposables, as a sum of S(nu)'s if 
    ## possible. We first check to see if the decomposition matrix for x is
    ## stored somewhere, and if not we try to calculate what we need. If we
    ## can't do this we return false.
    S:=function(x,silent) local y, c, d, mu, specht;     # P() -> S()
      if x=false or x=0*x then return x; 
      elif x.parts=[[]] then return H.operations.New("S",[x.coeffs[1]],[[]]);
      fi;
      d:=H.operations.KnownDecompositionMatrix(Sum(x.parts[1]));
      if d<>false then 
        y:=d.operations.P.S(d,x); 
        if y<>false then return y; fi;
      fi;
    
      ## since that didn't work, we use the LLT algorithm when 
      ## IsBound(H.Pq)
      if IsBound(H.Pq) then
        if not H.IsSpecht or ForAll(x.parts, c->IsERegular(H.e,c)) then
          return Sum([1..Length(x.parts)],c->
                   x.coeffs[c]*Specialized(Pq(x.parts[c])));
        fi;
      fi;
      
      d:=Concatenation(H.HeckeRing,"S");
      mu:=1; specht:=H.operations.New("S",[0],[[]]);
      while specht<>false and mu<=Length(x.parts) do
        if IsSimpleModule(H,ConjugatePartition(x.parts[mu])) then 
          specht:=specht+H.operations.New("S",[1],[x.parts[mu]]);
        elif EWeight(x.H.e,x.parts[mu])=1 then ## wrap e-hooks onto c 
          c:=H.operations.New("S",[1],[ECore(x.H.e, x.parts[mu])])
                   * H.operations.Omega("S",x.H.e); 
          y:=Position(c.parts, x.parts[mu]);
          if y=1 then specht:=specht+H.operations.New("S",[1],[c.parts[y]]);
	  elif IsInt(y) then specht:=specht+H.operations.New("S",[1,1],
                                          [c.parts[y-1],c.parts[y]]);
	  else specht:=false;
          fi;
        elif IsBound(H.operations.P.(d)) then 
          specht:=specht+x.H.operations.P.(d)(x.parts[mu]); 
        else specht:=false; 
        fi;
        mu:=mu+1;
      od;
      if specht<>false then return specht; 
      elif not silent then 
        Print("# P(<x>), unable to rewrite <x> as a sum of projectives\n");
      fi;
      return false;
    end, 
     
    P:=function(x,silent) return x; end,                   # P() -> P()

    D:=function(x,silent)                                  # P() -> D()
      x:=x.operations.S(x,silent);
      if x=false then return x;
      else return x.operations.D(x,silent); 
      fi;
    end
  );    ## H.operations.P()
    
  H.operations.D:=rec(
    \+:=function(a,b) return AddModules(a,b); end,
    \-:=function(a,b) 
      if a=false or b=false then return false;
      else 
        b:=Copy(b); b.coeffs:=-b.coeffs;
        return AddModules(a,b);

      fi;
    end,
    \*:=function(a,b) local x, nx;
      if not( IsRec(a) and IsBound(a.module) ) then
        return H.operations.New(b.module,a*b.coeffs,b.parts);
      elif a.module=b.module then
        x:=MultiplyModules(a.operations.S(a,false),b.operations.S(b,false));
        nx:=x.operations.D(x,true);
        if nx<>false then return nx; else return x; fi;
      else
        return MultiplyModules(a.operations.S(a,false),
                               b.operations.S(b,false));
      fi;
    end,
    \/:=H.operations.P.\/,

    Print:=function(x)
      if x.H.IsSpecht then PrintModule("D",x);
      else PrintModule("F",x);
      fi;
    end,
    
    Coefficient:=Coefficient,
    PositiveCoefficients:=PositiveCoefficients,
    IntegralCoefficients:=IntegralCoefficients,
    InnerProduct:=InnerProduct,
    
    InducedModule:=function(x,list) local nx;
      x:=H.operations.S.InducedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.D(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    RestrictedModule:=function(x,list) local nx;
      x:=H.operations.S.RestrictedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.D(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    SInducedModule:=function(x,list) local nx;
      x:=H.operations.S.SInducedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.D(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    SRestrictedModule:=function(x,list) local x,nx;
      x:=H.operations.S.SRestrictedModule(x.operations.S(x,false),list);
      if x=false or x=0*x then return x; fi;
      nx:=x.operations.D(x,false);
      if nx<>false then return nx; else return x; fi;
    end,
    
    ## Next the functions S(), P(), and D().
    
    #F Writes D(mu) as a sum of S(nu)'s if possible. We first check to see
    ## if the decomposition matrix for Sum(mu) is stored in the library, and
    ## then try to calculate it directly. If we are unable to do this either
    ## we return false.
    S:=function(x,silent) local c, d, y, a;                # D() -> S()
      if x=false or x=0*x then return x;
      elif x.parts=[[]] then return H.operations.New("S",[x.coeffs[1]],[[]]);
      fi;
    
      ## look for the decomposition matrix
      d:=H.operations.KnownDecompositionMatrix(Sum(x.parts[1]));
      if d<>false then 
        y:=d.operations.D.S(d,x); 
        if y<>false then return y; fi;
      fi;
    
      ## since that didn't work, we use the LLT algorithm when IsBound(H.Pq)
      if IsBound(H.Pq) and H.IsSpecht then
        return Sum([1..Length(x.parts)],
                     c->x.coeffs[c]*Specialized(Dq(x.parts[c])));
      fi;
    
      ## Next, see if we can calculate it.
      d:=Concatenation(H.HeckeRing, "S");
      if IsBound(H.operations.D.(d)) then 
        return H.operations.D.(d)(x,silent); 
      fi;
    
      ## Finally, hope only e-weights<2 are involved.
      c:=1; d:=true; y:=H.operations.New("S",[0],[[]]);
      while d and c<=Length(x.parts) do
        if IsSimpleModule(H, x.parts[c]) then 
          y:=y+H.operations.New("S",[x.coeffs[c]],[x.parts[c]]);
        elif IsERegular(x.H.e, x.parts[c]) and EWeight(x.H.e,x.parts[c])=1 
        then ## wrap e-hooks onto c
          a:=H.operations.New("S",[1],[ECore(x.H.e,x.parts[c])])
               * H.operations.Omega("S",x.H.e); 
          a.parts:=a.parts{[1..Position(a.parts, x.parts[c])]};
          a.coeffs:=a.coeffs{[1..Length(a.parts)]}*(-1)^(1+Length(a.parts));
          y:=y+a;
        else d:=false;
        fi;
        c:=c+1;
      od;
      if d<>false then return y; 
      elif not silent then
        Print("# Unable to calculate D(mu)\n");
      fi;
      return false;
    end,  # H.operations.D.S()
                                   # D() -> P()
    P:=function(x,silent) 
      x:=x.operations.S(x,silent);
      if x=false then return x;
      else return x.operations.P(x,silent);
      fi;
    end,
    
    D:=function(x,silent) return x; end                     # D() -> D()
  );  ## H.operations.D

  #################################################################
  ## Next, the functions for dealing with decomposition matrices ##
  #################################################################

  ## Finally, we keep a copy of SpechtDirectory in H so that we have a
  ## chance of finding new decomposition matrices when it changes.
  H.info.SpechtDirectory:=SpechtDirectory;

  ## This record will hold any decomposition matrices which Specht()
  ## (or rather its derivatives) read in. This used to be a public 
  ## record; it is now private because q-Schur algebra matrices and
  ## Hecke algebra matrices might need to coexist. 

  ## Keep ourselves honest when inducing decomposition matrices
  BUG:=function(arg) local a;
     PrintTo("*errout*",
             "\n ***********************************************************************",
             "\n *** Either the decomposition matrix that you have entered is incorrect",
             "\n *** or you may have uncovered a bug in SPECHT's function ",
             "\n ***          ", arg[2],"() - report #", arg[1], 
	     "\n *** If you are confident that your decomposition matrix is correct then", 
	     "\n *** please e-mail all possible details to ", SPECHT.Email,
             "\n ***********************************************************************",
	     "\n"
     );
     for a in arg{[3..Length(arg)]} do
       PrintTo("*errout*", a);
     od;
     Error("\n");
  end;

  ## Pretty printing (and TeXing) of a (decomposition) matrix.
  ## d=decomposition matrix record, TeX=true of false
  ## Actually, d can be any record having .d=matrix, .labels=[strings],
  ## row, and cols components (this is also used by KappaMatrix for example).
  ##   tex=0 normal printing
  ##   tex=1 LaTeX output
  PrintDecompositionMatrix:=function(d, tex)
    local rows, cols, r, c, col, len, endBit, sep, M, label, rowlabel,
          spacestr, dotstr, PrintFn, x, i,j,texstr;
  
    ## catch the empty case
    if d.d=[] then Print(d.d); return; fi;
       
    ## might have to fix up the ordering before printing d
    rows:=Copy(d.rows); 
    cols:=Copy(d.cols); 
    if H.Ordering=Lexicographic then
      rows:=rows{[Length(rows),Length(rows)-1..1]};
      cols:=cols{[Length(cols),Length(cols)-1..1]};
    else
      Sort(rows, H.Ordering);
      Sort(cols, H.Ordering);
    fi;
    rows:=List(rows, r->Position(d.rows,r));
    cols:=List(cols, c->Position(d.cols,c));

    rowlabel:=List(d.rows, LabelPartition);

    if tex=1 then # print tex output
      PrintFn:=a->PrintToString(texstr,TeX(a)); 
                    ## PrintFn() allows us to tex() matrix elements (which
                    ## is necessary for crystallized decomposition matices).
      texstr:=SPrint("$$\\begin{array}{l|*{", Length(d.cols)+1,"}{r}}\n");
      sep:="&";
      endBit:=function(i) PrintToString(texstr,"\\\\\n"); end;
    
      ## gangely work around to tex 1^10 properly as 1^{10} etc.
      label:=function(i) local str, bad, l; 
        bad:=Filtered(Collected(d.rows[i]),l->l[2]>9);
        if bad=[] then PrintToString(texstr,rowlabel[i]);
        else # assume no conflicts as 1^10 and 1^101
          str:=Copy(rowlabel[i]); 
          IsString(str);   ## seems to be necessary...
          for l in bad do
            str:=Replace(str,
                   Concatenation(String(l[1]),"^",String(l[2])),
                   Concatenation(String(l[1]),"^{",String(l[2]),"}"));
          od;
          PrintToString(texstr,str);
        fi;
        PrintToString(texstr,"&");
      end;
    else
      PrintFn:=x->Print(String(x,len)); 
      if tex=0 then sep:=" "; else sep:="#";fi;
      endBit:=function(i)if i<>rows[Length(rows)] then Print("\n"); fi; end;

      M:=-Maximum( List(rows, r->Length(rowlabel[r])) );
      label:=function(i) Print(String(rowlabel[i],M),"| ");end;
    
      ## used to be able to print the dimensions at the end of the row.
      # if false then 
      #   endBit:=function(i) Print(" ", String(d.dim[i],-10),"\n");end;
      # fi;
    fi;
    
    ## Find out how wide the columns have to be (very expensive for 
    ## crystallized matrices - also slightly incorrect as String(<poly>) 
    ## returns such wonders as (2)*v rather than 2*v).
    if tex<>2 then
      if tex=1 then len:=0; 
      else 
        len:=1;
        for i in d.d do
          if i.coeffs<>[] then
            len:=Maximum(len,Maximum(List(i.coeffs,j->Length(String(j)))));
          fi;
        od;
      fi;
      spacestr:=String("",len); 
      dotstr:=String(".",len);
    else
      spacestr:="#";
      dotstr:=".";
    fi;
    col:=0;
    for r in rows do
      label(r); 
      if d.rows[r] in d.cols then col:=col+1; fi;
      for c in [1..Length(cols)] do
        if IsBound(d.d[cols[c]]) and r in d.d[cols[c]].parts then
          PrintFn(d.d[cols[c]].coeffs[Position(d.d[cols[c]].parts,r)]);  
          if c<>cols[1] then PrintFn(sep); fi;
        elif c<=col then PrintFn(dotstr); PrintFn(sep);
        else PrintFn(spacestr); PrintFn(sep);
        fi;
      od;
      endBit(r);
    od;
    if tex=1 then 
      PrintToString(texstr,"\\end{array}$$\n"); 
      return texstr;
    fi;
  end;   ## PrintDecompositionMatrix()

  CrystalMatrixOps:=rec(name:="CrystalMatrixOps",
    \*:=function(x,y)
       if x.cols<>y.rows then return false; fi;
       return DecompositionMatrixMatrix(x.H,
                MatrixDecompositionMatrix(x)*MatrixDecompositionMatrix(y),
                x.rows[Length(x.rows)][1]);
    end,

    \=:=function ( x, y )
      return IsDecompositionMatrix(x) and IsDecompositionMatrix(y) and
           x.d=y.d and x.cols=y.cols and x.rows=y.rows;
      end,
    Print  :=function(d) PrintDecompositionMatrix(d,0); end,
    TeX    :=function(d) return PrintDecompositionMatrix(d,1); end,

    #P This function adds the column for Px to the decomposition matrix <d>.
    ## if <checking>=true then Px is checked against its image under the
    ## Mullineux map; if this image is not there then we also insert it
    ## into <d>.
    AddIndecomposable:=function(d,Px,checking) 
      local c, mPx, mx, coeffs, parts, rp, r;
      if Px=false or Px.parts=[] or Px.parts=[[]] or Px=0*Px then return false; fi;
      c:=Position(d.cols, Px.parts[Length(Px.parts)]);
      if checking then
        ## first look to see if <Px> already exists in <d>
        if IsBound(d.d[c]) then ## Px already exists
          SpechtInfo("# AddIndecomposable: overwriting old value of P(",
                TightStringList(Px.parts[Length(Px.parts)]),") in <d>\n");
          d.inverse:=[];      # just in case these were computed
          d.dimensions:=[];
        fi;
        ## now looks at the image of <Px> under Mullineux
        if (Px.H.IsSpecht or IsERegular(Px.H.e,Px.parts[Length(Px.parts)])) 
        and Px.parts[Length(Px.parts)]<>ConjugatePartition(Px.parts[1]) 
	      and not IsBound(Px.H.PartitionLength) then
          mPx:=MullineuxMap(Px);
	        mx:=Position(d.cols,mPx.parts[Length(Px.parts)]);
          if mx<>false then
	          if IsBound(d.d[mx]) and d.P(d,mPx.parts[Length(Px.parts)]) <> mPx then
              SpechtInfo("# AddIndecomposable(<d>,<Px>), WARNING: P(",
                  TightStringList(Px.parts[Length(Px.parts)]), ") and P(",
                  TightStringList(mPx.parts[Length(Px.parts)]),
                  ") in <d> are incompatible\n");
            else
              SpechtInfo("# AddIndecomposable(<d>,<Px>): adding MullineuxMap(<Px>) ",
                  "to <d>\n");
              d.operations.AddIndecomposable(d,mPx,false);
            fi;
          fi;
        fi;
      fi;      # end of check
      coeffs:=[]; parts:=[];
      for r in [1..Length(Px.coeffs)] do
	      rp:=Position(d.rows,Px.parts[r]);
	      if IsInt(rp) then 
	        Add(coeffs, Px.coeffs[r]);
	        Add(parts, rp);
	      fi;
      od;
      SpechtInfo("Adding P(",Px.parts[Length(Px.parts)],")\n");
      d.d[c]:=rec(coeffs:=coeffs, parts:=parts );
    end, # AddIndecomposable

    ## Used by SaveDecomposition matrix to update CrystalMatrices[]
    Store:=function(d) local mu;
      for mu in d.cols do SPECHT.PutCached(H,"Pq",d.P(d,mu));od;
    end,

    ## and now the converstion functions X.Y() writes a sum of X(mu)'s as
    ## a sum of Y(mu)'s.
    S:=rec(
      S:=function(d,x) return x; end,       # S() -> S()

      P:=function(d,x) local nx, r, c, P, S;      # S() -> P()
        if x=false or x=0*x then return x; fi;
        if d.IsDecompositionMatrix then P:="P"; S:="S";
        else P:="Pq"; S:="Sq"; 
        fi;

        nx:=x.H.operations.New(P,[0],[[]]);
        while x<>false and x<>0*x do
          c:=Position(d.cols,x.parts[Length(x.parts)]);
          if c=false or not IsBound(d.d[c]) then return false; fi;
          nx:=nx+x.H.operations.New(P,[x.coeffs[Length(x.parts)]],[d.cols[c]]);
          x:=x+H.operations.New(S,-x.coeffs[Length(x.parts)]*d.d[c].coeffs,
                                List(d.d[c].parts,r->d.rows[r]));
        od;
        return nx;
      end,

      D:=function(d,x) local nx, y, r, rr, c, D, core;   # S() -> D()
        if x=false or x=0*x then return x; fi;
        if d.IsDecompositionMatrix then D:="D"; else D:="Dq"; fi;

        nx:=x.H.operations.New(D,[0],[[]]);
        for y in [1..Length(x.parts)] do
          r:=Position(d.rows, x.parts[y]);
          if r=false then return false; fi;
          core:=ECore(x.H.e,x.parts[y]);
          c:=Length(d.cols);
          while c>0 and d.cols[c]>=x.parts[y] do
            if IsBound(d.d[c]) then
              rr:=Position(d.d[c].parts,r);
              if rr<>false then nx:=nx+H.operations.New(D,
                             [x.coeffs[y]*d.d[c].coeffs[rr]],[d.cols[c]]);
              fi;
            elif ECore(x.H.e,d.cols[c])=core then return false;
            fi;
            c:=c-1;
          od;
        od;
        return nx;
      end
    ),  #S()
    
    P:=rec(
      S:=function(d,x) local S, nx, y, r, c;     # P() -> S()
        if x=false or x=0*x then return x; fi;
        if d.IsDecompositionMatrix then S:="S"; else S:="Sq"; fi;

        nx:=H.operations.New(S,[0],[[]]);
        for y in [1..Length(x.parts)] do
          c:=Position(d.cols,x.parts[y]);
          if c=false or not IsBound(d.d[c]) then return false; fi;
          nx:=nx+H.operations.New(S,x.coeffs[y]*d.d[c].coeffs,
                                  List(d.d[c].parts,r->d.rows[r]));
        od;
        return nx;
      end,

      P:=function(d,x) return x; end,         # P() -> P()

      D:=function(d,x) return d.operations.S.D(d,d.operations.P.S(d,x)); end
    ),  # P()
    
    ## writes D(mu) as a linear combination of S(nu)'s if possible
    Invert:=function(d,mu) local c, S, D, inv, smu, l;
      if d.IsDecompositionMatrix then S:="S"; D:="D";
      else S:="Sq"; D:="Dq";
      fi;

      c:=Position(d.cols,mu);
      if c=false then return false; 
      elif IsBound(d.inverse[c]) then 
        return H.operations.New(S,d.inverse[c].coeffs,
                      List(d.inverse[c].parts,l->d.cols[l]));
      fi;

      inv:=H.operations.New(S,[1],[mu]);
      smu:=H.operations.New(D,[1],[mu])-d.operations.S.D(d,inv);

      while smu<>false and smu<>0*smu do
        inv:=inv+H.operations.New(S,[smu.coeffs[1]],[smu.parts[1]]);
        smu:=smu+d.operations.S.D(d,
                   H.operations.New(S,[-smu.coeffs[1]],[smu.parts[1]]));
      od;
      if smu=false then return false; fi;

      d.inverse[c]:=rec(coeffs:=inv.coeffs,  
                    parts:=List(inv.parts,l->Position(d.cols,l)));
      return inv;
    end,
     
    D:=rec(
      S:=function(d,x) local S, nx, y, c, inv;    # D() -> S()
        if x=false or x=0*x then return x; fi;
        if d.IsDecompositionMatrix then S:="S"; else S:="Sq"; fi;

        nx:=H.operations.New(S,[0],[[]]);
        for y in [1..Length(x.parts)] do
          c:=Position(d.cols,x.parts[y]);
          if c=false then return false; fi; 
          inv:=d.operations.Invert(d,x.parts[y]);
          if inv=false then return inv; 
          else nx:=nx+x.coeffs[y]*inv;
          fi;
        od;
        return nx;
      end,
      P:=function(d,x) 
        return d.operations.S.P(d,d.operations.D.S(d,x)); 
      end,
      D:=function(d,x) return x; end       # D() -> D()
    )  # D()

  ); ## end of crystallized matrix operations

  ## And now the operations for normal decomposition matrices (we still 
  ## have to add a few functions to CrystalMatrixOps).
  DecompositionMatrixOps:=Copy(CrystalMatrixOps);
  DecompositionMatrixOps.name:="DecompositionMatrixOps";

  ## Used by SaveDecompositionMatrix to update DecompositionMatrices[]
  DecompositionMatrixOps.Store:=Ignore;

  #F This function checks to see whether Px is indecomposable using the
  ## decomposition matrix d. The basic idea is to loop through all of the
  ## (e-regular, unless IsSpecht=false) projectives Py in Px such that Px-Py 
  ## has non-negative coefficients and to then apply the q-Schaper theorem 
  ## and induce simples, together with a few other tricks to decide 
  ## whether or not Py splits off. If yes, then we move on; if we can't
  ## decide (or don't know the value of Py), we spit the dummy and return
  ## false, printing a "message from our sponsor" explaining why we failed
  ## if called from outside InducedModule(). If we can count for all
  ## projectives then we return true. Note that in this case the value
  ## of Px may have changed, but we update the original value of px=Px
  ## before leaving (whether we return false or true).
  DecompositionMatrixOps.IsNewIndecomposable:=function(d,px,oldd,Mu) 
    local Px,nu,regs,Py,y,z,a,b,n,mu,m,M,Message;
  
    if px=false and px=0*px then return false; fi;
  
    if Mu=false then Message:=Ignore;
    else Message:=Print;
    fi;

    Px:=px; Py:=true;      # strip PIMs from the top of Px
    while Py<>false and Px<>0*Px do
      Py:=d.P(d,Px.parts[Length(Px.parts)]);
      if Py<>false then Px:=Px-Px.coeffs[Length(Px.parts)]*Py;fi;
      if not PositiveCoefficients(Px) then 
        BUG(11,"IsNewIndecomposable");
      fi;
    od;
    if Px=0*Px then
      Message("# This module is a sum of known indecomposables.\n");
      return false; 
    fi;
    if IntegralCoefficients(Px/Px.coeffs[Length(Px.parts)]) then 
      Px:=Px/Px.coeffs[Length(Px.parts)]; 
   fi;
    
    regs:=Obstructions(d,Px);
    if IsList(Mu) then regs:=Filtered(regs,mu->mu<Mu); fi;
    
    for y in regs do   ## loop through projectives that might split
    
      if H.IsSpecht and not IsBound(H.PartitionLength) 
	and MullineuxMap(H.e,ConjugatePartition(Px.parts[1]))<>Px.parts[Length(Px.parts)] then
        Py:=true;  ## strip any known indecomposables off the bottom of Px
        while Py<>false and Px<>0*Px do
          Py:=d.P(d,ConjugatePartition(Px.parts[1]));
          if Py<>false and IsERegular(Py.H.e,Py.parts[Length(Py.parts)]) then
            Px:=Px-Px.coeffs[1]*MullineuxMap(Py);
          else Py:=false;
          fi;
        od;
        if not PositiveCoefficients(Px) then BUG(8,"IsNewIndecomposable");fi;
      fi;

      m:=0;            ## lower and upper bounds on decomposition the number
      if Px=0*Px then M:=0;
      else M:=Coefficient(Px,y)/Px.coeffs[Length(Px.parts)]; 
      fi;
      if M<>0 then Py:=d.P(d,y); fi;     
      
      if not ( m=M or Px.parts[Length(Px.parts)]>=y ) then
        if Py=false then 
          Message("# The multiplicity of S(", TightStringList(y),
            ") in <Px> is zero; however, S(", TightStringList(y),
            ") is not known\n");
          px.coeffs:=Px.coeffs; px.parts:=Px.parts; return false;
        else 
          Px:=Px-Coefficient(Px,y)*Py; 
          if not PositiveCoefficients(Px) then BUG(1,"IsNewIndecomposable");fi;
          M:=0;
        fi;
      fi;
      if Px<>0*Px and IntegralCoefficients(Px/Px.coeffs[Length(Px.parts)]) then 
        Px:=Px/Px.coeffs[Length(Px.parts)]; 
      fi;
  
      ## remember that Px.coeffs[Length(Px.parts)] could be greater than 1
      if m<>M and (Coefficient(Px,y) mod Px.coeffs[Length(Px.parts)])<>0 then 
        if Py=false then 
          Message("# <Px> is not indecomposable, as at least ", 
            Coefficient(Px,y) mod Px.coeffs[Length(Px.parts)], " copies of P(",
            TightStringList(y), ") split off.\n# However, P(",
            TightStringList(y), ") is not known\n");
          px.coeffs:=Px.coeffs; px.parts:=Px.parts; return false;
        else
          ## this is at least projective; perhaps more still come off though
          Px:=Px-(Coefficient(Px,y) mod Px.coeffs[Length(Px.parts)])*Py; 
          if not PositiveCoefficients(Px) then BUG(2,"IsNewIndecomposable");fi;
          if IntegralCoefficients(Px/Px.coeffs[Length(Px.parts)]) then 
            Px:=Px/Px.coeffs[Length(Px.parts)];
          fi;
        fi;
      fi;
      
      ## At this point the coefficient of Sx in Px divides the coefficient of 
      ## Sy in Px. If Py splits off then it does so in multiples of 
      ## (Px:Sx)=Px.coeffs[Length(Px.parts)].
      if m<>M and (Py=false 
      or PositiveCoefficients(Px-Px.coeffs[Length(Px.parts)]*Py)) then
        ## use the q-Schaper theorem to test whether Sy is contained in Px
        M:=Minimum(M,InnerProduct(Px/Px.coeffs[Length(Px.parts)],
                                  Schaper(H, y)));
        if M=0 then # Schaper says NO!
          ## Px-(Px:Sy)Py is still projective so substract Py if it is 
          ## known. If Py=false (ie. not known), then at least we know 
          ## that Px is not indecomposable, even though we couldn't
          ## calculate Py.
          if Py=false then 
            Message("# The multiplicity of S(", TightStringList(y),
              ") in P(", TightStringList(Px.parts[Length(Px.parts)]), 
             ") is zero;\n#  however, P(", TightStringList(y), 
             ") is not known.\n");
            px.coeffs:=Px.coeffs; px.parts:=Px.parts; return false;
          else 
            Px:=Px-Coefficient(Px,y)*Py; 
            if not PositiveCoefficients(Px) then BUG(3,"IsNewIndecomposable"); 
            fi;
          fi;
        elif Px.coeffs[Length(Px.parts)]<>Coefficient(Px,y) then
          ## We know that (Px:Sy)>=m>0, but perhaps some Py's still split off
  
          m:=1;                         # up until now m was 0
          if m=M then
            if Coefficient(Px,y)<>m*Px.coeffs[Length(Px.parts)] then
              if Py=false then
                Message("# The multiplicity of S(", TightStringList(y),
                    ") in P(",TightStringList(Px.parts[Length(Px.parts)]), 
                    ") is ", m, "however, P(",TightStringList(y),
                    ") is not known.\n");
                px.coeffs:=Px.coeffs; px.parts:=Px.parts;
                return false;
              else
		Print("Px = ", Px,"\nPy = ", Py,"\n");
                Px:=Px-(Coefficient(Px,y)-Px.coeffs[Length(Px.parts)]*m)*Py;
                if not PositiveCoefficients(Px) or m>M then 
                  BUG(4,"IsNewIndecomposable",Px,m,M,d);
                fi;
              fi;
            fi;
          fi;
        else m:=1;   ## up until this point m was 0
        fi;          ## q-Schaper test
  
        if m<>M and SPECHT.DecompositionNumberOK then
           ## see if we can calculate this decomposition number (this uses
           ## row and column removal)
           a:=DecompositionNumber(Px.H, y, Px.parts[Length(Px.parts)]);
           if a<>false then
             if Px.coeffs[Length(Px.parts)]*a=Coefficient(Px,y) then m:=a; M:=a;
             elif Py<>false then
               ## precisely this many Py's come off
               Px:=Px-(Coefficient(Px,y)-Px.coeffs[Length(Px.parts)]*a)*Py; 
               m:=a; M:=a; # upper and lower bounds are equal
               if not PositiveCoefficients(Px) then BUG(5,"IsNewIndecomposable"); fi;
             fi;
           fi;
        fi;
   
        if m<>M and Py=false then ## nothing else we can do
          Message("# The multiplicity of S(", TightStringList(y),
                  ") in P(",TightStringList(Px.parts[Length(Px.parts)]),
                  ") is at least ", m, " and at most ", M, 
                  ";\n# however, P(",TightStringList(y),") is not known.\n");
          px.coeffs:=Px.coeffs; px.parts:=Px.parts; 
          return false;
        fi;

        if m<>M then
          ## Maybe the Mullineux map can lower M...
          M:=Minimum(M,Coefficient(Px,
          Py.parts[Length(Py.parts)]/Px.coeffs[Length(Px.parts)]));
          if Coefficient(Px,y)>M*Px.coeffs[Length(Px.parts)] then
            Px:=Px-(Coefficient(Px,y)-M*Px.coeffs[Length(Px.parts)])*Py;
	    if not PositiveCoefficients(Px) then BUG(6,"IsNewIndecomposable"); fi;
          fi;
          while m<M and not 
          PositiveCoefficients(Px-(Px.coeffs[Length(Px.parts)]*(M-m))*Py) do
            m:=m+1;
          od;
        fi;

        ## finally, we take a look at inducing the simples in oldd
        if m<>M and oldd<>false then
          if not IsBound(oldd.simples) then ## we have to first induce them
            oldd.simples:=[];
            if IsBound(oldd.ind) then       ## defined in InducedModule()
              for mu in [1..Length(oldd.cols)] do
                a:=oldd.operations.Invert(oldd,oldd.cols[mu]);
                if a<>false then
                  for n in [1..H.e] do  ## induce simples of oldd
                    z:=Sum([1..Length(a.coeffs)],b->a.coeffs[b]
                         *oldd.ind[Position(oldd.rows,a.parts[b])][n]);
                    if z<>0*z then Add(oldd.simples,z);fi;
                  od;
                fi;
              od;
            else   ## do everything by hand
              for mu in [1..Length(oldd.cols)] do
                a:=oldd.operations.Invert(oldd,oldd.cols[mu]);
                if a<>false then
                  for n in [0..H.e-1] do
                    z:=InducedModule(a,H.e,n);
                    if z<>0*z then Add(oldd.simples,z); fi;
                  od;
                fi;
              od;
            fi;
          fi;
          mu:=Length(oldd.simples);
          while mu >0 and m<M do
            z:=oldd.simples[mu];
            if y=regs[Length(regs)] 
            or Lexicographic(z.parts[1],Py.parts[Length(Py.parts)]) then 
              a:=InnerProduct(z,Py);
              if a<>0 then 
                b:=InnerProduct(z,Px)/Px.coeffs[Length(Px.parts)];
                m:=Maximum(m,M-Int(b/a));
              fi;
            fi;
            mu:=mu-1;
          od;
          if Coefficient(Px,y)>M*Px.coeffs[Length(Px.parts)] then
            Px:=Px-(Coefficient(Px,y)-M*Px.coeffs[Length(Px.parts)])*Py;
	    if not PositiveCoefficients(Px) then BUG(7,"IsNewIndecomposable"); fi;
          fi;
        fi;
        if m<>M then ## nothing else we can do
          px.coeffs:=Px.coeffs; px.parts:=Px.parts;
          Message("# The multiplicity of S(", TightStringList(y),
              ") in P(",TightStringList(Px.parts[Length(Px.parts)]),
              ") is at least ", m, " and at most ", M, ".\n");
          return false;
        fi;
      fi;
    od;
    if Px=0*Px then
      Message("# This module is a sum of known indecomposables.\n");
      return false;
    elif Px.coeffs[Length(Px.parts)]<>1 then BUG(9,"IsNewIndecomposable");
    else px.coeffs:=Px.coeffs; px.parts:=Px.parts; return true;
    fi; 
  end;  # IsNewIndecomposable

  ## Given a decomposition matrix induce it to find as many columns as
  ## possible of the next higher matrix using simple minded induction.
  ## Not as simple minded as it was originally, as it now tries to use
  ## Schaper's theorem [JM2] to break up troublesome projectives. The
  ## function looks deceptively simple because all of the work is now
  ## done by IsNewIndecomposable().
  ## Usage: InducedDecompositionMatrix(dn)
  ## in the second form new columns are added to d{n+1}.
  DecompositionMatrixOps.Induced:=function(arg)
    local d, newd, mu, nu, Px, Py, n,r;
  
    d:=arg[1];
    if not (IsBound(d.IsDecompositionMatrix) and d.IsDecompositionMatrix=true)
    then Error("InducedDecompositionMatrix(d): ",
                 "<d> must be a decomposition matrix.");
    fi;
  
    n:=Sum(d.rows[1])+1;
    if n>8 then                            ## print dots to let the user
      PrintTo("*stdout*","# Inducing.");   ## know something is happening.
    fi;
     
    nu:=Partitions(n);
    if not d.H.IsSpecht then 
      newd:=H.operations.NewDecompositionMatrix(nu, nu, true);
    else newd:=H.operations.NewDecompositionMatrix(nu,
                ERegularPartitions(H.e,n),true);
    fi;
    ## add any P(mu)'s with EWeight(mu)<=1 or P(mu)=S(mu) <=> S(mu')=D(mu')
    for mu in newd.cols do
      if EWeight(H.e,mu)<=1 then
        newd.operations.AddIndecomposable(newd, 
           H.operations.P.S(H.operations.New("P",[1],[mu]),true),false);
      elif IsSimpleModule(H,ConjugatePartition(mu)) then
        newd.operations.AddIndecomposable(newd, 
           H.operations.New("S",[1],[mu]),false);
      fi;
    od;
        
    ## next we r-induce all of the partitions in d so we can just add 
    ## them up as we need them later.
    ## (note that this InducedModule() is Specht()'s and not the generic one)
    d.ind:=List(d.rows, mu->List([0..d.H.e-1],
              r->InducedModule(H.operations.New("S",[1],[mu]),H.e,r)));
    
    if n<9 then n:=Length(d.cols)+1; fi; ## fudge for user friendliness
  
    for mu in [1..Length(d.cols)] do
      if IsBound(d.d[mu]) then
        for r in [1..d.H.e] do   ## really the e-residues; see ind above
          ## Here we calculate InducedModule(P(mu),H.e,r).
          Px:=Sum([1..Length(d.d[mu].parts)],
                     nu->d.d[mu].coeffs[nu]*d.ind[d.d[mu].parts[nu]][r]); 
          if d.operations.IsNewIndecomposable(newd,Px,d,false) then
            if IsERegular(Px.H.e,Px.parts[Length(Px.parts)]) and
	      not IsBound(d.H.PartitionLength) then   
              # can apply MullineuxMap
              nu:=ConjugatePartition(Px.parts[1]);
              if nu<>MullineuxMap(d.H.e,Px.parts[Length(Px.parts)]) then
                ## wrong image under the Mullineux map
                BUG(10, "Induce", "nu = ", nu, ", Px = ", Px);
              else   ## place the Mullineux image of Px as well
                newd.operations.AddIndecomposable(newd,MullineuxMap(Px),false);
              fi; 
            fi;
            newd.operations.AddIndecomposable(newd,Px,false);
          fi;
        od;
        if mu mod n = 0 then PrintTo("*stdout*",".");fi;
      fi;
    od;
    Unbind(d.ind); Unbind(d.simples); ## maybe we should leave these.
  
    if n>8 then Print("\n"); fi;
    MissingIndecomposables(newd);
    return newd;
  end; # DecompositionMatrixOps.Induce
      
  ## Finally, we can define the creation function for decomposition matrices
  ## (note that NewDM() does not add the partition labels to the decomp.
  ## matrix; this used to be done here but now happens in PrintDM() because
  ## crystallized matrices may never be printed and this operation is
  ## expensive).
  ## **NOTE: we assume when extracting entries from d that d.rows is
  ## ordered lexicographically. If this is not the case then addition
  ## will not work properly.
  H.operations.NewDecompositionMatrix:=function(rows, cols, decompmat) local mu;
    if IsBound(H.PartitionLength) then
      rows:=Filtered(rows, mu->Length(mu)<=H.PartitionLength);
      cols:=Filtered(cols, mu->Length(mu)<=H.PartitionLength);
    fi;
    if decompmat then 
      return rec(d:=[],      # matrix entries
             rows:=rows, # matrix rows
             cols:=cols, # matrix cols
             inverse:=[], #inverse  matrix
             dimensions:=[], pdimensions:=[], ## simple and proj. dim.
             IsDecompositionMatrix:=true,  # not crystallized
             H:=H, operations:=DecompositionMatrixOps,
             P:=function(d,mu)     ## a lazy helper
               return d.operations.P.S(d, d.H.operations.New("P",[1],[mu]));
             end);
    else
      return rec(d:=[],      # matrix entries
             rows:=rows, # matrix rows
             cols:=cols, # matrix cols
             inverse:=[], #inverse  matrix
             dimensions:=[], pdimensions:=[], ## simple and proj. dim.
             IsDecompositionMatrix:=false,   # crystallized
             H:=H, operations:=CrystalMatrixOps,
             P:=function(d,mu)     ## a lazy helper
               return d.operations.P.S(d, d.H.operations.New("P",[1],[mu]));
             end);
    fi;
  end;   ## end of H.operations.NewDecompositonMatrix
      
  #U CartanMatrix(H,n)......................................................
  #U CartanMatrix(d)........................................................
  #M Return the Cartan matrix of H(n), or the decomposition matrix <d>.     
  CartanMatrix:=function(arg) local d, dmat, c;
    if Length(arg)=1 and IsDecompositionMatrix(arg[1]) then
      d:=arg[1];
    elif Length(arg)=2 then
      d:=DecompositionMatrix(arg[1],arg[2]);
    else Error("usage: CartanMatrix(H,n), or CartanMatrix(d)");
    fi;
    dmat:=MatrixDecompositionMatrix(d);
    c:=DecompositionMatrixMatrix(d.H,TransposedMat(dmat)*dmat,d.rows,d.cols);
    c.matname:="Cartan Matrix";
    return c;
  end;  

  ## This function will read the file "n" if n is a string; otherwise
  ## it will look for the relevant file for H(Sym_n). If crystal=true
  ## this will be a crystal decomposition matrix, otherwise it will be
  ## a 'normal' decomposition matrix. 
  H.operations.ReadDecompositionMatrix:=function(H,n, crystal)
    local msg, file, M, d, c, parts, coeffs, p, x, r, cm, rm;
  
    if crystal then 
      if IsString(n) then file:=n;
      else file:=Concatenation("e",String(H.e),"crys.",String(n));
      fi;
      msg:="ReadCrystalMatrix-";
    else 
      msg:="ReadDecompositionMatrix-";
      if IsString(n) then file:=n;
      else file:=Concatenation(H.HeckeRing,".",String(n));
      fi;
    fi;

    A_Specht_Decomposition_Matrix:=false;  ## just in case.
    d:=false;
  
    if not ReadPath("",file,"",Concatenation(msg,"CurrentDirectory")) then
      if SpechtDirectory="" or not 
      ReadPath(SpechtDirectory,file,"",Concatenation(msg,"SpechtDirectory"))
      then
        ReadPath(H.info.Library,file,"",Concatenation(msg,"SpechtLibrary"));
      fi;
    fi;
  
    if A_Specht_Decomposition_Matrix<>false then   ## extract matrix from M
      M:=A_Specht_Decomposition_Matrix;
      A_Specht_Decomposition_Matrix:=false;
      r:=Set(M.rows); c:=Set(M.cols);
      if H.IsSpecht and r=c then
        d:=H.operations.NewDecompositionMatrix(r,
              Filtered(c,x->IsERegular(H.e,x)),not IsBound(M.crystal));
      elif H.IsSpecht then
        d:=H.operations.NewDecompositionMatrix(r,c,not IsBound(M.crystal));
      else
        d:=H.operations.NewDecompositionMatrix(r,r,not IsBound(M.crystal));
      fi;
      if IsSet(M.rows) and IsSet(M.cols) then ## new format
        if IsBound(M.matname) then d.matname:=M.matname; fi;
        for c in [1..Length(d.cols)] do
          cm:=Position(M.cols, d.cols[c]);
          if cm<>false and IsBound(M.d[cm]) then
            x:=M.d[cm];
            parts:=[]; coeffs:=[];
            for rm in [1..Length(x)/2] do
              r:=Position(d.rows,M.rows[x[rm+Length(x)/2]]);
              if r<>false then
                Add(parts,r);
                if IsInt(x[rm]) then Add(coeffs,x[rm]);
                else
                  p:=Copy(v);
                  p.valuation:=x[rm][1];
                  p.coefficients:=x[rm]{[2..Length(x[rm])]};
                  Add(coeffs,p);
                fi;
              fi;
            od;
            if parts<>[] then   ## paranoia
              SortParallel(parts,coeffs);
              d.d[c]:=rec(parts:=parts,coeffs:=coeffs);
           fi;
          fi;
        od;
      else  ## old format
        d.d:=List(c, r->rec(coeffs:=[], parts:=[]));
        ## next, we unravel the decomposition matrix
        for rm in [1..Length(M.rows)] do
          r:=Position(d.rows,M.rows[rm]);
          if r<>false then
            x:=1;
            while x<Length(M.d[rm]) do
              c:=Position(d.cols,M.cols[M.d[rm][x]]);
              if c<>false then
                Add(d.d[c].coeffs, M.d[rm][x+1]);
                Add(d.d[c].parts, r);
              fi;
              x:=x+2;
            od;
          fi;
        od;
        for c in [1..Length(d.d)] do
          if d.d[c].parts=[] then Unbind(d.d[c]);
          else SortParallel(d.d[c].parts, d.d[c].coeffs);
          fi;
        od;
      fi;
    fi;
    if crystal and d<>false then d.operations.Store(d); fi;
    return d;
  end;

  ## Look up the decomposition matrix in the library files and in the
  ## internal lists and return it if it is known.
  ## NOTE: this function does not use the crystal basis to calculate the
  ## decomposition matrix because it is called by the various conversion
  ## functions X->Y which will only need a small part of the matrix in 
  ## general. The function FindDecompositionMatrix() also uses the crystal
  ## basis.
  H.operations.KnownDecompositionMatrix:=function(n) local d, x, r, c;
    d:=H.operations.ReadDecompositionMatrix(H,n,false);

    ## next we look for crystal matrices
    if d=false and H.IsSpecht and IsBound(H.Pq) then
      d:=H.operations.ReadDecompositionMatrix(H,n,true);
      if d<>false then d:=Specialized(d); fi;
    fi;

    if d=false and n<2*H.e then   
      ## decomposition matrix can be calculated
      r:=Partitions(n); 
      if H.IsSpecht then c:=ERegularPartitions(H.e,n);
      else c:=r;
      fi;
      d:=H.operations.NewDecompositionMatrix(r, c, true);

      for x in [1..Length(d.cols)] do
        if IsECore(H,d.cols[x]) then    ## an e-core
          d.operations.AddIndecomposable(d,
          H.operations.New("S",[1],[d.cols[x]]),false);
        elif IsSimpleModule(H,d.cols[x]) then ## first e-weight 1 partition
          c:=H.operations.New( "S",[1],[ECore(H,d.cols[x])] )
                   *H.operations.Omega("S",H.e);
          for r in [2..Length(c.parts)] do
            d.operations.AddIndecomposable(d,
                  H.operations.New("S",[1,1],c.parts{[r-1,r]}),false);
          od;
          if not H.IsSpecht then 
            d.operations.AddIndecomposable(d,
            H.operations.New("S",[1],[c.parts[1]]),false);
          fi;
        fi;
      od;
    fi;
    return d;
  end;  # end of KnownDecompositionMatrix

  ## almost identical to KnownDecompositionMatrix except that if this function
  ## fails then the crystalized decomposition matrix is calculated.
  H.operations.FindDecompositionMatrix:=function(n) local d,c;
    d:=H.operations.KnownDecompositionMatrix(n);
    if d=false and H.IsSpecht and IsBound(H.Pq) then
      d:=H.operations.NewDecompositionMatrix(
              Partitions(n),ERegularPartitions(H.e,n),false);
      for c in [1..Length(d.cols)] do
        if not IsBound(d.d[c]) then
          d.operations.AddIndecomposable(d,Pq(d.cols[c]),false);
        fi;
      od;
      d:=Specialized(d);
    fi;
    return d;
  end;
 
  #########################################################################
  ## Next, for fields of characteristic 0, we implement the LLT algorithm.
  ## The actual LLT algorithm is contained in the function Pq (should 
  ## really be called Pv...), but there are also functions Sq and Dq as
  ## well. These functions work as follows:
  ##   Pq(mu) -> sum_nu d_{nu,mu} S(nu)
  ##   Sq(mu) -> sum_nu d_{mu,nu} D(nu)
  ##   Dq(mu) -> sum_nu d_{mu,nu}^{-1} S(nu)
  ## The later two functions will call Pq() as needed. The "modules" x
  ## returned by these functions have x.module equal to "Pq", "Sq" and
  ## "Dq" respectively to distinguish them from the specialized versions.
  ## Accordingly we need H.operations.Xq for X = S, P, and D; these are
  ## defined after Pq(), Sq(), and Dq() (which they make use of).

  if H.p<>0 then v:=1;    ## just in case
  else                    ## define crystal basis operations

    ## first make sure that we have an indeterminant
    H.info.Indeterminate:=Indeterminate(Integers);
    v:=H.info.Indeterminate;
    if not IsBound(v.name) then H.info.Indeterminate.name:="v"; fi;

    ## next define the handling function H.Pq() - cut from Handler()
    HSC:=function(module,arg) local mu, z;
      mu:=Flat(arg);
      if not ForAll(mu,z->IsInt(z)) then 
        Error("usage: H.", module, "(<mu1,mu2,...>)\n");
      fi;

      z:=Copy(mu);
      Sort(mu, function(a,b) return a>b;end); # non-increasing
      if mu<>z then
        Print("## ", module,"(mu), warning <mu> is not a partition.\n");
      fi;
      if Length(mu)>0 and mu[Length(mu)]<0 then
        Error("## B", module,"(mu): <mu> contains negative parts.\n");
      fi;
      z:=Position(mu,0);
      if z<>false then mu:=mu{[1..z-1]}; fi;  ## remove any zeros from mu
  
      if module<>"Sq" then
        if H.IsSpecht and not IsERegular(H.e,mu) then
          Error("Pq(mu): <mu>=[",TightStringList(mu), 
                "] must be ", H.e,"-regular\n\n");
        else return Pq(mu);
        fi;
      else return H.operations.New(module,[1],[mu]);
      fi;
    end;

    H.Sq:=function(arg) return HSC("Sq", arg); end;
    H.Pq:=function(arg) return HSC("Pq", arg); end;

    # recrusively compute the LLT A(mu) function. Note that we can not
    # going through the handler function HSC for reasons of efficiency...
    # at some point I'd like to decouple all of the other functions to make
    # everything faster...
    # We assume that mu is regular...
    H.Aq:=function(mu) local removables, L, i, nu, r;
      if mu=[] then return H.operations.New("Sq",[v^0],[[]]); fi;
      a:=SPECHT.GetCached(H,"Aq",mu);
      if a<> false then return a; fi;
      removables:=StronglyNormalNodes(H.e,mu);
      L:=Maximum(List(removables,Length));
	    i:=First([1..H.e], s->Length(removables[s])=L);
      nu:=Copy(mu);
      for r in removables[i] do
        nu[r]:=nu[r]-1;
      od;
      if nu[Length(nu)]=0  then Unbind( nu[Length(nu)] ); fi;
      a:=qSInducedModule(H.Aq(nu), Length(removables[i]), (i-1) mod H.e);
      SPECHT.PutCached(H,"Aq",a);
      return a;
    end;

   #######################################################################
   ## Next we define qInduction and qRestriction functions

   # add n nodes of residue r to the partition y from the i-th row down
   # here exp is the exponent of the indeterminate
   qsinduced:=function(y, n, r, i, exp) local ny, j, z;
     ny:=[];
     for j in [i..Length(y)-n+1] do
       if y[j]>0 and r=(y[j]-j) mod H.e and (j=Length(y) or y[j]>y[j+1])
       then exp:=exp-n;               ## removeable node of residue r
       elif r=(y[j] - j + 1) mod H.e then
         if j=1 or y[j] < y[j-1] then ## addable node of residue r
           z:=Copy(y);
           z[j]:=z[j] + 1; # only one node of residue r can be added
           if n=1 then Add(ny, [exp,z] );   # no more nodes to add
           else Append(ny, qsinduced(z, n-1, r, j+1, exp));
           fi;
           exp:=exp+n;
         fi;
       fi;
     od;
     return ny;
   end;

   ## notice that we can't pull the trick to induce all residues at once
   ## that we used in InducedModule() etc. as we have to keep track of the
   ## number of addable and removable nodes of each residue. Rather than
   ## do this we just call this function as many times as necessary.
   qSInducedModule:=function(x,s,r) local coeffs, parts, y, z;
      coeffs:=[]; parts:=[];
      for y in [1..Length(x.parts)] do
        if r=( -Length(x.parts[y]) mod H.e) then # add a node to the bottom
          z:=Copy(x.parts[y]);
          Add(z,0);
          for z in qsinduced(z,s,r,1,0)  do
            if z[1]=0 then Add(coeffs, x.coeffs[y]);
            else Add(coeffs, x.coeffs[y]*v^z[1]);
            fi;
            if z[2][Length(z[2])]=0 then Unbind(z[2][Length(z[2])]); fi;
            Add(parts, z[2]);
          od;
        else
          for z in qsinduced(x.parts[y],s,r,1,0) do
            if z[1]=0 then Add(coeffs, x.coeffs[y]);
            else Add(coeffs, x.coeffs[y]*v^z[1]);
            fi;
            Add(parts, z[2]);
          od;
        fi;
      od;

      if coeffs=[] then return H.operations.New(x.module,[0],[[]]);
      else return H.operations.Collect(x.module, coeffs, parts);
      fi;
    end;  # qSInduce

    qrestricted:=function(y, n, r, i, exp) local ny, j, z;
      ny:=[];
      for j in [i,i-1..n] do
        if y[j]>0 and r=(y[j]+1-j) mod H.e and (j=1 or y[j]<y[j-1]) then
           exp:=exp-n;                 ## an addable node of residue r
        elif r=(y[j] - j) mod H.e then   ## removeable node of residue r
          if j=Length(y) or y[j] > y[j+1] then
            z:=Copy(y);
            z[j]:=z[j]-1;
            if z[j]=0 then Unbind(z[j]); fi;
            if n=1 then Add(ny, [exp,z]); # no mode nodes to remove
            else Append(ny, qrestricted(z, n-1, r, j-1, exp));
            fi;
            exp:=exp+n;
          fi;
        fi;
      od;
      return ny;
    end;
  
    ## string-restriction: remove m r's from each partition in x
    ## ** should allow restricting s-times also
    qSRestrictedModule:=function(x, s, r) local coeffs, parts, z, y, i, e, exp;
      e:=x.H.e;
      coeffs:=[]; parts:=[];
      if s=0 then return x.H.operations.New("Sq",[1],[[]]); fi;
      for y in [1..Length(x.parts)] do
        if -Length(x.parts[y]) mod H.e = r then exp:=-s; else exp:=0; fi;
        for z in qrestricted(x.parts[y], s, r, Length(x.parts[y]), exp) do
          if z[1]=0 then Add(coeffs, x.coeffs[y]); 
          else Add(coeffs, x.coeffs[y]*v^z[1]);
          fi;
          Add(parts, z[2]);
        od;
      od;
      if parts=[] then return x.H.operations.New("Sq",[0],[[]]);
      else return x.H.operations.Collect("Sq", coeffs, parts);
      fi;
    end;  # qSRestrict

    #######################################################################
    ## Retrieves or calculates the crystal basis element Pq(mu)
    Pq:=function(mu) 
      local n, x, r, c, mucore, removables, L, i, Hplus, abacus, muplus, 
            xplus, nu, s, y, eweight;
      if mu=[] then return H.operations.New("Sq",[v^0],[[]]); fi;
      n:=Sum(mu);
      ## first we see if we have already calculated Pq(mu)
      x:=SPECHT.GetCached(H,"Pq",mu);
      if x<>false then 
        return x; fi;
      mucore:=ECore(H.e,mu);
      eweight:=(Sum(mu)-Sum(mucore))/H.e;
      if mu=mucore or mu=[1..n]*0+1 then
        x:=H.operations.New("Sq",[v^0],[mu]);
      elif eweight=1 then
        x:=H.operations.New("Sq",[v^0],[ECore(H.e,mu)])
                   * H.operations.Omega("Sq",H.e);
        r:=Position(x.parts,mu);
        if r=1 then 
          x.parts:=x.parts{[1]};
          x.coeffs:=[v^0];
        else 
          x.parts:=x.parts{[r-1,r]};
          x.coeffs:=[v,v^0];
        fi;
      else
	      removables:=StronglyNormalNodes(H.e,mu);
	      L:=Maximum(List(removables,Length));
	      i:=First([1..H.e], s->Length(removables[s])=L);
	      if L=0 then # this must be a singular partition so we add an
	                  # empty runner and calculate with e+1...
	
	        # calculate using the empty runner trick from the J-M equating # paper
	        # first, deifne the "next" Hecke algebra - which is cached 
	        Hplus:=Specht(H.e+1);
	        # L is the length of the longest partition in the same
	        # block as mu
	        L:=Length(mucore)+n-Sum(mucore);
	        # add an empty runner to the abacus for mu
	        abacus:=EAbacusRunners(H.e,mu,L);
	        Add(abacus,[]);
	        muplus:=PartitionAbacusRunners(abacus);
          xplus:=Hplus.Pq(muplus);
          x:=H.operations.New("Sq",[v^0],[mu]);
	        for i in [1..Length(xplus.parts)-1] do
	          if Length(xplus.parts[i])<=L then
	            abacus:=EAbacusRunners(Hplus.e,xplus.parts[i],L);
	            if abacus[Hplus.e]=[] then # partition of form nu^+
	              Unbind(abacus[Hplus.e]);
	              x:=x+H.operations.New("Sq",[xplus.coeffs[i]],
	                                         [PartitionAbacusRunners(abacus)]);
	            fi;
	          fi;
	        od;
        else  ## we calculate Pq(mu) recursively using LLT
          ## don't want to change the original mu
          nu:=Copy(mu);
	        for r in removables[i] do
	          nu[r]:=nu[r]-1;
	        od;
          if nu[Length(nu)]=0  then Unbind( nu[Length(nu)] ); fi;
          x:=qSInducedModule(Pq(nu), Length(removables[i]), (i-1) mod H.e);
          ## if H.PartitionLength has been set limit calculations to
          ## partitions of that length.
          if IsBound(H.PartitionLength) then 
            s:=Filtered([1..Length(x.parts)], 
                          r->Length(x.parts[r])<=H.PartitionLength);
            x.coeffs:=x.coeffs{s};
            x.parts:=x.parts{s};
          fi;
          n:=1;
          while n<Length(x.parts) do
            if x.coeffs[Length(x.parts)-n].valuation>0 then n:=n+1;
            else
              r := Copy( x.coeffs[Length(x.parts)-n] );
              nu:=x.parts[Length(x.parts)-n];
              if Length(r.coefficients) < 1-r.coefficients then
                Append(r.coefficients,
                  List([1..Length(r.coefficients)-1-r.valuation], i->0));
              fi;
              r.coefficients := r.coefficients{[1..1-r.valuation]};
              Append(r.coefficients, Reversed(r.coefficients{[1..-r.valuation]}));
              x := x-r*Pq(nu);
              if nu in x.parts then n:= n+1; fi;
            fi;
          od;
          r:=Filtered([1..Length(x.coeffs)],s->s<>0*s);
          x.coeffs:=x.coeffs{r};
          x.coeffs:=x.coeffs{r};
        fi;
      fi;
      ## having found x we add it to the cache
      ## (might be [] if IsBound(H.PartitionLength) )
      if x.parts<>[] then 
        SPECHT.PutCached(H,"Pq",x); 
        ## When mu is regular we also add the Mullineux image of Pq(mu) to the cache
        ## (see LLT Theorem 7.2)
        #if not IsBound(H.PartitionLength) and IsERegular(H,mu) then
        #  if ConjugatePartition(x.parts[1])<>mu then
        #    y:=ShallowCopy(x);
        #    ## not self-image under MullineuxMap
        #    ## (but only when H.PartitionLength is not set)
        #    y.coeffs:=List(x.coeffs*v^0, i->v^eweight*Value(i,v^-1));     ## v -> v^-1
        #    y.parts:=List(x.parts, ConjugatePartition);
        #    SortParallel(y.parts,y.coeffs);
        #    SPECHT.PutPq(H,y);
        #  fi;
        #fi;
      fi;
      return x;
    end;

    ## Strictly speaking the functions Sq() and Dq() are superfluous as
    ## these functions could be carried out by the decomposition matrix
    ## operations; however the point is that the crystal matrices are
    ## allowed to have missing columns, and if any needed columns are 
    ## missing then these functions calculate them on the fly via using Pq().

    ## writes S(mu), from the Fock space, as a sum of D(nu)
    Sq:=function(mu) local n, parts, x, r, Pr, c;
      n:=Sum(mu);
      parts:=PartitionsInSameBlock(H,mu);
      if H.IsSpecht then 
        parts:=Filtered(parts,nu->Dominates(nu,mu) and IsERegular(H,nu));
      else
        parts:=Filtered(parts,nu->Dominates(nu,mu));
      fi;

      x:=H.operations.New("Dq",[0],[[]]);
      r:=Length(parts);
      while r>0 do
        Pr:=Pq(parts[r]);
        c:=Position(Pr.parts,mu);
        if c<>false then x:=x+H.operations.New("Dq",[Pr.coeffs[c]],[parts[r]]); fi;
        r:=r-1;
      od;
      return x;
    end;

    ## write D(mu), from the Fock space, as a sum of S(nu)
    Dq:=function(mu) local inv, x, c;
      inv:=H.operations.New("Sq",[1],[mu]);
      x:=H.operations.New("Dq",[1],[mu])-Sq(mu);
      while x<>0*x do
        c:=Length(x.parts);
        inv:=inv+H.operations.New("Sq",[x.coeffs[c]],[x.parts[c]]);
        x:=x-x.coeffs[c]*Sq(x.parts[c]);
      od;
      return inv;
    end;

    ## Now that we have our functions Pq(), Sq(), and Dq() we define the
    ## operations records for the Fock space elements. The basic operation 
    ## set is the same as for 'ordinary' modules

    H.operations.Sq:=Copy(H.operations.S);
    H.operations.Pq:=Copy(H.operations.P);
    H.operations.Dq:=Copy(H.operations.D);
   
    ## however there are a few obvious differences

    H.operations.Sq.InducedModule:=function(x, list) local r;
      if list=[] then return Sum([0..H.e-1],r->qSInducedModule(x,1,r));
      elif H.e=0 then
        Error("Induce, r-induction is not defined when e=0.");
      elif ForAny(list,r-> r>=H.e or r<0) then
        Error("Induce, r-induction is defined only when 0<=r<e.\n");
      else 
        for r in list do   ## we could do slightly better here
          x:= qSInducedModule(x,1,r);
        od;
        return x;
      fi;
    end;
    H.operations.Pq.InducedModule:=function(x,list)
      return H.operations.Sq.P(
        H.operations.Sq.InducedModule(x.operations.S(x,false),list),false);
    end;
    H.operations.Dq.InducedModule:=function(x,list)
      return H.operations.Sq.D(
        H.operations.Sq.InducedModule(x.operations.S(x,false),list),false);
    end;

    H.operations.Sq.SInducedModule:=function(x, list) local r;
      if Length(list)=1 then 
        list:=list[1];
        if list=0 then return H.operations.New("Sq",[1],[[]]); fi;
        while list > 0 do
          x:=Sum([0..H.e-1],r->qSInducedModule(x,1,r));
          list:=list-1;
        od;
        return x;
      elif H.e=0 then
        Error("SInduce, r-induction is not defined when e=0.");
      elif list[2]>H.e or list[2]<0 then
        Error("SInduce, r-induction is defined only when 0<=r<e.\n");
      else return qSInducedModule(x, list[1], list[2]);
      fi;
    end;
    H.operations.Pq.SInducedModule:=function(x,list)
      return H.operations.Sq.P(
        H.operations.Sq.SInducedModule(x.operations.S(x,false),list),false);
    end;
    H.operations.Dq.SInducedModule:=function(x,list)
      return H.operations.Sq.D(
        H.operations.Sq.SInducedModule(x.operations.S(x,false),list),false);
    end;

    H.operations.Sq.RestrictedModule:=function(x, list) local r;
      if list=[] then return Sum([0..H.e-1],r->qSRestrictedModule(x,1,r));
      elif H.e=0 then
        Error("Restrict, r-restriction is not defined when e=0.");
     elif ForAny(list,r-> r>=H.e or r<0) then
        Error("Restrict, r-restriction is defined only when 0<=r<e.\n");
      else
        for r in list do   ## we could do slightly better here
          x:= qSRestrictedModule(x,1,r);
        od;
        return x;
      fi;
    end;
    H.operations.Pq.RestrictedModule:=function(x,list)
      return H.operations.Sq.P(
        H.operations.Sq.RestrictedModule(x.operations.S(x,false),list),false);
    end;
    H.operations.Dq.RestrictedModule:=function(x,list)
      return H.operations.Sq.D(
        H.operations.Sq.RestrictedModule(x.operations.S(x,false),list),false);
    end;

    H.operations.Sq.SRestrictedModule:=function(x, list) local r;
      if Length(list)=1 then
        list:=list[1];
        if list=0 then return H.operations.New("Sq",[1],[[]]); fi;
        while list>0 do
          x:=Sum([0..H.e-1],r->qSRestrictedModule(x,1,r));
          list:=list-1;
        od;
        return x;
      elif H.e=0 then
        Error("SRestrict, r-restriction is not defined when e=0.");
      elif list[2]>H.e or list[2]<0 then
        Error("SRestrict, r-restriction is defined only when 0<=r<e.\n");
      else return qSRestrictedModule(x, list[1], list[2]);
      fi;
    end;
    H.operations.Pq.SRestrictedModule:=function(x,list)
      return H.operations.Sq.P(
        H.operations.Sq.SRestrictedModule(x.operations.S(x,false),list),false);
    end;
    H.operations.Dq.SRestrictedModule:=function(x,list)
      return H.operations.Sq.D(
        H.operations.Sq.SRestrictedModule(x.operations.S(x,false),list),false);
    end;

    ## Specialization taking Xq -> X
    H.operations.Sq.Specialized:=function(x,a) local coeffs, c;
      coeffs:=List(x.coeffs*v^0,c->Value(c,a));
      c:=Filtered([1..Length(coeffs)],c->c<>c*0);
      if c<>[] then return H.operations.New(x.module{[1]},coeffs{c},x.parts{c});
      else return H.operations.New(x.module{[1]},[0],[[]]);
      fi;
    end;
    H.operations.Pq.Specialized:=H.operations.Sq.Specialized;
    H.operations.Dq.Specialized:=H.operations.Sq.Specialized;

    ## This also needs to be done for crystal matrices (this wasn't
    ## put in above because normal decomposition matrices don't
    ## specialise.
    CrystalMatrixOps.Specialized:=function(d,a) local sd, c, p, coeffs;
      sd:=H.operations.NewDecompositionMatrix(d.rows,d.cols,true);
      for c in [1..Length(d.cols)] do
        if IsBound(d.d[c]) then
          sd.d[c]:=rec();
          coeffs:=List(d.d[c].coeffs*v^0,p->Value(p,a));
          p:=List(coeffs,p->p<>0);
          if true in p then
            sd.d[c]:=rec(coeffs:=ListBlist(coeffs,p),
                          parts:=ListBlist(d.d[c].parts,p) );
          else sd.d[c]:=rec(coeffs:=[0], parts:=[ [] ] );
          fi;
        fi;
        if IsBound(d.inverse[c]) then
          coeffs:=List(d.inverse[c].coeffs*v^0,p->Value(p,a));
          p:=List(coeffs,p->p<>0);
          if true in p then
            sd.inverse[c]:=rec(coeffs:=ListBlist(coeffs,p),
                            parts:=ListBlist(d.inverse[c].parts,p) );
          else sd.d[c]:=rec(coeffs:=[0], parts:=[ [] ] );
          fi;
        fi;
      od;
      return sd;
    end;

    PositiveCrystalCoefficients:=function(x) local c, p;
      return ForAll(x.coeffs, p->( IsInt(p) and p>=0 ) or
                                ForAll(p.coefficients, c->c>=0) );
    end;
    H.operations.Sq.PositiveCoefficients:=PositiveCrystalCoefficients;
    H.operations.Pq.PositiveCoefficients:=PositiveCrystalCoefficients;
    H.operations.Dq.PositiveCoefficients:=PositiveCrystalCoefficients;

    IntegralCrystalCoefficients:=function(x) local c, p;
      return ForAll(x.coeffs, p-> ( IsInt(p) and p>=0 ) or
                                ForAll(p.coefficients, c->IsInt(c)) );
    end;
    H.operations.Sq.IntegralCoefficients:=IntegralCrystalCoefficients;
    H.operations.Pq.IntegralCoefficients:=IntegralCrystalCoefficients;
    H.operations.Dq.IntegralCoefficients:=IntegralCrystalCoefficients;

    ## Finally, change the various conversion functions X()->Y();
    ## in fact, we only have to change the four non-trivial ones:
    ##   P() <-> S() <-> D().

    H.operations.Sq.P:=function(x,silent) local proj;
      if x=false or x=0*x then return x;
      elif x.parts=[[]] then return H.operations.New("Pq",[x.coeffs[1]],[[]]);
      fi;

      proj:=H.operations.New("Pq",[0],[[]]);
      while x<>0*x and PositiveCrystalCoefficients(x) do
        proj:=proj+H.operations.New("Pq",[x.coeffs[Length(x.parts)]], 
                                         [x.parts[Length(x.parts)]]);
        x:=x-x.coeffs[Length(x.parts)]*Pq(x.parts[Length(x.parts)]);
      od;
      if x=0*x then return proj;
      elif not silent then
        Print("# P(<x>), unable to rewrite <x> as a sum of projectives\n");
      fi;
      return false;
    end;

    H.operations.Sq.D:=function(x,silent) local mu;
      if x=false or x=0*x then return x;
      elif x.parts=[[]] then return H.operations.New("Dq",[x.coeffs[1]],[[]]);
      fi;
      return Sum([1..Length(x.parts)], mu->x.coeffs[mu]*Sq(x.parts[mu]) );
    end;

    H.operations.Pq.S:=function(x,silent) local mu;
      if x=false or x=0*x then return x;
      elif x.parts=[[]] then return H.operations.New("Sq",[x.coeffs[1]],[[]]);
      fi;
     
      return Sum([1..Length(x.parts)], mu->x.coeffs[mu]*Pq(x.parts[mu]) );
    end;

    H.operations.Dq.S:=function(x,silent) local mu;
      if x=false or x=0*x then return x;
      elif x.parts=[[]] then return H.operations.New("Sq",[x.coeffs[1]],[[]]);
      fi;

      return Sum([1..Length(x.coeffs)], mu->x.coeffs[mu]*Dq(x.parts[mu]) );
    end;

  fi;           ## end of if-statement for crystal basis functions

  ########################### Special cases ##############################

  ## Finally we add functions to H.operations.X() which calculate PIMS etc
  ## in various known special cases (at present these are all for e=2 and
  ## p=0 and so these functions are now redundant because of the
  ## LLT-Grojnowski algorithm. They are still included, although not
  ## automatically accessed because they are more efficient than the
  ## crystal basis algorithm for large n (being explicit), and because I
  ## don't want to lose the code. Note that the conversions functions
  ## all make calls of the form H.operations.X.e#p#Y(); these refer to
  ## the functions below (but are never called is IsBound(H.Pq).

  if H.e=2 and H.p=0 then

    ## This function writes S(mu), for mu e-regular, as a sum of
    ## simples in the case when e=2 and mu has at most four parts.
    ## Called via D(linear comb. S(nu)). This is quite horrible; we just
    ## have a list of the simples in each Specht module.
    H.operations.S.e2p0D:=function(mu) local mumod, m;
  
      if Length(mu)>4 or (Length(mu)>2 and not IsERegular(2,mu)) then
        return false;
      fi;
  
      mumod:=List(mu, m->m mod 2);
      if Length(mu)=1 then return H.operations.Collect("D",1,mu);
      elif Length(mu)=2 then
        if (mu[1] + mu[2]) mod 2<>0 then 
          return H.operations.Collect("D",1,mu); 
        else 
          return H.operations.Collect("D",[1,1],[ mu,[mu[1]+1, mu[2]-1] ]);
        fi;
      elif Length(mu)=3 then
        if mumod=[0,0,0] or mumod=[1,1,1] then
          return H.operations.Collect("D",[1,1,1,1], 
            [ [mu[1]+2,mu[2],mu[3]-2],
              [mu[1]+1,mu[2]-1,mu[3]],[mu[1],mu[2]+1,mu[3]-1], mu]);
        elif mumod=[0,0,1] or mumod=[1,1,0] then
          return H.operations.Collect("D",[1,1,1,1],
           [ [mu[1]+1,mu[2]+1,mu[3]-2],[mu[1]+1,mu [2],mu[3]-1],
             [mu[1]+1,mu[2]-1,mu[3]],mu] );
        elif mumod=[0,1,0] or mumod=[1,0,1] then 
          return H.operations.Collect("D",1,mu);
        elif mumod=[0,1,1] or mumod=[1,0,0] then
          return H.operations.Collect("D",[1,1,1,1], 
            [ [mu[1]+2,mu[2]-1,mu[3]-1],[mu[1]+1,mu[2],mu[3]-1],
              [mu[1],mu[2]+1,mu[3]-1], mu] );
        fi;
      elif Length(mu)=4 then
        if mumod=[0,0,0,0] or mumod=[1,1,1,1] then
          return H.operations.Collect("D", 
            [1,1,2,1,1,1,1,1,2,1,1,1,1,1,1,1],
            [ [mu[1]+3,mu[2]+1,mu[3]-1,mu[4]-3],
              [mu[1]+3,mu[2]-1,mu[3]-1,mu[4]-1],
              [mu[1]+2,mu[2],mu[3],mu[4]-2],[mu[1]+2,mu[2]+2,mu[3]-2,mu[4]-2],
              [mu[1]+2,mu[2],mu[3]-1,mu[4]-1],[mu[1]+2,mu[2],mu[3]-2,mu[4]],
              [mu[1]+1,mu[2]+1,mu[3]+1,mu[4]-3],
              [mu[1]+1,mu[2]+1,mu[3],mu[4]-2],
              [mu[1]+1,mu[2]+1,mu[3]-1,mu[4]-1],[mu[1]+1,mu[2],mu[3],mu[4]-1],
              [mu[1]+1,mu[2]-1,mu[3]+1,mu[4]-1],[mu[1]+1,mu[2]-1,mu[3],mu[4]],
              [mu[1],mu[2]+2,mu[3],mu[4]-2], [mu[1],mu[2]+1,mu[3]-1,mu[4]],
              [mu[1],mu[2],mu[3]+1,mu[4]-1], mu] );
        elif mumod=[0,0,0,1] or mumod=[1,1,1,0] then
          m:=H.operations.Collect("D",[1,1,1,1,1,1,1,1], 
              [ [mu[1]+1,mu[2]+2,mu[3]-1,mu[4]-2],
              [mu[1]+1,mu[2]+1,mu[3],mu[4]-2],[mu[1]+1,mu[2]+1,mu[3]-2,mu[4]],
              [mu[1]+1,mu[2],mu[3]-1,mu[4]],[mu[1]+1,mu[2],mu[3]+1,mu[4]-2],
              [mu[1]+1,mu[2],mu[3],mu[4]-1],
              [mu[1]+1,mu[2]-1,mu[3],mu[4]], mu] );
          if mu{[3,4]}=[2,1] then 
            return m+H.operations.Collect("D",1,[mu[1]+2,mu[2]+1]);
          else return m;
          fi;
        elif mumod=[0,0,1,0] or mumod=[1,1,0,1] then
          m:=H.operations.Collect("D",[1,1,1,1,1,1,1,1],
             [ [mu[1]+1,mu[2]+2,mu[3]-1,mu[4]-2],
               [mu[1]+1,mu[2]+1,mu[3],mu[4]-2],
               [mu[1]+1,mu[2]+1,mu[3]-2,mu[4]],
               [mu[1]+1,mu[2],mu[3]-1,mu[4]],[mu[1]+1,mu[2],mu[3]+1,mu[4]-2],
               [mu[1]+1,mu[2],mu[3],mu[4]-1],[mu[1]+1,mu[2]-1,mu[3],mu[4]],
               mu] );
          if mu{[3,4]}=[2,1] then 
            return m+H.operations.Collect("D",1,[mu[1]+1,mu[2]+2]);
          else return m;
          fi;    
        elif mumod=[0,0,1,1] or mumod=[1,1,0,0] then
          return H.operations.Collect("D",[1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1],
            [ [mu[1]+3,mu[2],mu[3],mu[4]-3], [mu[1]+3,mu[2],mu[3]-2,mu[4]-1],
              [mu[1]+2,mu[2]-1,mu[3]+1,mu[4]-2],
              [mu[1]+2,mu[2]+1,mu[3]-1,mu[4]-2],
              [mu[1]+2,mu[2]-1,mu[3],mu[4]-1],[mu[1]+2,mu[2]-1,mu[3]-1,mu[4]],
              [mu[1]+1,mu[2]+2,mu[3],mu[4]-3], 
              [mu[1]+1,mu[2]+2,mu[3]-2,mu[4]-1],
              [mu[1]+1,mu[2]+1,mu[3]-1,mu[4]-1],
              [mu[1]+1,mu[2],mu[3]+1,mu[4]-2], [mu[1]+1,mu[2],mu[3]-1,mu[4]], 
              [mu[1]+1,mu[2],mu[3],mu[4]-1], [mu[1],mu[2]+1,mu[3]+1,mu[4]-2], 
              [mu[1],mu[2]+1,mu[3],mu[4]-1], [mu[1],mu[2]+1,mu[3]-1,mu[4]],
               mu]  );
        elif mumod=[0,1,0,0] or mumod=[1,0,1,1] then
          return H.operations.Collect("D",[1,1,1,1,1,1,1,1],
            [ [mu[1]+2,mu[2]+1,mu[3]-2,mu[4]-1],
              [mu[1]+2,mu[2],mu[3]-1,mu[4]-1],[mu[1]+2,mu[2]-1,mu[3],mu[4]-1],
              [mu[1]+1,mu[2],mu[3],mu[4]-1],[mu[1],mu[2]+2,mu[3]-1,mu[4]-1],
              [mu[1],mu[2]+1,mu[3],mu[4]-1],[mu[1],mu[2],mu[3]+1,mu[4]-1],
               mu] );
        elif mumod=[0,1,0,1] or mumod=[1,0,1,0] then 
          return H.operations.Collect("D",1,mu);
        elif mumod=[0,1,1,0] or mumod=[1,0,0,1] then 
          return H.operations.Collect("D",[1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1],
            [ [mu[1]+3,mu[2],mu[3],mu[4]-3], [mu[1]+3,mu[2],mu[3]-2,mu[4]-1],
              [mu[1]+2,mu[2]-1,mu[3]+1,mu[4]-2], 
              [mu[1]+2,mu[2]+1,mu[3]-1,mu[4]-2],
              [mu[1]+2,mu[2]-1,mu[3],mu[4]-1],
              [mu[1]+2,mu[2]-1,mu[3]-1,mu[4]],
              [mu[1]+1,mu[2]+2,mu[3]-2,mu[4]-1],
              [mu[1]+1,mu[2]+1,mu[3]-1,mu[4]-1],
              [mu[1]+1,mu[2],mu[3]+1,mu[4]-2],
              [mu[1]+1,mu[2],mu[3],mu[4]-1],[mu[1]+1,mu[2],mu[3]-1,mu[4]],
              [mu[1],mu[2]+1,mu[3]+1,mu[4]-2],[mu[1],mu[2]+1,mu[3],mu[4]-1],
              [mu[1],mu[2]+1,mu[3]-1,mu[4]],[mu[1],mu[2],mu[3]+1,mu[4]-1], 
               mu] );
        elif mumod=[0,1,1,1] or mumod=[1,0,0,0] then
          return H.operations.Collect("D",[1,1,1,1,1,1,1,1], 
           [ [mu[1],mu[2],mu[3]+1,mu[4]-1],
             [mu[1]+2,mu[2]+1,mu[3]-1,mu[4]-2],[mu[1]+2,mu[2],mu[3],mu[4]-2],
             [mu[1]+2,mu[2]-1,mu[3]-1,mu[4]],[mu[1]+1,mu[2],mu[3]-1,mu[4]],
             [mu[1],mu[2]+2,mu[3],mu[4]-2],[mu[1],mu[2]+1,mu[3]-1,mu[4]], 
              mu] );
        fi;
      fi;
    end;     ## e2p0D()

    ## here arg=lambda is a "staircase" partition (see [JM1]), and this
    ## function computes the corresponding PIM. 
    H.operations.P.e2p0SingleRegular:=function(arg) 
      local lambda,c,mu,x,pim,i,q;
  
      lambda:=Flat(arg);
      if not IsERegular(2,lambda) then return false; fi;
      c:=ECore(2,lambda);
      if c=lambda then return H.operations.New("S",[1],[c]); fi;
      x:=EQuotient(2,lambda);
      mu:=[];
      for i in [1..Length(x)] do
        if x[i]<>[] then
          if mu=[] then
            mu:=x[i];
            if i=1 then
              q:=function(a) return [a[1], ConjugatePartition(a[2])]; end;
            else
              q:=function(a) return [ConjugatePartition(a[1]), a[2]]; end;
            fi;
          else
            Error("SingleReg, lambda is not a suitable partition\n");
          fi;
        fi;
      od;
  
      mu:=Collected(InverseLittlewoodRichardsonRule(mu));
      pim:=H.operations.New("S",[0],[[]]);
      for x in mu do
        pim:=pim + x[2]
          * H.operations.New("S",[1],[CombineEQuotientECore(q(x[1]),c)]);
      od;
      return pim;
    end;  # H.operations.P.e2p0SingleRegular
  
    ## given a singular partition mu this expands S(mu) into a sum of regular
    ## Specht modules S(nu). e=2 only though... 
    ## ***** undocumented and unused
    H.operations.P.e2p0ExpandSingular:=function(arg) 
      local Twos, Vert, x, y, ymu, i;
  
      Vert:=function(mu) local mud, vmu, i, v;
        mud:=ConjugatePartition(mu);
        vmu:=[];
        repeat
          i:=Length(mud);
          if mud[i] >=2 then
            mud[i]:=mud[i] - 2;
            v:=1;
          else v:=0;
          fi;
          while i > 1 do
            i:=i - 1;
            if mud[i] - 2 >=mud[i+1] then
              mud[i]:=mud[i] - 2;
              v:=v + 1;
            fi;
          od;
          if v > 0 then Add(vmu, v); fi;
        until v=0;
        return [ConjugatePartition(mud), vmu];
      end;
  
      Twos:=function(t)
        return H.operations.New("S",[(-1)^t],[[2*t]])
           +Sum([0..t-1], i->H.operations.New("S",[(-1)^i],[[t+i,t-i]]) );
      end;
  
      if Length(arg)=1 and IsSpecht(arg[1]) then x:=arg[1];
      else
        y:=Flat(arg);
        ymu:=Vert(y);
        x:=H.operations.New("S",[1],[y]) 
            - H.operations.New("S",[1],[ymu[1]])*Product(ymu[2], i->Twos(i));
      fi;
  
      repeat
        y:=1;
        while y<=Length(x.parts)  and IsERegular(x.H.e,x.parts[y]) do
          y:=y+1;
        od;
        if y<=Length(x.parts) then
          ymu:=Vert(x.parts[y]);
          x:=x - H.operations.New("S",[x.coeffs[y]],[ymu[1]])
                    * Product(ymu[2], i->Twos(i));
        fi;
      until y>Length(x.parts);
      return x;
    end; # H.e2p0ExpandSingular 
  
    ## (Theorem) formulae for the PIMS of the two part partitions when e=2
    H.operations.P.e2p0S_two:=function(k, l)
      local pim, x, y, a, n, k2, l2, ltop, K, L, lambda,
            pimfn, eveneven, Oddeven, OddOdd;
  
      Oddeven:=function(x,y,a)
        if y=0 then
          if x mod 2=k2 then pim:=pim+H.operations.New("S",[1],
                                    [H.operations.DoubleHook(n,x,y,a)]); 
          fi;
        else
          if x mod 2=l2 then
            if y mod 2=l2 and a mod 2=l2 then
              pim:=pim
                +H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
            fi;
          elif x mod 2=k2 and y mod 2=k2 then
            if a mod 2=l2 then 
              pim:=pim
                +H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]); 
            fi;
          else 
            pim:=pim
              +H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
          fi;
        fi;
      end;
  
      OddOdd:=function(x,y,a)
        if y=0 then 
          pim:=pim+H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
        else
          if x mod 2<>y mod 2 then
            if a mod 2=x mod 2 then 
              pim:=pim
                +H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]); 
            fi;
          elif x<>y or x mod 2=0 then 
            pim:=pim 
              + H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
          fi;
        fi;
      end;
  
      eveneven:=function(x,y,a)
        if y=0 then 
          pim:=pim+H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
        else
          if a=k-x and a=l+1-y then lambda:=0;
          else
            if a=k-x or a=l+1-y or a=l-x then lambda:=1;
            else lambda:=2;
            fi;
            if x mod 2<>y mod 2 then
              if a mod 2=x mod 2 then
                pim:=pim+lambda       
                * H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
              fi;
            else pim:=pim+lambda
                * H.operations.New("S",[1],[H.operations.DoubleHook(n,x,y,a)]);
            fi;
          fi;
        fi;
      end;
  
      n:=k + l;
      pim:=H.operations.New("S",[0],[[]]);
  
      k2:=k mod 2; l2:=l mod 2;
      if k2<>l2 then
        pimfn:=Oddeven;
        ltop:=l;
        L:=function(x,y) return Maximum(l-x,0);end;
        K:=function(x,y) return Minimum(k-x,l-y);end;
      elif l2=0 then
        pimfn:=eveneven;
        ltop:=l+1;
        L:=function(x,y) return Maximum(l-x,0);end;
        K:=function(x,y) return Minimum(k-x,l+1-y);end;
      else
        pimfn:=OddOdd;
        ltop:=l + 1;
        L:=function(x,y)
          if y mod 2=0 then return Maximum(l-x,0);
          else return Maximum(l+1-x,0);
          fi;
        end;
        K:=function(x,y)
          if   x mod 2=0 and y mod 2=0 then return Minimum(k-1-x,l+1-y);
          elif x mod 2=0 and y mod 2=1 then return Minimum(k-1-x,l-y);
          elif x mod 2=1 and y mod 2=0 then return Minimum(k-x,l+1-y);
          else return Minimum(k-x,l-y);
          fi;
        end;
      fi;
  
      for x in [2..k] do
        if x > l then pimfn(x,0,0); fi;
        for y in [2..Minimum(ltop,x)] do
          for a in [L(x,y)..K(x,y)] do
            pimfn(x,y,a);
          od;
        od;
      od;
      return pim;
    end;  # P.e2p0S_two
  
    H.operations.P.e2p0S_three:=function(k,l,m) local th, mum;
      mum:=[k mod 2, l mod 2, m mod 2];
  
      if mum=[0,0,0] then
        th:=RestrictedModule(SRestrictedModule(SRestrictedModule(
           H.operations.P.e2p0SingleRegular(k+1,l+2,m+3,2,1),H.e,5,0),
             H.e,3,1),H.e,0);
        if m>2 then
        ## In [JM1] there is an ambiguity in the projectives of this sort.
        ## The following use of Schaper's theorem seems to resolve it.
        mum:=th-SRestrictedModule(
                  H.operations.P.e2p0SingleRegular(k-2,l-1,m,3,2,1),H.e,3,1);
          if H.operations.S.InnerProduct(mum,Schaper(th.H, 
                        H.operations.Hook(k+l+m,k-2,m,m,3)))=0 then
            Print("# Warning: the projective P(",k,",",l,",",m,") is partly ",
                  "conjectural.\n");
          fi;
        fi;
        return th;
      elif mum=[0,0,1] then
       return SRestrictedModule(SRestrictedModule(
         H.operations.P.e2p0SingleRegular(k+1,l+2,m+2,2,1),H.e,5,0),H.e,3,1);
      elif mum=[0,1,0] then
        return RestrictedModule(SRestrictedModule(
          H.operations.P.e2p0SingleRegular(k+1,l+1,m+1,2,1),H.e,5,0),H.e,1);
      elif mum=[0,1,1] then
        return SRestrictedModule(SRestrictedModule(
          H.operations.P.e2p0SingleRegular(k+1,l+1,m+2,2,1),H.e,5,0),H.e,2,1);
      elif mum=[1,0,0] then
        return RestrictedModule(H.operations.P.e2p0S_three(k,l,m+1),H.e,0);
        th:=RestrictedModule(RestrictedModule(SRestrictedModule(
         H.operations.P.e2p0SingleRegular(k,l,m+1,2,1),H.e,2,0),H.e,1),H.e,0);
        return th;
      elif mum=[1,0,1] then
        if m=1 then return H.operations.P.e2p0SingleRegular(k,l,m);
        else
          th:=RestrictedModule(SRestrictedModule(
               H.operations.P.e2p0SingleRegular(k,l,m,2,1),H.e,2,0),H.e,1);
          if m>2 then
            if m > 4 then 
              th:=th-RestrictedModule(
                       H.operations.P.e2p0SingleRegular(k,l,m-2,2,1),H.e,0);
            fi;
            if l > m+2 then
              th:=th-RestrictedModule(
                       H.operations.P.e2p0SingleRegular(k,l-2,m,2,1),H.e,0);
            fi;
            if k > l+2 then
              th:=th-RestrictedModule(
                       H.operations.P.e2p0SingleRegular(k-2,l,m,2,1),H.e,0);
            fi;
          fi;
        fi;
        return th;
      elif mum=[1,1,0] then
        th:=SRestrictedModule(H.operations.P.e2p0S_three(k,l+1,m+1),H.e,2,0);
        if m > 2 then
          th:=th - SRestrictedModule(
                  H.operations.P.e2p0SingleRegular(k,l+1,m-1,2,1),H.e,3,0);
        fi;
        if l<>m + 1 then
          th:=th - SRestrictedModule(
                    H.operations.P.e2p0SingleRegular(k,l-1,m+1,2,1),H.e,3,0);
        fi;
        if k<>l + 2 then
          th:=th-SRestrictedModule(
                   H.operations.P.e2p0SingleRegular(k-2,l+1,m+1,2,1),H.e,3,0);
        fi;
        return th;
      else # mum=[1,1,1]
        th:=InducedModule(H.operations.P.e2p0S_three(k-1,l,m),H.e,0);
        if m > 3 then
          return th - SRestrictedModule(SInducedModule(
          H.operations.P.e2p0SingleRegular(k-3,l-2,m-1,3,2,1),H.e,2,0),H.e,2,1);
        else return th;
        fi;
      fi;
    end; # H.operations.P.e2p0S_three
  
    # the case e=2: P->S
    H.operations.P.e2p0S:=function(mu) local n, i, j;
      if not IsERegular(2,mu) then return false; fi;
  
      n:=Sum(mu);
      if Length(mu)=1 then  # mu=(n)
        if n mod 2=0 then 
          return H.operations.Collect("S",List([1..n], i->1), 
                       List([1..n],i->H.operations.Hook(n,n-i+1)));
        else 
          return H.operations.Collect("S",List([0..(n-1)/2],i->1), 
                           List([0..(n-1)/2],i->H.operations.Hook(n,n-2*i)));
        fi;
      elif Length(mu)=2 then return H.operations.P.e2p0S_two(mu[1], mu[2]);
      elif Length(mu)=3 then 
        return H.operations.P.e2p0S_three(mu[1],mu[2],mu[3]);
      fi;
      
      ## If we have a "staircase partition", or a restriction of one then
      ## [JM1] tells us how to construct P(mu). This first bit gets most
      ## of them, and we get the remainder by sending ourselves back into
      ## this case.
      i:=Length(mu);
      if mu[i]=1 then
        while mu[i-1]=mu[i]+1 do i:=i-1; od;  # note, mu not a 2-core => i>1
        if Length(mu)+4>=2*i and IsSimpleModule(H, mu{[1..i-1]}) then
          if (Length(mu) mod 2)=(mu[1] mod 2) then
            if Length(mu)+3>=2*i 
              then return H.operations.P.e2p0SingleRegular(mu);
            else i:=1;
            fi;
          fi;
          mu:=Copy(mu);
          for j in [i..Length(mu)] do 
            mu[j]:=mu[j]+1; 
          od;
          Add(mu,1);
          return SRestrictedModule(H.operations.P.e2p0SingleRegular(mu),H.e,
                   Length(mu)-i+1,(Length(mu)+1) mod 2);
        fi;
      elif mu[i]=2 then  ## the remaining case of 'staircase' partitions
        Add(mu,1); return RestrictedModule(H.operations.P.e2p0S(mu),H.e,
                                    [(Length(mu)+1)mod 2]);
      fi;
      return false;
    end;  # P.e2p0S()
  fi;

  return H;
end;  ## end of Specht()

#M The record returned by Schur() is essentially identical to that returned
#M by Specht(). The only differences are superficial; namely, the functions
#M are called H.W, H.P, and H.F rather than H.S, H.P, and H.D and 
#M H.IsSpecht is false rather than true. Even though these names are 
#M different the value of S.X(x).module is still "S", "P", or "D" and 
#M H.operations() contains ## records S, P, and D respectively. 
#M Consequently if you want to create either Specht() or Schur() 'modules'
#M in a function it is best not to do this via S.X(*) but instead to use 
#M the more explicit H.operations.New(*) if you want to these functions to
#M  work for both the Hecke and q-Schur algebras.
Schur:=function(arg) local schur;
  schur:=Specht(arg);
  schur.IsSpecht:=false;
  schur.W:=schur.S; #Unbind(schur.S);
  if IsBound(schur.Sq) then schur.Wq:=schur.Sq; fi; #Unbind(schur.Sq); fi;
  schur.F:=schur.D; #Unbind(schur.D);
  return schur;
end;

